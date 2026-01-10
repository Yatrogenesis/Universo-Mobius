"""
API Routes
==========

FastAPI routes for VMS audio processing services.
"""

from fastapi import APIRouter, UploadFile, File, HTTPException, Query, BackgroundTasks
from fastapi.responses import StreamingResponse, JSONResponse
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
import numpy as np
import io
import base64
import uuid
from datetime import datetime


router = APIRouter()


# =============================================================================
# MODELS
# =============================================================================

class CleaningOptions(BaseModel):
    """Audio cleaning configuration options."""
    noise_floor_percentile: float = Field(25, ge=1, le=50)
    harmonic_weight: float = Field(0.6, ge=0, le=1)
    gain_floor: float = Field(0.15, ge=0, le=0.5)
    frame_size: int = Field(2048, ge=512, le=8192)
    hop_size: int = Field(512, ge=128, le=2048)


class CleaningResult(BaseModel):
    """Audio cleaning result."""
    job_id: str
    status: str
    input_snr_db: Optional[float] = None
    output_snr_db: Optional[float] = None
    improvement_db: Optional[float] = None
    processing_time_ms: Optional[float] = None
    download_url: Optional[str] = None


class AnalysisResult(BaseModel):
    """Signal analysis result."""
    job_id: str
    duration_s: float
    sample_rate: int
    n_channels: int
    peak_amplitude: float
    rms_level_db: float
    spectral_centroid_hz: float
    is_hexagonal: bool
    hexagonal_confidence: float


class SpectrogramData(BaseModel):
    """Spectrogram data for visualization."""
    magnitude: List[List[float]]  # 2D array
    frequencies: List[float]
    times: List[float]
    sample_rate: int


class SignalBlock(BaseModel):
    """Signal block definition."""
    id: int
    freq_range: List[float]  # [min, max] Hz
    time_range: List[float]  # [start, end] seconds
    threshold_db: float
    label: str = ""
    is_signal: bool = True


class FilterConfig(BaseModel):
    """Filter configuration."""
    blocks: List[SignalBlock] = []
    preset: Optional[str] = None  # 'voice', 'music', 'hum'
    soft_edges: bool = True


# =============================================================================
# JOB STORAGE (In production, use Redis/database)
# =============================================================================

_jobs: Dict[str, Dict[str, Any]] = {}
_audio_cache: Dict[str, bytes] = {}


# =============================================================================
# HEALTH CHECK
# =============================================================================

@router.get("/", tags=["Health"])
async def root():
    """API root - health check."""
    return {
        "service": "VMS Audio Processing API",
        "version": "0.1.0",
        "status": "healthy",
        "timestamp": datetime.utcnow().isoformat()
    }


@router.get("/health", tags=["Health"])
async def health_check():
    """Detailed health check."""
    return {
        "status": "healthy",
        "components": {
            "api": "healthy",
            "vms_engine": "available",
            "discriminator": "available"
        },
        "timestamp": datetime.utcnow().isoformat()
    }


# =============================================================================
# AUDIO CLEANING
# =============================================================================

@router.post("/clean", response_model=CleaningResult, tags=["Audio Processing"])
async def clean_audio(
    file: UploadFile = File(...),
    background_tasks: BackgroundTasks = None,
    noise_floor_percentile: float = Query(25, ge=1, le=50),
    harmonic_weight: float = Query(0.6, ge=0, le=1),
):
    """
    Clean audio file using VMS algorithm.

    Removes noise while preserving voice characteristics using
    algorithms derived from gravitational wave detection.
    """
    # Validate file type
    if not file.filename.endswith(('.wav', '.mp3', '.flac', '.ogg')):
        raise HTTPException(400, "Unsupported audio format. Use WAV, MP3, FLAC, or OGG.")

    # Generate job ID
    job_id = str(uuid.uuid4())[:8]

    try:
        # Read audio file
        audio_bytes = await file.read()

        # Process audio
        result = await _process_audio_cleaning(
            audio_bytes,
            file.filename,
            noise_floor_percentile,
            harmonic_weight
        )

        # Store result
        _jobs[job_id] = result
        if result.get('audio_bytes'):
            _audio_cache[job_id] = result['audio_bytes']

        return CleaningResult(
            job_id=job_id,
            status="completed",
            input_snr_db=result.get('input_snr'),
            output_snr_db=result.get('output_snr'),
            improvement_db=result.get('improvement'),
            processing_time_ms=result.get('processing_time_ms'),
            download_url=f"/download/{job_id}"
        )

    except Exception as e:
        raise HTTPException(500, f"Processing error: {str(e)}")


async def _process_audio_cleaning(
    audio_bytes: bytes,
    filename: str,
    noise_floor_percentile: float,
    harmonic_weight: float
) -> Dict[str, Any]:
    """Process audio cleaning asynchronously."""
    import time
    from scipy.io import wavfile
    import io

    start_time = time.time()

    # Load audio
    audio_io = io.BytesIO(audio_bytes)

    try:
        sr, audio = wavfile.read(audio_io)
    except Exception:
        raise HTTPException(400, "Could not read audio file. Ensure it's a valid WAV file.")

    # Convert to float
    if audio.dtype == np.int16:
        audio = audio.astype(np.float32) / 32768.0
    elif audio.dtype == np.int32:
        audio = audio.astype(np.float32) / 2147483648.0

    # Convert to mono
    if len(audio.shape) > 1:
        audio = np.mean(audio, axis=1)

    # Process with VMS
    from ..vms_engine import VMSAudioCleaner

    cleaner = VMSAudioCleaner(
        sample_rate=sr,
        noise_floor_percentile=noise_floor_percentile,
        harmonic_weight=harmonic_weight
    )

    result = cleaner.clean_audio(audio)

    # Convert back to WAV bytes
    output_io = io.BytesIO()
    output_audio = (result.audio * 32767).astype(np.int16)
    wavfile.write(output_io, sr, output_audio)
    output_bytes = output_io.getvalue()

    processing_time = (time.time() - start_time) * 1000

    return {
        'input_snr': result.input_snr,
        'output_snr': result.output_snr,
        'improvement': result.improvement,
        'processing_time_ms': processing_time,
        'audio_bytes': output_bytes
    }


@router.get("/download/{job_id}", tags=["Audio Processing"])
async def download_cleaned_audio(job_id: str):
    """Download cleaned audio file."""
    if job_id not in _audio_cache:
        raise HTTPException(404, "Job not found or audio expired")

    audio_bytes = _audio_cache[job_id]

    return StreamingResponse(
        io.BytesIO(audio_bytes),
        media_type="audio/wav",
        headers={
            "Content-Disposition": f"attachment; filename=cleaned_{job_id}.wav"
        }
    )


# =============================================================================
# SIGNAL ANALYSIS
# =============================================================================

@router.post("/analyze", response_model=AnalysisResult, tags=["Analysis"])
async def analyze_audio(file: UploadFile = File(...)):
    """
    Analyze audio signal properties.

    Returns spectral characteristics and hexagonal pattern detection.
    """
    # Read audio
    audio_bytes = await file.read()

    try:
        from scipy.io import wavfile
        audio_io = io.BytesIO(audio_bytes)
        sr, audio = wavfile.read(audio_io)

        if audio.dtype == np.int16:
            audio = audio.astype(np.float32) / 32768.0

        if len(audio.shape) > 1:
            n_channels = audio.shape[1]
            audio_mono = np.mean(audio, axis=1)
        else:
            n_channels = 1
            audio_mono = audio

        # Basic analysis
        duration = len(audio_mono) / sr
        peak_amp = np.max(np.abs(audio_mono))
        rms = np.sqrt(np.mean(audio_mono ** 2))
        rms_db = 20 * np.log10(max(rms, 1e-10))

        # Spectral centroid
        spectrum = np.abs(np.fft.rfft(audio_mono[:sr]))  # First second
        freqs = np.fft.rfftfreq(sr, 1/sr)
        spectral_centroid = np.sum(freqs * spectrum) / max(np.sum(spectrum), 1e-10)

        # Hexagonal analysis
        from ..octh_core import FrequencyRatioAnalyzer
        analyzer = FrequencyRatioAnalyzer()
        hex_result = analyzer.analyze(spectrum, freqs, n_peaks=10)

        job_id = str(uuid.uuid4())[:8]

        return AnalysisResult(
            job_id=job_id,
            duration_s=duration,
            sample_rate=sr,
            n_channels=n_channels,
            peak_amplitude=float(peak_amp),
            rms_level_db=float(rms_db),
            spectral_centroid_hz=float(spectral_centroid),
            is_hexagonal=hex_result['is_hexagonal'],
            hexagonal_confidence=hex_result['hexagonal_fraction']
        )

    except Exception as e:
        raise HTTPException(500, f"Analysis error: {str(e)}")


@router.post("/spectrogram", response_model=SpectrogramData, tags=["Analysis"])
async def get_spectrogram(
    file: UploadFile = File(...),
    frame_size: int = Query(2048, ge=512, le=8192),
    hop_size: int = Query(512, ge=128, le=2048),
    max_freq: float = Query(8000, ge=100, le=22050)
):
    """
    Compute spectrogram for visualization.

    Returns frequency-time-amplitude data for 3D plotting.
    """
    audio_bytes = await file.read()

    try:
        from scipy.io import wavfile
        audio_io = io.BytesIO(audio_bytes)
        sr, audio = wavfile.read(audio_io)

        if audio.dtype == np.int16:
            audio = audio.astype(np.float32) / 32768.0

        if len(audio.shape) > 1:
            audio = np.mean(audio, axis=1)

        # Compute spectrogram
        from ..discriminator_3d import create_3d_spectrogram
        spec = create_3d_spectrogram(
            audio, sr, frame_size=frame_size, hop_size=hop_size
        )

        # Filter by max frequency
        freq_mask = spec.frequencies <= max_freq
        filtered_freqs = spec.frequencies[freq_mask]
        filtered_mag = spec.magnitude[freq_mask, :]

        return SpectrogramData(
            magnitude=filtered_mag.tolist(),
            frequencies=filtered_freqs.tolist(),
            times=spec.times.tolist(),
            sample_rate=sr
        )

    except Exception as e:
        raise HTTPException(500, f"Spectrogram error: {str(e)}")


# =============================================================================
# SIGNAL DISCRIMINATION
# =============================================================================

@router.post("/discriminate", tags=["Discrimination"])
async def discriminate_signal(
    file: UploadFile = File(...),
    config: FilterConfig = None
):
    """
    Apply signal discrimination using filter blocks.

    Separates signal from noise based on defined regions.
    """
    audio_bytes = await file.read()

    try:
        from scipy.io import wavfile
        audio_io = io.BytesIO(audio_bytes)
        sr, audio = wavfile.read(audio_io)

        if audio.dtype == np.int16:
            audio = audio.astype(np.float32) / 32768.0

        if len(audio.shape) > 1:
            audio = np.mean(audio, axis=1)

        # Create spectrogram
        from ..discriminator_3d import (
            create_3d_spectrogram,
            SignalDiscriminator,
            FilterSculptor
        )

        spec = create_3d_spectrogram(audio, sr)
        discriminator = SignalDiscriminator(spec)

        # Auto-detect or use provided config
        if config and config.blocks:
            sculptor = FilterSculptor()
            sculptor.initialize(spec.frequencies, spec.times)
            for block in config.blocks:
                sculptor.add_block(
                    time_range=tuple(block.time_range),
                    freq_range=tuple(block.freq_range),
                    threshold_db=block.threshold_db,
                    label=block.label,
                    is_signal=block.is_signal
                )
            mask = sculptor.generate_mask(spec.magnitude, soft_edges=config.soft_edges)
        elif config and config.preset:
            sculptor = FilterSculptor()
            sculptor.initialize(spec.frequencies, spec.times)
            if config.preset == 'voice':
                sculptor.create_preset_voice()
            elif config.preset == 'music':
                sculptor.create_preset_music()
            mask = sculptor.generate_mask(spec.magnitude)
        else:
            # Auto-detect
            discriminator.detect_signal_blocks()
            discriminator.detect_noise_regions()
            mask = discriminator.create_filter_mask()

        # Reconstruct audio
        filtered_audio = discriminator.reconstruct_audio(mask, audio)

        # Convert to WAV
        output_io = io.BytesIO()
        output_audio = (filtered_audio * 32767).astype(np.int16)
        wavfile.write(output_io, sr, output_audio)

        job_id = str(uuid.uuid4())[:8]
        _audio_cache[job_id] = output_io.getvalue()

        return {
            "job_id": job_id,
            "status": "completed",
            "download_url": f"/download/{job_id}",
            "blocks_detected": len(discriminator.signal_blocks),
            "noise_regions": len(discriminator.noise_regions)
        }

    except Exception as e:
        raise HTTPException(500, f"Discrimination error: {str(e)}")


# =============================================================================
# HEXAGONAL ANALYSIS
# =============================================================================

@router.post("/hexagonal", tags=["OCTH Analysis"])
async def hexagonal_analysis(file: UploadFile = File(...)):
    """
    Perform OCTH hexagonal frequency analysis.

    Detects hexagonal frequency patterns in the signal,
    similar to gravitational wave analysis.
    """
    audio_bytes = await file.read()

    try:
        from scipy.io import wavfile
        audio_io = io.BytesIO(audio_bytes)
        sr, audio = wavfile.read(audio_io)

        if audio.dtype == np.int16:
            audio = audio.astype(np.float32) / 32768.0

        if len(audio.shape) > 1:
            audio = np.mean(audio, axis=1)

        # Compute spectrum
        spectrum = np.abs(np.fft.rfft(audio))
        freqs = np.fft.rfftfreq(len(audio), 1/sr)

        # Hexagonal analysis
        from ..octh_core import (
            FrequencyRatioAnalyzer,
            HexagonalTransform,
            calculate_binomial_significance,
            HEXAGONAL_RATIOS
        )

        # Ratio analysis
        analyzer = FrequencyRatioAnalyzer()
        ratio_result = analyzer.analyze(spectrum, freqs)

        # Transform analysis
        transform = HexagonalTransform()
        f0 = transform.detect_fundamental(spectrum, freqs)
        modes = transform.find_modes(spectrum, freqs, f0)

        # Statistical significance
        if ratio_result['total_ratios'] > 0:
            sig_result = calculate_binomial_significance(
                ratio_result['hexagonal_count'],
                ratio_result['total_ratios'],
                null_probability=0.5
            )
        else:
            sig_result = None

        return {
            "fundamental_frequency_hz": float(f0),
            "hexagonal_ratios": HEXAGONAL_RATIOS,
            "detected_modes": [
                {
                    "frequency": m.frequency,
                    "ratio": m.ratio_value,
                    "expected_ratio": m.expected_ratio,
                    "deviation": m.deviation,
                    "snr": m.snr,
                    "confidence": m.confidence
                }
                for m in modes
            ],
            "ratio_analysis": {
                "hexagonal_count": ratio_result['hexagonal_count'],
                "total_ratios": ratio_result['total_ratios'],
                "hexagonal_fraction": ratio_result['hexagonal_fraction'],
                "is_hexagonal": ratio_result['is_hexagonal']
            },
            "statistical_significance": {
                "sigma": sig_result.sigma if sig_result else 0,
                "p_value": sig_result.p_value if sig_result else 1.0,
                "z_score": sig_result.z_score if sig_result else 0
            } if sig_result else None
        }

    except Exception as e:
        raise HTTPException(500, f"Hexagonal analysis error: {str(e)}")
