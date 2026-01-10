"""
VMS Audio Cleaner - Core Implementation
=======================================

Primary audio cleaning class with full processing pipeline.
Adapted from LIGO gravitational wave signal extraction algorithms.
"""

import numpy as np
from scipy import signal
from scipy.io import wavfile
from dataclasses import dataclass
from typing import Optional, Callable, Tuple, List
import warnings

warnings.filterwarnings('ignore')

# Optional high-quality audio I/O
try:
    import soundfile as sf
    SOUNDFILE_AVAILABLE = True
except ImportError:
    SOUNDFILE_AVAILABLE = False


@dataclass
class CleaningResult:
    """Result of audio cleaning operation."""
    audio: np.ndarray           # Cleaned audio
    sample_rate: int            # Sample rate
    input_snr: float            # Estimated input SNR (dB)
    output_snr: float           # Estimated output SNR (dB)
    improvement: float          # SNR improvement (dB)
    correlation: float          # Correlation with estimated clean signal
    processing_time: float      # Time taken (seconds)
    frames_processed: int       # Number of frames processed


class VMSAudioCleaner:
    """
    Vibrational Mode Separator for Audio Cleaning.

    This algorithm extracts coherent voice modes from noisy audio
    using techniques derived from gravitational wave analysis.

    The key innovation is using hexagonal frequency relationships
    to better preserve voice harmonics while suppressing noise.
    """

    # Voice frequency ranges
    VOICE_FUNDAMENTAL = (80, 300)    # F0 range (fundamental)
    VOICE_HARMONICS = (300, 4000)    # Harmonics and formants
    VOICE_FULL = (80, 8000)          # Full voice spectrum

    # Harmonic ratios (natural + OCTH-enhanced)
    VOICE_RATIOS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    ENHANCED_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512,
        noise_floor_percentile: float = 25,
        harmonic_weight: float = 0.6,
        spectral_smoothing: float = 0.4,
        gain_floor: float = 0.15
    ):
        """
        Initialize VMS Audio Cleaner.

        Args:
            sample_rate: Audio sample rate in Hz
            frame_size: FFT frame size (larger = better frequency resolution)
            hop_size: Samples between frames (smaller = smoother output)
            noise_floor_percentile: Percentile for noise estimation
            harmonic_weight: Weight for harmonic preservation (0-1)
            spectral_smoothing: Smoothing factor for spectral mask
            gain_floor: Minimum gain to prevent complete silence
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size
        self.noise_floor_percentile = noise_floor_percentile
        self.harmonic_weight = harmonic_weight
        self.spectral_smoothing = spectral_smoothing
        self.gain_floor = gain_floor

        # Precompute frequency bins (for rfft output)
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Analysis/synthesis windows (sqrt-Hann for WOLA)
        self.analysis_window = np.sqrt(np.hanning(frame_size))
        self.synthesis_window = np.sqrt(np.hanning(frame_size))
        self.wola_window = np.hanning(frame_size)  # For normalization

        # Noise profile (updated during processing)
        self.noise_profile: Optional[np.ndarray] = None

    def estimate_noise_profile(
        self,
        audio: np.ndarray,
        noise_duration: float = 0.5
    ) -> np.ndarray:
        """
        Estimate noise profile from initial segment of audio.

        In real-time mode, this would be continuously updated.

        Args:
            audio: Input audio signal
            noise_duration: Duration of noise-only segment (seconds)

        Returns:
            Estimated noise spectrum
        """
        n_samples = int(noise_duration * self.sample_rate)
        noise_segment = audio[:min(n_samples, len(audio))]

        # Compute average spectrum of noise
        n_frames = max(1, (len(noise_segment) - self.frame_size) // self.hop_size)

        spectra = []
        for i in range(n_frames):
            start = i * self.hop_size
            frame = noise_segment[start:start + self.frame_size]
            if len(frame) < self.frame_size:
                frame = np.pad(frame, (0, self.frame_size - len(frame)))

            windowed = frame * self.analysis_window
            spectrum = np.abs(np.fft.rfft(windowed))
            spectra.append(spectrum)

        if spectra:
            self.noise_profile = np.percentile(
                spectra, self.noise_floor_percentile, axis=0
            )
        else:
            self.noise_profile = np.zeros(self.frame_size // 2 + 1)

        return self.noise_profile

    def detect_fundamental(self, spectrum: np.ndarray) -> float:
        """
        Detect fundamental frequency (F0) using peak detection.

        Args:
            spectrum: Amplitude spectrum

        Returns:
            Detected F0 in Hz
        """
        f_min, f_max = self.VOICE_FUNDAMENTAL
        freq_resolution = self.sample_rate / self.frame_size

        idx_min = int(f_min / freq_resolution)
        idx_max = min(int(f_max / freq_resolution), len(spectrum) - 1)

        if idx_max <= idx_min or idx_max >= len(spectrum):
            return 150.0  # Default F0

        voice_spectrum = spectrum[idx_min:idx_max + 1]
        if len(voice_spectrum) == 0:
            return 150.0

        peak_idx = np.argmax(voice_spectrum)
        f0 = (idx_min + peak_idx) * freq_resolution

        return max(f_min, min(f_max, f0))

    def create_harmonic_mask(
        self,
        spectrum: np.ndarray,
        f0: float
    ) -> np.ndarray:
        """
        Create a mask that preserves harmonic frequencies.

        Uses OCTH-inspired frequency relationships for better
        voice naturalness.

        Args:
            spectrum: Amplitude spectrum
            f0: Fundamental frequency

        Returns:
            Harmonic preservation mask
        """
        mask = np.zeros(len(spectrum))
        freq_resolution = self.sample_rate / self.frame_size

        # Combine natural harmonics with OCTH-enhanced ratios
        all_ratios = sorted(set(self.VOICE_RATIOS + self.ENHANCED_RATIOS))

        for ratio in all_ratios:
            f_harmonic = f0 * ratio

            if f_harmonic > self.sample_rate / 2:
                break

            idx = int(f_harmonic / freq_resolution)
            if idx >= len(mask):
                break

            # Width proportional to frequency
            width = max(2, int(f_harmonic / 50))

            for i in range(max(0, idx - width), min(len(mask), idx + width + 1)):
                distance = abs(i - idx)
                mask[i] = max(mask[i], np.exp(-distance**2 / (2 * (width/2)**2)))

        return mask

    def spectral_subtraction(
        self,
        spectrum: np.ndarray,
        phase: np.ndarray,
        noise_estimate: np.ndarray
    ) -> np.ndarray:
        """
        Advanced spectral subtraction with harmonic preservation.

        Uses Wiener filtering approach for more natural sound.

        Args:
            spectrum: Input amplitude spectrum
            phase: Input phase spectrum
            noise_estimate: Estimated noise spectrum

        Returns:
            Cleaned amplitude spectrum
        """
        # Avoid division by zero
        spectrum_safe = np.maximum(spectrum, 1e-10)
        noise_safe = np.maximum(noise_estimate, 1e-10)

        # Power-domain Wiener filtering
        signal_power = spectrum_safe ** 2
        noise_power = noise_safe ** 2

        # Wiener gain
        gain_sq = np.maximum(signal_power - noise_power, 0) / signal_power
        gain = np.sqrt(gain_sq)
        gain = np.maximum(gain, self.gain_floor)

        # Detect F0 and create harmonic mask
        f0 = self.detect_fundamental(spectrum)
        harmonic_mask = self.create_harmonic_mask(spectrum, f0)

        # Apply gain
        subtracted = spectrum * gain

        # Boost harmonics that might have been over-suppressed
        harmonic_boost = 1 + (1 - gain) * harmonic_mask * self.harmonic_weight
        final_spectrum = subtracted * harmonic_boost

        # Don't exceed original
        final_spectrum = np.minimum(final_spectrum, spectrum)

        return final_spectrum

    def process_frame(self, frame: np.ndarray) -> np.ndarray:
        """
        Process a single audio frame.

        Args:
            frame: Input frame (frame_size samples)

        Returns:
            Processed frame (windowed for overlap-add)
        """
        if len(frame) < self.frame_size:
            frame = np.pad(frame, (0, self.frame_size - len(frame)))

        # Analysis: window and FFT
        windowed = frame * self.analysis_window
        fft_data = np.fft.rfft(windowed)
        spectrum = np.abs(fft_data)
        phase = np.angle(fft_data)

        # Get noise estimate
        if self.noise_profile is None:
            noise_estimate = np.full_like(
                spectrum,
                np.percentile(spectrum, self.noise_floor_percentile)
            )
        else:
            noise_estimate = self.noise_profile

        # Clean spectrum
        cleaned_spectrum = self.spectral_subtraction(spectrum, phase, noise_estimate)

        # Reconstruct with original phase
        reconstructed = cleaned_spectrum * np.exp(1j * phase)

        # Inverse FFT
        cleaned_frame = np.fft.irfft(reconstructed, n=self.frame_size)

        # Synthesis window
        cleaned_frame *= self.synthesis_window

        return cleaned_frame

    def clean_audio(
        self,
        audio: np.ndarray,
        noise_sample: Optional[np.ndarray] = None,
        progress_callback: Optional[Callable[[float], None]] = None
    ) -> CleaningResult:
        """
        Clean a complete audio signal.

        Args:
            audio: Input audio (mono, normalized to [-1, 1])
            noise_sample: Optional noise-only sample for better estimation
            progress_callback: Optional callback(progress_percent)

        Returns:
            CleaningResult with cleaned audio and metrics
        """
        import time
        start_time = time.time()

        # Normalize input
        audio = np.asarray(audio, dtype=np.float32)
        original_max = np.max(np.abs(audio))
        if original_max > 1.0:
            audio = audio / original_max

        # Estimate noise profile
        if noise_sample is not None:
            self.estimate_noise_profile(noise_sample)
        else:
            self.estimate_noise_profile(audio, noise_duration=0.5)

        # Output buffer with overlap-add
        output = np.zeros(len(audio) + self.frame_size)
        window_sum = np.zeros(len(audio) + self.frame_size)

        # Process frames
        n_frames = (len(audio) - self.frame_size) // self.hop_size + 1

        for i in range(n_frames):
            start = i * self.hop_size
            frame = audio[start:start + self.frame_size]

            cleaned_frame = self.process_frame(frame)

            # Overlap-add
            output[start:start + self.frame_size] += cleaned_frame
            window_sum[start:start + self.frame_size] += self.wola_window

            if progress_callback and i % 100 == 0:
                progress_callback(100 * i / n_frames)

        # WOLA normalization
        min_window_sum = np.max(window_sum) * 0.5
        valid_mask = window_sum >= min_window_sum

        output_normalized = np.zeros_like(output)
        output_normalized[valid_mask] = output[valid_mask] / window_sum[valid_mask]

        # Handle edges
        edge_mask = ~valid_mask & (window_sum > 1e-8)
        if np.any(edge_mask):
            expected_sum = np.max(window_sum)
            output_normalized[edge_mask] = output[edge_mask] / expected_sum

        # Trim to original length
        output_final = output_normalized[:len(audio)]

        # Prevent clipping
        output_max = np.max(np.abs(output_final))
        if output_max > 0.95:
            output_final = output_final / output_max * 0.95

        processing_time = time.time() - start_time

        # Calculate metrics (skip edges)
        margin = min(10000, len(audio) // 10)
        start_idx = margin
        end_idx = len(audio) - margin

        if end_idx > start_idx:
            input_segment = audio[start_idx:end_idx]
            output_segment = output_final[start_idx:end_idx]

            # Rough SNR estimation
            input_power = np.mean(input_segment ** 2)
            output_power = np.mean(output_segment ** 2)
            noise_reduction = np.mean((input_segment - output_segment) ** 2)

            input_snr = 10 * np.log10(input_power / max(noise_reduction, 1e-10))
            output_snr = 10 * np.log10(output_power / max(noise_reduction, 1e-10))

            correlation = np.corrcoef(input_segment, output_segment)[0, 1]
        else:
            input_snr = 0.0
            output_snr = 0.0
            correlation = 1.0

        return CleaningResult(
            audio=output_final,
            sample_rate=self.sample_rate,
            input_snr=input_snr,
            output_snr=output_snr,
            improvement=output_snr - input_snr,
            correlation=correlation,
            processing_time=processing_time,
            frames_processed=n_frames
        )

    @staticmethod
    def load_audio(filepath: str) -> Tuple[np.ndarray, int]:
        """Load audio file, returns (audio, sample_rate)."""
        if SOUNDFILE_AVAILABLE:
            audio, sr = sf.read(filepath)
        else:
            sr, audio = wavfile.read(filepath)
            if audio.dtype == np.int16:
                audio = audio.astype(np.float32) / 32768.0
            elif audio.dtype == np.int32:
                audio = audio.astype(np.float32) / 2147483648.0

        # Convert to mono if stereo
        if len(audio.shape) > 1:
            audio = np.mean(audio, axis=1)

        return audio, sr

    @staticmethod
    def save_audio(filepath: str, audio: np.ndarray, sample_rate: int):
        """Save audio to file."""
        if SOUNDFILE_AVAILABLE:
            sf.write(filepath, audio, sample_rate)
        else:
            audio_int = (audio * 32767).astype(np.int16)
            wavfile.write(filepath, sample_rate, audio_int)
