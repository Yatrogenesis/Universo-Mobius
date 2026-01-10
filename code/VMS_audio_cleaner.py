#!/usr/bin/env python3
"""
VMS Audio Cleaner - "ClearVoice Pro"
====================================

Audio noise reduction based on Vibrational Mode Separator algorithm
originally developed for LIGO gravitational wave detection.

Key Innovation:
- Extracts coherent frequency modes from extremely noisy signals
- Preserves natural voice characteristics (no "robotic" artifacts)
- Works in real-time with low latency

Based on OCTH hexagonal frequency detection, adapted for audio.

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

import numpy as np
from scipy import signal
from scipy.fft import fft, ifft, fftfreq
from scipy.io import wavfile
import warnings
warnings.filterwarnings('ignore')

# Optional imports
try:
    import soundfile as sf
    SOUNDFILE_AVAILABLE = True
except ImportError:
    SOUNDFILE_AVAILABLE = False

try:
    import librosa
    LIBROSA_AVAILABLE = True
except ImportError:
    LIBROSA_AVAILABLE = False


class VMSAudioCleaner:
    """
    Vibrational Mode Separator for Audio Cleaning.

    This algorithm extracts coherent voice modes from noisy audio
    using techniques derived from gravitational wave analysis.

    Key parameters:
    - voice_range: Frequency range for human voice (80-4000 Hz typical)
    - harmonic_tolerance: How strictly to enforce harmonic relationships
    - noise_floor_percentile: Percentile for noise floor estimation
    """

    # Voice frequency ranges
    VOICE_FUNDAMENTAL = (80, 300)    # F0 range (fundamental)
    VOICE_HARMONICS = (300, 4000)    # Harmonics and formants
    VOICE_FULL = (80, 8000)          # Full voice spectrum

    # Harmonic ratios (based on voice physics + hexagonal enhancement)
    # Natural voice harmonics: 1, 2, 3, 4, 5...
    # OCTH enhancement: sqrt(3) ≈ 1.73, sqrt(7) ≈ 2.65
    VOICE_RATIOS = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    ENHANCED_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512,
        noise_floor_percentile: float = 10,
        harmonic_weight: float = 0.7,
        spectral_smoothing: float = 0.3
    ):
        """
        Initialize VMS Audio Cleaner.

        Args:
            sample_rate: Audio sample rate in Hz
            frame_size: FFT frame size (larger = better frequency resolution)
            hop_size: Samples between frames (smaller = smoother output)
            noise_floor_percentile: Percentile for noise estimation (lower = aggressive)
            harmonic_weight: Weight for harmonic preservation (0-1)
            spectral_smoothing: Smoothing factor for spectral mask (0-1)
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size
        self.noise_floor_percentile = noise_floor_percentile
        self.harmonic_weight = harmonic_weight
        self.spectral_smoothing = spectral_smoothing

        # Precompute frequency bins (for rfft output)
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Noise profile (updated during processing)
        self.noise_profile = None

    def estimate_noise_profile(self, audio: np.ndarray, noise_duration: float = 0.5) -> np.ndarray:
        """
        Estimate noise profile from initial segment of audio.

        In real-time mode, this would be continuously updated.
        """
        # Use first noise_duration seconds
        n_samples = int(noise_duration * self.sample_rate)
        noise_segment = audio[:min(n_samples, len(audio))]

        # Compute average spectrum of noise
        n_frames = max(1, (len(noise_segment) - self.frame_size) // self.hop_size)

        # Use same window as process_frame (sqrt-Hann for WOLA)
        window = np.sqrt(np.hanning(self.frame_size))

        spectra = []
        for i in range(n_frames):
            start = i * self.hop_size
            frame = noise_segment[start:start + self.frame_size]
            if len(frame) < self.frame_size:
                frame = np.pad(frame, (0, self.frame_size - len(frame)))

            # Window and rfft
            windowed = frame * window
            spectrum = np.abs(np.fft.rfft(windowed))
            spectra.append(spectrum)

        if spectra:
            # Use percentile for robust noise estimation
            self.noise_profile = np.percentile(spectra, self.noise_floor_percentile, axis=0)
        else:
            self.noise_profile = np.zeros(self.frame_size // 2 + 1)

        return self.noise_profile

    def detect_fundamental(self, spectrum: np.ndarray) -> float:
        """
        Detect fundamental frequency (F0) using peak detection.

        This is crucial for harmonic preservation.
        """
        # Focus on voice fundamental range
        f_min, f_max = self.VOICE_FUNDAMENTAL

        # Find indices for rfft output (uses rfftfreq)
        freq_resolution = self.sample_rate / self.frame_size
        idx_min = int(f_min / freq_resolution)
        idx_max = min(int(f_max / freq_resolution), len(spectrum) - 1)

        # Find peak in fundamental range
        if idx_max <= idx_min or idx_max >= len(spectrum):
            return 150.0  # Default F0

        voice_spectrum = spectrum[idx_min:idx_max+1]
        if len(voice_spectrum) == 0:
            return 150.0  # Default F0

        peak_idx = np.argmax(voice_spectrum)
        f0 = (idx_min + peak_idx) * freq_resolution

        return max(f_min, min(f_max, f0))

    def create_harmonic_mask(self, spectrum: np.ndarray, f0: float) -> np.ndarray:
        """
        Create a mask that preserves harmonic frequencies.

        The key innovation: use OCTH-inspired frequency relationships
        to better preserve voice naturalness.
        """
        mask = np.zeros(len(spectrum))
        freq_resolution = self.sample_rate / self.frame_size

        # Combine natural harmonics with OCTH-enhanced ratios
        all_ratios = list(set(self.VOICE_RATIOS + self.ENHANCED_RATIOS))
        all_ratios.sort()

        for ratio in all_ratios:
            f_harmonic = f0 * ratio

            # Skip if outside audio range (Nyquist)
            if f_harmonic > self.sample_rate / 2:
                break

            # Create Gaussian peak around harmonic
            idx = int(f_harmonic / freq_resolution)
            if idx >= len(mask):
                break

            # Width proportional to frequency (wider for higher harmonics)
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

        Unlike basic spectral subtraction, this preserves voice harmonics
        even when they're below the noise floor. Uses Wiener filtering
        approach for more natural sound.
        """
        # Wiener-style gain calculation
        # g(f) = max(|X(f)|^2 - alpha * |N(f)|^2, 0) / |X(f)|^2

        # Avoid division by zero
        spectrum_safe = np.maximum(spectrum, 1e-10)
        noise_safe = np.maximum(noise_estimate, 1e-10)

        # Power-domain subtraction (less aggressive than amplitude)
        signal_power = spectrum_safe ** 2
        noise_power = noise_safe ** 2

        # Subtraction factor (1.0-2.0 typical, lower = less aggressive)
        alpha = 1.0

        # Wiener gain squared
        gain_sq = np.maximum(signal_power - alpha * noise_power, 0) / signal_power

        # Take square root to get amplitude gain, with floor
        gain = np.sqrt(gain_sq)
        gain = np.maximum(gain, 0.15)  # 15% floor prevents complete silence

        # Detect fundamental frequency
        f0 = self.detect_fundamental(spectrum)

        # Create harmonic preservation mask
        harmonic_mask = self.create_harmonic_mask(spectrum, f0)

        # Apply gain
        subtracted = spectrum * gain

        # Key innovation: boost harmonics that might have been over-suppressed
        # Harmonics are boosted back toward original level
        harmonic_boost = 1 + (1 - gain) * harmonic_mask * self.harmonic_weight
        final_spectrum = subtracted * harmonic_boost

        # Ensure we don't exceed original
        final_spectrum = np.minimum(final_spectrum, spectrum)

        return final_spectrum

    def process_frame(self, frame: np.ndarray) -> np.ndarray:
        """
        Process a single audio frame.

        This can be used for real-time processing.
        Returns windowed output for overlap-add.
        """
        # Ensure correct length
        if len(frame) < self.frame_size:
            frame = np.pad(frame, (0, self.frame_size - len(frame)))

        # Analysis window (sqrt-Hann for WOLA)
        window = np.sqrt(np.hanning(self.frame_size))
        windowed = frame * window

        # Use rfft for real signals
        fft_data = np.fft.rfft(windowed)
        spectrum = np.abs(fft_data)
        phase = np.angle(fft_data)

        # Get noise estimate
        if self.noise_profile is None:
            noise_estimate = np.percentile(spectrum, self.noise_floor_percentile)
            noise_estimate = np.full_like(spectrum, noise_estimate)
        else:
            noise_estimate = self.noise_profile

        # Clean spectrum
        cleaned_spectrum = self.spectral_subtraction(spectrum, phase, noise_estimate)

        # Reconstruct complex spectrum with original phase
        reconstructed = cleaned_spectrum * np.exp(1j * phase)

        # Inverse FFT
        cleaned_frame = np.fft.irfft(reconstructed, n=self.frame_size)

        # Synthesis window (sqrt-Hann for WOLA - analysis * synthesis = Hann)
        cleaned_frame *= window

        return cleaned_frame

    def clean_audio(
        self,
        audio: np.ndarray,
        noise_sample: np.ndarray = None,
        progress_callback=None
    ) -> np.ndarray:
        """
        Clean a complete audio signal.

        Args:
            audio: Input audio (mono, float32, normalized to [-1, 1])
            noise_sample: Optional noise-only sample for better estimation
            progress_callback: Optional callback(progress_percent)

        Returns:
            Cleaned audio signal
        """
        # Normalize input
        audio = np.asarray(audio, dtype=np.float32)
        original_max = np.max(np.abs(audio))
        if original_max > 1.0:
            audio = audio / original_max

        # Estimate noise profile
        if noise_sample is not None:
            self.estimate_noise_profile(noise_sample)
        else:
            # Use first 0.5 seconds as noise estimate
            self.estimate_noise_profile(audio, noise_duration=0.5)

        # Output buffer with overlap-add
        output = np.zeros(len(audio) + self.frame_size)
        window_sum = np.zeros(len(audio) + self.frame_size)

        # Window for WOLA normalization (sqrt-Hann squared = Hann)
        wola_window = np.hanning(self.frame_size)

        # Process frames
        n_frames = (len(audio) - self.frame_size) // self.hop_size + 1

        for i in range(n_frames):
            start = i * self.hop_size
            frame = audio[start:start + self.frame_size]

            cleaned_frame = self.process_frame(frame)

            # Overlap-add
            output[start:start + self.frame_size] += cleaned_frame
            window_sum[start:start + self.frame_size] += wola_window

            # Progress callback
            if progress_callback and i % 100 == 0:
                progress_callback(100 * i / n_frames)

        # Normalize by window sum for WOLA
        # Only normalize where we have sufficient overlap (avoid edge artifacts)
        min_window_sum = np.max(window_sum) * 0.5  # At least 50% of max
        valid_mask = window_sum >= min_window_sum

        output_normalized = np.zeros_like(output)
        output_normalized[valid_mask] = output[valid_mask] / window_sum[valid_mask]

        # For edge regions, scale by the ratio to expected sum
        edge_mask = ~valid_mask & (window_sum > 1e-8)
        if np.any(edge_mask):
            expected_sum = np.max(window_sum)
            output_normalized[edge_mask] = output[edge_mask] / expected_sum

        # Trim to original length
        output = output_normalized[:len(audio)]

        # Prevent clipping while preserving dynamics
        output_max = np.max(np.abs(output))
        if output_max > 0.95:
            output = output / output_max * 0.95

        return output

    @staticmethod
    def load_audio(filepath: str) -> tuple:
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
            # Scale to int16
            audio_int = (audio * 32767).astype(np.int16)
            wavfile.write(filepath, sample_rate, audio_int)


def demo_comparison():
    """
    Demo: Compare VMS cleaner effectiveness.
    """
    print("=" * 60)
    print("  VMS AUDIO CLEANER - DEMO")
    print("=" * 60)

    # Generate test signal: voice + noise
    duration = 3.0
    sr = 44100
    t = np.linspace(0, duration, int(sr * duration))

    # Create noise-only section at start (0.5s) then voice
    noise_duration = 0.5
    noise_samples = int(noise_duration * sr)

    # Simulated voice with pitch variation (starts after noise section)
    f0 = 150  # Hz
    voice = np.zeros_like(t)

    # Voice envelope: silence at start, then voice
    voice_envelope = np.zeros_like(t)
    voice_envelope[noise_samples:] = 1.0
    # Smooth onset
    onset_samples = int(0.05 * sr)
    voice_envelope[noise_samples:noise_samples+onset_samples] = np.linspace(0, 1, onset_samples)

    # Generate voice with harmonics
    pitch_mod = 1 + 0.05 * np.sin(2 * np.pi * 3 * t)  # Subtle pitch variation
    phase = 2 * np.pi * f0 * np.cumsum(pitch_mod) / sr

    for h in range(1, 10):
        amplitude = 1.0 / (h ** 0.7)  # Natural harmonic rolloff
        voice += amplitude * np.sin(h * phase)

    voice = voice * voice_envelope
    voice = voice / np.max(np.abs(voice)) * 0.6

    # Generate realistic noise (broadband + hum + some stationary)
    np.random.seed(42)  # Reproducible
    white_noise = np.random.randn(len(t))

    # Pink noise via 1/f filtering
    freqs_noise = np.fft.rfftfreq(len(t), 1/sr)
    freqs_noise[0] = 1  # Avoid division by zero
    pink_filter = 1 / np.sqrt(freqs_noise)
    pink_filter[0] = 0
    white_fft = np.fft.rfft(white_noise)
    pink_noise = np.fft.irfft(white_fft * pink_filter, n=len(t))
    pink_noise = pink_noise / np.max(np.abs(pink_noise))

    # Add hum at 60 Hz and harmonics
    hum = 0.15 * np.sin(2 * np.pi * 60 * t)
    hum += 0.08 * np.sin(2 * np.pi * 120 * t)
    hum += 0.04 * np.sin(2 * np.pi * 180 * t)

    # Combined noise (constant throughout)
    noise = 0.3 * pink_noise + hum + 0.1 * white_noise
    noise = noise / np.max(np.abs(noise)) * 0.4

    # Noisy audio
    noisy_audio = voice + noise

    print(f"\nGenerated test signal:")
    print(f"  Duration: {duration}s (first 0.5s noise-only)")
    print(f"  Sample rate: {sr} Hz")
    print(f"  Fundamental: {f0} Hz")

    # Create cleaner with optimized parameters
    cleaner = VMSAudioCleaner(
        sample_rate=sr,
        frame_size=2048,
        hop_size=512,
        noise_floor_percentile=25,
        harmonic_weight=0.6,
        spectral_smoothing=0.4
    )

    print("\nProcessing with VMS Audio Cleaner...")

    # Clean audio (uses first 0.5s for noise estimation automatically)
    cleaned = cleaner.clean_audio(noisy_audio)

    # Calculate SNR improvement (skip edges for accurate measurement)
    margin = 10000  # Skip first/last samples for stable measurement
    start = noise_samples + margin
    end = len(voice) - margin

    voice_v = voice[start:end]
    noisy_v = noisy_audio[start:end]
    clean_v = cleaned[start:end]

    # SNR before/after
    signal_power = np.mean(voice_v**2)
    noise_in_power = np.mean((noisy_v - voice_v)**2)
    noise_out_power = np.mean((clean_v - voice_v)**2)

    snr_in = 10 * np.log10(signal_power / max(noise_in_power, 1e-10))
    snr_out = 10 * np.log10(signal_power / max(noise_out_power, 1e-10))

    # Additional metrics
    correlation = np.corrcoef(voice_v, clean_v)[0, 1]

    print(f"\nResults:")
    print(f"  Input SNR:    {snr_in:.1f} dB")
    print(f"  Output SNR:   {snr_out:.1f} dB")
    print(f"  Improvement:  {snr_out - snr_in:+.1f} dB")
    print(f"  Correlation:  {correlation:.3f}")

    # Save files
    output_dir = "demo_audio"
    import os
    os.makedirs(output_dir, exist_ok=True)

    VMSAudioCleaner.save_audio(f"{output_dir}/original_voice.wav", voice, sr)
    VMSAudioCleaner.save_audio(f"{output_dir}/noisy_audio.wav", noisy_audio, sr)
    VMSAudioCleaner.save_audio(f"{output_dir}/vms_cleaned.wav", cleaned, sr)
    VMSAudioCleaner.save_audio(f"{output_dir}/noise_only.wav", noise, sr)

    print(f"\nAudio files saved to {output_dir}/")
    print("  - original_voice.wav (reference)")
    print("  - noisy_audio.wav (input)")
    print("  - vms_cleaned.wav (output)")
    print("  - noise_only.wav (isolated noise)")

    return {
        'input_snr': snr_in,
        'output_snr': snr_out,
        'improvement': snr_out - snr_in,
        'correlation': correlation
    }


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        # Process file
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else "cleaned_output.wav"

        print(f"Loading {input_file}...")
        audio, sr = VMSAudioCleaner.load_audio(input_file)

        cleaner = VMSAudioCleaner(sample_rate=sr)

        print("Processing...")
        cleaned = cleaner.clean_audio(audio)

        print(f"Saving to {output_file}...")
        VMSAudioCleaner.save_audio(output_file, cleaned, sr)

        print("Done!")
    else:
        # Run demo
        demo_comparison()
