"""
VMS Spectral Analysis Tools
===========================

Low-level spectral processing functions for VMS engine.
These are the building blocks for audio cleaning and
signal analysis.
"""

import numpy as np
from typing import Tuple, Optional, List
from dataclasses import dataclass


@dataclass
class SpectralFrame:
    """Represents a single spectral frame."""
    magnitude: np.ndarray      # Amplitude spectrum
    phase: np.ndarray          # Phase spectrum
    frequencies: np.ndarray    # Frequency bins
    time: float                # Time position in seconds
    sample_rate: int           # Sample rate


class SpectralAnalyzer:
    """
    Spectral analysis toolkit for VMS engine.

    Provides methods for spectral decomposition, peak detection,
    and frequency analysis.
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512,
        window_type: str = 'hann'
    ):
        """
        Initialize spectral analyzer.

        Args:
            sample_rate: Audio sample rate
            frame_size: FFT frame size
            hop_size: Samples between frames
            window_type: Window function ('hann', 'hamming', 'blackman')
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size

        # Create window
        if window_type == 'hann':
            self.window = np.hanning(frame_size)
        elif window_type == 'hamming':
            self.window = np.hamming(frame_size)
        elif window_type == 'blackman':
            self.window = np.blackman(frame_size)
        else:
            self.window = np.ones(frame_size)

        # Frequency bins
        self.frequencies = np.fft.rfftfreq(frame_size, 1/sample_rate)

    def analyze_frame(
        self,
        frame: np.ndarray,
        time_pos: float = 0.0
    ) -> SpectralFrame:
        """
        Analyze a single audio frame.

        Args:
            frame: Audio samples (frame_size)
            time_pos: Time position in seconds

        Returns:
            SpectralFrame with magnitude and phase
        """
        if len(frame) < self.frame_size:
            frame = np.pad(frame, (0, self.frame_size - len(frame)))

        windowed = frame[:self.frame_size] * self.window
        fft_data = np.fft.rfft(windowed)

        return SpectralFrame(
            magnitude=np.abs(fft_data),
            phase=np.angle(fft_data),
            frequencies=self.frequencies,
            time=time_pos,
            sample_rate=self.sample_rate
        )

    def analyze_audio(self, audio: np.ndarray) -> List[SpectralFrame]:
        """
        Analyze complete audio signal.

        Args:
            audio: Audio signal

        Returns:
            List of SpectralFrame objects
        """
        frames = []
        n_frames = (len(audio) - self.frame_size) // self.hop_size + 1

        for i in range(n_frames):
            start = i * self.hop_size
            frame = audio[start:start + self.frame_size]
            time_pos = start / self.sample_rate

            frames.append(self.analyze_frame(frame, time_pos))

        return frames

    def get_spectrogram(
        self,
        audio: np.ndarray,
        db_scale: bool = True
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute spectrogram.

        Args:
            audio: Audio signal
            db_scale: Convert to dB scale

        Returns:
            (spectrogram, frequencies, times)
        """
        frames = self.analyze_audio(audio)

        if not frames:
            return np.array([[]]), self.frequencies, np.array([])

        spectrogram = np.vstack([f.magnitude for f in frames]).T
        times = np.array([f.time for f in frames])

        if db_scale:
            spectrogram = 20 * np.log10(np.maximum(spectrogram, 1e-10))

        return spectrogram, self.frequencies, times

    def find_peaks(
        self,
        spectrum: np.ndarray,
        n_peaks: int = 10,
        min_height: float = None,
        min_separation_hz: float = 50.0
    ) -> List[Tuple[float, float]]:
        """
        Find spectral peaks.

        Args:
            spectrum: Amplitude spectrum
            n_peaks: Maximum number of peaks
            min_height: Minimum peak height (auto if None)
            min_separation_hz: Minimum Hz between peaks

        Returns:
            List of (frequency, amplitude) tuples
        """
        from scipy.signal import find_peaks as scipy_peaks

        if min_height is None:
            min_height = np.median(spectrum) * 2

        freq_res = self.sample_rate / self.frame_size
        distance = int(min_separation_hz / freq_res)

        peak_idx, _ = scipy_peaks(
            spectrum,
            height=min_height,
            distance=max(1, distance)
        )

        if len(peak_idx) == 0:
            return []

        # Sort by amplitude
        amplitudes = spectrum[peak_idx]
        sorted_idx = np.argsort(amplitudes)[::-1][:n_peaks]

        return [
            (self.frequencies[peak_idx[i]], spectrum[peak_idx[i]])
            for i in sorted_idx
        ]


def wiener_filter(
    spectrum: np.ndarray,
    noise_spectrum: np.ndarray,
    alpha: float = 1.0,
    gain_floor: float = 0.1
) -> np.ndarray:
    """
    Apply Wiener filter to spectrum.

    Wiener filter: G(f) = max(|S|^2 - α|N|^2, 0) / |S|^2

    Args:
        spectrum: Input amplitude spectrum
        noise_spectrum: Noise amplitude spectrum
        alpha: Noise subtraction factor (1.0 = full subtraction)
        gain_floor: Minimum gain to prevent complete silence

    Returns:
        Filtered spectrum
    """
    spectrum_safe = np.maximum(spectrum, 1e-10)
    noise_safe = np.maximum(noise_spectrum, 1e-10)

    signal_power = spectrum_safe ** 2
    noise_power = noise_safe ** 2

    # Wiener gain
    gain_sq = np.maximum(signal_power - alpha * noise_power, 0) / signal_power
    gain = np.sqrt(gain_sq)
    gain = np.maximum(gain, gain_floor)

    return spectrum * gain


def harmonic_preserving_filter(
    spectrum: np.ndarray,
    frequencies: np.ndarray,
    f0: float,
    harmonic_ratios: List[float] = None,
    harmonic_weight: float = 0.5,
    bandwidth_factor: float = 0.05
) -> np.ndarray:
    """
    Create harmonic-preserving filter mask.

    Boosts frequencies at harmonic positions relative to f0.

    Args:
        spectrum: Input spectrum (for reference, not modified)
        frequencies: Frequency array
        f0: Fundamental frequency
        harmonic_ratios: Ratios to preserve (default: 1-8)
        harmonic_weight: Weight for harmonic boost (0-1)
        bandwidth_factor: Fraction of frequency for bandwidth

    Returns:
        Harmonic boost mask (multiply with gain)
    """
    if harmonic_ratios is None:
        harmonic_ratios = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

    mask = np.ones(len(spectrum))

    for ratio in harmonic_ratios:
        f_harm = f0 * ratio

        if f_harm > frequencies[-1]:
            break

        # Gaussian centered on harmonic
        sigma = f_harm * bandwidth_factor
        gaussian = np.exp(-0.5 * ((frequencies - f_harm) / sigma) ** 2)

        # Add to mask
        mask += gaussian * harmonic_weight

    return mask


def spectral_gate(
    spectrum: np.ndarray,
    threshold: float,
    attack_ratio: float = 0.5,
    hold_bins: int = 2
) -> np.ndarray:
    """
    Apply spectral gate (noise gate in frequency domain).

    Args:
        spectrum: Input amplitude spectrum
        threshold: Gate threshold
        attack_ratio: Soft-knee ratio (0=hard, 1=soft)
        hold_bins: Bins to hold gate open

    Returns:
        Gated spectrum
    """
    gain = np.ones(len(spectrum))

    # Hard gate positions
    below_threshold = spectrum < threshold

    if attack_ratio > 0:
        # Soft knee
        for i in range(len(spectrum)):
            if below_threshold[i]:
                ratio = spectrum[i] / threshold
                gain[i] = ratio ** (1 / attack_ratio)
    else:
        # Hard gate
        gain[below_threshold] = 0.0

    # Hold: don't gate if neighbors are loud
    if hold_bins > 0:
        for i in range(len(spectrum)):
            if gain[i] < 1.0:
                # Check neighbors
                start = max(0, i - hold_bins)
                end = min(len(spectrum), i + hold_bins + 1)
                if np.any(gain[start:end] == 1.0):
                    gain[i] = 1.0

    return spectrum * gain


def estimate_noise_floor(
    spectrogram: np.ndarray,
    percentile: float = 10,
    smoothing: int = 5
) -> np.ndarray:
    """
    Estimate noise floor from spectrogram.

    Args:
        spectrogram: 2D spectrogram (freq x time)
        percentile: Percentile for noise estimation
        smoothing: Smoothing window size

    Returns:
        Noise floor estimate per frequency bin
    """
    # Take percentile across time
    noise_floor = np.percentile(spectrogram, percentile, axis=1)

    # Smooth
    if smoothing > 1:
        from scipy.ndimage import uniform_filter1d
        noise_floor = uniform_filter1d(noise_floor, smoothing)

    return noise_floor


def detect_voice_activity(
    spectrum: np.ndarray,
    frequencies: np.ndarray,
    voice_range: Tuple[float, float] = (80, 4000),
    threshold_db: float = -30
) -> bool:
    """
    Simple voice activity detection.

    Args:
        spectrum: Amplitude spectrum
        frequencies: Frequency array
        voice_range: Voice frequency range
        threshold_db: Detection threshold in dB

    Returns:
        True if voice activity detected
    """
    # Mask for voice frequencies
    mask = (frequencies >= voice_range[0]) & (frequencies <= voice_range[1])

    if not np.any(mask):
        return False

    # Energy in voice band
    voice_energy = np.mean(spectrum[mask] ** 2)

    # Total energy
    total_energy = np.mean(spectrum ** 2)

    if total_energy < 1e-10:
        return False

    # Voice band ratio
    ratio_db = 10 * np.log10(voice_energy / total_energy)

    return ratio_db > threshold_db
