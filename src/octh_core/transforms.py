"""
OCTH Transforms and Frequency Analysis
=======================================

Core algorithms for detecting hexagonal frequency patterns in signals.
These transforms are the foundation for both gravitational wave
analysis and audio processing.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from .constants import HEXAGONAL_RATIOS, DEFAULT_TOLERANCE


@dataclass
class HexagonalMode:
    """Represents a detected hexagonal frequency mode."""
    frequency: float          # Hz
    ratio_index: int          # Which hexagonal ratio (0-5)
    ratio_value: float        # Actual ratio to fundamental
    expected_ratio: float     # Expected hexagonal ratio
    deviation: float          # Fractional deviation from expected
    amplitude: float          # Signal amplitude at this frequency
    snr: float               # Local signal-to-noise ratio
    confidence: float        # Detection confidence (0-1)


class HexagonalTransform:
    """
    Transform for detecting hexagonal frequency patterns.

    The hexagonal transform projects a frequency spectrum onto
    the space of hexagonal ratios, revealing coherent structures
    that emerge from hexagonal spacetime geometry.
    """

    def __init__(
        self,
        ratios: List[float] = None,
        tolerance: float = DEFAULT_TOLERANCE,
        min_snr: float = 3.0
    ):
        """
        Initialize hexagonal transform.

        Args:
            ratios: List of expected frequency ratios (default: HEXAGONAL_RATIOS)
            tolerance: Fractional tolerance for ratio matching
            min_snr: Minimum SNR for mode detection
        """
        self.ratios = ratios if ratios is not None else HEXAGONAL_RATIOS
        self.tolerance = tolerance
        self.min_snr = min_snr

    def detect_fundamental(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        freq_range: Tuple[float, float] = (10, 500)
    ) -> float:
        """
        Detect the fundamental frequency (F0) in a spectrum.

        Uses peak detection with harmonic verification.

        Args:
            spectrum: Amplitude spectrum
            freqs: Corresponding frequencies
            freq_range: Range to search for fundamental

        Returns:
            Detected fundamental frequency in Hz
        """
        # Mask for frequency range
        mask = (freqs >= freq_range[0]) & (freqs <= freq_range[1])
        if not np.any(mask):
            return freq_range[0]  # Default

        freq_masked = freqs[mask]
        spec_masked = spectrum[mask]

        # Find prominent peaks
        from scipy.signal import find_peaks
        peaks, properties = find_peaks(
            spec_masked,
            height=np.median(spec_masked) * 2,
            distance=5
        )

        if len(peaks) == 0:
            # Fallback to simple max
            return freq_masked[np.argmax(spec_masked)]

        # Score each peak by harmonic support
        # Prefer lower frequencies as fundamentals (they generate harmonics)
        candidates = []

        for peak_idx in peaks:
            f0_candidate = freq_masked[peak_idx]
            harmonic_score = self._harmonic_score(spectrum, freqs, f0_candidate)
            peak_amplitude = spec_masked[peak_idx]
            candidates.append((f0_candidate, harmonic_score, peak_amplitude))

        if not candidates:
            return freq_masked[np.argmax(spec_masked)]

        # Sort by: 1) harmonic score (desc), 2) frequency (asc - prefer lower)
        # Among peaks with similar harmonic scores, prefer lower frequency
        max_score = max(c[1] for c in candidates)
        threshold = max_score * 0.7  # Within 70% of best score

        good_candidates = [c for c in candidates if c[1] >= threshold]

        # Among good candidates, pick the lowest frequency (fundamental)
        good_candidates.sort(key=lambda x: x[0])  # Sort by frequency ascending

        return good_candidates[0][0]

    def _harmonic_score(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        f0: float,
        n_harmonics: int = 8
    ) -> float:
        """Score a candidate F0 by harmonic energy."""
        score = 0.0
        for h in range(1, n_harmonics + 1):
            f_harm = f0 * h
            idx = np.argmin(np.abs(freqs - f_harm))
            if idx < len(spectrum):
                score += spectrum[idx] / h  # Weight by 1/harmonic number
        return score

    def find_modes(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        f0: float = None,
        noise_floor: np.ndarray = None
    ) -> List[HexagonalMode]:
        """
        Find hexagonal frequency modes in a spectrum.

        Args:
            spectrum: Amplitude spectrum
            freqs: Corresponding frequencies
            f0: Fundamental frequency (auto-detected if None)
            noise_floor: Noise spectrum for SNR calculation

        Returns:
            List of detected HexagonalMode objects
        """
        if f0 is None:
            f0 = self.detect_fundamental(spectrum, freqs)

        if noise_floor is None:
            noise_floor = np.percentile(spectrum, 25)
            if isinstance(noise_floor, (int, float)):
                noise_floor = np.full_like(spectrum, noise_floor)

        modes = []

        for i, expected_ratio in enumerate(self.ratios):
            expected_freq = f0 * expected_ratio

            # Search window around expected frequency
            window = expected_freq * self.tolerance
            mask = (freqs >= expected_freq - window) & (freqs <= expected_freq + window)

            if not np.any(mask):
                continue

            # Find peak in window
            freq_window = freqs[mask]
            spec_window = spectrum[mask]
            noise_window = noise_floor[mask] if len(noise_floor) > 1 else noise_floor

            peak_idx = np.argmax(spec_window)
            peak_freq = freq_window[peak_idx]
            peak_amp = spec_window[peak_idx]

            # Calculate SNR
            local_noise = np.median(noise_window) if hasattr(noise_window, '__len__') else noise_window
            snr = peak_amp / max(local_noise, 1e-10)

            if snr < self.min_snr:
                continue

            # Calculate actual ratio and deviation
            actual_ratio = peak_freq / f0
            deviation = abs(actual_ratio - expected_ratio) / expected_ratio

            # Confidence based on SNR and deviation
            confidence = min(1.0, (snr / self.min_snr) * (1 - deviation / self.tolerance))

            mode = HexagonalMode(
                frequency=peak_freq,
                ratio_index=i,
                ratio_value=actual_ratio,
                expected_ratio=expected_ratio,
                deviation=deviation,
                amplitude=peak_amp,
                snr=snr,
                confidence=confidence
            )
            modes.append(mode)

        return modes

    def transform(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        f0: float = None
    ) -> np.ndarray:
        """
        Project spectrum onto hexagonal ratio space.

        Returns a vector of amplitudes at each hexagonal ratio.

        Args:
            spectrum: Amplitude spectrum
            freqs: Corresponding frequencies
            f0: Fundamental frequency

        Returns:
            Array of shape (len(ratios),) with amplitudes
        """
        if f0 is None:
            f0 = self.detect_fundamental(spectrum, freqs)

        result = np.zeros(len(self.ratios))

        for i, ratio in enumerate(self.ratios):
            expected_freq = f0 * ratio

            # Gaussian weighting around expected frequency
            sigma = expected_freq * self.tolerance
            weights = np.exp(-0.5 * ((freqs - expected_freq) / sigma) ** 2)

            # Weighted sum
            result[i] = np.sum(spectrum * weights) / max(np.sum(weights), 1e-10)

        return result

    def inverse_transform(
        self,
        hex_amplitudes: np.ndarray,
        freqs: np.ndarray,
        f0: float
    ) -> np.ndarray:
        """
        Reconstruct spectrum from hexagonal amplitudes.

        This creates a "hexagonally filtered" version of the spectrum.

        Args:
            hex_amplitudes: Amplitudes at each hexagonal ratio
            freqs: Target frequency array
            f0: Fundamental frequency

        Returns:
            Reconstructed spectrum
        """
        spectrum = np.zeros_like(freqs)

        for i, (ratio, amp) in enumerate(zip(self.ratios, hex_amplitudes)):
            expected_freq = f0 * ratio

            # Gaussian peak at each hexagonal frequency
            sigma = expected_freq * self.tolerance
            peak = amp * np.exp(-0.5 * ((freqs - expected_freq) / sigma) ** 2)
            spectrum += peak

        return spectrum


class FrequencyRatioAnalyzer:
    """
    Analyze frequency ratios in spectral data.

    Specialized for detecting and validating hexagonal patterns
    in gravitational wave and audio signals.
    """

    def __init__(self, tolerance: float = DEFAULT_TOLERANCE):
        self.tolerance = tolerance
        self.hexagonal_ratios = HEXAGONAL_RATIOS

    def extract_peaks(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        n_peaks: int = 10,
        min_separation: float = 10.0
    ) -> List[Tuple[float, float]]:
        """
        Extract prominent frequency peaks from spectrum.

        Args:
            spectrum: Amplitude spectrum
            freqs: Frequency array
            n_peaks: Maximum number of peaks to extract
            min_separation: Minimum Hz between peaks

        Returns:
            List of (frequency, amplitude) tuples
        """
        from scipy.signal import find_peaks

        # Find peaks
        freq_resolution = freqs[1] - freqs[0] if len(freqs) > 1 else 1.0
        distance = int(min_separation / freq_resolution)

        peak_indices, properties = find_peaks(
            spectrum,
            height=np.median(spectrum) * 2,
            distance=max(1, distance)
        )

        if len(peak_indices) == 0:
            return []

        # Sort by amplitude and take top n
        amplitudes = spectrum[peak_indices]
        sorted_idx = np.argsort(amplitudes)[::-1][:n_peaks]

        peaks = [
            (freqs[peak_indices[i]], spectrum[peak_indices[i]])
            for i in sorted_idx
        ]

        return peaks

    def compute_all_ratios(
        self,
        peaks: List[Tuple[float, float]]
    ) -> List[Tuple[float, float, float, float]]:
        """
        Compute all pairwise frequency ratios.

        Args:
            peaks: List of (frequency, amplitude) tuples

        Returns:
            List of (freq1, freq2, ratio, combined_amplitude) tuples
        """
        ratios = []

        for i, (f1, a1) in enumerate(peaks):
            for j, (f2, a2) in enumerate(peaks):
                if j <= i:
                    continue

                # Always ratio > 1
                if f1 > f2:
                    ratio = f1 / f2
                else:
                    ratio = f2 / f1

                combined_amp = np.sqrt(a1 * a2)
                ratios.append((f1, f2, ratio, combined_amp))

        return ratios

    def classify_ratio(
        self,
        ratio: float
    ) -> Tuple[int, float, bool]:
        """
        Classify a frequency ratio as hexagonal or not.

        Args:
            ratio: Frequency ratio to classify

        Returns:
            Tuple of (closest_ratio_index, deviation, is_hexagonal)
        """
        # Find closest hexagonal ratio
        deviations = [
            abs(ratio - hex_ratio) / hex_ratio
            for hex_ratio in self.hexagonal_ratios
        ]

        closest_idx = np.argmin(deviations)
        min_deviation = deviations[closest_idx]
        is_hexagonal = min_deviation <= self.tolerance

        return closest_idx, min_deviation, is_hexagonal

    def analyze(
        self,
        spectrum: np.ndarray,
        freqs: np.ndarray,
        n_peaks: int = 10
    ) -> Dict:
        """
        Complete hexagonal ratio analysis of a spectrum.

        Args:
            spectrum: Amplitude spectrum
            freqs: Frequency array
            n_peaks: Number of peaks to analyze

        Returns:
            Dictionary with analysis results
        """
        # Extract peaks
        peaks = self.extract_peaks(spectrum, freqs, n_peaks)

        if len(peaks) < 2:
            return {
                "peaks": peaks,
                "ratios": [],
                "hexagonal_count": 0,
                "total_ratios": 0,
                "hexagonal_fraction": 0.0,
                "is_hexagonal": False
            }

        # Compute ratios
        all_ratios = self.compute_all_ratios(peaks)

        # Classify each ratio
        hexagonal_count = 0
        classified_ratios = []

        for f1, f2, ratio, amp in all_ratios:
            idx, dev, is_hex = self.classify_ratio(ratio)
            classified_ratios.append({
                "freq1": f1,
                "freq2": f2,
                "ratio": ratio,
                "amplitude": amp,
                "closest_hex_index": idx,
                "closest_hex_ratio": self.hexagonal_ratios[idx],
                "deviation": dev,
                "is_hexagonal": is_hex
            })
            if is_hex:
                hexagonal_count += 1

        total = len(all_ratios)
        hex_fraction = hexagonal_count / total if total > 0 else 0.0

        # Is the spectrum hexagonal? (majority of ratios match)
        is_hexagonal = hex_fraction > 0.5

        return {
            "peaks": peaks,
            "ratios": classified_ratios,
            "hexagonal_count": hexagonal_count,
            "total_ratios": total,
            "hexagonal_fraction": hex_fraction,
            "is_hexagonal": is_hexagonal
        }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def compute_ratio_proximity(
    ratio: float,
    target_ratios: List[float] = None,
    tolerance: float = DEFAULT_TOLERANCE
) -> Tuple[int, float, bool]:
    """
    Check if a ratio is close to any hexagonal ratio.

    Args:
        ratio: Ratio to check
        target_ratios: Target ratios (default: HEXAGONAL_RATIOS)
        tolerance: Matching tolerance

    Returns:
        (closest_index, deviation, is_match)
    """
    if target_ratios is None:
        target_ratios = HEXAGONAL_RATIOS

    deviations = [abs(ratio - t) / t for t in target_ratios]
    closest = np.argmin(deviations)

    return closest, deviations[closest], deviations[closest] <= tolerance


def find_hexagonal_modes(
    spectrum: np.ndarray,
    freqs: np.ndarray,
    f0: float = None,
    tolerance: float = DEFAULT_TOLERANCE,
    min_snr: float = 3.0
) -> List[HexagonalMode]:
    """
    Convenience function to find hexagonal modes in a spectrum.

    Args:
        spectrum: Amplitude spectrum
        freqs: Frequency array
        f0: Fundamental frequency (auto-detected if None)
        tolerance: Matching tolerance
        min_snr: Minimum SNR threshold

    Returns:
        List of detected HexagonalMode objects
    """
    transform = HexagonalTransform(tolerance=tolerance, min_snr=min_snr)
    return transform.find_modes(spectrum, freqs, f0)
