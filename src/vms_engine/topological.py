"""
VMS Topological Engine
======================

Audio cleaner using topological signal representation.
Treats signal as a 3D solid form that can be sculpted.

Key insight: In 3D (time × frequency × amplitude), coherent signal
forms connected "mountains" while noise is scattered "dust".
We can separate them topologically rather than by thresholds.

Methods combined:
1. Spectral Subtraction (Boll)
2. Wiener Filter (MMSE)
3. Harmonic Enhancement (VMS)
4. Topological Filtering (connected components)

The output blends all methods, weighted by local signal confidence.

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
from typing import Optional, Tuple, Dict, List
from dataclasses import dataclass
from scipy import ndimage
from scipy import signal as scipy_signal


@dataclass
class TopologicalResult:
    """Result from topological cleaning."""
    audio: np.ndarray
    sample_rate: int
    components_found: int          # Number of signal "mountains"
    noise_regions_removed: int     # Number of noise regions filtered
    method_weights: Dict[str, float]  # Final blend weights
    snr_improvement_db: float


class VMSTopologicalCleaner:
    """
    Topological audio cleaner.

    Represents the spectrogram as a 3D landscape where:
    - X = time
    - Y = frequency
    - Z = amplitude (dB)

    Signal forms connected mountain ranges.
    Noise forms isolated peaks or flat plains.

    We identify and preserve the connected mountains (signal)
    while removing isolated noise.
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512
    ):
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size

        # Windows
        self.window = np.hanning(frame_size)
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Hexagonal ratios for harmonic detection
        self.hex_ratios = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def _compute_spectrogram(self, audio: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute spectrogram with phase."""
        n_frames = (len(audio) - self.frame_size) // self.hop_size + 1
        n_freqs = self.frame_size // 2 + 1

        magnitude = np.zeros((n_freqs, n_frames))
        phase = np.zeros((n_freqs, n_frames), dtype=complex)

        for i in range(n_frames):
            start = i * self.hop_size
            frame = audio[start:start + self.frame_size]
            if len(frame) < self.frame_size:
                frame = np.pad(frame, (0, self.frame_size - len(frame)))

            windowed = frame * self.window
            spec = np.fft.rfft(windowed)
            magnitude[:, i] = np.abs(spec)
            phase[:, i] = spec / (np.abs(spec) + 1e-10)

        times = np.arange(n_frames) * self.hop_size / self.sample_rate
        return magnitude, phase, times

    def _estimate_noise_spectrum(self, magnitude: np.ndarray, noise_frames: int = None) -> np.ndarray:
        """Estimate noise from initial frames."""
        if noise_frames is None:
            noise_frames = int(0.5 * self.sample_rate / self.hop_size)

        noise_frames = min(noise_frames, magnitude.shape[1])
        return np.median(magnitude[:, :noise_frames], axis=1)

    # =========================================================================
    # CLEANING METHODS (all combined)
    # =========================================================================

    def _spectral_subtraction(
        self,
        magnitude: np.ndarray,
        noise_spectrum: np.ndarray,
        alpha: float = 1.0
    ) -> np.ndarray:
        """Classic Boll spectral subtraction."""
        noise_2d = noise_spectrum[:, np.newaxis]
        cleaned = np.maximum(magnitude - alpha * noise_2d, 0.1 * magnitude)
        return cleaned

    def _wiener_filter(
        self,
        magnitude: np.ndarray,
        noise_spectrum: np.ndarray
    ) -> np.ndarray:
        """MMSE Wiener filter."""
        noise_power = noise_spectrum[:, np.newaxis] ** 2
        signal_power = magnitude ** 2

        gain = np.maximum(signal_power - noise_power, 0) / (signal_power + 1e-10)
        gain = np.sqrt(gain)
        gain = np.maximum(gain, 0.1)

        return magnitude * gain

    def _harmonic_enhancement(
        self,
        magnitude: np.ndarray,
        f0: float = 150.0
    ) -> np.ndarray:
        """Enhance harmonic structure using hexagonal ratios."""
        enhanced = magnitude.copy()
        freq_res = self.sample_rate / self.frame_size

        # Create harmonic mask
        mask = np.zeros(magnitude.shape[0])
        all_ratios = list(range(1, 9)) + self.hex_ratios

        for ratio in all_ratios:
            f_harm = f0 * ratio
            if f_harm > self.sample_rate / 2:
                break

            idx = int(f_harm / freq_res)
            if idx < len(mask):
                width = max(2, int(f_harm / 100))
                for i in range(max(0, idx - width), min(len(mask), idx + width + 1)):
                    dist = abs(i - idx)
                    mask[i] = max(mask[i], np.exp(-dist**2 / (width**2)))

        # Apply mask as boost
        mask_2d = mask[:, np.newaxis]
        enhanced = magnitude * (1 + 0.5 * mask_2d)

        return enhanced

    def _topological_filter(
        self,
        magnitude: np.ndarray,
        threshold_percentile: float = 30,
        min_component_size: int = 50
    ) -> Tuple[np.ndarray, int, int]:
        """
        Filter based on connected component topology.

        Signal forms large connected regions (mountains).
        Noise forms small isolated regions (dust).
        """
        # Convert to dB for thresholding
        mag_db = 20 * np.log10(magnitude + 1e-10)

        # Threshold to binary
        threshold = np.percentile(mag_db, threshold_percentile)
        binary = mag_db > threshold

        # Find connected components
        labeled, n_components = ndimage.label(binary)

        # Measure component sizes
        component_sizes = ndimage.sum(binary, labeled, range(1, n_components + 1))

        # Keep only large components (signal mountains)
        large_components = []
        small_removed = 0

        for i, size in enumerate(component_sizes):
            if size >= min_component_size:
                large_components.append(i + 1)
            else:
                small_removed += 1

        # Create mask keeping only large components
        mask = np.zeros_like(labeled, dtype=float)
        for comp_id in large_components:
            mask[labeled == comp_id] = 1.0

        # Smooth mask edges
        mask = ndimage.gaussian_filter(mask.astype(float), sigma=1)

        # Apply mask with soft floor
        filtered = magnitude * (0.1 + 0.9 * mask)

        return filtered, len(large_components), small_removed

    def _detect_f0(self, magnitude: np.ndarray) -> float:
        """Detect fundamental frequency."""
        avg_spectrum = np.mean(magnitude, axis=1)
        freq_res = self.sample_rate / self.frame_size

        # Search voice range
        idx_min = int(60 / freq_res)
        idx_max = min(int(400 / freq_res), len(avg_spectrum) - 1)

        if idx_max <= idx_min:
            return 150.0

        voice_spec = avg_spectrum[idx_min:idx_max + 1]
        peak_idx = np.argmax(voice_spec)

        return (idx_min + peak_idx) * freq_res

    def _blend_methods(
        self,
        results: Dict[str, np.ndarray],
        magnitude: np.ndarray,
        noise_spectrum: np.ndarray
    ) -> Tuple[np.ndarray, Dict[str, float]]:
        """
        Blend all cleaning methods based on local confidence.

        High SNR regions: prefer aggressive methods (spectral sub)
        Low SNR regions: prefer conservative methods (Wiener)
        Harmonic regions: boost harmonic enhancement weight
        """
        n_freqs, n_frames = magnitude.shape

        # Calculate local SNR
        noise_2d = noise_spectrum[:, np.newaxis]
        local_snr = magnitude / (noise_2d + 1e-10)
        local_snr_normalized = np.clip(local_snr / 3, 0, 1)  # Normalize to 0-1

        # Weights that vary spatially
        # High SNR → more aggressive
        # Low SNR → more conservative
        w_spectral = local_snr_normalized * 0.4  # 0 to 0.4
        w_wiener = (1 - local_snr_normalized) * 0.4  # 0.4 to 0
        w_harmonic = 0.1  # Constant small boost
        w_topo = 0.1  # Constant small contribution

        # Normalize weights
        w_total = w_spectral + w_wiener + w_harmonic + w_topo
        w_spectral /= w_total
        w_wiener /= w_total
        w_harmonic /= (w_total / 0.1)
        w_topo /= (w_total / 0.1)

        # Blend
        blended = (
            results['spectral'] * w_spectral +
            results['wiener'] * w_wiener +
            results['harmonic'] * 0.1 +
            results['topological'] * 0.1
        )

        # Compute average weights for reporting
        avg_weights = {
            'spectral_sub': float(np.mean(w_spectral)),
            'wiener': float(np.mean(w_wiener)),
            'harmonic': 0.1,
            'topological': 0.1
        }

        return blended, avg_weights

    def _reconstruct_audio(
        self,
        magnitude: np.ndarray,
        phase: np.ndarray,
        length: int
    ) -> np.ndarray:
        """Reconstruct audio from magnitude and phase."""
        n_frames = magnitude.shape[1]

        output = np.zeros(length + self.frame_size)
        window_sum = np.zeros(length + self.frame_size)

        for i in range(n_frames):
            start = i * self.hop_size

            # Reconstruct complex spectrum
            spec = magnitude[:, i] * phase[:, i]

            # IFFT
            frame = np.fft.irfft(spec, n=self.frame_size)

            # Overlap-add with synthesis window
            output[start:start + self.frame_size] += frame * self.window
            window_sum[start:start + self.frame_size] += self.window ** 2

        # Normalize
        window_sum = np.maximum(window_sum, 1e-8)
        output = output / window_sum

        return output[:length]

    def clean_audio(
        self,
        audio: np.ndarray,
        noise_sample: np.ndarray = None
    ) -> TopologicalResult:
        """
        Clean audio using topological multi-method approach.

        Args:
            audio: Input audio (mono)
            noise_sample: Optional noise-only sample

        Returns:
            TopologicalResult with cleaned audio and metrics
        """
        audio = np.asarray(audio, dtype=np.float32)
        original_length = len(audio)

        # Compute spectrogram
        magnitude, phase, times = self._compute_spectrogram(audio)

        # Estimate noise
        if noise_sample is not None:
            noise_mag, _, _ = self._compute_spectrogram(noise_sample)
            noise_spectrum = np.median(noise_mag, axis=1)
        else:
            noise_spectrum = self._estimate_noise_spectrum(magnitude)

        # Detect F0
        f0 = self._detect_f0(magnitude)

        # Apply all methods
        results = {}
        results['spectral'] = self._spectral_subtraction(magnitude, noise_spectrum)
        results['wiener'] = self._wiener_filter(magnitude, noise_spectrum)
        results['harmonic'] = self._harmonic_enhancement(magnitude, f0)
        results['topological'], n_components, n_removed = self._topological_filter(magnitude)

        # Blend methods
        blended, weights = self._blend_methods(results, magnitude, noise_spectrum)

        # Ensure we don't exceed original (prevent amplification artifacts)
        blended = np.minimum(blended, magnitude * 1.2)

        # Reconstruct audio
        cleaned = self._reconstruct_audio(blended, phase, original_length)

        # Prevent clipping
        if np.max(np.abs(cleaned)) > 0.95:
            cleaned = cleaned / np.max(np.abs(cleaned)) * 0.95

        # Calculate improvement
        # Use middle segment for measurement
        margin = min(10000, len(audio) // 10)
        if len(audio) > 2 * margin:
            input_seg = audio[margin:-margin]
            output_seg = cleaned[margin:-margin]

            input_power = np.mean(input_seg ** 2)
            output_power = np.mean(output_seg ** 2)

            # Rough SNR improvement estimate
            if output_power > 0:
                improvement = 10 * np.log10(output_power / (input_power + 1e-10))
            else:
                improvement = 0.0
        else:
            improvement = 0.0

        return TopologicalResult(
            audio=cleaned,
            sample_rate=self.sample_rate,
            components_found=n_components,
            noise_regions_removed=n_removed,
            method_weights=weights,
            snr_improvement_db=improvement
        )

    def visualize_topology(
        self,
        audio: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Generate data for 3D topology visualization.

        Returns data suitable for Plotly 3D surface plots.
        """
        magnitude, phase, times = self._compute_spectrogram(audio)
        mag_db = 20 * np.log10(magnitude + 1e-10)

        # Connected components
        threshold = np.percentile(mag_db, 30)
        binary = mag_db > threshold
        labeled, n_components = ndimage.label(binary)

        return {
            'magnitude_db': mag_db,
            'frequencies': self.freqs,
            'times': times,
            'components': labeled,
            'n_components': n_components,
            'threshold_db': threshold
        }
