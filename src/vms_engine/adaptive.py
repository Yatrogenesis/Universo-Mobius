"""
VMS Adaptive Engine
===================

Audio cleaner with fully interdependent parameters.
All settings derive from signal characteristics - no manual tuning needed.

The key insight: parameters should not be independent constants,
but functions of each other and of the signal itself.

Core dependencies:
- noise_estimate → derives gain_floor
- gain_floor → affects harmonic_weight
- harmonic_weight → modulates spectral_smoothing
- SNR_estimate → adjusts all of the above

This creates a self-regulating system that adapts to any input.

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
from typing import Optional, Tuple
from dataclasses import dataclass


@dataclass
class AdaptiveState:
    """Dynamic state of the adaptive cleaner."""
    snr_estimate: float = 15.0          # Estimated input SNR (dB) - start optimistic
    noise_power: float = 0.01           # Current noise power estimate
    signal_power: float = 0.1           # Current signal power estimate
    fundamental_freq: float = 150.0     # Detected F0 (Hz)
    harmonic_strength: float = 0.3      # How harmonic is the signal (0-1)
    stationarity: float = 0.7           # How stationary is the noise (0-1) - assume stationary

    # Derived parameters (all depend on above)
    @property
    def gain_floor(self) -> float:
        """
        Minimum gain - depends on SNR AND stationarity.
        Low SNR + stationary noise → can still be aggressive
        Low SNR + non-stationary → be conservative
        """
        # Base floor from SNR (0.02 to 0.15)
        snr_clamped = max(0, min(30, self.snr_estimate))
        base_floor = 0.02 + 0.13 * (1 - snr_clamped / 30)

        # Adjust by stationarity - stationary noise = lower floor ok
        stationarity_factor = 0.5 + 0.5 * self.stationarity

        return base_floor * stationarity_factor

    @property
    def harmonic_weight(self) -> float:
        """
        Harmonic preservation weight - depends on harmonic_strength and SNR.
        More harmonic signal → preserve more harmonics.
        Lower SNR → be more careful with harmonics.
        """
        base_weight = self.harmonic_strength * 0.8
        snr_factor = min(1.0, self.snr_estimate / 15)
        return base_weight * (0.5 + 0.5 * snr_factor)

    @property
    def noise_subtraction_alpha(self) -> float:
        """
        Noise subtraction aggressiveness - depends on stationarity AND SNR.
        Stationary noise + good SNR → can subtract aggressively.
        Non-stationary or low SNR → be conservative.
        """
        # Base from stationarity (0.8 to 1.5)
        base_alpha = 0.8 + 0.7 * self.stationarity

        # Boost when SNR is good (confident separation)
        snr_factor = min(1.5, 1 + self.snr_estimate / 20)

        return base_alpha * snr_factor

    @property
    def spectral_smoothing(self) -> float:
        """
        Smoothing factor - depends on SNR and stationarity.
        Low SNR + non-stationary → more smoothing.
        High SNR + stationary → less smoothing.
        """
        return 0.2 + 0.4 * (1 - self.stationarity) * (1 - min(1, self.snr_estimate / 20))


class VMSAdaptiveCleaner:
    """
    Fully adaptive VMS audio cleaner.

    All parameters are interdependent and derive from signal analysis.
    The system self-regulates based on:
    1. Estimated SNR
    2. Signal harmonicity
    3. Noise stationarity
    4. Spectral characteristics

    No manual tuning required - just feed it audio.
    """

    # Hexagonal ratios for enhanced harmonic detection
    HEXAGONAL_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512
    ):
        """
        Initialize adaptive cleaner.

        Only signal-agnostic parameters are set here.
        Everything else adapts to the input.
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size

        # Frequency bins
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Windows (WOLA)
        self.analysis_window = np.sqrt(np.hanning(frame_size))
        self.synthesis_window = np.sqrt(np.hanning(frame_size))
        self.wola_window = np.hanning(frame_size)

        # Adaptive state
        self.state = AdaptiveState()

        # Running estimates (exponential moving averages)
        self._noise_spectrum_ema: Optional[np.ndarray] = None
        self._signal_spectrum_ema: Optional[np.ndarray] = None
        self._ema_alpha = 0.1  # EMA smoothing factor

    def _update_snr_estimate(self, spectrum: np.ndarray, noise_spectrum: np.ndarray):
        """Update SNR estimate from current frame."""
        signal_power = np.mean(spectrum ** 2)
        noise_power = np.mean(noise_spectrum ** 2)

        # Calculate instantaneous SNR
        if noise_power > 1e-10:
            instant_snr = 10 * np.log10(max(signal_power / noise_power, 0.01))
            instant_snr = max(-10, min(40, instant_snr))
        else:
            instant_snr = 20.0

        # Voice activity: check if signal is significantly above noise
        # AND has harmonic structure
        voice_active = (
            signal_power > noise_power * 1.5 and
            self.state.harmonic_strength > 0.4
        )

        if voice_active:
            # Voice detected - update SNR more aggressively toward actual value
            alpha = 0.2
            target_snr = instant_snr
        else:
            # No voice - slowly decay SNR estimate
            alpha = 0.02
            target_snr = max(5.0, instant_snr)  # Don't go below 5 dB estimate

        self.state.snr_estimate = (1 - alpha) * self.state.snr_estimate + alpha * target_snr
        self.state.signal_power = signal_power
        self.state.noise_power = noise_power

    def _update_harmonic_strength(self, spectrum: np.ndarray, f0: float):
        """
        Measure how harmonic the signal is.
        Harmonic signals have energy concentrated at integer multiples of F0.
        """
        if f0 < 50:
            self.state.harmonic_strength = 0.3
            return

        freq_res = self.sample_rate / self.frame_size
        total_energy = np.sum(spectrum ** 2)

        if total_energy < 1e-10:
            self.state.harmonic_strength = 0.3
            return

        # Energy at harmonics
        harmonic_energy = 0
        for h in range(1, 10):
            f_harm = f0 * h
            if f_harm > self.sample_rate / 2:
                break

            idx = int(f_harm / freq_res)
            if idx < len(spectrum):
                # Sum energy in small window around harmonic
                start = max(0, idx - 2)
                end = min(len(spectrum), idx + 3)
                harmonic_energy += np.sum(spectrum[start:end] ** 2)

        ratio = harmonic_energy / total_energy
        # Smooth update
        self.state.harmonic_strength = 0.8 * self.state.harmonic_strength + 0.2 * min(1.0, ratio)

    def _update_stationarity(self, spectrum: np.ndarray):
        """
        Measure noise stationarity by comparing to running average.
        Stationary noise has consistent spectral shape.
        """
        if self._noise_spectrum_ema is None:
            self._noise_spectrum_ema = spectrum.copy()
            self.state.stationarity = 0.5
            return

        # Compare current to EMA
        diff = np.abs(spectrum - self._noise_spectrum_ema)
        relative_diff = np.mean(diff) / (np.mean(self._noise_spectrum_ema) + 1e-10)

        # Low difference = stationary
        stationarity = 1.0 / (1.0 + relative_diff * 10)

        # Update EMA
        self._noise_spectrum_ema = (
            (1 - self._ema_alpha) * self._noise_spectrum_ema +
            self._ema_alpha * spectrum
        )

        # Smooth update
        self.state.stationarity = 0.9 * self.state.stationarity + 0.1 * stationarity

    def _detect_fundamental(self, spectrum: np.ndarray) -> float:
        """Detect F0 with harmonic verification."""
        freq_res = self.sample_rate / self.frame_size

        # Search in voice range
        f_min, f_max = 60, 400
        idx_min = int(f_min / freq_res)
        idx_max = min(int(f_max / freq_res), len(spectrum) - 1)

        if idx_max <= idx_min:
            return self.state.fundamental_freq

        # Find peaks
        voice_spec = spectrum[idx_min:idx_max + 1]
        if len(voice_spec) < 3:
            return self.state.fundamental_freq

        # Score candidates by harmonic support
        best_f0 = self.state.fundamental_freq
        best_score = 0

        for i in range(len(voice_spec)):
            if voice_spec[i] < np.median(voice_spec) * 1.5:
                continue

            f0_candidate = (idx_min + i) * freq_res
            score = self._harmonic_score(spectrum, f0_candidate)

            if score > best_score:
                best_score = score
                best_f0 = f0_candidate

        # Smooth F0 tracking
        self.state.fundamental_freq = 0.7 * self.state.fundamental_freq + 0.3 * best_f0
        return self.state.fundamental_freq

    def _harmonic_score(self, spectrum: np.ndarray, f0: float) -> float:
        """Score F0 candidate by harmonic energy."""
        freq_res = self.sample_rate / self.frame_size
        score = 0

        for h in range(1, 8):
            f_harm = f0 * h
            if f_harm > self.sample_rate / 2:
                break

            idx = int(f_harm / freq_res)
            if idx < len(spectrum):
                score += spectrum[idx] / h  # Weight by 1/harmonic

        return score

    def _create_adaptive_harmonic_mask(self, spectrum: np.ndarray, f0: float) -> np.ndarray:
        """
        Create harmonic mask with adaptive width.
        Width depends on SNR - lower SNR needs wider protection.
        """
        mask = np.zeros(len(spectrum))
        freq_res = self.sample_rate / self.frame_size

        # Width adapts to SNR
        base_width = max(2, int(10 / (1 + self.state.snr_estimate / 10)))

        # Combine natural + hexagonal harmonics
        all_ratios = sorted(set(
            list(range(1, 9)) +  # Natural: 1, 2, 3, ..., 8
            [r for r in self.HEXAGONAL_RATIOS if r <= 8]  # Hexagonal
        ))

        for ratio in all_ratios:
            f_harm = f0 * ratio
            if f_harm > self.sample_rate / 2:
                break

            idx = int(f_harm / freq_res)
            if idx >= len(mask):
                break

            # Width proportional to frequency
            width = base_width + int(f_harm / 200)

            for i in range(max(0, idx - width), min(len(mask), idx + width + 1)):
                dist = abs(i - idx)
                gauss = np.exp(-dist ** 2 / (2 * (width / 2) ** 2))
                mask[i] = max(mask[i], gauss)

        return mask

    def _adaptive_spectral_subtraction(
        self,
        spectrum: np.ndarray,
        noise_spectrum: np.ndarray
    ) -> np.ndarray:
        """
        Spectral subtraction with all parameters derived from state.
        """
        # Get adaptive parameters from state
        alpha = self.state.noise_subtraction_alpha
        gain_floor = self.state.gain_floor
        harmonic_weight = self.state.harmonic_weight

        # Wiener filtering
        spectrum_safe = np.maximum(spectrum, 1e-10)
        noise_safe = np.maximum(noise_spectrum, 1e-10)

        signal_power = spectrum_safe ** 2
        noise_power = noise_safe ** 2

        # Adaptive gain
        gain_sq = np.maximum(signal_power - alpha * noise_power, 0) / signal_power
        gain = np.sqrt(gain_sq)
        gain = np.maximum(gain, gain_floor)

        # Apply base gain
        cleaned = spectrum * gain

        # Harmonic preservation (only if signal is harmonic)
        if self.state.harmonic_strength > 0.3:
            f0 = self.state.fundamental_freq
            harmonic_mask = self._create_adaptive_harmonic_mask(spectrum, f0)

            # Boost harmonics that were over-suppressed
            harmonic_boost = 1 + (1 - gain) * harmonic_mask * harmonic_weight
            cleaned = cleaned * harmonic_boost

            # Don't exceed original
            cleaned = np.minimum(cleaned, spectrum)

        return cleaned

    def process_frame(self, frame: np.ndarray, noise_frame: np.ndarray = None) -> np.ndarray:
        """
        Process single frame with full adaptation.

        All parameters update based on this frame's characteristics.
        """
        if len(frame) < self.frame_size:
            frame = np.pad(frame, (0, self.frame_size - len(frame)))

        # Analysis
        windowed = frame * self.analysis_window
        fft_data = np.fft.rfft(windowed)
        spectrum = np.abs(fft_data)
        phase = np.angle(fft_data)

        # Get noise estimate
        if noise_frame is not None:
            noise_windowed = noise_frame * self.analysis_window
            noise_spectrum = np.abs(np.fft.rfft(noise_windowed))
        elif self._noise_spectrum_ema is not None:
            noise_spectrum = self._noise_spectrum_ema
        else:
            # Bootstrap: use lower percentile of current frame
            noise_spectrum = np.percentile(spectrum, 20) * np.ones_like(spectrum)

        # Update all adaptive state (ORDER MATTERS!)
        # 1. First detect fundamental and harmonic strength
        f0 = self._detect_fundamental(spectrum)
        self._update_harmonic_strength(spectrum, f0)
        # 2. Then update SNR (uses harmonic_strength for voice detection)
        self._update_snr_estimate(spectrum, noise_spectrum)
        # 3. Finally stationarity
        self._update_stationarity(spectrum)

        # Clean with adaptive parameters
        cleaned_spectrum = self._adaptive_spectral_subtraction(spectrum, noise_spectrum)

        # Reconstruct
        reconstructed = cleaned_spectrum * np.exp(1j * phase)
        cleaned_frame = np.fft.irfft(reconstructed, n=self.frame_size)

        return cleaned_frame * self.synthesis_window

    def clean_audio(
        self,
        audio: np.ndarray,
        noise_sample: np.ndarray = None
    ) -> Tuple[np.ndarray, AdaptiveState]:
        """
        Clean complete audio with full adaptation.

        Returns cleaned audio and final adaptive state.
        """
        audio = np.asarray(audio, dtype=np.float32)

        # Initialize noise estimate from first segment
        if noise_sample is not None:
            init_segment = noise_sample
        else:
            init_segment = audio[:int(0.5 * self.sample_rate)]

        # Bootstrap noise spectrum
        n_init_frames = max(1, (len(init_segment) - self.frame_size) // self.hop_size)
        noise_spectra = []

        for i in range(n_init_frames):
            start = i * self.hop_size
            frame = init_segment[start:start + self.frame_size]
            if len(frame) < self.frame_size:
                frame = np.pad(frame, (0, self.frame_size - len(frame)))

            windowed = frame * self.analysis_window
            spectrum = np.abs(np.fft.rfft(windowed))
            noise_spectra.append(spectrum)

        if noise_spectra:
            self._noise_spectrum_ema = np.median(noise_spectra, axis=0)

        # Process all frames
        output = np.zeros(len(audio) + self.frame_size)
        window_sum = np.zeros(len(audio) + self.frame_size)

        n_frames = (len(audio) - self.frame_size) // self.hop_size + 1

        for i in range(n_frames):
            start = i * self.hop_size
            frame = audio[start:start + self.frame_size]

            cleaned_frame = self.process_frame(frame)

            output[start:start + self.frame_size] += cleaned_frame
            window_sum[start:start + self.frame_size] += self.wola_window

        # Normalize
        min_sum = np.max(window_sum) * 0.5
        valid = window_sum >= min_sum

        result = np.zeros_like(output)
        result[valid] = output[valid] / window_sum[valid]

        # Handle edges
        edge = ~valid & (window_sum > 1e-8)
        if np.any(edge):
            result[edge] = output[edge] / np.max(window_sum)

        result = result[:len(audio)]

        # Prevent clipping
        if np.max(np.abs(result)) > 0.95:
            result = result / np.max(np.abs(result)) * 0.95

        return result, self.state

    def get_state_summary(self) -> dict:
        """Get current adaptive state as dictionary."""
        return {
            'snr_estimate_db': self.state.snr_estimate,
            'fundamental_freq_hz': self.state.fundamental_freq,
            'harmonic_strength': self.state.harmonic_strength,
            'stationarity': self.state.stationarity,
            'derived': {
                'gain_floor': self.state.gain_floor,
                'harmonic_weight': self.state.harmonic_weight,
                'noise_subtraction_alpha': self.state.noise_subtraction_alpha,
                'spectral_smoothing': self.state.spectral_smoothing,
            }
        }
