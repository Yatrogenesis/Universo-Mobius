"""
GeoVMS Laser Microphone Signal Recovery
========================================

Extract voice/audio from laser microphone measurements of surface vibrations.

How laser microphones work:
1. Laser beam directed at glass window/surface
2. Sound waves cause micro-displacements (nm to um)
3. Reflected beam shows Doppler shift / interferometric pattern
4. Demodulate to recover acoustic signal

Challenges:
- Extremely low SNR (environmental vibrations dominate)
- Non-linear response (surface curvature, angle)
- Multi-path reflections from room interior
- Security/counter-surveillance applications

Why VMS Topological approach excels here:
- Voice modes form coherent patterns in solid medium
- Building vibrations are broad, incoherent
- Can separate direct sound from structural resonances
- Topological "mountains" = coherent voice, "dust" = noise

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
from typing import Optional, Tuple, Dict, List
from dataclasses import dataclass
from scipy import signal as scipy_signal
from scipy import ndimage


@dataclass
class LaserMicResult:
    """Result from laser microphone demodulation."""
    audio: np.ndarray
    sample_rate: int
    snr_improvement_db: float
    voice_activity_ratio: float
    structural_modes_removed: int
    dominant_voice_freq: float


class LaserMicrophoneProcessor:
    """
    Process laser microphone signals using topological methods.

    The key insight is that vibrations travel through solid materials
    (glass, walls) with coherent modal structure. Voice creates
    specific patterns that can be separated from building noise.

    Laser microphone raw signal typically contains:
    1. Voice-induced surface displacement (target signal)
    2. HVAC vibrations (10-100 Hz, periodic)
    3. Traffic/footsteps (1-20 Hz, impulsive)
    4. Electronic noise (50/60 Hz harmonics)
    5. Wind on window (broadband)

    VMS topological separation:
    - Voice: Connected frequency bands (harmonics)
    - HVAC: Narrow, periodic, isolated
    - Traffic: Low-frequency blobs
    - Wind: Scattered, no structure

    Example:
        proc = LaserMicrophoneProcessor(sample_rate=44100)
        result = proc.recover_voice(raw_laser_signal)
        # result.audio contains recovered voice
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512,
        voice_range: Tuple[float, float] = (80, 4000),
        window_resonance_range: Tuple[float, float] = (100, 500)
    ):
        """
        Initialize processor.

        Args:
            sample_rate: Sample rate in Hz
            frame_size: STFT frame size
            hop_size: STFT hop size
            voice_range: Voice frequency range (Hz)
            window_resonance_range: Expected glass resonance range
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size
        self.voice_range = voice_range
        self.window_resonance_range = window_resonance_range

        # Windows
        self.window = np.hanning(frame_size)
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Voice harmonic ratios (fundamental + overtones)
        self.voice_ratios = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]

        # OCTH hexagonal ratios for enhanced detection
        self.hex_ratios = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def _compute_spectrogram(
        self,
        audio: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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

    def _remove_powerline_harmonics(
        self,
        magnitude: np.ndarray,
        powerline_freq: float = 60.0,
        n_harmonics: int = 10,
        width_hz: float = 2.0
    ) -> np.ndarray:
        """Remove powerline interference (50/60 Hz and harmonics)."""
        cleaned = magnitude.copy()
        freq_res = self.sample_rate / self.frame_size

        for h in range(1, n_harmonics + 1):
            f_harm = powerline_freq * h
            idx = int(f_harm / freq_res)
            width_bins = int(width_hz / freq_res)

            if idx < len(self.freqs):
                start = max(0, idx - width_bins)
                end = min(len(self.freqs), idx + width_bins + 1)

                # Interpolate across notch
                if start > 0 and end < len(self.freqs):
                    left_val = cleaned[start - 1, :]
                    right_val = cleaned[end, :]
                    for i in range(start, end):
                        alpha = (i - start) / (end - start)
                        cleaned[i, :] = (1 - alpha) * left_val + alpha * right_val

        return cleaned

    def _identify_structural_modes(
        self,
        magnitude: np.ndarray
    ) -> Tuple[np.ndarray, List[float]]:
        """
        Identify and remove structural resonance modes.

        Glass windows have characteristic resonances based on:
        - Size
        - Thickness
        - Mounting
        - Temperature

        These appear as persistent horizontal lines in spectrogram.
        """
        # Average spectrum (persistent = structural)
        avg_spectrum = np.mean(magnitude, axis=1)

        # Find peaks in average (structural modes)
        peaks, properties = scipy_signal.find_peaks(
            avg_spectrum,
            height=np.mean(avg_spectrum) * 2,
            distance=10
        )

        # Filter to window resonance range
        freq_res = self.sample_rate / self.frame_size
        structural_freqs = []
        mask = np.ones_like(magnitude)

        for peak in peaks:
            f = peak * freq_res
            if self.window_resonance_range[0] <= f <= self.window_resonance_range[1]:
                structural_freqs.append(f)
                # Create notch filter
                width = max(3, int(f / 50))  # Wider at higher freq
                start = max(0, peak - width)
                end = min(len(self.freqs), peak + width + 1)
                mask[start:end, :] = 0.1  # Don't fully remove

        cleaned = magnitude * mask
        return cleaned, structural_freqs

    def _detect_voice_f0(
        self,
        magnitude: np.ndarray
    ) -> float:
        """Detect fundamental frequency of voice."""
        freq_res = self.sample_rate / self.frame_size

        # Voice range
        idx_min = int(self.voice_range[0] / freq_res)
        idx_max = int(min(400, self.voice_range[1]) / freq_res)

        # Time-averaged spectrum
        avg_spectrum = np.mean(magnitude, axis=1)

        if idx_max <= idx_min:
            return 150.0

        voice_spec = avg_spectrum[idx_min:idx_max + 1]
        peak_idx = np.argmax(voice_spec)

        return (idx_min + peak_idx) * freq_res

    def _create_voice_harmonic_mask(
        self,
        f0: float,
        n_freqs: int
    ) -> np.ndarray:
        """Create mask preserving voice harmonics."""
        mask = np.zeros(n_freqs)
        freq_res = self.sample_rate / self.frame_size

        # Combine standard + hexagonal ratios
        all_ratios = sorted(set(self.voice_ratios + self.hex_ratios))

        for ratio in all_ratios:
            f_harm = f0 * ratio
            if f_harm > self.voice_range[1]:
                break

            idx = int(f_harm / freq_res)
            if idx >= n_freqs:
                break

            # Gaussian window around harmonic
            width = max(3, int(f_harm / 100))
            for i in range(max(0, idx - width), min(n_freqs, idx + width + 1)):
                dist = abs(i - idx)
                mask[i] = max(mask[i], np.exp(-dist**2 / (width**2)))

        return mask

    def _topological_voice_separation(
        self,
        magnitude: np.ndarray,
        min_component_size: int = 30
    ) -> Tuple[np.ndarray, int]:
        """
        Separate voice using topological connected components.

        Voice in solid medium forms coherent "ridges" connecting
        fundamental to harmonics. Noise is scattered.
        """
        # Convert to dB
        mag_db = 20 * np.log10(magnitude + 1e-10)

        # Adaptive threshold
        threshold = np.percentile(mag_db, 40)
        binary = mag_db > threshold

        # Voice frequency mask
        freq_res = self.sample_rate / self.frame_size
        freq_mask = np.zeros(len(self.freqs), dtype=bool)
        for i, f in enumerate(self.freqs):
            if self.voice_range[0] <= f <= self.voice_range[1]:
                freq_mask[i] = True

        binary = binary & freq_mask[:, np.newaxis]

        # Find connected components
        labeled, n_components = ndimage.label(binary)

        # Keep large components (voice) remove small (noise)
        output_mask = np.zeros_like(magnitude)
        large_count = 0

        for i in range(1, n_components + 1):
            component = labeled == i
            size = np.sum(component)

            if size >= min_component_size:
                large_count += 1
                output_mask[component] = 1.0

        # Smooth mask edges
        output_mask = ndimage.gaussian_filter(output_mask.astype(float), sigma=1)

        # Apply with soft floor
        separated = magnitude * (0.05 + 0.95 * output_mask)

        return separated, n_components - large_count

    def _voice_activity_detection(
        self,
        magnitude: np.ndarray,
        threshold_ratio: float = 3.0
    ) -> np.ndarray:
        """
        Detect frames with voice activity.

        Returns binary mask for time frames.
        """
        freq_res = self.sample_rate / self.frame_size

        # Voice band energy
        voice_idx_min = int(self.voice_range[0] / freq_res)
        voice_idx_max = int(self.voice_range[1] / freq_res)

        voice_energy = np.sum(magnitude[voice_idx_min:voice_idx_max, :], axis=0)

        # Noise band energy (very low freq)
        noise_energy = np.sum(magnitude[:voice_idx_min, :], axis=0)

        # Ratio
        ratio = voice_energy / (noise_energy + 1e-10)

        # Binary decision
        vad = ratio > threshold_ratio

        # Smooth (remove isolated detections)
        vad = ndimage.binary_opening(vad, structure=np.ones(3))
        vad = ndimage.binary_closing(vad, structure=np.ones(5))

        return vad

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

            spec = magnitude[:, i] * phase[:, i]
            frame = np.fft.irfft(spec, n=self.frame_size)

            output[start:start + self.frame_size] += frame * self.window
            window_sum[start:start + self.frame_size] += self.window ** 2

        window_sum = np.maximum(window_sum, 1e-8)
        output = output / window_sum

        return output[:length]

    def recover_voice(
        self,
        laser_signal: np.ndarray,
        powerline_freq: float = 60.0,
        aggressive: bool = False
    ) -> LaserMicResult:
        """
        Recover voice from laser microphone signal.

        Args:
            laser_signal: Raw laser microphone signal
            powerline_freq: Power line frequency (50 or 60 Hz)
            aggressive: Use more aggressive noise removal

        Returns:
            LaserMicResult with recovered audio
        """
        signal = np.asarray(laser_signal, dtype=np.float64)
        original_length = len(signal)

        # Normalize
        signal = signal - np.mean(signal)
        if np.std(signal) > 0:
            signal = signal / np.std(signal)

        # Compute spectrogram
        magnitude, phase, times = self._compute_spectrogram(signal)

        # Step 1: Remove powerline harmonics
        magnitude = self._remove_powerline_harmonics(
            magnitude, powerline_freq
        )

        # Step 2: Identify and attenuate structural modes
        magnitude, structural_freqs = self._identify_structural_modes(magnitude)

        # Step 3: Detect voice F0
        f0 = self._detect_voice_f0(magnitude)

        # Step 4: Create voice harmonic mask
        voice_mask = self._create_voice_harmonic_mask(f0, len(self.freqs))

        # Step 5: Topological voice separation
        separated, noise_components = self._topological_voice_separation(
            magnitude,
            min_component_size=50 if aggressive else 30
        )

        # Step 6: Blend topological + harmonic
        harmonic_enhanced = separated * (1 + voice_mask[:, np.newaxis] * 0.5)

        # Step 7: Voice activity detection
        vad = self._voice_activity_detection(magnitude)
        vad_2d = vad[np.newaxis, :].repeat(len(self.freqs), axis=0)

        # Apply VAD (softer in non-speech regions)
        final_magnitude = harmonic_enhanced * (0.1 + 0.9 * vad_2d)

        # Reconstruct
        recovered = self._reconstruct_audio(final_magnitude, phase, original_length)

        # Prevent clipping
        if np.max(np.abs(recovered)) > 0.95:
            recovered = recovered / np.max(np.abs(recovered)) * 0.95

        # Calculate metrics
        input_power = np.mean(signal ** 2)
        output_power = np.mean(recovered ** 2)
        snr_improvement = 10 * np.log10(output_power / (input_power + 1e-10))

        voice_ratio = np.mean(vad.astype(float))

        return LaserMicResult(
            audio=recovered.astype(np.float32),
            sample_rate=self.sample_rate,
            snr_improvement_db=snr_improvement,
            voice_activity_ratio=voice_ratio,
            structural_modes_removed=len(structural_freqs),
            dominant_voice_freq=f0
        )


def demo_laser_mic():
    """Demonstrate laser microphone processing."""
    print("=" * 60)
    print("GeoVMS Laser Microphone Signal Recovery Demo")
    print("=" * 60)

    sample_rate = 44100
    duration = 5  # seconds
    t = np.linspace(0, duration, int(duration * sample_rate))

    # Simulate voice (target signal)
    f0 = 150  # Fundamental frequency
    voice = np.zeros_like(t)
    for h in range(1, 8):
        voice += np.sin(2 * np.pi * f0 * h * t) / h

    # Add formants (voice characteristic)
    formant1 = scipy_signal.lfilter(*scipy_signal.butter(2, [500/(sample_rate/2), 700/(sample_rate/2)], 'band'), voice)
    voice = voice + formant1 * 0.5

    # Simulate speech envelope (not continuous)
    envelope = np.zeros_like(t)
    for start in np.arange(0.5, duration - 0.3, 0.8):
        envelope += np.exp(-((t - start) / 0.2)**2)
    voice = voice * np.minimum(envelope, 1)

    # Structural resonances (window glass ~200 Hz)
    structural = np.sin(2 * np.pi * 180 * t) * 0.5
    structural += np.sin(2 * np.pi * 350 * t) * 0.3

    # HVAC noise (low frequency)
    hvac = np.sin(2 * np.pi * 20 * t) * 0.3

    # Power line interference
    powerline = np.sin(2 * np.pi * 60 * t) * 0.2
    powerline += np.sin(2 * np.pi * 120 * t) * 0.1

    # Wind noise (broadband)
    wind = np.random.randn(len(t)) * 0.1
    wind = scipy_signal.lfilter(*scipy_signal.butter(2, 100/(sample_rate/2), 'low'), wind)

    # Combine (very noisy!)
    raw_signal = voice * 0.1 + structural + hvac + powerline + wind

    print(f"Simulated laser mic signal: {duration}s at {sample_rate} Hz")
    print(f"Voice at F0={f0} Hz, heavily contaminated with:")
    print(f"  - Structural resonances (180, 350 Hz)")
    print(f"  - HVAC (20 Hz)")
    print(f"  - Powerline (60, 120 Hz)")
    print(f"  - Wind noise")

    # Process
    proc = LaserMicrophoneProcessor(sample_rate=sample_rate)
    result = proc.recover_voice(raw_signal)

    print(f"\nRecovery Results:")
    print(f"  SNR improvement: {result.snr_improvement_db:.1f} dB")
    print(f"  Voice activity: {result.voice_activity_ratio*100:.1f}%")
    print(f"  Structural modes removed: {result.structural_modes_removed}")
    print(f"  Detected F0: {result.dominant_voice_freq:.1f} Hz (actual: {f0})")

    # Correlation with original voice
    min_len = min(len(voice), len(result.audio))
    correlation = np.corrcoef(voice[:min_len], result.audio[:min_len])[0, 1]
    print(f"  Voice correlation: {correlation:.4f}")


if __name__ == "__main__":
    demo_laser_mic()
