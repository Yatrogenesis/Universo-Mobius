"""
VMS Real-Time Processor
=======================

Low-latency audio processing for real-time applications.
Designed for integration with call centers, video conferencing,
and live streaming.
"""

import numpy as np
from typing import Optional, Callable
from collections import deque
from dataclasses import dataclass
import threading


@dataclass
class RealtimeStats:
    """Statistics for real-time processing."""
    frames_processed: int = 0
    total_latency_ms: float = 0.0
    avg_processing_time_ms: float = 0.0
    buffer_underruns: int = 0
    buffer_overruns: int = 0


class RealtimeProcessor:
    """
    Low-latency real-time audio processor.

    Uses circular buffers and optimized processing for
    minimal latency audio cleaning.

    Typical latency: 20-50ms depending on frame size.
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 1024,  # Smaller for lower latency
        hop_size: int = 256,
        noise_estimation_frames: int = 50,
        noise_floor_percentile: float = 25,
        harmonic_weight: float = 0.6,
        gain_floor: float = 0.15
    ):
        """
        Initialize real-time processor.

        Args:
            sample_rate: Audio sample rate
            frame_size: FFT frame size
            hop_size: Samples between frames
            noise_estimation_frames: Frames for rolling noise estimate
            noise_floor_percentile: Percentile for noise floor
            harmonic_weight: Harmonic preservation weight
            gain_floor: Minimum gain
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size
        self.noise_estimation_frames = noise_estimation_frames
        self.noise_floor_percentile = noise_floor_percentile
        self.harmonic_weight = harmonic_weight
        self.gain_floor = gain_floor

        # Windows
        self.analysis_window = np.sqrt(np.hanning(frame_size))
        self.synthesis_window = np.sqrt(np.hanning(frame_size))
        self.wola_window = np.hanning(frame_size)

        # Frequency bins
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

        # Circular buffers
        self.input_buffer = np.zeros(frame_size * 4)
        self.output_buffer = np.zeros(frame_size * 4)
        self.window_sum_buffer = np.zeros(frame_size * 4)
        self.buffer_pos = 0

        # Rolling noise estimation
        self.noise_spectra = deque(maxlen=noise_estimation_frames)
        self.noise_profile: Optional[np.ndarray] = None

        # Stats
        self.stats = RealtimeStats()

        # Thread safety
        self._lock = threading.Lock()

        # Voice parameters
        self.voice_fundamental = (80, 300)
        self.voice_ratios = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        self.enhanced_ratios = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]

    def reset(self):
        """Reset processor state."""
        with self._lock:
            self.input_buffer.fill(0)
            self.output_buffer.fill(0)
            self.window_sum_buffer.fill(0)
            self.buffer_pos = 0
            self.noise_spectra.clear()
            self.noise_profile = None
            self.stats = RealtimeStats()

    def _update_noise_profile(self, spectrum: np.ndarray):
        """Update rolling noise estimate."""
        self.noise_spectra.append(spectrum.copy())

        if len(self.noise_spectra) >= 10:
            stacked = np.vstack(list(self.noise_spectra))
            self.noise_profile = np.percentile(
                stacked, self.noise_floor_percentile, axis=0
            )

    def _detect_fundamental(self, spectrum: np.ndarray) -> float:
        """Fast F0 detection."""
        f_min, f_max = self.voice_fundamental
        freq_res = self.sample_rate / self.frame_size

        idx_min = int(f_min / freq_res)
        idx_max = min(int(f_max / freq_res), len(spectrum) - 1)

        if idx_max <= idx_min:
            return 150.0

        voice_spec = spectrum[idx_min:idx_max + 1]
        peak_idx = np.argmax(voice_spec)

        return (idx_min + peak_idx) * freq_res

    def _create_harmonic_mask(self, spectrum_len: int, f0: float) -> np.ndarray:
        """Create harmonic preservation mask."""
        mask = np.zeros(spectrum_len)
        freq_res = self.sample_rate / self.frame_size

        all_ratios = sorted(set(self.voice_ratios + self.enhanced_ratios))

        for ratio in all_ratios:
            f_harm = f0 * ratio
            if f_harm > self.sample_rate / 2:
                break

            idx = int(f_harm / freq_res)
            if idx >= spectrum_len:
                break

            width = max(2, int(f_harm / 100))  # Narrower for speed

            start = max(0, idx - width)
            end = min(spectrum_len, idx + width + 1)

            for i in range(start, end):
                dist = abs(i - idx)
                mask[i] = max(mask[i], np.exp(-dist**2 / (width**2)))

        return mask

    def process_chunk(self, chunk: np.ndarray) -> np.ndarray:
        """
        Process an audio chunk.

        Args:
            chunk: Input audio samples

        Returns:
            Processed audio samples (may be shorter due to latency)
        """
        import time
        start_time = time.perf_counter()

        with self._lock:
            chunk = np.asarray(chunk, dtype=np.float32)
            chunk_len = len(chunk)

            # Add to input buffer
            buf_len = len(self.input_buffer)
            new_pos = self.buffer_pos + chunk_len

            if new_pos <= buf_len:
                self.input_buffer[self.buffer_pos:new_pos] = chunk
            else:
                # Wrap around - shift buffer
                shift = new_pos - buf_len
                self.input_buffer[:-shift] = self.input_buffer[shift:]
                self.output_buffer[:-shift] = self.output_buffer[shift:]
                self.window_sum_buffer[:-shift] = self.window_sum_buffer[shift:]
                self.input_buffer[-chunk_len:] = chunk
                new_pos = buf_len
                self.buffer_pos = buf_len - chunk_len
                self.stats.buffer_overruns += 1

            # Process available frames
            frames_available = (new_pos - self.frame_size) // self.hop_size
            output_samples = []

            for _ in range(frames_available):
                frame_start = self.buffer_pos
                frame_end = frame_start + self.frame_size

                # Bounds check
                if frame_end > len(self.input_buffer):
                    break

                frame = self.input_buffer[frame_start:frame_end]

                # Process frame
                processed = self._process_frame(frame)

                # Overlap-add to output
                self.output_buffer[frame_start:frame_end] += processed
                self.window_sum_buffer[frame_start:frame_end] += self.wola_window

                self.buffer_pos += self.hop_size
                self.stats.frames_processed += 1

            # Extract output
            output_len = min(chunk_len, self.buffer_pos - self.frame_size // 2)
            if output_len > 0:
                start = max(0, self.buffer_pos - output_len - self.frame_size // 2)
                end = start + output_len

                # Normalize by window sum
                output = np.zeros(output_len)
                for i in range(output_len):
                    idx = start + i
                    if self.window_sum_buffer[idx] > 1e-8:
                        output[i] = self.output_buffer[idx] / self.window_sum_buffer[idx]

                output_samples = output

            self.buffer_pos = new_pos

        # Update stats
        processing_time = (time.perf_counter() - start_time) * 1000
        self.stats.avg_processing_time_ms = (
            self.stats.avg_processing_time_ms * 0.95 + processing_time * 0.05
        )
        self.stats.total_latency_ms = (self.frame_size / self.sample_rate) * 1000

        return np.array(output_samples) if len(output_samples) > 0 else np.array([])

    def _process_frame(self, frame: np.ndarray) -> np.ndarray:
        """Process single frame."""
        if len(frame) < self.frame_size:
            frame = np.pad(frame, (0, self.frame_size - len(frame)))

        # Analysis
        windowed = frame * self.analysis_window
        fft_data = np.fft.rfft(windowed)
        spectrum = np.abs(fft_data)
        phase = np.angle(fft_data)

        # Update noise estimate
        self._update_noise_profile(spectrum)

        # Get noise estimate
        if self.noise_profile is None:
            noise_estimate = np.full_like(spectrum, np.median(spectrum) * 0.5)
        else:
            noise_estimate = self.noise_profile

        # Wiener filtering
        spectrum_safe = np.maximum(spectrum, 1e-10)
        noise_safe = np.maximum(noise_estimate, 1e-10)

        signal_power = spectrum_safe ** 2
        noise_power = noise_safe ** 2

        gain_sq = np.maximum(signal_power - noise_power, 0) / signal_power
        gain = np.sqrt(gain_sq)
        gain = np.maximum(gain, self.gain_floor)

        # Harmonic preservation
        f0 = self._detect_fundamental(spectrum)
        harmonic_mask = self._create_harmonic_mask(len(spectrum), f0)

        subtracted = spectrum * gain
        harmonic_boost = 1 + (1 - gain) * harmonic_mask * self.harmonic_weight
        cleaned_spectrum = np.minimum(subtracted * harmonic_boost, spectrum)

        # Reconstruct
        reconstructed = cleaned_spectrum * np.exp(1j * phase)
        cleaned_frame = np.fft.irfft(reconstructed, n=self.frame_size)

        return cleaned_frame * self.synthesis_window

    def get_latency_ms(self) -> float:
        """Get current processing latency in milliseconds."""
        return (self.frame_size / self.sample_rate) * 1000

    def get_stats(self) -> RealtimeStats:
        """Get processing statistics."""
        return self.stats


class StreamingCleaner:
    """
    High-level streaming interface for VMS audio cleaning.

    Suitable for integration with audio frameworks like
    PyAudio, sounddevice, or web audio APIs.
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        channels: int = 1,
        **processor_kwargs
    ):
        """
        Initialize streaming cleaner.

        Args:
            sample_rate: Audio sample rate
            channels: Number of audio channels (1=mono, 2=stereo)
            **processor_kwargs: Arguments passed to RealtimeProcessor
        """
        self.sample_rate = sample_rate
        self.channels = channels

        # Create processor for each channel
        self.processors = [
            RealtimeProcessor(sample_rate=sample_rate, **processor_kwargs)
            for _ in range(channels)
        ]

    def process(self, audio: np.ndarray) -> np.ndarray:
        """
        Process audio (mono or stereo).

        Args:
            audio: Input audio, shape (samples,) or (samples, channels)

        Returns:
            Processed audio, same shape as input
        """
        if audio.ndim == 1:
            # Mono
            return self.processors[0].process_chunk(audio)

        else:
            # Multi-channel
            output = np.zeros_like(audio)
            for ch in range(min(audio.shape[1], self.channels)):
                processed = self.processors[ch].process_chunk(audio[:, ch])
                # Handle case where processor returns shorter array during warmup
                if len(processed) == 0:
                    continue  # Leave zeros
                elif len(processed) < audio.shape[0]:
                    output[:len(processed), ch] = processed
                else:
                    output[:, ch] = processed[:audio.shape[0]]
            return output

    def reset(self):
        """Reset all processors."""
        for p in self.processors:
            p.reset()

    def get_latency_ms(self) -> float:
        """Get processing latency."""
        return self.processors[0].get_latency_ms()
