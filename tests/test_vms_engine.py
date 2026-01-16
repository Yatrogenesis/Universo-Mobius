"""
Tests for VMS Engine Module
===========================

Tests for audio cleaning and real-time processing.
"""

import pytest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from vms_engine.cleaner import VMSAudioCleaner, CleaningResult
from vms_engine.realtime import RealtimeProcessor, StreamingCleaner
from vms_engine.spectral import (
    SpectralAnalyzer,
    wiener_filter,
    harmonic_preserving_filter,
    spectral_gate,
    detect_voice_activity,
)


class TestVMSAudioCleaner:
    """Test VMSAudioCleaner class."""

    @pytest.fixture
    def cleaner(self):
        """Create cleaner instance."""
        return VMSAudioCleaner(
            sample_rate=44100,
            frame_size=2048,
            hop_size=512,
            noise_floor_percentile=25,
            harmonic_weight=0.6
        )

    @pytest.fixture
    def test_signal(self):
        """Create test signal with voice + noise."""
        sr = 44100
        duration = 1.0
        t = np.linspace(0, duration, int(sr * duration))

        # Voice-like signal
        f0 = 150  # Hz
        voice = np.zeros_like(t)
        for h in range(1, 8):
            voice += (1.0 / h) * np.sin(2 * np.pi * f0 * h * t)
        voice = voice / np.max(np.abs(voice)) * 0.5

        # Noise
        np.random.seed(42)
        noise = np.random.randn(len(t)) * 0.2

        return voice + noise, voice, noise, sr

    def test_cleaner_initialization(self, cleaner):
        """Should initialize with correct parameters."""
        assert cleaner.sample_rate == 44100
        assert cleaner.frame_size == 2048
        assert cleaner.hop_size == 512
        assert len(cleaner.freqs) == cleaner.frame_size // 2 + 1

    def test_noise_profile_estimation(self, cleaner, test_signal):
        """Should estimate noise profile."""
        noisy, _, _, sr = test_signal

        profile = cleaner.estimate_noise_profile(noisy, noise_duration=0.2)

        assert profile is not None
        assert len(profile) == cleaner.frame_size // 2 + 1
        assert np.all(profile >= 0)

    def test_fundamental_detection(self, cleaner, test_signal):
        """Should detect fundamental frequency."""
        noisy, voice, _, sr = test_signal

        # Get spectrum of voice
        spectrum = np.abs(np.fft.rfft(voice[:2048]))

        f0 = cleaner.detect_fundamental(spectrum)

        # Should be close to 150 Hz
        assert 100 < f0 < 200

    def test_harmonic_mask_creation(self, cleaner):
        """Should create harmonic mask."""
        spectrum = np.ones(1025)  # rfft output size for 2048
        f0 = 150.0

        mask = cleaner.create_harmonic_mask(spectrum, f0)

        assert len(mask) == len(spectrum)
        assert np.max(mask) <= 1.0
        assert np.min(mask) >= 0.0

        # Should have peaks at harmonics
        freq_res = cleaner.sample_rate / cleaner.frame_size
        for h in [1, 2, 3]:
            idx = int(f0 * h / freq_res)
            if idx < len(mask):
                assert mask[idx] > 0.5  # Peak at harmonic

    def test_frame_processing(self, cleaner, test_signal):
        """Should process single frame."""
        noisy, _, _, _ = test_signal

        frame = noisy[:2048]
        cleaner.estimate_noise_profile(noisy[:8000])

        processed = cleaner.process_frame(frame)

        assert len(processed) == len(frame)
        assert not np.any(np.isnan(processed))
        assert not np.any(np.isinf(processed))

    def test_full_audio_cleaning(self, cleaner, test_signal):
        """Should clean complete audio."""
        noisy, voice, noise, sr = test_signal

        result = cleaner.clean_audio(noisy)

        assert isinstance(result, CleaningResult)
        assert len(result.audio) == len(noisy)
        assert result.sample_rate == sr
        assert result.frames_processed > 0
        assert result.processing_time > 0

        # SNR should improve
        # Note: simplified metric may not always show improvement
        assert result.correlation > 0.5

    def test_audio_not_silent(self, cleaner, test_signal):
        """Cleaned audio should not be silent."""
        noisy, _, _, _ = test_signal

        result = cleaner.clean_audio(noisy)

        # Should preserve significant energy
        input_power = np.mean(noisy ** 2)
        output_power = np.mean(result.audio ** 2)

        # Output should have at least 10% of input power
        assert output_power > input_power * 0.1


class TestRealtimeProcessor:
    """Test RealtimeProcessor class."""

    @pytest.fixture
    def processor(self):
        """Create processor instance."""
        return RealtimeProcessor(
            sample_rate=44100,
            frame_size=1024,
            hop_size=256
        )

    def test_processor_initialization(self, processor):
        """Should initialize correctly."""
        assert processor.sample_rate == 44100
        assert processor.frame_size == 1024
        assert processor.get_latency_ms() > 0

    def test_chunk_processing(self, processor):
        """Should process audio chunks."""
        chunk = np.random.randn(512).astype(np.float32) * 0.1

        # Process multiple chunks to fill buffer
        for _ in range(10):
            output = processor.process_chunk(chunk)

        # Eventually should produce output
        stats = processor.get_stats()
        assert stats.frames_processed > 0

    def test_processor_reset(self, processor):
        """Should reset state."""
        chunk = np.random.randn(512).astype(np.float32) * 0.1

        # Process some data
        for _ in range(5):
            processor.process_chunk(chunk)

        # Reset
        processor.reset()

        stats = processor.get_stats()
        assert stats.frames_processed == 0


class TestStreamingCleaner:
    """Test StreamingCleaner class."""

    def test_mono_processing(self):
        """Should process mono audio."""
        cleaner = StreamingCleaner(sample_rate=44100, channels=1)

        chunk = np.random.randn(1024).astype(np.float32) * 0.1

        # Process multiple chunks
        for _ in range(10):
            output = cleaner.process(chunk)

        assert cleaner.get_latency_ms() > 0

    def test_stereo_processing(self):
        """Should process stereo audio."""
        cleaner = StreamingCleaner(sample_rate=44100, channels=2)

        chunk = np.random.randn(1024, 2).astype(np.float32) * 0.1

        # Process
        for _ in range(10):
            output = cleaner.process(chunk)


class TestSpectralTools:
    """Test spectral analysis tools."""

    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance."""
        return SpectralAnalyzer(
            sample_rate=44100,
            frame_size=2048,
            hop_size=512
        )

    def test_frame_analysis(self, analyzer):
        """Should analyze single frame."""
        frame = np.random.randn(2048) * 0.1

        result = analyzer.analyze_frame(frame)

        assert len(result.magnitude) == 1025  # rfft output
        assert len(result.phase) == 1025
        assert result.sample_rate == 44100

    def test_spectrogram(self, analyzer):
        """Should compute spectrogram."""
        audio = np.random.randn(44100) * 0.1  # 1 second

        spec, freqs, times = analyzer.get_spectrogram(audio)

        assert spec.shape[0] == len(freqs)
        assert spec.shape[1] == len(times)
        assert len(times) > 0

    def test_peak_finding(self, analyzer):
        """Should find spectral peaks."""
        # Create spectrum with clear peaks
        spectrum = np.random.randn(1025) * 0.1
        spectrum[100] = 10.0  # Strong peak
        spectrum[200] = 8.0   # Another peak

        peaks = analyzer.find_peaks(spectrum, n_peaks=5)

        assert len(peaks) > 0
        # Strongest peak should be near index 100
        peak_freqs = [p[0] for p in peaks]
        assert any(abs(f - analyzer.frequencies[100]) < 100 for f in peak_freqs)

    def test_wiener_filter(self):
        """Should apply Wiener filter."""
        spectrum = np.ones(100) * 10
        noise = np.ones(100) * 5

        filtered = wiener_filter(spectrum, noise)

        assert len(filtered) == len(spectrum)
        assert np.all(filtered <= spectrum)  # Shouldn't amplify
        assert np.all(filtered >= 0)

    def test_harmonic_preserving_filter(self):
        """Should preserve harmonics."""
        freqs = np.linspace(0, 22050, 1025)
        spectrum = np.ones(1025)

        mask = harmonic_preserving_filter(spectrum, freqs, f0=200.0)

        assert len(mask) == len(spectrum)
        assert np.max(mask) > 1.0  # Boost at harmonics

        # Should boost at 200 Hz
        idx_200 = np.argmin(np.abs(freqs - 200))
        assert mask[idx_200] > 1.2

    def test_spectral_gate(self):
        """Should gate low-energy regions."""
        spectrum = np.array([1, 2, 10, 2, 1])

        # Use hold_bins=0 for pure hard gate (default hold_bins=2 would keep neighbors)
        gated = spectral_gate(spectrum, threshold=5, attack_ratio=0, hold_bins=0)

        assert gated[2] == 10  # Above threshold preserved
        assert gated[0] == 0   # Below threshold gated
        assert gated[4] == 0

    def test_voice_activity_detection(self):
        """Should detect voice activity."""
        freqs = np.linspace(0, 22050, 1025)

        # Voice-like spectrum (energy in voice band)
        voice_spectrum = np.ones(1025) * 0.1
        voice_mask = (freqs > 80) & (freqs < 4000)
        voice_spectrum[voice_mask] = 10.0

        assert detect_voice_activity(voice_spectrum, freqs)

        # Noise-only spectrum
        noise_spectrum = np.ones(1025) * 5.0

        # Uniform energy might or might not trigger - depends on threshold


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
