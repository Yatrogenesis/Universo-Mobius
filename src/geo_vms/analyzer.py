"""
GeoVMS Seismic Topology Analyzer
================================

Apply VMS topological methods to seismic data.

Key insight: Seismic waves maintain coherent structure in solid media.
Unlike acoustic waves in air, seismic waves form distinct "arrivals"
that appear as connected ridges in time-frequency space.

Methods:
1. Topological phase picking (P, S, surface waves)
2. Coda wave analysis (scattering structure)
3. Source characterization from spectral topology
4. Multi-station coherence analysis

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
from typing import Optional, Tuple, Dict, List, Any
from dataclasses import dataclass
from scipy import ndimage
from scipy import signal as scipy_signal


@dataclass
class PhaseArrival:
    """A seismic phase arrival."""
    phase_name: str      # P, S, PP, SS, PcP, ScS, Love, Rayleigh, etc.
    arrival_time: float  # Seconds from trace start
    confidence: float    # 0-1
    frequency_band: Tuple[float, float]  # Hz
    component_id: int    # Topological component ID


@dataclass
class SeismicTopologyResult:
    """Result from topological seismic analysis."""
    phase_arrivals: List[PhaseArrival]
    n_coherent_components: int
    n_scattered_components: int
    coda_decay_time: float  # Seconds
    dominant_frequency: float  # Hz
    spectral_width: float
    topology_map: np.ndarray  # Labeled connected components
    signal_strength_db: float


class SeismicTopologyAnalyzer:
    """
    Topological analyzer for seismic signals.

    Seismic waves are ideal for topological analysis because:
    1. P-waves and S-waves form distinct time-frequency ridges
    2. Surface waves (Love, Rayleigh) have characteristic dispersion
    3. Scattered coda can be separated from direct arrivals
    4. Multi-pathing creates connected but distinct structures

    Unlike audio noise reduction (where we separate signal/noise),
    seismic analysis separates DIFFERENT TYPES of signal:
    - Direct arrivals vs. scattered coda
    - Body waves vs. surface waves
    - Source effects vs. path effects

    Example:
        analyzer = SeismicTopologyAnalyzer(sample_rate=40.0)
        result = analyzer.analyze(seismic_trace)
        for arrival in result.phase_arrivals:
            print(f"{arrival.phase_name} at {arrival.arrival_time:.1f}s")
    """

    def __init__(
        self,
        sample_rate: float = 40.0,
        frame_size: int = 256,
        hop_size: int = 64,
        freq_bands: Dict[str, Tuple[float, float]] = None
    ):
        """
        Initialize analyzer.

        Args:
            sample_rate: Sample rate in Hz
            frame_size: STFT frame size
            hop_size: STFT hop size
            freq_bands: Frequency bands for different wave types
        """
        self.sample_rate = sample_rate
        self.frame_size = frame_size
        self.hop_size = hop_size

        # Default seismic frequency bands
        if freq_bands is None:
            freq_bands = {
                "teleseismic_p": (0.5, 2.0),    # Teleseismic P-waves
                "teleseismic_s": (0.05, 0.5),   # Teleseismic S-waves
                "local_p": (1.0, 10.0),         # Local P-waves
                "local_s": (0.5, 5.0),          # Local S-waves
                "surface": (0.01, 0.1),         # Surface waves
                "microseism": (0.1, 0.3),       # Ocean microseism
                "high_freq": (5.0, 20.0),       # Local/near sources
            }
        self.freq_bands = freq_bands

        # Analysis window
        self.window = np.hanning(frame_size)
        self.freqs = np.fft.rfftfreq(frame_size, 1/sample_rate)

    def _compute_spectrogram(
        self,
        data: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute spectrogram."""
        n_frames = (len(data) - self.frame_size) // self.hop_size + 1
        n_freqs = self.frame_size // 2 + 1

        magnitude = np.zeros((n_freqs, n_frames))
        phase = np.zeros((n_freqs, n_frames))

        for i in range(n_frames):
            start = i * self.hop_size
            frame = data[start:start + self.frame_size]
            if len(frame) < self.frame_size:
                frame = np.pad(frame, (0, self.frame_size - len(frame)))

            windowed = frame * self.window
            spec = np.fft.rfft(windowed)
            magnitude[:, i] = np.abs(spec)
            phase[:, i] = np.angle(spec)

        times = np.arange(n_frames) * self.hop_size / self.sample_rate
        return magnitude, phase, times

    def _find_connected_components(
        self,
        magnitude: np.ndarray,
        threshold_percentile: float = 50,
        min_size: int = 20
    ) -> Tuple[np.ndarray, int, List[Dict]]:
        """
        Find connected components (topological features) in spectrogram.

        Returns:
            labeled: Array with component labels
            n_components: Number of large components
            component_info: List of component properties
        """
        # Convert to dB
        mag_db = 20 * np.log10(magnitude + 1e-10)

        # Dynamic threshold based on local statistics
        threshold = np.percentile(mag_db, threshold_percentile)
        binary = mag_db > threshold

        # Find connected components
        labeled, n_total = ndimage.label(binary)

        # Measure properties of each component
        component_info = []
        n_large = 0

        for i in range(1, n_total + 1):
            mask = labeled == i
            size = np.sum(mask)

            if size >= min_size:
                n_large += 1

                # Find component extent
                rows, cols = np.where(mask)
                freq_indices = rows
                time_indices = cols

                # Calculate centroid and extent
                freq_centroid = np.mean(freq_indices)
                time_centroid = np.mean(time_indices)
                freq_extent = np.max(freq_indices) - np.min(freq_indices)
                time_extent = np.max(time_indices) - np.min(time_indices)

                # Average amplitude in component
                avg_amplitude = np.mean(mag_db[mask])

                # Convert to physical units
                freq_hz = freq_centroid * self.sample_rate / self.frame_size
                time_sec = time_centroid * self.hop_size / self.sample_rate

                component_info.append({
                    'id': i,
                    'size': size,
                    'freq_centroid_hz': freq_hz,
                    'time_centroid_sec': time_sec,
                    'freq_extent_hz': freq_extent * self.sample_rate / self.frame_size,
                    'time_extent_sec': time_extent * self.hop_size / self.sample_rate,
                    'amplitude_db': avg_amplitude,
                    'aspect_ratio': time_extent / (freq_extent + 1)  # Time/Freq
                })

        return labeled, n_large, component_info

    def _classify_components(
        self,
        components: List[Dict]
    ) -> List[PhaseArrival]:
        """
        Classify topological components as seismic phases.

        Uses:
        1. Frequency band (P vs S vs surface)
        2. Aspect ratio (impulsive vs dispersive)
        3. Arrival time (early = P, later = S, latest = surface)
        4. Amplitude pattern
        """
        arrivals = []

        # Sort by arrival time
        sorted_components = sorted(components, key=lambda c: c['time_centroid_sec'])

        for i, comp in enumerate(sorted_components):
            freq = comp['freq_centroid_hz']
            aspect = comp['aspect_ratio']
            time = comp['time_centroid_sec']
            amp = comp['amplitude_db']

            # Classification logic
            if freq > 1.0 and aspect < 5:
                # High frequency, short duration → P-wave
                phase_name = "P" if i == 0 else "PP"
                confidence = 0.8

            elif 0.2 < freq < 2.0 and 2 < aspect < 20:
                # Medium frequency, moderate duration → S-wave
                phase_name = "S" if i <= 1 else "SS"
                confidence = 0.7

            elif freq < 0.2 and aspect > 20:
                # Low frequency, long duration → Surface wave
                if aspect > 50:
                    phase_name = "Rayleigh"
                else:
                    phase_name = "Love"
                confidence = 0.6

            elif freq < 0.3 and 0.1 < freq:
                # Microseism band
                phase_name = "Microseism"
                confidence = 0.5

            else:
                phase_name = f"Unknown_{i}"
                confidence = 0.3

            # Determine frequency band
            freq_band = (
                max(0, freq - comp['freq_extent_hz']/2),
                freq + comp['freq_extent_hz']/2
            )

            arrivals.append(PhaseArrival(
                phase_name=phase_name,
                arrival_time=time,
                confidence=confidence,
                frequency_band=freq_band,
                component_id=comp['id']
            ))

        return arrivals

    def _analyze_coda(
        self,
        data: np.ndarray,
        s_arrival: float = None
    ) -> Tuple[float, float]:
        """
        Analyze coda (scattered wave train after S-wave).

        Returns:
            decay_time: Coda Q decay time constant
            coda_q: Quality factor
        """
        if s_arrival is None:
            s_arrival = len(data) * 0.3 / self.sample_rate

        # Get coda portion (after S-wave)
        s_sample = int(s_arrival * self.sample_rate)
        coda = data[s_sample:]

        if len(coda) < 100:
            return 0.0, 0.0

        # Compute envelope
        analytic = scipy_signal.hilbert(coda)
        envelope = np.abs(analytic)

        # Smooth envelope
        window_size = min(100, len(envelope) // 10)
        if window_size > 0:
            envelope = np.convolve(envelope, np.ones(window_size)/window_size, mode='same')

        # Find decay time (time to drop to 1/e)
        max_idx = np.argmax(envelope)
        max_val = envelope[max_idx]
        target = max_val / np.e

        decay_idx = None
        for i in range(max_idx, len(envelope)):
            if envelope[i] < target:
                decay_idx = i
                break

        if decay_idx is None:
            decay_time = (len(envelope) - max_idx) / self.sample_rate
        else:
            decay_time = (decay_idx - max_idx) / self.sample_rate

        # Coda Q (quality factor)
        # For frequency f: Q_c = 2 * pi * f * decay_time
        dominant_freq = 1.0  # Assume 1 Hz dominant
        coda_q = 2 * np.pi * dominant_freq * decay_time

        return decay_time, coda_q

    def analyze(
        self,
        data: np.ndarray,
        threshold_percentile: float = 50,
        min_component_size: int = 20
    ) -> SeismicTopologyResult:
        """
        Perform complete topological analysis of seismic trace.

        Args:
            data: Seismic waveform data
            threshold_percentile: Threshold for component detection
            min_component_size: Minimum pixels for a component

        Returns:
            SeismicTopologyResult with arrivals and statistics
        """
        data = np.asarray(data, dtype=np.float64)

        # Detrend and normalize
        data = data - np.mean(data)
        if np.std(data) > 0:
            data = data / np.std(data)

        # Compute spectrogram
        magnitude, phase, times = self._compute_spectrogram(data)

        # Find topological components
        labeled, n_large, components = self._find_connected_components(
            magnitude,
            threshold_percentile=threshold_percentile,
            min_size=min_component_size
        )

        # Count scattered (small) components
        total_components = np.max(labeled)
        n_scattered = total_components - n_large

        # Classify components as seismic phases
        arrivals = self._classify_components(components)

        # Find S-wave for coda analysis
        s_arrival = None
        for arr in arrivals:
            if arr.phase_name == "S":
                s_arrival = arr.arrival_time
                break

        # Analyze coda
        coda_decay, _ = self._analyze_coda(data, s_arrival)

        # Dominant frequency
        avg_spectrum = np.mean(magnitude, axis=1)
        dom_idx = np.argmax(avg_spectrum)
        dominant_freq = self.freqs[dom_idx]

        # Spectral width (bandwidth)
        spectrum_norm = avg_spectrum / np.sum(avg_spectrum)
        spectral_centroid = np.sum(self.freqs * spectrum_norm)
        spectral_variance = np.sum((self.freqs - spectral_centroid)**2 * spectrum_norm)
        spectral_width = np.sqrt(spectral_variance)

        # Signal strength
        signal_db = 20 * np.log10(np.max(magnitude) + 1e-10)

        return SeismicTopologyResult(
            phase_arrivals=arrivals,
            n_coherent_components=n_large,
            n_scattered_components=n_scattered,
            coda_decay_time=coda_decay,
            dominant_frequency=dominant_freq,
            spectral_width=spectral_width,
            topology_map=labeled,
            signal_strength_db=signal_db
        )

    def compare_traces(
        self,
        trace1: np.ndarray,
        trace2: np.ndarray
    ) -> Dict[str, float]:
        """
        Compare topological structure of two traces.

        Useful for:
        - Comparing earthquake doublets
        - Matching repeating earthquakes
        - Detecting similar source mechanisms
        """
        result1 = self.analyze(trace1)
        result2 = self.analyze(trace2)

        # Compare arrival patterns
        phases1 = set(a.phase_name for a in result1.phase_arrivals)
        phases2 = set(a.phase_name for a in result2.phase_arrivals)
        phase_overlap = len(phases1 & phases2) / max(len(phases1 | phases2), 1)

        # Compare spectral properties
        freq_diff = abs(result1.dominant_frequency - result2.dominant_frequency)
        freq_similarity = 1 / (1 + freq_diff)

        # Compare coda
        coda_diff = abs(result1.coda_decay_time - result2.coda_decay_time)
        coda_similarity = 1 / (1 + coda_diff/10)

        # Overall similarity
        overall = (phase_overlap + freq_similarity + coda_similarity) / 3

        return {
            'phase_similarity': phase_overlap,
            'spectral_similarity': freq_similarity,
            'coda_similarity': coda_similarity,
            'overall_similarity': overall
        }


def demo_analyzer():
    """Demonstrate seismic topology analysis."""
    print("=" * 60)
    print("GeoVMS Seismic Topology Analyzer Demo")
    print("=" * 60)

    # Create synthetic seismogram
    sample_rate = 40.0
    duration = 300  # 5 minutes
    t = np.linspace(0, duration, int(duration * sample_rate))

    # Synthetic seismic arrivals
    # P-wave at 30s (high freq, impulsive)
    p_arrival = 30
    p_wave = np.exp(-((t - p_arrival) / 2)**2) * np.sin(2 * np.pi * 3 * t)

    # S-wave at 50s (lower freq, longer duration)
    s_arrival = 50
    s_wave = np.exp(-((t - s_arrival) / 5)**2) * np.sin(2 * np.pi * 1 * t)

    # Surface wave at 100s (low freq, dispersive)
    surf_arrival = 100
    chirp = np.linspace(0.05, 0.2, len(t))
    surf_wave = np.exp(-((t - surf_arrival) / 30)**2) * np.sin(2 * np.pi * np.cumsum(chirp) / sample_rate)

    # Coda (scattered waves)
    coda_start = 60
    coda = np.exp(-(t - coda_start) / 100) * (t > coda_start) * np.random.randn(len(t)) * 0.1

    # Noise
    noise = np.random.randn(len(t)) * 0.05

    # Combined signal
    seismogram = p_wave + s_wave * 2 + surf_wave * 1.5 + coda + noise

    print(f"\nSynthetic seismogram: {duration}s at {sample_rate} Hz")
    print(f"P-wave at {p_arrival}s, S-wave at {s_arrival}s, Surface at {surf_arrival}s")

    # Analyze
    analyzer = SeismicTopologyAnalyzer(sample_rate=sample_rate)
    result = analyzer.analyze(seismogram)

    print(f"\nTopological Analysis Results:")
    print(f"  Coherent components: {result.n_coherent_components}")
    print(f"  Scattered components: {result.n_scattered_components}")
    print(f"  Coda decay time: {result.coda_decay_time:.1f} s")
    print(f"  Dominant frequency: {result.dominant_frequency:.2f} Hz")
    print(f"  Signal strength: {result.signal_strength_db:.1f} dB")

    print(f"\nPhase Arrivals:")
    for arrival in result.phase_arrivals:
        print(f"  {arrival.phase_name}: {arrival.arrival_time:.1f}s "
              f"({arrival.frequency_band[0]:.2f}-{arrival.frequency_band[1]:.2f} Hz) "
              f"conf={arrival.confidence:.2f}")


if __name__ == "__main__":
    demo_analyzer()
