"""
3D Spectrogram Visualizer
=========================

Core visualization components for 3D signal analysis.
Uses Plotly for interactive WebGL rendering.
"""

import numpy as np
from typing import Tuple, List, Optional, Dict, Any
from dataclasses import dataclass


@dataclass
class Spectrogram3D:
    """3D spectrogram data container."""
    magnitude: np.ndarray       # Shape: (n_freq, n_time)
    frequencies: np.ndarray     # Hz
    times: np.ndarray           # seconds
    sample_rate: int
    frame_size: int
    hop_size: int

    @property
    def duration(self) -> float:
        """Audio duration in seconds."""
        return self.times[-1] if len(self.times) > 0 else 0.0

    @property
    def freq_resolution(self) -> float:
        """Frequency resolution in Hz."""
        return self.sample_rate / self.frame_size

    @property
    def time_resolution(self) -> float:
        """Time resolution in seconds."""
        return self.hop_size / self.sample_rate


def create_3d_spectrogram(
    audio: np.ndarray,
    sample_rate: int = 44100,
    frame_size: int = 2048,
    hop_size: int = 512,
    db_scale: bool = True,
    db_floor: float = -80.0
) -> Spectrogram3D:
    """
    Create 3D spectrogram from audio.

    Args:
        audio: Audio signal (mono)
        sample_rate: Sample rate in Hz
        frame_size: FFT frame size
        hop_size: Samples between frames
        db_scale: Convert to dB scale
        db_floor: Minimum dB value

    Returns:
        Spectrogram3D object
    """
    # Window
    window = np.hanning(frame_size)

    # Frequencies
    frequencies = np.fft.rfftfreq(frame_size, 1/sample_rate)

    # Compute spectrogram
    n_frames = (len(audio) - frame_size) // hop_size + 1
    spectrogram = np.zeros((len(frequencies), n_frames))
    times = np.zeros(n_frames)

    for i in range(n_frames):
        start = i * hop_size
        frame = audio[start:start + frame_size]

        if len(frame) < frame_size:
            frame = np.pad(frame, (0, frame_size - len(frame)))

        windowed = frame * window
        spectrum = np.abs(np.fft.rfft(windowed))
        spectrogram[:, i] = spectrum
        times[i] = start / sample_rate

    # Convert to dB if requested
    if db_scale:
        spectrogram = 20 * np.log10(np.maximum(spectrogram, 1e-10))
        spectrogram = np.maximum(spectrogram, db_floor)

    return Spectrogram3D(
        magnitude=spectrogram,
        frequencies=frequencies,
        times=times,
        sample_rate=sample_rate,
        frame_size=frame_size,
        hop_size=hop_size
    )


class SignalDiscriminator:
    """
    Interactive 3D signal discriminator.

    Provides tools for visualizing and separating signal from noise
    using volumetric 3D representation.
    """

    def __init__(self, spectrogram: Spectrogram3D):
        """
        Initialize discriminator.

        Args:
            spectrogram: Input spectrogram data
        """
        self.spectrogram = spectrogram
        self.signal_mask = np.ones_like(spectrogram.magnitude, dtype=bool)
        self.noise_regions: List[Dict[str, Any]] = []
        self.signal_blocks: List[Dict[str, Any]] = []

    def detect_signal_blocks(
        self,
        threshold_db: float = -40.0,
        min_freq: float = 80.0,
        max_freq: float = 8000.0,
        min_duration: float = 0.05
    ) -> List[Dict[str, Any]]:
        """
        Automatically detect signal blocks (regions above threshold).

        Args:
            threshold_db: Threshold in dB for signal detection
            min_freq: Minimum frequency to consider
            max_freq: Maximum frequency to consider
            min_duration: Minimum block duration in seconds

        Returns:
            List of signal block definitions
        """
        from scipy import ndimage

        spec = self.spectrogram

        # Create binary mask
        freq_mask = (spec.frequencies >= min_freq) & (spec.frequencies <= max_freq)
        signal_binary = spec.magnitude > threshold_db

        # Apply frequency mask
        signal_binary[~freq_mask, :] = False

        # Label connected components
        labeled, n_features = ndimage.label(signal_binary)

        self.signal_blocks = []

        for i in range(1, n_features + 1):
            # Find bounding box
            positions = np.where(labeled == i)

            if len(positions[0]) == 0:
                continue

            freq_min_idx = positions[0].min()
            freq_max_idx = positions[0].max()
            time_min_idx = positions[1].min()
            time_max_idx = positions[1].max()

            # Check duration
            duration = (time_max_idx - time_min_idx) * spec.time_resolution
            if duration < min_duration:
                continue

            block = {
                'id': i,
                'freq_range': (
                    spec.frequencies[freq_min_idx],
                    spec.frequencies[min(freq_max_idx, len(spec.frequencies)-1)]
                ),
                'time_range': (
                    spec.times[time_min_idx],
                    spec.times[min(time_max_idx, len(spec.times)-1)]
                ),
                'peak_amplitude': spec.magnitude[positions].max(),
                'mean_amplitude': spec.magnitude[positions].mean(),
                'size': len(positions[0]),
                'duration': duration,
                'mask_indices': (positions[0], positions[1])
            }
            self.signal_blocks.append(block)

        return self.signal_blocks

    def detect_noise_regions(
        self,
        threshold_db: float = -60.0,
        hum_frequencies: List[float] = [50, 60, 100, 120, 180]
    ) -> List[Dict[str, Any]]:
        """
        Detect noise regions (low energy + known noise frequencies).

        Args:
            threshold_db: Threshold for noise floor
            hum_frequencies: Known hum frequencies to flag

        Returns:
            List of noise region definitions
        """
        spec = self.spectrogram

        self.noise_regions = []

        # Find low-energy regions
        noise_binary = spec.magnitude < threshold_db

        # Also flag hum frequencies
        for hum_freq in hum_frequencies:
            # Find closest frequency bin
            idx = np.argmin(np.abs(spec.frequencies - hum_freq))

            # Check if this frequency is consistently present
            hum_power = spec.magnitude[idx, :].mean()

            if hum_power > threshold_db:  # Hum is present
                region = {
                    'type': 'hum',
                    'frequency': hum_freq,
                    'freq_idx': idx,
                    'mean_power': hum_power,
                    'description': f'{hum_freq}Hz hum'
                }
                self.noise_regions.append(region)

        # Find broadband noise floor
        noise_floor = np.percentile(spec.magnitude, 10, axis=1)
        noise_region = {
            'type': 'floor',
            'spectrum': noise_floor,
            'mean_level': noise_floor.mean(),
            'description': 'Noise floor'
        }
        self.noise_regions.append(noise_region)

        return self.noise_regions

    def create_filter_mask(
        self,
        include_blocks: List[int] = None,
        exclude_hum: bool = True,
        noise_gate_db: float = -50.0
    ) -> np.ndarray:
        """
        Create filter mask based on detected regions.

        Args:
            include_blocks: Block IDs to include (None = all)
            exclude_hum: Remove hum frequencies
            noise_gate_db: Gate threshold for noise floor

        Returns:
            2D mask array (same shape as spectrogram)
        """
        spec = self.spectrogram
        mask = np.ones_like(spec.magnitude)

        # Start with noise gate
        mask[spec.magnitude < noise_gate_db] = 0.0

        # Remove hum
        if exclude_hum:
            for region in self.noise_regions:
                if region.get('type') == 'hum':
                    idx = region['freq_idx']
                    # Notch out hum and harmonics
                    for mult in [1, 2, 3]:
                        notch_idx = idx * mult
                        if notch_idx < len(spec.frequencies):
                            # Soft notch
                            width = 3
                            for i in range(max(0, notch_idx-width), min(len(spec.frequencies), notch_idx+width+1)):
                                mask[i, :] *= 0.1

        # If specific blocks selected, only keep those
        if include_blocks is not None:
            block_mask = np.zeros_like(mask)
            for block in self.signal_blocks:
                if block['id'] in include_blocks:
                    idx = block['mask_indices']
                    block_mask[idx] = 1.0
            mask *= block_mask

        self.signal_mask = mask > 0.5
        return mask

    def create_plotly_figure(
        self,
        show_blocks: bool = True,
        show_noise: bool = True,
        colorscale: str = 'Viridis',
        opacity: float = 0.7
    ) -> 'plotly.graph_objects.Figure':
        """
        Create interactive Plotly 3D figure.

        Args:
            show_blocks: Highlight detected signal blocks
            show_noise: Highlight noise regions
            colorscale: Plotly colorscale name
            opacity: Surface opacity

        Returns:
            Plotly Figure object
        """
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        spec = self.spectrogram

        # Create meshgrid for surface
        T, F = np.meshgrid(spec.times, spec.frequencies)

        # Main figure with subplots
        fig = make_subplots(
            rows=2, cols=2,
            specs=[
                [{'type': 'surface', 'rowspan': 2}, {'type': 'xy'}],
                [None, {'type': 'xy'}]
            ],
            subplot_titles=(
                '3D Spectrogram',
                'Time-Frequency View',
                'Frequency Spectrum'
            ),
            column_widths=[0.6, 0.4]
        )

        # 3D Surface
        fig.add_trace(
            go.Surface(
                x=T,
                y=F,
                z=spec.magnitude,
                colorscale=colorscale,
                opacity=opacity,
                name='Spectrogram',
                showscale=True,
                colorbar=dict(title='dB', x=0.45)
            ),
            row=1, col=1
        )

        # Add signal block markers
        if show_blocks and self.signal_blocks:
            for block in self.signal_blocks[:10]:  # Limit for performance
                t_range = block['time_range']
                f_range = block['freq_range']

                # Block boundary
                t_center = (t_range[0] + t_range[1]) / 2
                f_center = (f_range[0] + f_range[1]) / 2
                z_height = block['peak_amplitude']

                fig.add_trace(
                    go.Scatter3d(
                        x=[t_center],
                        y=[f_center],
                        z=[z_height + 5],
                        mode='markers+text',
                        marker=dict(size=8, color='lime', symbol='diamond'),
                        text=[f'Block {block["id"]}'],
                        textposition='top center',
                        name=f'Signal Block {block["id"]}',
                        showlegend=True
                    ),
                    row=1, col=1
                )

        # 2D heatmap view
        fig.add_trace(
            go.Heatmap(
                x=spec.times,
                y=spec.frequencies,
                z=spec.magnitude,
                colorscale=colorscale,
                showscale=False,
                name='2D View'
            ),
            row=1, col=2
        )

        # Average spectrum
        avg_spectrum = spec.magnitude.mean(axis=1)
        fig.add_trace(
            go.Scatter(
                x=spec.frequencies,
                y=avg_spectrum,
                mode='lines',
                name='Avg Spectrum',
                line=dict(color='cyan')
            ),
            row=2, col=2
        )

        # Mark noise floor
        if show_noise:
            for region in self.noise_regions:
                if region.get('type') == 'floor':
                    floor = region['spectrum']
                    fig.add_trace(
                        go.Scatter(
                            x=spec.frequencies,
                            y=floor,
                            mode='lines',
                            name='Noise Floor',
                            line=dict(color='red', dash='dash')
                        ),
                        row=2, col=2
                    )

        # Layout
        fig.update_layout(
            title='VMS Signal Discriminator',
            scene=dict(
                xaxis_title='Time (s)',
                yaxis_title='Frequency (Hz)',
                zaxis_title='Amplitude (dB)',
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=0.8)
                )
            ),
            height=800,
            showlegend=True
        )

        # Update 2D axes
        fig.update_xaxes(title_text='Time (s)', row=1, col=2)
        fig.update_yaxes(title_text='Frequency (Hz)', row=1, col=2)
        fig.update_xaxes(title_text='Frequency (Hz)', row=2, col=2)
        fig.update_yaxes(title_text='Amplitude (dB)', row=2, col=2)

        return fig

    def export_filtered_spectrogram(
        self,
        mask: np.ndarray = None
    ) -> np.ndarray:
        """
        Export spectrogram with filter mask applied.

        Args:
            mask: Filter mask (uses self.signal_mask if None)

        Returns:
            Filtered spectrogram
        """
        if mask is None:
            mask = self.signal_mask.astype(float)

        return self.spectrogram.magnitude * mask

    def reconstruct_audio(
        self,
        mask: np.ndarray = None,
        original_audio: np.ndarray = None
    ) -> np.ndarray:
        """
        Reconstruct audio from filtered spectrogram.

        Note: This is a simplified reconstruction. For best quality,
        use the VMS engine with the mask as a guide.

        Args:
            mask: Filter mask
            original_audio: Original audio for phase information

        Returns:
            Reconstructed audio
        """
        if mask is None:
            mask = self.signal_mask.astype(float)

        spec = self.spectrogram

        # Need original audio for phase
        if original_audio is None:
            raise ValueError("Original audio required for reconstruction")

        # Compute original STFT for phase
        window = np.hanning(spec.frame_size)
        n_frames = spec.magnitude.shape[1]

        # Reconstruct using Griffin-Lim with magnitude mask
        filtered_mag = spec.magnitude * mask

        # Convert back from dB
        filtered_linear = 10 ** (filtered_mag / 20)

        # Get phase from original
        output = np.zeros(len(original_audio))
        window_sum = np.zeros(len(original_audio))

        for i in range(n_frames):
            start = i * spec.hop_size
            end = start + spec.frame_size

            if end > len(original_audio):
                break

            # Get original phase
            frame = original_audio[start:end] * window
            orig_fft = np.fft.rfft(frame)
            phase = np.angle(orig_fft)

            # Reconstruct with filtered magnitude
            filtered_fft = filtered_linear[:, i] * np.exp(1j * phase)
            filtered_frame = np.fft.irfft(filtered_fft, n=spec.frame_size)

            # Overlap-add
            output[start:end] += filtered_frame * window
            window_sum[start:end] += window ** 2

        # Normalize
        window_sum = np.maximum(window_sum, 1e-8)
        output /= window_sum

        return output
