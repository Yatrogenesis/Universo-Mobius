"""
Filter Sculptor
===============

Interactive tools for "sculpting" filter masks in 3D.
Allows users to visually define signal regions by
manipulating 3D volumes.
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass, field


@dataclass
class SignalBlock:
    """
    Represents a 3D block of signal.

    A signal block is a rectangular region in the time-frequency
    plane with an amplitude threshold.
    """
    id: int
    time_start: float       # seconds
    time_end: float         # seconds
    freq_start: float       # Hz
    freq_end: float         # Hz
    threshold_db: float     # Minimum amplitude in dB
    label: str = ""         # Optional label
    color: str = "green"    # Display color
    is_signal: bool = True  # True for signal, False for noise

    @property
    def duration(self) -> float:
        return self.time_end - self.time_start

    @property
    def bandwidth(self) -> float:
        return self.freq_end - self.freq_start

    def contains_point(
        self,
        time: float,
        freq: float,
        amplitude: float
    ) -> bool:
        """Check if a point falls within this block."""
        return (
            self.time_start <= time <= self.time_end and
            self.freq_start <= freq <= self.freq_end and
            amplitude >= self.threshold_db
        )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'id': self.id,
            'time_range': [self.time_start, self.time_end],
            'freq_range': [self.freq_start, self.freq_end],
            'threshold_db': self.threshold_db,
            'label': self.label,
            'color': self.color,
            'is_signal': self.is_signal
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'SignalBlock':
        """Create from dictionary."""
        return cls(
            id=data['id'],
            time_start=data['time_range'][0],
            time_end=data['time_range'][1],
            freq_start=data['freq_range'][0],
            freq_end=data['freq_range'][1],
            threshold_db=data['threshold_db'],
            label=data.get('label', ''),
            color=data.get('color', 'green'),
            is_signal=data.get('is_signal', True)
        )


@dataclass
class FilterSculptor:
    """
    Interactive filter mask sculptor.

    Allows creation and manipulation of signal blocks to
    define custom filter masks.
    """
    blocks: List[SignalBlock] = field(default_factory=list)
    n_freq_bins: int = 0
    n_time_bins: int = 0
    frequencies: np.ndarray = field(default_factory=lambda: np.array([]))
    times: np.ndarray = field(default_factory=lambda: np.array([]))
    _next_id: int = 1

    def initialize(
        self,
        frequencies: np.ndarray,
        times: np.ndarray
    ):
        """
        Initialize sculptor with spectrogram dimensions.

        Args:
            frequencies: Frequency array
            times: Time array
        """
        self.frequencies = frequencies
        self.times = times
        self.n_freq_bins = len(frequencies)
        self.n_time_bins = len(times)

    def add_block(
        self,
        time_range: Tuple[float, float],
        freq_range: Tuple[float, float],
        threshold_db: float = -50.0,
        label: str = "",
        is_signal: bool = True
    ) -> SignalBlock:
        """
        Add a new signal block.

        Args:
            time_range: (start, end) in seconds
            freq_range: (start, end) in Hz
            threshold_db: Amplitude threshold
            label: Optional label
            is_signal: True for signal, False for noise

        Returns:
            Created SignalBlock
        """
        block = SignalBlock(
            id=self._next_id,
            time_start=time_range[0],
            time_end=time_range[1],
            freq_start=freq_range[0],
            freq_end=freq_range[1],
            threshold_db=threshold_db,
            label=label,
            color='green' if is_signal else 'red',
            is_signal=is_signal
        )
        self.blocks.append(block)
        self._next_id += 1
        return block

    def remove_block(self, block_id: int) -> bool:
        """Remove a block by ID."""
        for i, block in enumerate(self.blocks):
            if block.id == block_id:
                self.blocks.pop(i)
                return True
        return False

    def update_block(
        self,
        block_id: int,
        **kwargs
    ) -> Optional[SignalBlock]:
        """Update block properties."""
        for block in self.blocks:
            if block.id == block_id:
                for key, value in kwargs.items():
                    if hasattr(block, key):
                        setattr(block, key, value)
                return block
        return None

    def get_block(self, block_id: int) -> Optional[SignalBlock]:
        """Get block by ID."""
        for block in self.blocks:
            if block.id == block_id:
                return block
        return None

    def generate_mask(
        self,
        spectrogram: np.ndarray = None,
        soft_edges: bool = True,
        edge_width: int = 3
    ) -> np.ndarray:
        """
        Generate filter mask from all blocks.

        Args:
            spectrogram: Optional spectrogram for amplitude thresholding
            soft_edges: Apply soft edges to transitions
            edge_width: Width of soft edge in bins

        Returns:
            2D mask array (n_freq, n_time)
        """
        mask = np.zeros((self.n_freq_bins, self.n_time_bins))

        # Process each block
        for block in self.blocks:
            if not block.is_signal:
                continue  # Noise blocks handled separately

            # Find bin indices
            freq_mask = (
                (self.frequencies >= block.freq_start) &
                (self.frequencies <= block.freq_end)
            )
            time_mask = (
                (self.times >= block.time_start) &
                (self.times <= block.time_end)
            )

            # Create block mask
            block_mask = np.outer(freq_mask.astype(float), time_mask.astype(float))

            # Apply amplitude threshold if spectrogram provided
            if spectrogram is not None:
                block_mask *= (spectrogram >= block.threshold_db)

            mask = np.maximum(mask, block_mask)

        # Apply noise blocks (subtract from mask)
        for block in self.blocks:
            if block.is_signal:
                continue

            freq_mask = (
                (self.frequencies >= block.freq_start) &
                (self.frequencies <= block.freq_end)
            )
            time_mask = (
                (self.times >= block.time_start) &
                (self.times <= block.time_end)
            )

            noise_region = np.outer(freq_mask.astype(float), time_mask.astype(float))
            mask *= (1 - noise_region)

        # Soft edges
        if soft_edges and edge_width > 0:
            from scipy.ndimage import gaussian_filter
            mask = gaussian_filter(mask, sigma=edge_width / 2)

        return np.clip(mask, 0, 1)

    def create_preset_voice(
        self,
        threshold_db: float = -40.0
    ) -> List[SignalBlock]:
        """
        Create preset blocks for voice isolation.

        Returns:
            List of created blocks
        """
        self.blocks.clear()
        self._next_id = 1

        # Fundamental voice range
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(80, 300),
            threshold_db=threshold_db,
            label="Voice Fundamental"
        )

        # Formants
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(300, 3500),
            threshold_db=threshold_db - 10,
            label="Voice Formants"
        )

        # Sibilants
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(3500, 8000),
            threshold_db=threshold_db - 20,
            label="Sibilants"
        )

        # Remove hum
        for hum_freq in [50, 60]:
            self.add_block(
                time_range=(self.times[0], self.times[-1]),
                freq_range=(hum_freq - 5, hum_freq + 5),
                threshold_db=-100,
                label=f"{hum_freq}Hz Hum",
                is_signal=False
            )

        return self.blocks

    def create_preset_music(
        self,
        threshold_db: float = -50.0
    ) -> List[SignalBlock]:
        """
        Create preset blocks for music isolation.

        Returns:
            List of created blocks
        """
        self.blocks.clear()
        self._next_id = 1

        # Bass
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(20, 250),
            threshold_db=threshold_db,
            label="Bass"
        )

        # Mid range
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(250, 4000),
            threshold_db=threshold_db - 5,
            label="Mid Range"
        )

        # High frequencies
        self.add_block(
            time_range=(self.times[0], self.times[-1]),
            freq_range=(4000, 16000),
            threshold_db=threshold_db - 15,
            label="High Frequencies"
        )

        return self.blocks

    def export_config(self) -> Dict[str, Any]:
        """Export sculptor configuration."""
        return {
            'blocks': [b.to_dict() for b in self.blocks],
            'n_freq_bins': self.n_freq_bins,
            'n_time_bins': self.n_time_bins
        }

    def import_config(self, config: Dict[str, Any]):
        """Import sculptor configuration."""
        self.blocks = [SignalBlock.from_dict(b) for b in config.get('blocks', [])]
        self._next_id = max([b.id for b in self.blocks], default=0) + 1


def create_filter_mask(
    spectrogram: np.ndarray,
    frequencies: np.ndarray,
    times: np.ndarray,
    signal_regions: List[Dict[str, Any]] = None,
    noise_regions: List[Dict[str, Any]] = None,
    threshold_db: float = -50.0,
    soft_edges: bool = True
) -> np.ndarray:
    """
    Convenience function to create a filter mask.

    Args:
        spectrogram: 2D spectrogram (n_freq, n_time)
        frequencies: Frequency array
        times: Time array
        signal_regions: List of signal region definitions
        noise_regions: List of noise region definitions
        threshold_db: Default amplitude threshold
        soft_edges: Apply soft edge transitions

    Returns:
        2D filter mask
    """
    sculptor = FilterSculptor()
    sculptor.initialize(frequencies, times)

    # Add signal regions
    if signal_regions:
        for region in signal_regions:
            sculptor.add_block(
                time_range=region.get('time_range', (times[0], times[-1])),
                freq_range=region.get('freq_range', (frequencies[0], frequencies[-1])),
                threshold_db=region.get('threshold_db', threshold_db),
                label=region.get('label', ''),
                is_signal=True
            )

    # Add noise regions
    if noise_regions:
        for region in noise_regions:
            sculptor.add_block(
                time_range=region.get('time_range', (times[0], times[-1])),
                freq_range=region.get('freq_range', (0, 0)),
                threshold_db=-100,
                label=region.get('label', ''),
                is_signal=False
            )

    return sculptor.generate_mask(spectrogram, soft_edges=soft_edges)
