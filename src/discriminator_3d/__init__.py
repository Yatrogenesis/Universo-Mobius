"""
3D Signal Discriminator
=======================

Interactive 3D visualization and manipulation tool for signal analysis.
Allows visual separation of signal from noise using volumetric
representation and interactive filters.

Key Features:
- 3D spectrogram visualization (time × frequency × amplitude)
- Interactive "sculpting" of signal blocks
- Real-time filter preview
- Export filtered signal

Based on the concept of treating signal as "solid blocks" that can
be carved out from surrounding noise.

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

from .visualizer import (
    Spectrogram3D,
    SignalDiscriminator,
    create_3d_spectrogram,
)

from .sculptor import (
    FilterSculptor,
    SignalBlock,
    create_filter_mask,
)

from .app import (
    create_dash_app,
    run_discriminator,
)

__version__ = "0.1.0"
__all__ = [
    "Spectrogram3D",
    "SignalDiscriminator",
    "create_3d_spectrogram",
    "FilterSculptor",
    "SignalBlock",
    "create_filter_mask",
    "create_dash_app",
    "run_discriminator",
]
