"""
VMS Engine - Vibrational Mode Separator
=======================================

High-performance audio noise reduction using algorithms derived
from gravitational wave signal processing.

Key Features:
- Extracts coherent voice modes from extremely noisy signals (SNR < 0 dB)
- Preserves natural voice characteristics (no "robotic" artifacts)
- Real-time capable with low latency (<50ms)
- Based on OCTH hexagonal frequency detection

Products:
- ClearVoice Pro: Consumer audio cleaning
- HexaMonitor: Industrial vibration analysis

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

from .cleaner import VMSAudioCleaner, CleaningResult
from .realtime import RealtimeProcessor
from .spectral import (
    SpectralAnalyzer,
    wiener_filter,
    harmonic_preserving_filter,
)

__version__ = "0.1.0"
__all__ = [
    "VMSAudioCleaner",
    "CleaningResult",
    "RealtimeProcessor",
    "SpectralAnalyzer",
    "wiener_filter",
    "harmonic_preserving_filter",
]
