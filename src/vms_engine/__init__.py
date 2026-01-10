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

Versions:
- V1 (VMSAudioCleaner): Manual parameter tuning
- V2 (VMSAdaptiveCleaner): Fully adaptive, interdependent parameters

Products:
- ClearVoice Pro: Consumer audio cleaning
- HexaMonitor: Industrial vibration analysis

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

# V1: Original cleaner with manual parameters
from .cleaner import VMSAudioCleaner, CleaningResult

# V2: Adaptive cleaner with interdependent parameters
from .adaptive import VMSAdaptiveCleaner, AdaptiveState

# Real-time processing
from .realtime import RealtimeProcessor

# Spectral tools
from .spectral import (
    SpectralAnalyzer,
    wiener_filter,
    harmonic_preserving_filter,
)

# V3: Topological cleaner (combines all methods)
from .topological import VMSTopologicalCleaner, TopologicalResult

__version__ = "0.3.0"
__all__ = [
    # V1: Manual parameters
    "VMSAudioCleaner",
    "CleaningResult",
    # V2: Adaptive parameters
    "VMSAdaptiveCleaner",
    "AdaptiveState",
    # V3: Topological (RECOMMENDED)
    "VMSTopologicalCleaner",
    "TopologicalResult",
    # Real-time
    "RealtimeProcessor",
    # Tools
    "SpectralAnalyzer",
    "wiener_filter",
    "harmonic_preserving_filter",
]
