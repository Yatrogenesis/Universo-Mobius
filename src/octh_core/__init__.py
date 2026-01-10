"""
OCTH Core - Ontología del Campo Tensorial Hexagonal
====================================================

Core mathematical foundations for hexagonal spacetime theory.
Provides constants, transforms, and frequency analysis primitives.

Key concepts:
- Hexagonal frequency ratios (1, sqrt(3), 2, sqrt(7), 3, sqrt(12))
- Tensor field equations for curved spacetime
- Signal resonance detection algorithms

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

from .constants import (
    HEXAGONAL_RATIOS,
    GOLDEN_RATIO,
    PLANCK_UNITS,
    TOLERANCE_SIGMA,
    CONFIDENCE_LEVELS,
)

from .transforms import (
    HexagonalTransform,
    FrequencyRatioAnalyzer,
    compute_ratio_proximity,
    find_hexagonal_modes,
)

from .statistics import (
    calculate_binomial_significance,
    calculate_combined_significance,
    monte_carlo_null,
)

__version__ = "0.1.0"
__all__ = [
    "HEXAGONAL_RATIOS",
    "GOLDEN_RATIO",
    "PLANCK_UNITS",
    "TOLERANCE_SIGMA",
    "CONFIDENCE_LEVELS",
    "HexagonalTransform",
    "FrequencyRatioAnalyzer",
    "compute_ratio_proximity",
    "find_hexagonal_modes",
    "calculate_binomial_significance",
    "calculate_combined_significance",
    "monte_carlo_null",
]
