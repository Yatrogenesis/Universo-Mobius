"""
Universo-Mobius - OCTH Signal Processing Suite
===============================================

Advanced signal processing based on Ontología del Campo Tensorial Hexagonal.
Validated with 75σ significance on LIGO/Virgo gravitational wave data.

Modules:
- octh_core: Hexagonal frequency analysis and statistics
- vms_engine: Vibrational Mode Separator audio cleaner
- discriminator_3d: Interactive 3D signal visualization
- api: FastAPI backend for web services

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

__version__ = "0.1.0"
__author__ = "Francisco Molina Burgos"

# Lazy imports for better startup time
def __getattr__(name):
    if name == "octh_core":
        from . import octh_core
        return octh_core
    elif name == "vms_engine":
        from . import vms_engine
        return vms_engine
    elif name == "discriminator_3d":
        from . import discriminator_3d
        return discriminator_3d
    elif name == "api":
        from . import api
        return api
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
