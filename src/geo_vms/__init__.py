"""
GeoVMS - Seismic Topological Analysis
======================================

Application of VMS topological signal processing to seismic data
for Earth and planetary interior modeling.

Key Insight: Seismic waves in solid media preserve topological
structure better than sound in air. The VMS "mountain detection"
approach is ideal for identifying coherent wave arrivals.

Applications:
1. Earth interior tomography
2. Lunar interior modeling (Apollo data)
3. Earthquake source characterization
4. Laser microphone signal recovery

Data Sources:
- IRIS/NSF SAGE: https://ds.iris.edu/
- Apollo PSE: https://pds-geosciences.wustl.edu/missions/apollo/
- GSN Real-time: SeedLink protocol

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

from .fetcher import SeismicFetcher, LunarFetcher
from .analyzer import SeismicTopologyAnalyzer
from .tomography import TopologicalTomography
from .laser_mic import LaserMicrophoneProcessor, LaserMicResult

__version__ = "0.1.0"
__all__ = [
    # Data fetching
    "SeismicFetcher",
    "LunarFetcher",
    # Topology analysis
    "SeismicTopologyAnalyzer",
    "TopologicalTomography",
    # Laser microphone
    "LaserMicrophoneProcessor",
    "LaserMicResult",
]
