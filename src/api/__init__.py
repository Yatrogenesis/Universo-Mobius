"""
VMS API - FastAPI Backend
=========================

REST API for VMS audio processing and signal discrimination.
Provides endpoints for audio cleaning, analysis, and real-time
processing.

Author: Francisco Molina Burgos
Date: January 2026
License: Proprietary - Patent Pending
"""

from .main import create_app, app
from .routes import router

__version__ = "0.1.0"
__all__ = ["create_app", "app", "router"]
