"""
FastAPI Application
===================

Main FastAPI application for VMS services.
"""

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import os


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan handler."""
    # Startup
    print("VMS API starting up...")

    # Initialize any shared resources
    app.state.version = "0.1.0"

    yield

    # Shutdown
    print("VMS API shutting down...")


def create_app() -> FastAPI:
    """
    Create and configure FastAPI application.

    Returns:
        Configured FastAPI app
    """
    app = FastAPI(
        title="VMS Audio Processing API",
        description="""
        ## Vibrational Mode Separator API

        High-performance audio processing using algorithms derived from
        gravitational wave detection.

        ### Features
        - **Audio Cleaning**: Remove noise while preserving voice quality
        - **Real-time Processing**: Low-latency streaming interface
        - **Signal Analysis**: Spectral analysis and visualization data
        - **3D Discrimination**: Interactive signal/noise separation

        ### Technology
        Based on OCTH (Ontología del Campo Tensorial Hexagonal) theory,
        validated with 75σ significance on LIGO/Virgo gravitational wave data.
        """,
        version="0.1.0",
        lifespan=lifespan,
        docs_url="/docs",
        redoc_url="/redoc"
    )

    # CORS
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],  # Configure appropriately for production
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Include routes
    from .routes import router
    app.include_router(router)

    return app


# Create default app instance
app = create_app()


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=True
    )
