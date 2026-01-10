#!/usr/bin/env python3
"""
Universo-Mobius Setup
=====================

Installation script for OCTH signal processing suite.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
long_description = ""
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")

# Core dependencies
install_requires = [
    "numpy>=1.21.0,<2.0.0",
    "scipy>=1.7.0",
    "soundfile>=0.10.0",
]

# Optional dependencies
extras_require = {
    "api": [
        "fastapi>=0.100.0",
        "uvicorn[standard]>=0.22.0",
        "python-multipart>=0.0.6",
        "pydantic>=2.0.0",
    ],
    "visualization": [
        "plotly>=5.10.0",
        "dash>=2.10.0",
    ],
    "dev": [
        "pytest>=7.0.0",
        "pytest-cov>=4.0.0",
        "pytest-asyncio>=0.21.0",
        "black>=23.0.0",
        "isort>=5.12.0",
        "mypy>=1.0.0",
    ],
    "full": [
        # All optional dependencies
        "fastapi>=0.100.0",
        "uvicorn[standard]>=0.22.0",
        "python-multipart>=0.0.6",
        "pydantic>=2.0.0",
        "plotly>=5.10.0",
        "dash>=2.10.0",
    ],
}

setup(
    name="universo-mobius",
    version="0.1.0",
    author="Francisco Molina Burgos",
    author_email="francisco@example.com",  # Update with real email
    description="OCTH Signal Processing Suite - Audio cleaning based on gravitational wave algorithms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yatrogenesis/universo-mobius",  # Update with real URL
    license="Proprietary",

    # Package configuration
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.9",

    # Dependencies
    install_requires=install_requires,
    extras_require=extras_require,

    # Entry points
    entry_points={
        "console_scripts": [
            "vms-clean=vms_engine.cleaner:main",
            "vms-discriminator=discriminator_3d.app:run_discriminator",
            "vms-api=api.main:main",
        ],
    },

    # Classifiers
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Multimedia :: Sound/Audio :: Analysis",
        "Topic :: Multimedia :: Sound/Audio :: Speech",
    ],

    # Keywords
    keywords=[
        "audio processing",
        "noise reduction",
        "voice cleaning",
        "spectral analysis",
        "gravitational waves",
        "hexagonal",
        "OCTH",
        "signal processing",
    ],

    # Include package data
    include_package_data=True,
    zip_safe=False,
)
