# MASTER ROADMAP - Universo Möbius / OCTH Project
## Complete Development Blueprint

**Author:** Francisco Molina-Burgos
**Email:** fmolina@avermex.com
**Affiliation:** Avermex Research Division, Mérida, Yucatán, México
**Repository:** https://github.com/Yatrogenesis/Universo-Mobius
**Date:** January 2026
**Version:** 1.0

---

## TABLE OF CONTENTS

1. [Executive Summary](#1-executive-summary)
2. [Project Architecture](#2-project-architecture)
3. [Repository Structure](#3-repository-structure)
4. [Component Specifications](#4-component-specifications)
5. [Development Roadmap](#5-development-roadmap)
6. [Technical Stack](#6-technical-stack)
7. [API Design](#7-api-design)
8. [Deployment Strategy](#8-deployment-strategy)
9. [Continuation Protocol](#9-continuation-protocol)

---

## 1. EXECUTIVE SUMMARY

### What We Have (Completed)

| Component | Status | Evidence |
|-----------|--------|----------|
| OCTH Theory | ✅ Complete | Mathematical formalization, 6 axioms |
| GWTC-3 Analysis | ✅ Complete | 80/80 events (100%), 75σ significance |
| 13 Dataset Validation | ✅ Complete | Combined >>10σ |
| VMS Audio Cleaner | ✅ Working | +4.8dB improvement, 0.96 correlation |
| 3D Visualizations | ✅ Complete | GWTC-3 + Discriminator concept |

### OCTH Scientific Predictions (January 2026)

| Prediction | Status | Significance | Code |
|------------|--------|--------------|------|
| Rotation Curves | ✅ Verified | >>10σ | `OCTH_rotation_curves.py` |
| Tully-Fisher | ✅ Verified | >>10σ | `OCTH_tully_fisher.py` |
| CMB Anti-correlation | ✅ Verified | 5.2σ | `OCTH_cmb_anticorrelation.py` |
| Hubble Tension | ✅ Verified | 100% match | `OCTH_hubble_tension_v2.py` |
| Galaxy Clusters | ⚠️ Partial | ~0.7x ratio | `OCTH_galaxy_clusters.py` |
| Solar System | ⚠️ Partial | 10^-8 diff. | `OCTH_solar_system.py` |
| CMB Spectrum | ⏳ Pending | Needs CLASS | `OCTH_cmb_spectrum.py` |

**Combined statistical significance:** p < 10⁻⁸⁹

### What We're Building (This Roadmap)

1. **VMS Audio Cleaner** - Open-source audio cleaning tool
2. **OCTH Analyzer** - Scientific analysis platform
3. **3D Signal Discriminator** - Interactive signal processing tool

### Key Metrics

- **Scientific:** 75σ detection, p < 10⁻³⁰⁰
- **Technical:** Real-time processing < 50ms latency

---

## 2. PROJECT ARCHITECTURE

### High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    UNIVERSO MÖBIUS PLATFORM                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │   OCTH CORE  │  │  VMS ENGINE  │  │ 3D VISUALIZER│          │
│  │   (Theory)   │  │   (Signal)   │  │  (Interface) │          │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘          │
│         │                 │                 │                   │
│         └────────────┬────┴────────────────┘                   │
│                      │                                          │
│              ┌───────▼───────┐                                  │
│              │  UNIFIED API  │                                  │
│              └───────┬───────┘                                  │
│                      │                                          │
│    ┌─────────────────┼─────────────────┐                       │
│    │                 │                 │                        │
│    ▼                 ▼                 ▼                        │
│ ┌──────┐       ┌──────────┐     ┌───────────┐                  │
│ │ CLI  │       │ Web App  │     │   API     │                  │
│ │ Tool │       │ (React)  │     │ (FastAPI) │                  │
│ └──────┘       └──────────┘     └───────────┘                  │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Module Dependency Graph

```
OCTH_CORE (Foundation)
    │
    ├── hexagonal_math.py      # Frequency ratios, tensors
    ├── temporal_permeability.py # Ψ(z) calculations
    └── topology.py            # Möbius structure
         │
         ▼
VMS_ENGINE (Signal Processing)
    │
    ├── spectral_analysis.py   # FFT, STFT
    ├── noise_estimation.py    # Profile detection
    ├── harmonic_preservation.py # OCTH-enhanced
    └── realtime_processor.py  # Streaming
         │
         ▼
DISCRIMINATOR_3D (Visualization + Interaction)
    │
    ├── spectrogram_3d.py      # 3D surface generation
    ├── block_detection.py     # Solid signal regions
    ├── filter_sculptor.py     # Interactive mask creation
    └── export_filter.py       # Generate filter coefficients
         │
         ▼
APPLICATIONS
    │
    ├── clearvoice_pro/        # Audio cleaning product
    ├── octh_analyzer/         # Scientific tool
    └── hexamonitor/           # Industrial IoT
```

---

## 3. REPOSITORY STRUCTURE

### Current Structure (Existing)

```
Universo-Mobius/
├── README.md
├── MASTER_ROADMAP.md              # THIS FILE
├── COMPREHENSIVE_VALIDATION_REPORT.md
├── KNOWN_LIMITATIONS.md
├── REPRODUCIBILITY.md
│
├── paper/
│   ├── OCTH_Nature_Article.tex
│   ├── OCTH_Amplification_Mechanism.md
│   ├── OCTH_Tabla_Unica_Verdad.md
│   └── bibliography.bib
│
├── code/
│   ├── GWTC3_raw_pipeline.py      # GW analysis
│   ├── VMS_audio_cleaner.py       # Audio cleaner (WORKING)
│   ├── OCTH_3D_visualization.py   # 3D viz (TOY + OCTH)
│   ├── OCTH_3D_real_data.py       # 3D viz (real data)
│   ├── test_asesino_A_inyeccion_GR.py
│   └── test_asesino_B_null_sky.py
│
├── data/
│   ├── planck/
│   └── ligo/
│
├── results/
│   └── gwtc3/
│
├── figures/
│   ├── 3d_visualization/
│   └── 3d_real_data/
│
└── demo_audio/
    ├── original_voice.wav
    ├── noisy_audio.wav
    └── vms_cleaned.wav
```

### Target Structure (To Build)

```
Universo-Mobius/
├── ... (existing files)
│
├── src/                           # NEW: Main source code
│   ├── __init__.py
│   │
│   ├── octh_core/                 # OCTH mathematical foundation
│   │   ├── __init__.py
│   │   ├── constants.py           # PHI, SQRT3, SQRT7, HEX_RATIOS
│   │   ├── hexagonal.py           # Hexagonal math operations
│   │   ├── permeability.py        # Ψ(z) temporal permeability
│   │   ├── tensor_field.py        # Field calculations
│   │   └── topology.py            # Möbius topology
│   │
│   ├── vms_engine/                # Signal processing engine
│   │   ├── __init__.py
│   │   ├── audio_cleaner.py       # Main cleaner class
│   │   ├── noise_profiler.py      # Noise estimation
│   │   ├── harmonic_detector.py   # F0 + harmonics
│   │   ├── spectral_filter.py     # Wiener + OCTH enhancement
│   │   ├── realtime.py            # Streaming processor
│   │   └── io_utils.py            # Load/save audio
│   │
│   ├── discriminator_3d/          # 3D visualization + interaction
│   │   ├── __init__.py
│   │   ├── spectrogram.py         # 3D spectrogram generation
│   │   ├── block_detector.py      # Signal block identification
│   │   ├── filter_sculptor.py     # Interactive filter design
│   │   ├── visualizer.py          # Plotly/Matplotlib rendering
│   │   └── export.py              # Export filter coefficients
│   │
│   ├── analyzers/                 # Analysis tools
│   │   ├── __init__.py
│   │   ├── gwtc3_analyzer.py      # Gravitational wave analysis
│   │   ├── audio_analyzer.py      # Audio signal analysis
│   │   └── industrial_analyzer.py # Vibration/IoT analysis
│   │
│   └── api/                       # API layer
│       ├── __init__.py
│       ├── main.py                # FastAPI app
│       ├── routes/
│       │   ├── audio.py           # /api/audio/*
│       │   ├── analysis.py        # /api/analysis/*
│       │   └── visualize.py       # /api/visualize/*
│       └── schemas.py             # Pydantic models
│
├── apps/                          # Applications
│   ├── clearvoice_cli/            # CLI audio cleaner
│   │   ├── __init__.py
│   │   └── main.py
│   │
│   ├── clearvoice_web/            # Web application (React)
│   │   ├── package.json
│   │   ├── src/
│   │   │   ├── App.tsx
│   │   │   ├── components/
│   │   │   │   ├── AudioUploader.tsx
│   │   │   │   ├── Spectrogram3D.tsx
│   │   │   │   ├── FilterSculptor.tsx
│   │   │   │   └── AudioPlayer.tsx
│   │   │   └── hooks/
│   │   │       └── useAudioProcessor.ts
│   │   └── public/
│   │
│   └── discriminator_app/         # Standalone 3D discriminator
│       ├── __init__.py
│       └── main.py                # Plotly Dash app
│
├── tests/
│   ├── test_octh_core.py
│   ├── test_vms_engine.py
│   ├── test_discriminator.py
│   └── test_api.py
│
├── docker/
│   ├── Dockerfile
│   ├── docker-compose.yml
│   └── nginx.conf
│
├── scripts/
│   ├── setup.sh
│   ├── run_demo.py
│   └── benchmark.py
│
├── docs/
│   ├── API.md
│   ├── INSTALLATION.md
│   ├── USER_GUIDE.md
│   └── DEVELOPER_GUIDE.md
│
├── requirements.txt
├── setup.py
├── pyproject.toml
└── Makefile
```

---

## 4. COMPONENT SPECIFICATIONS

### 4.1 OCTH Core Module

**Purpose:** Mathematical foundation for hexagonal spacetime calculations

```python
# src/octh_core/constants.py
PHI = 1.618033988749895      # Golden ratio
SQRT3 = 1.7320508075688772   # √3
SQRT7 = 2.6457513110645907   # √7
SQRT12 = 3.4641016151377544  # √12

HEX_RATIOS = [1.0, SQRT3, 2.0, SQRT7, 3.0, SQRT12]

# Planck units
L_PLANCK = 1.616255e-35      # meters
T_PLANCK = 5.391247e-44      # seconds
```

```python
# src/octh_core/hexagonal.py
class HexagonalField:
    """Hexagonal tensor field calculations."""

    def __init__(self, resolution=100):
        self.resolution = resolution

    def compute_field(self, psi, kappa, phi):
        """
        Compute field intensity F(Ψ, κ, φ).

        Parameters:
            psi: Temporal permeability
            kappa: Hexagonal curvature
            phi: Topological phase

        Returns:
            Field intensity with hexagonal symmetry
        """
        pass

    def detect_hexagonal_pattern(self, frequencies, tolerance=0.05):
        """
        Detect if frequency ratios match hexagonal pattern.

        Returns:
            hex_score: 0-1 (1 = perfect hexagonal)
            matched_ratios: list of detected ratios
        """
        pass
```

### 4.2 VMS Engine Module

**Purpose:** Vibrational Mode Separator for signal processing

```python
# src/vms_engine/audio_cleaner.py
class VMSAudioCleaner:
    """
    Main audio cleaning class using OCTH-enhanced spectral processing.

    Key features:
    - Wiener filtering with harmonic preservation
    - OCTH ratio enhancement (√3, √7 harmonics)
    - Real-time capable (< 50ms latency)
    """

    def __init__(
        self,
        sample_rate: int = 44100,
        frame_size: int = 2048,
        hop_size: int = 512,
        noise_floor_percentile: float = 25,
        harmonic_weight: float = 0.6
    ):
        pass

    def clean_audio(self, audio, noise_sample=None) -> np.ndarray:
        """Clean complete audio signal."""
        pass

    def process_frame(self, frame) -> np.ndarray:
        """Process single frame (for real-time)."""
        pass

    def estimate_noise_profile(self, audio, duration=0.5):
        """Estimate noise from audio segment."""
        pass
```

```python
# src/vms_engine/realtime.py
class RealtimeProcessor:
    """
    Streaming audio processor for real-time applications.

    Usage:
        processor = RealtimeProcessor(sample_rate=44100)
        for chunk in audio_stream:
            cleaned = processor.process_chunk(chunk)
            output_stream.write(cleaned)
    """

    def __init__(self, sample_rate, chunk_size=1024):
        self.cleaner = VMSAudioCleaner(sample_rate=sample_rate)
        self.buffer = CircularBuffer(size=chunk_size * 4)

    def process_chunk(self, chunk) -> np.ndarray:
        """Process audio chunk with overlap."""
        pass

    def get_latency_ms(self) -> float:
        """Return current processing latency."""
        pass
```

### 4.3 3D Discriminator Module

**Purpose:** Interactive 3D signal/noise separation

```python
# src/discriminator_3d/spectrogram.py
class Spectrogram3D:
    """
    Generate and manipulate 3D spectrograms.

    Axes:
    - X: Time
    - Y: Frequency
    - Z: Amplitude (dB)
    """

    def __init__(self, sample_rate, nperseg=512):
        self.sr = sample_rate
        self.nperseg = nperseg

    def compute(self, audio) -> tuple:
        """Return (times, frequencies, amplitudes_db)."""
        pass

    def to_mesh(self) -> dict:
        """Convert to mesh format for 3D rendering."""
        pass
```

```python
# src/discriminator_3d/block_detector.py
class BlockDetector:
    """
    Detect solid signal blocks in 3D spectrogram.

    A "block" is a contiguous region of high-energy signal
    that can be extracted or masked.
    """

    def __init__(self, threshold_db=-40):
        self.threshold = threshold_db

    def detect_blocks(self, spectrogram) -> list:
        """
        Find signal blocks above threshold.

        Returns:
            List of Block objects with:
            - time_range: (t_start, t_end)
            - freq_range: (f_start, f_end)
            - amplitude_range: (a_min, a_max)
            - volume: 3D volume of block
        """
        pass

    def detect_harmonics(self, spectrogram, f0) -> list:
        """Detect harmonic pillars for given fundamental."""
        pass
```

```python
# src/discriminator_3d/filter_sculptor.py
class FilterSculptor:
    """
    Interactive 3D filter design.

    User defines regions in 3D space to create filters:
    - Include regions: signal to preserve
    - Exclude regions: noise to remove

    Outputs filter coefficients for VMSEngine.
    """

    def __init__(self, spectrogram):
        self.spec = spectrogram
        self.include_regions = []
        self.exclude_regions = []

    def add_include_region(self, t_range, f_range, a_range=None):
        """Add region to preserve."""
        pass

    def add_exclude_region(self, t_range, f_range, a_range=None):
        """Add region to remove."""
        pass

    def generate_mask_3d(self) -> np.ndarray:
        """Generate 3D mask from defined regions."""
        pass

    def export_filter(self, format='numpy') -> dict:
        """
        Export filter coefficients.

        Formats: 'numpy', 'json', 'vms' (for VMSEngine)
        """
        pass
```

```python
# src/discriminator_3d/visualizer.py
class Visualizer3D:
    """
    3D visualization using Plotly for interactivity.

    Features:
    - Rotate/zoom/pan
    - Click to select regions
    - Real-time filter preview
    """

    def __init__(self, backend='plotly'):
        self.backend = backend

    def render_spectrogram(self, spec, blocks=None) -> go.Figure:
        """Render 3D spectrogram with optional blocks highlighted."""
        pass

    def render_filter_sculptor(self, sculptor) -> go.Figure:
        """Render interactive filter sculptor interface."""
        pass

    def render_comparison(self, original, filtered) -> go.Figure:
        """Side-by-side before/after comparison."""
        pass
```

### 4.4 API Specifications

```python
# src/api/main.py
from fastapi import FastAPI, UploadFile, File
from pydantic import BaseModel

app = FastAPI(
    title="ClearVoice Pro API",
    description="OCTH-enhanced audio processing",
    version="1.0.0"
)

# Endpoints:
# POST /api/audio/clean          - Clean uploaded audio
# POST /api/audio/analyze        - Analyze audio (no cleaning)
# GET  /api/audio/{job_id}       - Get processing result
#
# POST /api/visualize/spectrogram - Generate 3D spectrogram
# POST /api/visualize/blocks      - Detect signal blocks
#
# POST /api/filter/create        - Create filter from regions
# POST /api/filter/apply         - Apply filter to audio
```

---

## 5. DEVELOPMENT ROADMAP

### Phase 1: Foundation (Week 1-2)
**Status: 80% Complete**

- [x] OCTH mathematical formalization
- [x] GWTC-3 analysis pipeline
- [x] VMS Audio Cleaner (basic)
- [x] 3D visualizations (static)
- [ ] Refactor code into `src/` structure
- [ ] Add unit tests
- [ ] Create `setup.py` / `pyproject.toml`

### Phase 2: Core Engine (Week 3-4)

- [ ] **octh_core module**
  - [ ] `constants.py` - All OCTH constants
  - [ ] `hexagonal.py` - Hexagonal field calculations
  - [ ] `permeability.py` - Ψ(z) functions
  - [ ] Unit tests

- [ ] **vms_engine module**
  - [ ] Refactor `VMS_audio_cleaner.py` into module
  - [ ] Add `realtime.py` for streaming
  - [ ] Optimize for < 50ms latency
  - [ ] Benchmark tests

### Phase 3: 3D Discriminator (Week 5-6)

- [ ] **discriminator_3d module**
  - [ ] `spectrogram.py` - 3D spectrogram generation
  - [ ] `block_detector.py` - Signal block detection
  - [ ] `filter_sculptor.py` - Interactive filter design
  - [ ] `visualizer.py` - Plotly integration
  - [ ] `export.py` - Filter export

- [ ] **Interactive Demo**
  - [ ] Plotly Dash application
  - [ ] Load audio file
  - [ ] Generate 3D spectrogram
  - [ ] Click to define regions
  - [ ] Preview filtered result
  - [ ] Export cleaned audio

### Phase 4: API & Web App (Week 7-8)

- [ ] **FastAPI Backend**
  - [ ] Audio endpoints
  - [ ] Visualization endpoints
  - [ ] Filter endpoints
  - [ ] WebSocket for real-time

- [ ] **React Frontend**
  - [ ] Audio uploader component
  - [ ] 3D spectrogram viewer (Three.js/Plotly)
  - [ ] Filter sculptor interface
  - [ ] Audio player with A/B comparison

### Phase 5: Productization (Week 9-10)

- [ ] **ClearVoice Pro CLI**
  - [ ] `clearvoice clean input.wav output.wav`
  - [ ] `clearvoice analyze input.wav`
  - [ ] `clearvoice interactive input.wav`

- [ ] **Docker Deployment**
  - [ ] Dockerfile
  - [ ] docker-compose for full stack
  - [ ] nginx configuration

- [ ] **Documentation**
  - [ ] API documentation
  - [ ] User guide
  - [ ] Developer guide

### Phase 6: Launch (Week 11-12)

- [ ] **Beta Testing**
  - [ ] 10 podcasters beta test
  - [ ] Collect feedback
  - [ ] Fix issues

- [ ] **Marketing**
  - [ ] Landing page (clearvoicepro.ai)
  - [ ] Demo videos
  - [ ] Pitch deck

- [ ] **Launch**
  - [ ] GitHub release
  - [ ] Zenodo DOI
  - [ ] Documentation complete

---

## 6. TECHNICAL STACK

### Backend

| Component | Technology | Reason |
|-----------|------------|--------|
| Language | Python 3.10+ | Scientific computing ecosystem |
| API | FastAPI | Async, fast, OpenAPI docs |
| Signal Processing | NumPy, SciPy | Industry standard |
| Audio I/O | soundfile, librosa | Format support |
| 3D Visualization | Plotly | Interactive, web-ready |
| Task Queue | Celery + Redis | Background processing |
| Database | PostgreSQL | User data, jobs |
| Cache | Redis | Session, results |

### Frontend

| Component | Technology | Reason |
|-----------|------------|--------|
| Framework | React 18 + TypeScript | Type safety, ecosystem |
| 3D Rendering | Three.js / react-three-fiber | WebGL performance |
| State | Zustand | Simple, performant |
| Styling | Tailwind CSS | Rapid development |
| Audio | Web Audio API | Browser-native |

### Infrastructure

| Component | Technology | Reason |
|-----------|------------|--------|
| Container | Docker | Reproducibility |
| Orchestration | Docker Compose / K8s | Scaling |
| CI/CD | GitHub Actions | Integration |
| Hosting | Vercel (FE) + Railway (BE) | Easy deployment |
| CDN | Cloudflare | Performance |

### Dependencies (requirements.txt)

```
# Core
numpy>=1.24.0
scipy>=1.10.0
matplotlib>=3.7.0

# Audio
soundfile>=0.12.0
librosa>=0.10.0

# API
fastapi>=0.100.0
uvicorn>=0.22.0
python-multipart>=0.0.6
pydantic>=2.0.0

# 3D Visualization
plotly>=5.15.0
dash>=2.11.0

# Database
sqlalchemy>=2.0.0
asyncpg>=0.28.0
redis>=4.6.0

# Testing
pytest>=7.4.0
pytest-asyncio>=0.21.0
httpx>=0.24.0

# Development
black>=23.0.0
ruff>=0.0.280
mypy>=1.4.0
```

---

## 7. API DESIGN

### REST Endpoints

```yaml
# Audio Processing
POST /api/v1/audio/upload:
  description: Upload audio file for processing
  request: multipart/form-data (audio file)
  response: { job_id: string, status: "queued" }

POST /api/v1/audio/clean:
  description: Clean audio with VMS engine
  request: { job_id: string, settings: CleanSettings }
  response: { job_id: string, status: "processing" }

GET /api/v1/audio/{job_id}:
  description: Get job status and result
  response: { status: string, result_url?: string, progress?: number }

GET /api/v1/audio/{job_id}/download:
  description: Download processed audio
  response: audio/wav

# 3D Visualization
POST /api/v1/visualize/spectrogram:
  description: Generate 3D spectrogram data
  request: { job_id: string, settings: SpectrogramSettings }
  response: { mesh_data: MeshData, blocks: Block[] }

# Filter Design
POST /api/v1/filter/create:
  description: Create filter from 3D regions
  request: { regions: Region[], type: "include" | "exclude" }
  response: { filter_id: string, coefficients: FilterCoeffs }

POST /api/v1/filter/apply:
  description: Apply filter to audio
  request: { job_id: string, filter_id: string }
  response: { job_id: string, status: "processing" }

# Analysis
POST /api/v1/analyze/hexagonal:
  description: Analyze signal for hexagonal patterns
  request: { job_id: string }
  response: { hex_score: number, ratios: number[], classification: string }
```

### WebSocket Events

```yaml
# Real-time processing updates
WS /api/v1/ws/{job_id}:
  events:
    - progress: { percent: number, stage: string }
    - spectrogram_update: { partial_data: MeshData }
    - complete: { result_url: string }
    - error: { message: string }
```

### Data Models

```python
# src/api/schemas.py
from pydantic import BaseModel
from typing import Optional, List

class CleanSettings(BaseModel):
    noise_reduction: float = 0.7      # 0-1
    harmonic_weight: float = 0.6      # 0-1
    use_octh_enhancement: bool = True
    preserve_dynamics: bool = True

class SpectrogramSettings(BaseModel):
    nperseg: int = 512
    freq_min: float = 0
    freq_max: float = 8000
    db_min: float = -80
    db_max: float = 0

class Region(BaseModel):
    time_start: float
    time_end: float
    freq_start: float
    freq_end: float
    amplitude_min: Optional[float] = None
    amplitude_max: Optional[float] = None

class Block(BaseModel):
    id: str
    region: Region
    volume: float
    energy: float
    classification: str  # "signal", "noise", "harmonic"

class HexagonalAnalysis(BaseModel):
    hex_score: float           # 0-1
    detected_ratios: List[float]
    expected_ratios: List[float]
    classification: str        # "HEXAGONAL" or "STANDARD"
    confidence: float
```

---

## 8. DEPLOYMENT STRATEGY

### Development Environment

```bash
# Clone and setup
git clone https://github.com/Yatrogenesis/Universo-Mobius.git
cd Universo-Mobius
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# Run tests
pytest

# Run demo
python scripts/run_demo.py

# Start API (dev)
uvicorn src.api.main:app --reload
```

### Docker Deployment

```dockerfile
# docker/Dockerfile
FROM python:3.10-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/
COPY scripts/ ./scripts/

EXPOSE 8000
CMD ["uvicorn", "src.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

```yaml
# docker/docker-compose.yml
version: '3.8'

services:
  api:
    build: .
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql://user:pass@db:5432/clearvoice
      - REDIS_URL=redis://redis:6379
    depends_on:
      - db
      - redis

  worker:
    build: .
    command: celery -A src.api.worker worker -l info
    depends_on:
      - redis

  db:
    image: postgres:15
    environment:
      - POSTGRES_USER=user
      - POSTGRES_PASSWORD=pass
      - POSTGRES_DB=clearvoice
    volumes:
      - postgres_data:/var/lib/postgresql/data

  redis:
    image: redis:7-alpine

  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx.conf:/etc/nginx/nginx.conf

volumes:
  postgres_data:
```

### CI/CD Pipeline

```yaml
# .github/workflows/ci.yml
name: CI

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install -e ".[dev]"
      - run: pytest --cov=src
      - run: ruff check src/
      - run: mypy src/

  build:
    needs: test
    runs-on: ubuntu-latest
    if: github.ref == 'refs/heads/main'
    steps:
      - uses: actions/checkout@v3
      - uses: docker/build-push-action@v4
        with:
          push: true
          tags: clearvoicepro/api:latest
```

---

## 9. CONTINUATION PROTOCOL

### If Session Ends - How to Continue

1. **Read this file first:**
   ```
   /Users/yatrogenesis/Desktop/Universo-Mobius/MASTER_ROADMAP.md
   ```

2. **Check current status:**
   ```bash
   cd /Users/yatrogenesis/Desktop/Universo-Mobius
   git status
   git log --oneline -10
   ```

3. **Key files to understand project:**
   - `MASTER_ROADMAP.md` - This roadmap (start here)
   - `COMPREHENSIVE_VALIDATION_REPORT.md` - Scientific results
   - `code/VMS_audio_cleaner.py` - Working audio cleaner
   - `code/OCTH_3D_real_data.py` - 3D visualizations

4. **Working demos to verify:**
   ```bash
   # Test audio cleaner
   python code/VMS_audio_cleaner.py

   # Test 3D visualization
   python code/OCTH_3D_real_data.py
   ```

5. **Next immediate tasks (in order):**
   - [ ] Create `src/` directory structure
   - [ ] Move code into proper modules
   - [ ] Add `setup.py` and `pyproject.toml`
   - [ ] Create interactive 3D discriminator with Plotly Dash
   - [ ] Add real-time audio processing

### Critical Context

- **Author:** Francisco Molina-Burgos (fmolina@avermex.com)
- **Affiliation:** Avermex Research Division, Mérida, Yucatán, México
- **Unique value:** Algorithms from gravitational wave detection (75σ)
- **Tech foundation:** VMS (Vibrational Mode Separator)
- **License:** AGPL-3.0

### Commands to Resume Work

```bash
# Quick start
cd /Users/yatrogenesis/Desktop/Universo-Mobius
source venv/bin/activate 2>/dev/null || python -m venv venv && source venv/bin/activate
pip install -q numpy scipy matplotlib plotly

# Run demo to verify everything works
python -c "
import sys
sys.path.insert(0, 'code')
from VMS_audio_cleaner import VMSAudioCleaner
print('VMS Audio Cleaner: OK')
"

# Check repo status
git status
```

---

## APPENDIX: Quick Reference

### OCTH Hexagonal Ratios

| Ratio | Value | Physical Meaning |
|-------|-------|------------------|
| 1 | 1.000 | Fundamental |
| √3 | 1.732 | Hexagonal first mode |
| 2 | 2.000 | First harmonic |
| √7 | 2.646 | Hexagonal second mode |
| 3 | 3.000 | Second harmonic |
| √12 | 3.464 | Hexagonal third mode |

### Key Metrics

| Metric | Value | Source |
|--------|-------|--------|
| GW events hexagonal | 80/80 (100%) | GWTC-3 analysis |
| Statistical significance | 75σ | Z-score |
| Combined p-value | < 10⁻³⁰⁰ | Fisher's method |
| Audio SNR improvement | +4.8 dB | VMS cleaner test |
| Correlation preservation | 0.96 | VMS cleaner test |

### File Locations

| What | Where |
|------|-------|
| This roadmap | `MASTER_ROADMAP.md` |
| Audio cleaner | `code/VMS_audio_cleaner.py` |
| 3D visualizations | `code/OCTH_3D_real_data.py` |
| GW analysis | `code/GWTC3_raw_pipeline.py` |
| Demo audio | `demo_audio/` |
| Figures | `figures/3d_real_data/` |

---

*Last updated: January 2026*
*Document version: 1.0*
*Status: ACTIVE DEVELOPMENT*
