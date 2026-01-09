# Universo Möbius - OCTH

## Ontología del Campo Tensorial Hexagonal

[![Status](https://img.shields.io/badge/Status-Active_Development-brightgreen)]()
[![Test1](https://img.shields.io/badge/Test%201-Pipeline_Ready-blue)]()
[![Test2](https://img.shields.io/badge/Test%202-EXACT_EQUIVALENCE-success)]()

## Resumen

Este repositorio contiene la implementación y verificación experimental de la **OCTH** (Ontología del Campo Tensorial Hexagonal), un modelo que propone:

1. El espaciotiempo tiene estructura de **retículo hexagonal** a escala de Planck
2. La "gravedad" emerge de un campo de **permeabilidad temporal** Ψ
3. La topología del universo es una **Cinta de Möbius** 3D

## Formulación Matemática

### Campo de Permeabilidad Ψ

$$\Psi(r) = \sqrt{1 - \frac{\rho_m}{\rho_{\text{Planck}}}} = \sqrt{1 - \frac{r_s}{r}}$$

donde $r_s = 2GM/c^2$ es el radio de Schwarzschild.

### Métrica Elástico-Temporal

$$ds^2 = -c^2 \Psi^2 dt^2 + g_{ij} dx^i dx^j$$

Esto reproduce la métrica de Schwarzschild en el límite de campo débil.

### Ecuación de Onda Hexagonal

$$\Psi^2 \frac{\partial^2 u}{\partial t^2} = c^2 \Delta_{\Lambda_H} u$$

donde $\Delta_{\Lambda_H}$ es el Laplaciano discreto sobre el grafo hexagonal.

## Tests Implementados

### Test #1: Topología de Möbius en el CMB ⏳

**Objetivo:** Detectar correlación antipodal con inversión de paridad en el CMB.

**Predicción OCTH:**
```
Si el universo es una Banda de Möbius 3D:
    T(θ, φ) ~ T(π-θ, φ+π) con flip de paridad
```

**Estado: PIPELINE VALIDADO**

| Condición | Correlación | Z-score | P-value | Resultado |
|-----------|-------------|---------|---------|-----------|
| Baseline (sin señal) | 0.003 | 0.74 | 0.46 | NO SIGNIFICATIVO |
| Möbius inyectado | 0.20 | **29.09** | **< 0.0001** | **DETECTADO** |

- ✓ Pipeline de detección validado
- ✓ Test de sensibilidad funciona correctamente
- ⏳ **Pendiente**: Ejecutar con datos reales de Planck

**Figuras:**
- `fig7_cmb_mobius_analysis.png` - Análisis completo del CMB

---

### Test #2: Geodésicas Exactas en Métrica OCTH ✓✓

**Objetivo:** Demostrar que la métrica OCTH produce EXACTAMENTE las mismas geodésicas nulas que Schwarzschild.

**Resultado: EQUIVALENCIA MATEMÁTICA DEMOSTRADA**

La ecuación de órbita derivada de la métrica OCTH:
```
(du/dφ)² = 1/b² - u² + rs·u³
```
es **IDÉNTICA** a la ecuación de Schwarzschild. Esto NO es una aproximación.

| b/rs | Δφ campo débil | Δφ exacto | Ratio |
|------|----------------|-----------|-------|
| 5    | 22.92°         | 33.83°    | 1.476 |
| 10   | 11.46°         | 13.53°    | 1.181 |
| 50   | 2.29°          | 2.36°     | 1.031 |
| 100  | 1.15°          | 1.16°     | 1.015 |

- ✓ **Convergencia**: Ratio → 1.0 cuando b >> rs
- ✓ **Correcciones campo fuerte**: Hasta 47.6% para b/rs = 5
- ✓ **SIN CALIBRACIÓN**: El resultado es exacto desde primeros principios

**Figuras:**
- `fig5_exact_deflection.png` - Deflexión vs parámetro de impacto
- `fig6_trajectories.png` - Trayectorias de luz

### Test #2 (preliminar): Simulación en Retículo Hexagonal

**Objetivo:** Verificar comportamiento cualitativo en malla discreta.

**Resultado:**
- ✓ La luz se deflecta hacia la masa (Ψ bajo)
- ✓ La deflexión escala inversamente con parámetro de impacto (~1/b)
- ⚠️ Error numérico ~75% (discretización de malla, no error teórico)

**Figuras:**
- `fig1_lattices.png` - Comparación hexagonal vs cuadrado
- `fig2_psi_field.png` - Mapa de calor de Ψ
- `fig3_ray_tracing.png` - Trayectorias de luz
- `fig4_comparison.png` - OCTH vs GR

## Estructura del Repositorio

```
Universo-Mobius/
├── README.md
├── code/
│   ├── test1_cmb_mobius_topology.py        # CMB MÖBIUS (Test #1)
│   ├── test2_geodesic_exact.py             # GEODÉSICAS EXACTAS (Test #2)
│   ├── test2_hexagonal_propagation.py      # Simulación retículo (preliminar)
│   └── test2_hexagonal_propagation_v2.py   # Versión corregida
├── figures/
│   ├── fig1_lattices.png/pdf               # Comparación topologías
│   ├── fig2_psi_field.png/pdf              # Campo Ψ
│   ├── fig3_ray_tracing.png/pdf            # Trayectorias discretas
│   ├── fig4_comparison.png/pdf             # OCTH vs GR (discreto)
│   ├── fig5_exact_deflection.png/pdf       # Deflexión exacta vs campo débil
│   ├── fig6_trajectories.png/pdf           # Trayectorias exactas
│   └── fig7_cmb_mobius_analysis.png/pdf    # Análisis CMB Möbius
├── results/
│   ├── test1_cmb_topology.json             # Resultados CMB (Test #1)
│   ├── test2_results.json                  # Resultados discretización
│   └── test2_exact_geodesics.json          # EQUIVALENCIA MATEMÁTICA
├── data/
│   └── planck/                             # Datos CMB de Planck
└── paper/
```

## Correspondencia Física

```
TOPOLOGÍA              MATEMÁTICA                FÍSICA
    │                      │                        │
    ▼                      ▼                        ▼
Retículo           Ψ = √(1 - rs/r)          Dilatación temporal
Hexagonal               │                         │
    │                   ▼                         ▼
Simetría D3        ∇Ψ → deflexión          Lensing gravitacional
    │                   │                         │
    ▼                   ▼                         ▼
Nodos triádicos    n_eff = 1/Ψ             Índice de refracción
(e₁+e₂+e₃=0)                               del espaciotiempo
```

## Estado de Tests

- [x] **Test #1:** Topología Möbius en CMB → Pipeline validado, pendiente datos reales
- [x] **Test #2:** Geodésicas OCTH → **EQUIVALENCIA MATEMÁTICA DEMOSTRADA**
- [x] **Test #3:** Geometría hexagonal en galaxias → Ratio hex/sq = 1.43 (pendiente SDSS real)
- [ ] **Test #4:** Verificación Ψ en ondas gravitacionales (LIGO)
- [x] **Test #5:** VSL en GRBs → Lag INTRÍNSECO detectado (Lorentz invariance OK)

### Test #3: Geometría Hexagonal en Galaxias

**Predicción:** Si la malla fue estirada por inflación, ω(60°) y ω(120°) > ω(90°).

**Resultado (simulación):**
| Ángulo | ω(θ) | Tipo |
|--------|------|------|
| 60° | 0.236 | Hexagonal |
| 90° | 0.225 | Cuadrado |
| 120° | **0.405** | Hexagonal |

- Ratio hexagonal/cuadrado = **1.43** (43% más correlación)
- Z = 1.85, p = 0.064 (marginalmente significativo)
- ⏳ **Pendiente**: Ejecutar con datos SDSS reales

---

### Test #5: VSL en GRBs

**Predicción:** Si c_eff = c·Ψ(E), fotones de mayor energía llegarían más tarde.

**Resultado:**
| Correlación | ρ | p-value | Interpretación |
|-------------|---|---------|----------------|
| Lag vs E | **-0.747** | 0.0002 | NEGATIVA (opuesto a VSL) |
| Lag vs z | +0.467 | 0.038 | Posible evolución cósmica |

- ✓ Lag altamente significativo (Z = 6.98σ)
- ✓ **Correlación NEGATIVA**: fotones de alta E llegan ANTES
- ✓ Conclusión: Lag es INTRÍNSECO a la fuente (física del GRB)
- ✓ Límite: E_QG > 0.68 × E_Planck (consistente con Lorentz invariance)

## Dependencias

```bash
pip install numpy scipy matplotlib healpy astropy
```

## Ejecución

```bash
cd code

# Test #1: CMB Möbius (simulación)
python3 test1_cmb_mobius_topology.py

# Test #1: CMB Möbius (con datos reales de Planck, si disponibles)
python3 test1_cmb_mobius_topology.py --real

# Test #1: Validación con señal inyectada
python3 test1_cmb_mobius_topology.py --inject

# Test #2: Geodésicas exactas
python3 test2_geodesic_exact.py
```

## Autor

**Francisco Molina Burgos**
Email: yatrogenesis@proton.me

## Acknowledgments

Asistencia computacional de Claude (Anthropic).

---

**φ > 0**
