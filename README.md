# Universo Möbius - OCTH

## Ontología del Campo Tensorial Hexagonal

[![Status](https://img.shields.io/badge/Status-In_Development-yellow)]()
[![Test2](https://img.shields.io/badge/Test%202-Qualitative_✓-green)]()

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

### Test #2: Propagación en Retículo Hexagonal ✓

**Objetivo:** Demostrar que ondas en malla hexagonal con Ψ variable reproducen lensing gravitacional.

**Resultado:**
- ✓ La luz se deflecta hacia la masa (Ψ bajo)
- ✓ La deflexión escala inversamente con parámetro de impacto (~1/b)
- ⚠️ Error cuantitativo ~75% vs GR (requiere calibración)

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
│   ├── test2_hexagonal_propagation.py      # Test inicial
│   └── test2_hexagonal_propagation_v2.py   # Versión corregida
├── figures/
│   ├── fig1_lattices.png/pdf
│   ├── fig2_psi_field.png/pdf
│   ├── fig3_ray_tracing.png/pdf
│   └── fig4_comparison.png/pdf
├── data/
├── results/
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

## Próximos Tests

- [ ] **Test #1:** Topología Möbius en CMB (Planck data)
- [x] **Test #2:** Propagación hexagonal (COMPLETADO - cualitativo)
- [ ] **Test #3:** Geometría hexagonal en distribución de galaxias (SDSS)
- [ ] **Test #4:** Verificación Ψ en ondas gravitacionales (LIGO)
- [ ] **Test #5:** VSL en GRBs (Fermi-LAT)

## Dependencias

```bash
pip install numpy scipy matplotlib
```

## Ejecución

```bash
cd code
python3 test2_hexagonal_propagation_v2.py
```

## Autor

**Francisco Molina Burgos**
Email: yatrogenesis@proton.me

## Acknowledgments

Asistencia computacional de Claude (Anthropic).

---

**φ > 0**
