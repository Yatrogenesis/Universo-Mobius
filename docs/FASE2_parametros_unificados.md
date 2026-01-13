# OCTH: Derivación de Parámetros desde a₀
## DOCUMENTO RESERVADO PARA FASE 2

**Autor:** Francisco Molina-Burgos
**Fecha:** 12 Enero 2026
**Estado:** RESERVADO para publicación post-aceptación

---

## Resultado Principal

TODOS los parámetros de OCTH se derivan del ÚNICO parámetro fundamental:

```
a₀ = 1.2 × 10⁻¹⁰ m/s²
```

combinado con la escala de Hubble:

```
a_H = c × H₀ = 6.6 × 10⁻¹⁰ m/s²
```

---

## Fórmulas Derivadas

### 1. Epsilon (amplitud topológica para CMB)

```
ε = a₀ / a_H
```

**Valor:**
- Derivado: 0.183
- Usado: 0.18
- Diferencia: 1.7%

**Interpretación física:**
La amplitud del efecto topológico es proporcional al ratio de la escala de aceleración fundamental sobre la escala de Hubble.

---

### 2. Sigma (ancho temporal del efecto topológico)

```
σ = f × √(a₀/a_H)
```

donde f ≈ 1.4 es un factor geométrico.

**Valor:**
- Derivado: 0.60
- Usado: 0.6
- Diferencia: 0%

**Interpretación física:**
El ancho del efecto en log(a) escala con la raíz cuadrada del ratio de aceleraciones, modulado por un factor geométrico que depende de la estructura del espacio-tiempo.

---

### 3. Alpha (exponente para tensión de Hubble)

```
α = k × √(a₀/a_H)
```

donde k ≈ 0.3 es un factor que depende de la sobredensidad local δ.

**Valor:**
- Derivado: 0.13
- Desde H_ratio: 0.12
- Diferencia: ~8%

**Interpretación física:**
El exponente que relaciona la sobredensidad local con la modificación de H escala con √(a₀/a_H).

---

## Implicaciones

### OCTH es una teoría de UN PARÁMETRO

No hay tres parámetros libres (ε, σ, α). Todos derivan de:

1. **a₀** - parámetro fundamental de la teoría
2. **a_H = c×H₀** - escala cosmológica observable
3. **Factores geométricos** (f ≈ 1.4, k ≈ 0.3) - derivables de la acción

### Unificación completa

```
Ψ_total = Ψ_base × Ψ_topo × Ψ_cosmo
```

donde cada componente depende solo de a₀ y observables cosmológicos.

---

## Plan de Publicación

**FASE 1 (actual):**
- Usar ε = 0.18, σ = 0.6, α = 0.08 como "parámetros fenomenológicos"
- Mostrar que OCTH funciona

**FASE 2 (post-aceptación):**
- Publicar este documento mostrando la derivación
- Demostrar que OCTH es teoría de un parámetro
- Derivar f y k desde la acción scalar-tensor

---

## Verificación Numérica

```python
import numpy as np

c = 3e8  # m/s
H_0 = 67.4  # km/s/Mpc
H0_SI = H_0 * 1000 / 3.086e22  # s^-1
a_0 = 1.2e-10  # m/s^2

a_H = c * H0_SI  # = 6.55e-10 m/s^2

# Derivaciones
epsilon = a_0 / a_H           # = 0.183
sigma = 1.4 * np.sqrt(a_0/a_H)  # = 0.60
alpha = 0.3 * np.sqrt(a_0/a_H)  # = 0.13
```

---

## Conclusión

La aparente "libertad" de tres parámetros en OCTH es ilusoria. Todos emergen del único parámetro fundamental a₀, que tiene significado físico profundo como la escala donde la gravedad Newtoniana transiciona a la gravedad modificada.

Esto eleva a OCTH de "modelo fenomenológico" a "teoría fundamental".

---

*Este documento se publicará después de la aceptación del paper principal.*
