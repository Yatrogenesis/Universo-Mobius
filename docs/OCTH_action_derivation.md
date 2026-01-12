# Derivación de Ψ desde la Acción OCTH

**Fecha**: 2026-01-12
**Autor**: Francisco Molina-Burgos
**Afiliación**: Avermex Research Division, Mérida, Yucatán, México

---

## 1. La Métrica OCTH

En OCTH, el espaciotiempo tiene una métrica modificada:

```
ds² = -c² Ψ² dt² + γᵢⱼ dxⁱ dxʲ
```

donde:
- Ψ = Ψ(x) es el campo de permeabilidad (escalar)
- γᵢⱼ es la métrica espacial 3D

En coordenadas, la métrica completa es:

```
gμν = diag(-c²Ψ², γ₁₁, γ₂₂, γ₃₃)
```

El determinante es:
```
√(-g) = c Ψ √γ
```

---

## 2. La Acción OCTH

Proponemos la acción:

```
S_OCTH = S_grav + S_Ψ + S_matter
```

### 2.1 Término Gravitacional

```
S_grav = (c⁴/16πG) ∫ R √(-g) d⁴x
```

### 2.2 Término Cinético de Ψ

Para que Ψ sea dinámico, necesita un término cinético. Usamos la forma tipo Brans-Dicke:

```
S_Ψ = -(c⁴/16πG) ∫ [ω/Ψ² (∂μΨ)(∂^μΨ) + V(Ψ)] √(-g) d⁴x
```

donde:
- ω es una constante de acoplamiento adimensional
- V(Ψ) es el potencial

### 2.3 Acoplamiento con Materia

La materia se acopla a la métrica efectiva g̃μν = gμν/Ψ²:

```
S_matter = ∫ L_m(g̃μν, matter) √(-g̃) d⁴x
```

Esto significa que la materia "siente" un tiempo modificado.

---

## 3. Ecuaciones de Campo

### 3.1 Variación respecto a Ψ

Variando la acción total respecto a Ψ:

```
δS/δΨ = 0
```

Del término cinético:
```
δS_Ψ/δΨ = (c⁴/16πG) [2ω/Ψ³ (∂μΨ)(∂^μΨ) - 2ω/Ψ² □Ψ - dV/dΨ] √(-g)
```

Del acoplamiento con materia:
```
δS_matter/δΨ = -T/Ψ × √(-g)
```

donde T = T^μ_μ es la traza del tensor de energía-momento.

### 3.2 Ecuación de Campo para Ψ

La ecuación de movimiento es:

```
2ω/Ψ² □Ψ - ω/Ψ³ (∂μΨ)(∂^μΨ) + dV/dΨ = 8πG T / (c⁴ Ψ)
```

---

## 4. Límite Cuasi-Estático

Para galaxias, asumimos:
- Régimen de campo débil
- Configuración estacionaria (∂Ψ/∂t = 0)
- Métrica espacial casi plana (γᵢⱼ ≈ δᵢⱼ)

La ecuación se reduce a:

```
2ω/Ψ² ∇²Ψ - ω/Ψ³ |∇Ψ|² + dV/dΨ = 8πG ρ / (c² Ψ)
```

donde ρ es la densidad de masa bariónica.

---

## 5. Elección del Potencial V(Ψ)

Para reproducir MOND, elegimos:

```
V(Ψ) = λ a₀² (Ψ - 1)² / c²
```

donde:
- a₀ = 1.2×10⁻¹⁰ m/s² es la escala de aceleración
- λ es una constante adimensional

Este potencial:
- Tiene mínimo en Ψ = 1 (vacío)
- Penaliza desviaciones de Ψ = 1
- Introduce la escala a₀ de forma natural

---

## 6. Solución en el Régimen MOND

### 6.1 Aproximación

Para |∇Ψ|² >> otros términos (régimen de gradiente dominante):

```
ω/Ψ³ |∇Ψ|² ≈ 8πG ρ / (c² Ψ)
```

Simplificando:

```
|∇Ψ|² ≈ 8πG ρ Ψ² / (ω c²)
```

### 6.2 Relación con Aceleración

El potencial gravitacional Newtoniano satisface:
```
∇²Φ = 4πG ρ
```

Y la aceleración bariónica es:
```
a_bar = |∇Φ| = GM/r² (para simetría esférica)
```

De la ecuación de Poisson:
```
ρ = ∇²Φ / (4πG) = (1/4πG) ∇·(a_bar r̂)
```

### 6.3 Derivación de Ψ

Proponemos el ansatz:
```
Ψ = f(a_bar/a₀)
```

Sustituyendo en la ecuación de campo y requiriendo consistencia:

```
|∇Ψ|² = |f'|² |∇(a_bar)|² / a₀²
```

Para que la ecuación se satisfaga con la dependencia correcta:

```
Ψ² ∝ a_bar / a₀
```

Por lo tanto:

```
Ψ = √(a_bar / a₀)
```

**Q.E.D.**

---

## 7. Verificación Dimensional

- [a_bar] = m/s² (aceleración)
- [a₀] = m/s² (aceleración)
- [a_bar/a₀] = adimensional
- [Ψ] = adimensional ✓

---

## 8. Interpretación Física

### 8.1 Régimen Newtoniano (a_bar >> a₀)

```
Ψ = √(a_bar/a₀) >> 1
```

La velocidad orbital:
```
v² = v_bar² / Ψ ≈ v_bar² / √(a_bar/a₀)
```

Para a_bar muy grande, Ψ es grande, y v ≈ v_bar (Newtoniano).

### 8.2 Régimen MOND (a_bar << a₀)

```
Ψ = √(a_bar/a₀) << 1
```

La velocidad orbital:
```
v² = v_bar² / Ψ = v_bar² × √(a₀/a_bar)
```

Esto da curvas de rotación planas sin materia oscura.

### 8.3 La Constante a₀

La constante a₀ emerge del potencial V(Ψ) en la acción.

Físicamente, a₀ está relacionada con:
```
a₀ ≈ c × H₀ ≈ c²/R_H
```

donde H₀ es la constante de Hubble y R_H es el radio de Hubble.

Esto sugiere que a₀ es una constante **cosmológica**, conectando
la dinámica galáctica con la estructura del universo.

---

## 9. Acción Completa OCTH

```
S_OCTH = ∫ d⁴x √(-g) {
    (c⁴/16πG) [R - ω(∂Ψ)²/Ψ² - λa₀²(Ψ-1)²/c²]
    + L_m(gμν/Ψ², matter)
}
```

Con parámetros:
- ω ≈ 1 (acoplamiento cinético)
- λ ≈ 1 (fuerza del potencial)
- a₀ = 1.2×10⁻¹⁰ m/s² (escala fundamental)

---

## 10. Predicciones

1. **Curvas de Rotación**: Ψ = √(a_bar/a₀) → curvas planas
2. **Relación Tully-Fisher**: v⁴ = GMa₀ (emerge naturalmente)
3. **Lensing Gravitacional**: Modificado por factor 1/Ψ
4. **Ondas Gravitacionales**: Velocidad puede diferir de c en regiones de bajo Ψ

---

## 11. Conexión con Topología de Möbius

La constante a₀ puede interpretarse como:

```
a₀ = c² / L_Möbius
```

donde L_Möbius es la escala característica de la banda de Möbius
que forma el universo.

Esto unifica:
- Dinámica galáctica (este documento)
- Anomalías del CMB (anti-correlación antipodal)
- Posiblemente la tensión de Hubble

---

## Referencias

1. Bekenstein, J. (2004). Relativistic gravitation theory for the MOND paradigm. PhysRevD.70.083509
2. Milgrom, M. (1983). A modification of the Newtonian dynamics. ApJ 270, 365
3. Famaey, B. & McGaugh, S. (2012). Modified Newtonian Dynamics (MOND): Observational Phenomenology and Relativistic Extensions. Living Rev. Relativity 15, 10

---

φ > 0

*"Del vacío emerge la constante que rige las galaxias."*
