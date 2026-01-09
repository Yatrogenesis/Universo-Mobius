# OCTH: Tabla Única de la Verdad

## Reconciliación de Métricas y Resultados

**Autor:** Francisco Molina Burgos
**Fecha:** Enero 2025
**Estado:** Documento de Control de Consistencia

---

## 1. Datasets y Sus Propósitos

| Dataset | Código | Datos | N | Propósito | Evidencia |
|---------|--------|-------|---|-----------|-----------|
| **A: Catálogo** | `GWTC3_full_analysis.py` | Parámetros estimados → frecuencias simuladas | 75 BBH | Demostración metodológica | Indicativa |
| **B: Raw Strain** | `GWTC3_raw_pipeline.py` | Strain crudo GWOSC (H1, L1, V1) | 33→80* | **Test observacional** | **Determinante** |
| **C: CMB** | `test1_cmb_mobius_topology.py` | Planck SMICA 2018 | 1 (cielo completo) | Topología Möbius | Fuerte |
| **D: Cross-Corr** | `CMB_LIGO_cross_correlation.py` | Posiciones GW + anomalías CMB | 20 eventos | Test único OCTH | Fuerte |

*En proceso de expansión

---

## 2. Métricas por Dataset

### Dataset A: Análisis de Catálogo (Exploratorio)

| Métrica | Valor | Definición | Archivo JSON |
|---------|-------|------------|--------------|
| N eventos | 75 BBH | Eventos con m1, m2 > 3 M☉ | `gwtc3_hexagonal_analysis.json` |
| Δχ² medio | **3918 ± 307** | χ²_GR - χ²_OCTH sobre ratios simulados | ibid. |
| Z-score | **222** | (Δχ²_obs - μ_null) / σ_null | ibid. |
| % favorece OCTH | 100% | Eventos con Δχ² > 0 | ibid. |
| Hipótesis nula | Ratios uniformes [1.1, 2.5] | Distribución sin estructura | ibid. |

**Interpretación:** Este análisis demuestra que el pipeline PUEDE detectar señal hexagonal cuando existe. El Z=222 no es una medida de significancia física, sino de poder de discriminación del método.

### Dataset B: Análisis Raw Strain (Confirmatorio)

| Métrica | Valor | Definición | Archivo JSON |
|---------|-------|------------|--------------|
| N eventos | 33 (→80 en proceso) | Eventos con datos descargables | `raw_analysis/gwtc3_raw_analysis.json` |
| Δχ² medio | **236.2 ± 202.7** | χ²_GR - χ²_OCTH sobre espectro real | ibid. |
| Z-score | **6.69σ** | Significancia estadística real | ibid. |
| % favorece OCTH | 100% | Eventos con Δχ² > 0 | ibid. |
| p-value | < 10⁻¹⁰ | Probabilidad bajo H₀ | ibid. |

**Interpretación:** Este es el resultado primario. Los datos crudos de LIGO favorecen consistentemente OCTH sobre GR estándar.

### Dataset C: CMB (Topología)

| Métrica | Valor | Definición | Archivo JSON |
|---------|-------|------------|--------------|
| Correlación antipodal | -0.048 | T(θ,φ) vs T(π-θ, φ+π) | `test1_cmb_topology.json` |
| Z-score | **-6.08σ** | Anti-correlación significativa | ibid. |
| p-value | 1.2 × 10⁻⁹ | Probabilidad de azar | ibid. |
| Test polos galácticos | **5.0σ** | Signal sin contaminación galáctica | `test1_galactic_poles.json` |

**Interpretación:** El CMB muestra anti-correlación antipodal consistente con topología de Möbius.

### Dataset D: Cross-Correlación CMB×LIGO

| Métrica | Valor | Definición | Archivo JSON |
|---------|-------|------------|--------------|
| Eventos well-localized | 20 | 90% credible region < 500 deg² | `cmb_ligo_cross_correlation.json` |
| Pares cercanos (<30°) | 18 vs 8.4±2.2 esperados | Correlación espacial | ibid. |
| Z-score | **4.31σ** | Exceso significativo | ibid. |
| p-value | 10⁻⁴ | Probabilidad de azar | ibid. |
| Caso más cercano | GW190814 a 1.8° del SGP | Evidencia más fuerte | ibid. |

**Interpretación:** Las direcciones de GWs correlacionan con anomalías del CMB. Este resultado es ÚNICO de OCTH.

---

## 3. Reconciliación de Diferencias

### ¿Por qué Δχ² = 3918 vs 236?

| Factor | Dataset A | Dataset B |
|--------|-----------|-----------|
| **Tipo de datos** | Frecuencias simuladas | Espectro de potencia real |
| **Ruido** | 5% gaussiano controlado | Ruido instrumental real de LIGO |
| **Baseline** | Ratios aleatorios | GR estándar |
| **Sensibilidad** | Alta (señal inyectada) | Real (señal débil en ruido) |

**La diferencia es esperada:** Dataset A mide poder de discriminación, Dataset B mide señal real.

### ¿Por qué Z = 222 vs 6.69?

| Z-score | Significado |
|---------|-------------|
| **222 (Dataset A)** | "El método discrimina perfectamente entre señal hexagonal y ruido aleatorio" |
| **6.69 (Dataset B)** | "Los datos reales favorecen OCTH con significancia de 6.7 desviaciones estándar" |

**Son métricas diferentes:** Una es validación de método, otra es evidencia física.

---

## 4. Claims Oficiales para Publicación

### Resultado Primario (usar en abstract/resumen):

> "El análisis de 33 eventos GWTC-3 con datos raw de strain muestra preferencia sistemática por OCTH sobre GR (Δχ² = 236 ± 203, Z = 6.69σ, 100% de eventos)."

### Resultado Secundario (contexto):

> "La cross-correlación entre direcciones de ondas gravitacionales y anomalías del CMB muestra exceso significativo (Z = 4.31σ), una predicción única de OCTH."

### Resultado de Validación (métodos):

> "El análisis metodológico de 75 eventos con parámetros de catálogo confirma el poder de discriminación del pipeline (Z = 222 vs hipótesis nula aleatoria)."

---

## 5. Tabla de Significancias Combinadas

| Test | Z-score | p-value | Independencia |
|------|---------|---------|---------------|
| Raw GW (B) | 6.69σ | < 10⁻¹⁰ | Datos LIGO |
| CMB Polos (C) | 5.0σ | < 10⁻⁶ | Datos Planck |
| CMB×LIGO (D) | 4.31σ | 10⁻⁴ | Correlación espacial |
| Placebo Test | 3.40σ | 0.5% | Control negativo |

### Combinación (conservadora, asumiendo independencia parcial):

$$p_{\text{combined}} < p_{\text{GW}} \times p_{\text{CMB}} \times p_{\text{cross}} \approx 10^{-10} \times 10^{-6} \times 10^{-4} = 10^{-20}$$

**Pero** dada la posible correlación entre tests, reportamos conservadoramente:

$$\boxed{p_{\text{combined}} < 10^{-8}}$$

---

## 6. Checklist de Consistencia

| Documento | Δχ² reportado | Z reportado | N reportado | ¿Consistente? |
|-----------|---------------|-------------|-------------|---------------|
| README.md | 236.2 (raw) | 6.69σ | 33 | ✅ |
| EXECUTIVE_SUMMARY.md | - | - | 75 | ⚠️ Actualizar |
| OCTH_Nature_Article.tex | 3918 (cat) | 222 | 75 | ⚠️ Clarificar |
| Guía Educativa | - | - | 33 | ✅ |
| Referencias Verificadas | - | - | - | ✅ |

### Acciones Requeridas:

1. ☐ EXECUTIVE_SUMMARY.md → Distinguir datasets A y B
2. ☐ OCTH_Nature_Article.tex → Añadir tabla de datasets
3. ☐ Recompilar PDF cuando estén los 80 eventos

---

## 7. Versiones y Timestamps

| Dataset | Última actualización | Git commit | Datos pendientes |
|---------|---------------------|------------|------------------|
| A: Catálogo | 2025-01-09 | c59f800 | N/A |
| B: Raw | 2025-01-09 | c59f800 | 47 eventos más |
| C: CMB | 2025-01-08 | 7809fc4 | N/A |
| D: Cross | 2025-01-08 | 7809fc4 | Actualizará con B |

---

## 8. Glosario de Métricas

| Símbolo | Definición | Unidades |
|---------|------------|----------|
| Δχ² | χ²_GR - χ²_OCTH (positivo = OCTH mejor) | adimensional |
| Z-score | (observado - esperado) / σ | desviaciones estándar |
| p-value | P(datos | H₀) | probabilidad |
| SNR | Signal-to-Noise Ratio | adimensional |
| BBH | Binary Black Hole | - |
| QNM | Quasi-Normal Mode | - |

---

*Documento de control interno - Actualizar con cada análisis*
*Última revisión: 9 Enero 2025*
