# Predicciones OCTH del gist (2026-01-09) contra GWTC-5.0 / O4b: resultado

Fecha: 2026-11-01. Prerregistro: PREREGISTRO_LIGO_O4B.md (sha256 951c6a84..., commit c7a4c6b), con enmiendas A1-A3 y LT2b registradas antes de calcularlas.
Ceguera: GWTC-5.0 se publicó en GWOSC el 2026-05-22 y en Zenodo el 2026-05-26, después del gist (una sola versión, 2026-01-09T23:46:35Z).

| Prueba | Predicción del gist | Resultado O4b | Veredicto |
|---|---|---|---|
| LT1 cielo (sección 2.4) | >= 60 % de eventos bien localizados a < 30 grados de 4 posiciones | 8/40 = 20 %; isótropo f0 = 25 %; Z = -0.72, p = 0.82 | Falla. Se cumple el criterio de falsación 4.4 del gist (Z < 2) |
| LT3 ringdown medido por LVK (resumen) | Frecuencias de ringdown en f1 x {1, sqrt3, 2, sqrt7}, f1 = 2200/M | 0/55 ajustes con f_220 en bandas OCTH (marco fuente y detector). f_220 entre 110 y 210 Hz; el modo OCTH más alto está entre 60 y 96 Hz | Falla |
| LT2 potencia en bandas OCTH (prerregistrada) | Señal en >= 60 % de eventos | 24/35 con exceso; banda ancha 30-300 Hz 33/35; GR f_220 12/35 | No discriminante: mide la presencia de la señal CBC, no la estructura hexagonal (defecto de diseño de LT2) |
| LT2b especificidad (exploratoria) | Las razones sqrt3 : 2 : sqrt7 concentran potencia | 0/35 eventos con p < 0.05 contra 500 peines aleatorios; mediana p = 0.47; Fisher p = 0.72 | Sin evidencia de que las razones hexagonales sean especiales |

Datos: Zenodo 20348005 (mapas de cielo y tabla resumen de PE, md5 verificado), 21454847 (TGR, RD y QNMRF, md5 verificado) y strain O4b a 4 kHz de GWOSC.
Código: lt1_cielo.py, lt2_potencia.py, lt2_agregar.py, lt2b_especificidad.py, lt3_ringdown.py. Resultados: lt1_resultado.json, lt2_resultado.json, lt2b_resultado.json, lt3_resultado.json.

Limitaciones:
- El nulo de LT1 supone isotropía; la red de detectores no es isótropa.
- El GPS de GWOSC tiene resolución de 0.1 s.
- LT3 cubre solo los 5 eventos con análisis de ringdown publicado.
- El analista no estaba ciego (fue el mismo que hizo la auditoría).

## LT4 corregida (posfusión) y LT5 (verificador), 2026-11-01 20:4x

Criterios fijados antes de ver los datos (enmiendas B1 y B2).

- Verificador LT5 (63 unidades evento-detector, 2520 inyecciones por SNR): detección por exceso fuera de la fuente de 21.7 % (SNR 4), 46.5 % (6), 70.4 % (8) y 95.7 % (12). Falsas alarmas en ruido puro: 6.19 % (78/1260; binomial contra 5 %, p = 0.03).
- La regla prerregistrada exige >= 90 % de detección a algún SNR <= 12 (se cumple, 95.7 %) Y falsas alarmas <= 5 % (no se cumple, 6.19 %). Veredicto formal: LT4 NO INFORMATIVA.
- Datos de LT4 (35 eventos, posfusión [t_c + 0.05, t_c + 0.30] s): 2/35 eventos con p < 0.05 (5.7 %; binomial contra 5 %, p = 0.53; Fisher global p = 0.59). El evento más fuerte, GW250114 (SNR 77), da p = 0.13. Secundaria de peines: 10/35 (no decide; ese nulo está mal calibrado, B2).
- Lectura, no veredicto: la tasa observada es la del azar, y unas falsas alarmas algo altas facilitarían encontrar positivos, no esconderlos. Los datos no muestran ningún exceso posfusión en las frecuencias OCTH. Por la regla prerregistrada, la prueba no se declara concluyente.

## LT6 confirmatoria (O4a, detector recalibrado), 2026-11-01 21:2x

- Criterio de calidad: CUMPLE.
  - Detección a SNR 12: 91.9 % (>= 90 %). A SNR 8: 67.6 %.
  - Falsas alarmas: 4.76 % (785/16500; binomial contra 5 %, p = 0.93).
- Datos: 18 eventos de O4a (33 unidades evento-detector; 3 descartadas por calidad de datos alrededor del evento). 0/18 eventos con p < 0.05; Fisher global p = 0.92. El más fuerte, GW230814_230901 (SNR 42), da p = 0.46.
- Veredicto prerregistrado: DESCARTADO. No hay potencia posfusión en las frecuencias OCTH f1 x {1, sqrt3, 2, sqrt7}. Un hexágono con SNR >= 12 por evento se habría detectado en más del 90 % de los casos.
- Junto con O4b (cielo 20 %, ringdown LVK 0/55, razones no especiales 0/35, posfusión 2/35 con verificador marginal): la predicción hexagonal no se sostiene en ningún régimen probado.
