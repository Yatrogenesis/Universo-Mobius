# Prueba ciega de las predicciones OCTH para O4 con GWTC-5.0 (O4b): prerregistro

Fecha: 2026-11-01. Aprobado por Frank en la misma instrucción que el prerregistro topológico ("registra el preregistro ... y corre").

## Ceguera

- El gist "PHG O4 Predictions - Timestamped January 9, 2026" (gist 0ddbb1a15441ca2e6c998ae7b1690163, archivo OCTH_O4_Predictions.md) tiene una sola versión, b05c2577, creada el 2026-01-09T23:46:35Z. No fue modificado.
- GWTC-4.0 (O4a) se publicó antes del gist. NO es prueba ciega y no se usa para decidir.
- GWTC-5.0 (O4b, 161 eventos) se agregó a GWOSC el 2026-05-22 (date_added de los parámetros) y los productos de PE salieron en Zenodo el 2026-05-26. Ambos son posteriores al gist: es una prueba ciega de verdad.
- Este análisis lo hace quien escribió la auditoría (sin ceguera del analista). Por eso las reglas se fijan aquí antes de abrir mapas de cielo o strain.

## Convención de masa

- Primaria: M = masa total en el marco de la fuente (`total_mass_source`, mediana de PE), como en MASS_CATALOG del repositorio (GW150914 = 65.3).
- f1 = 2200 / M Hz. Modos: f1 x {1, sqrt 3, 2, sqrt 7}.
- Secundaria (solo reporte): masa en el marco del detector, M (1 + z).

## LT1. Agrupamiento en el cielo (gist, sección 2.4)

- Predicción: al menos 60 % de los eventos bien localizados (región creíble del 90 % < 500 deg²) caen a menos de 30 grados de una de estas posiciones:
  - (RA 12h51m, Dec -27)
  - (RA 3h15m, Dec -19)
  - (RA 11h, Dec +25)
  - (l 220, b -20)
- Muestra: eventos de GWTC-5.0 con mapa de cielo de PE en Archived_Skymaps y área del 90 % < 500 deg². Si un evento tiene varios mapas, se usa el preferido o, si no se puede determinar, el primero en orden alfabético.
- Medición primaria: si el píxel de máxima probabilidad está a menos de 30 grados de alguna posición.
- Medición secundaria: si la probabilidad posterior dentro de la unión de los 4 casquetes es mayor que 0.5.
- Nulo: f0 = área de la unión de los casquetes / 4 pi, calculada con HEALPix a NSIDE 256. Prueba binomial de una cola.
- Decisión:
  - la predicción se cumple si la fracción observada es >= 0.60;
  - el criterio de falsación del gist (sección 4.4) se cumple si Z < 2 respecto a f0.
- Limitación declarada: la red de detectores no tiene sensibilidad isótropa. f0 es una aproximación; se reporta como tal.

## LT2. Potencia en las frecuencias hexagonales (gist, secciones 2.1, 2.2 y 4.1)

- Muestra: eventos de GWTC-5.0 con BBH (`mass_2_source` > 3), SNR de red >= 12 y strain público de H1 y/o L1 a 4 kHz.
- Ventana en la fuente: [t_c - 0.1 s, t_c + 0.4 s], con ventana de Tukey (alpha = 0.2). t_c es el GPS del evento en GWOSC.
- Blanqueo: PSD por Welch (segmentos de 4 s, 50 % de traslape) sobre ±64 s alrededor del evento, excluyendo ±4 s.
- Estadístico por detector: S = media sobre los 4 modos de la potencia blanqueada promediada en la banda f_n (1 ± 0.05). Los modos fuera de [20, 1000] Hz se descartan; si no queda ninguno, el evento se excluye.
- Nulo por evento: el mismo S en 200 ventanas fuera de la fuente, de 0.5 s cada una, tomadas del mismo archivo de strain a más de 8 s del evento y sin segmentos marcados como inválidos. p_evento = (1 + #{S_off >= S_on}) / 201. Si hay dos detectores, se combinan por Fisher.
- Resultado combinado: fracción de eventos con p < 0.05 (esperado 5 % bajo el nulo) y Fisher sobre todos los eventos.
- Control positivo: la misma maquinaria en la banda de GR f_220 (1 ± 0.15), con f_220 = 32.3 kHz x F(a_f) / M_f,det y F(a) = 1.5251 - 1.1568 (1 - a)^0.1292 (ajuste de Berti, Cardoso y Will 2006, gr-qc/0512160), F(0.69) = 0.5307, espín final fijo a 0.69. Debe dar exceso en los eventos con más SNR. Si no lo da, el pipeline no es sensible y LT2 se declara no concluyente.
- Decisión:
  - criterio del gist 4.1: OCTH queda falsada si menos del 60 % de los eventos muestra la señal (aquí, p < 0.05);
  - criterio del gist 4.2 (independencia de masa): se reporta.

## LT3. Modos de ringdown medidos por LVK (gist, resumen: "ringdown frequencies will exhibit ratios of 1 : sqrt 3 : 2 : sqrt 7")

- Fuente: Zenodo 21454847, "Data release for GWTC-5.0: Tests of General Relativity", parte RD (ringdown).
- Para cada evento con frecuencias de QNM medidas: se compara la mediana y el intervalo del 90 % de f_220 contra el conjunto OCTH f1 x {1, sqrt 3, 2, sqrt 7} (± 5 %). Si hay modos medidos, también sus cocientes.
- Decisión: si f_220 medido cae fuera de todas las bandas OCTH en la mayoría de los eventos, la predicción del resumen no se cumple.
- Si el formato del archivo no permite extraer f_220 sin ambigüedad, LT3 se reporta como no ejecutado.

## Semilla y reproducibilidad

- Semilla 20261101 para elegir las ventanas fuera de la fuente.
- Hashes SHA-256 de todas las entradas descargadas y de las salidas.

## Enmiendas (antes de ver resultados agregados)

- A1 (2026-11-01): el GPS de GWOSC tiene resolución de 0.1 s; t_c se toma tal cual. La ventana de 0.5 s absorbe ese error.
- A2 (2026-11-01): filtro paso-alto Butterworth de orden 8 a 15 Hz, de fase cero, sobre el strain completo antes de cortar ventanas. Se agregó tras inspeccionar UN evento (GW240413_022019) durante la depuración. La mediana fuera de la fuente variaba entre 4 y 89 por fuga espectral del ruido sísmico, sin relación con la señal. Las ventanas fuera de la fuente exigen datos válidos a ±4 s (antes ±1 s), por los transitorios del filtro en los huecos.
- A3: diagnóstico secundario NO prerregistrado, potencia de banda ancha 30-300 Hz, como control positivo adicional de sensibilidad a la señal de CBC. No entra en la decisión.

## Resultado LT2 prerregistrado y prueba exploratoria LT2b (registrada antes de correrla)

- LT2, tal como se prerregistró: exceso en las bandas OCTH (masa de la fuente) en 24/35 eventos. El diagnóstico de banda ancha 30-300 Hz da exceso en 33/35, y el control GR f_220 en 12/35. El estadístico no es específico: cualquier señal de CBC deposita potencia en esas bandas durante [t_c - 0.1, t_c + 0.4] s. El "cumplimiento" del criterio 4.1 del gist no distingue OCTH de GR. Es un defecto del diseño de LT2, no evidencia.
- LT2b (EXPLORATORIA, post hoc respecto de LT2, definida sin haber calculado este estadístico):
  - En la ventana de la fuente y con la misma PSD y el mismo filtro, se compara S en el peine {1, sqrt3, 2, sqrt7} x f1 contra 500 peines {1, r2, r3, r4}, con r_i uniformes e independientes en [1.3, 2.9] (semilla 20261102).
  - p_spec = (1 + #{S_ctrl >= S_OCTH}) / 501, combinado por Fisher entre detectores.
  - Si las razones hexagonales no son especiales, p_spec es uniforme: unos 5 % de eventos con p < 0.05.
  - Variante sin la fundamental: el peine {sqrt3, 2, sqrt7} contra {r2, r3, r4}.
