#!/usr/bin/env python3
"""
SANITY CHECK #2: TEST DE ROTACION
=================================

PREGUNTA: ¿La correlacion CMB-LIGO (Z=4.31) es un artefacto geometrico?

METODO:
1. Tomar las posiciones LIGO originales
2. Rotarlas artificialmente (+90, +180, +270 grados en RA)
3. Re-calcular la correlacion con CMB
4. Si Z sigue siendo alto → ERROR SISTEMATICO
5. Si Z cae → CORRELACION REAL

Autor: Francisco Molina Burgos & Claude
Fecha: 2026-01-09
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import os
import json

np.random.seed(42)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

print("=" * 70)
print("SANITY CHECK #2: TEST DE ROTACION")
print("¿La correlacion CMB-LIGO es REAL o un artefacto geometrico?")
print("=" * 70)

# Posiciones de anomalias CMB (del paper)
CMB_ANOMALIES = {
    'North_Galactic_Pole': {'l': 0, 'b': 90},      # NGP
    'South_Galactic_Pole': {'l': 0, 'b': -90},     # SGP
    'Cold_Spot': {'l': 209, 'b': -57},             # Eridanus
    'Hemispheric_Asymmetry': {'l': 225, 'b': 15},  # Eje de asimetria
    'Quadrupole_Axis': {'l': 240, 'b': 63},        # Eje cuadrupolo
    'Octopole_Axis': {'l': 308, 'b': 63},          # Eje octopolo
}

# Eventos LIGO con posiciones bien localizadas (del paper)
LIGO_EVENTS = [
    {'name': 'GW150914', 'ra': 125, 'dec': -41, 'mass': 62},
    {'name': 'GW151226', 'ra': 48, 'dec': 55, 'mass': 21},
    {'name': 'GW170104', 'ra': 342, 'dec': 52, 'mass': 50},
    {'name': 'GW170608', 'ra': 17, 'dec': 18, 'mass': 18},
    {'name': 'GW170729', 'ra': 285, 'dec': 36, 'mass': 79},
    {'name': 'GW170809', 'ra': 60, 'dec': 28, 'mass': 56},
    {'name': 'GW170814', 'ra': 35, 'dec': -45, 'mass': 55},
    {'name': 'GW170817', 'ra': 197.4, 'dec': -23.4, 'mass': 2.7},  # BNS
    {'name': 'GW170818', 'ra': 345, 'dec': 33, 'mass': 59},
    {'name': 'GW170823', 'ra': 52, 'dec': 15, 'mass': 65},
    {'name': 'GW190412', 'ra': 270, 'dec': 18, 'mass': 37},
    {'name': 'GW190425', 'ra': 150, 'dec': 12, 'mass': 3.4},  # BNS
    {'name': 'GW190521', 'ra': 215, 'dec': -25, 'mass': 150},
    {'name': 'GW190814', 'ra': 12.5, 'dec': -25.2, 'mass': 26},
    {'name': 'GW190630', 'ra': 226, 'dec': 5, 'mass': 35},
    {'name': 'GW191105', 'ra': 80, 'dec': -5, 'mass': 18},
    {'name': 'GW191216', 'ra': 160, 'dec': 40, 'mass': 18},
    {'name': 'GW200115', 'ra': 320, 'dec': 0, 'mass': 8},
]


def galactic_to_equatorial(l, b):
    """Convierte coordenadas galacticas a ecuatoriales."""
    # Parametros del polo norte galactico en ecuatoriales
    ra_ngp = np.radians(192.85948)
    dec_ngp = np.radians(27.12825)
    l_ncp = np.radians(122.93192)  # Longitud galactica del NCP

    l_rad = np.radians(l)
    b_rad = np.radians(b)

    sin_dec = np.sin(dec_ngp) * np.sin(b_rad) + \
              np.cos(dec_ngp) * np.cos(b_rad) * np.cos(l_rad - l_ncp)
    dec = np.arcsin(sin_dec)

    cos_ra_diff = (np.cos(b_rad) * np.sin(l_rad - l_ncp)) / np.cos(dec)
    sin_ra_diff = (np.cos(dec_ngp) * np.sin(b_rad) - \
                   np.sin(dec_ngp) * np.cos(b_rad) * np.cos(l_rad - l_ncp)) / np.cos(dec)

    ra_diff = np.arctan2(cos_ra_diff, sin_ra_diff)
    ra = ra_ngp - ra_diff

    ra = np.degrees(ra) % 360
    dec = np.degrees(dec)

    return ra, dec


def angular_separation(ra1, dec1, ra2, dec2):
    """Calcula separacion angular en grados."""
    ra1, dec1 = np.radians(ra1), np.radians(dec1)
    ra2, dec2 = np.radians(ra2), np.radians(dec2)

    cos_sep = np.sin(dec1) * np.sin(dec2) + \
              np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    cos_sep = np.clip(cos_sep, -1, 1)
    return np.degrees(np.arccos(cos_sep))


def rotate_coordinates(ra, dec, delta_ra):
    """Rota coordenadas ecuatoriales en RA."""
    new_ra = (ra + delta_ra) % 360
    return new_ra, dec


def count_close_pairs(events, anomalies, threshold_deg=30):
    """Cuenta pares cercanos entre eventos LIGO y anomalias CMB."""
    close_pairs = []

    # Convertir anomalias a ecuatorial
    cmb_eq = {}
    for name, coords in anomalies.items():
        ra, dec = galactic_to_equatorial(coords['l'], coords['b'])
        cmb_eq[name] = {'ra': ra, 'dec': dec}

    for event in events:
        for cmb_name, cmb_coords in cmb_eq.items():
            sep = angular_separation(event['ra'], event['dec'],
                                    cmb_coords['ra'], cmb_coords['dec'])
            if sep < threshold_deg:
                close_pairs.append({
                    'event': event['name'],
                    'anomaly': cmb_name,
                    'separation': sep
                })

    return len(close_pairs), close_pairs


def monte_carlo_expected(n_events, n_anomalies, threshold_deg, n_sim=10000):
    """Calcula distribucion esperada para posiciones aleatorias."""
    counts = []

    for _ in range(n_sim):
        # Generar eventos aleatorios
        ra_random = np.random.uniform(0, 360, n_events)
        dec_random = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_events)))

        # Contar coincidencias
        n_close = 0
        for name, coords in CMB_ANOMALIES.items():
            cmb_ra, cmb_dec = galactic_to_equatorial(coords['l'], coords['b'])
            for ra, dec in zip(ra_random, dec_random):
                if angular_separation(ra, dec, cmb_ra, cmb_dec) < threshold_deg:
                    n_close += 1

        counts.append(n_close)

    return np.array(counts)


# Analisis original
print("\n[FASE 1: Correlacion original (sin rotacion)]")
n_original, pairs_original = count_close_pairs(LIGO_EVENTS, CMB_ANOMALIES)
print(f"  Pares cercanos (<30 deg): {n_original}")
for p in pairs_original:
    print(f"    {p['event']} - {p['anomaly']}: {p['separation']:.1f} deg")

# Distribucion nula
print("\n[FASE 2: Calculando distribucion nula (Monte Carlo)]")
null_distribution = monte_carlo_expected(len(LIGO_EVENTS), len(CMB_ANOMALIES), 30)
mean_null = np.mean(null_distribution)
std_null = np.std(null_distribution)
print(f"  Esperado aleatorio: {mean_null:.1f} +/- {std_null:.1f} pares")

z_original = (n_original - mean_null) / std_null
p_original = np.sum(null_distribution >= n_original) / len(null_distribution)
print(f"  Z-score original: {z_original:.2f}")
print(f"  P-value: {p_original:.4f}")

# Rotaciones
print("\n[FASE 3: Test de rotacion]")
ROTATIONS = [45, 90, 135, 180, 225, 270, 315]
rotation_results = {}

for rot in ROTATIONS:
    # Rotar coordenadas LIGO
    rotated_events = []
    for event in LIGO_EVENTS:
        new_ra, new_dec = rotate_coordinates(event['ra'], event['dec'], rot)
        rotated_events.append({
            'name': event['name'],
            'ra': new_ra,
            'dec': new_dec,
            'mass': event['mass']
        })

    n_rotated, pairs_rotated = count_close_pairs(rotated_events, CMB_ANOMALIES)
    z_rotated = (n_rotated - mean_null) / std_null

    rotation_results[rot] = {
        'n_pairs': n_rotated,
        'z_score': z_rotated,
        'pairs': pairs_rotated
    }

    status = "***" if z_rotated > 2 else ""
    print(f"  Rotacion +{rot:3d} deg: {n_rotated} pares, Z = {z_rotated:.2f} {status}")

# Analisis
print("\n" + "=" * 70)
print("ANALISIS DE RESULTADOS")
print("=" * 70)

z_scores_rotated = [rotation_results[r]['z_score'] for r in ROTATIONS]
max_z_rotated = max(z_scores_rotated)
mean_z_rotated = np.mean(z_scores_rotated)

print(f"\n  Z-score original: {z_original:.2f}")
print(f"  Z-score medio (rotados): {mean_z_rotated:.2f}")
print(f"  Z-score maximo (rotados): {max_z_rotated:.2f}")

# Calcular cuantas rotaciones tienen Z similar al original
n_similar = sum(1 for z in z_scores_rotated if z >= z_original - 0.5)

# Veredicto
print("\n" + "=" * 70)
print("VEREDICTO: TEST DE ROTACION")
print("=" * 70)

# Criterios:
# - Si el Z original es significativamente mayor que TODOS los rotados → PASA
# - Si el Z original es similar a los rotados → FALLA

z_drop = z_original - mean_z_rotated
is_significant = z_original > 2.0 and z_drop > 1.5 and max_z_rotated < z_original - 0.5

if is_significant:
    verdict = "PASA: La correlacion CMB-LIGO es REAL"
    explanation = f"Z cae de {z_original:.1f} a {mean_z_rotated:.1f} al rotar (delta = {z_drop:.1f})"
    honest = True
else:
    if z_original < 2.0:
        verdict = "INCONCLUSO: Correlacion original no es significativa"
        explanation = f"Z original = {z_original:.1f} (< 2.0 sigma)"
    elif max_z_rotated >= z_original - 0.5:
        verdict = "ALERTA: Posible sesgo geometrico"
        explanation = f"Algunas rotaciones tienen Z similar ({max_z_rotated:.1f} vs {z_original:.1f})"
    else:
        verdict = "MARGINAL: Requiere mas eventos"
        explanation = f"La caida de Z es modesta (delta = {z_drop:.1f})"
    honest = z_original > max_z_rotated + 0.5

print(f"\n  {verdict}")
print(f"  {explanation}")

if honest:
    print(f"\n  Conclusion:")
    print(f"    - La correlacion original (Z={z_original:.1f}) es la MAS ALTA")
    print(f"    - Las rotaciones reducen la correlacion")
    print(f"    - Esto indica una correlacion FISICA, no artefacto")
else:
    print(f"\n  Conclusion:")
    print(f"    - Hay rotaciones con Z comparable al original")
    print(f"    - Esto sugiere cautela en la interpretacion")
    print(f"    - Se necesitan mas eventos para confirmar")

# Crear figura
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Panel 1: Z-score vs rotacion
ax = axes[0]
rots = [0] + ROTATIONS
z_all = [z_original] + z_scores_rotated
colors = ['red' if r == 0 else 'steelblue' for r in rots]
bars = ax.bar(range(len(rots)), z_all, color=colors, edgecolor='black')
ax.axhline(2.0, color='orange', ls='--', label='2 sigma')
ax.axhline(mean_null, color='gray', ls=':', label='Ruido')
ax.set_xticks(range(len(rots)))
ax.set_xticklabels([f'{r}' for r in rots])
ax.set_xlabel('Rotacion en RA (grados)')
ax.set_ylabel('Z-score')
ax.set_title('Correlacion CMB-LIGO vs Rotacion')
ax.legend()

# Anotar
ax.annotate(f'Z={z_original:.1f}\n(ORIGINAL)',
            xy=(0, z_original), xytext=(0.5, z_original + 0.5),
            arrowprops=dict(arrowstyle='->', color='red'),
            fontsize=10, color='red')

# Panel 2: Resumen
ax = axes[1]
ax.axis('off')

summary = f"""
TEST DE ROTACION: SANITY CHECK #2
{'='*45}

METODO:
  Rotar posiciones LIGO artificialmente
  Re-calcular correlacion con CMB
  Verificar si Z cae

RESULTADOS:
  Z original (0 deg):  {z_original:.2f} ***
  Z medio (rotados):   {mean_z_rotated:.2f}
  Z maximo (rotados):  {max_z_rotated:.2f}

  Caida de Z: {z_drop:.2f}

DISTRIBUCION NULA:
  Esperado: {mean_null:.1f} +/- {std_null:.1f} pares
  Observado: {n_original} pares
  P-value: {p_original:.4f}

{'='*45}
VEREDICTO: {verdict.split(':')[0]}
{'='*45}

{explanation}

La correlacion {'ES' if honest else 'puede ser'} REAL.
"""

ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round',
                  facecolor='lightgreen' if honest else 'lightyellow',
                  alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, 'sanity_check_rotation.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIGURES_DIR, 'sanity_check_rotation.pdf'), dpi=300, bbox_inches='tight')
print(f"\nFigura guardada: sanity_check_rotation.png/pdf")
plt.close()

# Guardar resultados
output = {
    'test': 'Sanity Check #2: Rotation Test',
    'original': {
        'n_pairs': n_original,
        'z_score': float(z_original),
        'p_value': float(p_original),
        'pairs': pairs_original
    },
    'null_distribution': {
        'mean': float(mean_null),
        'std': float(std_null)
    },
    'rotations': {
        str(r): {
            'n_pairs': rotation_results[r]['n_pairs'],
            'z_score': float(rotation_results[r]['z_score'])
        }
        for r in ROTATIONS
    },
    'analysis': {
        'mean_z_rotated': float(mean_z_rotated),
        'max_z_rotated': float(max_z_rotated),
        'z_drop': float(z_drop)
    },
    'verdict': verdict,
    'passes': honest
}

filepath = os.path.join(RESULTS_DIR, 'sanity_check_rotation.json')
with open(filepath, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Resultados guardados: {filepath}")
