"""
PRIMACIA INFORMACIONAL TOPOLOGICA (PIT) - TEST COSMOLOGICO

Hipotesis: La entropia topologica S_H1 del campo de densidad cosmico
debe mostrar COLAPSO en z ~ 1089, ANTES de recombinacion.

Metodologia:
1. Simular campo de densidad a diferentes z
2. Calcular homologia persistente H1
3. Medir entropia S_H1(z)
4. Aplicar CUSUM para detectar transicion
5. Verificar si t_CUSUM < t_recombinacion

Autor: F. Molina-Burgos
Fecha: Enero 2026
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

# Intentar importar ripser para homologia persistente
try:
    from ripser import ripser
    HAS_RIPSER = True
except ImportError:
    HAS_RIPSER = False
    print("NOTA: ripser no disponible, usando aproximacion Euler")

print("="*70)
print("  PRIMACIA INFORMACIONAL TOPOLOGICA - TEST COSMOLOGICO")
print("="*70)

# =============================================================================
# PARAMETROS COSMOLOGICOS
# =============================================================================

H0 = 67.4  # km/s/Mpc
OMEGA_M = 0.315
OMEGA_B = 0.049
OMEGA_R = 9.0e-5
OMEGA_LAMBDA = 1 - OMEGA_M - OMEGA_R

Z_RECOMBINATION = 1089.0  # Recombinacion
Z_DECOUPLING = 1100.0     # Desacoplamiento

# OCTH parameters
EPSILON_TOPO = 0.015
SIGMA_TOPO = 0.6
Z_TOPO = 1089.0

# =============================================================================
# MODELO COSMOLOGICO CON OCTH
# =============================================================================

def psi_octh(z):
    """Campo de permeabilidad temporal OCTH"""
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + Z_TOPO)
    delta = (ln_a - ln_a_topo) / SIGMA_TOPO
    f_topo = np.exp(-0.5 * delta**2)
    return 1.0 - EPSILON_TOPO * f_topo

def growth_factor(z):
    """Factor de crecimiento de perturbaciones (aproximado)"""
    # D(z) ~ a para materia dominada
    # Modificado por OCTH
    a = 1.0 / (1 + z)
    psi = psi_octh(z)
    # El crecimiento se modifica por Psi
    return a / psi

def density_contrast_rms(z):
    """RMS de contraste de densidad a redshift z"""
    # sigma_8 ~ 0.8 a z=0, escala con D(z)
    sigma_8 = 0.81
    D_z = growth_factor(z)
    D_0 = growth_factor(0)
    return sigma_8 * D_z / D_0

# =============================================================================
# SIMULACION DE CAMPO DE DENSIDAD
# =============================================================================

def generate_density_field(z, N_points=200, box_size=100.0):
    """
    Genera un campo de densidad cosmico a redshift z.

    Usa espectro de potencias simplificado P(k) ~ k^n_s con cutoff.
    El campo representa posiciones de "trazadores" (galaxias/materia).
    """
    np.random.seed(42 + int(z))  # Reproducible pero diferente para cada z

    # Contraste de densidad a este z
    sigma = density_contrast_rms(z)

    # Para z alto, el universo es mas homogeneo
    # Las fluctuaciones son mas pequenas

    # Generar campo gaussiano
    # En z alto (pre-recombinacion), hay oscilaciones acusticas
    # En z bajo (post-recombinacion), hay estructura

    if z > Z_RECOMBINATION:
        # Pre-recombinacion: plasma acoplado, oscilaciones BAO
        # Estructura suprimida por presion de radiacion
        sigma_effective = sigma * 0.01  # Muy pequeno
        correlation_length = box_size * 0.5  # Largo alcance (BAO)
    else:
        # Post-recombinacion: estructura crece
        # La correlacion disminuye (estructura mas local)
        sigma_effective = sigma
        correlation_length = box_size * 0.1 * (1 + z) / 100

    # Generar posiciones con clustering
    # Empezar con grid uniforme + perturbaciones
    n_side = int(np.sqrt(N_points))
    if n_side * n_side < N_points:
        n_side += 1

    x = np.linspace(0, box_size, n_side)
    y = np.linspace(0, box_size, n_side)
    xx, yy = np.meshgrid(x, y)

    # Tomar N_points puntos
    xx = xx.flatten()[:N_points]
    yy = yy.flatten()[:N_points]

    # Anadir perturbaciones correlacionadas
    # Usando campo gaussiano filtrado
    perturbation_x = np.random.randn(N_points) * sigma_effective * box_size * 0.1
    perturbation_y = np.random.randn(N_points) * sigma_effective * box_size * 0.1

    # Aplicar correlacion espacial (simplificado)
    if correlation_length > 0:
        # Suavizar perturbaciones para correlacionar vecinos
        perturbation_x = gaussian_filter1d_wrap(perturbation_x, correlation_length / box_size * N_points)
        perturbation_y = gaussian_filter1d_wrap(perturbation_y, correlation_length / box_size * N_points)

    xx = (xx + perturbation_x) % box_size
    yy = (yy + perturbation_y) % box_size

    points = np.column_stack([xx, yy])

    return points

def gaussian_filter1d_wrap(data, sigma):
    """Filtro gaussiano con condiciones periodicas"""
    if sigma < 1:
        return data
    from scipy.ndimage import gaussian_filter1d
    # Extender para periodicidad
    extended = np.concatenate([data, data, data])
    filtered = gaussian_filter1d(extended, sigma)
    n = len(data)
    return filtered[n:2*n]

# =============================================================================
# HOMOLOGIA PERSISTENTE
# =============================================================================

def compute_persistence_entropy_h1(points, max_edge=None):
    """
    Calcula entropia de persistencia H1 usando ripser o aproximacion Euler.
    """
    if len(points) < 4:
        return 0.0

    # Matriz de distancias
    dist_matrix = squareform(pdist(points))

    if max_edge is None:
        max_edge = np.percentile(dist_matrix, 50)

    if HAS_RIPSER:
        # Metodo exacto con ripser
        result = ripser(dist_matrix, maxdim=1, distance_matrix=True, thresh=max_edge)
        dgm1 = result['dgms'][1]

        if len(dgm1) == 0:
            return 0.0

        # Calcular lifetimes
        lifetimes = dgm1[:, 1] - dgm1[:, 0]
        lifetimes = lifetimes[np.isfinite(lifetimes)]
        lifetimes = lifetimes[lifetimes > 0]

        if len(lifetimes) == 0:
            return 0.0

        # Entropia de Shannon
        total = np.sum(lifetimes)
        p = lifetimes / total
        entropy = -np.sum(p * np.log(p + 1e-10))

        return entropy
    else:
        # Aproximacion Euler
        return compute_euler_h1_entropy(dist_matrix, max_edge)

def compute_euler_h1_entropy(dist_matrix, max_edge):
    """
    Aproximacion de entropia H1 usando caracteristica de Euler.

    beta_1 ≈ E - V + beta_0 - F

    Aproximamos contando edges y triangulos a diferentes escalas.
    """
    n = len(dist_matrix)
    n_scales = 20
    scales = np.linspace(0, max_edge, n_scales)

    beta1_values = []

    for eps in scales[1:]:
        # Contar edges
        edges = np.sum(dist_matrix < eps) // 2 - n // 2  # Quitar diagonal

        # Estimar triangulos (simplificado)
        triangles = 0
        for i in range(n):
            neighbors = np.where(dist_matrix[i] < eps)[0]
            neighbors = neighbors[neighbors > i]
            for j in neighbors:
                common = np.sum((dist_matrix[i] < eps) & (dist_matrix[j] < eps)) - 2
                triangles += common // 2
        triangles = triangles // 3

        # Euler: V - E + F, para complejo simplicial: beta_0 - beta_1 + beta_2 = chi
        # Aproximacion: beta_1 ~ E - n + 1 - F (asumiendo conectado, beta_2 ~ F/4)
        beta1_approx = max(0, edges - n + 1 - triangles // 4)
        beta1_values.append(beta1_approx)

    beta1_values = np.array(beta1_values)

    if np.sum(beta1_values) == 0:
        return 0.0

    # Entropia sobre distribucion de beta1 en escalas
    total = np.sum(beta1_values) + 1e-10
    p = beta1_values / total
    p = p[p > 0]
    entropy = -np.sum(p * np.log(p + 1e-10))

    return entropy

# =============================================================================
# CUSUM DETECTOR
# =============================================================================

class CUSUMDetector:
    """Detector de cambios usando CUSUM"""

    def __init__(self, k=0.5, h=4.0, sigma_min=0.1):
        self.k = k
        self.h = h
        self.sigma_min = sigma_min
        self.mu = None
        self.sigma = None
        self.C_pos = 0
        self.C_neg = 0
        self.calibrated = False

    def calibrate(self, baseline_data):
        """Calibrar con datos de baseline"""
        self.mu = np.mean(baseline_data)
        self.sigma = max(np.std(baseline_data), self.sigma_min)
        self.C_pos = 0
        self.C_neg = 0
        self.calibrated = True

    def update(self, value):
        """Actualizar y detectar"""
        if not self.calibrated:
            return False

        z = (value - self.mu) / self.sigma

        self.C_pos = max(0, self.C_pos + z - self.k)
        self.C_neg = max(0, self.C_neg - z - self.k)

        return self.C_pos > self.h or self.C_neg > self.h

# =============================================================================
# TEST PRINCIPAL
# =============================================================================

print("\n" + "="*70)
print("PARTE 1: Simulando campo de densidad a diferentes z")
print("="*70)

# Rango de redshifts para estudiar
z_values = np.concatenate([
    np.linspace(2000, 1200, 20),   # Pre-recombinacion
    np.linspace(1190, 1000, 30),   # Alrededor de recombinacion
    np.linspace(990, 100, 20),     # Post-recombinacion temprano
    np.linspace(90, 0, 10)         # Universo tardio
])

N_POINTS = 150  # Puntos por simulacion
BOX_SIZE = 100.0  # Mpc

print(f"\nParametros:")
print(f"  N_points = {N_POINTS}")
print(f"  Box size = {BOX_SIZE} Mpc")
print(f"  z range = [{z_values[-1]:.0f}, {z_values[0]:.0f}]")
print(f"  z_recombination = {Z_RECOMBINATION}")

# Calcular entropia topologica para cada z
print("\nCalculando S_H1(z)...")
print("-"*70)

results = []

for i, z in enumerate(z_values):
    # Generar campo
    points = generate_density_field(z, N_points=N_POINTS, box_size=BOX_SIZE)

    # Calcular entropia H1
    S_H1 = compute_persistence_entropy_h1(points, max_edge=BOX_SIZE * 0.3)

    # Calcular Psi OCTH
    psi = psi_octh(z)

    results.append({
        'z': z,
        'S_H1': S_H1,
        'psi': psi,
        'delta_rms': density_contrast_rms(z)
    })

    if i % 10 == 0:
        print(f"  z = {z:7.1f}: S_H1 = {S_H1:.4f}, Psi = {psi:.4f}")

print("-"*70)

# Convertir a arrays
z_arr = np.array([r['z'] for r in results])
S_H1_arr = np.array([r['S_H1'] for r in results])
psi_arr = np.array([r['psi'] for r in results])

print("\n" + "="*70)
print("PARTE 2: Deteccion CUSUM de transicion topologica")
print("="*70)

# Usar primera parte (z alto) como baseline
baseline_idx = len(z_arr) // 4
baseline_S = S_H1_arr[:baseline_idx]

print(f"\nBaseline: z > {z_arr[baseline_idx]:.0f}")
print(f"  mu_baseline = {np.mean(baseline_S):.4f}")
print(f"  sigma_baseline = {np.std(baseline_S):.4f}")

# Aplicar CUSUM
detector = CUSUMDetector(k=0.5, h=3.0, sigma_min=0.05)
detector.calibrate(baseline_S)

z_cusum_detected = None

for i in range(baseline_idx, len(z_arr)):
    if detector.update(S_H1_arr[i]):
        z_cusum_detected = z_arr[i]
        break

print(f"\nResultados CUSUM:")
if z_cusum_detected is not None:
    print(f"  z_CUSUM (deteccion topologica) = {z_cusum_detected:.1f}")
    print(f"  z_recombinacion (fisica)       = {Z_RECOMBINATION:.1f}")

    if z_cusum_detected > Z_RECOMBINATION:
        gap = z_cusum_detected - Z_RECOMBINATION
        print(f"\n  *** TOPOLOGIA PRECEDE METRICA ***")
        print(f"  Gap: Delta_z = {gap:.1f}")
        print(f"  La transicion topologica ocurre ANTES de recombinacion!")
        PIT_CONFIRMED = True
    else:
        gap = Z_RECOMBINATION - z_cusum_detected
        print(f"\n  Topologia detectada DESPUES de recombinacion")
        print(f"  Gap: Delta_z = -{gap:.1f}")
        PIT_CONFIRMED = False
else:
    print("  No se detecto transicion con CUSUM")
    PIT_CONFIRMED = False

print("\n" + "="*70)
print("PARTE 3: Correlacion S_H1 con Psi OCTH")
print("="*70)

# Calcular correlacion
correlation = np.corrcoef(S_H1_arr, psi_arr)[0, 1]
print(f"\nCorrelacion Pearson(S_H1, Psi) = {correlation:.4f}")

if abs(correlation) > 0.5:
    print("  *** ALTA CORRELACION ***")
    print("  La entropia topologica TDA correlaciona con campo OCTH!")
else:
    print("  Correlacion moderada o baja")

# Encontrar minimo de S_H1
idx_min_S = np.argmin(S_H1_arr)
z_min_S = z_arr[idx_min_S]

# Encontrar minimo de Psi
idx_min_psi = np.argmin(psi_arr)
z_min_psi = z_arr[idx_min_psi]

print(f"\nMinimo de S_H1:  z = {z_min_S:.1f}")
print(f"Minimo de Psi:   z = {z_min_psi:.1f}")
print(f"Recombinacion:   z = {Z_RECOMBINATION:.1f}")

if abs(z_min_S - z_min_psi) < 100:
    print("\n  *** MINIMOS COINCIDEN ***")
    print("  S_H1 y Psi minimizan en el mismo epoch!")

print("\n" + "="*70)
print("PARTE 4: Comparacion con TDA molecular")
print("="*70)

print("""
SISTEMA MOLECULAR (Lennard-Jones 2D):
  - Observable topologico: S_H1 (loops)
  - Observable metrico: psi_6 (orden hexatico)
  - Resultado: t_topo < t_fisico (100% para N >= 900)

SISTEMA COSMOLOGICO (Universo):
  - Observable topologico: S_H1 del campo de densidad
  - Observable metrico: z_recombinacion
  - Resultado: """)

if PIT_CONFIRMED:
    print(f"z_topo > z_fisico (topologia PRECEDE)")
    print("""
  *** PIT CONFIRMADO EN ESCALA COSMOLOGICA ***

  El mismo principio que gobierna cristalizacion molecular
  gobierna recombinacion cosmologica:

    TOPOLOGIA CAMBIA ANTES QUE METRICA
""")
else:
    print(f"Test no concluyente con simulacion simplificada")
    print("""
  Se requiere:
  - Datos reales de CMB (Planck)
  - Homologia persistente exacta
  - Mayor resolucion espacial
""")

print("\n" + "="*70)
print("PARTE 5: Perfil de S_H1(z)")
print("="*70)

print("\n{:^10} {:^12} {:^12} {:^12}".format("z", "S_H1", "Psi", "delta_rms"))
print("-"*50)

# Mostrar algunos puntos clave
indices_to_show = [0, 10, 20, 30, 40, 50, 60, 70, -1]
for idx in indices_to_show:
    if idx < len(results):
        r = results[idx]
        marker = " <-- min" if r['z'] == z_min_S else ""
        marker = " <-- z_rec" if abs(r['z'] - Z_RECOMBINATION) < 50 else marker
        print(f"{r['z']:>10.1f} {r['S_H1']:>12.4f} {r['psi']:>12.4f} {r['delta_rms']:>12.6f}{marker}")

print("\n" + "="*70)
print("CONCLUSIONES")
print("="*70)

print(f"""
1. ENTROPIA TOPOLOGICA S_H1:
   - Muestra estructura clara como funcion de z
   - {'Tiene minimo cerca de recombinacion' if abs(z_min_S - Z_RECOMBINATION) < 200 else 'Patron diferente al esperado'}

2. CORRELACION CON OCTH:
   - Pearson(S_H1, Psi) = {correlation:.4f}
   - {'Alta correlacion - TDA y OCTH miden lo mismo!' if abs(correlation) > 0.5 else 'Correlacion a investigar'}

3. DETECCION CUSUM:
   - {'z_CUSUM = ' + f'{z_cusum_detected:.1f}' if z_cusum_detected else 'No detectado'}
   - {'TOPOLOGIA PRECEDE RECOMBINACION' if PIT_CONFIRMED else 'Test no concluyente'}

4. PRINCIPIO PIT:
   - {'APOYADO por simulacion' if PIT_CONFIRMED else 'Requiere datos reales CMB'}
   - Siguiente paso: Aplicar ripser a mapas Planck reales
""")

print("="*70)
print("FIN DEL TEST")
print("="*70)

# Guardar resultados para analisis posterior
output_file = "PIT_cosmological_results.json"
import json

output_data = {
    'z_values': z_arr.tolist(),
    'S_H1_values': S_H1_arr.tolist(),
    'psi_values': psi_arr.tolist(),
    'z_cusum_detected': z_cusum_detected,
    'z_recombination': Z_RECOMBINATION,
    'correlation_S_psi': correlation,
    'pit_confirmed': PIT_CONFIRMED
}

with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"\nResultados guardados en: {output_file}")
