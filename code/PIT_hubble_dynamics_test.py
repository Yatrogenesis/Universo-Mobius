"""
PIT TEST v2: TDA sobre DINAMICA de Hubble

En vez de campo espacial, aplicamos TDA a la TRAYECTORIA H(z)
usando delay embedding (reconstruccion de Takens).

Esto es analogo a:
- Brusselator: TDA sobre trayectoria (X,Y)
- Kuramoto: TDA sobre fases theta(t)
- FitzHugh-Nagumo: TDA sobre (v,w)

Para cosmologia:
- Trayectoria: H(z), dH/dz, d2H/dz2
- O: rho(z), P(z), w(z)

Autor: F. Molina-Burgos
Fecha: Enero 2026
"""

import numpy as np
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

try:
    from ripser import ripser
    HAS_RIPSER = True
except ImportError:
    HAS_RIPSER = False

print("="*70)
print("  PIT TEST v2: TDA sobre DINAMICA de Hubble H(z)")
print("="*70)

# =============================================================================
# PARAMETROS
# =============================================================================

H0 = 67.4
OMEGA_M = 0.315
OMEGA_R = 9.0e-5
OMEGA_LAMBDA = 1 - OMEGA_M - OMEGA_R

# OCTH
EPSILON = 0.015
SIGMA = 0.6
Z_TOPO = 1089.0

Z_REC = 1089.0

# =============================================================================
# MODELO COSMOLOGICO
# =============================================================================

def psi_octh(z):
    """Campo Psi OCTH"""
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + Z_TOPO)
    delta = (ln_a - ln_a_topo) / SIGMA
    return 1.0 - EPSILON * np.exp(-0.5 * delta**2)

def H_lcdm(z):
    """Hubble parameter LCDM"""
    opz = 1 + z
    return H0 * np.sqrt(OMEGA_R * opz**4 + OMEGA_M * opz**3 + OMEGA_LAMBDA)

def H_octh(z):
    """Hubble parameter OCTH"""
    return H_lcdm(z) / psi_octh(z)

def dH_dz_lcdm(z):
    """Derivada de H respecto a z (LCDM)"""
    opz = 1 + z
    H = H_lcdm(z)
    dE2_dz = 4 * OMEGA_R * opz**3 + 3 * OMEGA_M * opz**2
    return H0**2 * dE2_dz / (2 * H)

def dH_dz_octh(z):
    """Derivada de H respecto a z (OCTH)"""
    # d(H/psi)/dz = (dH/dz * psi - H * dpsi/dz) / psi^2
    H = H_lcdm(z)
    psi = psi_octh(z)
    dH = dH_dz_lcdm(z)

    # dpsi/dz
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + Z_TOPO)
    delta = (ln_a - ln_a_topo) / SIGMA
    d_delta_dz = 1 / (SIGMA * (1 + z))
    dpsi_dz = EPSILON * delta * np.exp(-0.5 * delta**2) * d_delta_dz

    return (dH * psi - H * dpsi_dz) / psi**2

def q_deceleration(z, model='lcdm'):
    """Parametro de desaceleracion q = -a*a''/a'^2"""
    if model == 'lcdm':
        H = H_lcdm(z)
        dH = dH_dz_lcdm(z)
    else:
        H = H_octh(z)
        dH = dH_dz_octh(z)

    return (1 + z) * dH / H - 1

# =============================================================================
# DELAY EMBEDDING (Reconstruccion de Takens)
# =============================================================================

def delay_embedding(signal, dim=3, tau=1):
    """
    Reconstruccion de Takens del atractor.

    signal: serie temporal
    dim: dimension de embedding
    tau: delay (en indices)

    Retorna: array de puntos en espacio de fase reconstruido
    """
    n = len(signal)
    m = n - (dim - 1) * tau

    if m <= 0:
        return np.array([])

    embedded = np.zeros((m, dim))
    for i in range(dim):
        embedded[:, i] = signal[i*tau:i*tau + m]

    return embedded

# =============================================================================
# PERSISTENCIA
# =============================================================================

def persistence_entropy_h1(points, max_edge=None):
    """Entropia H1 de puntos"""
    if len(points) < 4:
        return 0.0

    dm = squareform(pdist(points))

    if max_edge is None:
        max_edge = np.percentile(dm, 90)

    if HAS_RIPSER:
        result = ripser(dm, maxdim=1, distance_matrix=True, thresh=max_edge)
        dgm1 = result['dgms'][1]

        if len(dgm1) == 0:
            return 0.0

        lifetimes = dgm1[:, 1] - dgm1[:, 0]
        lifetimes = lifetimes[np.isfinite(lifetimes) & (lifetimes > 0)]

        if len(lifetimes) == 0:
            return 0.0

        total = np.sum(lifetimes)
        p = lifetimes / total
        return -np.sum(p * np.log(p + 1e-10))
    else:
        # Aproximacion simple
        return euler_h1_approx(dm, max_edge)

def euler_h1_approx(dm, max_edge):
    """Aproximacion Euler de beta_1"""
    n = len(dm)
    scales = np.linspace(0.1 * max_edge, max_edge, 10)

    beta1_vals = []
    for eps in scales:
        edges = np.sum(dm < eps) // 2 - n // 2
        # Muy simplificado
        beta1 = max(0, edges - n + 1)
        beta1_vals.append(beta1)

    beta1_vals = np.array(beta1_vals, dtype=float)
    if np.sum(beta1_vals) == 0:
        return 0.0

    total = np.sum(beta1_vals) + 1e-10
    p = beta1_vals / total
    p = p[p > 0]
    return -np.sum(p * np.log(p + 1e-10))

# =============================================================================
# CUSUM
# =============================================================================

class CUSUM:
    def __init__(self, k=0.5, h=4.0, sigma_min=0.01):
        self.k = k
        self.h = h
        self.sigma_min = sigma_min
        self.mu = None
        self.sigma = None
        self.C = 0

    def calibrate(self, data):
        self.mu = np.mean(data)
        self.sigma = max(np.std(data), self.sigma_min)
        self.C = 0

    def update(self, x):
        if self.mu is None:
            return False
        z = (x - self.mu) / self.sigma
        self.C = max(0, self.C + abs(z) - self.k)
        return self.C > self.h

# =============================================================================
# TEST PRINCIPAL
# =============================================================================

print("\n" + "="*70)
print("Generando trayectorias H(z)")
print("="*70)

# Alta resolucion alrededor de recombinacion
z_all = np.concatenate([
    np.linspace(3000, 1500, 100),
    np.linspace(1499, 1000, 200),  # Alta resolucion cerca de z_rec
    np.linspace(999, 100, 100),
    np.linspace(99, 0, 50)
])

# Calcular H y derivadas
H_lcdm_arr = np.array([H_lcdm(z) for z in z_all])
H_octh_arr = np.array([H_octh(z) for z in z_all])
psi_arr = np.array([psi_octh(z) for z in z_all])
q_lcdm_arr = np.array([q_deceleration(z, 'lcdm') for z in z_all])
q_octh_arr = np.array([q_deceleration(z, 'octh') for z in z_all])

# Diferencia relativa
delta_H = (H_octh_arr - H_lcdm_arr) / H_lcdm_arr * 100

print(f"\nRango de z: [{z_all[-1]:.0f}, {z_all[0]:.0f}]")
print(f"Puntos: {len(z_all)}")

# =============================================================================
# ANALISIS TDA: Embedding de H(z), dH/dz
# =============================================================================

print("\n" + "="*70)
print("Analisis TDA de trayectoria (H, dH/dz)")
print("="*70)

# Crear embedding 2D: (H, dH/dz)
dH_lcdm = np.gradient(H_lcdm_arr, z_all)
dH_octh = np.gradient(H_octh_arr, z_all)

# Normalizar para TDA
def normalize(arr):
    return (arr - np.mean(arr)) / (np.std(arr) + 1e-10)

# Trayectoria LCDM en espacio de fase
traj_lcdm = np.column_stack([
    normalize(H_lcdm_arr),
    normalize(dH_lcdm)
])

# Trayectoria OCTH en espacio de fase
traj_octh = np.column_stack([
    normalize(H_octh_arr),
    normalize(dH_octh)
])

# Diferencia de trayectorias
traj_diff = traj_octh - traj_lcdm

print("\nCalculando entropia topologica en ventanas deslizantes...")

# Ventana deslizante para calcular S_H1(z)
window_size = 30
step = 5

S_H1_lcdm = []
S_H1_octh = []
S_H1_diff = []
z_windows = []

for i in range(0, len(z_all) - window_size, step):
    window_lcdm = traj_lcdm[i:i+window_size]
    window_octh = traj_octh[i:i+window_size]
    window_diff = traj_diff[i:i+window_size]

    s_lcdm = persistence_entropy_h1(window_lcdm)
    s_octh = persistence_entropy_h1(window_octh)
    s_diff = persistence_entropy_h1(window_diff)

    S_H1_lcdm.append(s_lcdm)
    S_H1_octh.append(s_octh)
    S_H1_diff.append(s_diff)
    z_windows.append(z_all[i + window_size // 2])

S_H1_lcdm = np.array(S_H1_lcdm)
S_H1_octh = np.array(S_H1_octh)
S_H1_diff = np.array(S_H1_diff)
z_windows = np.array(z_windows)

print(f"\nVentanas calculadas: {len(z_windows)}")

# =============================================================================
# DETECCION CUSUM
# =============================================================================

print("\n" + "="*70)
print("Deteccion CUSUM de transicion topologica")
print("="*70)

# Usar z > 2000 como baseline (bien antes de recombinacion)
baseline_mask = z_windows > 2000
baseline_idx = np.sum(baseline_mask)

if baseline_idx > 5:
    baseline_S = S_H1_diff[baseline_mask]

    detector = CUSUM(k=0.3, h=3.0, sigma_min=0.01)
    detector.calibrate(baseline_S)

    print(f"\nBaseline (z > 2000):")
    print(f"  mu = {detector.mu:.4f}")
    print(f"  sigma = {detector.sigma:.4f}")

    z_detected = None
    for i in range(baseline_idx, len(z_windows)):
        if detector.update(S_H1_diff[i]):
            z_detected = z_windows[i]
            break

    print(f"\nResultados:")
    if z_detected:
        print(f"  z_CUSUM (topologia) = {z_detected:.1f}")
        print(f"  z_rec (fisica)      = {Z_REC:.1f}")

        if z_detected > Z_REC:
            print(f"\n  *** TOPOLOGIA PRECEDE FISICA por Delta_z = {z_detected - Z_REC:.1f} ***")
            PIT_CONFIRMED = True
        else:
            print(f"\n  Topologia detectada despues de fisica")
            PIT_CONFIRMED = False
    else:
        print("  No se detecto transicion CUSUM")
        PIT_CONFIRMED = False
else:
    print("  Baseline insuficiente")
    PIT_CONFIRMED = False

# =============================================================================
# ANALISIS DE MINIMOS
# =============================================================================

print("\n" + "="*70)
print("Analisis de extremos")
print("="*70)

# Minimo de Psi
idx_psi_min = np.argmin(psi_arr)
z_psi_min = z_all[idx_psi_min]

# Maximo de delta_H
idx_delta_max = np.argmax(np.abs(delta_H))
z_delta_max = z_all[idx_delta_max]

# Maximo de S_H1_diff (maxima "actividad topologica")
if len(S_H1_diff) > 0:
    idx_S_max = np.argmax(S_H1_diff)
    z_S_max = z_windows[idx_S_max]
else:
    z_S_max = None

print(f"\nPuntos criticos:")
print(f"  z(Psi_min)     = {z_psi_min:.1f}")
print(f"  z(Delta_H_max) = {z_delta_max:.1f}")
if z_S_max:
    print(f"  z(S_H1_max)    = {z_S_max:.1f}")
print(f"  z_recombinacion = {Z_REC:.1f}")

# Verificar si S_H1 pico ANTES de recombinacion
if z_S_max and z_S_max > Z_REC:
    print(f"\n  *** S_H1 PICO ANTES DE RECOMBINACION ***")
    print(f"  La actividad topologica precede la transicion fisica!")
    PEAK_PRECEDES = True
else:
    PEAK_PRECEDES = False

# =============================================================================
# PERFIL DETALLADO
# =============================================================================

print("\n" + "="*70)
print("Perfil S_H1(z) alrededor de recombinacion")
print("="*70)

print("\n{:>10} {:>12} {:>12} {:>12} {:>10}".format(
    "z", "S_H1(LCDM)", "S_H1(OCTH)", "S_H1(diff)", "Delta_H%"))
print("-"*60)

# Filtrar alrededor de recombinacion
mask_rec = (z_windows > 800) & (z_windows < 1500)
indices_rec = np.where(mask_rec)[0]

for i in indices_rec[::3]:  # Cada 3er punto
    marker = ""
    if abs(z_windows[i] - Z_REC) < 50:
        marker = " <-- z_rec"
    if z_S_max and abs(z_windows[i] - z_S_max) < 30:
        marker = " <-- S_max"

    print(f"{z_windows[i]:>10.1f} {S_H1_lcdm[i]:>12.4f} {S_H1_octh[i]:>12.4f} {S_H1_diff[i]:>12.4f} {delta_H[i*step]:>10.3f}{marker}")

# =============================================================================
# CONCLUSION
# =============================================================================

print("\n" + "="*70)
print("CONCLUSIONES")
print("="*70)

print(f"""
RESULTADOS DEL TEST PIT v2:

1. TRAYECTORIA EN ESPACIO DE FASE:
   - LCDM: trayectoria suave
   - OCTH: perturbacion gaussiana en z ~ {Z_TOPO}

2. ENTROPIA TOPOLOGICA S_H1:
   - Mide "complejidad" de la trayectoria
   - S_H1(diff) captura desviacion OCTH vs LCDM

3. DETECCION CUSUM:
   - {'z_CUSUM = ' + f'{z_detected:.1f}' if z_detected else 'No detectado claramente'}
   - {'PRECEDE z_rec!' if PIT_CONFIRMED else 'Test no concluyente'}

4. PICO DE ACTIVIDAD TOPOLOGICA:
   - {'z_peak = ' + f'{z_S_max:.1f}' if z_S_max else 'No identificado'}
   - {'ANTES de recombinacion!' if PEAK_PRECEDES else 'Patron a investigar'}

5. INTERPRETACION PIT:
""")

if PIT_CONFIRMED or PEAK_PRECEDES:
    print("""   *** EVIDENCIA A FAVOR DE PIT ***

   La desviacion topologica (S_H1 de diferencia OCTH-LCDM)
   muestra actividad ANTES de z_recombinacion.

   Esto sugiere que el cambio topologico (Mobius/Psi)
   PRECEDE y posiblemente CAUSA la transicion fisica.
""")
else:
    print("""   Test no concluyente con esta metodologia.

   Posibles mejoras:
   - Usar ripser exacto (no aproximacion Euler)
   - Mayor resolucion en z
   - Embedding de mayor dimension
   - Datos observacionales reales
""")

print("="*70)

# Guardar resultados
import json
output = {
    'z_windows': z_windows.tolist(),
    'S_H1_lcdm': S_H1_lcdm.tolist(),
    'S_H1_octh': S_H1_octh.tolist(),
    'S_H1_diff': S_H1_diff.tolist(),
    'z_detected': z_detected if 'z_detected' in dir() else None,
    'z_S_max': z_S_max,
    'z_rec': Z_REC,
    'pit_confirmed': PIT_CONFIRMED if 'PIT_CONFIRMED' in dir() else False,
    'peak_precedes': PEAK_PRECEDES if 'PEAK_PRECEDES' in dir() else False
}

with open('PIT_hubble_dynamics_results.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nResultados guardados en PIT_hubble_dynamics_results.json")
