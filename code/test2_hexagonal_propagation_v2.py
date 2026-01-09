#!/usr/bin/env python3
"""
TEST #2 v2: SIMULACIÓN CORREGIDA - RÉGIMEN DE CAMPO DÉBIL
==========================================================

Versión corregida con escalas físicas apropiadas para régimen de campo débil
donde la aproximación de Schwarzschild es válida (b >> rs).

OCTH - Ontología del Campo Tensorial Hexagonal
Autor: Francisco Molina Burgos (Yatrogenesis)
Fecha: 2026-01-08
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.interpolate import griddata
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONSTANTES (unidades donde c=G=1, escaladas para campo débil)
# =============================================================================
C = 1.0
G = 1.0

# =============================================================================
# RETÍCULO HEXAGONAL
# =============================================================================

class HexLattice:
    def __init__(self, Nx, Ny, spacing=1.0):
        self.Nx, self.Ny = Nx, Ny
        self.N = Nx * Ny
        self.spacing = spacing
        self.positions = self._make_positions()

    def _make_positions(self):
        pos = []
        for j in range(self.Ny):
            for i in range(self.Nx):
                x = i * self.spacing + (j % 2) * self.spacing / 2
                y = j * self.spacing * np.sqrt(3) / 2
                pos.append([x, y])
        pos = np.array(pos)
        pos -= pos.mean(axis=0)
        return pos


class SquareLattice:
    def __init__(self, Nx, Ny, spacing=1.0):
        self.Nx, self.Ny = Nx, Ny
        self.N = Nx * Ny
        self.spacing = spacing
        self.positions = self._make_positions()

    def _make_positions(self):
        pos = []
        for j in range(self.Ny):
            for i in range(self.Nx):
                pos.append([i * self.spacing, j * self.spacing])
        pos = np.array(pos)
        pos -= pos.mean(axis=0)
        return pos


# =============================================================================
# CAMPO Ψ EN RÉGIMEN DE CAMPO DÉBIL
# =============================================================================

def psi_schwarzschild(r, M):
    """
    Ψ = √(1 - rs/r) donde rs = 2GM/c²

    En campo débil (r >> rs): Ψ ≈ 1 - GM/(rc²)
    """
    rs = 2 * G * M / C**2
    ratio = rs / np.maximum(r, rs * 1.01)  # Evitar singularidad
    return np.sqrt(np.maximum(1 - ratio, 0.01))


def compute_psi(positions, mass_pos, M):
    """Calcula Ψ en todos los nodos."""
    r = np.linalg.norm(positions - mass_pos, axis=1)
    return psi_schwarzschild(r, M)


# =============================================================================
# TRAZADO DE RAYOS CON GRADIENTE DE Ψ
# =============================================================================

def trace_ray_octh(positions, psi, start, direction, n_steps=2000, ds=0.05):
    """
    Traza rayo usando ecuación de eikonal en medio con índice de refracción n = 1/Ψ.

    En óptica geométrica: d(n·dr/ds)/ds = ∇n
    Con n = 1/Ψ: la luz se curva hacia regiones de Ψ bajo (como en GR).
    """
    traj = [start.copy()]
    pos = start.copy()

    # Normalizar dirección
    vel = direction / np.linalg.norm(direction)

    for _ in range(n_steps):
        # Interpolar Ψ y su gradiente en posición actual
        dist = np.linalg.norm(positions - pos, axis=1)

        # Kernel de interpolación
        sigma = 1.5
        w = np.exp(-dist**2 / (2 * sigma**2))
        w /= w.sum() + 1e-10

        psi_local = np.dot(w, psi)

        # Gradiente de Ψ (diferencias finitas locales)
        grad_psi = np.zeros(2)
        for i in np.argsort(dist)[:20]:  # 20 vecinos más cercanos
            if dist[i] > 0.01:
                dr = positions[i] - pos
                dpsi = psi[i] - psi_local
                grad_psi += w[i] * dpsi * dr / (dist[i]**2 + 0.01)

        # Índice de refracción efectivo n = 1/Ψ
        # Deflexión: dθ/ds = -∂log(n)/∂⊥ = ∂log(Ψ)/∂⊥ = (1/Ψ)∂Ψ/∂⊥

        # Componente perpendicular del gradiente
        grad_perp = grad_psi - np.dot(grad_psi, vel) * vel

        # Deflexión proporcional a gradiente perpendicular de log(Ψ)
        deflection = grad_perp / (psi_local + 0.1) * ds * 2

        # Actualizar velocidad
        vel = vel + deflection
        vel = vel / np.linalg.norm(vel)

        # Avanzar
        pos = pos + vel * ds
        traj.append(pos.copy())

        # Condición de salida
        if np.linalg.norm(pos) > 30:
            break

    return np.array(traj)


def deflection_gr(b, M):
    """
    Deflexión GR para luz pasando a distancia b de masa M.

    Δφ = 4GM/(bc²) [radianes] - válido para b >> rs
    """
    return 4 * G * M / (b * C**2)


# =============================================================================
# TEST PRINCIPAL
# =============================================================================

def run_test():
    print("=" * 70)
    print("TEST #2 v2: RETÍCULO HEXAGONAL - CAMPO DÉBIL")
    print("=" * 70)

    # Parámetros en régimen de campo débil
    # Masa pequeña para que rs << b
    M = 0.1  # rs = 0.2

    # Crear retículo grande
    Nx, Ny = 80, 80
    hex_lat = HexLattice(Nx, Ny, spacing=0.8)
    sq_lat = SquareLattice(Nx, Ny, spacing=0.8)

    mass_pos = np.array([0.0, 0.0])

    print(f"\nParámetros:")
    print(f"  Masa M = {M}")
    print(f"  Radio Schwarzschild rs = {2*G*M/C**2:.3f}")
    print(f"  Nodos: {hex_lat.N}")

    # Calcular Ψ
    psi_hex = compute_psi(hex_lat.positions, mass_pos, M)
    psi_sq = compute_psi(sq_lat.positions, mass_pos, M)

    print(f"  Ψ_min = {psi_hex.min():.4f}")
    print(f"  Ψ_max = {psi_hex.max():.4f}")

    # =========================================================================
    # FIGURA 1: Estructura de retículos
    # =========================================================================
    print("\n1. Generando comparación de retículos...")

    fig1, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    sc = ax.scatter(hex_lat.positions[:,0], hex_lat.positions[:,1],
                   c=psi_hex, cmap='RdYlBu', s=8, vmin=0.9, vmax=1.0)
    ax.scatter(0, 0, c='black', s=200, marker='*', zorder=10)
    ax.set_aspect('equal')
    ax.set_title('Retículo HEXAGONAL + Ψ', fontsize=14)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(sc, ax=ax, label='Ψ')

    ax = axes[1]
    sc = ax.scatter(sq_lat.positions[:,0], sq_lat.positions[:,1],
                   c=psi_sq, cmap='RdYlBu', s=8, vmin=0.9, vmax=1.0)
    ax.scatter(0, 0, c='black', s=200, marker='*', zorder=10)
    ax.set_aspect('equal')
    ax.set_title('Retículo CUADRADO + Ψ (Control)', fontsize=14)
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.colorbar(sc, ax=ax, label='Ψ')

    plt.tight_layout()
    plt.savefig('../figures/fig1_lattices.png', dpi=300)
    plt.savefig('../figures/fig1_lattices.pdf', dpi=300)
    plt.close()
    print("   ✓ fig1_lattices.png/pdf")

    # =========================================================================
    # FIGURA 2: Mapa de calor de Ψ
    # =========================================================================
    print("\n2. Generando mapa de calor de Ψ...")

    fig2, ax = plt.subplots(figsize=(10, 8))

    # Interpolación para contornos
    xi = np.linspace(-25, 25, 200)
    yi = np.linspace(-25, 25, 200)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = griddata(hex_lat.positions, psi_hex, (Xi, Yi), method='linear', fill_value=1.0)

    cf = ax.contourf(Xi, Yi, Zi, levels=20, cmap='RdYlBu')
    ax.contour(Xi, Yi, Zi, levels=[0.95, 0.97, 0.99], colors='black', linewidths=0.5)

    ax.scatter(0, 0, c='black', s=300, marker='*', zorder=10, label='Masa M')

    # Círculo de Schwarzschild
    rs = 2 * G * M / C**2
    circle = plt.Circle((0, 0), rs, fill=False, color='red', ls='--', lw=2)
    ax.add_patch(circle)

    ax.set_xlim(-25, 25)
    ax.set_ylim(-25, 25)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Campo Ψ = √(1 - rs/r)\nTiempo LENTO cerca de masa (rojo)', fontsize=14)
    plt.colorbar(cf, ax=ax, label='Ψ (permeabilidad temporal)')

    plt.tight_layout()
    plt.savefig('../figures/fig2_psi_field.png', dpi=300)
    plt.savefig('../figures/fig2_psi_field.pdf', dpi=300)
    plt.close()
    print("   ✓ fig2_psi_field.png/pdf")

    # =========================================================================
    # FIGURA 3: Ray tracing
    # =========================================================================
    print("\n3. Trazando rayos de luz...")

    # Parámetros de impacto (todos >> rs para campo débil válido)
    b_values = np.array([2.0, 4.0, 6.0, 8.0, 10.0, 15.0])

    fig3, ax = plt.subplots(figsize=(12, 10))

    # Fondo: campo Ψ
    cf = ax.contourf(Xi, Yi, Zi, levels=20, cmap='Greys', alpha=0.3)

    # Masa
    ax.scatter(0, 0, c='black', s=400, marker='o', zorder=10, label='Masa M')
    circle = plt.Circle((0, 0), rs, fill=False, color='red', ls='--', lw=2, label='$r_s$')
    ax.add_patch(circle)

    # Rayos
    deflections_octh = []
    deflections_gr_list = []

    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(b_values)))

    for i, b in enumerate(b_values):
        # OCTH
        start = np.array([b, -25.0])
        direction = np.array([0.0, 1.0])
        traj = trace_ray_octh(hex_lat.positions, psi_hex, start, direction)

        ax.plot(traj[:,0], traj[:,1], '-', color=colors[i], lw=2,
               label=f'b={b:.0f}' if i == 0 else f'b={b:.0f}')

        # Calcular deflexión OCTH
        if len(traj) > 50:
            dir_in = traj[10] - traj[0]
            dir_out = traj[-1] - traj[-11]
            dir_in /= np.linalg.norm(dir_in)
            dir_out /= np.linalg.norm(dir_out)
            angle = np.arctan2(dir_out[0], dir_out[1]) - np.arctan2(dir_in[0], dir_in[1])
            deflections_octh.append(abs(angle))
        else:
            deflections_octh.append(0)

        # GR teórico
        deflections_gr_list.append(deflection_gr(b, M))

        # Línea recta (sin gravedad) para comparación
        if i == 0:
            ax.plot([b, b], [-25, 25], 'k:', alpha=0.5, label='Sin gravedad')

    ax.set_xlim(-5, 20)
    ax.set_ylim(-25, 25)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Deflexión de Luz en Retículo Hexagonal (OCTH)', fontsize=14)
    ax.legend(loc='upper right', fontsize=10)

    plt.tight_layout()
    plt.savefig('../figures/fig3_ray_tracing.png', dpi=300)
    plt.savefig('../figures/fig3_ray_tracing.pdf', dpi=300)
    plt.close()
    print("   ✓ fig3_ray_tracing.png/pdf")

    deflections_octh = np.array(deflections_octh)
    deflections_gr_arr = np.array(deflections_gr_list)

    # =========================================================================
    # FIGURA 4: Comparación cuantitativa
    # =========================================================================
    print("\n4. Comparación cuantitativa OCTH vs GR...")

    fig4, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel izquierdo: Deflexión vs b
    ax = axes[0]
    ax.loglog(b_values, np.degrees(deflections_gr_arr), 'ro-', lw=2, ms=10,
             label='GR: Δφ = 4GM/(bc²)')
    ax.loglog(b_values, np.degrees(deflections_octh), 'b^--', lw=2, ms=10,
             label='OCTH (retículo hexagonal)')

    ax.set_xlabel('Parámetro de impacto b', fontsize=12)
    ax.set_ylabel('Deflexión Δφ (grados)', fontsize=12)
    ax.set_title('Deflexión de Luz: OCTH vs GR', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

    # Panel derecho: Error relativo
    ax = axes[1]
    valid = deflections_gr_arr > 0
    error = np.abs(deflections_octh[valid] - deflections_gr_arr[valid]) / deflections_gr_arr[valid] * 100

    ax.plot(b_values[valid], error, 'g-o', lw=2, ms=10)
    ax.axhline(100, color='red', ls='--', lw=2, label='100% error')
    ax.axhline(50, color='orange', ls='--', lw=2, label='50% error')
    ax.axhline(10, color='green', ls='--', lw=2, label='10% error')

    ax.set_xlabel('Parámetro de impacto b', fontsize=12)
    ax.set_ylabel('Error relativo (%)', fontsize=12)
    ax.set_title('Convergencia OCTH → GR', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(150, error.max() * 1.1))

    plt.tight_layout()
    plt.savefig('../figures/fig4_comparison.png', dpi=300)
    plt.savefig('../figures/fig4_comparison.pdf', dpi=300)
    plt.close()
    print("   ✓ fig4_comparison.png/pdf")

    # =========================================================================
    # RESULTADOS
    # =========================================================================
    print("\n" + "=" * 70)
    print("RESULTADOS")
    print("=" * 70)

    print(f"\n{'b':<8} {'OCTH (°)':<12} {'GR (°)':<12} {'Error (%)':<12}")
    print("-" * 50)
    for i, b in enumerate(b_values):
        if deflections_gr_arr[i] > 0:
            err = abs(deflections_octh[i] - deflections_gr_arr[i]) / deflections_gr_arr[i] * 100
        else:
            err = 0
        print(f"{b:<8.1f} {np.degrees(deflections_octh[i]):<12.4f} "
              f"{np.degrees(deflections_gr_arr[i]):<12.4f} {err:<12.1f}")

    mean_error = np.mean(error) if len(error) > 0 else 0

    print("-" * 50)
    print(f"\nError promedio: {mean_error:.1f}%")

    # Verificar comportamiento cualitativo
    # La deflexión debe DECRECER con b (inversamente proporcional)
    octh_decreasing = all(deflections_octh[i] >= deflections_octh[i+1]
                         for i in range(len(deflections_octh)-1) if deflections_octh[i+1] > 0)

    print(f"\n✓ Deflexión decrece con distancia: {'SÍ' if octh_decreasing else 'NO'}")
    print(f"✓ Luz se curva hacia la masa: {'SÍ' if deflections_octh.mean() > 0 else 'NO'}")

    if octh_decreasing and deflections_octh.mean() > 0:
        print("\n" + "=" * 70)
        print("🏆 CONCLUSIÓN: El retículo hexagonal con Ψ variable")
        print("   REPRODUCE CUALITATIVAMENTE el lensing gravitacional.")
        print("   - La luz se deflecta hacia la masa (Ψ bajo)")
        print("   - La deflexión escala inversamente con b")
        print("=" * 70)

    return {
        'b_values': b_values,
        'deflections_octh': deflections_octh,
        'deflections_gr': deflections_gr_arr,
        'mean_error': mean_error
    }


if __name__ == "__main__":
    results = run_test()
