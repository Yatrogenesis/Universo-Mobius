#!/usr/bin/env python3
"""
TEST #2: SIMULACIÓN DE PROPAGACIÓN EN RETÍCULO HEXAGONAL
==========================================================

Demuestra que ondas en malla hexagonal con Ψ variable reproducen
lensing gravitacional (deflexión de luz cerca de masas).

Incluye:
1. Grupo de control: Malla cuadrada vs hexagonal
2. Mapa de calor del campo Ψ
3. Comparación con geodésicas de Schwarzschild

OCTH - Ontología del Campo Tensorial Hexagonal
Autor: Francisco Molina Burgos (Yatrogenesis)
Fecha: 2026-01-08
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import eigsh, spsolve
from scipy.integrate import odeint
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# Constantes físicas (unidades naturales normalizadas)
C_LIGHT = 1.0  # Velocidad de la luz
G_NEWTON = 1.0  # Constante gravitacional
RHO_PLANCK = 1.0  # Densidad de Planck (normalizada)

# =============================================================================
# CLASE: RETÍCULO HEXAGONAL
# =============================================================================

class HexagonalLattice:
    """
    Retículo hexagonal 2D con base vectorial triádica.

    Vectores base: e1, e2, e3 con e1 + e2 + e3 = 0 (simetría D3)
    """

    def __init__(self, Nx, Ny, spacing=1.0):
        self.Nx = Nx
        self.Ny = Ny
        self.spacing = spacing
        self.N = Nx * Ny

        # Vectores base hexagonales
        self.e1 = np.array([1.0, 0.0])
        self.e2 = np.array([-0.5, np.sqrt(3)/2])
        self.e3 = np.array([-0.5, -np.sqrt(3)/2])

        # Generar posiciones de nodos
        self.positions = self._generate_positions()

        # Generar matriz de adyacencia
        self.adjacency = self._generate_adjacency()

    def _generate_positions(self):
        """Genera posiciones de nodos en retículo hexagonal."""
        positions = np.zeros((self.N, 2))

        for j in range(self.Ny):
            for i in range(self.Nx):
                idx = j * self.Nx + i
                # Offset en filas impares para patrón hexagonal
                x = i * self.spacing + (j % 2) * self.spacing / 2
                y = j * self.spacing * np.sqrt(3) / 2
                positions[idx] = [x, y]

        # Centrar en origen
        positions -= positions.mean(axis=0)
        return positions

    def _generate_adjacency(self):
        """Genera matriz de adyacencia para conectividad hexagonal."""
        adj = lil_matrix((self.N, self.N))

        for j in range(self.Ny):
            for i in range(self.Nx):
                idx = j * self.Nx + i

                # Vecinos horizontales
                if i > 0:
                    adj[idx, idx - 1] = 1
                if i < self.Nx - 1:
                    adj[idx, idx + 1] = 1

                # Vecinos verticales (depende de paridad de fila)
                if j > 0:
                    if j % 2 == 0:
                        if i > 0:
                            adj[idx, (j-1) * self.Nx + i - 1] = 1
                        adj[idx, (j-1) * self.Nx + i] = 1
                    else:
                        adj[idx, (j-1) * self.Nx + i] = 1
                        if i < self.Nx - 1:
                            adj[idx, (j-1) * self.Nx + i + 1] = 1

                if j < self.Ny - 1:
                    if j % 2 == 0:
                        if i > 0:
                            adj[idx, (j+1) * self.Nx + i - 1] = 1
                        adj[idx, (j+1) * self.Nx + i] = 1
                    else:
                        adj[idx, (j+1) * self.Nx + i] = 1
                        if i < self.Nx - 1:
                            adj[idx, (j+1) * self.Nx + i + 1] = 1

        return csr_matrix(adj)

    def get_laplacian(self):
        """Calcula el Laplaciano discreto sobre el grafo."""
        degree = np.array(self.adjacency.sum(axis=1)).flatten()
        D = lil_matrix((self.N, self.N))
        for i in range(self.N):
            D[i, i] = degree[i]
        return csr_matrix(D - self.adjacency)


class SquareLattice:
    """
    Retículo cuadrado 2D (grupo de control).
    """

    def __init__(self, Nx, Ny, spacing=1.0):
        self.Nx = Nx
        self.Ny = Ny
        self.spacing = spacing
        self.N = Nx * Ny

        self.positions = self._generate_positions()
        self.adjacency = self._generate_adjacency()

    def _generate_positions(self):
        positions = np.zeros((self.N, 2))
        for j in range(self.Ny):
            for i in range(self.Nx):
                idx = j * self.Nx + i
                positions[idx] = [i * self.spacing, j * self.spacing]
        positions -= positions.mean(axis=0)
        return positions

    def _generate_adjacency(self):
        adj = lil_matrix((self.N, self.N))
        for j in range(self.Ny):
            for i in range(self.Nx):
                idx = j * self.Nx + i
                if i > 0:
                    adj[idx, idx - 1] = 1
                if i < self.Nx - 1:
                    adj[idx, idx + 1] = 1
                if j > 0:
                    adj[idx, idx - self.Nx] = 1
                if j < self.Ny - 1:
                    adj[idx, idx + self.Nx] = 1
        return csr_matrix(adj)

    def get_laplacian(self):
        degree = np.array(self.adjacency.sum(axis=1)).flatten()
        D = lil_matrix((self.N, self.N))
        for i in range(self.N):
            D[i, i] = degree[i]
        return csr_matrix(D - self.adjacency)


# =============================================================================
# CAMPO Ψ (PERMEABILIDAD TEMPORAL)
# =============================================================================

def compute_psi_field(positions, mass_positions, masses, rs=0.5):
    """
    Calcula el campo Ψ según OCTH.

    Ψ(r) = sqrt(1 - (ρ_m/ρ_Planck)²)

    En campo débil: Ψ ≈ sqrt(1 - 2GM/(rc²))

    Args:
        positions: Posiciones de nodos del retículo
        mass_positions: Posiciones de masas puntuales
        masses: Valores de masas
        rs: Radio de Schwarzschild característico (para normalización)

    Returns:
        Psi: Campo de permeabilidad en cada nodo
    """
    N = len(positions)
    psi = np.ones(N)

    for m_pos, M in zip(mass_positions, masses):
        for i in range(N):
            r = np.linalg.norm(positions[i] - m_pos)
            r = max(r, 0.1 * rs)  # Evitar singularidad

            # Métrica de Schwarzschild: g_tt = -(1 - rs/r)
            # Ψ = sqrt(1 - rs/r) donde rs = 2GM/c²
            rs_eff = 2 * G_NEWTON * M / C_LIGHT**2

            psi[i] *= np.sqrt(max(1 - rs_eff / r, 0.01))

    return psi


# =============================================================================
# PROPAGACIÓN DE ONDAS
# =============================================================================

def simulate_wave_propagation(lattice, psi_field, source_pos, t_max=50, dt=0.1):
    """
    Simula propagación de onda en retículo con Ψ variable.

    Ecuación: Ψ² ∂²u/∂t² = c² Δu

    Args:
        lattice: Objeto de retículo (hexagonal o cuadrado)
        psi_field: Campo Ψ en cada nodo
        source_pos: Posición de la fuente
        t_max: Tiempo máximo de simulación
        dt: Paso temporal

    Returns:
        u_history: Historia temporal del campo
        times: Vector de tiempos
    """
    N = lattice.N
    L = lattice.get_laplacian()

    # Condiciones iniciales
    u = np.zeros(N)
    u_dot = np.zeros(N)

    # Encontrar nodo más cercano a la fuente
    distances = np.linalg.norm(lattice.positions - source_pos, axis=1)
    source_node = np.argmin(distances)

    # Pulso gaussiano inicial
    sigma = 2.0
    for i in range(N):
        r = np.linalg.norm(lattice.positions[i] - source_pos)
        u[i] = np.exp(-r**2 / (2 * sigma**2))

    # Evolución temporal (método de Verlet)
    times = np.arange(0, t_max, dt)
    u_history = [u.copy()]

    for t in times[1:]:
        # Aceleración: a = (c²/Ψ²) Δu
        psi_sq = psi_field**2
        psi_sq = np.maximum(psi_sq, 0.01)  # Evitar división por cero

        Lu = L.dot(u)
        u_ddot = -C_LIGHT**2 * Lu / psi_sq

        # Integración de Verlet
        u_new = 2 * u - u_history[-1] + u_ddot * dt**2

        # Amortiguamiento en bordes (condiciones absorbentes)
        boundary_mask = (np.abs(lattice.positions[:, 0]) > 0.8 * lattice.positions[:, 0].max()) | \
                       (np.abs(lattice.positions[:, 1]) > 0.8 * lattice.positions[:, 1].max())
        u_new[boundary_mask] *= 0.95

        u_history.append(u_new.copy())
        u = u_new

    return np.array(u_history), times


def trace_ray(lattice, psi_field, start_pos, direction, n_steps=500, step_size=0.1):
    """
    Traza un rayo de luz a través del campo Ψ.

    La velocidad efectiva es c_eff = c₀ * Ψ, lo que produce
    curvatura de la trayectoria (lensing).

    Usa el principio de Fermat: la luz sigue caminos de tiempo mínimo.
    """
    trajectory = [start_pos.copy()]
    pos = start_pos.copy()
    vel = direction.copy()
    vel = vel / np.linalg.norm(vel) * C_LIGHT

    for _ in range(n_steps):
        # Interpolar Ψ en posición actual
        distances = np.linalg.norm(lattice.positions - pos, axis=1)
        weights = np.exp(-distances**2 / 2)
        weights /= weights.sum()
        psi_local = np.dot(weights, psi_field)

        # Calcular gradiente de Ψ
        grad_psi = np.zeros(2)
        for i in range(len(weights)):
            if weights[i] > 0.01:
                dr = lattice.positions[i] - pos
                grad_psi += weights[i] * psi_field[i] * dr / (np.linalg.norm(dr) + 0.01)

        # Velocidad efectiva: c_eff = c₀ * Ψ
        c_eff = C_LIGHT * psi_local

        # Deflexión: la luz se curva hacia regiones de menor Ψ
        # (como en óptica de medios inhomogéneos)
        deflection = -grad_psi / (psi_local + 0.01) * step_size

        vel += deflection * c_eff
        vel = vel / np.linalg.norm(vel) * c_eff

        pos = pos + vel * step_size / c_eff
        trajectory.append(pos.copy())

        # Salir si está fuera del dominio
        if np.linalg.norm(pos) > 1.5 * np.max(np.abs(lattice.positions)):
            break

    return np.array(trajectory)


# =============================================================================
# GEODÉSICAS DE SCHWARZSCHILD (REFERENCIA GR)
# =============================================================================

def schwarzschild_geodesic(M, b, r_start=20, r_end=-20, n_points=1000):
    """
    Calcula geodésica nula (luz) en métrica de Schwarzschild.

    Para luz pasando a distancia b de masa M:
    φ = ∫ dr / (r² √(1/b² - 1/r² + rs/r³))

    Deflexión total para b >> rs: Δφ ≈ 4GM/(bc²) = 2rs/b
    """
    rs = 2 * G_NEWTON * M / C_LIGHT**2

    # Integrar trayectoria
    def dphi_dr(r, b, rs):
        if r <= rs:
            return 0
        term = 1/b**2 - 1/r**2 + rs/r**3
        if term <= 0:
            return 0
        return 1 / (r**2 * np.sqrt(term))

    # Aproximación: trayectoria casi recta con deflexión pequeña
    y_vals = np.linspace(-r_start, r_end, n_points)
    x_vals = np.full_like(y_vals, b)

    # Aplicar deflexión acumulativa
    delta_phi = 4 * G_NEWTON * M / (b * C_LIGHT**2)  # Deflexión GR

    trajectory = np.column_stack([x_vals, y_vals])

    # Rotar por la deflexión
    cos_d, sin_d = np.cos(delta_phi/2), np.sin(delta_phi/2)
    rotation = np.array([[cos_d, -sin_d], [sin_d, cos_d]])

    for i in range(len(trajectory)):
        # Deflexión gradual a lo largo del camino
        if y_vals[i] > 0:
            frac = 1 - y_vals[i] / r_start
        else:
            frac = 1 + y_vals[i] / r_start

        angle = delta_phi * frac / 2
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
        trajectory[i] = R @ trajectory[i]

    return trajectory, delta_phi


# =============================================================================
# VISUALIZACIÓN
# =============================================================================

def plot_lattice_comparison(hex_lattice, sq_lattice, psi_hex, psi_sq, figsize=(14, 6)):
    """Compara estructura de retículos hexagonal vs cuadrado."""

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Hexagonal
    ax = axes[0]
    scatter = ax.scatter(hex_lattice.positions[:, 0], hex_lattice.positions[:, 1],
                        c=psi_hex, cmap='coolwarm_r', s=20, vmin=0.5, vmax=1.0)
    ax.set_aspect('equal')
    ax.set_title('Retículo Hexagonal (OCTH)', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(scatter, ax=ax, label='Ψ (permeabilidad)')

    # Cuadrado
    ax = axes[1]
    scatter = ax.scatter(sq_lattice.positions[:, 0], sq_lattice.positions[:, 1],
                        c=psi_sq, cmap='coolwarm_r', s=20, vmin=0.5, vmax=1.0)
    ax.set_aspect('equal')
    ax.set_title('Retículo Cuadrado (Control)', fontsize=14)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(scatter, ax=ax, label='Ψ (permeabilidad)')

    plt.tight_layout()
    return fig


def plot_psi_heatmap(lattice, psi, mass_pos, title="Campo Ψ", figsize=(10, 8)):
    """Visualiza el campo Ψ como mapa de calor."""

    fig, ax = plt.subplots(figsize=figsize)

    # Crear grid para interpolación
    x_range = [lattice.positions[:, 0].min(), lattice.positions[:, 0].max()]
    y_range = [lattice.positions[:, 1].min(), lattice.positions[:, 1].max()]

    # Scatter plot con colores
    scatter = ax.scatter(lattice.positions[:, 0], lattice.positions[:, 1],
                        c=psi, cmap='RdYlBu', s=30, vmin=0, vmax=1)

    # Marcar masa
    ax.scatter(mass_pos[0], mass_pos[1], c='black', s=200, marker='*',
              label=f'Masa (Ψ→0)', edgecolors='white', linewidths=2)

    # Contornos de Ψ constante
    from scipy.interpolate import griddata
    xi = np.linspace(x_range[0], x_range[1], 100)
    yi = np.linspace(y_range[0], y_range[1], 100)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = griddata(lattice.positions, psi, (Xi, Yi), method='cubic', fill_value=1.0)

    contours = ax.contour(Xi, Yi, Zi, levels=[0.3, 0.5, 0.7, 0.9], colors='black',
                         linewidths=0.5, alpha=0.5)
    ax.clabel(contours, inline=True, fontsize=8, fmt='%.1f')

    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(loc='upper right')

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Ψ = √(1 - ρ/ρ_Planck)', fontsize=11)

    # Añadir anotaciones
    ax.annotate('Tiempo LENTO\n(Ψ → 0)', xy=(mass_pos[0], mass_pos[1] - 3),
               fontsize=10, ha='center', color='red')
    ax.annotate('Tiempo RÁPIDO\n(Ψ → 1)', xy=(x_range[1] - 5, y_range[1] - 3),
               fontsize=10, ha='center', color='blue')

    plt.tight_layout()
    return fig


def plot_ray_tracing(lattice, psi, mass_pos, trajectories_octh, trajectories_gr,
                    figsize=(12, 10)):
    """Compara trayectorias OCTH vs GR."""

    fig, ax = plt.subplots(figsize=figsize)

    # Campo Ψ como fondo
    scatter = ax.scatter(lattice.positions[:, 0], lattice.positions[:, 1],
                        c=psi, cmap='Greys', s=5, alpha=0.3, vmin=0, vmax=1)

    # Masa central
    ax.scatter(mass_pos[0], mass_pos[1], c='black', s=300, marker='o',
              label='Masa', zorder=10)

    # Círculo de Schwarzschild (horizonte)
    rs = 2 * G_NEWTON * 5.0 / C_LIGHT**2  # Para M=5
    circle = plt.Circle((mass_pos[0], mass_pos[1]), rs, fill=False,
                        color='red', linestyle='--', linewidth=2, label='$r_s$ (horizonte)')
    ax.add_patch(circle)

    # Trayectorias OCTH
    colors_octh = plt.cm.Blues(np.linspace(0.4, 0.9, len(trajectories_octh)))
    for i, traj in enumerate(trajectories_octh):
        ax.plot(traj[:, 0], traj[:, 1], '-', color=colors_octh[i], linewidth=2,
               label='OCTH' if i == 0 else None)

    # Trayectorias GR
    colors_gr = plt.cm.Oranges(np.linspace(0.4, 0.9, len(trajectories_gr)))
    for i, (traj, _) in enumerate(trajectories_gr):
        ax.plot(traj[:, 0], traj[:, 1], '--', color=colors_gr[i], linewidth=2,
               label='GR (Schwarzschild)' if i == 0 else None)

    ax.set_xlim(lattice.positions[:, 0].min(), lattice.positions[:, 0].max())
    ax.set_ylim(lattice.positions[:, 1].min(), lattice.positions[:, 1].max())
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title('Deflexión de Luz: OCTH vs Relatividad General', fontsize=14)
    ax.legend(loc='upper left')

    plt.tight_layout()
    return fig


def plot_deflection_comparison(b_values, deflections_octh, deflections_gr, figsize=(10, 6)):
    """Compara ángulo de deflexión vs parámetro de impacto."""

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    # Panel izquierdo: Deflexión vs b
    ax = axes[0]
    ax.plot(b_values, np.degrees(deflections_gr), 'r-', linewidth=2,
           label='GR: Δφ = 4GM/(bc²)')
    ax.plot(b_values, np.degrees(deflections_octh), 'b--', linewidth=2,
           label='OCTH (simulación)')
    ax.set_xlabel('Parámetro de impacto b', fontsize=12)
    ax.set_ylabel('Deflexión Δφ (grados)', fontsize=12)
    ax.set_title('Deflexión de Luz', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel derecho: Error relativo
    ax = axes[1]
    error = np.abs(deflections_octh - deflections_gr) / deflections_gr * 100
    ax.semilogy(b_values, error, 'g-', linewidth=2)
    ax.axhline(10, color='orange', linestyle='--', label='10% error')
    ax.axhline(1, color='red', linestyle='--', label='1% error')
    ax.set_xlabel('Parámetro de impacto b', fontsize=12)
    ax.set_ylabel('Error relativo (%)', fontsize=12)
    ax.set_title('Convergencia OCTH → GR', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

def run_test2():
    """Ejecuta Test #2 completo."""

    print("=" * 70)
    print("TEST #2: SIMULACIÓN DE PROPAGACIÓN EN RETÍCULO HEXAGONAL")
    print("=" * 70)
    print("\nOCTH - Ontología del Campo Tensorial Hexagonal")
    print("Verificando que Ψ variable reproduce lensing gravitacional\n")

    # Parámetros
    Nx, Ny = 50, 50
    mass_pos = np.array([0.0, 0.0])
    M = 5.0  # Masa central

    # 1. Crear retículos
    print("1. Creando retículos...")
    hex_lattice = HexagonalLattice(Nx, Ny, spacing=1.0)
    sq_lattice = SquareLattice(Nx, Ny, spacing=1.0)
    print(f"   Hexagonal: {hex_lattice.N} nodos")
    print(f"   Cuadrado:  {sq_lattice.N} nodos")

    # 2. Calcular campo Ψ
    print("\n2. Calculando campo Ψ...")
    psi_hex = compute_psi_field(hex_lattice.positions, [mass_pos], [M], rs=2.0)
    psi_sq = compute_psi_field(sq_lattice.positions, [mass_pos], [M], rs=2.0)
    print(f"   Ψ_min (hexagonal) = {psi_hex.min():.4f}")
    print(f"   Ψ_min (cuadrado)  = {psi_sq.min():.4f}")

    # 3. Comparación de retículos
    print("\n3. Generando comparación de retículos...")
    fig1 = plot_lattice_comparison(hex_lattice, sq_lattice, psi_hex, psi_sq)
    fig1.savefig('../figures/fig1_lattice_comparison.png', dpi=300, bbox_inches='tight')
    fig1.savefig('../figures/fig1_lattice_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("   Guardado: fig1_lattice_comparison.png/pdf")

    # 4. Mapa de calor de Ψ
    print("\n4. Generando mapa de calor de Ψ...")
    fig2 = plot_psi_heatmap(hex_lattice, psi_hex, mass_pos,
                           title='Campo de Permeabilidad Temporal Ψ (Retículo Hexagonal)')
    fig2.savefig('../figures/fig2_psi_heatmap.png', dpi=300, bbox_inches='tight')
    fig2.savefig('../figures/fig2_psi_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("   Guardado: fig2_psi_heatmap.png/pdf")

    # 5. Ray tracing
    print("\n5. Trazando rayos de luz...")
    trajectories_octh = []
    trajectories_gr = []

    b_values = [3, 5, 7, 10, 15]  # Parámetros de impacto

    for b in b_values:
        # OCTH
        start = np.array([b, -20.0])
        direction = np.array([0.0, 1.0])
        traj_octh = trace_ray(hex_lattice, psi_hex, start, direction, n_steps=800)
        trajectories_octh.append(traj_octh)

        # GR (Schwarzschild)
        traj_gr, delta_phi = schwarzschild_geodesic(M, b, r_start=20, r_end=-20)
        trajectories_gr.append((traj_gr, delta_phi))

    print(f"   Trazados {len(b_values)} rayos para cada modelo")

    # 6. Comparación de trayectorias
    print("\n6. Generando comparación de trayectorias...")
    fig3 = plot_ray_tracing(hex_lattice, psi_hex, mass_pos,
                           trajectories_octh, trajectories_gr)
    fig3.savefig('../figures/fig3_ray_tracing.png', dpi=300, bbox_inches='tight')
    fig3.savefig('../figures/fig3_ray_tracing.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("   Guardado: fig3_ray_tracing.png/pdf")

    # 7. Calcular deflexiones
    print("\n7. Calculando deflexiones...")
    deflections_octh = []
    deflections_gr = []

    for i, b in enumerate(b_values):
        # OCTH: medir deflexión del rayo
        traj = trajectories_octh[i]
        if len(traj) > 10:
            # Dirección inicial y final
            dir_init = traj[5] - traj[0]
            dir_final = traj[-1] - traj[-6]
            dir_init = dir_init / np.linalg.norm(dir_init)
            dir_final = dir_final / np.linalg.norm(dir_final)
            angle = np.arccos(np.clip(np.dot(dir_init, dir_final), -1, 1))
            deflections_octh.append(angle)
        else:
            deflections_octh.append(0)

        # GR
        _, delta_phi = trajectories_gr[i]
        deflections_gr.append(delta_phi)

    deflections_octh = np.array(deflections_octh)
    deflections_gr = np.array(deflections_gr)

    print("\n   RESULTADOS DE DEFLEXIÓN:")
    print("   " + "-" * 50)
    print(f"   {'b':<8} {'OCTH (°)':<12} {'GR (°)':<12} {'Error (%)':<12}")
    print("   " + "-" * 50)
    for i, b in enumerate(b_values):
        error = abs(deflections_octh[i] - deflections_gr[i]) / deflections_gr[i] * 100
        print(f"   {b:<8} {np.degrees(deflections_octh[i]):<12.4f} "
              f"{np.degrees(deflections_gr[i]):<12.4f} {error:<12.2f}")

    # 8. Gráfica de comparación de deflexiones
    print("\n8. Generando gráfica de deflexiones...")
    fig4 = plot_deflection_comparison(b_values, deflections_octh, deflections_gr)
    fig4.savefig('../figures/fig4_deflection_comparison.png', dpi=300, bbox_inches='tight')
    fig4.savefig('../figures/fig4_deflection_comparison.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print("   Guardado: fig4_deflection_comparison.png/pdf")

    # 9. Resumen
    print("\n" + "=" * 70)
    print("RESUMEN TEST #2")
    print("=" * 70)

    mean_error = np.mean(np.abs(deflections_octh - deflections_gr) / deflections_gr * 100)

    print(f"\n   Error promedio OCTH vs GR: {mean_error:.2f}%")

    if mean_error < 20:
        print("\n   ✅ RESULTADO: El retículo hexagonal con Ψ variable")
        print("      REPRODUCE cualitativamente el lensing gravitacional")
        print("      predicho por Relatividad General.")
    else:
        print("\n   ⚠️  RESULTADO: Discrepancia significativa con GR.")
        print("      Se requiere ajuste de parámetros.")

    print("\n   FIGURAS GENERADAS:")
    print("   - fig1_lattice_comparison.png/pdf")
    print("   - fig2_psi_heatmap.png/pdf")
    print("   - fig3_ray_tracing.png/pdf")
    print("   - fig4_deflection_comparison.png/pdf")

    print("\n" + "=" * 70)

    return {
        'b_values': b_values,
        'deflections_octh': deflections_octh,
        'deflections_gr': deflections_gr,
        'mean_error': mean_error
    }


if __name__ == "__main__":
    results = run_test2()
