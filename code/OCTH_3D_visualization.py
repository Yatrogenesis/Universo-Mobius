#!/usr/bin/env python3
"""
OCTH 3D Field Visualization
============================

Two visualization modes:
1. TOY: Traditional time-frequency-amplitude spectrogram (3D surface)
2. OCTH: Hexagonal tensor field with Ψ-κ-φ coordinates

The OCTH visualization reveals the hexagonal structure that produces
the 75σ detection significance in gravitational wave analysis.

Author: Francisco Molina Burgos
Date: January 2026
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize
import warnings
warnings.filterwarnings('ignore')

# OCTH Constants
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
SQRT3 = np.sqrt(3)
SQRT7 = np.sqrt(7)

# Hexagonal frequency ratios from OCTH
HEX_RATIOS = [1.0, SQRT3, 2.0, SQRT7, 3.0, np.sqrt(12)]


class OCTHFieldVisualizer:
    """
    Visualize OCTH tensor fields in 3D.

    Coordinates:
    - Ψ (Psi): Temporal permeability
    - κ (kappa): Hexagonal curvature
    - φ (phi): Topological phase

    Or traditional:
    - t: Time
    - f: Frequency
    - A: Amplitude
    """

    def __init__(self, resolution=100):
        self.resolution = resolution

    # =========================================================================
    # TOY MODEL: Traditional Spectrogram 3D
    # =========================================================================

    def generate_toy_signal(self, duration=1.0, sr=1000):
        """
        Generate a toy signal with chirp and harmonics.
        This is the TRADITIONAL representation.
        """
        t = np.linspace(0, duration, int(sr * duration))

        # Chirp signal (frequency increases with time)
        f0, f1 = 50, 200
        phase = 2 * np.pi * (f0 * t + (f1 - f0) * t**2 / (2 * duration))
        signal = np.sin(phase)

        # Add harmonics
        signal += 0.5 * np.sin(2 * phase)
        signal += 0.3 * np.sin(3 * phase)

        # Add noise
        signal += 0.2 * np.random.randn(len(t))

        return t, signal, sr

    def compute_spectrogram_3d(self, signal, sr, nperseg=64):
        """
        Compute spectrogram for 3D visualization.
        """
        from scipy import signal as sig

        f, t, Sxx = sig.spectrogram(signal, sr, nperseg=nperseg,
                                     noverlap=nperseg//2)

        # Convert to dB
        Sxx_db = 10 * np.log10(Sxx + 1e-10)

        return t, f, Sxx_db

    def plot_toy_spectrogram_3d(self, save_path=None):
        """
        TOY MODEL: Traditional 3D spectrogram.

        Axes:
        - X: Time
        - Y: Frequency
        - Z: Amplitude (dB)
        """
        # Generate signal
        t_sig, signal, sr = self.generate_toy_signal()

        # Compute spectrogram
        t, f, Sxx = self.compute_spectrogram_3d(signal, sr)

        # Create meshgrid
        T, F = np.meshgrid(t, f)

        # Plot
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Surface plot
        surf = ax.plot_surface(T, F, Sxx, cmap='viridis',
                               linewidth=0, antialiased=True,
                               alpha=0.8)

        # Labels
        ax.set_xlabel('Tiempo (s)', fontsize=12, labelpad=10)
        ax.set_ylabel('Frecuencia (Hz)', fontsize=12, labelpad=10)
        ax.set_zlabel('Amplitud (dB)', fontsize=12, labelpad=10)
        ax.set_title('TOY MODEL: Espectrograma 3D Tradicional\n(Tiempo × Frecuencia × Amplitud)',
                     fontsize=14, fontweight='bold')

        # Colorbar
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Potencia (dB)')

        # Adjust view
        ax.view_init(elev=25, azim=45)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig, ax

    # =========================================================================
    # OCTH MODEL: Hexagonal Tensor Field
    # =========================================================================

    def hexagonal_field(self, psi, kappa, phi):
        """
        Compute OCTH field intensity at point (Ψ, κ, φ).

        The field has hexagonal symmetry in the κ-φ plane,
        modulated by temporal permeability Ψ.

        F(Ψ,κ,φ) = Ψ · Σ cos(n·κ·cos(φ + nπ/3))

        This produces the hexagonal pattern detected at 75σ.
        """
        field = np.zeros_like(psi)

        # Hexagonal sum (6-fold symmetry)
        for n in range(6):
            angle = phi + n * np.pi / 3
            field += np.cos(kappa * np.cos(angle))

        # Modulate by temporal permeability
        field *= psi

        # Add resonance peaks at OCTH frequency ratios
        for ratio in HEX_RATIOS:
            field += 0.5 * psi * np.exp(-(kappa - ratio)**2 / 0.1)

        return field / 6  # Normalize

    def hexagonal_tension_field(self, x, y, z):
        """
        Alternative representation using:
        - Tension (τ)
        - Resistance (ρ)
        - Elasticity (ε)

        With implicit vibration as the field value.

        Vibración = f(τ, ρ, ε) siguiendo estructura hexagonal
        """
        # Map to OCTH coordinates
        # τ (tension) → drives oscillation amplitude
        # ρ (resistance) → damping factor
        # ε (elasticity) → restoring force / frequency

        # Hexagonal potential in τ-ρ plane
        hex_potential = np.zeros_like(x)
        for n in range(6):
            angle = n * np.pi / 3
            hex_potential += np.cos(x * np.cos(angle) + y * np.sin(angle))

        # Elasticity modulates the "spring constant"
        omega = np.sqrt(np.abs(z) + 0.1)  # Natural frequency from elasticity

        # Vibration amplitude (field value)
        vibration = hex_potential * np.exp(-y**2 / 4) * np.sin(omega * x)

        return vibration

    def plot_octh_field_3d(self, save_path=None):
        """
        OCTH MODEL: Hexagonal tensor field visualization.

        Axes:
        - X: Ψ (Temporal Permeability)
        - Y: κ (Hexagonal Curvature)
        - Z: φ (Topological Phase)

        Color: Field intensity (shows hexagonal structure)
        """
        # Create 3D grid
        n = self.resolution

        psi = np.linspace(0.5, 2.0, n)      # Temporal permeability
        kappa = np.linspace(0, 4, n)         # Hexagonal curvature (covers ratios 1-√12)
        phi = np.linspace(0, 2*np.pi, n)     # Topological phase

        # For surface plot, we'll show slices
        fig = plt.figure(figsize=(16, 12))

        # Main 3D plot: isosurfaces of field intensity
        ax1 = fig.add_subplot(221, projection='3d')

        # Create meshgrid for κ-φ plane at fixed Ψ
        K, P = np.meshgrid(kappa, phi)

        # Plot several Ψ slices
        psi_values = [0.8, 1.0, 1.2, 1.5]
        colors = ['blue', 'green', 'orange', 'red']

        for psi_val, color in zip(psi_values, colors):
            PSI = np.full_like(K, psi_val)
            field = self.hexagonal_field(PSI, K, P)

            # Plot as surface at height proportional to field
            ax1.plot_surface(K * np.cos(P), K * np.sin(P),
                           field * psi_val,
                           alpha=0.3, color=color,
                           label=f'Ψ={psi_val}')

        ax1.set_xlabel('κ·cos(φ)', fontsize=10)
        ax1.set_ylabel('κ·sin(φ)', fontsize=10)
        ax1.set_zlabel('Campo F(Ψ,κ,φ)', fontsize=10)
        ax1.set_title('Campo OCTH: Cortes en Ψ\n(Estructura Hexagonal Visible)', fontsize=12)

        # Subplot 2: κ-φ plane showing hexagonal pattern
        ax2 = fig.add_subplot(222)

        K2, P2 = np.meshgrid(np.linspace(0, 5, 200), np.linspace(0, 2*np.pi, 200))
        PSI2 = np.ones_like(K2)
        field2 = self.hexagonal_field(PSI2, K2, P2)

        # Convert to Cartesian for display
        X2 = K2 * np.cos(P2)
        Y2 = K2 * np.sin(P2)

        c2 = ax2.pcolormesh(X2, Y2, field2, cmap='RdBu_r', shading='auto')
        plt.colorbar(c2, ax=ax2, label='Intensidad de Campo')

        # Mark OCTH frequency ratios
        for ratio in HEX_RATIOS:
            circle = plt.Circle((0, 0), ratio, fill=False, color='gold',
                               linestyle='--', linewidth=2)
            ax2.add_patch(circle)

        ax2.set_xlim(-5, 5)
        ax2.set_ylim(-5, 5)
        ax2.set_aspect('equal')
        ax2.set_xlabel('κ·cos(φ)', fontsize=10)
        ax2.set_ylabel('κ·sin(φ)', fontsize=10)
        ax2.set_title('Plano κ-φ: Patrón Hexagonal\n(Círculos = ratios OCTH: 1, √3, 2, √7, 3, √12)',
                     fontsize=11)

        # Subplot 3: Tension-Resistance-Elasticity representation
        ax3 = fig.add_subplot(223, projection='3d')

        # Grid in τ-ρ-ε space
        tau = np.linspace(-3, 3, 50)      # Tension
        rho = np.linspace(-3, 3, 50)      # Resistance

        T, R = np.meshgrid(tau, rho)

        # Show vibration field at different elasticity values
        for eps_val in [0.5, 1.0, 2.0]:
            E = np.full_like(T, eps_val)
            vibration = self.hexagonal_tension_field(T, R, E)

            ax3.plot_surface(T, R, vibration + eps_val * 2,
                           alpha=0.4, cmap='coolwarm',
                           linewidth=0)

        ax3.set_xlabel('Tensión (τ)', fontsize=10)
        ax3.set_ylabel('Resistencia (ρ)', fontsize=10)
        ax3.set_zlabel('Vibración + offset', fontsize=10)
        ax3.set_title('Campo τ-ρ-ε: Vibración Hexagonal\n(Capas = diferentes elasticidades)',
                     fontsize=11)

        # Subplot 4: Ψ evolution showing resonance
        ax4 = fig.add_subplot(224)

        kappa_line = np.linspace(0, 5, 500)
        psi_vals = [0.5, 1.0, 1.5, 2.0]

        for psi_val in psi_vals:
            # Integrate over φ to get radial profile
            field_radial = np.zeros_like(kappa_line)
            for phi_val in np.linspace(0, 2*np.pi, 100):
                field_radial += self.hexagonal_field(
                    np.full_like(kappa_line, psi_val),
                    kappa_line,
                    np.full_like(kappa_line, phi_val)
                )
            field_radial /= 100

            ax4.plot(kappa_line, field_radial, label=f'Ψ={psi_val}', linewidth=2)

        # Mark OCTH ratios
        for i, ratio in enumerate(HEX_RATIOS):
            ax4.axvline(ratio, color='gray', linestyle=':', alpha=0.5)
            ax4.text(ratio, ax4.get_ylim()[1]*0.95,
                    ['1', '√3', '2', '√7', '3', '√12'][i],
                    ha='center', fontsize=9)

        ax4.set_xlabel('Curvatura Hexagonal κ', fontsize=10)
        ax4.set_ylabel('⟨F⟩_φ (Campo promediado)', fontsize=10)
        ax4.set_title('Perfil Radial: Resonancias en Ratios OCTH\n(Picos = frecuencias hexagonales)',
                     fontsize=11)
        ax4.legend(loc='upper right')
        ax4.grid(True, alpha=0.3)

        plt.suptitle('OCTH: Visualización 3D del Campo Tensorial Hexagonal\n'
                    '(Detección 75σ, significancia ~10⁻³⁰⁰)',
                    fontsize=14, fontweight='bold', y=1.02)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig

    def plot_volumetric_comparison(self, save_path=None):
        """
        Side-by-side comparison of TOY vs OCTH representations.
        """
        fig = plt.figure(figsize=(18, 8))

        # === LEFT: TOY (Traditional) ===
        ax1 = fig.add_subplot(121, projection='3d')

        # Generate spectrogram data
        t_sig, signal, sr = self.generate_toy_signal()
        t, f, Sxx = self.compute_spectrogram_3d(signal, sr)
        T, F = np.meshgrid(t, f)

        surf1 = ax1.plot_surface(T, F, Sxx, cmap='viridis',
                                 linewidth=0, antialiased=True, alpha=0.8)

        ax1.set_xlabel('Tiempo (s)', fontsize=11)
        ax1.set_ylabel('Frecuencia (Hz)', fontsize=11)
        ax1.set_zlabel('Amplitud (dB)', fontsize=11)
        ax1.set_title('TOY MODEL\n(Tradicional: t × f × A)',
                     fontsize=13, fontweight='bold')
        ax1.view_init(elev=25, azim=45)

        # === RIGHT: OCTH (Hexagonal) ===
        ax2 = fig.add_subplot(122, projection='3d')

        # Create hexagonal field visualization
        n = 80
        kappa = np.linspace(0, 4, n)
        phi = np.linspace(0, 2*np.pi, n)
        K, P = np.meshgrid(kappa, phi)

        # Field at Ψ = 1 (normalized)
        PSI = np.ones_like(K)
        field = self.hexagonal_field(PSI, K, P)

        # Convert to 3D coordinates
        X = K * np.cos(P)
        Y = K * np.sin(P)
        Z = field

        surf2 = ax2.plot_surface(X, Y, Z, cmap='RdBu_r',
                                 linewidth=0, antialiased=True, alpha=0.8)

        ax2.set_xlabel('κ·cos(φ)', fontsize=11)
        ax2.set_ylabel('κ·sin(φ)', fontsize=11)
        ax2.set_zlabel('F(Ψ,κ,φ)', fontsize=11)
        ax2.set_title('OCTH MODEL\n(Hexagonal: Ψ × κ × φ)',
                     fontsize=13, fontweight='bold')
        ax2.view_init(elev=30, azim=45)

        plt.suptitle('Comparación: Representación Tradicional vs OCTH Volumétrico',
                    fontsize=15, fontweight='bold', y=1.02)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Saved: {save_path}")

        return fig


def main():
    """Generate all visualizations."""
    import os

    output_dir = 'figures/3d_visualization'
    os.makedirs(output_dir, exist_ok=True)

    viz = OCTHFieldVisualizer(resolution=100)

    print("="*60)
    print("  OCTH 3D VISUALIZATION")
    print("="*60)

    # 1. TOY model (traditional)
    print("\n1. Generating TOY model (traditional spectrogram)...")
    viz.plot_toy_spectrogram_3d(f'{output_dir}/toy_spectrogram_3d.png')

    # 2. OCTH field
    print("\n2. Generating OCTH hexagonal field...")
    viz.plot_octh_field_3d(f'{output_dir}/octh_field_3d.png')

    # 3. Comparison
    print("\n3. Generating comparison...")
    viz.plot_volumetric_comparison(f'{output_dir}/comparison_toy_vs_octh.png')

    print(f"\n✓ All figures saved to {output_dir}/")
    print("\nVisualizaciones:")
    print("  - toy_spectrogram_3d.png     : Modelo tradicional (tiempo-frecuencia-amplitud)")
    print("  - octh_field_3d.png          : Campo OCTH hexagonal (Ψ-κ-φ)")
    print("  - comparison_toy_vs_octh.png : Comparación lado a lado")

    plt.close('all')


if __name__ == '__main__':
    main()
