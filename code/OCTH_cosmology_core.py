#!/usr/bin/env python3
"""
OCTH COSMOLOGY CORE - Modulo base para calculos cosmologicos
=============================================================

Este modulo provee las funciones fundamentales para calcular
cantidades cosmologicas en OCTH (Ontologia del Campo Tensorial Hexagonal).

Funciones principales:
  - Psi(z): Campo de permeabilidad temporal
  - H(z): Parametro de Hubble
  - r_s(z): Horizonte acustico (sound horizon)
  - D_C(z): Distancia comoving
  - D_A(z): Distancia angular
  - D_V(z): Distancia volumetrica (BAO)

Author: F. Molina-Burgos
Date: 13 January 2026
Version: 2.0 (corregido)
"""

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from dataclasses import dataclass
from typing import Optional, Tuple, Callable

# =============================================================================
# CONSTANTES FISICAS
# =============================================================================

C_KMS = 2.998e5           # Velocidad de luz en km/s
C_MPS = 2.998e8           # Velocidad de luz en m/s
G_SI = 6.674e-11          # Constante gravitacional m^3/kg/s^2
HBAR = 1.055e-34          # Constante de Planck reducida J*s
K_B = 1.381e-23           # Constante de Boltzmann J/K
T_CMB = 2.7255            # Temperatura CMB en K

# =============================================================================
# PARAMETROS COSMOLOGICOS (Planck 2018 baseline)
# =============================================================================

@dataclass
class CosmologyParams:
    """Parametros cosmologicos."""
    H_0: float = 67.4           # km/s/Mpc
    Omega_b: float = 0.0493     # Densidad barionica
    Omega_c: float = 0.264      # Densidad CDM (LCDM)
    Omega_m: float = 0.315      # Densidad materia total
    Omega_r: float = 9.0e-5     # Densidad radiacion
    Omega_Lambda: float = 0.685 # Densidad energia oscura
    T_CMB: float = 2.7255       # Temperatura CMB K
    N_eff: float = 3.046        # Numero efectivo neutrinos

    @property
    def h(self) -> float:
        return self.H_0 / 100.0

    @property
    def Omega_b_h2(self) -> float:
        return self.Omega_b * self.h**2

    @property
    def Omega_m_h2(self) -> float:
        return self.Omega_m * self.h**2

    @property
    def Omega_gamma(self) -> float:
        """Densidad de fotones."""
        return 2.47e-5 / self.h**2


@dataclass
class OCTHParams:
    """Parametros del modelo OCTH."""
    a_0: float = 1.2e-10        # Aceleracion caracteristica m/s^2
    sigma: float = 0.6          # Anchura del efecto topologico
    z_topo: float = 1089        # Redshift del pico topologico

    @property
    def epsilon(self) -> float:
        """epsilon = a_0 / a_H donde a_H = c * H_0."""
        H0_SI = 67.4 / 3.086e19  # H_0 en s^-1
        a_H = C_MPS * H0_SI
        return self.a_0 / a_H


# Instancias por defecto
PLANCK2018 = CosmologyParams()
OCTH_DEFAULT = OCTHParams()


# =============================================================================
# CAMPO DE PERMEABILIDAD TEMPORAL Psi(z)
# =============================================================================

def Psi(z: float, params: OCTHParams = OCTH_DEFAULT) -> float:
    """
    Campo de permeabilidad temporal Psi(z).

    En OCTH, Psi modifica la metrica:
        ds^2 = -c^2 * Psi^2 * dt^2 + g_ij dx^i dx^j

    Psi = 1 - epsilon * f_topo(z)

    donde f_topo es un perfil Gaussiano centrado en z_topo.

    Args:
        z: Redshift
        params: Parametros OCTH

    Returns:
        Valor de Psi (0 < Psi <= 1)
    """
    ln_a = -np.log(1 + z)
    ln_a_topo = -np.log(1 + params.z_topo)

    # Perfil Gaussiano en ln(a)
    f_topo = np.exp(-0.5 * ((ln_a - ln_a_topo) / params.sigma)**2)

    return 1.0 - params.epsilon * f_topo


def Psi_array(z_array: np.ndarray, params: OCTHParams = OCTH_DEFAULT) -> np.ndarray:
    """Version vectorizada de Psi."""
    return np.array([Psi(z, params) for z in z_array])


# =============================================================================
# PARAMETRO DE HUBBLE H(z)
# =============================================================================

def H_LCDM(z: float, cosmo: CosmologyParams = PLANCK2018) -> float:
    """
    Parametro de Hubble en LCDM estandar.

    H(z) = H_0 * E(z)
    E(z)^2 = Omega_r*(1+z)^4 + Omega_m*(1+z)^3 + Omega_Lambda

    Args:
        z: Redshift
        cosmo: Parametros cosmologicos

    Returns:
        H(z) en km/s/Mpc
    """
    E2 = (cosmo.Omega_r * (1+z)**4 +
          cosmo.Omega_m * (1+z)**3 +
          cosmo.Omega_Lambda)
    return cosmo.H_0 * np.sqrt(E2)


def H_OCTH(z: float,
           cosmo: CosmologyParams = PLANCK2018,
           octh: OCTHParams = OCTH_DEFAULT) -> float:
    """
    Parametro de Hubble en OCTH.

    En OCTH, el potencial gravitacional es mas profundo cuando Psi < 1,
    lo que equivale a una densidad de materia efectiva mayor:

        Omega_m_eff = Omega_m / Psi^2

    Args:
        z: Redshift
        cosmo: Parametros cosmologicos
        octh: Parametros OCTH

    Returns:
        H(z) en km/s/Mpc
    """
    psi = Psi(z, octh)
    Omega_m_eff = cosmo.Omega_m / psi**2

    E2 = (cosmo.Omega_r * (1+z)**4 +
          Omega_m_eff * (1+z)**3 +
          cosmo.Omega_Lambda)
    return cosmo.H_0 * np.sqrt(E2)


# =============================================================================
# VELOCIDAD DEL SONIDO Y BARYON LOADING
# =============================================================================

def R_baryon(z: float, cosmo: CosmologyParams = PLANCK2018) -> float:
    """
    Baryon-to-photon momentum density ratio.

    R = 3 * rho_b / (4 * rho_gamma)
      = (3/4) * (Omega_b / Omega_gamma) / (1+z)

    NOTA: En OCTH, R NO se modifica porque depende de densidades,
    no del potencial gravitacional.

    Args:
        z: Redshift
        cosmo: Parametros cosmologicos

    Returns:
        R(z)
    """
    return 0.75 * (cosmo.Omega_b / cosmo.Omega_gamma) / (1 + z)


def c_s(z: float, cosmo: CosmologyParams = PLANCK2018) -> float:
    """
    Velocidad del sonido en el fluido barion-foton.

    c_s = c / sqrt(3 * (1 + R))

    NOTA: En OCTH, c_s NO cambia directamente porque R depende
    de densidades, no del potencial. Solo H(z) cambia.

    Args:
        z: Redshift
        cosmo: Parametros cosmologicos

    Returns:
        c_s(z) en km/s
    """
    R = R_baryon(z, cosmo)
    return C_KMS / np.sqrt(3.0 * (1.0 + R))


# =============================================================================
# DRAG EPOCH (z_drag)
# =============================================================================

def z_drag_EH98(cosmo: CosmologyParams = PLANCK2018) -> float:
    """
    Redshift del drag epoch usando formula de Eisenstein & Hu (1998).

    Args:
        cosmo: Parametros cosmologicos

    Returns:
        z_drag
    """
    Omega_b_h2 = cosmo.Omega_b_h2
    Omega_m_h2 = cosmo.Omega_m_h2

    b1 = 0.313 * Omega_m_h2**(-0.419) * (1 + 0.607 * Omega_m_h2**0.674)
    b2 = 0.238 * Omega_m_h2**0.223

    z_d = (1291 * Omega_m_h2**0.251 /
           (1 + 0.659 * Omega_m_h2**0.828) *
           (1 + b1 * Omega_b_h2**b2))

    return z_d


def z_star_EH98(cosmo: CosmologyParams = PLANCK2018) -> float:
    """
    Redshift de la superficie de ultima dispersion (decoupling).
    Formula de Eisenstein & Hu (1998).

    Args:
        cosmo: Parametros cosmologicos

    Returns:
        z_star (z de decoupling)
    """
    Omega_m_h2 = cosmo.Omega_m_h2
    Omega_b_h2 = cosmo.Omega_b_h2

    g1 = 0.0783 * Omega_b_h2**(-0.238) / (1 + 39.5 * Omega_b_h2**0.763)
    g2 = 0.560 / (1 + 21.1 * Omega_b_h2**1.81)

    z_star = 1048 * (1 + 0.00124 * Omega_b_h2**(-0.738)) * (1 + g1 * Omega_m_h2**g2)

    return z_star


# =============================================================================
# HORIZONTE ACUSTICO r_s (SOUND HORIZON)
# =============================================================================

def r_s(z_final: float,
        H_func: Callable[[float], float],
        cosmo: CosmologyParams = PLANCK2018,
        z_max: float = 1e6) -> float:
    """
    Calcula el horizonte acustico r_s(z).

    r_s = integral_z^inf c_s(z') / H(z') dz'

    Esta es la distancia comoving que el sonido puede viajar
    desde el Big Bang hasta el redshift z_final.

    Args:
        z_final: Redshift final (tipicamente z_drag)
        H_func: Funcion H(z) a usar (H_LCDM o H_OCTH)
        cosmo: Parametros cosmologicos
        z_max: Limite superior de integracion

    Returns:
        r_s en Mpc (comoving)
    """
    def integrand(z):
        return c_s(z, cosmo) / H_func(z)

    result, error = quad(integrand, z_final, z_max, limit=500)
    return result


def r_s_LCDM(cosmo: CosmologyParams = PLANCK2018) -> float:
    """Horizonte acustico en LCDM."""
    z_d = z_drag_EH98(cosmo)
    return r_s(z_d, lambda z: H_LCDM(z, cosmo), cosmo)


def r_s_OCTH(cosmo: CosmologyParams = PLANCK2018,
             octh: OCTHParams = OCTH_DEFAULT) -> float:
    """Horizonte acustico en OCTH."""
    z_d = z_drag_EH98(cosmo)
    return r_s(z_d, lambda z: H_OCTH(z, cosmo, octh), cosmo)


# =============================================================================
# DISTANCIAS COSMOLOGICAS
# =============================================================================

def D_C(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia comoving D_C(z).

    D_C = integral_0^z c/H(z') dz'

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_C en Mpc (comoving)
    """
    def integrand(zp):
        return C_KMS / H_func(zp)

    result, error = quad(integrand, 0, z, limit=200)
    return result


def D_M(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia comoving transversa D_M(z).
    Para universo plano, D_M = D_C.

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_M en Mpc
    """
    return D_C(z, H_func)


def D_A(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia angular D_A(z).

    D_A = D_M / (1+z)

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_A en Mpc (proper)
    """
    return D_M(z, H_func) / (1 + z)


def D_L(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia de luminosidad D_L(z).

    D_L = D_M * (1+z)

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_L en Mpc
    """
    return D_M(z, H_func) * (1 + z)


def D_H(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia de Hubble D_H(z).

    D_H = c / H(z)

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_H en Mpc
    """
    return C_KMS / H_func(z)


def D_V(z: float,
        H_func: Callable[[float], float]) -> float:
    """
    Distancia volumetrica D_V(z) para BAO.

    D_V = [z * D_H * D_M^2]^(1/3)

    Args:
        z: Redshift
        H_func: Funcion H(z)

    Returns:
        D_V en Mpc
    """
    d_m = D_M(z, H_func)
    d_h = D_H(z, H_func)
    return (z * d_h * d_m**2)**(1.0/3.0)


# =============================================================================
# OBSERVABLES BAO
# =============================================================================

def DV_rd(z: float,
          H_func: Callable[[float], float],
          r_d: float) -> float:
    """
    D_V(z) / r_d - observable BAO isotropico.

    Args:
        z: Redshift
        H_func: Funcion H(z)
        r_d: Sound horizon at drag epoch (Mpc)

    Returns:
        D_V/r_d (adimensional)
    """
    return D_V(z, H_func) / r_d


def DM_rd(z: float,
          H_func: Callable[[float], float],
          r_d: float) -> float:
    """
    D_M(z) / r_d - observable BAO transverso.

    Args:
        z: Redshift
        H_func: Funcion H(z)
        r_d: Sound horizon (Mpc)

    Returns:
        D_M/r_d (adimensional)
    """
    return D_M(z, H_func) / r_d


def DH_rd(z: float,
          H_func: Callable[[float], float],
          r_d: float) -> float:
    """
    D_H(z) / r_d - observable BAO radial.

    Args:
        z: Redshift
        H_func: Funcion H(z)
        r_d: Sound horizon (Mpc)

    Returns:
        D_H/r_d (adimensional)
    """
    return D_H(z, H_func) / r_d


# =============================================================================
# OBSERVABLE CMB: theta_s
# =============================================================================

def theta_s(cosmo: CosmologyParams = PLANCK2018,
            H_func: Optional[Callable[[float], float]] = None,
            r_d: Optional[float] = None) -> float:
    """
    Angulo subtendido por el horizonte acustico en el CMB.

    theta_s = r_s(z_drag) / D_A(z_star)

    Este es el OBSERVABLE DIRECTO del CMB, no r_s.

    Args:
        cosmo: Parametros cosmologicos
        H_func: Funcion H(z). Si None, usa H_LCDM.
        r_d: Sound horizon. Si None, lo calcula.

    Returns:
        theta_s en radianes
    """
    if H_func is None:
        H_func = lambda z: H_LCDM(z, cosmo)

    if r_d is None:
        z_d = z_drag_EH98(cosmo)
        r_d = r_s(z_d, H_func, cosmo)

    z_star = z_star_EH98(cosmo)
    d_a = D_C(z_star, H_func)  # Usamos D_C (comoving) para comparar con r_s (comoving)

    return r_d / d_a


# =============================================================================
# FUNCIONES DE CONVENIENCIA
# =============================================================================

def compute_all_LCDM(cosmo: CosmologyParams = PLANCK2018) -> dict:
    """Calcula todas las cantidades relevantes para LCDM."""
    H_func = lambda z: H_LCDM(z, cosmo)

    z_d = z_drag_EH98(cosmo)
    z_s = z_star_EH98(cosmo)
    r_d = r_s(z_d, H_func, cosmo)

    return {
        'model': 'LCDM',
        'z_drag': z_d,
        'z_star': z_s,
        'r_s': r_d,
        'theta_s': theta_s(cosmo, H_func, r_d),
        'D_C_star': D_C(z_s, H_func),
        'H_0': cosmo.H_0,
        'Omega_m': cosmo.Omega_m,
    }


def compute_all_OCTH(cosmo: CosmologyParams = PLANCK2018,
                     octh: OCTHParams = OCTH_DEFAULT) -> dict:
    """Calcula todas las cantidades relevantes para OCTH."""
    H_func = lambda z: H_OCTH(z, cosmo, octh)

    z_d = z_drag_EH98(cosmo)
    z_s = z_star_EH98(cosmo)
    r_d = r_s(z_d, H_func, cosmo)

    return {
        'model': 'OCTH',
        'z_drag': z_d,
        'z_star': z_s,
        'r_s': r_d,
        'theta_s': theta_s(cosmo, H_func, r_d),
        'D_C_star': D_C(z_s, H_func),
        'H_0': cosmo.H_0,
        'Omega_m': cosmo.Omega_m,
        'epsilon': octh.epsilon,
        'sigma': octh.sigma,
        'Psi_at_drag': Psi(z_d, octh),
    }


# =============================================================================
# MAIN - Test basico
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("OCTH COSMOLOGY CORE - Test basico")
    print("=" * 70)

    lcdm = compute_all_LCDM()
    octh = compute_all_OCTH()

    print(f"\n{'Cantidad':<20} {'LCDM':>15} {'OCTH':>15} {'Diff':>10}")
    print("-" * 65)
    print(f"{'z_drag':<20} {lcdm['z_drag']:>15.2f} {octh['z_drag']:>15.2f} {'-':>10}")
    print(f"{'z_star':<20} {lcdm['z_star']:>15.2f} {octh['z_star']:>15.2f} {'-':>10}")
    print(f"{'r_s (Mpc)':<20} {lcdm['r_s']:>15.2f} {octh['r_s']:>15.2f} {100*(octh['r_s']/lcdm['r_s']-1):>+9.1f}%")
    print(f"{'D_C(z*) (Mpc)':<20} {lcdm['D_C_star']:>15.2f} {octh['D_C_star']:>15.2f} {100*(octh['D_C_star']/lcdm['D_C_star']-1):>+9.1f}%")
    print(f"{'theta_s (rad)':<20} {lcdm['theta_s']:>15.6f} {octh['theta_s']:>15.6f} {100*(octh['theta_s']/lcdm['theta_s']-1):>+9.1f}%")
    print(f"{'theta_s (arcmin)':<20} {np.degrees(lcdm['theta_s'])*60:>15.2f} {np.degrees(octh['theta_s'])*60:>15.2f} {'-':>10}")

    print(f"\nPsi(z_drag) = {octh['Psi_at_drag']:.4f}")
    print(f"epsilon = {octh['epsilon']:.4f}")

    # Comparacion con Planck
    theta_planck = 0.0104110  # rad
    print(f"\nComparacion con Planck (theta_s = {theta_planck:.6f} rad):")
    print(f"  LCDM:  {100*(lcdm['theta_s']/theta_planck - 1):+.2f}% off")
    print(f"  OCTH:  {100*(octh['theta_s']/theta_planck - 1):+.2f}% off")

    print("\n" + "=" * 70)
    print("Test completado")
    print("=" * 70)
