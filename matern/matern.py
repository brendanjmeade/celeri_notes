import numpy as np
from scipy import special


def spectral_matern(
    f: np.ndarray | float,
    sigma: float,
    nu: float,
    rho: float,
    n: int = 2,
    *,
    unnormalized: bool = False,
) -> np.ndarray | float:
    """Spectral density S(f) for a Matérn GP in n dimensions.

    If unnormalized=False (default, canonical form; requires ν>0):
        S(f) = σ² * [2^n π^{n/2} Γ(ν + n/2) (2ν)^ν] / [Γ(ν) ρ^{2ν}]
               * ((2ν/ρ²) + (2π f)²)^(-(ν + n/2))

    If unnormalized=True (extended form; allows any real ν):
        Use λ² = 2*max(ν, 0.5)/ρ² inside the frequency term and a simple σ² scale:
            S(f) ∝ σ² * (λ² + (2π f)²)^(-(ν + n/2))
        S(f) = σ² * [2^n π^{n/2} Γ(ν + n/2) (2ν)^ν] / [Γ(ν) ρ^{2ν}]
               * ((2ν/ρ²) + (2π f)²)^(-(ν + n/2))

    If unnormalized=True (extended form; allows any real ν):
        Use λ² = 2*max(ν, 0.5)/ρ² inside the frequency term and a simple σ² scale:
            S(f) ∝ σ² * (λ² + (2π f)²)^(-(ν + n/2))
    """
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if rho <= 0:
        raise ValueError("rho must be positive")
    if n < 1:
        raise ValueError("n must be at least 1")

    f_is_scalar = np.isscalar(f)
    f = np.atleast_1d(f).astype(float)

    omega2 = (2.0 * np.pi * f) ** 2
    exponent = -(nu + n / 2.0)

    if not unnormalized:
        if nu <= 0:
            raise ValueError("nu must be positive for normalized Matérn; set unnormalized=True to extend.")
        # log-prefactor:
        # log σ² + n log 2 + (n/2) log π + log Γ(ν + n/2) + ν log(2ν) - log Γ(ν) - 2ν log ρ
        log_pref = (
            2.0 * np.log(sigma)
            + n * np.log(2.0)
            + (n / 2.0) * np.log(np.pi)
            + special.gammaln(nu + n / 2.0)
            + nu * np.log(2.0 * nu)
            - special.gammaln(nu)
            - 2.0 * nu * np.log(rho)
        )
        lam2 = 2.0 * nu / (rho * rho)  # > 0
        log_base = np.log(lam2 + omega2)
        logS = log_pref + exponent * log_base
    else:
        # Extended form that works when nu is negative and matches at nu=0.5
        lam2 = 2.0 * max(nu, 0.5) / (rho * rho)
        log_pref = (
            2.0 * np.log(sigma)
            + n * np.log(2.0)
            + ((n - 1) / 2.0) * np.log(np.pi)
            + special.gammaln((n + 1) / 2.0)
            - np.log(rho)
        )
        log_base = np.log(lam2 + omega2)
        logS = log_pref + exponent * log_base

    S = np.exp(logS)
    if f_is_scalar:
        return float(S[0])
    return S
