"""
OCTH Statistical Analysis
=========================

Statistical tools for validating hexagonal patterns against null hypotheses.
Includes binomial tests, Monte Carlo simulations, and combined significance
calculations.
"""

import numpy as np
from scipy import stats
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass


@dataclass
class SignificanceResult:
    """Result of a significance calculation."""
    sigma: float              # Significance in standard deviations
    p_value: float            # Two-tailed p-value
    z_score: float            # Z-score from binomial test
    observed: int             # Observed count
    expected: float           # Expected count under null
    total: int                # Total trials
    null_probability: float   # Probability under null hypothesis
    method: str               # Method used for calculation


def calculate_binomial_significance(
    successes: int,
    trials: int,
    null_probability: float = 0.5,
    alternative: str = 'greater'
) -> SignificanceResult:
    """
    Calculate significance of observed success rate vs null hypothesis.

    For OCTH: If hexagonal patterns are random, we expect ~50% of
    frequency ratios to match hexagonal ratios by chance.

    Args:
        successes: Number of observed successes (hexagonal matches)
        trials: Total number of trials
        null_probability: Expected probability under null (default 0.5)
        alternative: 'greater', 'less', or 'two-sided'

    Returns:
        SignificanceResult with statistical measures
    """
    if trials == 0:
        return SignificanceResult(
            sigma=0.0, p_value=1.0, z_score=0.0,
            observed=0, expected=0.0, total=0,
            null_probability=null_probability, method='binomial'
        )

    # Expected value under null
    expected = trials * null_probability

    # Binomial test
    result = stats.binomtest(successes, trials, null_probability, alternative=alternative)
    p_value = result.pvalue

    # Z-score approximation (for large n)
    std = np.sqrt(trials * null_probability * (1 - null_probability))
    if std > 0:
        z_score = (successes - expected) / std
    else:
        z_score = 0.0

    # Convert to sigma (one-tailed to two-tailed adjustment)
    if p_value > 0:
        sigma = stats.norm.ppf(1 - p_value / 2)
    else:
        sigma = np.inf

    return SignificanceResult(
        sigma=abs(sigma),
        p_value=p_value,
        z_score=z_score,
        observed=successes,
        expected=expected,
        total=trials,
        null_probability=null_probability,
        method='binomial'
    )


def calculate_combined_significance(
    results: List[SignificanceResult],
    method: str = 'fisher'
) -> SignificanceResult:
    """
    Combine multiple significance results into one.

    Methods:
    - 'fisher': Fisher's method (combines p-values)
    - 'stouffer': Stouffer's method (combines z-scores)
    - 'weighted': Weighted Stouffer (by sample size)

    Args:
        results: List of SignificanceResult objects
        method: Combination method

    Returns:
        Combined SignificanceResult
    """
    if not results:
        return SignificanceResult(
            sigma=0.0, p_value=1.0, z_score=0.0,
            observed=0, expected=0.0, total=0,
            null_probability=0.5, method=f'combined_{method}'
        )

    # Aggregate statistics
    total_observed = sum(r.observed for r in results)
    total_expected = sum(r.expected for r in results)
    total_trials = sum(r.total for r in results)

    if method == 'fisher':
        # Fisher's method: chi-squared from -2 * sum(log(p))
        p_values = [r.p_value for r in results if r.p_value > 0]
        if not p_values:
            combined_p = 0.0
        else:
            chi2_stat = -2 * sum(np.log(p) for p in p_values)
            df = 2 * len(p_values)
            combined_p = 1 - stats.chi2.cdf(chi2_stat, df)

        combined_z = stats.norm.ppf(1 - combined_p / 2) if combined_p > 0 else np.inf

    elif method == 'stouffer':
        # Stouffer's method: sum z-scores / sqrt(n)
        z_scores = [r.z_score for r in results]
        combined_z = sum(z_scores) / np.sqrt(len(z_scores))
        combined_p = 2 * (1 - stats.norm.cdf(abs(combined_z)))

    elif method == 'weighted':
        # Weighted Stouffer: weight by sqrt(sample size)
        weights = [np.sqrt(r.total) for r in results]
        z_scores = [r.z_score for r in results]

        combined_z = sum(w * z for w, z in zip(weights, z_scores)) / np.sqrt(sum(w**2 for w in weights))
        combined_p = 2 * (1 - stats.norm.cdf(abs(combined_z)))

    else:
        raise ValueError(f"Unknown combination method: {method}")

    # Average null probability
    avg_null = np.mean([r.null_probability for r in results])

    return SignificanceResult(
        sigma=abs(combined_z),
        p_value=combined_p,
        z_score=combined_z,
        observed=total_observed,
        expected=total_expected,
        total=total_trials,
        null_probability=avg_null,
        method=f'combined_{method}'
    )


def monte_carlo_null(
    observed_statistic: float,
    null_generator: callable,
    statistic_function: callable,
    n_iterations: int = 10000,
    seed: int = None
) -> Dict:
    """
    Monte Carlo estimation of null distribution.

    Args:
        observed_statistic: The observed test statistic
        null_generator: Function that generates null hypothesis data
        statistic_function: Function to compute statistic from data
        n_iterations: Number of Monte Carlo iterations
        seed: Random seed for reproducibility

    Returns:
        Dictionary with p-value, null distribution, and percentile
    """
    if seed is not None:
        np.random.seed(seed)

    null_statistics = []

    for _ in range(n_iterations):
        null_data = null_generator()
        null_stat = statistic_function(null_data)
        null_statistics.append(null_stat)

    null_statistics = np.array(null_statistics)

    # One-tailed p-value (probability of observing value >= observed)
    p_value = np.mean(null_statistics >= observed_statistic)

    # Percentile of observed in null distribution
    percentile = np.mean(null_statistics <= observed_statistic) * 100

    # Null distribution statistics
    null_mean = np.mean(null_statistics)
    null_std = np.std(null_statistics)

    # Z-score of observed
    z_score = (observed_statistic - null_mean) / max(null_std, 1e-10)

    return {
        "observed": observed_statistic,
        "p_value": max(p_value, 1 / n_iterations),  # Floor at 1/N
        "percentile": percentile,
        "z_score": z_score,
        "null_mean": null_mean,
        "null_std": null_std,
        "null_distribution": null_statistics,
        "n_iterations": n_iterations
    }


def hexagonal_null_generator(
    n_peaks: int,
    freq_range: Tuple[float, float] = (10, 1000),
    seed: int = None
) -> callable:
    """
    Create a null generator for hexagonal ratio testing.

    Under null hypothesis, peak frequencies are uniformly distributed,
    and ratios should not preferentially match hexagonal values.

    Args:
        n_peaks: Number of frequency peaks
        freq_range: Range of frequencies
        seed: Random seed

    Returns:
        Generator function
    """
    rng = np.random.default_rng(seed)

    def generator():
        # Random frequencies (uniform in log space is more realistic)
        log_freqs = rng.uniform(
            np.log10(freq_range[0]),
            np.log10(freq_range[1]),
            n_peaks
        )
        return 10 ** log_freqs

    return generator


def hexagonal_statistic(
    frequencies: np.ndarray,
    hexagonal_ratios: List[float] = None,
    tolerance: float = 0.03
) -> float:
    """
    Compute hexagonal matching statistic from frequency list.

    Returns fraction of pairwise ratios that match hexagonal ratios.

    Args:
        frequencies: Array of frequencies
        hexagonal_ratios: Target ratios (default from constants)
        tolerance: Matching tolerance

    Returns:
        Fraction of hexagonal matches (0-1)
    """
    if hexagonal_ratios is None:
        from .constants import HEXAGONAL_RATIOS
        hexagonal_ratios = HEXAGONAL_RATIOS

    if len(frequencies) < 2:
        return 0.0

    matches = 0
    total = 0

    for i, f1 in enumerate(frequencies):
        for j, f2 in enumerate(frequencies):
            if j <= i:
                continue

            # Compute ratio (always > 1)
            ratio = max(f1, f2) / min(f1, f2)

            # Check if matches any hexagonal ratio
            for hex_ratio in hexagonal_ratios:
                if abs(ratio - hex_ratio) / hex_ratio <= tolerance:
                    matches += 1
                    break

            total += 1

    return matches / total if total > 0 else 0.0


def validate_gwtc3_result(
    hexagonal_events: int = 80,
    total_events: int = 80,
    null_probability: float = 0.5
) -> SignificanceResult:
    """
    Validate the GWTC-3 hexagonal detection result.

    This is the flagship result: 80/80 events show hexagonal patterns.

    Args:
        hexagonal_events: Events with hexagonal patterns (default: 80)
        total_events: Total events analyzed (default: 80)
        null_probability: Expected rate under null (default: 0.5)

    Returns:
        SignificanceResult for the GWTC-3 analysis
    """
    result = calculate_binomial_significance(
        hexagonal_events,
        total_events,
        null_probability,
        alternative='greater'
    )

    # Add context
    print(f"\n{'='*60}")
    print("GWTC-3 HEXAGONAL DETECTION VALIDATION")
    print(f"{'='*60}")
    print(f"Observed: {hexagonal_events}/{total_events} events ({100*hexagonal_events/total_events:.1f}%)")
    print(f"Expected under null: {total_events * null_probability:.0f} events ({100*null_probability:.0f}%)")
    print(f"Z-score: {result.z_score:.1f}")
    print(f"Sigma significance: {result.sigma:.1f}σ")
    print(f"P-value: {result.p_value:.2e}")
    print(f"{'='*60}")

    return result
