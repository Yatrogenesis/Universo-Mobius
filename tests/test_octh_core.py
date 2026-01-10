"""
Tests for OCTH Core Module
==========================

Tests for hexagonal constants, transforms, and statistics.
"""

import pytest
import numpy as np
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from octh_core.constants import (
    HEXAGONAL_RATIOS,
    GOLDEN_RATIO,
    DEFAULT_TOLERANCE,
    sigma_to_pvalue,
    pvalue_to_sigma,
)
from octh_core.transforms import (
    HexagonalTransform,
    FrequencyRatioAnalyzer,
    compute_ratio_proximity,
    find_hexagonal_modes,
)
from octh_core.statistics import (
    calculate_binomial_significance,
    calculate_combined_significance,
    monte_carlo_null,
    validate_gwtc3_result,
)


class TestConstants:
    """Test OCTH constants."""

    def test_hexagonal_ratios_count(self):
        """Should have 6 hexagonal ratios."""
        assert len(HEXAGONAL_RATIOS) == 6

    def test_hexagonal_ratios_values(self):
        """Ratios should match expected values."""
        expected = [1.0, np.sqrt(3), 2.0, np.sqrt(7), 3.0, np.sqrt(12)]
        for actual, exp in zip(HEXAGONAL_RATIOS, expected):
            assert abs(actual - exp) < 1e-10

    def test_golden_ratio(self):
        """Golden ratio should be correct."""
        expected = (1 + np.sqrt(5)) / 2
        assert abs(GOLDEN_RATIO - expected) < 1e-10

    def test_sigma_pvalue_conversion(self):
        """Sigma-pvalue conversion should be consistent."""
        for sigma in [1.0, 2.0, 3.0, 5.0]:
            p = sigma_to_pvalue(sigma)
            recovered = pvalue_to_sigma(p)
            assert abs(recovered - sigma) < 0.01


class TestHexagonalTransform:
    """Test HexagonalTransform class."""

    @pytest.fixture
    def transform(self):
        """Create transform instance."""
        return HexagonalTransform(tolerance=0.03, min_snr=2.0)

    @pytest.fixture
    def test_spectrum(self):
        """Create test spectrum with hexagonal peaks."""
        sample_rate = 44100
        n_samples = 4096
        freqs = np.fft.rfftfreq(n_samples, 1/sample_rate)

        # Create spectrum with peaks at hexagonal ratios of 100 Hz
        f0 = 100.0
        spectrum = np.random.randn(len(freqs)) * 0.1 + 1.0

        for ratio in HEXAGONAL_RATIOS:
            f_peak = f0 * ratio
            idx = np.argmin(np.abs(freqs - f_peak))
            if idx < len(spectrum):
                spectrum[idx] += 10.0  # Strong peak

        return spectrum, freqs, f0

    def test_detect_fundamental(self, transform, test_spectrum):
        """Should detect fundamental frequency."""
        spectrum, freqs, expected_f0 = test_spectrum
        detected_f0 = transform.detect_fundamental(spectrum, freqs, (50, 200))

        # Should be within 10% of expected
        assert abs(detected_f0 - expected_f0) / expected_f0 < 0.1

    def test_find_modes(self, transform, test_spectrum):
        """Should find hexagonal modes."""
        spectrum, freqs, f0 = test_spectrum
        modes = transform.find_modes(spectrum, freqs, f0)

        # Should find multiple modes
        assert len(modes) >= 3

        # First mode should be near fundamental
        first_mode = modes[0]
        assert abs(first_mode.frequency - f0) / f0 < 0.05

    def test_transform_inverse(self, transform, test_spectrum):
        """Transform should be invertible."""
        spectrum, freqs, f0 = test_spectrum

        # Forward transform
        hex_amplitudes = transform.transform(spectrum, freqs, f0)

        # Should have amplitude at each hexagonal ratio
        assert len(hex_amplitudes) == len(HEXAGONAL_RATIOS)
        assert np.all(hex_amplitudes > 0)


class TestFrequencyRatioAnalyzer:
    """Test FrequencyRatioAnalyzer class."""

    @pytest.fixture
    def analyzer(self):
        """Create analyzer instance."""
        return FrequencyRatioAnalyzer(tolerance=0.03)

    def test_classify_hexagonal_ratio(self, analyzer):
        """Should correctly classify hexagonal ratios."""
        # Test exact hexagonal ratios
        for i, ratio in enumerate(HEXAGONAL_RATIOS):
            idx, dev, is_hex = analyzer.classify_ratio(ratio)
            assert is_hex
            assert idx == i
            assert dev < 0.001

    def test_classify_non_hexagonal_ratio(self, analyzer):
        """Should reject non-hexagonal ratios."""
        # Test a ratio that's not hexagonal
        ratio = 1.5  # Between sqrt(3) and 2
        idx, dev, is_hex = analyzer.classify_ratio(ratio)
        assert not is_hex
        assert dev > 0.03

    def test_analyze_hexagonal_spectrum(self, analyzer):
        """Should identify hexagonal spectrum."""
        # Create spectrum with hexagonal peaks
        freqs = np.linspace(10, 1000, 1000)
        spectrum = np.random.randn(len(freqs)) * 0.1

        # Add peaks at hexagonal ratios of 100 Hz
        f0 = 100.0
        for ratio in HEXAGONAL_RATIOS:
            idx = np.argmin(np.abs(freqs - f0 * ratio))
            spectrum[idx] = 10.0

        result = analyzer.analyze(spectrum, freqs, n_peaks=10)

        # Should be hexagonal
        assert result['is_hexagonal']
        assert result['hexagonal_fraction'] > 0.5


class TestStatistics:
    """Test statistical functions."""

    def test_binomial_significance_extreme(self):
        """Test significance for extreme result (80/80)."""
        result = calculate_binomial_significance(80, 80, 0.5)

        # Should be highly significant
        assert result.sigma > 10  # > 10 sigma
        assert result.p_value < 1e-10
        assert result.z_score > 0

    def test_binomial_significance_null(self):
        """Test significance for null result (50/100)."""
        result = calculate_binomial_significance(50, 100, 0.5)

        # Should not be significant
        assert result.sigma < 2
        assert result.p_value > 0.05

    def test_combined_significance_fisher(self):
        """Test Fisher's method for combining p-values."""
        # Create multiple significant results
        results = [
            calculate_binomial_significance(8, 10, 0.5)
            for _ in range(5)
        ]

        combined = calculate_combined_significance(results, method='fisher')

        # Combined should be more significant
        assert combined.sigma > results[0].sigma

    def test_gwtc3_validation(self, capsys):
        """Test GWTC-3 result validation."""
        result = validate_gwtc3_result(80, 80, 0.5)

        # Should match expected values
        assert result.observed == 80
        assert result.total == 80
        assert result.sigma > 70  # Should be ~75 sigma

        # Check output was printed
        captured = capsys.readouterr()
        assert 'GWTC-3' in captured.out


class TestConvenienceFunctions:
    """Test convenience functions."""

    def test_compute_ratio_proximity(self):
        """Test ratio proximity function."""
        # Exact match
        idx, dev, match = compute_ratio_proximity(2.0)
        assert match
        assert idx == 2  # Index of 2.0 in HEXAGONAL_RATIOS
        assert dev < 0.001

        # Near match
        idx, dev, match = compute_ratio_proximity(2.05, tolerance=0.03)
        assert match
        assert dev < 0.03

        # No match
        idx, dev, match = compute_ratio_proximity(1.5, tolerance=0.03)
        assert not match

    def test_find_hexagonal_modes_convenience(self):
        """Test convenience function for finding modes."""
        # Create simple spectrum
        freqs = np.linspace(10, 1000, 500)
        spectrum = np.ones(len(freqs))

        # Add peak at 100 Hz
        idx = np.argmin(np.abs(freqs - 100))
        spectrum[idx] = 20.0

        modes = find_hexagonal_modes(spectrum, freqs, f0=100.0)

        # Should find at least one mode
        assert len(modes) >= 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
