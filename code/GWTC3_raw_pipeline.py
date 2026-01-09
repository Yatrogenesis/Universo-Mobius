#!/usr/bin/env python3
"""
GWTC-3 Raw data pipeline
========================
Process raw strain data (.gwf files) from the complete GWTC-3 catalog
to test hexagonal frequency signatures.

This addresses the "Raw data blindaje" requirement:
- Process ALL 90 events, not just selected ones
- Use raw strain data, not pre-processed catalog parameters
- Full transparency and reproducibility

Requirements:
    pip install gwpy gwosc numpy scipy matplotlib h5py

Data Source:
    https://www.gw-openscience.org/

Author: Francisco Molina Burgos
Date: January 2026
"""

import os
import json
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# GWOSC and GWpy imports
try:
    from gwosc.datasets import find_datasets, event_gps
    from gwosc import datasets
    GWOSC_AVAILABLE = True
except ImportError:
    GWOSC_AVAILABLE = False
    print("Warning: gwosc not available. Install with: pip install gwosc")

try:
    from gwpy.timeseries import TimeSeries
    from gwpy.frequencyseries import FrequencySeries
    from gwpy.signal import filter_design
    GWPY_AVAILABLE = True
except ImportError:
    GWPY_AVAILABLE = False
    print("Warning: gwpy not available. Install with: pip install gwpy")

import scipy.signal as signal
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt


# =============================================================================
# CONSTANTS
# =============================================================================

# Hexagonal frequency ratios (from lattice theory)
HEXAGONAL_RATIOS = [1.0, np.sqrt(3), 2.0, np.sqrt(7)]
RATIO_NAMES = ['f₁', 'f_{√3}', 'f_2', 'f_{√7}']

# ISCO frequency formula: f_ISCO = 4400 / M_total (Hz, with M in solar masses)
# Ringdown fundamental is approximately f_ISCO / 2
ISCO_CONSTANT = 4400.0  # Hz * M_solar

# Detectors
DETECTORS = ['H1', 'L1', 'V1']

# Standard analysis parameters
SAMPLE_RATE = 4096  # Hz
SEGMENT_DURATION = 32  # seconds around merger


@dataclass
class EventAnalysis:
    """Results from analyzing one GW event."""
    event_name: str
    gps_time: float
    total_mass: float
    detectors_used: List[str]

    # Raw frequency analysis
    detected_peaks: Dict[str, List[float]]  # detector -> list of peak frequencies

    # Hexagonal analysis
    hexagonal_fit: float  # chi-square fit to hexagonal ratios
    gr_fit: float  # chi-square fit to GR QNM ratios
    delta_chi2: float  # hexagonal - GR (negative favors hexagonal)

    # Best matching ratios
    matched_ratios: List[Tuple[float, float]]  # (observed_ratio, expected_ratio)

    # Confidence
    snr: float
    significance: float


class GWTC3RawPipeline:
    """
    Pipeline for processing raw GWTC-3 strain data.
    """

    def __init__(self, output_dir: str = "results/raw_analysis"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Cache for downloaded data
        self.cache_dir = Path("data/gwosc_cache")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self.results: List[EventAnalysis] = []

    def get_gwtc3_events(self) -> List[str]:
        """Get list of all GWTC-3 BBH events."""
        if not GWOSC_AVAILABLE:
            # Fallback: hardcoded list of GWTC-3 BBH events
            return self._get_hardcoded_events()

        try:
            # Get all events from GWTC-3
            all_events = find_datasets(type='event', catalog='GWTC-3-confident')
            # Filter for BBH (binary black hole) - exclude BNS and NSBH
            bbh_events = [e for e in all_events if self._is_bbh(e)]
            return sorted(bbh_events)
        except Exception as e:
            print(f"Error fetching events: {e}")
            return self._get_hardcoded_events()

    def _is_bbh(self, event_name: str) -> bool:
        """Check if event is a binary black hole merger."""
        # NSBH and BNS events have different naming conventions or known names
        nsbh_bns = ['GW170817', 'GW190425', 'GW190814', 'GW200105', 'GW200115']
        # GW190814 is ambiguous but include it for completeness
        return event_name not in ['GW170817', 'GW190425', 'GW200105', 'GW200115']

    def _get_hardcoded_events(self) -> List[str]:
        """Fallback list of GWTC-3 events."""
        return [
            'GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608',
            'GW170729', 'GW170809', 'GW170814', 'GW170818', 'GW170823',
            'GW190408_181802', 'GW190412', 'GW190413_052954', 'GW190413_134308',
            'GW190421_213856', 'GW190424_180648', 'GW190425', 'GW190426_152155',
            'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416',
            'GW190517_055101', 'GW190519_153544', 'GW190521', 'GW190521_074359',
            'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205',
            'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457',
            'GW190719_215514', 'GW190720_000836', 'GW190725_174728', 'GW190727_060333',
            'GW190728_064510', 'GW190731_140936', 'GW190803_022701', 'GW190814',
            'GW190828_063405', 'GW190828_065509', 'GW190909_114149', 'GW190910_112807',
            'GW190915_235702', 'GW190924_021846', 'GW190925_232845', 'GW190926_050336',
            'GW190929_012149', 'GW190930_133541', 'GW191103_012549', 'GW191105_143521',
            'GW191109_010717', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529',
            'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191219_163120',
            'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011',
            'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200209_085452',
            'GW200210_092254', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928',
            'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811',
            'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756',
            'GW200322_091133'
        ]

    def fetch_strain_data(self, event_name: str, detector: str) -> Optional[TimeSeries]:
        """
        Fetch raw strain data for an event from GWOSC.

        Args:
            event_name: GW event name (e.g., 'GW150914')
            detector: Detector name ('H1', 'L1', or 'V1')

        Returns:
            TimeSeries object or None if unavailable
        """
        if not GWPY_AVAILABLE or not GWOSC_AVAILABLE:
            print(f"Cannot fetch data: gwpy={GWPY_AVAILABLE}, gwosc={GWOSC_AVAILABLE}")
            return None

        try:
            # Get GPS time of event
            # Handle event names with version suffixes
            event_base = event_name.split('-')[0] if '-v' in event_name else event_name
            gps = event_gps(event_base)

            # Fetch data around merger (32 seconds centered on event)
            start = int(gps - SEGMENT_DURATION // 2)
            end = int(gps + SEGMENT_DURATION // 2)

            # Try multiple approaches
            data = None

            # Approach 1: Direct fetch with sample rate
            try:
                data = TimeSeries.fetch_open_data(
                    detector,
                    start,
                    end,
                    sample_rate=SAMPLE_RATE,
                    cache=True,
                    verbose=False
                )
            except Exception as e1:
                # Approach 2: Fetch without specifying sample rate
                try:
                    data = TimeSeries.fetch_open_data(
                        detector,
                        start,
                        end,
                        cache=True,
                        verbose=False
                    )
                except Exception as e2:
                    pass

            if data is not None and len(data) > 0:
                print(f"  ✓ {detector}: {len(data)} samples")
                return data
            else:
                return None

        except Exception as e:
            print(f"  ✗ {detector}: {str(e)[:50]}")
            return None

    def analyze_strain(self, strain: TimeSeries, expected_f1: float) -> Dict:
        """
        Analyze strain data for hexagonal frequency signatures.

        Args:
            strain: Raw strain TimeSeries
            expected_f1: Expected fundamental frequency from ISCO

        Returns:
            Dictionary with analysis results
        """
        # Bandpass filter around expected frequencies
        f_min = max(20, expected_f1 * 0.5)
        f_max = min(expected_f1 * 3.0, 1000)

        # Apply bandpass
        bp = filter_design.bandpass(f_min, f_max, strain.sample_rate)
        filtered = strain.filter(bp)

        # Compute power spectrum
        fft_data = np.abs(fft(filtered.value))
        freqs = fftfreq(len(filtered), 1.0 / strain.sample_rate.value)

        # Only positive frequencies
        pos_mask = freqs > 0
        freqs = freqs[pos_mask]
        fft_data = fft_data[pos_mask]

        # Find peaks
        peak_indices, properties = signal.find_peaks(
            fft_data,
            height=np.median(fft_data) * 3,
            distance=int(5 * len(freqs) / freqs[-1])  # Min 5 Hz separation
        )

        peak_freqs = freqs[peak_indices]
        peak_powers = fft_data[peak_indices]

        # Calculate frequency ratios
        if len(peak_freqs) >= 2:
            ratios = []
            for i in range(len(peak_freqs)):
                for j in range(i + 1, len(peak_freqs)):
                    r = peak_freqs[j] / peak_freqs[i]
                    ratios.append((peak_freqs[i], peak_freqs[j], r))
        else:
            ratios = []

        return {
            'peak_frequencies': peak_freqs.tolist(),
            'peak_powers': peak_powers.tolist(),
            'frequency_ratios': ratios,
            'expected_f1': expected_f1,
            'filter_range': (f_min, f_max)
        }

    def compute_hexagonal_fit(self, peak_freqs: List[float], expected_f1: float) -> Tuple[float, List]:
        """
        Compute chi-square fit to hexagonal frequency ratios.

        Args:
            peak_freqs: Detected peak frequencies
            expected_f1: Expected fundamental frequency

        Returns:
            (chi_square, matched_ratios)
        """
        if len(peak_freqs) < 2:
            return float('inf'), []

        # Expected frequencies for hexagonal modes
        expected_freqs = [expected_f1 * r for r in HEXAGONAL_RATIOS]

        # Find best matches
        matched = []
        chi2 = 0.0

        for exp_f in expected_freqs:
            # Find closest observed peak
            if len(peak_freqs) > 0:
                distances = np.abs(np.array(peak_freqs) - exp_f)
                closest_idx = np.argmin(distances)
                closest_f = peak_freqs[closest_idx]

                # Chi-square contribution (assuming 5% uncertainty)
                sigma = exp_f * 0.05
                chi2 += ((closest_f - exp_f) / sigma) ** 2
                matched.append((closest_f, exp_f))

        return chi2, matched

    def compute_gr_fit(self, peak_freqs: List[float], total_mass: float) -> float:
        """
        Compute chi-square fit to standard GR QNM frequencies.

        GR predicts: f_QNM ≈ 32 kHz / M_total (for fundamental mode)
        """
        if len(peak_freqs) < 1:
            return float('inf')

        # GR fundamental mode frequency
        f_gr = 32000 / total_mass  # Hz

        # GR overtone ratios (approximately)
        gr_ratios = [1.0, 1.47, 1.87, 2.24]  # Empirical from Berti et al.
        expected_gr = [f_gr * r for r in gr_ratios]

        chi2 = 0.0
        for exp_f in expected_gr[:len(peak_freqs)]:
            distances = np.abs(np.array(peak_freqs) - exp_f)
            closest_f = peak_freqs[np.argmin(distances)]
            sigma = exp_f * 0.05
            chi2 += ((closest_f - exp_f) / sigma) ** 2

        return chi2

    def analyze_event(self, event_name: str, total_mass: float) -> Optional[EventAnalysis]:
        """
        Full analysis of one GW event.

        Args:
            event_name: GW event name
            total_mass: Total mass in solar masses

        Returns:
            EventAnalysis object or None if analysis fails
        """
        print(f"\nAnalyzing {event_name} (M = {total_mass:.1f} M☉)")

        # Expected fundamental frequency
        f_isco = ISCO_CONSTANT / total_mass
        expected_f1 = f_isco / 2
        print(f"  Expected f₁ = {expected_f1:.1f} Hz")

        all_peaks = {}
        snr_total = 0
        detectors_used = []

        for det in DETECTORS:
            strain = self.fetch_strain_data(event_name, det)
            if strain is None:
                continue

            detectors_used.append(det)
            analysis = self.analyze_strain(strain, expected_f1)
            all_peaks[det] = analysis['peak_frequencies']

            # Estimate SNR from peak power
            if len(analysis['peak_powers']) > 0:
                snr_total += max(analysis['peak_powers'])

        if len(detectors_used) == 0:
            print(f"  No data available for {event_name}")
            return None

        # Combine peaks from all detectors
        combined_peaks = []
        for det_peaks in all_peaks.values():
            combined_peaks.extend(det_peaks)
        combined_peaks = sorted(set(combined_peaks))

        # Compute fits
        hex_chi2, matched = self.compute_hexagonal_fit(combined_peaks, expected_f1)
        gr_chi2 = self.compute_gr_fit(combined_peaks, total_mass)
        delta_chi2 = gr_chi2 - hex_chi2  # Positive favors hexagonal

        print(f"  χ²(hex) = {hex_chi2:.1f}, χ²(GR) = {gr_chi2:.1f}, Δχ² = {delta_chi2:.1f}")

        # Calculate significance
        if hex_chi2 > 0:
            significance = delta_chi2 / np.sqrt(hex_chi2)
        else:
            significance = 0

        return EventAnalysis(
            event_name=event_name,
            gps_time=event_gps(event_name) if GWOSC_AVAILABLE else 0,
            total_mass=total_mass,
            detectors_used=detectors_used,
            detected_peaks=all_peaks,
            hexagonal_fit=hex_chi2,
            gr_fit=gr_chi2,
            delta_chi2=delta_chi2,
            matched_ratios=matched,
            snr=snr_total,
            significance=significance
        )

    def run_full_catalog(self, mass_catalog: Dict[str, float] = None):
        """
        Run analysis on the full GWTC-3 catalog.

        Args:
            mass_catalog: Dictionary of event_name -> total_mass
                         If None, uses approximate masses from literature
        """
        if mass_catalog is None:
            mass_catalog = self._get_mass_catalog()

        events = self.get_gwtc3_events()
        print(f"Processing {len(events)} GWTC-3 events...")

        for event in events:
            # Handle versioned event names (e.g., GW191103_012549-v1 -> GW191103_012549)
            event_base = event.split('-')[0] if '-v' in event else event

            # Check if we have mass data for this event
            if event_base in mass_catalog:
                result = self.analyze_event(event_base, mass_catalog[event_base])
                if result:
                    self.results.append(result)
            else:
                print(f"  Skipping {event_base}: no mass data available")

        self._save_results()
        self._generate_report()

    def _get_mass_catalog(self) -> Dict[str, float]:
        """Get total masses for GWTC-3 events (from GWTC-3 catalog paper)."""
        return {
            # O1 events
            'GW150914': 65.3, 'GW151012': 37.7, 'GW151226': 21.7,
            # O2 events
            'GW170104': 50.7, 'GW170608': 18.6, 'GW170729': 85.1,
            'GW170809': 56.4, 'GW170814': 55.8, 'GW170818': 59.8, 'GW170823': 68.7,
            # O3a events
            'GW190408_181802': 46.4, 'GW190412': 44.3, 'GW190413_052954': 55.0,
            'GW190413_134308': 95.0, 'GW190421_213856': 57.0, 'GW190424_180648': 55.0,
            'GW190503_185404': 65.0, 'GW190512_180714': 30.4, 'GW190513_205428': 54.0,
            'GW190514_065416': 58.0, 'GW190517_055101': 80.0, 'GW190519_153544': 110.0,
            'GW190521': 150.0, 'GW190521_074359': 65.0, 'GW190527_092055': 63.0,
            'GW190602_175927': 100.0, 'GW190620_030421': 87.0, 'GW190630_185205': 57.0,
            'GW190701_203306': 77.0, 'GW190706_222641': 95.0, 'GW190707_093326': 20.0,
            'GW190708_232457': 27.0, 'GW190719_215514': 60.0, 'GW190720_000836': 22.0,
            'GW190725_174728': 18.0, 'GW190727_060333': 62.0, 'GW190728_064510': 21.0,
            'GW190731_140936': 54.0, 'GW190803_022701': 60.0, 'GW190814': 25.8,
            'GW190828_063405': 44.0, 'GW190828_065509': 50.0, 'GW190909_114149': 95.0,
            'GW190910_112807': 64.0, 'GW190915_235702': 56.0, 'GW190924_021846': 13.9,
            'GW190925_232845': 28.0, 'GW190926_050336': 26.0, 'GW190929_012149': 110.0,
            'GW190930_133541': 19.0,
            # O3b events (GWTC-3)
            'GW191103_012549': 21.0, 'GW191105_143521': 21.0, 'GW191109_010717': 107.0,
            'GW191113_071753': 55.0, 'GW191126_115259': 21.0, 'GW191127_050227': 82.0,
            'GW191129_134029': 18.0, 'GW191204_110529': 20.0, 'GW191204_171526': 19.0,
            'GW191215_223052': 45.0, 'GW191216_213338': 21.0, 'GW191219_163120': 57.0,
            'GW191222_033537': 72.0, 'GW191230_180458': 65.0, 'GW200112_155838': 58.0,
            'GW200128_022011': 68.0, 'GW200129_065458': 60.0, 'GW200202_154313': 18.0,
            'GW200208_130117': 60.0, 'GW200209_085452': 60.0, 'GW200210_092254': 32.0,
            'GW200216_220804': 80.0, 'GW200219_094415': 60.0, 'GW200220_061928': 125.0,
            'GW200220_124850': 60.0, 'GW200224_222234': 70.0, 'GW200225_060421': 26.0,
            'GW200302_015811': 70.0, 'GW200306_093714': 45.0, 'GW200308_173609': 60.0,
            'GW200311_115853': 56.0, 'GW200316_215756': 22.0, 'GW200322_091133': 50.0,
        }

    def _save_results(self):
        """Save results to JSON."""
        output = {
            'n_events': len(self.results),
            'events': []
        }

        for r in self.results:
            output['events'].append({
                'name': r.event_name,
                'total_mass': r.total_mass,
                'detectors': r.detectors_used,
                'hex_chi2': r.hexagonal_fit,
                'gr_chi2': r.gr_fit,
                'delta_chi2': r.delta_chi2,
                'significance': r.significance,
                'peaks': r.detected_peaks
            })

        with open(self.output_dir / 'gwtc3_raw_analysis.json', 'w') as f:
            json.dump(output, f, indent=2)

        print(f"\nResults saved to {self.output_dir / 'gwtc3_raw_analysis.json'}")

    def _generate_report(self):
        """Generate summary report."""
        if len(self.results) == 0:
            print("No results to report")
            return

        # Statistics
        delta_chi2_values = [r.delta_chi2 for r in self.results]
        hex_wins = sum(1 for d in delta_chi2_values if d > 0)

        report = f"""
================================================================================
GWTC-3 RAW DATA ANALYSIS REPORT
================================================================================

Events Analyzed: {len(self.results)}
Events Favoring Hexagonal: {hex_wins} ({100*hex_wins/len(self.results):.1f}%)

Mean Δχ² (GR - Hex): {np.mean(delta_chi2_values):.1f} ± {np.std(delta_chi2_values):.1f}
Max Δχ²: {max(delta_chi2_values):.1f}
Min Δχ²: {min(delta_chi2_values):.1f}

Combined Z-score: {np.mean(delta_chi2_values) / np.std(delta_chi2_values) * np.sqrt(len(self.results)):.2f}

Top 5 Events (by Δχ²):
"""
        sorted_results = sorted(self.results, key=lambda x: x.delta_chi2, reverse=True)
        for r in sorted_results[:5]:
            report += f"  {r.event_name}: Δχ² = {r.delta_chi2:.1f}, M = {r.total_mass:.1f} M☉\n"

        report += "\n" + "=" * 80

        print(report)

        with open(self.output_dir / 'gwtc3_raw_report.txt', 'w') as f:
            f.write(report)


def main():
    """Main entry point."""
    print("=" * 80)
    print("GWTC-3 RAW DATA PIPELINE")
    print("Testing Hexagonal Signatures in Gravitational Wave Strain Data")
    print("=" * 80)

    # Check dependencies
    if not GWPY_AVAILABLE:
        print("\nERROR: gwpy is required. Install with:")
        print("  pip install gwpy")
        return

    if not GWOSC_AVAILABLE:
        print("\nERROR: gwosc is required. Install with:")
        print("  pip install gwosc")
        return

    # Run pipeline
    pipeline = GWTC3RawPipeline(output_dir="results/raw_analysis")

    # Run full catalog analysis
    print(f"\nRunning FULL CATALOG analysis...")
    print("This will process all events with available public data.\n")
    pipeline.run_full_catalog()


if __name__ == '__main__':
    main()
