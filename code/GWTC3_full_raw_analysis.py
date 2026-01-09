#!/usr/bin/env python3
"""
GWTC-3 Full Raw Analysis Pipeline (No gwpy dependency)
======================================================
Downloads and analyzes ALL GWTC-3 events using GWOSC data.
Uses h5py + scipy instead of gwpy to avoid C compilation issues.

Author: Francisco Molina Burgos & Claude
Date: January 2025
"""

import os
import sys
import json
import requests
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, asdict
import h5py
from scipy import signal
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import warnings
warnings.filterwarnings('ignore')

# GWOSC imports
from gwosc.datasets import find_datasets, event_gps
from gwosc.locate import get_urls

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "gwosc_raw"
RESULTS_DIR = BASE_DIR / "results" / "full_raw_analysis"
FIGURES_DIR = BASE_DIR / "figures" / "full_raw"

DATA_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Analysis parameters
SAMPLE_RATE = 4096  # Hz
SEGMENT_DURATION = 32  # seconds around merger
DETECTORS = ['H1', 'L1', 'V1']

# Hexagonal frequency ratios
HEXAGONAL_RATIOS = np.array([1.0, np.sqrt(3), 2.0, np.sqrt(7)])
RATIO_NAMES = ['f1', 'f_sqrt3', 'f_2', 'f_sqrt7']

# ISCO constant
ISCO_CONSTANT = 4400.0  # Hz * M_solar


# =============================================================================
# GWTC-3 MASS CATALOG (from LIGO publications)
# =============================================================================

MASS_CATALOG = {
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
    # O3b events
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


@dataclass
class EventResult:
    """Result from analyzing one event."""
    name: str
    total_mass: float
    gps_time: float
    detectors_used: List[str]
    hex_chi2: float
    gr_chi2: float
    delta_chi2: float
    significance: float
    favors_hexagonal: bool
    n_peaks: Dict[str, int]
    peak_frequencies: Dict[str, List[float]]
    status: str  # 'success', 'partial', 'failed'
    error: Optional[str] = None


class GWTC3FullAnalyzer:
    """Full GWTC-3 raw data analyzer."""

    def __init__(self):
        self.results: List[EventResult] = []
        self.failed_downloads: List[str] = []

    def get_all_events(self) -> List[str]:
        """Get list of all GWTC-3 BBH events with mass data."""
        return list(MASS_CATALOG.keys())

    def download_event_data(self, event_name: str) -> Dict[str, Optional[Path]]:
        """
        Download strain data for an event from GWOSC.

        Returns dict of detector -> file path (or None if unavailable)
        """
        files = {}

        # Try different catalogs in order
        catalogs = [
            'GWTC-3-confident',
            'GWTC-2.1-confident',
            'GWTC-1-confident',
            'O3_Discovery_Papers',
            'O2_4KHZ_R1',
            'O1_4KHZ_R1'
        ]

        event_data = None
        for catalog in catalogs:
            try:
                api_url = f"https://gwosc.org/eventapi/json/{catalog}/{event_name}/"
                resp = requests.get(api_url, timeout=30)
                if resp.status_code == 200:
                    data = resp.json()
                    if 'events' in data:
                        # Find the event key (may have version suffix)
                        for key in data['events']:
                            if event_name in key:
                                event_data = data['events'][key]
                                break
                    if event_data:
                        break
            except:
                continue

        if not event_data:
            print(f"    [!] Event {event_name} not found in any catalog")
            return files

        strain_list = event_data.get('strain', [])
        if not strain_list:
            print(f"    [!] No strain data for {event_name}")
            return files

        for detector in DETECTORS:
            cache_file = DATA_DIR / f"{event_name}_{detector}.hdf5"

            # Check cache
            if cache_file.exists() and cache_file.stat().st_size > 1000:
                files[detector] = cache_file
                continue

            # Find 4KHz HDF5 file for this detector
            strain_url = None
            for s in strain_list:
                if (s.get('detector') == detector and
                    s.get('sampling_rate') == 4096 and
                    s.get('url', '').endswith('.hdf5')):
                    strain_url = s['url']
                    break

            if not strain_url:
                continue

            try:
                print(f"    Downloading {detector}...", end=" ", flush=True)
                response = requests.get(strain_url, timeout=120)
                response.raise_for_status()

                with open(cache_file, 'wb') as f:
                    f.write(response.content)

                print(f"{len(response.content)//1024}KB")
                files[detector] = cache_file

            except Exception as e:
                print(f"FAILED: {e}")

        return files

    def read_strain_data(self, filepath: Path) -> Tuple[Optional[np.ndarray], float]:
        """
        Read strain data from HDF5 file.

        Returns (strain_array, sample_rate) or (None, 0) if failed
        """
        try:
            with h5py.File(filepath, 'r') as f:
                # GWOSC HDF5 structure
                if 'strain' in f:
                    strain_group = f['strain']
                    if 'Strain' in strain_group:
                        strain = strain_group['Strain'][:]
                        # Get sample rate from attributes
                        dt = strain_group['Strain'].attrs.get('Xspacing', 1/4096)
                        sample_rate = 1.0 / dt
                        return strain, sample_rate

                # Alternative structure
                for key in f.keys():
                    if 'strain' in key.lower():
                        data = f[key][:]
                        return data, SAMPLE_RATE

            return None, 0

        except Exception as e:
            return None, 0

    def analyze_strain(self, strain: np.ndarray, sample_rate: float,
                      expected_f1: float) -> Dict:
        """
        Analyze strain data for hexagonal frequency signatures.
        """
        # Bandpass filter
        f_min = max(20, expected_f1 * 0.3)
        f_max = min(expected_f1 * 4.0, sample_rate / 2 - 10)

        # Design bandpass filter
        nyq = sample_rate / 2
        low = f_min / nyq
        high = f_max / nyq

        try:
            b, a = signal.butter(4, [low, high], btype='band')
            filtered = signal.filtfilt(b, a, strain)
        except:
            filtered = strain

        # Compute power spectrum
        n = len(filtered)
        fft_data = np.abs(fft(filtered))[:n//2]
        freqs = fftfreq(n, 1/sample_rate)[:n//2]

        # Find peaks
        # Normalize spectrum
        fft_norm = fft_data / np.max(fft_data)

        # Find peaks above threshold
        peak_indices, properties = signal.find_peaks(
            fft_norm,
            height=0.05,  # 5% of max
            distance=int(5 * n / sample_rate),  # Min 5 Hz separation
            prominence=0.02
        )

        # Filter to relevant frequency range
        mask = (freqs[peak_indices] >= f_min) & (freqs[peak_indices] <= f_max)
        peak_freqs = freqs[peak_indices][mask]
        peak_powers = fft_data[peak_indices][mask]

        # Sort by power
        sort_idx = np.argsort(peak_powers)[::-1]
        peak_freqs = peak_freqs[sort_idx][:50]  # Top 50 peaks

        return {
            'peak_frequencies': peak_freqs.tolist(),
            'n_peaks': len(peak_freqs),
            'expected_f1': expected_f1
        }

    def compute_hexagonal_chi2(self, peak_freqs: List[float],
                               expected_f1: float) -> float:
        """Compute chi-square fit to hexagonal ratios."""
        if len(peak_freqs) < 2:
            return float('inf')

        expected = expected_f1 * HEXAGONAL_RATIOS
        chi2 = 0.0

        for exp_f in expected:
            if len(peak_freqs) > 0:
                distances = np.abs(np.array(peak_freqs) - exp_f)
                closest = np.min(distances)
                sigma = exp_f * 0.1  # 10% tolerance
                chi2 += (closest / sigma) ** 2

        return chi2

    def compute_gr_chi2(self, peak_freqs: List[float], total_mass: float) -> float:
        """Compute chi-square fit to GR QNM predictions."""
        if len(peak_freqs) < 2:
            return float('inf')

        # GR fundamental QNM frequency (Berti et al.)
        f_qnm = 32000 / total_mass
        gr_ratios = np.array([1.0, 1.47, 1.87, 2.24])  # QNM overtones
        expected = f_qnm * gr_ratios

        chi2 = 0.0
        for exp_f in expected:
            if len(peak_freqs) > 0:
                distances = np.abs(np.array(peak_freqs) - exp_f)
                closest = np.min(distances)
                sigma = exp_f * 0.1
                chi2 += (closest / sigma) ** 2

        return chi2

    def analyze_event(self, event_name: str) -> Optional[EventResult]:
        """Full analysis of one event."""
        if event_name not in MASS_CATALOG:
            return None

        total_mass = MASS_CATALOG[event_name]
        expected_f1 = (ISCO_CONSTANT / total_mass) / 2

        print(f"\n[{event_name}] M={total_mass:.1f} Msol, f1_exp={expected_f1:.1f} Hz")

        # Download data
        files = self.download_event_data(event_name)

        if not files:
            print(f"    [FAIL] No data available")
            return EventResult(
                name=event_name,
                total_mass=total_mass,
                gps_time=0,
                detectors_used=[],
                hex_chi2=float('inf'),
                gr_chi2=float('inf'),
                delta_chi2=0,
                significance=0,
                favors_hexagonal=False,
                n_peaks={},
                peak_frequencies={},
                status='failed',
                error='No data available'
            )

        # Analyze each detector
        all_peaks = {}
        n_peaks = {}
        detectors_used = []

        try:
            gps = event_gps(event_name)
        except:
            gps = 0

        for detector, filepath in files.items():
            strain, sr = self.read_strain_data(filepath)

            if strain is None:
                print(f"    [{detector}] Failed to read data")
                continue

            print(f"    [{detector}] {len(strain)} samples @ {sr:.0f} Hz")

            analysis = self.analyze_strain(strain, sr, expected_f1)
            all_peaks[detector] = analysis['peak_frequencies']
            n_peaks[detector] = analysis['n_peaks']
            detectors_used.append(detector)

        if not detectors_used:
            return EventResult(
                name=event_name,
                total_mass=total_mass,
                gps_time=gps,
                detectors_used=[],
                hex_chi2=float('inf'),
                gr_chi2=float('inf'),
                delta_chi2=0,
                significance=0,
                favors_hexagonal=False,
                n_peaks={},
                peak_frequencies={},
                status='failed',
                error='Could not read any detector data'
            )

        # Combine peaks from all detectors
        combined_peaks = []
        for peaks in all_peaks.values():
            combined_peaks.extend(peaks)
        combined_peaks = sorted(set(combined_peaks))

        # Compute fits
        hex_chi2 = self.compute_hexagonal_chi2(combined_peaks, expected_f1)
        gr_chi2 = self.compute_gr_chi2(combined_peaks, total_mass)
        delta_chi2 = gr_chi2 - hex_chi2  # Positive = hexagonal better

        # Significance
        if hex_chi2 > 0 and hex_chi2 != float('inf'):
            significance = delta_chi2 / np.sqrt(hex_chi2)
        else:
            significance = 0

        favors_hex = delta_chi2 > 0

        status = 'success' if len(detectors_used) >= 2 else 'partial'

        print(f"    chi2_hex={hex_chi2:.1f}, chi2_GR={gr_chi2:.1f}, Delta_chi2={delta_chi2:.1f}")
        print(f"    {'--> HEXAGONAL' if favors_hex else '--> GR'}")

        return EventResult(
            name=event_name,
            total_mass=total_mass,
            gps_time=gps,
            detectors_used=detectors_used,
            hex_chi2=hex_chi2,
            gr_chi2=gr_chi2,
            delta_chi2=delta_chi2,
            significance=significance,
            favors_hexagonal=favors_hex,
            n_peaks=n_peaks,
            peak_frequencies=all_peaks,
            status=status
        )

    def run_full_analysis(self):
        """Run analysis on all events."""
        events = self.get_all_events()
        total = len(events)

        print("=" * 70)
        print(f"GWTC-3 FULL RAW ANALYSIS: {total} events")
        print("=" * 70)

        for i, event in enumerate(events, 1):
            print(f"\n[{i}/{total}]", end="")
            result = self.analyze_event(event)
            if result:
                self.results.append(result)

            # Save intermediate results every 10 events
            if i % 10 == 0:
                self._save_results()

        self._save_results()
        self._generate_report()
        self._generate_figures()

    def _save_results(self):
        """Save results to JSON."""
        output = {
            'n_events': len(self.results),
            'n_success': sum(1 for r in self.results if r.status == 'success'),
            'n_partial': sum(1 for r in self.results if r.status == 'partial'),
            'n_failed': sum(1 for r in self.results if r.status == 'failed'),
            'events': [asdict(r) for r in self.results]
        }

        with open(RESULTS_DIR / 'gwtc3_full_raw.json', 'w') as f:
            json.dump(output, f, indent=2, default=str)

    def _generate_report(self):
        """Generate summary report."""
        successful = [r for r in self.results if r.status in ['success', 'partial']]

        if not successful:
            print("\nNo successful analyses to report.")
            return

        hex_wins = sum(1 for r in successful if r.favors_hexagonal)
        delta_values = [r.delta_chi2 for r in successful if r.delta_chi2 != float('inf')]

        report = f"""
================================================================================
GWTC-3 FULL RAW DATA ANALYSIS REPORT
================================================================================

SUMMARY
-------
Total Events Attempted: {len(self.results)}
Successful Analyses: {len(successful)}
  - Full (2+ detectors): {sum(1 for r in successful if r.status == 'success')}
  - Partial (1 detector): {sum(1 for r in successful if r.status == 'partial')}
Failed: {sum(1 for r in self.results if r.status == 'failed')}

RESULTS
-------
Events Favoring Hexagonal: {hex_wins}/{len(successful)} ({100*hex_wins/len(successful):.1f}%)
Events Favoring GR: {len(successful) - hex_wins}/{len(successful)} ({100*(len(successful)-hex_wins)/len(successful):.1f}%)

Mean Delta_chi2 (GR - Hex): {np.mean(delta_values):.1f} ± {np.std(delta_values):.1f}
Median Delta_chi2: {np.median(delta_values):.1f}
Max Delta_chi2: {max(delta_values):.1f}
Min Delta_chi2: {min(delta_values):.1f}

Combined Z-score: {np.mean(delta_values) / (np.std(delta_values) / np.sqrt(len(delta_values))):.2f}

BY OBSERVING RUN
----------------
"""
        # Group by run
        o1 = [r for r in successful if r.name.startswith('GW15')]
        o2 = [r for r in successful if r.name.startswith('GW17')]
        o3a = [r for r in successful if r.name.startswith('GW19') and int(r.name[4:6]) < 11]
        o3b = [r for r in successful if r.name.startswith('GW19') and int(r.name[4:6]) >= 11 or r.name.startswith('GW20')]

        for name, group in [('O1', o1), ('O2', o2), ('O3a', o3a), ('O3b', o3b)]:
            if group:
                hex_w = sum(1 for r in group if r.favors_hexagonal)
                report += f"{name}: {len(group)} events, {hex_w} favor hex ({100*hex_w/len(group):.0f}%)\n"

        report += f"""
TOP 10 EVENTS (by Delta_chi2)
----------------------
"""
        sorted_results = sorted(successful, key=lambda x: x.delta_chi2 if x.delta_chi2 != float('inf') else -999, reverse=True)
        for r in sorted_results[:10]:
            report += f"  {r.name}: Delta_chi2 = {r.delta_chi2:.1f}, M = {r.total_mass:.1f} Msol, {r.status}\n"

        report += f"""
KEY EVENTS
----------
"""
        key_events = ['GW150914', 'GW170817', 'GW190521', 'GW190814', 'GW170814']
        for name in key_events:
            r = next((x for x in self.results if x.name == name), None)
            if r:
                status = "HEXAGONAL" if r.favors_hexagonal else "GR"
                report += f"  {name}: {status}, Delta_chi2 = {r.delta_chi2:.1f}\n"
            else:
                report += f"  {name}: NOT ANALYZED\n"

        report += "\n" + "=" * 80

        print(report)

        with open(RESULTS_DIR / 'gwtc3_full_report.txt', 'w') as f:
            f.write(report)

    def _generate_figures(self):
        """Generate summary figures."""
        successful = [r for r in self.results if r.status in ['success', 'partial'] and r.delta_chi2 != float('inf')]

        if len(successful) < 5:
            print("Not enough data for figures")
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # 1. Delta_chi2 distribution
        ax = axes[0, 0]
        delta_values = [r.delta_chi2 for r in successful]
        ax.hist(delta_values, bins=30, edgecolor='black', alpha=0.7)
        ax.axvline(0, color='red', linestyle='--', label='Delta_chi2=0 (neutral)')
        ax.axvline(np.mean(delta_values), color='green', linestyle='-', label=f'Mean={np.mean(delta_values):.1f}')
        ax.set_xlabel('Delta_chi2 (GR - Hexagonal)')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of Delta_chi2 (Positive = Favors Hexagonal)')
        ax.legend()

        # 2. Delta_chi2 vs Mass
        ax = axes[0, 1]
        masses = [r.total_mass for r in successful]
        ax.scatter(masses, delta_values, alpha=0.6, c=['green' if d > 0 else 'red' for d in delta_values])
        ax.axhline(0, color='black', linestyle='--')
        ax.set_xlabel('Total Mass (Msol)')
        ax.set_ylabel('Delta_chi2 (GR - Hexagonal)')
        ax.set_title('Delta_chi2 vs Total Mass')

        # 3. Results by observing run
        ax = axes[1, 0]
        runs = {'O1': [], 'O2': [], 'O3a': [], 'O3b': []}
        for r in successful:
            if r.name.startswith('GW15'):
                runs['O1'].append(r.favors_hexagonal)
            elif r.name.startswith('GW17'):
                runs['O2'].append(r.favors_hexagonal)
            elif r.name.startswith('GW19') and int(r.name[4:6]) < 11:
                runs['O3a'].append(r.favors_hexagonal)
            else:
                runs['O3b'].append(r.favors_hexagonal)

        x = list(runs.keys())
        hex_counts = [sum(v) for v in runs.values()]
        total_counts = [len(v) for v in runs.values()]

        ax.bar(x, hex_counts, label='Favor Hexagonal', color='green', alpha=0.7)
        ax.bar(x, [t - h for t, h in zip(total_counts, hex_counts)], bottom=hex_counts,
               label='Favor GR', color='red', alpha=0.7)
        ax.set_ylabel('Number of Events')
        ax.set_title('Results by Observing Run')
        ax.legend()

        # 4. Cumulative significance
        ax = axes[1, 1]
        sorted_delta = sorted(delta_values, reverse=True)
        cumulative_z = []
        for i in range(1, len(sorted_delta) + 1):
            subset = sorted_delta[:i]
            if np.std(subset) > 0:
                z = np.mean(subset) / (np.std(subset) / np.sqrt(i))
            else:
                z = 0
            cumulative_z.append(z)
        ax.plot(range(1, len(sorted_delta) + 1), cumulative_z, 'b-', linewidth=2)
        ax.axhline(3, color='orange', linestyle='--', label='3σ')
        ax.axhline(5, color='red', linestyle='--', label='5σ')
        ax.set_xlabel('Number of Events (sorted by Delta_chi2)')
        ax.set_ylabel('Cumulative Z-score')
        ax.set_title('Cumulative Statistical Significance')
        ax.legend()

        plt.tight_layout()
        plt.savefig(FIGURES_DIR / 'gwtc3_full_summary.png', dpi=150)
        plt.savefig(FIGURES_DIR / 'gwtc3_full_summary.pdf')
        plt.close()

        print(f"\nFigures saved to {FIGURES_DIR}")


def main():
    print("=" * 70)
    print("GWTC-3 FULL RAW DATA ANALYSIS")
    print("Testing Hexagonal vs GR Frequency Signatures")
    print("=" * 70)

    analyzer = GWTC3FullAnalyzer()
    analyzer.run_full_analysis()


if __name__ == '__main__':
    main()
