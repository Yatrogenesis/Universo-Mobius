"""
GeoVMS Data Fetchers
====================

Fetch seismic data from public sources:
- IRIS/NSF SAGE (Earth seismology)
- NASA PDS (Apollo lunar seismometer)
- GSN real-time (SeedLink)

Requires: obspy, requests
"""

import numpy as np
from typing import Optional, List, Tuple, Dict, Any
from dataclasses import dataclass
from datetime import datetime, timedelta
import os
import requests
from pathlib import Path


@dataclass
class SeismicTrace:
    """A single seismic trace (one component, one station)."""
    data: np.ndarray
    sample_rate: float
    station: str
    network: str
    channel: str  # e.g., BHZ (broadband, high-gain, vertical)
    location: str
    start_time: datetime
    end_time: datetime
    latitude: Optional[float] = None
    longitude: Optional[float] = None
    elevation: Optional[float] = None

    @property
    def duration_seconds(self) -> float:
        return len(self.data) / self.sample_rate

    @property
    def npts(self) -> int:
        return len(self.data)


@dataclass
class SeismicEvent:
    """An earthquake or seismic event."""
    event_id: str
    origin_time: datetime
    latitude: float
    longitude: float
    depth_km: float
    magnitude: float
    magnitude_type: str  # ML, Mw, mb, etc.
    region: str
    traces: List[SeismicTrace] = None


class SeismicFetcher:
    """
    Fetch seismic data from IRIS/FDSN web services.

    Uses FDSN standard web services:
    - fdsnws-station: Station metadata
    - fdsnws-dataselect: Waveform data
    - fdsnws-event: Earthquake catalog

    Example:
        fetcher = SeismicFetcher()
        event = fetcher.get_event_data(
            event_id="usp000hvnu",  # 2011 Tohoku earthquake
            stations=["IU.ANMO", "II.PFO"],
            window_minutes=60
        )
    """

    # FDSN web service endpoints
    IRIS_BASE = "https://service.iris.edu"
    STATION_URL = f"{IRIS_BASE}/fdsnws/station/1/query"
    DATASELECT_URL = f"{IRIS_BASE}/fdsnws/dataselect/1/query"
    EVENT_URL = f"{IRIS_BASE}/fdsnws/event/1/query"

    def __init__(self, cache_dir: str = None):
        """
        Initialize fetcher.

        Args:
            cache_dir: Directory to cache downloaded data
        """
        if cache_dir is None:
            cache_dir = os.path.expanduser("~/.geovms_cache")
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def search_events(
        self,
        start_time: datetime,
        end_time: datetime,
        min_magnitude: float = 5.0,
        max_magnitude: float = 10.0,
        min_latitude: float = -90,
        max_latitude: float = 90,
        min_longitude: float = -180,
        max_longitude: float = 180,
        limit: int = 100
    ) -> List[Dict[str, Any]]:
        """
        Search earthquake catalog.

        Returns list of event dictionaries.
        """
        params = {
            "starttime": start_time.isoformat(),
            "endtime": end_time.isoformat(),
            "minmagnitude": min_magnitude,
            "maxmagnitude": max_magnitude,
            "minlatitude": min_latitude,
            "maxlatitude": max_latitude,
            "minlongitude": min_longitude,
            "maxlongitude": max_longitude,
            "limit": limit,
            "format": "text",
            "orderby": "magnitude"
        }

        response = requests.get(self.EVENT_URL, params=params)
        response.raise_for_status()

        events = []
        lines = response.text.strip().split('\n')

        # Skip header line
        for line in lines[1:]:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split('|')
            if len(parts) >= 11:
                events.append({
                    'event_id': parts[0],
                    'time': parts[1],
                    'latitude': float(parts[2]),
                    'longitude': float(parts[3]),
                    'depth_km': float(parts[4]),
                    'magnitude': float(parts[10]),
                    'magnitude_type': parts[9],
                    'region': parts[12] if len(parts) > 12 else ""
                })

        return events

    def get_stations(
        self,
        network: str = "*",
        station: str = "*",
        channel: str = "BH*",
        min_latitude: float = None,
        max_latitude: float = None,
        min_longitude: float = None,
        max_longitude: float = None,
    ) -> List[Dict[str, Any]]:
        """Get available seismic stations."""
        params = {
            "network": network,
            "station": station,
            "channel": channel,
            "level": "station",
            "format": "text"
        }

        if min_latitude is not None:
            params["minlatitude"] = min_latitude
            params["maxlatitude"] = max_latitude
            params["minlongitude"] = min_longitude
            params["maxlongitude"] = max_longitude

        response = requests.get(self.STATION_URL, params=params)
        response.raise_for_status()

        stations = []
        lines = response.text.strip().split('\n')

        for line in lines[1:]:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split('|')
            if len(parts) >= 6:
                stations.append({
                    'network': parts[0],
                    'station': parts[1],
                    'latitude': float(parts[2]),
                    'longitude': float(parts[3]),
                    'elevation': float(parts[4]),
                    'sitename': parts[5]
                })

        return stations

    def fetch_waveforms(
        self,
        network: str,
        station: str,
        location: str,
        channel: str,
        start_time: datetime,
        end_time: datetime
    ) -> Optional[SeismicTrace]:
        """
        Fetch waveform data for a single channel.

        Returns SeismicTrace or None if no data available.
        """
        params = {
            "network": network,
            "station": station,
            "location": location,
            "channel": channel,
            "starttime": start_time.isoformat(),
            "endtime": end_time.isoformat()
        }

        try:
            response = requests.get(self.DATASELECT_URL, params=params)

            if response.status_code == 204:
                # No data available
                return None

            response.raise_for_status()

            # Parse miniSEED data
            # This requires obspy - for now, return raw bytes
            # In production, use: obspy.read(io.BytesIO(response.content))

            # Placeholder: return synthetic data
            # TODO: Integrate ObsPy for real miniSEED parsing
            duration = (end_time - start_time).total_seconds()
            sample_rate = 40.0  # Typical broadband rate
            npts = int(duration * sample_rate)

            return SeismicTrace(
                data=np.zeros(npts),  # Placeholder
                sample_rate=sample_rate,
                station=station,
                network=network,
                channel=channel,
                location=location,
                start_time=start_time,
                end_time=end_time
            )

        except requests.exceptions.RequestException as e:
            print(f"Error fetching {network}.{station}.{location}.{channel}: {e}")
            return None

    def fetch_event_waveforms(
        self,
        event_time: datetime,
        latitude: float,
        longitude: float,
        networks: List[str] = ["IU", "II"],
        channels: List[str] = ["BHZ", "BHN", "BHE"],
        pre_seconds: float = 60,
        post_seconds: float = 600,
        max_stations: int = 20
    ) -> List[SeismicTrace]:
        """
        Fetch waveforms for an event from multiple stations.

        Args:
            event_time: Origin time of the event
            latitude, longitude: Event location
            networks: Networks to search (IU=IRIS/USGS, II=IRIS)
            channels: Channels to request
            pre_seconds: Seconds before event
            post_seconds: Seconds after event
            max_stations: Maximum number of stations

        Returns:
            List of SeismicTrace objects
        """
        start = event_time - timedelta(seconds=pre_seconds)
        end = event_time + timedelta(seconds=post_seconds)

        traces = []

        for network in networks:
            stations = self.get_stations(network=network)[:max_stations]

            for sta_info in stations:
                for channel in channels:
                    trace = self.fetch_waveforms(
                        network=network,
                        station=sta_info['station'],
                        location="*",
                        channel=channel,
                        start_time=start,
                        end_time=end
                    )

                    if trace is not None:
                        trace.latitude = sta_info['latitude']
                        trace.longitude = sta_info['longitude']
                        trace.elevation = sta_info['elevation']
                        traces.append(trace)

        return traces


class LunarFetcher:
    """
    Fetch Apollo lunar seismic data from NASA PDS.

    Data from Apollo 11, 12, 14, 15, 16 Passive Seismic Experiments (PSE).
    Available in SEED format from PDS Geosciences Node.

    The Moon "rings like a bell" - seismic waves propagate for hours
    due to the dry, fractured regolith scattering waves.

    Example:
        fetcher = LunarFetcher()
        events = fetcher.get_event_catalog()
        traces = fetcher.fetch_moonquake("A14_1972_01_04")
    """

    # NASA PDS endpoints
    PDS_BASE = "https://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-apollo_pse"
    CATALOG_URL = "https://pds-geosciences.wustl.edu/lunar/urn-nasa-pds-apollo_seismic_event_catalog"

    # Apollo station locations on the Moon
    STATIONS = {
        "S11": {"lat": 0.6732, "lon": 23.4732, "mission": "Apollo 11"},
        "S12": {"lat": -3.0128, "lon": -23.4219, "mission": "Apollo 12"},
        "S14": {"lat": -3.6453, "lon": -17.4714, "mission": "Apollo 14"},
        "S15": {"lat": 26.1322, "lon": 3.6339, "mission": "Apollo 15"},
        "S16": {"lat": -8.9734, "lon": 15.5011, "mission": "Apollo 16"},
    }

    # Event types
    EVENT_TYPES = {
        "A": "Artificial impact (rocket stage, etc.)",
        "M": "Meteorite impact",
        "SH": "Shallow moonquake",
        "DM": "Deep moonquake",
        "TH": "Thermal moonquake (day/night boundary)"
    }

    def __init__(self, cache_dir: str = None):
        if cache_dir is None:
            cache_dir = os.path.expanduser("~/.geovms_cache/lunar")
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_event_catalog(self) -> List[Dict[str, Any]]:
        """
        Get catalog of lunar seismic events.

        Returns list of events with:
        - Event ID
        - Event type (A=artificial, M=meteorite, SH=shallow, DM=deep, TH=thermal)
        - Date/time
        - Quality grade
        - Amplitude
        - Location (if determined)
        """
        # Catalog file path at PDS
        catalog_file = self.cache_dir / "event_catalog.csv"

        # For demo, return synthetic catalog based on published data
        # In production, download from PDS
        events = [
            {
                "event_id": "A12_1969_11_20",
                "type": "A",
                "description": "Apollo 12 LM ascent stage impact",
                "date": "1969-11-20T22:17:17",
                "station": "S12",
                "amplitude": "Large",
                "quality": "A",
                "lat": -3.94,
                "lon": -21.20
            },
            {
                "event_id": "DM_1972_01_04_001",
                "type": "DM",
                "description": "Deep moonquake cluster A1",
                "date": "1972-01-04T08:23:45",
                "station": "S14",
                "amplitude": "Medium",
                "quality": "B",
                "depth_km": 900
            },
            {
                "event_id": "M_1972_05_13_001",
                "type": "M",
                "description": "Large meteorite impact",
                "date": "1972-05-13T14:55:00",
                "station": "S12,S14,S15,S16",
                "amplitude": "Very Large",
                "quality": "A"
            },
            {
                "event_id": "SH_1973_03_13_001",
                "type": "SH",
                "description": "Shallow moonquake",
                "date": "1973-03-13T07:12:00",
                "station": "S15",
                "amplitude": "Small",
                "quality": "C"
            },
        ]

        return events

    def fetch_moonquake(
        self,
        event_id: str,
        channels: List[str] = ["MH1", "MH2", "MHZ", "SHZ"]
    ) -> List[SeismicTrace]:
        """
        Fetch waveforms for a moonquake event.

        Args:
            event_id: Event identifier from catalog
            channels: Seismometer channels
                - MH1, MH2: Mid-period horizontal
                - MHZ: Mid-period vertical
                - SHZ: Short-period vertical

        Returns:
            List of SeismicTrace objects

        Note: Lunar seismograms can extend for HOURS due to
        scattering in the fractured regolith.
        """
        # In production, download from PDS
        # For now, return synthetic characteristic lunar seismogram

        # Lunar seismograms have distinctive characteristics:
        # 1. Very slow rise (emergent P-wave)
        # 2. No clear S-wave arrival
        # 3. Long coda (hours of ringing)
        # 4. Spindle-shaped envelope

        sample_rate = 6.625  # Apollo PSE sample rate (mid-period)
        duration_hours = 2.0
        npts = int(duration_hours * 3600 * sample_rate)

        t = np.linspace(0, duration_hours * 3600, npts)

        # Characteristic lunar "spindle" waveform
        # Rise time ~ 10 minutes, decay time ~ 1 hour
        rise_time = 600  # 10 minutes
        decay_time = 3600  # 1 hour

        envelope = np.exp(-t / decay_time) - np.exp(-t / rise_time)
        envelope = envelope / np.max(envelope)

        # Add scattered arrivals (simulating regolith scattering)
        noise = np.random.randn(npts) * 0.1
        filtered_noise = np.convolve(noise, np.ones(100)/100, mode='same')

        signal = envelope * (np.sin(2 * np.pi * 0.5 * t) + filtered_noise)

        traces = []
        for channel in channels:
            trace = SeismicTrace(
                data=signal.astype(np.float32),
                sample_rate=sample_rate,
                station="S12",  # Example station
                network="XA",   # Apollo network code
                channel=channel,
                location="",
                start_time=datetime(1972, 1, 4, 8, 23, 45),
                end_time=datetime(1972, 1, 4, 10, 23, 45),
                latitude=self.STATIONS["S12"]["lat"],
                longitude=self.STATIONS["S12"]["lon"]
            )
            traces.append(trace)

        return traces

    def describe_lunar_interior(self) -> Dict[str, Any]:
        """
        Return known lunar interior structure from Apollo seismology.

        Based on analysis of artificial impacts and moonquakes.
        """
        return {
            "crust": {
                "thickness_km": 45,  # Average (variable: 30-60 km)
                "description": "Anorthositic upper crust, fractured"
            },
            "upper_mantle": {
                "depth_range_km": (45, 500),
                "composition": "Olivine-pyroxene",
                "state": "Solid"
            },
            "lower_mantle": {
                "depth_range_km": (500, 1100),
                "composition": "Possibly more Fe-rich",
                "state": "Partially molten (attenuating zone)"
            },
            "core": {
                "radius_km": 350,
                "description": "Small, possibly liquid outer + solid inner",
                "composition": "Iron-rich"
            },
            "notes": [
                "No global magnetic field (unlike Earth)",
                "Deep moonquakes at 700-1100 km depth",
                "Shallow moonquakes near surface (thermal?)",
                "Extreme scattering in upper 20 km (regolith)"
            ]
        }


def demo_fetch():
    """Demonstrate data fetching capabilities."""
    print("=" * 60)
    print("GeoVMS Seismic Data Fetcher Demo")
    print("=" * 60)

    # Earth seismic
    print("\n[1] Earth Seismic Data (IRIS)")
    fetcher = SeismicFetcher()

    print("Searching for recent M6+ earthquakes...")
    try:
        events = fetcher.search_events(
            start_time=datetime(2024, 1, 1),
            end_time=datetime(2024, 12, 31),
            min_magnitude=6.5,
            limit=5
        )
        print(f"Found {len(events)} events")
        for e in events[:3]:
            print(f"  - M{e['magnitude']:.1f} {e['region'][:40]}")
    except Exception as ex:
        print(f"  (Network error: {ex})")

    # Lunar seismic
    print("\n[2] Lunar Seismic Data (Apollo)")
    lunar = LunarFetcher()

    events = lunar.get_event_catalog()
    print(f"Apollo seismic catalog: {len(events)} events")
    for e in events:
        print(f"  - {e['event_id']}: {e['type']} - {e['description'][:40]}")

    print("\nLunar interior structure:")
    interior = lunar.describe_lunar_interior()
    print(f"  Crust: {interior['crust']['thickness_km']} km")
    print(f"  Core radius: {interior['core']['radius_km']} km")

    print("\nFetching synthetic moonquake waveform...")
    traces = lunar.fetch_moonquake("DM_1972_01_04_001")
    print(f"Got {len(traces)} traces")
    if traces:
        t = traces[0]
        print(f"  Duration: {t.duration_seconds/3600:.1f} hours")
        print(f"  Sample rate: {t.sample_rate} Hz")
        print(f"  Samples: {t.npts}")


if __name__ == "__main__":
    demo_fetch()
