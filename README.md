# pbg-bioreactordesign

Process-bigraph wrapper for [BioReactorDesign (BiRD)](https://github.com/NatLabRockies/BioReactorDesign) bioreactor simulation.

Implements a 0D (well-mixed) bioreactor model using BiRD's physics correlations:
- **Higbie penetration theory** for volumetric mass transfer coefficient (kLa)
- **Temperature-dependent Henry's law** for gas-liquid equilibrium
- **Wilke-Chang correlation** for liquid diffusivity
- **Monod kinetics** for microbial growth
- Supports **bubble column**, **stirred tank**, and **airlift** reactor types

## Installation

```bash
git clone <repo-url>
cd pbg-bioreactordesign
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

## Quick Start

```python
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter, gather_emitter_results
from pbg_bioreactordesign import BiRDReactorProcess, make_reactor_document

core = allocate_core()
core.register_link('BiRDReactorProcess', BiRDReactorProcess)
core.register_link('ram-emitter', RAMEmitter)

doc = make_reactor_document(
    reactor_type='bubble_column',
    volume_L=20.0,
    gas_flow_rate_Lpm=2.0,
    temperature_K=310.15,
    initial_biomass_gL=0.5,
    max_growth_rate_per_h=0.7,
    interval=0.5,
)

sim = Composite({'state': doc}, core=core)
sim.run(12.0)  # 12 hours

results = gather_emitter_results(sim)
for entry in results[('emitter',)]:
    print(f"t={entry['time']:.1f}h  DO={entry['dissolved_o2']:.2f} mg/L  X={entry['biomass']:.3f} g/L")
```

## API Reference

### BiRDReactorProcess

Time-driven Process. Time units: **hours**.

| Config Parameter | Type | Default | Description |
|---|---|---|---|
| `reactor_type` | string | `'bubble_column'` | `'bubble_column'`, `'stirred_tank'`, or `'airlift'` |
| `volume_L` | float | 20.0 | Reactor liquid volume (L) |
| `diameter_m` | float | 0.2 | Vessel inner diameter (m) |
| `liquid_height_m` | float | 0.5 | Liquid height (m) |
| `gas_flow_rate_Lpm` | float | 1.0 | Gas flow rate (L/min) |
| `temperature_K` | float | 298.15 | Temperature (K) |
| `pressure_atm` | float | 1.0 | Headspace pressure (atm) |
| `mean_bubble_diameter_mm` | float | 3.0 | Mean bubble diameter (mm) |
| `initial_biomass_gL` | float | 0.5 | Initial biomass (g/L) |
| `initial_do_mgL` | float | 8.0 | Initial dissolved O2 (mg/L) |
| `max_growth_rate_per_h` | float | 0.4 | Max specific growth rate (1/h) |
| `ks_oxygen_mgL` | float | 0.2 | Monod half-saturation for O2 (mg/L) |
| `yield_biomass_o2` | float | 1.2 | Biomass yield on O2 (g_X/g_O2) |
| `maintenance_coeff_per_h` | float | 0.01 | Maintenance coefficient (g_O2/(g_X·h)) |
| `respiratory_quotient` | float | 1.0 | RQ (mol CO2/mol O2) |
| `impeller_power_W` | float | 0.0 | Impeller power for stirred tanks (W) |

**Output ports** (all `overwrite[float]`):

| Port | Units | Description |
|---|---|---|
| `dissolved_o2` | mg/L | Dissolved oxygen concentration |
| `dissolved_co2` | mg/L | Dissolved CO2 concentration |
| `biomass` | g/L | Biomass concentration |
| `gas_holdup` | — | Gas volume fraction (0–1) |
| `kla_o2` | 1/h | O2 mass transfer coefficient |
| `kla_co2` | 1/h | CO2 mass transfer coefficient |
| `o2_saturation` | mg/L | O2 saturation concentration |
| `specific_growth_rate` | 1/h | Current Monod growth rate |
| `o2_uptake_rate` | mg/(L·h) | Oxygen uptake rate |
| `co2_evolution_rate` | mg/(L·h) | CO2 evolution rate |

## Architecture

```
BiRDReactorProcess (Process)
├── inputs:  gas_flow_rate_Lpm ← stores
├── outputs: dissolved_o2, biomass, kla_o2, ... → stores
└── internal: scipy ODE solver (RK45)
         ├── Higbie kLa (penetration theory)
         ├── Henry's law (T-dependent)
         ├── Wilke-Chang diffusivity
         └── Monod growth kinetics
```

## Demo

```bash
source .venv/bin/activate
python demo/demo_report.py
```

Generates `demo/report.html` — a self-contained interactive report with:
- 3D reactor geometry viewers (Three.js)
- Plotly time-series charts (DO, biomass, kLa, uptake rates)
- Bigraph architecture diagrams
- Interactive PBG document trees

Three configurations: E. coli bubble column (20L), CHO stirred tank (250L), yeast airlift (10,000L).

## Tests

```bash
pytest tests/ -v
```
