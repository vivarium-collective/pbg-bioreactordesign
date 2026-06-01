"""BioReactorDesign (BiRD) process-bigraph wrappers.

Three processes share one transport module (`transport.py`):

  - BiRDReactorProcess   — gas-liquid transport + internal Monod biomass ODE.
                           Standalone, single-process reactor. Owns biomass.
  - MonodCellProcess     — the biomass/consumption term factored out as a
                           cell-side engine conforming to the cell-side
                           interface contract (the trivial fixture for the
                           coupled studies; substitutable for dFBA / WCM).
  - BiRDTransportProcess — transport ONLY. Takes biomass / substrate /
                           dissolved-gas concentrations as read-only inputs
                           and emits only the gas-liquid transport
                           contribution as additive deltas. Owns no biomass.

The split lets a cell engine and the reactor compose at shared dissolved-
species stores: under process-bigraph update aggregation, the cell's exchange
deltas and the reactor's transport deltas merge additively (bare `float`
ports, never `overwrite[float]`) without either process owning the other's
state.

Time units: hours. Concentrations: mg/L (gas species) or g/L (biomass).
"""

import numpy as np
from scipy.integrate import solve_ivp
from process_bigraph import Process

# Re-export the transport physics from the shared module so existing imports
# (`from pbg_bioreactordesign.processes import henry_constant`, etc.) keep working.
from pbg_bioreactordesign.transport import (  # noqa: F401
    SPECIES_DATA,
    WC_PSI,
    WC_M,
    RHO_LIQ,
    G,
    water_viscosity,
    henry_constant,
    wilke_chang_diffusivity,
    saturation_concentration,
    higbie_kla,
    bubble_rise_velocity,
    gas_holdup,
    slip_velocity,
    superficial_gas_velocity,
    compute_transport_state,
    o2_transport_rate,
    co2_transport_rate,
)


def monod_kinetics(C_O2, X, cfg):
    """Monod biomass/consumption term (the cell side of the reactor ODE).

    Shared by BiRDReactorProcess (internal biomass) and MonodCellProcess
    (factored-out cell engine) so both compute identical kinetics.

    `cfg` is any mapping with: max_growth_rate_per_h, ks_oxygen_mgL,
    yield_biomass_o2, maintenance_coeff_per_h, respiratory_quotient.

    Returns dict: mu (1/h), q_o2 (mg_O2/(g_X·h)), OUR (mg/(L·h)), CER (mg/(L·h)).
    """
    C_O2_pos = max(C_O2, 0.0)
    mu = cfg['max_growth_rate_per_h'] * C_O2_pos / (
        cfg['ks_oxygen_mgL'] + C_O2_pos)

    # Specific O2 uptake rate (mg_O2 / (g_X · h))
    q_o2 = mu / max(cfg['yield_biomass_o2'], 1e-12) + cfg['maintenance_coeff_per_h']

    OUR = q_o2 * max(X, 0.0) * 1000.0  # mg/(L·h)
    RQ = cfg['respiratory_quotient']
    CER = RQ * OUR * (44.01 / 32.0)    # mg_CO2/(L·h)

    return {'mu': mu, 'q_o2': q_o2, 'OUR': OUR, 'CER': CER}


class BiRDReactorProcess(Process):
    """0D bioreactor: gas-liquid transport + internal Monod biomass ODE.

    Standalone, single-process reactor. Behavior is preserved exactly across
    the transport-module refactor: the transport coefficients come from
    `transport.compute_transport_state` and the kinetics from `monod_kinetics`,
    which are verbatim extractions of the pre-refactor inline computation.

    Time units: hours. All rates in 1/h, concentrations in mg/L or g/L.
    """

    config_schema = {
        # Reactor geometry
        'reactor_type': {'_type': 'string', '_default': 'bubble_column'},
        'volume_L': {'_type': 'float', '_default': 20.0},
        'diameter_m': {'_type': 'float', '_default': 0.2},
        'liquid_height_m': {'_type': 'float', '_default': 0.5},
        # Operating conditions
        'gas_flow_rate_Lpm': {'_type': 'float', '_default': 1.0},
        'temperature_K': {'_type': 'float', '_default': 298.15},
        'pressure_atm': {'_type': 'float', '_default': 1.0},
        'o2_fraction_inlet': {'_type': 'float', '_default': 0.21},
        'co2_fraction_inlet': {'_type': 'float', '_default': 0.0004},
        # Bubble properties
        'mean_bubble_diameter_mm': {'_type': 'float', '_default': 3.0},
        # Microbial kinetics
        'initial_biomass_gL': {'_type': 'float', '_default': 0.5},
        'initial_do_mgL': {'_type': 'float', '_default': 8.0},
        'initial_dco2_mgL': {'_type': 'float', '_default': 0.5},
        'max_growth_rate_per_h': {'_type': 'float', '_default': 0.4},
        'ks_oxygen_mgL': {'_type': 'float', '_default': 0.2},
        'yield_biomass_o2': {'_type': 'float', '_default': 1.2},
        'maintenance_coeff_per_h': {'_type': 'float', '_default': 0.01},
        'respiratory_quotient': {'_type': 'float', '_default': 1.0},
        # Stirred tank specific
        'impeller_power_W': {'_type': 'float', '_default': 0.0},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._state = None  # internal ODE state [C_O2, C_CO2, X]
        self._derived = {}  # derived quantities from last update

    def inputs(self):
        return {
            'gas_flow_rate_Lpm': 'float',
        }

    def outputs(self):
        return {
            'dissolved_o2': 'overwrite[float]',
            'dissolved_co2': 'overwrite[float]',
            'biomass': 'overwrite[float]',
            'gas_holdup': 'overwrite[float]',
            'kla_o2': 'overwrite[float]',
            'kla_co2': 'overwrite[float]',
            'bubble_diameter_mm': 'overwrite[float]',
            'o2_saturation': 'overwrite[float]',
            'co2_saturation': 'overwrite[float]',
            'specific_growth_rate': 'overwrite[float]',
            'o2_uptake_rate': 'overwrite[float]',
            'co2_evolution_rate': 'overwrite[float]',
            'superficial_gas_velocity': 'overwrite[float]',
        }

    def _compute_derived(self, C_O2, C_CO2, X, gas_flow_Lpm):
        """Combine shared transport coefficients + shared Monod kinetics.

        Identical in value to the pre-refactor inline computation: the
        transport half is `compute_transport_state`, the cell half is
        `monod_kinetics`.
        """
        t = compute_transport_state(self.config, gas_flow_Lpm)
        k = monod_kinetics(C_O2, X, self.config)
        return {**t, **k}

    def _ode_rhs(self, t, y, gas_flow_Lpm):
        """RHS of the ODE: y = [C_O2 (mg/L), C_CO2 (mg/L), X (g/L)], dy/dt per hour."""
        C_O2, C_CO2, X = y
        C_O2 = max(C_O2, 0.0)
        C_CO2 = max(C_CO2, 0.0)
        X = max(X, 0.0)

        d = self._compute_derived(C_O2, C_CO2, X, gas_flow_Lpm)

        # dC_O2/dt = kLa·(C* − C_O2) − OUR        (transport − consumption)
        dC_O2 = o2_transport_rate(d['kla_o2'], d['cstar_o2'], C_O2) - d['OUR']
        # dC_CO2/dt = CER + transport stripping
        dC_CO2 = d['CER'] + co2_transport_rate(d['kla_co2'], d['cstar_co2'], C_CO2)
        # dX/dt = mu·X
        dX = d['mu'] * X

        return [dC_O2, dC_CO2, dX]

    def initial_state(self):
        cfg = self.config
        C_O2 = cfg['initial_do_mgL']
        C_CO2 = cfg['initial_dco2_mgL']
        X = cfg['initial_biomass_gL']
        self._state = np.array([C_O2, C_CO2, X])

        d = self._compute_derived(C_O2, C_CO2, X, cfg['gas_flow_rate_Lpm'])
        self._derived = d

        return {
            'dissolved_o2': C_O2,
            'dissolved_co2': C_CO2,
            'biomass': X,
            'gas_holdup': d['alpha_gas'],
            'kla_o2': d['kla_o2'],
            'kla_co2': d['kla_co2'],
            'bubble_diameter_mm': cfg['mean_bubble_diameter_mm'],
            'o2_saturation': d['cstar_o2'],
            'co2_saturation': d['cstar_co2'],
            'specific_growth_rate': d['mu'],
            'o2_uptake_rate': d['OUR'],
            'co2_evolution_rate': d['CER'],
            'superficial_gas_velocity': d['Ug'],
        }

    def update(self, state, interval):
        if self._state is None:
            self.initial_state()

        cfg = self.config
        gas_flow = state.get('gas_flow_rate_Lpm', cfg['gas_flow_rate_Lpm'])

        sol = solve_ivp(
            self._ode_rhs,
            [0.0, interval],
            self._state.tolist(),
            method='RK45',
            args=(gas_flow,),
            rtol=1e-8,
            atol=1e-10,
            max_step=interval / 10.0,
        )

        if sol.success:
            self._state = np.array([
                max(sol.y[0, -1], 0.0),
                max(sol.y[1, -1], 0.0),
                max(sol.y[2, -1], 0.0),
            ])
        else:
            dy = self._ode_rhs(0, self._state.tolist(), gas_flow)
            self._state = np.maximum(self._state + np.array(dy) * interval, 0.0)

        C_O2, C_CO2, X = self._state
        d = self._compute_derived(C_O2, C_CO2, X, gas_flow)
        self._derived = d

        return {
            'dissolved_o2': float(C_O2),
            'dissolved_co2': float(C_CO2),
            'biomass': float(X),
            'gas_holdup': d['alpha_gas'],
            'kla_o2': d['kla_o2'],
            'kla_co2': d['kla_co2'],
            'bubble_diameter_mm': cfg['mean_bubble_diameter_mm'],
            'o2_saturation': d['cstar_o2'],
            'co2_saturation': d['cstar_co2'],
            'specific_growth_rate': d['mu'],
            'o2_uptake_rate': d['OUR'],
            'co2_evolution_rate': d['CER'],
            'superficial_gas_velocity': d['Ug'],
        }


class MonodCellProcess(Process):
    """Monod-growth cell engine — the biomass/consumption term factored out.

    Conforms (in simplified flat form) to the cell-side interface contract:
    it OWNS biomass and exposes the three contract outputs —
    ``cell_mass`` (g/L), ``external_exchange_fluxes`` (per-species specific
    rates), and ``growth_rate`` (1/h). It reads dissolved O2 from the reactor
    and emits its O2/CO2 exchange and its biomass growth as additive deltas so
    it composes with BiRDTransportProcess at the shared stores.

    Trivial by design: it is the conforming fixture that makes the contract
    concrete and the substitution target for dFBA / whole-cell engines.

    The full contract publishes fluxes under ``agents.*.metabolism...``; this
    flat fixture exposes them directly and would adapt at a coupler boundary
    for an agent-structured engine (see cell_side_interface_contract.md).
    """

    config_schema = {
        'initial_biomass_gL': {'_type': 'float', '_default': 0.5},
        'max_growth_rate_per_h': {'_type': 'float', '_default': 0.4},
        'ks_oxygen_mgL': {'_type': 'float', '_default': 0.2},
        'yield_biomass_o2': {'_type': 'float', '_default': 1.2},
        'maintenance_coeff_per_h': {'_type': 'float', '_default': 0.01},
        'respiratory_quotient': {'_type': 'float', '_default': 1.0},
        # When False, biomass is held at initial_biomass_gL (static high-density
        # fixture for the closed-loop liveness study; no growth dynamics needed).
        'growth_enabled': {'_type': 'boolean', '_default': True},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        self._X = None  # internal biomass state (g/L) — the cell owns biomass

    def inputs(self):
        # Reactor → engine: dissolved O2 drives growth + uptake.
        return {
            'dissolved_o2': 'float',
        }

    def outputs(self):
        return {
            # Additive deltas to shared stores (bare float → compose with transport).
            'biomass': 'float',
            'dissolved_o2': 'float',
            'dissolved_co2': 'float',
            # Cell-side contract readouts (diagnostics; overwrite each step).
            'cell_mass': 'overwrite[float]',
            'growth_rate': 'overwrite[float]',
            'external_exchange_fluxes': 'overwrite[map[float]]',
        }

    def _kinetics(self, C_O2):
        return monod_kinetics(C_O2, self._X, self.config)

    def initial_state(self):
        self._X = self.config['initial_biomass_gL']
        k = self._kinetics(self.config['initial_do_mgL'] if 'initial_do_mgL'
                           in self.config else 8.0)
        return {
            'biomass': 0.0,       # delta — store init comes from the composite
            'dissolved_o2': 0.0,
            'dissolved_co2': 0.0,
            'cell_mass': float(self._X),
            'growth_rate': k['mu'],
            'external_exchange_fluxes': {
                'OXYGEN-MOLECULE[p]': -k['q_o2'],
                'CARBON-DIOXIDE[p]': k['q_o2'] * self.config['respiratory_quotient']
                * (44.01 / 32.0),
            },
        }

    def update(self, state, interval):
        if self._X is None:
            self.initial_state()

        C_O2 = state.get('dissolved_o2', 0.0)
        k = self._kinetics(C_O2)

        # Biomass growth (owned by the cell). Held constant when growth disabled.
        dX = k['mu'] * max(self._X, 0.0) * interval if self.config['growth_enabled'] else 0.0
        self._X = max(self._X + dX, 0.0)

        # Exchange applied to shared dissolved stores as deltas. Cap O2 uptake
        # at the O2 actually present so the explicit coupling step can't drive
        # dissolved O2 negative; the exact d[O2]/dt = transport − consumption
        # holds whenever O2 is not limiting (the regime bird-02 validates, and
        # the cap binds less tightly as the coupling interval shrinks — bird-03).
        d_o2 = -min(k['OUR'] * interval, max(C_O2, 0.0))  # consumption (negative)
        d_co2 = k['CER'] * interval                       # evolution (positive)

        return {
            'biomass': dX,
            'dissolved_o2': d_o2,
            'dissolved_co2': d_co2,
            'cell_mass': float(self._X),
            'growth_rate': k['mu'],
            'external_exchange_fluxes': {
                'OXYGEN-MOLECULE[p]': -k['q_o2'],
                'CARBON-DIOXIDE[p]': k['q_o2'] * self.config['respiratory_quotient']
                * (44.01 / 32.0),
            },
        }


class BiRDTransportProcess(Process):
    """Transport-only reactor: emits ONLY the gas-liquid transport contribution.

    Biomass, substrate, and dissolved-gas concentrations are read-only INPUTS;
    the process owns no biomass equation and no substrate consumption. It emits
    the O2/CO2 transport deltas (kLa·(C*−C)) as additive `float` updates to the
    shared dissolved stores, so it composes with a cell engine that supplies the
    consumption side. It NEVER writes the biomass store (biomass is not an
    output), and emits zero glucose transport (no gas-phase glucose).

    This is the REACTOR side of the cell-side interface contract.
    """

    config_schema = {
        'reactor_type': {'_type': 'string', '_default': 'bubble_column'},
        'volume_L': {'_type': 'float', '_default': 20.0},
        'diameter_m': {'_type': 'float', '_default': 0.2},
        'liquid_height_m': {'_type': 'float', '_default': 0.5},
        'gas_flow_rate_Lpm': {'_type': 'float', '_default': 1.0},
        'temperature_K': {'_type': 'float', '_default': 298.15},
        'pressure_atm': {'_type': 'float', '_default': 1.0},
        'o2_fraction_inlet': {'_type': 'float', '_default': 0.21},
        'co2_fraction_inlet': {'_type': 'float', '_default': 0.0004},
        'mean_bubble_diameter_mm': {'_type': 'float', '_default': 3.0},
        'impeller_power_W': {'_type': 'float', '_default': 0.0},
    }

    def inputs(self):
        # Biomass / substrate are read-only inputs (never written back).
        return {
            'dissolved_o2': 'float',
            'dissolved_co2': 'float',
            'biomass': 'float',
            'glucose': 'float',
            'gas_flow_rate_Lpm': 'float',
        }

    def outputs(self):
        return {
            # Additive transport deltas to the shared dissolved stores.
            'dissolved_o2': 'float',
            'dissolved_co2': 'float',
            # Diagnostics (overwrite readouts). NOTE: no `biomass`, no `glucose`
            # output — biomass is read-only and glucose has no gas transport.
            'o2_transport_delta': 'overwrite[float]',
            'co2_transport_delta': 'overwrite[float]',
            'glucose_transport': 'overwrite[float]',
            'biomass_transport_delta': 'overwrite[float]',
            'kla_o2': 'overwrite[float]',
            'kla_co2': 'overwrite[float]',
            'o2_saturation': 'overwrite[float]',
            'co2_saturation': 'overwrite[float]',
            'gas_holdup': 'overwrite[float]',
            'superficial_gas_velocity': 'overwrite[float]',
        }

    def _transport(self, C_O2, C_CO2, gas_flow):
        t = compute_transport_state(self.config, gas_flow)
        d_o2_rate = o2_transport_rate(t['kla_o2'], t['cstar_o2'], C_O2)
        d_co2_rate = co2_transport_rate(t['kla_co2'], t['cstar_co2'], C_CO2)
        return t, d_o2_rate, d_co2_rate

    def initial_state(self):
        t, d_o2_rate, d_co2_rate = self._transport(
            8.0, 0.5, self.config['gas_flow_rate_Lpm'])
        return {
            'dissolved_o2': 0.0,
            'dissolved_co2': 0.0,
            'o2_transport_delta': 0.0,
            'co2_transport_delta': 0.0,
            'glucose_transport': 0.0,
            'biomass_transport_delta': 0.0,
            'kla_o2': t['kla_o2'],
            'kla_co2': t['kla_co2'],
            'o2_saturation': t['cstar_o2'],
            'co2_saturation': t['cstar_co2'],
            'gas_holdup': t['alpha_gas'],
            'superficial_gas_velocity': t['Ug'],
        }

    def update(self, state, interval):
        C_O2 = state.get('dissolved_o2', 0.0)
        C_CO2 = state.get('dissolved_co2', 0.0)
        gas_flow = state.get('gas_flow_rate_Lpm', self.config['gas_flow_rate_Lpm'])

        t, d_o2_rate, d_co2_rate = self._transport(C_O2, C_CO2, gas_flow)
        d_o2 = d_o2_rate * interval
        d_co2 = d_co2_rate * interval

        return {
            'dissolved_o2': d_o2,
            'dissolved_co2': d_co2,
            'o2_transport_delta': d_o2,
            'co2_transport_delta': d_co2,
            'glucose_transport': 0.0,
            'biomass_transport_delta': 0.0,
            'kla_o2': t['kla_o2'],
            'kla_co2': t['kla_co2'],
            'o2_saturation': t['cstar_o2'],
            'co2_saturation': t['cstar_co2'],
            'gas_holdup': t['alpha_gas'],
            'superficial_gas_velocity': t['Ug'],
        }
