"""Unit tests for BiRDReactorProcess."""

import math
import pytest
from process_bigraph import allocate_core

from viva_bioreactordesign.processes import (
    BiRDReactorProcess,
    water_viscosity,
    henry_constant,
    wilke_chang_diffusivity,
    saturation_concentration,
    higbie_kla,
    bubble_rise_velocity,
)


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('BiRDReactorProcess', BiRDReactorProcess)
    return c


# --- Helper function tests ---

def test_water_viscosity_25C():
    mu = water_viscosity(298.15)
    # Water viscosity at 25°C ~ 8.9e-4 Pa·s
    assert 5e-4 < mu < 1.5e-3


def test_water_viscosity_37C():
    mu = water_viscosity(310.15)
    # Water viscosity at 37°C ~ 6.9e-4 Pa·s
    assert 4e-4 < mu < 1.0e-3


def test_henry_constant_o2():
    He = henry_constant('O2', 298.15)
    assert He == pytest.approx(0.032, rel=0.01)


def test_henry_constant_temperature_dependence():
    He_cold = henry_constant('O2', 288.15)
    He_warm = henry_constant('O2', 308.15)
    # Dimensionless Henry constant (c_liq/c_gas) decreases with
    # temperature — gas is less soluble in warmer water
    assert He_cold > He_warm


def test_wilke_chang_o2():
    D = wilke_chang_diffusivity('O2', 298.15)
    # O2 diffusivity in water at 25°C ~ 2.1e-9 m²/s
    assert 1e-10 < D < 1e-8


def test_saturation_concentration_o2():
    cstar = saturation_concentration('O2', 298.15, 0.21)
    # O2 saturation at 25°C, air ~ 8.2 mg/L
    assert 6.0 < cstar < 10.0


def test_higbie_kla_positive():
    kla = higbie_kla('O2', 298.15, u_slip=0.2, d_bubble_m=0.003,
                     alpha_gas=0.05)
    assert kla > 0


def test_higbie_kla_zero_gas():
    kla = higbie_kla('O2', 298.15, u_slip=0.0, d_bubble_m=0.003,
                     alpha_gas=0.0)
    assert kla == 0.0


def test_bubble_rise_velocity():
    u = bubble_rise_velocity(0.003)  # 3mm bubble
    # Terminal velocity of 3mm bubble ~ 0.2-0.3 m/s
    assert 0.1 < u < 0.5


# --- Process tests ---

def test_instantiation(core):
    proc = BiRDReactorProcess(config={}, core=core)
    assert proc.config['reactor_type'] == 'bubble_column'
    assert proc.config['volume_L'] == 20.0


def test_instantiation_custom_config(core):
    proc = BiRDReactorProcess(
        config={'reactor_type': 'stirred_tank', 'volume_L': 250.0},
        core=core)
    assert proc.config['reactor_type'] == 'stirred_tank'
    assert proc.config['volume_L'] == 250.0


def test_initial_state(core):
    proc = BiRDReactorProcess(config={}, core=core)
    state = proc.initial_state()
    assert 'dissolved_o2' in state
    assert 'biomass' in state
    assert 'gas_holdup' in state
    assert 'kla_o2' in state
    assert state['dissolved_o2'] == pytest.approx(8.0)
    assert state['biomass'] == pytest.approx(0.5)


def test_outputs_schema(core):
    proc = BiRDReactorProcess(config={}, core=core)
    outputs = proc.outputs()
    expected = [
        'dissolved_o2', 'dissolved_co2', 'biomass', 'gas_holdup',
        'kla_o2', 'kla_co2', 'bubble_diameter_mm', 'o2_saturation',
        'co2_saturation', 'specific_growth_rate', 'o2_uptake_rate',
        'co2_evolution_rate', 'superficial_gas_velocity',
    ]
    for port in expected:
        assert port in outputs


def test_single_update(core):
    proc = BiRDReactorProcess(config={}, core=core)
    proc.initial_state()
    result = proc.update({}, interval=1.0)
    assert isinstance(result['dissolved_o2'], float)
    assert isinstance(result['biomass'], float)
    assert result['biomass'] > 0


def test_biomass_increases(core):
    """Biomass should increase under aerobic conditions."""
    proc = BiRDReactorProcess(config={
        'gas_flow_rate_Lpm': 2.0,
        'initial_biomass_gL': 0.5,
    }, core=core)
    state0 = proc.initial_state()
    X0 = state0['biomass']
    result = proc.update({}, interval=2.0)
    assert result['biomass'] > X0


def test_o2_responds_to_demand(core):
    """O2 should decrease when biomass consumes it."""
    proc = BiRDReactorProcess(config={
        'gas_flow_rate_Lpm': 0.1,  # low aeration
        'initial_biomass_gL': 5.0,  # high biomass
        'max_growth_rate_per_h': 0.5,
    }, core=core)
    state0 = proc.initial_state()
    result = proc.update({}, interval=0.5)
    # With high biomass and low aeration, DO should drop
    assert result['dissolved_o2'] < state0['dissolved_o2']


def test_multi_step_stability(core):
    """Multiple updates should not produce NaN or negative values."""
    proc = BiRDReactorProcess(config={
        'gas_flow_rate_Lpm': 1.0,
        'initial_biomass_gL': 0.5,
    }, core=core)
    proc.initial_state()
    for _ in range(20):
        result = proc.update({}, interval=0.5)
        assert not math.isnan(result['dissolved_o2'])
        assert not math.isnan(result['biomass'])
        assert result['dissolved_o2'] >= 0
        assert result['biomass'] >= 0


def test_stirred_tank_kla(core):
    """Stirred tank with impeller power should have positive kLa."""
    proc = BiRDReactorProcess(config={
        'reactor_type': 'stirred_tank',
        'impeller_power_W': 50.0,
        'gas_flow_rate_Lpm': 1.0,
        'volume_L': 250.0,
        'diameter_m': 0.6,
    }, core=core)
    state = proc.initial_state()
    assert state['kla_o2'] > 0


def test_airlift_reactor(core):
    """Airlift reactor should produce valid state."""
    proc = BiRDReactorProcess(config={
        'reactor_type': 'airlift',
        'volume_L': 10000.0,
        'diameter_m': 1.8,
        'gas_flow_rate_Lpm': 50.0,
    }, core=core)
    state = proc.initial_state()
    assert state['gas_holdup'] > 0
    assert state['kla_o2'] > 0


def test_config_defaults(core):
    proc = BiRDReactorProcess(config={}, core=core)
    assert proc.config['temperature_K'] == 298.15
    assert proc.config['pressure_atm'] == 1.0
    assert proc.config['mean_bubble_diameter_mm'] == 3.0
    assert proc.config['max_growth_rate_per_h'] == 0.4
    assert proc.config['respiratory_quotient'] == 1.0
