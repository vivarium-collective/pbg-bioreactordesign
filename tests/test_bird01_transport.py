"""bird-01 acceptance tests: the transport/biomass responsibility split.

One test per behavior in studies/bird-01-transport-process/study.yaml, plus a
coupled-loop smoke test that previews bird-02. These are workspace pytest tests
(direct process-bigraph), the runnable counterpart to the study's behavior_tests
stubs (which use the dashboard measure grammar).
"""

import math
import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter

from pbg_bioreactordesign.processes import (
    BiRDReactorProcess,
    BiRDTransportProcess,
    MonodCellProcess,
)
from pbg_bioreactordesign.transport import (
    compute_transport_state,
    o2_transport_rate,
    co2_transport_rate,
)
from pbg_bioreactordesign.composites import (
    make_reactor_document,
    make_coupled_document,
)

TRANSPORT_CFG = {
    'reactor_type': 'bubble_column',
    'volume_L': 20.0,
    'diameter_m': 0.2,
    'gas_flow_rate_Lpm': 1.0,
    'temperature_K': 298.15,
    'pressure_atm': 1.0,
    'o2_fraction_inlet': 0.21,
    'co2_fraction_inlet': 0.0004,
    'mean_bubble_diameter_mm': 3.0,
    'impeller_power_W': 0.0,
}


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('BiRDReactorProcess', BiRDReactorProcess)
    c.register_link('BiRDTransportProcess', BiRDTransportProcess)
    c.register_link('MonodCellProcess', MonodCellProcess)
    c.register_link('ram-emitter', RAMEmitter)
    return c


# --- Behavior 1: transport-process-emits-zero-biomass-delta ---

def test_transport_process_emits_zero_biomass_delta(core):
    proc = BiRDTransportProcess(config=TRANSPORT_CFG, core=core)
    # Structural guarantee: biomass is not an output port, so it can never be written.
    assert 'biomass' not in proc.outputs()
    # And biomass IS a read-only input.
    assert 'biomass' in proc.inputs()
    update = proc.update(
        {'dissolved_o2': 6.0, 'dissolved_co2': 0.5, 'biomass': 5.0,
         'glucose': 10.0, 'gas_flow_rate_Lpm': 1.0},
        1.0,
    )
    assert 'biomass' not in update
    assert update['biomass_transport_delta'] == 0.0


# --- Behavior 2: o2-transport-equals-kla-driving-force ---

def test_o2_transport_equals_kla_driving_force(core):
    proc = BiRDTransportProcess(config=TRANSPORT_CFG, core=core)
    C_O2, interval = 6.0, 1.0
    t = compute_transport_state(TRANSPORT_CFG, TRANSPORT_CFG['gas_flow_rate_Lpm'])
    expected = o2_transport_rate(t['kla_o2'], t['cstar_o2'], C_O2) * interval

    update = proc.update(
        {'dissolved_o2': C_O2, 'dissolved_co2': 0.5, 'biomass': 5.0,
         'glucose': 10.0, 'gas_flow_rate_Lpm': 1.0},
        interval,
    )
    assert update['dissolved_o2'] == pytest.approx(expected, rel=1e-9)
    assert update['o2_transport_delta'] == pytest.approx(expected, rel=1e-9)
    # Driving force is positive when below saturation → O2 transferred IN.
    assert expected > 0.0


# --- Behavior 3: glucose-transport-is-zero ---

def test_glucose_transport_is_zero(core):
    proc = BiRDTransportProcess(config=TRANSPORT_CFG, core=core)
    # No gas-phase glucose: glucose is a read-only input, never an output store.
    assert 'glucose' in proc.inputs()
    assert 'glucose' not in proc.outputs()
    update = proc.update(
        {'dissolved_o2': 6.0, 'dissolved_co2': 0.5, 'biomass': 5.0,
         'glucose': 10.0, 'gas_flow_rate_Lpm': 1.0},
        1.0,
    )
    assert update['glucose_transport'] == 0.0


# --- Behavior 4: co2-stripping-symmetric-with-o2 ---

def test_co2_stripping_symmetric_with_o2(core):
    proc = BiRDTransportProcess(config=TRANSPORT_CFG, core=core)
    C_CO2, interval = 5.0, 1.0  # above saturation → should strip OUT (negative)
    t = compute_transport_state(TRANSPORT_CFG, TRANSPORT_CFG['gas_flow_rate_Lpm'])
    expected = co2_transport_rate(t['kla_co2'], t['cstar_co2'], C_CO2) * interval

    update = proc.update(
        {'dissolved_o2': 6.0, 'dissolved_co2': C_CO2, 'biomass': 5.0,
         'glucose': 10.0, 'gas_flow_rate_Lpm': 1.0},
        interval,
    )
    assert update['dissolved_co2'] == pytest.approx(expected, rel=1e-9)
    # Same kLa·(C*−C) form, CO2 Henry constant, stripping sign (C_CO2 > C*).
    assert expected < 0.0


# --- Behavior 5: reactor-process-matches-prerefactor-baseline ---

def test_reactor_process_matches_prerefactor_baseline(core):
    # Pin the transport physics the refactor must not change.
    # cstar_o2 at 298.15 K, 1 atm air (0.21 atm O2) = 8.21 mg/L by construction.
    t = compute_transport_state(BiRDReactorProcess.config_schema and {
        k: v['_default'] for k, v in BiRDReactorProcess.config_schema.items()
    }, 1.0)
    assert t['cstar_o2'] == pytest.approx(8.21, rel=1e-6)
    assert t['kla_o2'] > 0.0

    # Standalone reactor still runs and produces a physically sensible trajectory:
    # dissolved O2 stays within (0, saturation], biomass grows.
    doc = make_reactor_document(interval=1.0)
    sim = Composite({'state': doc}, core=core)
    sim.run(5.0)
    stores = sim.state['stores']
    assert 0.0 < stores['dissolved_o2'] <= 8.21 + 1e-6
    assert stores['biomass'] > 0.5  # grew from the 0.5 g/L inoculum
    assert stores['o2_saturation'] == pytest.approx(8.21, rel=1e-6)


# --- Behavior 6: monod-cell-conforms-to-interface-contract ---

def test_monod_cell_conforms_to_interface_contract(core):
    proc = MonodCellProcess(config={'initial_biomass_gL': 1.0}, core=core)
    out = proc.outputs()
    # The three cell-side contract outputs are present.
    assert 'cell_mass' in out
    assert 'growth_rate' in out
    assert 'external_exchange_fluxes' in out

    proc.initial_state()
    update = proc.update({'dissolved_o2': 8.0}, 1.0)
    assert update['cell_mass'] > 0.0                       # fg/gL biomass
    assert 0.0 <= update['growth_rate'] <= 0.4 + 1e-9      # 1/h, ≤ mu_max
    fluxes = update['external_exchange_fluxes']
    assert isinstance(fluxes, dict)
    assert fluxes['OXYGEN-MOLECULE[p]'] < 0.0              # uptake (negative)
    assert fluxes['CARBON-DIOXIDE[p]'] > 0.0              # evolution (positive)


# --- bird-02 preview: closed-loop liveness (coupled transport + cell) ---

def test_coupled_loop_drops_o2_below_saturation(core):
    # Impactful biomass + small coupling interval → DO settles below saturation.
    doc = make_coupled_document(
        initial_biomass_gL=0.5, growth_enabled=False,
        initial_do_mgL=8.0, interval=0.02,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(3.0)
    do = sim.state['stores']['dissolved_o2']
    assert 0.0 <= do < 8.0          # dropped below saturation, stayed non-negative
    assert sim.state['stores']['biomass'] == pytest.approx(0.5)  # cell owns it; growth off


def test_coupled_loop_recovers_without_biomass(core):
    # Negligible biomass → no consumption → transport holds DO near saturation.
    doc = make_coupled_document(
        initial_biomass_gL=1e-6, growth_enabled=False,
        initial_do_mgL=8.0, interval=0.05,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(3.0)
    do = sim.state['stores']['dissolved_o2']
    assert do == pytest.approx(8.21, abs=0.2)  # pulled to ~saturation
