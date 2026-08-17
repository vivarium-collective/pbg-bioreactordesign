"""bird-02 acceptance tests: closed-loop liveness at impactful biomass.

The integration test the single-process reactor cannot provide — the
reactor→cell direction only fires when biomass is dense enough to perturb
reactor state. Verifies the three bird-02 acceptance criteria against the
coupled_reactor_cell composite (BiRDTransportProcess + MonodCellProcess).
"""

import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter, gather_emitter_results

from viva_bioreactordesign.processes import (
    BiRDTransportProcess,
    MonodCellProcess,
    monod_kinetics,
)
from viva_bioreactordesign.transport import compute_transport_state, o2_transport_rate
from viva_bioreactordesign.composites import make_coupled_document

TRANSPORT_CFG = {
    'reactor_type': 'bubble_column', 'volume_L': 20.0, 'diameter_m': 0.2,
    'gas_flow_rate_Lpm': 1.0, 'temperature_K': 298.15, 'pressure_atm': 1.0,
    'o2_fraction_inlet': 0.21, 'co2_fraction_inlet': 0.0004,
    'mean_bubble_diameter_mm': 3.0, 'impeller_power_W': 0.0,
}


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('BiRDTransportProcess', BiRDTransportProcess)
    c.register_link('MonodCellProcess', MonodCellProcess)
    c.register_link('ram-emitter', RAMEmitter)
    return c


# --- Behavior 1: dissolved-o2-drops-below-saturation-at-high-biomass ---

def test_do_drops_below_saturation_at_high_biomass(core):
    doc = make_coupled_document(
        initial_biomass_gL=5.0, growth_enabled=False,
        initial_do_mgL=8.0, interval=0.02,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(4.0)
    raw = gather_emitter_results(sim)[('emitter',)]
    do_final = raw[-1]['dissolved_o2']
    do_prev = raw[-2]['dissolved_o2']

    # Consumption outpaces transport → DO settles below saturation, non-negative.
    assert 0.0 <= do_final < 8.21
    assert do_final < 8.0
    # Steady state: last step barely moves.
    assert abs(do_final - do_prev) < 0.1


# --- Behavior 2: o2-mass-balance-closes ---

def test_o2_mass_balance_closes_process_math(core):
    """d[O2] = transport − consumption, computed independently."""
    cell_cfg = {'initial_biomass_gL': 0.1, 'growth_enabled': False}
    transport = BiRDTransportProcess(config=TRANSPORT_CFG, core=core)
    cell = MonodCellProcess(config=cell_cfg, core=core)
    cell.initial_state()

    C_O2, interval = 6.0, 0.01  # below saturation; low biomass → uptake cap inactive
    tu = transport.update(
        {'dissolved_o2': C_O2, 'dissolved_co2': 0.5, 'biomass': 0.1,
         'glucose': 10.0, 'gas_flow_rate_Lpm': 1.0}, interval)
    cu = cell.update({'dissolved_o2': C_O2}, interval)

    net = tu['dissolved_o2'] + cu['dissolved_o2']

    t = compute_transport_state(TRANSPORT_CFG, 1.0)
    k = monod_kinetics(C_O2, 0.1, {
        'max_growth_rate_per_h': 0.4, 'ks_oxygen_mgL': 0.2,
        'yield_biomass_o2': 1.2, 'maintenance_coeff_per_h': 0.01,
        'respiratory_quotient': 1.0})
    expected = (o2_transport_rate(t['kla_o2'], t['cstar_o2'], C_O2)
                - k['OUR']) * interval

    assert net == pytest.approx(expected, rel=1e-9)
    # Reported diagnostics close the balance exactly.
    assert net == pytest.approx(
        tu['o2_transport_delta'] + cu['o2_exchange_delta'], rel=1e-9)


def test_o2_mass_balance_closes_in_composite(core):
    """The shared store actually aggregates the two processes' deltas."""
    doc = make_coupled_document(
        initial_biomass_gL=0.1, growth_enabled=False,
        initial_do_mgL=6.0, interval=0.01,
    )
    sim = Composite({'state': doc}, core=core)
    do0 = sim.state['stores']['dissolved_o2']
    sim.run(0.01)  # one coupling step
    s = sim.state['stores']
    net_change = s['dissolved_o2'] - do0
    assert net_change == pytest.approx(
        s['o2_transport_delta'] + s['o2_exchange_delta'], abs=1e-9)


# --- Behavior 3: dissolved-o2-recovers-when-biomass-removed ---

def test_do_recovers_when_biomass_removed(core):
    # Negligible biomass → no consumption → transport restores DO to saturation.
    doc = make_coupled_document(
        initial_biomass_gL=1e-9, growth_enabled=False,
        initial_do_mgL=2.0, interval=0.05,  # start depleted, watch it recover
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(5.0)
    do = sim.state['stores']['dissolved_o2']
    assert do == pytest.approx(8.21, abs=0.2)  # returned to ~saturation
