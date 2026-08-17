"""bird-01 acceptance tests: the transport/biomass responsibility split.

One test per behavior in studies/bird-01-transport-process/study.yaml, plus a
coupled-loop smoke test that previews bird-02. These are workspace pytest tests
(direct process-bigraph), the runnable counterpart to the study's behavior_tests
stubs (which use the dashboard measure grammar).
"""

import math
import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter, gather_emitter_results

from viva_bioreactordesign.processes import (
    BiRDReactorProcess,
    BiRDTransportProcess,
    MonodCellProcess,
)
from viva_bioreactordesign.transport import (
    compute_transport_state,
    o2_transport_rate,
    co2_transport_rate,
    saturation_concentration,
)
from viva_bioreactordesign.composites import (
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

# Golden trajectory captured from the PRE-REFACTOR monolithic BiRDReactorProcess
# at commit b3fe651 (the parent of the transport-split commit dc3714f), running
# make_reactor_document(interval=1.0) at 298.15 K for 12 h. The transport-module
# refactor must reproduce these values to numerical tolerance — that is the whole
# point of the split. Same-machine the diff was exactly 0.0; rel=1e-6 leaves
# headroom for cross-platform RK45 last-bit drift while catching any real
# regression. To regenerate: run the b3fe651 process on the same document and
# re-read these (dissolved_o2, dissolved_co2, biomass) at each time.
#
# Steps start at t=3 deliberately: the t=0 emitter row reports zeros (the emitter
# reads the stores before the process's first write — a pre-existing emission-
# ordering quirk, identical in the old and new code, not a physics value).
PREREFACTOR_TRAJECTORY = {  # t (h): (dissolved_o2, dissolved_co2, biomass) mg/L, mg/L, g/L
    3.0:  (0.08788079558518498, 12.742840246798902, 0.8048602270351891),
    6.0:  (0.05554226888838669, 12.791390444160461, 1.0951615179760936),
    9.0:  (0.04016070617669189, 12.814464044517873, 1.3760488831432685),
    12.0: (0.031190157217811914, 12.827915021801916, 1.647455748884813),
}


def test_reactor_process_matches_prerefactor_baseline(core):
    # Pin the transport physics the refactor must not change.
    # cstar_o2 at 298.15 K, 1 atm air (0.21 atm O2) = 8.21 mg/L by construction.
    t = compute_transport_state(BiRDReactorProcess.config_schema and {
        k: v['_default'] for k, v in BiRDReactorProcess.config_schema.items()
    }, 1.0)
    assert t['cstar_o2'] == pytest.approx(8.21, rel=1e-6)
    assert t['kla_o2'] > 0.0

    # Standalone reactor reproduces the pinned pre-refactor trajectory exactly.
    doc = make_reactor_document(interval=1.0)
    sim = Composite({'state': doc}, core=core)
    sim.run(12.0)
    by_time = {round(r['time']): r for r in gather_emitter_results(sim)[('emitter',)]}

    for t_h, (do, dco2, x) in PREREFACTOR_TRAJECTORY.items():
        row = by_time[round(t_h)]
        assert row['dissolved_o2'] == pytest.approx(do, rel=1e-6), f't={t_h}h DO'
        assert row['dissolved_co2'] == pytest.approx(dco2, rel=1e-6), f't={t_h}h dCO2'
        assert row['biomass'] == pytest.approx(x, rel=1e-6), f't={t_h}h biomass'


# --- Regression: O2 saturation must DECREASE with temperature ---

def test_o2_saturation_decreases_with_temperature():
    """Gas solubility in water falls as T rises (van 't Hoff). c*(O2) at 310 K
    (37 C, the cell-coupled target) must be lower than at 298.15 K (25 C, the
    BiRD test condition), and land near the literature air-saturation value.

    Literature O2 saturation in water under 1 atm air:
        25 C / 298.15 K  ~ 8.2-8.4 mg/L
        37 C / 310.15 K  ~ 6.5-7.0 mg/L  (~20-25% lower)
    """
    p_o2 = 0.21  # atm O2 in air at 1 atm total
    c298 = saturation_concentration('O2', 298.15, p_o2)
    c310 = saturation_concentration('O2', 310.0, p_o2)

    # 298 K value preserved at its calibrated ~8.21 mg/L.
    assert c298 == pytest.approx(8.21, rel=1e-6)
    # Solubility falls with temperature (the bug had it rising).
    assert c310 < c298
    # And lands in the literature 37 C air-saturation band.
    assert 6.5 <= c310 <= 7.5
    # CO2 shares the same Henry/van 't Hoff form: it must fall with T too.
    co2_298 = saturation_concentration('CO2', 298.15, p_o2)
    co2_310 = saturation_concentration('CO2', 310.0, p_o2)
    assert co2_310 < co2_298


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
    assert update['cell_mass'] > 0.0                       # g/L biomass
    assert 0.0 <= update['growth_rate'] <= 0.4 + 1e-9      # 1/h, ≤ mu_max
    fluxes = update['external_exchange_fluxes']
    assert isinstance(fluxes, dict)
    assert fluxes['OXYGEN-MOLECULE[p]'] < 0.0              # uptake (negative)
    assert fluxes['CARBON-DIOXIDE[p]'] > 0.0              # evolution (positive)

    # Contract units: external_exchange_fluxes are molar specific rates
    # mmol/(gDW·h), NOT the mass-specific g_O2/(gDW·h) the Monod q_o2 carries.
    from viva_bioreactordesign.processes import monod_kinetics
    from viva_bioreactordesign.transport import SPECIES_DATA
    k = monod_kinetics(8.0, 1.0, proc.config)
    o2_molar = k['q_o2'] / SPECIES_DATA['O2']['MW'] * 1000.0
    assert fluxes['OXYGEN-MOLECULE[p]'] == pytest.approx(-o2_molar, rel=1e-12)
    # CO2 evolution = RQ · molar O2 uptake (RQ is a molar ratio).
    assert fluxes['CARBON-DIOXIDE[p]'] == pytest.approx(
        proc.config['respiratory_quotient'] * o2_molar, rel=1e-12)
    # And the molar O2 uptake lands in the biologically standard E. coli range.
    assert 1.0 < -fluxes['OXYGEN-MOLECULE[p]'] < 30.0      # mmol/(gDW·h)


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
