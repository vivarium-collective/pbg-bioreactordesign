"""bird-04 acceptance tests: selectable reactor geometry + stirred-tank kLa.

A geometry/correlation configuration field selects the kLa form, and under the
stirred-tank correlation kLa matches the published Van't Riet (1979) result for
the configured power input and superficial gas velocity.
"""

import pytest

from viva_bioreactordesign.transport import (
    compute_transport_state,
    vant_riet_kla,
    superficial_gas_velocity,
)

# A vigorously aerated, mechanically stirred vessel.
ST_CFG = {
    'reactor_type': 'stirred_tank',
    'volume_L': 250.0,
    'diameter_m': 0.6,
    'gas_flow_rate_Lpm': 200.0,
    'temperature_K': 298.15,
    'pressure_atm': 1.0,
    'o2_fraction_inlet': 0.21,
    'co2_fraction_inlet': 0.0004,
    'mean_bubble_diameter_mm': 3.0,
    'impeller_power_W': 500.0,
}


# --- Behavior: geometry-selectable-via-config ---

def test_geometry_selectable_via_config():
    # Geometry (reactor_type) changes the gas holdup correlation.
    bc = compute_transport_state(
        {**ST_CFG, 'reactor_type': 'bubble_column'}, ST_CFG['gas_flow_rate_Lpm'])
    st = compute_transport_state(
        {**ST_CFG, 'reactor_type': 'stirred_tank'}, ST_CFG['gas_flow_rate_Lpm'])
    assert bc['alpha_gas'] != st['alpha_gas']

    # The kLa correlation can be forced independently of geometry: Higbie vs
    # Van't Riet on the same vessel give different kLa.
    higbie = compute_transport_state(
        {**ST_CFG, 'kla_correlation': 'higbie'}, ST_CFG['gas_flow_rate_Lpm'])
    vant_riet = compute_transport_state(
        {**ST_CFG, 'kla_correlation': 'vant_riet'}, ST_CFG['gas_flow_rate_Lpm'])
    assert vant_riet['kla_o2'] != higbie['kla_o2']

    # The default 'auto' resolves the correlation FROM the geometry, so selecting
    # a stirred tank gets the stirred-tank kLa without a second knob (the
    # least-surprise behavior — reactor_type='stirred_tank' was previously
    # silently served Higbie kLa). 'auto' is also the schema default.
    st_default = compute_transport_state(ST_CFG, ST_CFG['gas_flow_rate_Lpm'])
    assert st_default['kla_o2'] == pytest.approx(vant_riet['kla_o2'], rel=1e-12)
    assert st_default['kla_o2'] != pytest.approx(higbie['kla_o2'], rel=1e-12)

    # A non-stirred geometry resolves 'auto' to Higbie (penetration theory).
    bc_default = compute_transport_state(
        {**ST_CFG, 'reactor_type': 'bubble_column'}, ST_CFG['gas_flow_rate_Lpm'])
    bc_higbie = compute_transport_state(
        {**ST_CFG, 'reactor_type': 'bubble_column', 'kla_correlation': 'higbie'},
        ST_CFG['gas_flow_rate_Lpm'])
    assert bc_default['kla_o2'] == pytest.approx(bc_higbie['kla_o2'], rel=1e-12)


# --- Behavior: kla-matches-stirred-tank-correlation ---

def test_kla_matches_stirred_tank_correlation():
    u_g = superficial_gas_velocity(ST_CFG['gas_flow_rate_Lpm'], ST_CFG['diameter_m'])
    p_per_v = ST_CFG['impeller_power_W'] / (ST_CFG['volume_L'] / 1000.0)

    # The helper reproduces the published Van't Riet coalescing-media formula.
    expected = 0.026 * p_per_v ** 0.4 * u_g ** 0.5 * 3600.0
    assert vant_riet_kla(p_per_v, u_g) == pytest.approx(expected, rel=1e-12)
    # Physically plausible for a vigorously aerated stirred tank.
    assert 50.0 < vant_riet_kla(p_per_v, u_g) < 1000.0

    # compute_transport_state under vant_riet uses exactly that correlation value.
    t = compute_transport_state(
        {**ST_CFG, 'kla_correlation': 'vant_riet'}, ST_CFG['gas_flow_rate_Lpm'])
    assert t['kla_o2'] == pytest.approx(vant_riet_kla(p_per_v, u_g), rel=1e-9)


def test_vant_riet_monotonic_in_power_and_aeration():
    # Sanity: kLa rises with both stirring power and gas throughput.
    base = vant_riet_kla(1000.0, 0.01)
    assert vant_riet_kla(2000.0, 0.01) > base   # more power → more transfer
    assert vant_riet_kla(1000.0, 0.02) > base   # more aeration → more transfer
    assert vant_riet_kla(0.0, 0.01) == 0.0      # no stirring → no correlation
