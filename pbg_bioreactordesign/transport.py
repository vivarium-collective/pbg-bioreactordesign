"""Shared gas-liquid transport physics for BiRD bioreactor processes.

Extracted from the monolithic BiRDReactorProcess so it can be consumed by
two processes without duplication:

  - BiRDReactorProcess     — transport + internal Monod biomass (standalone)
  - BiRDTransportProcess   — transport only, biomass-as-input (coupled)

All correlations are the same ones BiRD's CFD post-processing uses: Higbie
penetration theory for kLa, temperature-dependent Henry's law, Wilke-Chang
liquid diffusivity. Time units are hours; concentrations mg/L; rates 1/h.

These are pure functions (no process state) so the standalone reactor and the
transport-only process compute byte-identical transport coefficients.
"""

import math


# --- Species property data (from BiRD bird/constants/*.yaml) ---

SPECIES_DATA = {
    'O2': {
        'MW': 31.999,          # g/mol
        'H_298': 0.032,        # dimensionless Henry constant at 298.15 K
        'DH': 1700.0,          # K, temperature correction
        'WC_V': 25.6e-3,       # m³/kmol, Wilke-Chang molar volume
    },
    'CO2': {
        'MW': 44.010,
        'H_298': 0.83,
        'DH': 2400.0,
        'WC_V': 34.0e-3,
    },
    'N2': {
        'MW': 28.013,
        'H_298': 0.015,
        'DH': 1300.0,
        'WC_V': 31.2e-3,
    },
}

# Water properties
WC_PSI = 2.6        # Wilke-Chang association parameter for water
WC_M = 18.0         # kg/kmol, water molecular weight
RHO_LIQ = 1000.0    # kg/m³, liquid density
G = 9.81             # m/s², gravity


def water_viscosity(T):
    """Dynamic viscosity of water (Pa·s) as function of temperature (K).

    Correlation from BiRD bird/constants/h2o.yaml.
    """
    return 2.414e-5 * 10.0 ** (247.8 / (T - 140.0))


def henry_constant(species, T):
    """Temperature-dependent dimensionless Henry constant.

    He(T) = H_298 * exp(DH * (1/T - 1/298.15))
    From BiRD species YAML data.
    """
    sp = SPECIES_DATA[species]
    return sp['H_298'] * math.exp(sp['DH'] * (1.0 / T - 1.0 / 298.15))


def wilke_chang_diffusivity(species, T):
    """Liquid diffusivity (m²/s) via Wilke-Chang correlation.

    D = 1.173e-16 * (psi * M)^0.5 * T / mu / V^0.6
    From BiRD species property calculations.
    """
    mu = water_viscosity(T)
    V_b = SPECIES_DATA[species]['WC_V']
    return 1.173e-16 * (WC_PSI * WC_M) ** 0.5 * T / mu / V_b ** 0.6


def saturation_concentration(species, T, p_partial_atm):
    """Saturation concentration (mg/L) via Henry's law.

    Scales from reference saturation values with a Henry temperature
    correction. For O2 at 25°C / 1 atm air (0.21 atm O2): c* ≈ 8.21 mg/L.
    """
    ref_cstar = {
        'O2': 8.21 / 0.21,    # mg/L per atm O2 at 298.15 K
        'CO2': 1450.0,         # mg/L per atm CO2 at 298.15 K
        'N2': 18.0 / 0.78,    # mg/L per atm N2 at 298.15 K
    }
    sp = SPECIES_DATA[species]
    # Temperature correction: solubility decreases with temperature
    T_factor = math.exp(-sp['DH'] * (1.0 / T - 1.0 / 298.15))
    return ref_cstar.get(species, 10.0) * p_partial_atm * T_factor


def higbie_kla(species, T, u_slip, d_bubble_m, alpha_gas):
    """Volumetric mass transfer coefficient (1/h) via Higbie penetration theory.

    kL = 2 * sqrt(D * u_slip / (pi * d_b)); a = 6 * alpha_gas / d_b; kLa = kL*a*3600

    From BiRD postprocess/post_quantities.py compute_instantaneous_kla().
    """
    D = wilke_chang_diffusivity(species, T)
    if u_slip < 1e-12 or d_bubble_m < 1e-12 or alpha_gas < 1e-12:
        return 0.0
    kL = 2.0 * math.sqrt(D * u_slip / (math.pi * d_bubble_m))
    a = 6.0 * alpha_gas / d_bubble_m
    return kL * a * 3600.0  # convert s⁻¹ to h⁻¹


def bubble_rise_velocity(d_bubble_m, rho_gas=1.2):
    """Terminal rise velocity of a single bubble (m/s).

    Mendelson (1967) wave analogy: u_b = sqrt(2*sigma/(rho_L*d) + g*d/2)
    """
    sigma = 0.07  # N/m, surface tension of water
    d = max(d_bubble_m, 1e-6)
    return math.sqrt(2.0 * sigma / (RHO_LIQ * d) + G * d / 2.0)


def vant_riet_kla(power_per_volume, u_g, coalescing=True):
    """Stirred-tank volumetric mass transfer coefficient (1/h), Van't Riet (1979).

    The standard literature correlation for mechanically-stirred vessels:

        coalescing (water):       kLa [1/s] = 0.026 · (P/V)^0.4 · u_g^0.5
        non-coalescing (electrolyte): kLa [1/s] = 0.002 · (P/V)^0.7 · u_g^0.2

    P/V is gassed power per unit liquid volume (W/m³); u_g is the superficial
    gas velocity (m/s). Returns kLa in 1/h. Valid ~500–10000 W/m³.

    Ref: Van't Riet, K. (1979). Review of measuring methods and results in
    nonviscous gas-liquid mass transfer in stirred vessels. Ind. Eng. Chem.
    Process Des. Dev. 18(3), 357-364.
    """
    if power_per_volume <= 0.0 or u_g <= 0.0:
        return 0.0
    if coalescing:
        kla_per_s = 0.026 * power_per_volume ** 0.4 * u_g ** 0.5
    else:
        kla_per_s = 0.002 * power_per_volume ** 0.7 * u_g ** 0.2
    return kla_per_s * 3600.0


def gas_holdup(reactor_type, Ug, impeller_power_W, volume_L):
    """Gas holdup (-) from superficial gas velocity correlation.

    Selectable by reactor geometry. Extracted verbatim from the monolithic
    reactor's _compute_gas_holdup so the standalone reactor is unchanged.
    """
    if Ug < 1e-12:
        return 0.0
    if reactor_type == 'bubble_column':
        return min(0.6 * Ug ** 0.7, 0.45)
    elif reactor_type == 'stirred_tank':
        PV = impeller_power_W / max(volume_L / 1000.0, 1e-6)
        return min(0.4 * (PV ** 0.2) * (Ug ** 0.5), 0.40)
    else:  # airlift
        return min(0.5 * Ug ** 0.6, 0.35)


def slip_velocity(d_bubble_m, alpha_gas):
    """Slip velocity between gas and liquid (m/s), Richardson-Zaki swarm-corrected."""
    u_single = bubble_rise_velocity(d_bubble_m)
    n_rz = 2.39  # Richardson-Zaki exponent for intermediate Re
    return u_single * (1.0 - alpha_gas) ** n_rz


def superficial_gas_velocity(gas_flow_Lpm, diameter_m):
    """Superficial gas velocity Ug (m/s) for a cylindrical vessel."""
    cross_area = math.pi * (diameter_m / 2.0) ** 2
    Qg = gas_flow_Lpm / 60000.0  # L/min -> m³/s
    return Qg / max(cross_area, 1e-12)


def compute_transport_state(cfg, gas_flow_Lpm):
    """Compute the gas-liquid transport coefficients for the current geometry.

    This is the transport half of the monolithic reactor's _compute_derived —
    everything that does NOT depend on biomass or dissolved-species state.
    Consumed identically by BiRDReactorProcess (standalone) and
    BiRDTransportProcess (coupled), guaranteeing matched physics.

    `cfg` is any mapping with the transport config keys (reactor_type,
    volume_L, diameter_m, temperature_K, pressure_atm, o2_fraction_inlet,
    co2_fraction_inlet, mean_bubble_diameter_mm, impeller_power_W).

    Returns a dict: Ug, alpha_gas, u_slip, kla_o2, kla_co2, cstar_o2,
    cstar_co2, d_b.
    """
    T = cfg['temperature_K']
    P = cfg['pressure_atm']
    d_b = cfg['mean_bubble_diameter_mm'] / 1000.0  # m

    Ug = superficial_gas_velocity(gas_flow_Lpm, cfg['diameter_m'])
    alpha_gas = gas_holdup(
        cfg['reactor_type'], Ug, cfg['impeller_power_W'], cfg['volume_L'])
    u_slip = slip_velocity(d_b, alpha_gas)

    # kLa correlation is selectable (bird-04). Default 'higbie' (penetration
    # theory, all geometries) preserves pre-existing behavior; 'vant_riet'
    # selects the stirred-tank correlation grounded in P/V and u_g.
    correlation = cfg.get('kla_correlation', 'higbie')
    if correlation == 'vant_riet':
        p_per_v = cfg.get('impeller_power_W', 0.0) / max(cfg['volume_L'] / 1000.0, 1e-9)
        kla_o2 = vant_riet_kla(p_per_v, Ug)
        # CO2 kLa scales with diffusivity^0.5 (penetration-theory scaling),
        # consistent with how Higbie relates the two species.
        kla_co2 = kla_o2 * math.sqrt(
            wilke_chang_diffusivity('CO2', T) / wilke_chang_diffusivity('O2', T))
    else:
        kla_o2 = higbie_kla('O2', T, u_slip, d_b, alpha_gas)
        kla_co2 = higbie_kla('CO2', T, u_slip, d_b, alpha_gas)

    cstar_o2 = saturation_concentration('O2', T, P * cfg['o2_fraction_inlet'])
    cstar_co2 = saturation_concentration('CO2', T, P * cfg['co2_fraction_inlet'])

    return {
        'Ug': Ug,
        'alpha_gas': alpha_gas,
        'u_slip': u_slip,
        'kla_o2': kla_o2,
        'kla_co2': kla_co2,
        'cstar_o2': cstar_o2,
        'cstar_co2': cstar_co2,
        'd_b': d_b,
    }


def o2_transport_rate(kla_o2, cstar_o2, C_O2):
    """O2 gas-liquid transport contribution (mg/(L·h)): kLa·(C* − C). Positive = into liquid."""
    return kla_o2 * (cstar_o2 - C_O2)


def co2_transport_rate(kla_co2, cstar_co2, C_CO2):
    """CO2 transport contribution (mg/(L·h)). Symmetric counterpart of O2: −kLa·(C − C*).

    Negative when C_CO2 > C*_CO2 (stripping out of the liquid).
    """
    return -kla_co2 * (C_CO2 - cstar_co2)
