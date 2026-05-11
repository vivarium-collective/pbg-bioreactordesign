"""Pre-built composite document factories for BiRD bioreactor simulations.

Provides make_reactor_document() which returns a ready-to-run PBG document
with a BiRDReactorProcess wired to stores and a RAM emitter for time-series
data collection.
"""


def make_reactor_document(
    reactor_type='bubble_column',
    volume_L=20.0,
    diameter_m=0.2,
    liquid_height_m=0.5,
    gas_flow_rate_Lpm=1.0,
    temperature_K=298.15,
    pressure_atm=1.0,
    o2_fraction_inlet=0.21,
    co2_fraction_inlet=0.0004,
    mean_bubble_diameter_mm=3.0,
    initial_biomass_gL=0.5,
    initial_do_mgL=8.0,
    initial_dco2_mgL=0.5,
    max_growth_rate_per_h=0.4,
    ks_oxygen_mgL=0.2,
    yield_biomass_o2=1.2,
    maintenance_coeff_per_h=0.01,
    respiratory_quotient=1.0,
    impeller_power_W=0.0,
    interval=1.0,
):
    """Build a composite document for a bioreactor simulation.

    Returns a dict suitable for ``Composite({'state': doc}, core=core)``.

    Args:
        reactor_type: 'bubble_column', 'stirred_tank', or 'airlift'
        volume_L: Reactor liquid volume (L)
        diameter_m: Vessel inner diameter (m)
        liquid_height_m: Liquid height (m)
        gas_flow_rate_Lpm: Volumetric gas flow rate (L/min)
        temperature_K: Operating temperature (K)
        pressure_atm: Headspace pressure (atm)
        interval: Time between process updates (hours)
        ... (see BiRDReactorProcess config_schema for others)

    Returns:
        dict: PBG composite document with reactor, stores, and emitter.
    """
    config = {
        'reactor_type': reactor_type,
        'volume_L': volume_L,
        'diameter_m': diameter_m,
        'liquid_height_m': liquid_height_m,
        'gas_flow_rate_Lpm': gas_flow_rate_Lpm,
        'temperature_K': temperature_K,
        'pressure_atm': pressure_atm,
        'o2_fraction_inlet': o2_fraction_inlet,
        'co2_fraction_inlet': co2_fraction_inlet,
        'mean_bubble_diameter_mm': mean_bubble_diameter_mm,
        'initial_biomass_gL': initial_biomass_gL,
        'initial_do_mgL': initial_do_mgL,
        'initial_dco2_mgL': initial_dco2_mgL,
        'max_growth_rate_per_h': max_growth_rate_per_h,
        'ks_oxygen_mgL': ks_oxygen_mgL,
        'yield_biomass_o2': yield_biomass_o2,
        'maintenance_coeff_per_h': maintenance_coeff_per_h,
        'respiratory_quotient': respiratory_quotient,
        'impeller_power_W': impeller_power_W,
    }

    output_ports = [
        'dissolved_o2', 'dissolved_co2', 'biomass', 'gas_holdup',
        'kla_o2', 'kla_co2', 'bubble_diameter_mm', 'o2_saturation',
        'co2_saturation', 'specific_growth_rate', 'o2_uptake_rate',
        'co2_evolution_rate', 'superficial_gas_velocity',
    ]

    return {
        'reactor': {
            '_type': 'process',
            'address': 'local:BiRDReactorProcess',
            'config': config,
            'interval': interval,
            'inputs': {
                'gas_flow_rate_Lpm': ['stores', 'gas_flow_rate_Lpm'],
            },
            'outputs': {
                port: ['stores', port] for port in output_ports
            },
        },
        'stores': {
            'gas_flow_rate_Lpm': gas_flow_rate_Lpm,
        },
        'emitter': {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'config': {
                'emit': {
                    'dissolved_o2': 'float',
                    'dissolved_co2': 'float',
                    'biomass': 'float',
                    'gas_holdup': 'float',
                    'kla_o2': 'float',
                    'o2_uptake_rate': 'float',
                    'co2_evolution_rate': 'float',
                    'specific_growth_rate': 'float',
                    'time': 'float',
                },
            },
            'inputs': {
                'dissolved_o2': ['stores', 'dissolved_o2'],
                'dissolved_co2': ['stores', 'dissolved_co2'],
                'biomass': ['stores', 'biomass'],
                'gas_holdup': ['stores', 'gas_holdup'],
                'kla_o2': ['stores', 'kla_o2'],
                'o2_uptake_rate': ['stores', 'o2_uptake_rate'],
                'co2_evolution_rate': ['stores', 'co2_evolution_rate'],
                'specific_growth_rate': ['stores', 'specific_growth_rate'],
                'time': ['global_time'],
            },
        },
    }
