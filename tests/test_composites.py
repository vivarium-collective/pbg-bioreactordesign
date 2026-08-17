"""Integration tests for BiRD bioreactor composites."""

import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter, gather_emitter_results

from viva_bioreactordesign.processes import BiRDReactorProcess
from viva_bioreactordesign.composites import make_reactor_document


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('BiRDReactorProcess', BiRDReactorProcess)
    c.register_link('ram-emitter', RAMEmitter)
    return c


def test_composite_assembly(core):
    doc = make_reactor_document(interval=1.0)
    sim = Composite({'state': doc}, core=core)
    assert sim is not None


def test_composite_short_run(core):
    doc = make_reactor_document(interval=1.0)
    sim = Composite({'state': doc}, core=core)
    sim.run(3.0)
    stores = sim.state['stores']
    assert stores['dissolved_o2'] > 0
    assert stores['biomass'] > 0


def test_emitter_collects_timeseries(core):
    doc = make_reactor_document(interval=1.0)
    sim = Composite({'state': doc}, core=core)
    sim.run(5.0)
    raw = gather_emitter_results(sim)
    data = raw[('emitter',)]
    assert len(data) >= 2
    for entry in data:
        assert 'dissolved_o2' in entry
        assert 'biomass' in entry
        assert 'time' in entry


def test_factory_params_bubble_column(core):
    doc = make_reactor_document(
        reactor_type='bubble_column',
        volume_L=20.0,
        gas_flow_rate_Lpm=2.0,
        interval=0.5,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(2.0)
    stores = sim.state['stores']
    assert stores['gas_holdup'] > 0
    assert stores['kla_o2'] > 0


def test_factory_params_stirred_tank(core):
    doc = make_reactor_document(
        reactor_type='stirred_tank',
        volume_L=250.0,
        diameter_m=0.6,
        impeller_power_W=50.0,
        gas_flow_rate_Lpm=1.0,
        interval=1.0,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(3.0)
    stores = sim.state['stores']
    assert stores['biomass'] > 0


def test_factory_params_airlift(core):
    doc = make_reactor_document(
        reactor_type='airlift',
        volume_L=10000.0,
        diameter_m=1.8,
        gas_flow_rate_Lpm=50.0,
        interval=1.0,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(3.0)
    stores = sim.state['stores']
    assert stores['dissolved_o2'] >= 0
