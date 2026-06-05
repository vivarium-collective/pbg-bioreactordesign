"""bird-03 acceptance tests: coupling-interval tunability + convergence.

The coupling interval is a composite parameter (from Phase 1). These tests
confirm it is honored and that the coupled dissolved-O2 trajectory converges
as the interval shrinks — so the interval trades fidelity against cost without
changing the qualitative result.
"""

import pytest
from process_bigraph import Composite, allocate_core
from process_bigraph.emitter import RAMEmitter, gather_emitter_results

from pbg_bioreactordesign.processes import BiRDTransportProcess, MonodCellProcess
from pbg_bioreactordesign.composites import make_coupled_document


@pytest.fixture
def core():
    c = allocate_core()
    c.register_link('BiRDTransportProcess', BiRDTransportProcess)
    c.register_link('MonodCellProcess', MonodCellProcess)
    c.register_link('ram-emitter', RAMEmitter)
    return c


def _do_at(core, interval, t_final, biomass=0.1):
    """Dissolved O2 at t_final for a coupled run at the given coupling interval."""
    doc = make_coupled_document(
        initial_biomass_gL=biomass, growth_enabled=False,
        initial_do_mgL=8.0, interval=interval,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(t_final)
    return sim.state['stores']['dissolved_o2']


def _n_steps(core, interval, t_final=2.0):
    doc = make_coupled_document(
        initial_biomass_gL=0.05, growth_enabled=False, interval=interval)
    sim = Composite({'state': doc}, core=core)
    sim.run(t_final)
    return len(gather_emitter_results(sim)[('emitter',)])


# --- Behavior: coupling-interval-is-configurable ---

def test_coupling_interval_is_configurable(core):
    # The factory threads the interval onto both processes.
    doc = make_coupled_document(interval=0.25)
    assert doc['transport']['interval'] == 0.25
    assert doc['cell']['interval'] == 0.25

    # And runs honor it: a smaller interval takes more steps for the same horizon.
    assert _n_steps(core, 0.5) < _n_steps(core, 0.1)


# --- Behavior: coupled-trajectory-converges-as-interval-decreases ---

def test_coupled_trajectory_converges_as_interval_decreases(core):
    # Convergence vs a fine reference solution. DO at a fixed early-transient
    # time (t=0.4, an exact multiple of every interval) should approach the
    # reference as the coupling interval shrinks.
    #
    # Intervals are kept in the stable regime (kLa·interval < 1); a coarse
    # interval like 0.2 sits in the explicit-Euler overshoot regime for the
    # transport term and is qualitatively wrong — that's the interval-dependence
    # this study characterizes, and the reason to keep the interval small.
    t_final = 0.4
    ref = _do_at(core, 0.00625, t_final)
    intervals = [0.1, 0.05, 0.025, 0.0125]
    errors = [abs(_do_at(core, h, t_final) - ref) for h in intervals]

    # Each halving moves the trajectory strictly closer to the reference.
    assert errors[0] > errors[1] > errors[2] > errors[3]
    # The finest interval is converged to within tolerance.
    assert errors[-1] < 0.005
