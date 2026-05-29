# BiRD reactor transport refactor for cell-engine coupling

## Question

Can the BiRD 0D reactor physics be refactored so that biomass is an
**input** rather than internal state — letting an external cell engine
(Monod, dFBA, or a whole-cell model) own growth and substrate exchange
while the reactor process owns only gas–liquid transport — and does the
resulting coupled system exercise the reactor→cell feedback that a
single-process reactor cannot?

## Hypothesis

Extracting the shared transport physics (kLa, Higbie penetration, Henry's
law, Wilke–Chang diffusivity) into a module consumed by two processes —
`BiRDReactorProcess` (transport + internal Monod, standalone) and a new
`BiRDTransportProcess` (transport only, biomass-as-input) — preserves the
standalone reactor's behavior exactly while enabling a coupled
configuration in which dissolved-oxygen dynamics respond to externally
supplied biomass. At a biomass density high enough for consumption to
exceed transport supply, dissolved O2 drops below saturation; mass
balance closes across the cell and reactor contributions at the shared
stores.

## Background

The current `BiRDReactorProcess` owns biomass internally via a Monod ODE.
When coupled to a cell engine that *also* owns biomass (e.g. the v2ecoli
whole-cell model), the two ownerships conflict; today this is worked
around with a `bird_disable_internal_biomass` configuration flag that
runs the process with half of itself switched off. This is
interface-dishonest and breaks the substitutability symmetry with the
cell-side interface contract (`cell_side_interface_contract.md`): an
engine wired to the reactor is supposed to *provide* biomass and exchange
fluxes, while the reactor provides only the transport contribution.

The fix is a single-responsibility split. The reactor's transport
contribution for a 0D well-mixed vessel is `kLa·(C* − C)` for O2 (and
symmetric stripping for CO2); there is no gas-phase glucose, so the
glucose transport contribution is zero. Under process-bigraph
update-aggregation, the cell engine's exchange deltas and the reactor's
transport deltas merge at the same dissolved-species stores, so the two
processes compose without either owning the other's state.

This investigation is the upstream reactor-physics dependency for the
v2ecoli `multiscale-bioprocess` investigation: its study `mbp-03`
(cell↔reactor coupling) declares `pbg-bioreactor-transport-fork` as a
hard pipeline prerequisite and cannot enter its Build phase until the
`BiRDTransportProcess` lands here. Coupling-timestep and reactor-geometry
decisions raised by `mbp-03` are resolved within this investigation's
scope.

## Phase 1: Shared transport module + BiRDTransportProcess

**Objective.** Extract the gas–liquid transport physics into a shared
module and build two consumers: keep `BiRDReactorProcess` (transport +
internal Monod biomass ODE) for standalone use, and add
`BiRDTransportProcess` (transport only). The new process takes biomass,
substrate, and dissolved-gas concentrations as **input stores** and emits
**only** the transport contribution — no biomass equation, no substrate
consumption. The biomass/consumption term currently internal to
`BiRDReactorProcess` is factored into a small `MonodCellProcess` that
conforms to the cell-side interface contract; standalone
`BiRDReactorProcess` is then `MonodCellProcess` + transport composed, and
the same `MonodCellProcess` is the trivial cell-side fixture for the
coupled studies (Phases 2–4). Remove the `bird_disable_internal_biomass`
flag, which the refactor makes obsolete.

**Acceptance criteria.**
- `transport-process-emits-zero-biomass-delta` — `BiRDTransportProcess`
  never writes to the biomass store; biomass is read-only input.
- `o2-transport-equals-kla-driving-force` — the emitted dissolved-O2
  delta equals `kLa·(C*_O2 − C_O2)` for the configured operating point.
- `glucose-transport-is-zero` — the process emits no glucose transport
  contribution (no gas-phase glucose).
- `co2-stripping-symmetric-with-o2` — CO2 stripping uses the same
  `kLa·(C* − C)` form with the CO2 Henry constant and correct sign.
- `reactor-process-matches-prerefactor-baseline` — `BiRDReactorProcess`
  run standalone produces trajectories identical (within numerical
  tolerance) to the pre-refactor implementation.
- `monod-cell-conforms-to-interface-contract` — the extracted
  `MonodCellProcess` exposes the cell-side contract's output stores
  (`cell_mass`, `external_exchange_fluxes`, growth rate) with correct
  units, verified against the conformance fixture.

## Phase 2: Closed-loop liveness at impactful biomass density

**Objective.** Compose `BiRDTransportProcess` with the `MonodCellProcess`
from Phase 1 at a **static high biomass density** and confirm the
reactor→cell feedback fires: with consumption exceeding transport supply,
dissolved O2 drops below saturation, and O2 mass balance closes across the
cell-consumption and reactor-transport contributions at the shared store.
This is an integration test of the closed loop, distinct from Phase 1's
per-edge unit checks; it requires a density high enough to perturb reactor
state (a constant scale factor on a representative biomass is sufficient —
no growth dynamics needed at this phase).

**Acceptance criteria.**
- `dissolved-o2-drops-below-saturation-at-high-biomass` — at the high-
  density operating point, steady-state dissolved O2 settles below the
  saturation concentration.
- `o2-mass-balance-closes` — `d[O2]/dt` equals reactor transport minus
  cell consumption at every step, within tolerance.
- `dissolved-o2-recovers-when-biomass-removed` — setting biomass input to
  zero returns dissolved O2 to saturation, confirming the feedback is
  bidirectional and not a fixed offset.

## Phase 3: Coupling-interval tunability

**Objective.** Expose the reactor coupling/update interval as a tunable
parameter (it is currently fixed). Confirm the coupled trajectory is
stable and converges as the interval shrinks, so the interval can be
chosen to balance fidelity against cost.

**Acceptance criteria.**
- `coupling-interval-is-configurable` — the interval is a process
  parameter with a documented default; runs honor the configured value.
- `coupled-trajectory-converges-as-interval-decreases` — halving the
  interval changes the dissolved-O2 trajectory by less than a stated
  tolerance, demonstrating convergence rather than interval-dependence.

## Phase 4: Reactor geometry for the benchmark

**Objective.** Support a stirred-tank kLa correlation (the geometry of
the high-cell-density benchmark) alongside the existing bubble-column
form, selectable via configuration, so the coupled model can be
parameterized to a specific vessel.

**Acceptance criteria.**
- `geometry-selectable-via-config` — reactor geometry (stirred-tank vs
  bubble-column) is chosen by a configuration field.
- `kla-matches-stirred-tank-correlation` — under stirred-tank geometry,
  kLa matches the published correlation for the configured power input
  and superficial gas velocity within tolerance.

## References

- `cell_side_interface_contract.md` — the substitutability contract the
  reactor side must honor (biomass + exchange in, transport out). The
  `MonodCellProcess` built in Phase 1 is the trivial conforming engine.
- BiRD-toolbox physics lineage (kLa / Higbie / Henry / Wilke–Chang) —
  prior art for the transport correlations; see `pbg_bioreactordesign/`.
- v2ecoli `multiscale-bioprocess` investigation, study `mbp-03` — the
  downstream consumer that depends on `BiRDTransportProcess`.
