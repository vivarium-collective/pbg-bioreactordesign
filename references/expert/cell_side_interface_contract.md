> **Lifted** from the v2ecoli multiscale-bioprocess investigation
> (`references/expert/cell_side_interface_contract.md`). This workspace
> implements the **reactor side** of this contract (BiRDTransportProcess);
> the studies here test reactor conformance to it. Keep in sync with the
> v2ecoli source as the contract evolves.

# Cell-side interface contract (mbp investigation)

The substitutability story of the `multiscale-bioprocess` investigation —
inherited from the upstream multiscale-bioprocess roadmap (Phase 0 artifact)
— hinges on a **stable cell-side interface**: an engine wired under this
contract is drop-in interchangeable with any other engine that satisfies
it, while the reactor composite stays unchanged.

This document is **load-bearing**: the contract is what makes mbp-03's
"swap Monod placeholder → v2ecoli" a meaningful test rather than a
hand-wave, and what makes future engine substitutions (dFBA, surrogates,
WCM variants) tractable as their own studies rather than re-architectures.

## The contract

A **cell-side engine** is a PBG sub-composite or Process satisfying:

### Inputs (reactor → engine)

| Path | Type | Units | Source |
|---|---|---|---|
| `environment.external_concentrations.<molecule>` | `map[float]` | mM | Reactor liquid-phase concentrations; written by reactor / coupler each step |
| `reactor.volume_L` | `float` | L | Reactor liquid volume (for per-cell concentration units) |
| `reactor.temperature_K` | `float` | K | (optional) Reactor temperature; engines that depend on T re-read |
| `reactor.dissolved_o2` | `float` | mg/L → mM | Dissolved O2 (also surfaced via `environment.external_concentrations.OXYGEN-MOLECULE[p]`) |

### Outputs (engine → reactor)

| Path | Type | Units | Consumed by |
|---|---|---|---|
| `agents.*.metabolism.external_exchange_fluxes.<molecule>` | `array[float]` | mmol/(gDW·h) | Reactor / coupler aggregates across agents × biomass to compute dC/dt |
| `agents.*.listeners.mass.cell_mass` | `float` | fg | Population aggregator sums into `population.total_biomass_gDW` |
| `agents.*.listeners.mass.instantaneous_growth_rate` | `float` | 1/h | Diagnostics; coupler can optionally use to drive BiRD's O2 demand |

### Lifecycle hooks

The engine MUST (when applicable):

- Respect `Division` semantics — `cell_count` may grow over the run.
- Tolerate **time-varying** external concentrations (re-read each step;
  no caching of medium definition past one timestep).
- Run without modification under either a **static** or **driven**
  environment driver (mbp-01's `env_driver_mode` parameter).

### Population scale (added 2026-05-26)

Per chris_feedback_2026_05_26 §4 (Eran's adoption of representative-sampling
over literal-sum 2026-05-26): the cell-side interface treats `cells_per_agent`
as part of the contract.

The engine MUST EITHER:

- **Operate at literal scale** — every simulated agent represents exactly
  one real cell; `cells_per_agent` is implicitly 1.0. The
  PopulationAggregator simply sums per-cell mass. Appropriate when the
  simulated lineage is small AND the reactor density is small enough that
  the literal population is tractable. Default for v2ecoli baseline at
  single-cell inoculation.
- **OR declare an explicit `cells_per_agent` factor** — each simulated
  agent represents `cells_per_agent` real cells under the 0D well-mixed
  assumption (which makes the scaling mathematically equivalent to a
  literal-N simulation for any aggregate observable; see claim
  `well-mixed-representative-sampling-equivalence`). Required when the
  reactor density and the simulation budget would otherwise make literal-N
  intractable (e.g. Beulig 50–80 gDW/L → 10^12–10^13 cells/L). Cost: per-
  cell stochastic variance washes out at the bulk mean as k grows —
  correct biology at high density (population averages dominate), but
  obscures heterogeneity that may matter in regimes with bistable
  expression states.

The PopulationAggregator (mbp-02) consumes `cells_per_agent` and applies
the factor ONLY at the population-aggregate boundary; per-cell observables
are invariant under scaling (regression-guarded by
`per-cell-observables-invariant-under-scaling`). Engines that declare a
non-1.0 scaling factor MUST NOT apply it inside the cell-side state
themselves — the aggregator owns that math.

### Identifiers

External molecule keys use **EcoCyc-style compartment-tagged IDs** (matching
v2ecoli's existing convention; see `v2ecoli/AGENTS.md` for the prefix table):

| Molecule | Key |
|---|---|
| Glucose (periplasmic exchange) | `GLC[p]` |
| O2 (periplasmic) | `OXYGEN-MOLECULE[p]` |
| CO2 (periplasmic) | `CARBON-DIOXIDE[p]` |
| Ammonium | `AMMONIUM[p]` |
| Acetate | `ACET[p]` |

Engines that use other identifier schemes (e.g. BiGG IDs in iML1515)
MUST adapt at the coupler boundary, not inside the engine.

## Conformance

An engine **conforms** if all the input/output store paths above resolve
and types match. The mbp-03 spec PR defines a conformance test fixture
(`tests/test_cell_side_interface.py`) that exercises a minimal mock
engine against the contract.

## Engines tracked under this contract

| Engine | Status | First wired in | Notes |
|---|---|---|---|
| Monod placeholder | planned | mbp-01 | The first study's interface-baseline engine; intentionally trivial — its role is to make the contract concrete. |
| v2ecoli baseline (single cell + Division) | planned | mbp-03 | Already exposes the contract's surface via `media_update` / `exchange_data` / `metabolism.external_exchange_fluxes`; only the reactor-side coupler is new. |
| pbg-bioreactordesign `BiRDTransportProcess` | upstream-PR pending | mbp-03 | Per chris_feedback_2026_05_26 §5 + Eran's Option-B decision 2026-05-26: a transport-only sibling of `BiRDReactorProcess` (no internal biomass ODE by construction, NOT a flag-disabled variant). Tracked as candidate-future-study `pbg-bioreactor-transport-fork` in mbp-06; HARD PREREQ of mbp-03 entering Build. NOT a cell-side engine — it's the REACTOR side under this contract. |
| **pbg-oxidizeme** (OxidizeME ME-model) | **imported** 2026-05-22 | mbp-06 candidate-future | Wrapper at [vivarium-collective/pbg-oxidizeme](https://github.com/vivarium-collective/pbg-oxidizeme) (live demo: [vivarium-collective.github.io/pbg-oxidizeme](https://vivarium-collective.github.io/pbg-oxidizeme/)). Process-bigraph Step bridging Yang et al. 2019 OxidizeME / cobrame / qminospy. **Now installed in this workspace's venv** and registered in `workspace.yaml.imports`; `OxidizeMEStep` is auto-discovered via `allocate_core()`. Real solver gated by upstream stack install (Python 2.7 / proprietary qMINOS — see wrapper README). **THE comparator engine for the v2ecoli-vs-ME-model study in mbp-06.** |
| dFBA (iML1515) | candidate-future | mbp-06 / future | The upstream multiscale-bioprocess roadmap's Phase 3 engine. A comparator study under this contract. |
| Population-coarse-grained surrogate | candidate-future | mbp-06 / future | For high-density (50–80 gDW/L; Beulig 2025) where single-cell-with-Division is intractable. |
| v2ecoli `colony` composite | deferred | TBD | Multi-cell-population variant. Heavier; deferred per the user-confirmed scope. |

## Why this matters

1. **Substitutability is a discipline, not a hope.** Without a named
   contract, "swap the cell side" is implicit and brittle. With one,
   every new engine is reviewable against the same surface.

2. **Engine comparison becomes science.** Once two engines satisfy the
   contract on the same reactor composite under the same conditions,
   any behavioral divergence between them is a fact about the engines,
   not about wiring artifacts.

3. **High-density intractability becomes addressable.** v2ecoli's
   single-cell-with-Division population won't scale to 10^13 cells/L
   (Beulig 50–80 gDW/L). Under the contract, a coarse-grained surrogate
   engine that satisfies the same interface can be wired as a separate
   study; both surrogate and full-WCM produce comparable reactor traces.

## Relationship to upstream

The upstream multiscale-bioprocess roadmap names the substitutability
story in `docs/concepts/phase-based-development.md` ("engine substitution
as a phase boundary" — open question). This document is the v2ecoli-
workspace concrete instantiation of that contract; the upstream roadmap
remains the conceptual reference. When the upstream open question
resolves (e.g. with explicit interface specs), this document syncs.

## Open questions for the first impl PR (mbp-01)

- [ ] Should `agents.*.metabolism.external_exchange_fluxes` be the
      canonical output, or should the contract publish at a higher-level
      `engine.exchange_fluxes` store that the engine writes? (Trade-off:
      transparency vs. interface stability across engines that don't
      have `agents.*` structure — e.g. a coarse-grained surrogate.)
- [ ] How is the `metabolism.external_exchange_fluxes` aggregation across
      agents handled when the engine doesn't expose per-agent fluxes?
      (Likely answer: coupler always writes the aggregate to
      `engine.exchange_fluxes_aggregate`; engines with per-agent data
      populate that themselves.)
- [ ] Where does the conformance test fixture live —
      `v2ecoli/tests/`, the workspace's `tests/behavior/`, or a new
      `tests/interface/`?
