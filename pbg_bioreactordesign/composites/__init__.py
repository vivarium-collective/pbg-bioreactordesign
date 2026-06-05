"""BiRD bioreactor composite documents + composite-spec discovery.

Two flavors of composite construction live in this package:

1. **Hand-coded factory** — `make_reactor_document(...)` builds a PBG
   state-dict programmatically for callers that want full control over
   the reactor wiring. Re-exported from the legacy single-file API.

2. **Declarative `*.composite.yaml`** — sibling files in this directory
   follow the pbg-superpowers composite-spec convention.
   `build_composite()` loads one by name and instantiates
   `process_bigraph.Composite` with parameter substitution. The
   dashboard's composite explorer discovers these automatically once
   the package is installed in a workspace.

Both flavors are equivalent — pick the one that fits your use case.
"""
from __future__ import annotations
import re
from pathlib import Path
from typing import Any

import yaml
from process_bigraph import allocate_core
from process_bigraph.emitter import RAMEmitter

from pbg_bioreactordesign.processes import (
    BiRDReactorProcess,
    BiRDTransportProcess,
    MonodCellProcess,
)

# Re-export the hand-coded factories (standalone reactor + coupled reactor↔cell).
from pbg_bioreactordesign.composites._factory import (
    make_reactor_document,
    make_coupled_document,
)


def register_bioreactordesign(core=None):
    """Return a core with the BiRD processes, the RAM emitter, and the
    bioreactor Visualization registered.

    Registers all three processes so both the standalone reactor specs and the
    coupled reactor↔cell spec (BiRDTransportProcess + MonodCellProcess) build
    and are discoverable by the dashboard catalog."""
    if core is None:
        core = allocate_core()
    core.register_link('BiRDReactorProcess', BiRDReactorProcess)
    core.register_link('BiRDTransportProcess', BiRDTransportProcess)
    core.register_link('MonodCellProcess', MonodCellProcess)
    core.register_link('ram-emitter', RAMEmitter)
    # Register Visualization Steps so composites can wire them by name.
    from pbg_bioreactordesign.visualizations import BioreactorPlots
    core.register_link('BioreactorPlots', BioreactorPlots)
    return core


# ---------------------------------------------------------------------------
# Declarative composite-spec loader (*.composite.yaml)
# ---------------------------------------------------------------------------

_COMPOSITES_DIR = Path(__file__).parent

_FULL_PLACEHOLDER = re.compile(r"^\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}$")
_INLINE_PLACEHOLDER = re.compile(r"\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


def _cast(value: Any, declared_type: str | None) -> Any:
    if declared_type is None:
        return value
    if declared_type == "float":
        return float(value)
    if declared_type == "int":
        return int(value)
    if declared_type in ("string", "str"):
        return str(value)
    if declared_type == "bool":
        if isinstance(value, str):
            return value.strip().lower() in ("true", "1", "yes")
        return bool(value)
    return value


def _substitute(state: Any, params: dict, overrides: dict) -> Any:
    if isinstance(state, dict):
        return {k: _substitute(v, params, overrides) for k, v in state.items()}
    if isinstance(state, list):
        return [_substitute(v, params, overrides) for v in state]
    if isinstance(state, str):
        m = _FULL_PLACEHOLDER.match(state)
        if m:
            pname = m.group(1)
            pdef = params.get(pname, {})
            raw = overrides.get(pname, pdef.get("default"))
            return _cast(raw, pdef.get("type"))
        if _INLINE_PLACEHOLDER.search(state):
            return _INLINE_PLACEHOLDER.sub(
                lambda mm: str(overrides.get(mm.group(1), params.get(mm.group(1), {}).get("default", ""))),
                state,
            )
    return state


def list_composite_specs() -> list[str]:
    """Return short names of every `*.composite.yaml` shipped in this package."""
    out: list[str] = []
    for path in sorted(_COMPOSITES_DIR.glob("*.composite.yaml")):
        out.append(path.name[: -len(".composite.yaml")])
    return out


def load_composite_spec(name: str) -> dict:
    """Load and parse a named composite spec. `name` is the stem (no suffix)."""
    path = _COMPOSITES_DIR / f"{name}.composite.yaml"
    if not path.is_file():
        raise FileNotFoundError(f"composite spec not found: {path}")
    return yaml.safe_load(path.read_text())


def build_composite(name: str, *, overrides: dict | None = None, core=None):
    """Load a *.composite.yaml by name and instantiate process_bigraph.Composite.

    overrides: parameter overrides (keys must match spec.parameters)
    core:      optional pre-built core; otherwise register_bioreactordesign() is used
    """
    from process_bigraph import Composite

    spec = load_composite_spec(name)
    if not isinstance(spec, dict) or "state" not in spec or "name" not in spec:
        raise ValueError(f"composite '{name}' missing required keys (name, state)")

    if core is None:
        core = register_bioreactordesign()

    params = spec.get("parameters") or {}
    state = _substitute(spec.get("state") or {}, params, overrides or {})
    return Composite({"state": state}, core=core)


__all__ = [
    'make_reactor_document',
    'make_coupled_document',
    'register_bioreactordesign',
    'list_composite_specs',
    'load_composite_spec',
    'build_composite',
]
