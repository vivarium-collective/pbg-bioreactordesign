#!/usr/bin/env python3
"""Validate workspace.yaml + cross-references. Exit non-zero on failure."""
from __future__ import annotations
import hashlib
import importlib
import json
import re
import subprocess
import sys
from pathlib import Path

import yaml
from jsonschema import Draft7Validator, FormatChecker, ValidationError


def _try_get_registry(ws_root: Path, ws_data: dict) -> set | None:
    """Try to introspect build_core() and return a set of registered process names.

    Returns None if introspection fails (caller should warn, not fail).
    """
    package_name = ws_data.get("package_path") or ("pbg_" + ws_data.get("name", "").replace("-", "_"))
    try:
        py = sys.executable
        script = f"""
import json, sys
try:
    from {package_name}.core import build_core
    core = build_core()
    try:
        names = sorted(core.process_registry.list())
    except Exception:
        try:
            names = sorted(core.process_registry.registry.keys())
        except Exception:
            names = []
    print(json.dumps(names))
except Exception as e:
    print(json.dumps([]))
"""
        result = subprocess.run(
            [py, "-c", script],
            cwd=ws_root, capture_output=True, text=True, timeout=15,
        )
        if result.returncode != 0 or not result.stdout.strip():
            return None
        names = json.loads(result.stdout.strip().split("\n")[-1])
        if isinstance(names, list):
            return set(names)
        return None
    except Exception:
        return None


def _validate_simulation_processes(ws_data: dict, registry: set | None) -> list[str]:
    """Return error strings for simulation process refs not in registry.

    If registry is None (introspection failed), returns empty list (warning only).
    """
    if registry is None:
        return []
    errors = []
    for sim in ws_data.get("simulations", []) or []:
        if not isinstance(sim, dict):
            continue
        for proc in sim.get("processes", []) or []:
            if proc not in registry:
                errors.append(
                    f"simulation '{sim.get('name', '?')}' references missing process '{proc}'"
                )
    return errors


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


WS_ROOT = Path(__file__).resolve().parents[1]
WS_FILE = WS_ROOT / "workspace.yaml"


def _schema(name: str) -> dict:
    p = WS_ROOT / ".pbg" / "schemas" / name
    if not p.exists():
        sys.exit(f"missing schema at {p}; was workspace scaffolded?")
    return json.loads(p.read_text())


def _fail(msg: str) -> None:
    print(f"LINT FAIL: {msg}", file=sys.stderr)
    sys.exit(1)


def main() -> None:
    if not WS_FILE.exists():
        _fail(f"{WS_FILE} not found — run inside a workspace")
    ws = yaml.safe_load(WS_FILE.read_text())

    # Schema version guard: fail clearly for v1 workspaces.
    schema_ver = ws.get("schema_version")
    if schema_ver != 2:
        _fail(
            f"workspace.yaml is schema v{schema_ver}; "
            "run `python3 scripts/_migrate_v1_to_v2.py` to migrate to v2 before linting."
        )

    Draft7Validator(_schema("workspace.schema.json"), format_checker=FormatChecker()).validate(ws)

    # Datasets
    for d in ws.get("datasets", []):
        path, url, sha = d.get("path"), d.get("url"), d.get("sha256")
        if path is not None:
            full = WS_ROOT / path
            if not full.exists():
                _fail(f"dataset '{d['name']}' path missing: {path}")
            elif sha:
                actual = _sha256(full)
                if actual != sha:
                    _fail(f"dataset '{d['name']}' sha256 mismatch (recorded={sha[:16]}…, actual={actual[:16]}…)")
        elif url is not None:
            if not sha:
                _fail(f"dataset '{d['name']}' has url but no sha256")
        else:
            _fail(f"dataset '{d['name']}' has neither path nor url")

    # References cross-ref
    refs_yaml = WS_ROOT / "references" / "claims.yaml"
    bib = WS_ROOT / "references" / "papers.bib"
    if refs_yaml.exists() and bib.exists():
        claims = (yaml.safe_load(refs_yaml.read_text()) or {}).get("claims", {}) or {}
        bib_text = bib.read_text()
        bib_keys = set(re.findall(r"@\w+\{([A-Za-z0-9_:-]+),", bib_text))
        for claim, key in claims.items():
            if isinstance(key, list):
                for k in key:
                    if k not in bib_keys:
                        _fail(f"claim '{claim}' references missing bib key '{k}'")
            elif isinstance(key, str):
                if key not in bib_keys:
                    _fail(f"claim '{claim}' references missing bib key '{key}'")

    # Expert docs path existence + optional sha256 verification
    for doc in ws.get("expert_docs", []) or []:
        doc_path = doc.get("path", "")
        if doc_path:
            full = WS_ROOT / doc_path
            if not full.exists():
                _fail(
                    f"expert_docs entry '{doc.get('name', '?')}' path missing: {doc_path} "
                    f"(expected at {full})"
                )
            sha = doc.get("sha256")
            if sha:
                actual = _sha256(full)
                if actual != sha:
                    _fail(f"expert_docs '{doc.get('name', '?')}' sha256 mismatch (recorded={sha[:16]}…, actual={actual[:16]}…)")

    # References PDFs: always verify sha256 + bib_key cross-reference
    bib = WS_ROOT / "references" / "papers.bib"
    bib_keys: set = set()
    if bib.exists():
        bib_text = bib.read_text()
        bib_keys = set(re.findall(r"@\w+\{([A-Za-z0-9_:-]+),", bib_text))
    for pdf_entry in ws.get("references_pdfs", []) or []:
        bib_key = pdf_entry.get("bib_key", "")
        pdf_path = pdf_entry.get("path", "")
        sha = pdf_entry.get("sha256", "")
        if bib_key and bib_keys and bib_key not in bib_keys:
            _fail(f"references_pdfs entry '{bib_key}' has no matching entry in papers.bib")
        if pdf_path:
            full = WS_ROOT / pdf_path
            if not full.exists():
                _fail(f"references_pdfs '{bib_key}' path missing: {pdf_path}")
            elif sha:
                actual = _sha256(full)
                if actual != sha:
                    _fail(f"references_pdfs '{bib_key}' sha256 mismatch (recorded={sha[:16]}…, actual={actual[:16]}…)")

    # Collect registered simulation names.
    registered_sim_names: set = set()
    for sim in ws.get("simulations", []) or []:
        if not isinstance(sim, dict):
            continue
        sim_name = sim.get("name", "?")
        registered_sim_names.add(sim_name)
        # t_end > t_start
        t_start = sim.get("t_start")
        t_end = sim.get("t_end")
        if t_start is not None and t_end is not None:
            if float(t_end) <= float(t_start):
                _fail(f"simulation '{sim_name}' has t_end ({t_end}) <= t_start ({t_start})")

    # Visualization simulation references.
    for viz in ws.get("visualizations", []) or []:
        if not isinstance(viz, dict):
            continue
        viz_name = viz.get("name", "?")
        # simulation reference
        viz_sim = viz.get("simulation")
        if viz_sim and registered_sim_names and viz_sim not in registered_sim_names:
            _fail(f"visualization '{viz_name}' references missing simulation '{viz_sim}'")

    # Simulation process references (best-effort: warn if registry unavailable).
    has_process_refs = any(
        isinstance(s, dict) and s.get("processes")
        for s in ws.get("simulations", []) or []
    )
    if has_process_refs:
        registry = _try_get_registry(WS_ROOT, ws)
        if registry is None:
            print(
                "LINT WARNING: Could not introspect build_core() — "
                "simulation process references not validated.",
                file=sys.stderr,
            )
        else:
            for err in _validate_simulation_processes(ws, registry):
                _fail(err)

    print("workspace lint: OK")


if __name__ == "__main__":
    main()
