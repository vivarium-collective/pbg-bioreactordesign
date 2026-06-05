#!/usr/bin/env python3
"""Render the workspace dashboard.

After the dashboard runtime was extracted to ``vivarium-dashboard``, this
script just wraps :func:`vivarium_dashboard.lib.report.render_dashboard`.
Kept for backward compatibility with older `bash` callers and CI.
"""
from __future__ import annotations
import sys
from pathlib import Path


def main() -> int:
    try:
        from vivarium_dashboard.lib.report import render_dashboard
    except ImportError:
        print("ERROR: vivarium-dashboard is not installed.", file=sys.stderr)
        print("Install it into the workspace venv:", file=sys.stderr)
        print("    .venv/bin/pip install vivarium-dashboard", file=sys.stderr)
        return 2
    ws_root = Path.cwd()
    if not (ws_root / "workspace.yaml").is_file():
        print(f"ERROR: not a workspace (no workspace.yaml): {ws_root}", file=sys.stderr)
        return 1
    # Make the workspace's own package importable for build_core().
    ws_str = str(ws_root)
    if ws_str not in sys.path:
        sys.path.insert(0, ws_str)
    from vivarium_dashboard.lib._root import set_workspace_root
    set_workspace_root(ws_root)
    out = render_dashboard(ws_root, write_all=True)
    print(f"rendered {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
