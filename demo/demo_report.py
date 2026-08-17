"""Demo: BiRD multi-configuration bioreactor report.

Runs three E. coli fermentation scenarios at different scales and
aeration strategies, generates interactive 3D reactor geometry viewers
with Three.js, Plotly time-series charts, colored bigraph-viz
architecture diagrams, and navigatable PBG document trees — all in a
single self-contained HTML report.
"""

import json
import math
import os
import base64
import subprocess
import tempfile
import time as _time

import numpy as np
from process_bigraph import allocate_core

from viva_bioreactordesign.processes import BiRDReactorProcess
from viva_bioreactordesign.composites import make_reactor_document


# ── Simulation Configs ──────────────────────────────────────────────

CONFIGS = [
    {
        'id': 'ecoli_lab',
        'title': 'E. coli Lab-Scale Batch',
        'subtitle': '2 L bubble column — exponential growth, O\u2082 crash',
        'description': (
            'A benchtop bubble column (2 L) growing E. coli K-12 at 37\u00b0C in '
            'minimal medium. With moderate aeration (0.5 vvm), the culture '
            'initially enjoys excess dissolved oxygen. As biomass doubles every '
            '~60 min (μ_max = 0.7 h\u207b\u00b9), the oxygen uptake rate overtakes the '
            'mass transfer capacity, causing a characteristic DO crash. This is '
            'the classic aerobic E. coli batch profile seen in teaching labs.'
        ),
        'config': {
            'reactor_type': 'bubble_column',
            'volume_L': 2.0,
            'diameter_m': 0.1,
            'liquid_height_m': 0.25,
            'gas_flow_rate_Lpm': 1.0,
            'temperature_K': 310.15,
            'pressure_atm': 1.0,
            'mean_bubble_diameter_mm': 2.5,
            'initial_biomass_gL': 0.05,
            'initial_do_mgL': 7.0,
            'initial_dco2_mgL': 0.4,
            'max_growth_rate_per_h': 0.7,
            'ks_oxygen_mgL': 0.1,
            'yield_biomass_o2': 0.8,
            'maintenance_coeff_per_h': 0.02,
            'respiratory_quotient': 1.0,
        },
        'n_snapshots': 30,
        'total_time': 10.0,
        'color_scheme': 'indigo',
        'reactor_color': '#6366f1',
    },
    {
        'id': 'ecoli_str',
        'title': 'E. coli High-Cell-Density (Stirred Tank)',
        'subtitle': '50 L stirred tank — high OD, impeller-enhanced kLa',
        'description': (
            'A 50 L stirred-tank bioreactor growing E. coli BL21 for recombinant '
            'protein production at 37\u00b0C. The Rushton impeller (200 W) dramatically '
            'improves kLa, enabling much higher biomass than a bubble column alone. '
            'Starting from a higher inoculum (OD\u2248 1), the culture reaches high '
            'cell density while the combined agitation + sparging maintains '
            'dissolved oxygen above the critical threshold for aerobic metabolism.'
        ),
        'config': {
            'reactor_type': 'stirred_tank',
            'volume_L': 50.0,
            'diameter_m': 0.35,
            'liquid_height_m': 0.55,
            'gas_flow_rate_Lpm': 5.0,
            'temperature_K': 310.15,
            'pressure_atm': 1.0,
            'mean_bubble_diameter_mm': 2.0,
            'initial_biomass_gL': 0.5,
            'initial_do_mgL': 7.0,
            'initial_dco2_mgL': 0.5,
            'max_growth_rate_per_h': 0.65,
            'ks_oxygen_mgL': 0.1,
            'yield_biomass_o2': 0.9,
            'maintenance_coeff_per_h': 0.02,
            'respiratory_quotient': 1.0,
            'impeller_power_W': 200.0,
        },
        'n_snapshots': 30,
        'total_time': 14.0,
        'color_scheme': 'emerald',
        'reactor_color': '#10b981',
    },
    {
        'id': 'ecoli_ind',
        'title': 'E. coli Industrial Scale-Up',
        'subtitle': '10,000 L airlift — mass transfer bottleneck at scale',
        'description': (
            'An industrial airlift reactor (10,000 L) scaling up E. coli '
            'fermentation at 37\u00b0C. Despite 50 L/min of air, the large vessel '
            'diameter yields low superficial gas velocity and poor kLa. The '
            'culture quickly becomes oxygen-limited, demonstrating the classic '
            'scale-up problem: mass transfer that was adequate at bench scale '
            'becomes the bottleneck at production scale. Growth rate collapses '
            'as DO drops below the critical concentration.'
        ),
        'config': {
            'reactor_type': 'airlift',
            'volume_L': 10000.0,
            'diameter_m': 1.8,
            'liquid_height_m': 4.0,
            'gas_flow_rate_Lpm': 50.0,
            'temperature_K': 310.15,
            'pressure_atm': 1.0,
            'mean_bubble_diameter_mm': 4.0,
            'initial_biomass_gL': 0.5,
            'initial_do_mgL': 7.0,
            'initial_dco2_mgL': 0.5,
            'max_growth_rate_per_h': 0.65,
            'ks_oxygen_mgL': 0.1,
            'yield_biomass_o2': 0.8,
            'maintenance_coeff_per_h': 0.02,
            'respiratory_quotient': 1.05,
        },
        'n_snapshots': 30,
        'total_time': 16.0,
        'color_scheme': 'rose',
        'reactor_color': '#f43f5e',
    },
]

COLOR_SCHEMES = {
    'indigo': {'primary': '#6366f1', 'light': '#e0e7ff', 'dark': '#4338ca',
               'bg': '#eef2ff', 'accent': '#818cf8', 'text': '#312e81'},
    'emerald': {'primary': '#10b981', 'light': '#d1fae5', 'dark': '#059669',
                'bg': '#ecfdf5', 'accent': '#34d399', 'text': '#064e3b'},
    'rose': {'primary': '#f43f5e', 'light': '#ffe4e6', 'dark': '#e11d48',
             'bg': '#fff1f2', 'accent': '#fb7185', 'text': '#881337'},
}


# ── Reactor mesh generation ─────────────────────────────────────────

def generate_reactor_mesh(reactor_type, diameter, height, n_theta=32, n_z=20):
    """Generate a triangulated 3D reactor vessel mesh for visualization.

    Returns (vertices, faces) where vertices is list of [x,y,z] and
    faces is list of [i,j,k] triangle indices.
    """
    r = diameter / 2.0
    vertices = []
    faces = []

    # Cylinder shell
    for j in range(n_z + 1):
        z = j * height / n_z
        for i in range(n_theta):
            theta = 2.0 * math.pi * i / n_theta
            vertices.append([
                round(r * math.cos(theta), 6),
                round(r * math.sin(theta), 6),
                round(z, 6),
            ])

    # Shell faces
    for j in range(n_z):
        for i in range(n_theta):
            i1 = j * n_theta + i
            i2 = j * n_theta + (i + 1) % n_theta
            i3 = (j + 1) * n_theta + i
            i4 = (j + 1) * n_theta + (i + 1) % n_theta
            faces.append([i1, i3, i2])
            faces.append([i2, i3, i4])

    # Bottom cap
    bc = len(vertices)
    vertices.append([0.0, 0.0, 0.0])
    for i in range(n_theta):
        i1 = i
        i2 = (i + 1) % n_theta
        faces.append([bc, i2, i1])

    # Top cap (open — just a rim)
    tc = len(vertices)
    vertices.append([0.0, 0.0, round(height, 6)])
    for i in range(n_theta):
        i1 = n_z * n_theta + i
        i2 = n_z * n_theta + (i + 1) % n_theta
        faces.append([tc, i1, i2])

    # Sparger disk at z = 5% height
    sparger_r = r * 0.35
    sparger_z = round(0.05 * height, 6)
    n_sparger = 16
    sc = len(vertices)
    vertices.append([0.0, 0.0, sparger_z])
    for i in range(n_sparger):
        theta = 2.0 * math.pi * i / n_sparger
        vertices.append([
            round(sparger_r * math.cos(theta), 6),
            round(sparger_r * math.sin(theta), 6),
            sparger_z,
        ])
    for i in range(n_sparger):
        i1 = sc + 1 + i
        i2 = sc + 1 + (i + 1) % n_sparger
        faces.append([sc, i1, i2])

    # Draft tube for airlift reactors
    if reactor_type == 'airlift':
        dr = r * 0.55
        dz_bot = 0.15 * height
        dz_top = 0.85 * height
        n_draft = 16
        n_dz = 8
        d_base = len(vertices)
        for j in range(n_dz + 1):
            z = dz_bot + j * (dz_top - dz_bot) / n_dz
            for i in range(n_draft):
                theta = 2.0 * math.pi * i / n_draft
                vertices.append([
                    round(dr * math.cos(theta), 6),
                    round(dr * math.sin(theta), 6),
                    round(z, 6),
                ])
        for j in range(n_dz):
            for i in range(n_draft):
                i1 = d_base + j * n_draft + i
                i2 = d_base + j * n_draft + (i + 1) % n_draft
                i3 = d_base + (j + 1) * n_draft + i
                i4 = d_base + (j + 1) * n_draft + (i + 1) % n_draft
                faces.append([i1, i2, i3])
                faces.append([i2, i4, i3])

    # Impeller for stirred tanks
    if reactor_type == 'stirred_tank':
        imp_z = round(0.35 * height, 6)
        imp_r = r * 0.4
        shaft_r = r * 0.04
        n_blade = 6
        blade_w = r * 0.08
        # Shaft (small cylinder)
        n_sh = 8
        sh_base = len(vertices)
        for j in [0, 1]:
            z = imp_z if j == 0 else round(height * 0.95, 6)
            for i in range(n_sh):
                theta = 2.0 * math.pi * i / n_sh
                vertices.append([
                    round(shaft_r * math.cos(theta), 6),
                    round(shaft_r * math.sin(theta), 6),
                    round(z, 6),
                ])
        for i in range(n_sh):
            i1 = sh_base + i
            i2 = sh_base + (i + 1) % n_sh
            i3 = sh_base + n_sh + i
            i4 = sh_base + n_sh + (i + 1) % n_sh
            faces.append([i1, i3, i2])
            faces.append([i2, i3, i4])
        # Blades (flat quads)
        for b in range(n_blade):
            theta = 2.0 * math.pi * b / n_blade
            ct, st_ = math.cos(theta), math.sin(theta)
            # Normal perpendicular to radial direction
            nx, ny = -st_, ct
            bl_base = len(vertices)
            for ri in [shaft_r, imp_r]:
                for dz in [-blade_w, blade_w]:
                    vertices.append([
                        round(ri * ct, 6),
                        round(ri * st_, 6),
                        round(imp_z + dz, 6),
                    ])
            faces.append([bl_base, bl_base + 1, bl_base + 3])
            faces.append([bl_base, bl_base + 3, bl_base + 2])

    return vertices, faces


# ── Simulation ──────────────────────────────────────────────────────

def run_simulation(cfg_entry):
    """Run a single bioreactor simulation, returning snapshots and runtime."""
    core = allocate_core()
    core.register_link('BiRDReactorProcess', BiRDReactorProcess)

    t0 = _time.perf_counter()
    proc = BiRDReactorProcess(config=cfg_entry['config'], core=core)
    state0 = proc.initial_state()

    interval = cfg_entry['total_time'] / cfg_entry['n_snapshots']
    snapshots = [_snap(0.0, state0)]

    t = 0.0
    for i in range(cfg_entry['n_snapshots']):
        result = proc.update({}, interval=interval)
        t += interval
        if not result:
            break
        snapshots.append(_snap(round(t, 4), result))

    runtime = _time.perf_counter() - t0
    return snapshots, runtime


def _snap(t, s):
    return {
        'time': t,
        'dissolved_o2': s['dissolved_o2'],
        'dissolved_co2': s['dissolved_co2'],
        'biomass': s['biomass'],
        'gas_holdup': s['gas_holdup'],
        'kla_o2': s['kla_o2'],
        'kla_co2': s['kla_co2'],
        'o2_saturation': s['o2_saturation'],
        'specific_growth_rate': s['specific_growth_rate'],
        'o2_uptake_rate': s['o2_uptake_rate'],
        'co2_evolution_rate': s['co2_evolution_rate'],
        'superficial_gas_velocity': s['superficial_gas_velocity'],
    }


# ── Bigraph diagram ────────────────────────────────────────────────

def generate_bigraph_image(cfg_entry):
    """Generate a colored bigraph-viz PNG for the composite document."""
    from bigraph_viz import plot_bigraph

    doc = {
        'reactor': {
            '_type': 'process',
            'address': 'local:BiRDReactorProcess',
            'config': {k: v for k, v in cfg_entry['config'].items()
                       if k in ('reactor_type', 'volume_L', 'temperature_K')},
            'outputs': {
                'dissolved_o2': ['stores', 'dissolved_o2'],
                'biomass': ['stores', 'biomass'],
                'gas_holdup': ['stores', 'gas_holdup'],
                'kla_o2': ['stores', 'kla_o2'],
            },
        },
        'stores': {},
        'emitter': {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'inputs': {
                'dissolved_o2': ['stores', 'dissolved_o2'],
                'biomass': ['stores', 'biomass'],
                'time': ['global_time'],
            },
        },
    }

    node_colors = {
        ('reactor',): '#6366f1',
        ('emitter',): '#8b5cf6',
        ('stores',): '#e0e7ff',
    }

    outdir = tempfile.mkdtemp()
    plot_bigraph(
        state=doc,
        out_dir=outdir,
        filename='bigraph',
        file_format='png',
        remove_process_place_edges=True,
        rankdir='LR',
        node_fill_colors=node_colors,
        node_label_size='16pt',
        port_labels=False,
        dpi='150',
    )
    png_path = os.path.join(outdir, 'bigraph.png')
    with open(png_path, 'rb') as f:
        b64 = base64.b64encode(f.read()).decode()
    return f'data:image/png;base64,{b64}'


def build_pbg_document(cfg_entry):
    """Build the PBG composite document dict for display."""
    doc = make_reactor_document(**{
        k: v for k, v in cfg_entry['config'].items()
        if k in ('reactor_type', 'volume_L', 'diameter_m', 'liquid_height_m',
                 'gas_flow_rate_Lpm', 'temperature_K', 'mean_bubble_diameter_mm',
                 'initial_biomass_gL', 'max_growth_rate_per_h', 'impeller_power_W')
    }, interval=cfg_entry['total_time'] / cfg_entry['n_snapshots'])
    return doc


# ── HTML Generation ─────────────────────────────────────────────────

def generate_html(sim_results, output_path):
    """Generate self-contained HTML report."""

    sections_html = []
    all_js_data = {}

    for idx, (cfg, (snapshots, runtime)) in enumerate(sim_results):
        sid = cfg['id']
        cs = COLOR_SCHEMES[cfg['color_scheme']]
        rcfg = cfg['config']

        # Generate reactor mesh
        mesh_verts, mesh_faces = generate_reactor_mesh(
            rcfg['reactor_type'], rcfg['diameter_m'], rcfg['liquid_height_m'])

        # Time series extraction
        times = [s['time'] for s in snapshots]
        do2 = [s['dissolved_o2'] for s in snapshots]
        dco2 = [s['dissolved_co2'] for s in snapshots]
        biomass = [s['biomass'] for s in snapshots]
        holdup = [s['gas_holdup'] for s in snapshots]
        kla = [s['kla_o2'] for s in snapshots]
        mu = [s['specific_growth_rate'] for s in snapshots]
        our = [s['o2_uptake_rate'] for s in snapshots]
        cer = [s['co2_evolution_rate'] for s in snapshots]
        cstar = [s['o2_saturation'] for s in snapshots]

        # JS data — all variables available for 3D coloring and charts
        all_js_data[sid] = {
            'faces': mesh_faces,
            'vertices': mesh_verts,
            'n_snapshots': len(snapshots),
            'camera': [
                rcfg['diameter_m'] * 2.5,
                rcfg['diameter_m'] * 1.8,
                rcfg['liquid_height_m'] * 1.5,
            ],
            'height': rcfg['liquid_height_m'],
            # Color variables: each has values array and [min, max] range
            'color_vars': {
                'dissolved_o2': {
                    'values': do2, 'range': [min(do2), max(do2)],
                    'label': 'Dissolved O\u2082', 'unit': 'mg/L',
                    'high_color': 'good', # blue=high is good
                },
                'dissolved_co2': {
                    'values': dco2, 'range': [min(dco2), max(dco2)],
                    'label': 'Dissolved CO\u2082', 'unit': 'mg/L',
                    'high_color': 'bad', # red=high is bad (accumulation)
                },
                'biomass': {
                    'values': biomass, 'range': [min(biomass), max(biomass)],
                    'label': 'Biomass', 'unit': 'g/L',
                    'high_color': 'good',
                },
                'gas_holdup': {
                    'values': [h * 100 for h in holdup],
                    'range': [min(holdup) * 100, max(holdup) * 100],
                    'label': 'Gas Holdup', 'unit': '%',
                    'high_color': 'neutral',
                },
                'growth_rate': {
                    'values': mu, 'range': [min(mu), max(mu)],
                    'label': 'Growth Rate', 'unit': 'h\u207b\u00b9',
                    'high_color': 'good',
                },
                'kla_o2': {
                    'values': kla, 'range': [min(kla), max(kla)],
                    'label': 'kLa (O\u2082)', 'unit': 'h\u207b\u00b9',
                    'high_color': 'good',
                },
            },
            'charts': {
                'times': times,
                'dissolved_o2': do2,
                'dissolved_co2': dco2,
                'biomass': biomass,
                'gas_holdup': [h * 100 for h in holdup],
                'kla_o2': kla,
                'growth_rate': mu,
                'o2_uptake': our,
                'co2_evolution': cer,
                'o2_saturation': cstar,
            },
        }

        # Bigraph image
        print(f'  Generating bigraph diagram for {sid}...')
        bigraph_img = generate_bigraph_image(cfg)

        # PBG document
        pbg_doc = build_pbg_document(cfg)

        # Metrics
        do_final = do2[-1]
        x_final = biomass[-1]
        x_init = biomass[0]
        x_fold = x_final / x_init if x_init > 0 else 0
        kla_avg = sum(kla) / len(kla)
        holdup_pct = holdup[-1] * 100

        section = f"""
    <div class="sim-section" id="sim-{sid}">
      <div class="sim-header" style="border-left: 4px solid {cs['primary']};">
        <div class="sim-number" style="background:{cs['light']}; color:{cs['dark']};">{idx+1}</div>
        <div>
          <h2 class="sim-title">{cfg['title']}</h2>
          <p class="sim-subtitle">{cfg['subtitle']}</p>
        </div>
      </div>
      <p class="sim-description">{cfg['description']}</p>

      <div class="metrics-row">
        <div class="metric"><span class="metric-label">Reactor</span><span class="metric-value">{rcfg['reactor_type'].replace('_',' ').title()}</span></div>
        <div class="metric"><span class="metric-label">Volume</span><span class="metric-value">{rcfg['volume_L']:.0f} L</span></div>
        <div class="metric"><span class="metric-label">Final DO</span><span class="metric-value">{do_final:.2f} mg/L</span></div>
        <div class="metric"><span class="metric-label">Biomass</span><span class="metric-value">{x_fold:.1f}x</span><span class="metric-sub">{x_init:.2f} &rarr; {x_final:.2f} g/L</span></div>
        <div class="metric"><span class="metric-label">Avg kLa</span><span class="metric-value">{kla_avg:.1f} h\u207b\u00b9</span></div>
        <div class="metric"><span class="metric-label">Gas Holdup</span><span class="metric-value">{holdup_pct:.2f}%</span></div>
        <div class="metric"><span class="metric-label">Runtime</span><span class="metric-value">{runtime:.2f}s</span></div>
      </div>

      <h3 class="subsection-title">3D Reactor Geometry</h3>
      <div class="viewer-wrap">
        <div class="var-selector-row">
          <label class="var-label">Color by:</label>
          <select id="varsel-{sid}" class="var-select" onchange="changeColorVar('{sid}', this.value)" style="border-color:{cs['primary']}">
            <option value="dissolved_o2" selected>Dissolved O\u2082 (mg/L)</option>
            <option value="dissolved_co2">Dissolved CO\u2082 (mg/L)</option>
            <option value="biomass">Biomass (g/L)</option>
            <option value="gas_holdup">Gas Holdup (%)</option>
            <option value="growth_rate">Growth Rate (h\u207b\u00b9)</option>
            <option value="kla_o2">kLa O\u2082 (h\u207b\u00b9)</option>
          </select>
        </div>
        <div class="viewer-row">
          <canvas id="canvas-{sid}" class="mesh-canvas"></canvas>
          <div class="colorbar-wrap" id="cbar-{sid}">
            <div class="colorbar-label" id="cbar-max-{sid}">{max(do2):.1f}</div>
            <div class="colorbar-gradient" id="cbar-grad-{sid}"></div>
            <div class="colorbar-label" id="cbar-min-{sid}">{min(do2):.1f}</div>
            <div class="colorbar-unit" id="cbar-unit-{sid}">mg/L</div>
          </div>
        </div>
        <div class="state-readout" id="readout-{sid}"></div>
        <div class="slider-row">
          <button class="play-btn" id="playbtn-{sid}" onclick="togglePlay('{sid}')">&#9654;</button>
          <input type="range" id="slider-{sid}" min="0" max="{len(snapshots)-1}" value="0" class="time-slider" style="accent-color:{cs['primary']}">
          <span id="tval-{sid}" class="time-val">t = 0 h</span>
        </div>
        <div class="viewer-info">
          Reactor vessel colored by selected variable.
          <span style="color:#3b82f6">Blue = high</span> &rarr;
          <span style="color:#ef4444">Red = low</span>
          (inverted for CO\u2082: red = high accumulation).
          Drag to rotate. Scroll to zoom.
        </div>
      </div>

      <h3 class="subsection-title">Time-Series Dynamics</h3>
      <div class="charts-grid">
        <div id="chart-do2-{sid}" class="chart-box"></div>
        <div id="chart-co2-{sid}" class="chart-box"></div>
        <div id="chart-biomass-{sid}" class="chart-box"></div>
        <div id="chart-kla-{sid}" class="chart-box"></div>
        <div id="chart-rates-{sid}" class="chart-box"></div>
        <div id="chart-holdup-{sid}" class="chart-box"></div>
      </div>

      <h3 class="subsection-title">Architecture &amp; Document</h3>
      <div class="arch-grid">
        <div class="bigraph-card">
          <h4>Bigraph Architecture</h4>
          <div class="bigraph-img-wrap">
            <img src="{bigraph_img}" alt="Bigraph architecture diagram">
          </div>
        </div>
        <div class="json-card">
          <h4>PBG Composite Document</h4>
          <div class="json-tree" id="json-{sid}"></div>
        </div>
      </div>
    </div>"""
        sections_html.append(section)

        # Embed JSON doc for JS tree viewer
        all_js_data[sid]['pbg_doc'] = pbg_doc

    # Navigation
    nav_items = ''.join(
        f'<a href="#sim-{c["id"]}" class="nav-link" '
        f'style="border-bottom:2px solid {COLOR_SCHEMES[c["color_scheme"]]["primary"]}">'
        f'{c["title"]}</a>'
        for c in CONFIGS
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>BiRD Bioreactor Simulation Report</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
*{{margin:0;padding:0;box-sizing:border-box;}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;background:#fff;color:#1e293b;line-height:1.6;}}
.top-bar{{position:sticky;top:0;z-index:100;background:rgba(255,255,255,0.95);backdrop-filter:blur(8px);border-bottom:1px solid #e2e8f0;padding:0.6rem 2rem;display:flex;align-items:center;gap:1.5rem;}}
.top-title{{font-size:1.1rem;font-weight:700;color:#1e293b;white-space:nowrap;}}
.nav-links{{display:flex;gap:1rem;overflow-x:auto;}}
.nav-link{{font-size:0.82rem;color:#64748b;text-decoration:none;padding:0.3rem 0;white-space:nowrap;transition:color 0.2s;}}
.nav-link:hover{{color:#1e293b;}}
.container{{max-width:1200px;margin:0 auto;padding:2rem;}}
.hero{{text-align:center;padding:2.5rem 1rem 1.5rem;}}
.hero h1{{font-size:2rem;font-weight:800;background:linear-gradient(135deg,#6366f1,#10b981,#f43f5e);-webkit-background-clip:text;-webkit-text-fill-color:transparent;background-clip:text;}}
.hero p{{color:#64748b;font-size:0.95rem;margin-top:0.5rem;max-width:700px;margin-left:auto;margin-right:auto;}}
.sim-section{{margin:3rem 0;padding:2rem 0;border-top:1px solid #e2e8f0;}}
.sim-header{{display:flex;align-items:center;gap:1rem;padding-left:1rem;}}
.sim-number{{width:36px;height:36px;border-radius:50%;display:flex;align-items:center;justify-content:center;font-weight:700;font-size:1.1rem;flex-shrink:0;}}
.sim-title{{font-size:1.4rem;font-weight:700;}}
.sim-subtitle{{font-size:0.85rem;color:#64748b;}}
.sim-description{{margin:1rem 0 1.5rem;color:#475569;font-size:0.9rem;max-width:800px;}}
.metrics-row{{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:0.7rem;margin:1.5rem 0;}}
.metric{{background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;padding:0.7rem;text-align:center;}}
.metric-label{{display:block;font-size:0.7rem;color:#94a3b8;text-transform:uppercase;letter-spacing:0.05em;}}
.metric-value{{display:block;font-size:1.05rem;font-weight:700;color:#1e293b;margin-top:0.15rem;}}
.metric-sub{{display:block;font-size:0.65rem;color:#94a3b8;}}
.subsection-title{{font-size:1rem;font-weight:600;color:#334155;margin:2rem 0 0.8rem;}}
.viewer-wrap{{background:#f8fafc;border:1px solid #e2e8f0;border-radius:12px;padding:1rem;}}
.viewer-row{{display:flex;gap:0.8rem;align-items:stretch;}}
.mesh-canvas{{flex:1;height:400px;border-radius:8px;display:block;min-width:0;}}
.colorbar-wrap{{display:flex;flex-direction:column;align-items:center;justify-content:space-between;width:50px;padding:0.5rem 0;flex-shrink:0;}}
.colorbar-gradient{{width:18px;flex:1;margin:0.3rem 0;border-radius:4px;border:1px solid #e2e8f0;background:linear-gradient(to bottom,#3b82f6,#06b6d4,#22c55e,#eab308,#ef4444);}}
.colorbar-label{{font-size:0.65rem;color:#64748b;font-family:'SF Mono',Menlo,monospace;text-align:center;line-height:1.2;}}
.var-selector-row{{display:flex;align-items:center;gap:0.6rem;margin-bottom:0.6rem;}}
.var-label{{font-size:0.8rem;color:#64748b;font-weight:500;}}
.var-select{{font-size:0.8rem;padding:0.3rem 0.6rem;border-radius:6px;border:1.5px solid #cbd5e1;background:#fff;color:#334155;cursor:pointer;font-family:inherit;}}
.var-select:focus{{outline:none;border-color:#6366f1;box-shadow:0 0 0 2px rgba(99,102,241,0.15);}}
.state-readout{{display:flex;flex-wrap:wrap;justify-content:center;gap:0.5rem 1.2rem;font-size:0.78rem;color:#475569;padding:0.5rem 0 0.2rem;font-family:'SF Mono',Menlo,monospace;}}
.state-readout .rv{{white-space:nowrap;}}
.state-readout .rv-label{{color:#94a3b8;}}
.state-readout .rv-val{{font-weight:600;color:#1e293b;}}
.colorbar-unit{{font-size:0.6rem;color:#94a3b8;margin-top:0.2rem;white-space:nowrap;}}
.viewer-info{{text-align:center;font-size:0.75rem;color:#94a3b8;margin:0.5rem 0;}}
.slider-row{{display:flex;align-items:center;gap:0.8rem;padding:0.5rem 0.5rem 0;}}
.play-btn{{width:32px;height:32px;border-radius:50%;border:1px solid #cbd5e1;background:#fff;cursor:pointer;font-size:0.8rem;display:flex;align-items:center;justify-content:center;flex-shrink:0;}}
.play-btn:hover{{background:#f1f5f9;}}
.time-slider{{flex:1;height:4px;cursor:pointer;}}
.time-val{{font-size:0.8rem;color:#64748b;font-family:'SF Mono',Menlo,monospace;min-width:80px;text-align:right;}}
.charts-grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem;}}
.charts-grid .chart-box{{min-height:260px;}}
.chart-box{{background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:0.5rem;min-height:280px;}}
.arch-grid{{display:grid;grid-template-columns:1fr 1fr;gap:1rem;}}
.bigraph-card,.json-card{{background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:1.2rem;}}
.bigraph-card h4,.json-card h4{{font-size:0.85rem;color:#64748b;margin-bottom:0.8rem;}}
.bigraph-img-wrap{{text-align:center;}}
.bigraph-img-wrap img{{max-width:100%;height:auto;}}
.json-tree{{font-family:'SF Mono',Menlo,Monaco,monospace;font-size:0.75rem;line-height:1.5;max-height:400px;overflow:auto;}}
.jt-key{{color:#7c3aed;font-weight:600;}}
.jt-str{{color:#059669;}}
.jt-num{{color:#2563eb;}}
.jt-bool{{color:#d97706;}}
.jt-null{{color:#94a3b8;}}
.jt-bracket{{color:#94a3b8;}}
.jt-toggle{{cursor:pointer;user-select:none;margin-right:0.3rem;color:#94a3b8;font-size:0.65rem;}}
.jt-toggle:hover{{color:#475569;}}
.jt-collapsed{{display:none;}}
.footer{{text-align:center;padding:2rem 0;color:#94a3b8;font-size:0.8rem;border-top:1px solid #e2e8f0;margin-top:3rem;}}
@media(max-width:900px){{
  .charts-grid,.arch-grid{{grid-template-columns:1fr;}}
  .top-bar{{flex-wrap:wrap;padding:0.5rem 1rem;}}
  .viewer-row{{flex-direction:column;}}
  .colorbar-wrap{{flex-direction:row;width:100%;height:40px;}}
  .colorbar-gradient{{width:auto;height:14px;flex:1;background:linear-gradient(to right,#3b82f6,#06b6d4,#22c55e,#eab308,#ef4444);}}
}}
</style>
</head>
<body>

<div class="top-bar">
  <span class="top-title">BiRD Bioreactor Report</span>
  <div class="nav-links">{nav_items}</div>
</div>

<div class="container">
  <div class="hero">
    <h1>BioReactorDesign Simulation Report</h1>
    <p>Process-bigraph wrapper for BiRD: 0D bioreactor models with Higbie kLa,
    Monod kinetics, and temperature-dependent Henry's law. Three E. coli
    fermentation scenarios — from 2 L bench to 10,000 L production —
    demonstrating scale-up mass transfer challenges.</p>
  </div>

  {''.join(sections_html)}

  <div class="footer">
    Generated by <strong>viva-bioreactordesign</strong> &middot;
    Process-bigraph wrapper for BioReactorDesign (BiRD) &middot;
    <a href="https://github.com/NatLabRockies/BioReactorDesign" style="color:#6366f1">BiRD on GitHub</a>
  </div>
</div>

<script>
// ── Data ───────────────────────────────────────────────────────────
const DATA = {json.dumps(all_js_data)};

// ── JSON tree viewer ───────────────────────────────────────────────
function renderJson(obj, depth) {{
  if (depth === undefined) depth = 0;
  if (obj === null) return '<span class="jt-null">null</span>';
  if (typeof obj === 'boolean') return '<span class="jt-bool">' + obj + '</span>';
  if (typeof obj === 'number') return '<span class="jt-num">' + (Number.isFinite(obj) ? (Math.abs(obj) < 0.001 && obj !== 0 ? obj.toExponential(2) : (Number.isInteger(obj) ? obj : obj.toFixed(4))) : obj) + '</span>';
  if (typeof obj === 'string') return '<span class="jt-str">"' + obj.replace(/</g,'&lt;') + '"</span>';
  if (Array.isArray(obj)) {{
    if (obj.length === 0) return '<span class="jt-bracket">[]</span>';
    if (obj.length <= 5 && obj.every(function(x){{ return typeof x !== 'object' || x === null; }})) {{
      return '<span class="jt-bracket">[</span>' + obj.map(function(x){{ return renderJson(x, depth+1); }}).join(', ') + '<span class="jt-bracket">]</span>';
    }}
    var id = 'jt' + Math.random().toString(36).slice(2,9);
    var collapsed = depth >= 2;
    var h = '<span class="jt-toggle" onclick="toggleJt(this,\\''+id+'\\')">'+( collapsed ? '\\u25b6' : '\\u25bc')+'</span>';
    h += '<span class="jt-bracket">[</span> <span style="color:#94a3b8">' + obj.length + ' items</span>';
    h += '<div id="'+id+'"'+(collapsed?' class="jt-collapsed"':'')+' style="margin-left:1.2rem">';
    obj.forEach(function(v,i){{ h += '<div>' + renderJson(v,depth+1) + (i<obj.length-1?',':'') + '</div>'; }});
    h += '</div><span class="jt-bracket">]</span>';
    return h;
  }}
  if (typeof obj === 'object') {{
    var keys = Object.keys(obj);
    if (keys.length === 0) return '<span class="jt-bracket">{{}}</span>';
    var id = 'jt' + Math.random().toString(36).slice(2,9);
    var collapsed = depth >= 2;
    var h = '<span class="jt-toggle" onclick="toggleJt(this,\\''+id+'\\')">'+( collapsed ? '\\u25b6' : '\\u25bc')+'</span>';
    h += '<span class="jt-bracket">{{</span>';
    h += '<div id="'+id+'"'+(collapsed?' class="jt-collapsed"':'')+' style="margin-left:1.2rem">';
    keys.forEach(function(k,i){{ h += '<div><span class="jt-key">'+k+'</span>: '+renderJson(obj[k],depth+1)+(i<keys.length-1?',':'')+'</div>'; }});
    h += '</div><span class="jt-bracket">}}</span>';
    return h;
  }}
  return String(obj);
}}
function toggleJt(el, id) {{
  var c = document.getElementById(id);
  if (c.classList.contains('jt-collapsed')) {{
    c.classList.remove('jt-collapsed');
    el.innerHTML = '\\u25bc';
  }} else {{
    c.classList.add('jt-collapsed');
    el.innerHTML = '\\u25b6';
  }}
}}

// ── Three.js reactor viewers ───────────────────────────────────────
var viewers = {{}};
var playStates = {{}};

function colormap(t) {{
  // Blue -> Cyan -> Green -> Yellow -> Red (t=0 is blue, t=1 is red)
  var r, g, b;
  if (t < 0.25) {{
    var s = t / 0.25;
    r = 0.19; g = 0.07 + 0.63*s; b = 0.99 - 0.19*s;
  }} else if (t < 0.5) {{
    var s = (t-0.25)/0.25;
    r = 0.19 + 0.11*s; g = 0.70 + 0.15*s; b = 0.80 - 0.55*s;
  }} else if (t < 0.75) {{
    var s = (t-0.5)/0.25;
    r = 0.30 + 0.60*s; g = 0.85 - 0.10*s; b = 0.25 - 0.15*s;
  }} else {{
    var s = (t-0.75)/0.25;
    r = 0.90 + 0.10*s; g = 0.75 - 0.55*s; b = 0.10 - 0.05*s;
  }}
  return [r, g, b];
}}

// Track which variable is selected for coloring per simulation
var activeColorVar = {{}};

function getColor(val, vmin, vmax, invert) {{
  // invert=true: high value = blue (good), low = red
  // invert=false: high value = red (bad), low = blue
  var t = Math.max(0, Math.min(1, (val - vmin) / (vmax - vmin + 1e-12)));
  if (invert) t = 1.0 - t;
  return colormap(t);
}}

function initViewer(sid) {{
  var d = DATA[sid];
  var canvas = document.getElementById('canvas-' + sid);
  if (!canvas) return;
  var W = canvas.clientWidth, H = canvas.clientHeight || 400;
  var renderer = new THREE.WebGLRenderer({{canvas:canvas, antialias:true, alpha:true}});
  renderer.setSize(W, H);
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setClearColor(0xf8fafc);

  var scene = new THREE.Scene();
  var cam = new THREE.PerspectiveCamera(45, W/H, 0.001, 100);
  cam.position.set(d.camera[0], d.camera[1], d.camera[2]);
  cam.lookAt(0, 0, d.height / 2);

  var controls = new THREE.OrbitControls(cam, canvas);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.autoRotate = true;
  controls.autoRotateSpeed = 0.8;
  controls.target.set(0, 0, d.height / 2);

  scene.add(new THREE.AmbientLight(0xffffff, 0.5));
  var dl1 = new THREE.DirectionalLight(0xffffff, 0.7);
  dl1.position.set(3, 5, 4);
  scene.add(dl1);
  var dl2 = new THREE.DirectionalLight(0xcbd5e1, 0.4);
  dl2.position.set(-3, -2, -4);
  scene.add(dl2);

  var nv = d.vertices.length;
  var positions = new Float32Array(nv * 3);
  var colors = new Float32Array(nv * 3);

  for (var i = 0; i < nv; i++) {{
    positions[i*3]   = d.vertices[i][0];
    positions[i*3+1] = d.vertices[i][1];
    positions[i*3+2] = d.vertices[i][2];
  }}

  var geometry = new THREE.BufferGeometry();
  var indices = [];
  for (var i = 0; i < d.faces.length; i++) {{
    indices.push(d.faces[i][0], d.faces[i][1], d.faces[i][2]);
  }}
  geometry.setIndex(indices);
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

  var material = new THREE.MeshPhongMaterial({{
    vertexColors: true, side: THREE.DoubleSide,
    shininess: 40, specular: 0xdddddd, flatShading: false,
    transparent: true, opacity: 0.85
  }});
  scene.add(new THREE.Mesh(geometry, material));

  var wireMat = new THREE.MeshBasicMaterial({{
    color: 0x94a3b8, wireframe: true, transparent: true, opacity: 0.06
  }});
  scene.add(new THREE.Mesh(geometry, wireMat));

  activeColorVar[sid] = 'dissolved_o2';

  function updateColor(snapIdx) {{
    var varKey = activeColorVar[sid];
    var cv = d.color_vars[varKey];
    var val = cv.values[snapIdx];
    var invert = (cv.high_color === 'good');
    var c = getColor(val, cv.range[0], cv.range[1], invert);
    for (var i = 0; i < nv; i++) {{
      colors[i*3]   = c[0];
      colors[i*3+1] = c[1];
      colors[i*3+2] = c[2];
    }}
    geometry.attributes.color.needsUpdate = true;
    geometry.computeVertexNormals();
  }}
  updateColor(0);

  var slider = document.getElementById('slider-' + sid);
  var tval = document.getElementById('tval-' + sid);
  var readout = document.getElementById('readout-' + sid);

  function fmtVal(v) {{
    if (Math.abs(v) >= 100) return v.toFixed(1);
    if (Math.abs(v) >= 1) return v.toFixed(2);
    if (Math.abs(v) >= 0.01) return v.toFixed(3);
    return v.toExponential(2);
  }}

  function updateReadout(idx) {{
    updateColor(idx);
    tval.textContent = 't = ' + d.charts.times[idx].toFixed(1) + ' h';
    if (readout) {{
      var c = d.charts;
      readout.innerHTML =
        '<span class="rv"><span class="rv-label">DO </span><span class="rv-val">' + fmtVal(c.dissolved_o2[idx]) + '</span> mg/L</span>' +
        '<span class="rv"><span class="rv-label">CO\u2082 </span><span class="rv-val">' + fmtVal(c.dissolved_co2[idx]) + '</span> mg/L</span>' +
        '<span class="rv"><span class="rv-label">Biomass </span><span class="rv-val">' + fmtVal(c.biomass[idx]) + '</span> g/L</span>' +
        '<span class="rv"><span class="rv-label">\u03bc </span><span class="rv-val">' + fmtVal(c.growth_rate[idx]) + '</span> h\u207b\u00b9</span>' +
        '<span class="rv"><span class="rv-label">kLa </span><span class="rv-val">' + fmtVal(c.kla_o2[idx]) + '</span> h\u207b\u00b9</span>' +
        '<span class="rv"><span class="rv-label">Gas </span><span class="rv-val">' + fmtVal(c.gas_holdup[idx]) + '</span>%</span>';
    }}
  }}
  updateReadout(0);

  slider.addEventListener('input', function() {{
    updateReadout(parseInt(slider.value));
  }});

  function animate() {{
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, cam);
  }}
  animate();

  viewers[sid] = {{slider:slider, tval:tval, updateColor:updateColor, updateReadout:updateReadout, controls:controls}};
  playStates[sid] = {{playing:false, interval:null}};
}}

function changeColorVar(sid, varKey) {{
  var d = DATA[sid];
  var cv = d.color_vars[varKey];
  activeColorVar[sid] = varKey;
  // Update color bar labels
  var invert = (cv.high_color === 'good');
  document.getElementById('cbar-max-' + sid).textContent = (invert ? cv.range[1] : cv.range[0]).toFixed(2);
  document.getElementById('cbar-min-' + sid).textContent = (invert ? cv.range[0] : cv.range[1]).toFixed(2);
  document.getElementById('cbar-unit-' + sid).textContent = cv.unit;
  // Update gradient direction: for 'good' vars, blue on top (high=good)
  var grad = document.getElementById('cbar-grad-' + sid);
  if (invert) {{
    grad.style.background = 'linear-gradient(to bottom,#3b82f6,#06b6d4,#22c55e,#eab308,#ef4444)';
  }} else {{
    grad.style.background = 'linear-gradient(to bottom,#ef4444,#eab308,#22c55e,#06b6d4,#3b82f6)';
  }}
  // Re-color at current slider position
  var idx = parseInt(viewers[sid].slider.value);
  viewers[sid].updateColor(idx);
}}

function togglePlay(sid) {{
  var ps = playStates[sid];
  var v = viewers[sid];
  var d = DATA[sid];
  ps.playing = !ps.playing;
  var btn = document.getElementById('playbtn-' + sid);
  if (ps.playing) {{
    btn.innerHTML = '\\u23f8';
    v.controls.autoRotate = false;
    ps.interval = setInterval(function() {{
      var idx = parseInt(v.slider.value) + 1;
      if (idx >= d.n_snapshots + 1) idx = 0;
      v.slider.value = idx;
      v.updateReadout(idx);
    }}, 400);
  }} else {{
    btn.innerHTML = '\\u25b6';
    clearInterval(ps.interval);
    v.controls.autoRotate = true;
  }}
}}

// ── Plotly charts ──────────────────────────────────────────────────
function initCharts(sid) {{
  var c = DATA[sid].charts;
  var pLayout = {{
    paper_bgcolor:'#f8fafc', plot_bgcolor:'#f8fafc',
    font:{{color:'#64748b',family:'-apple-system,sans-serif',size:11}},
    margin:{{l:55,r:15,t:35,b:40}},
    xaxis:{{gridcolor:'#e2e8f0',title:{{text:'Time (h)',font:{{size:10}}}}}},
    yaxis:{{gridcolor:'#e2e8f0'}},
  }};
  var pCfg = {{responsive:true, displayModeBar:false}};

  // Chart 1: Dissolved O2 + saturation
  Plotly.newPlot('chart-do2-'+sid, [
    {{x:c.times, y:c.dissolved_o2, type:'scatter', mode:'lines+markers',
      line:{{color:'#3b82f6',width:2}}, marker:{{size:3}}, name:'DO (mg/L)'}},
    {{x:c.times, y:c.o2_saturation, type:'scatter', mode:'lines',
      line:{{color:'#93c5fd',width:1.5,dash:'dash'}}, name:'Saturation'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'mg/L',font:{{size:10}}}}}},
    title:{{text:'Dissolved O\\u2082',font:{{size:13}}}},
    legend:{{x:0.7,y:0.95,font:{{size:9}}}},
    showlegend:true
  }}), pCfg);

  // Chart 2: Dissolved CO2
  Plotly.newPlot('chart-co2-'+sid, [
    {{x:c.times, y:c.dissolved_co2, type:'scatter', mode:'lines+markers',
      line:{{color:'#f97316',width:2}}, marker:{{size:3}}, name:'DCO\\u2082',
      fill:'tozeroy', fillcolor:'rgba(249,115,22,0.06)'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'mg/L',font:{{size:10}}}}}},
    title:{{text:'Dissolved CO\\u2082',font:{{size:13}}}},
    showlegend:false
  }}), pCfg);

  // Chart 3: Biomass
  Plotly.newPlot('chart-biomass-'+sid, [
    {{x:c.times, y:c.biomass, type:'scatter', mode:'lines+markers',
      line:{{color:'#10b981',width:2}}, marker:{{size:3}}, name:'Biomass',
      fill:'tozeroy', fillcolor:'rgba(16,185,129,0.06)'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'g/L',font:{{size:10}}}}}},
    title:{{text:'Biomass Growth',font:{{size:13}}}},
    showlegend:false
  }}), pCfg);

  // Chart 3: kLa
  Plotly.newPlot('chart-kla-'+sid, [
    {{x:c.times, y:c.kla_o2, type:'scatter', mode:'lines+markers',
      line:{{color:'#8b5cf6',width:2}}, marker:{{size:3}}, name:'kLa O\\u2082'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'h\\u207b\\u00b9',font:{{size:10}}}}}},
    title:{{text:'kLa (O\\u2082 Mass Transfer)',font:{{size:13}}}},
    showlegend:false
  }}), pCfg);

  // Chart 4: O2 uptake + CO2 evolution rates
  Plotly.newPlot('chart-rates-'+sid, [
    {{x:c.times, y:c.o2_uptake, type:'scatter', mode:'lines+markers',
      line:{{color:'#f43f5e',width:2}}, marker:{{size:3}}, name:'OUR'}},
    {{x:c.times, y:c.co2_evolution, type:'scatter', mode:'lines+markers',
      line:{{color:'#f97316',width:2}}, marker:{{size:3}}, name:'CER'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'mg/(L\\u00b7h)',font:{{size:10}}}}}},
    title:{{text:'O\\u2082 Uptake & CO\\u2082 Evolution',font:{{size:13}}}},
    legend:{{x:0.65,y:0.95,font:{{size:9}}}},
    showlegend:true
  }}), pCfg);

  // Chart 6: Gas Holdup + Growth Rate (dual y-axis)
  Plotly.newPlot('chart-holdup-'+sid, [
    {{x:c.times, y:c.gas_holdup, type:'scatter', mode:'lines+markers',
      line:{{color:'#0ea5e9',width:2}}, marker:{{size:3}}, name:'Gas Holdup (%)'}},
    {{x:c.times, y:c.growth_rate, type:'scatter', mode:'lines+markers',
      line:{{color:'#a855f7',width:2,dash:'dot'}}, marker:{{size:3}}, name:'\\u03bc (h\\u207b\\u00b9)', yaxis:'y2'}},
  ], Object.assign({{}}, pLayout, {{
    yaxis:{{gridcolor:'#e2e8f0',title:{{text:'%',font:{{size:10}}}}}},
    yaxis2:{{overlaying:'y',side:'right',gridcolor:'transparent',title:{{text:'h\\u207b\\u00b9',font:{{size:10}}}}}},
    title:{{text:'Gas Holdup & Growth Rate',font:{{size:13}}}},
    legend:{{x:0.5,y:0.95,font:{{size:9}}}},
    showlegend:true
  }}), pCfg);
}}

// ── Init all ───────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', function() {{
  var sids = {json.dumps([c['id'] for c in CONFIGS])};
  sids.forEach(function(sid) {{
    initViewer(sid);
    initCharts(sid);
    var el = document.getElementById('json-'+sid);
    if (el && DATA[sid].pbg_doc) {{
      el.innerHTML = renderJson(DATA[sid].pbg_doc, 0);
    }}
  }});
}});
</script>
</body>
</html>"""

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(html)
    print(f'Report written to {output_path}')


# ── Main ────────────────────────────────────────────────────────────

def run_demo():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, 'report.html')

    print('BiRD Bioreactor Simulation Report')
    print('=' * 50)

    sim_results = []
    for cfg in CONFIGS:
        print(f'\nRunning: {cfg["title"]}...')
        snapshots, runtime = run_simulation(cfg)
        print(f'  {len(snapshots)} snapshots in {runtime:.2f}s')
        sim_results.append((cfg, (snapshots, runtime)))

    print('\nGenerating HTML report...')
    generate_html(sim_results, output_path)

    # Open in Safari
    print('Opening report in Safari...')
    subprocess.run(['open', '-a', 'Safari', output_path])


if __name__ == '__main__':
    run_demo()
