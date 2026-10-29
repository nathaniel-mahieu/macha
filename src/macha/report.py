"""Self-contained HTML explorers (plotly figures + tables) with light/dark themes.

Figures are built in Python with the light palette; the page swaps every
light token for its selected dark step when the viewer's theme is dark, so
both modes are designed rather than inverted.
"""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

PLOTLY_JS = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@4.1.1/plotly.min.js"

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SERIES_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]
STATUS = {"good": "#0ca30c", "warning": "#fab219", "serious": "#ec835a", "critical": "#d03b3b"}
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
LIGHT = dict(surface="#fcfcfb", page="#f9f9f7", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9",
             axis="#c3c2b7", wash="#f0efec")
DARK = dict(surface="#1a1a19", page="#0d0d0d", ink="#f0efec", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a",
            axis="#383835", wash="#383835")
FONT = 'system-ui, -apple-system, "Segoe UI", sans-serif'

# light -> dark token map used by the page script
DARK_MAP = {**dict(zip(SERIES, SERIES_DARK)), **{LIGHT[k]: DARK[k] for k in LIGHT if LIGHT[k] != DARK[k]}}

METHOD_COLORS = {}  # stable method -> slot assignment so colour follows the entity across pages


def method_color(name: str, order=("macha", "greedy", "slice", "rect", "density", "global", "dtw", "none",
                                   "sspp", "cwt", "als_zero", "als_locf", "rolling_min", "threshold",
                                   "threshold+corr")) -> str:
    if name not in METHOD_COLORS:
        idx = order.index(name) if name in order else len(METHOD_COLORS)
        METHOD_COLORS[name] = SERIES[idx % len(SERIES)]
    return METHOD_COLORS[name]


def style(fig, height: int = 380, legend: bool = True, title: str | None = None):
    """Apply the house chart style: hairline solid grids, ink-token text, quiet chrome."""
    ax = dict(gridcolor=LIGHT["grid"], linecolor=LIGHT["axis"], zerolinecolor=LIGHT["axis"], gridwidth=1,
              linewidth=1, ticks="outside", tickcolor=LIGHT["axis"], ticklen=4,
              tickfont=dict(color=LIGHT["muted"], size=11), title=dict(font=dict(color=LIGHT["ink2"], size=12)))
    fig.update_layout(
        template="none", height=height, paper_bgcolor=LIGHT["surface"], plot_bgcolor=LIGHT["surface"],
        font=dict(family=FONT, color=LIGHT["ink"], size=12), margin=dict(l=64, r=20, t=48 if title else 24, b=52),
        showlegend=legend, legend=dict(orientation="h", y=1.02, yanchor="bottom", x=0, font=dict(color=LIGHT["ink2"]),
                                       bgcolor="rgba(0,0,0,0)"),
        hoverlabel=dict(bgcolor=LIGHT["surface"], bordercolor=LIGHT["axis"], font=dict(color=LIGHT["ink"],
                                                                                         family=FONT)),
        title=dict(text=title, font=dict(size=14, color=LIGHT["ink"]), x=0, xanchor="left") if title else None,
    )
    fig.update_xaxes(**ax)
    fig.update_yaxes(**ax)
    for ann in fig.layout.annotations or []:
        ann.font.color = LIGHT["ink2"]
    return fig


def _trace_color(tr):
    """(colour, swatch kind) of a trace for the HTML legend."""
    mode = getattr(tr, "mode", None) or ""
    for attr, kind in (("line", "line"), ("marker", "dot")):
        obj = getattr(tr, attr, None)
        c = getattr(obj, "color", None) if obj is not None else None
        if isinstance(c, str):
            if attr == "line" and ("lines" in mode or tr.type in ("scatter",) and mode == ""):
                return c, "line"
            if attr == "marker":
                return c, "bar" if tr.type in ("bar", "histogram") else "dot"
    obj = getattr(tr, "line", None)
    c = getattr(obj, "color", None) if obj is not None else None
    return (c, "line") if isinstance(c, str) else (None, None)


def _fmt(v, digits=3):
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return "–"
        if abs(v) >= 1e4 or (abs(v) < 1e-3 and v != 0):
            return f"{v:.2e}"
        return f"{v:.{digits}f}"
    if isinstance(v, (int, np.integer)):
        return f"{v:,}"
    return html.escape(str(v))


class Page:
    def __init__(self, title: str, subtitle: str = "", nav: list | None = None):
        self.title = title
        self.subtitle = subtitle
        self.blocks: list[str] = []
        self.figs: list[tuple[str, str]] = []
        self.nav = nav or []

    def section(self, title: str, text: str = ""):
        slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
        self.blocks.append(f'<section id="{slug}"><h2>{html.escape(title)}</h2>{text}</section>')

    def text(self, html_text: str):
        self.blocks.append(f'<div class="prose">{html_text}</div>')

    def tiles(self, items):
        """items: list of (label, value, note)."""
        cells = "".join(f'<div class="tile"><div class="tlabel">{html.escape(a)}</div><div class="tval">{b}</div>'
                        f'<div class="tnote">{html.escape(c)}</div></div>' for a, b, c in items)
        self.blocks.append(f'<div class="tiles">{cells}</div>')

    def figure(self, fig, caption: str = ""):
        """Add a figure; its title and legend are rendered as HTML above the plot so they wrap at any width."""
        fid = f"fig{len(self.figs)}"
        title = fig.layout.title.text if fig.layout.title and fig.layout.title.text else ""
        chips, seen = [], set()
        if fig.layout.showlegend is not False:
            for tr in fig.data:
                name = tr.name
                if not name or tr.showlegend is False or name in seen:
                    continue
                seen.add(name)
                color, kind = _trace_color(tr)
                if color is None:
                    continue
                chips.append(f'<button type="button" class="chip" data-fig="{fid}" data-name="{html.escape(name)}">'
                             f'<span class="sw {kind}" data-color="{color}" style="background:{color}"></span>'
                             f'{html.escape(name)}</button>')
        has_menu = bool(fig.layout.updatemenus)
        fig.update_layout(title=None, showlegend=False, margin=dict(t=44 if has_menu else 12))
        self.figs.append((fid, fig.to_json()))
        head = ""
        if title or chips:
            legend = f'<div class="legend">{"".join(chips)}</div>' if len(chips) >= 2 else ""
            head = f'<div class="fighead"><div class="ftitle">{title}</div>{legend}</div>'
        cap = f"<figcaption>{caption}</figcaption>" if caption else ""
        self.blocks.append(f'<figure>{head}<div class="plot" id="{fid}"></div>{cap}</figure>')

    def table(self, df: pd.DataFrame, caption: str = "", better: dict | None = None, digits: int = 3,
              group_col: str | None = None):
        """``better``: {column: "max" | "min"}; the best value per group (or overall) is set in bold."""
        better = better or {}
        best = {}
        groups = df.groupby(group_col).groups if group_col else {None: df.index}
        for g, idx in groups.items():
            for c, how in better.items():
                if c in df:
                    vals = pd.to_numeric(df.loc[idx, c], errors="coerce")
                    if vals.notna().any():
                        best.setdefault(c, set()).add(vals.idxmax() if how == "max" else vals.idxmin())
        head = "".join(f"<th>{html.escape(str(c))}{' ↑' if better.get(c) == 'max' else ' ↓' if better.get(c) == 'min' else ''}</th>"
                       for c in df.columns)
        rows = []
        for i, r in df.iterrows():
            tds = []
            for c in df.columns:
                classes = []
                if isinstance(r[c], (int, float, np.integer, np.floating)):
                    classes.append("num")
                if i in best.get(c, ()):
                    classes.append("best")
                tds.append(f'<td class="{" ".join(classes)}">{_fmt(r[c], digits)}</td>')
            rows.append("<tr>" + "".join(tds) + "</tr>")
        cap = f"<caption>{caption}</caption>" if caption else ""
        self.blocks.append(f'<div class="tablewrap"><table>{cap}<thead><tr>{head}</tr></thead>'
                           f'<tbody>{"".join(rows)}</tbody></table></div>')

    def write(self, path: str | Path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        nav = " · ".join(f'<a href="{html.escape(h)}">{html.escape(t)}</a>' for t, h in self.nav)
        figs = ",".join(f'["{fid}", {spec}]' for fid, spec in self.figs)
        doc = TEMPLATE.format(title=html.escape(self.title), subtitle=self.subtitle, nav=nav,
                              body="\n".join(self.blocks), figs=figs, plotly=PLOTLY_JS,
                              darkmap=json.dumps(DARK_MAP), font=FONT, **{f"L_{k}": v for k, v in LIGHT.items()},
                              **{f"D_{k}": v for k, v in DARK.items()})
        path.write_text(doc)
        return path


TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script src="{plotly}"></script>
<style>
:root {{ color-scheme: light; --surface:{L_surface}; --page:{L_page}; --ink:{L_ink}; --ink2:{L_ink2}; --muted:{L_muted};
  --grid:{L_grid}; --axis:{L_axis}; --wash:{L_wash}; --ring: rgba(11,11,11,0.10); --accent:#2a78d6; }}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{ color-scheme: dark; --surface:{D_surface};
  --page:{D_page}; --ink:{D_ink}; --ink2:{D_ink2}; --muted:{D_muted}; --grid:{D_grid}; --axis:{D_axis}; --wash:{D_wash};
  --ring: rgba(255,255,255,0.10); --accent:#3987e5; }} }}
:root[data-theme="dark"] {{ color-scheme: dark; --surface:{D_surface}; --page:{D_page}; --ink:{D_ink}; --ink2:{D_ink2};
  --muted:{D_muted}; --grid:{D_grid}; --axis:{D_axis}; --wash:{D_wash}; --ring: rgba(255,255,255,0.10); --accent:#3987e5; }}
* {{ box-sizing: border-box; }}
body {{ margin:0; background:var(--page); color:var(--ink); font-family:{font}; font-size:15px; line-height:1.55; }}
main {{ max-width:1180px; margin:0 auto; padding:28px 16px 64px; }}
header {{ display:flex; justify-content:space-between; gap:16px; align-items:flex-start; flex-wrap:wrap; }}
h1 {{ font-size:26px; margin:0 0 4px; letter-spacing:-0.01em; }}
.sub {{ color:var(--ink2); margin:0; max-width:820px; }}
nav {{ font-size:13px; color:var(--muted); margin:10px 0 0; }}
nav a, .prose a {{ color:var(--accent); text-decoration:none; }}
button#theme {{ background:var(--surface); color:var(--ink2); border:1px solid var(--ring); border-radius:8px;
  padding:6px 12px; font:inherit; font-size:13px; cursor:pointer; }}
section h2 {{ font-size:19px; margin:40px 0 6px; padding-top:14px; border-top:1px solid var(--grid); }}
section p, .prose p, .prose li {{ color:var(--ink2); max-width:860px; }}
.prose code, section code {{ background:var(--wash); padding:1px 5px; border-radius:4px; font-size:13px; }}
figure {{ margin:14px 0 22px; background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:8px; }}
figcaption {{ color:var(--muted); font-size:13px; padding:4px 8px 4px; }}
.fighead {{ padding:6px 8px 2px; }}
.ftitle {{ font-weight:600; font-size:14px; color:var(--ink); }}
.legend {{ display:flex; flex-wrap:wrap; gap:4px 14px; margin-top:6px; }}
.chip {{ display:inline-flex; align-items:center; gap:6px; background:none; border:0; padding:2px 0; font:inherit;
  font-size:12px; color:var(--ink2); cursor:pointer; }}
.chip.off {{ opacity:0.4; }}
.sw {{ display:inline-block; flex:none; }}
.sw.dot {{ width:9px; height:9px; border-radius:50%; }}
.sw.line {{ width:16px; height:2px; border-radius:1px; }}
.sw.bar {{ width:10px; height:10px; border-radius:2px; }}
.plot {{ width:100%; min-height:200px; }}
.tiles {{ display:grid; grid-template-columns:repeat(auto-fit, minmax(170px, 1fr)); gap:10px; margin:16px 0; }}
.tile {{ background:var(--surface); border:1px solid var(--ring); border-radius:10px; padding:12px 14px; }}
.tlabel {{ font-size:12px; color:var(--muted); }}
.tval {{ font-size:26px; font-weight:600; margin:2px 0; }}
.tnote {{ font-size:12px; color:var(--ink2); }}
.tablewrap {{ overflow-x:auto; margin:12px 0 20px; background:var(--surface); border:1px solid var(--ring); border-radius:10px; }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
caption {{ text-align:left; color:var(--muted); padding:8px 12px; caption-side:top; }}
th, td {{ padding:6px 12px; border-bottom:1px solid var(--grid); text-align:left; white-space:nowrap; }}
th {{ color:var(--ink2); font-weight:600; }}
td.num, td.best.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
td.best {{ font-weight:700; }}
@media (max-width:640px) {{ h1 {{ font-size:22px; }} .tval {{ font-size:22px; }} }}
</style></head>
<body><main>
<header><div><h1>{title}</h1><p class="sub">{subtitle}</p><nav>{nav}</nav></div>
<button id="theme" type="button" aria-label="Toggle light or dark theme">Theme</button></header>
{body}
</main>
<script>
const FIGS = [{figs}];
const DARK = {darkmap};
function isDark() {{
  const t = document.documentElement.dataset.theme;
  if (t) return t === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
}}
function swap(x) {{
  if (typeof x === "string") return DARK[x.toLowerCase()] || x;
  if (Array.isArray(x)) return x.map(swap);
  if (x && typeof x === "object") {{ const o = {{}}; for (const k in x) o[k] = swap(x[k]); return o; }}
  return x;
}}
function render() {{
  const dark = isDark();
  for (const [id, spec] of FIGS) {{
    const el = document.getElementById(id);
    if (!el) continue;
    const s = dark ? swap(spec) : spec;
    Plotly.react(el, s.data, s.layout, {{responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d"]}});
  }}
  paintChips();
}}
function paintChips() {{
  const dark = isDark();
  for (const sw of document.querySelectorAll(".sw")) {{
    const c = sw.dataset.color; sw.style.background = dark ? (DARK[c.toLowerCase()] || c) : c;
  }}
}}
for (const chip of document.querySelectorAll(".chip")) {{
  chip.addEventListener("click", () => {{
    const el = document.getElementById(chip.dataset.fig);
    const name = chip.dataset.name;
    const off = !chip.classList.contains("off");
    chip.classList.toggle("off", off);
    const idx = el.data.map((t, i) => t.name === name ? i : -1).filter(i => i >= 0);
    // traces hidden by a menu stay hidden; the chip only toggles what the menu currently shows
    const vis = idx.map(i => off ? (el.data[i].visible === false ? false : "legendonly")
                               : (el.data[i].visible === "legendonly" ? true : el.data[i].visible));
    Plotly.restyle(el, {{visible: vis}}, idx);
  }});
}}
document.getElementById("theme").addEventListener("click", () => {{
  document.documentElement.dataset.theme = isDark() ? "light" : "dark"; render();
}});
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
render();
</script></body></html>
"""


def binned_hist(x, bins=50, range_=None, density=True):
    """Pre-binned histogram as plotly bar kwargs (keeps pages small; raw values are not embedded)."""

    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    h, e = np.histogram(x, bins=bins, range=range_, density=density)
    return dict(x=0.5 * (e[1:] + e[:-1]), y=h, width=np.diff(e) * 0.92)
