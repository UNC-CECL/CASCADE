#!/usr/bin/env python3
"""
A plain-language pamphlet that walks through BRIE's alongshore step on six sections, one step at a time.

    python brie_sand_pamphlet.py   ->  brie_sand_pamphlet.html

The same six-section headland as 2-the-code (x = 1000, 1000, 900, 900, 900, 1000 m; 500 m sections). It follows
BRIE's three lines (angle, rate, update) for the first half of one year and then checks the books at every
boundary. Static HTML with inline drawings; no script runs in the page. The three diffusivities are BRIE's
coast_diff values under option A waves at 79, 90 and 101 degrees (see 6-pamphlet in the README).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# --- CONFIG ------------------------------------------------------------------
OUT = HERE / "brie_sand_pamphlet.html"
X = [1000, 1000, 900, 900, 900, 1000]      # shoreline position, m (BRIE counts landward: smaller = further out)
DY = 500.0                                  # section length, m
DT = 1.0                                    # years
COAST_DIFF = {79: 56934.823, 90: -119.294, 101: 9214.415}   # brie._coast_diff under option A
DEPTH = 19.8                                # shoreface depth, m: 1 m of shoreline on one section = DY*DEPTH m3
FULL_YEAR_LOSS = 15.0                       # Brie.update() on this coast, both halves of the year, m

# -----------------------------------------------------------------------------
N = len(X)
THETA = [180 * math.atan2(X[(i + 1) % N] - X[i], DY) / math.pi for i in range(N)]
RATE = [max(0.0, COAST_DIFF[round(90 - t)] * DT / 2 / DY ** 2) for t in THETA]
NORTH = [RATE[i] * (X[(i + 1) % N] - X[i]) for i in range(N)]   # change in x booked with the north neighbour
SOUTH = [RATE[i] * (X[i - 1] - X[i]) for i in range(N)]         # change in x booked with the south neighbour
GAIN = [-(NORTH[i] + SOUTH[i]) for i in range(N)]               # beach gained, m (minus: x counts landward)
TOTAL = sum(GAIN)


def f(v, nd=2, sign=True):
    """Number with a real minus sign."""
    s = f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"
    return s.replace("-", "−")


# --- drawings ----------------------------------------------------------------
L, SW, BASE = 40, 93, 214


def top(x):
    return BASE - 54 - (1000 - x) * 0.9


def cx(i):
    return L + SW * i + SW / 2


def coast_svg(label, extra="", show_x=True):
    parts = [f'<svg viewBox="0 0 640 250" role="img" aria-label="{label}">',
             '<rect class="sea" x="0" y="0" width="640" height="250"/>',
             '<text class="tag" x="14" y="24">SEA</text>',
             '<text class="tag" x="626" y="24" text-anchor="end">NORTH &#8594;</text>',
             '<text class="tag" x="14" y="44">&#8592; south</text>']
    for i, x in enumerate(X):
        y = top(x)
        parts.append(f'<rect class="sand" x="{L + SW * i}" y="{y:.0f}" width="{SW}" height="{BASE - y:.0f}"/>')
        if show_x:
            parts.append(f'<text class="xval" x="{cx(i):.0f}" y="{y + 20:.0f}" text-anchor="middle">{x} m</text>')
        parts.append(f'<text class="sec" x="{cx(i):.0f}" y="{BASE + 26}" text-anchor="middle">{i}</text>')
    for i in range(1, N):
        parts.append(f'<line class="edge" x1="{L + SW * i}" y1="{min(top(X[i - 1]), top(X[i])):.0f}" '
                     f'x2="{L + SW * i}" y2="{BASE}"/>')
    parts.append(f'<text class="tag" x="{L - 8}" y="{BASE + 26}" text-anchor="end">section</text>')
    parts.append(extra)
    parts.append('</svg>')
    return "".join(parts)


def angle_extra():
    out = []
    for i in range(N - 1):
        x1, y1, x2, y2 = cx(i), top(X[i]), cx(i + 1), top(X[i + 1])
        hot = abs(THETA[i]) > 0.01
        cls = "look hot" if hot else "look"
        out.append(f'<line class="{cls}" x1="{x1:.0f}" y1="{y1:.0f}" x2="{x2:.0f}" y2="{y2:.0f}" '
                   f'marker-end="url(#{"ah" if hot else "a"})"/>')
        out.append(f'<circle class="dot" cx="{x1:.0f}" cy="{y1:.0f}" r="4"/>')
        ty = min(y1, y2) - 14
        out.append(f'<text class="ang{" hotT" if hot else ""}" x="{(x1 + x2) / 2:.0f}" y="{ty:.0f}" '
                   f'text-anchor="middle">{f(THETA[i], 1, sign=abs(THETA[i]) > 0.01)}&#176;</text>')
    out.append(f'<circle class="dot" cx="{cx(N - 1):.0f}" cy="{top(X[-1]):.0f}" r="4"/>')
    out.append(f'<text class="ang" x="{cx(N - 1) + 10:.0f}" y="{top(X[-1]) - 14:.0f}" text-anchor="middle">'
               f'0&#176;*</text>')
    defs = ('<defs><marker id="a" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
            'orient="auto"><path class="mk" d="M0,0 L10,5 L0,10 z"/></marker>'
            '<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" '
            'orient="auto"><path class="mkh" d="M0,0 L10,5 L0,10 z"/></marker></defs>')
    return defs + "".join(out)


def move_extra():
    out = []
    for i in range(N):
        y = top(X[i])
        if RATE[i] == 0:
            out.append(f'<text class="zero" x="{cx(i):.0f}" y="{y - 12:.0f}" text-anchor="middle">r = 0</text>')
            continue
        gain = GAIN[i] > 0
        cls = "made" if gain else "lost"
        y0, y1 = (y - 6, y - 40) if gain else (y - 40, y - 6)
        out.append(f'<line class="mv {cls}" x1="{cx(i):.0f}" y1="{y0:.0f}" x2="{cx(i):.0f}" y2="{y1:.0f}" '
                   f'marker-end="url(#m{cls})"/>')
        out.append(f'<text class="mvT {cls}T" x="{cx(i):.0f}" y="{y - 50:.0f}" text-anchor="middle">'
                   f'{f(GAIN[i])} m</text>')
    defs = ''.join(f'<marker id="m{c}" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="5" markerHeight="5" '
                   f'orient="auto"><path class="mk{c}" d="M0,0 L10,5 L0,10 z"/></marker>' for c in ("made", "lost"))
    return f'<defs>{defs}</defs>' + "".join(out)


# --- tables ------------------------------------------------------------------
def rate_rows():
    rows = []
    for i in range(N):
        k = round(90 - THETA[i])
        r = RATE[i]
        chip = '<span class="chip on">has a rate</span>' if r > 0 else '<span class="chip">no rate</span>'
        rows.append(f'<tr><td>{i}</td><td>{f(THETA[i], 1, sign=abs(THETA[i]) > 0.01)}&#176;</td>'
                    f'<td>{f(COAST_DIFF[k], 0, sign=False)}</td><td>{r:.4f}</td><td>{chip}</td></tr>')
    return "".join(rows)


def move_rows():
    rows = []
    for i in range(N):
        n_gain, s_gain = -NORTH[i], -SOUTH[i]
        cls = "made" if GAIN[i] > 1e-9 else ("lost" if GAIN[i] < -1e-9 else "")
        rows.append(f'<tr class="{cls}"><td>{i}</td><td>{RATE[i]:.4f}</td>'
                    f'<td>{RATE[i]:.4f} &#215; {f(-(X[(i + 1) % N] - X[i]), 0, sign=False)} = {f(n_gain + 0.0)}</td>'
                    f'<td>{RATE[i]:.4f} &#215; {f(-(X[i - 1] - X[i]), 0, sign=False)} = {f(s_gain + 0.0)}</td>'
                    f'<td><b>{f(GAIN[i] + 0.0)}</b></td></tr>')
    return "".join(rows)


def ledger(shared):
    """One card per boundary: what the south section and the north section each write down."""
    cards = []
    for b in range(N - 1):
        s, n = b, b + 1
        south_says = -NORTH[s]                       # sand the south section gains across this boundary
        r_b = RATE[s]                                # the boundary's own rate (its angle is the south section's)
        north_says = -(r_b if shared else RATE[n]) * (X[s] - X[n])
        left = south_says + north_says               # should be zero
        if abs(left) < 1e-9:
            verdict, cls = "balanced", "ok"
        elif left > 0:
            verdict, cls = f"{f(left)} m of sand appears", "made"
        else:
            verdict, cls = f"{f(-left, sign=False)} m of sand vanishes", "lost"
        rate_note = (f"both use r = {r_b:.4f}" if shared else
                     f"{s} uses r = {RATE[s]:.4f}, {n} uses r = {RATE[n]:.4f}")
        cards.append(
            f'<div class="card {cls}"><div class="bname">boundary {s}&#8202;|&#8202;{n}</div>'
            f'<div class="says"><span>Section {s} writes</span><b>{f(south_says + 0.0)} m</b></div>'
            f'<div class="says"><span>Section {n} writes</span><b>{f(north_says + 0.0)} m</b></div>'
            f'<div class="rn">{rate_note}</div><div class="verdict">{verdict}</div></div>')
    return "".join(cards)


# --- page --------------------------------------------------------------------
CSS = """
/* Layout: a single-column folded leaflet; each numbered panel is one step, a drawing above plain text. */
:root{
  --paper:#F3F6F6; --panel:#FFFFFF; --ink:#16303A; --muted:#5A6E75; --rule:#C9D5D6;
  --sea:#DCEBEE; --sand:#E2C487; --sand-edge:#B8964F; --accent:#1F6F86;
  --made:#2F7D45; --made-wash:#E2F1E5; --lost:#B23A25; --lost-wash:#F8E3DE; --code-bg:#EAF0F1;
  --display:"Familjen Grotesk", "Helvetica Neue", Arial, sans-serif;
  --body:"Source Serif 4", Georgia, "Times New Roman", serif;
  --mono:"JetBrains Mono", Consolas, "Courier New", monospace;
}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
  --paper:#0F181B; --panel:#152226; --ink:#E2ECEC; --muted:#9AB0B5; --rule:#2B3C41;
  --sea:#17333B; --sand:#9C8150; --sand-edge:#C9A766; --accent:#6BC0D4;
  --made:#6CC387; --made-wash:#183326; --lost:#EE8270; --lost-wash:#3A1D18; --code-bg:#1B2B30; color-scheme:dark}}
:root[data-theme="dark"]{
  --paper:#0F181B; --panel:#152226; --ink:#E2ECEC; --muted:#9AB0B5; --rule:#2B3C41;
  --sea:#17333B; --sand:#9C8150; --sand-edge:#C9A766; --accent:#6BC0D4;
  --made:#6CC387; --made-wash:#183326; --lost:#EE8270; --lost-wash:#3A1D18; --code-bg:#1B2B30; color-scheme:dark}
body{background:var(--paper);color:var(--ink);font:17px/1.6 var(--body);padding-inline:16px;padding-block:28px 64px}
main{max-width:780px;margin:0 auto;display:flex;flex-direction:column;gap:0}
h1,h2,h3{font-family:var(--display);text-wrap:balance;line-height:1.15;margin:0}
h1{font-size:clamp(2rem,6vw,3.1rem);letter-spacing:-.01em}
h2{font-size:1.5rem}
p{margin:0;max-width:65ch}
.cover{padding-block:8px 28px;display:flex;flex-direction:column;gap:14px}
.eyebrow{font:600 .78rem var(--display);letter-spacing:.12em;text-transform:uppercase;color:var(--accent)}
.lede{font-size:1.15rem;color:var(--muted)}
.quote{border-left:3px solid var(--sand-edge);padding:4px 0 4px 16px;font-style:italic}
.panel{background:var(--panel);border:1px solid var(--rule);border-radius:6px;padding:26px 26px 28px;
  display:flex;flex-direction:column;gap:14px;position:relative}
.fold{height:26px;border-left:2px dashed var(--rule);margin-left:40px}
.step{display:flex;align-items:baseline;gap:12px}
.num{font:700 .8rem var(--display);letter-spacing:.1em;color:var(--accent);text-transform:uppercase;white-space:nowrap}
figure{margin:4px 0 0;display:flex;flex-direction:column;gap:6px}
figure svg{width:100%;height:auto;display:block;border-radius:4px}
figcaption{font-size:.85rem;color:var(--muted)}
.sea{fill:var(--sea)} .sand{fill:var(--sand);stroke:var(--sand-edge);stroke-width:1}
.edge{stroke:var(--sand-edge);stroke-width:1;stroke-dasharray:3 3}
.tag{font:600 11px var(--display);letter-spacing:.1em;fill:var(--muted)}
.xval{font:600 13px var(--mono);fill:var(--ink)}
.sec{font:700 15px var(--display);fill:var(--ink)}
.look{stroke:var(--muted);stroke-width:1.6;fill:none} .look.hot{stroke:var(--accent);stroke-width:2.6}
.dot{fill:var(--ink)} .mk{fill:var(--muted)} .mkh{fill:var(--accent)}
.ang{font:600 13px var(--mono);fill:var(--muted)} .hotT{fill:var(--accent);font-weight:700}
.zero{font:500 12px var(--mono);fill:var(--muted)}
.mv{stroke-width:3} .made{stroke:var(--made)} .lost{stroke:var(--lost)}
.mkmade{fill:var(--made)} .mklost{fill:var(--lost)}
.mvT{font:700 14px var(--mono)} .madeT{fill:var(--made)} .lostT{fill:var(--lost)}
code,.code{font-family:var(--mono);font-size:.86em}
.code{display:block;background:var(--code-bg);border-radius:4px;padding:12px 14px;overflow-x:auto;white-space:pre;line-height:1.55}
.code .c{color:var(--muted)}
.tw{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:.92rem;font-variant-numeric:tabular-nums}
th{font:600 .74rem var(--display);letter-spacing:.08em;text-transform:uppercase;color:var(--muted);text-align:left;
  border-bottom:1px solid var(--rule);padding:6px 10px 6px 0}
td{padding:6px 10px 6px 0;border-bottom:1px solid var(--rule);font-family:var(--mono);font-size:.86rem;white-space:nowrap}
tr.made td{background:var(--made-wash)} tr.lost td{background:var(--lost-wash)}
.chip{font:600 .72rem var(--display);letter-spacing:.04em;padding:2px 8px;border-radius:99px;border:1px solid var(--rule);color:var(--muted)}
.chip.on{border-color:var(--accent);color:var(--accent)}
.note{font-size:.92rem;color:var(--muted)}
.key{background:var(--code-bg);border-radius:4px;padding:12px 16px}
.ledger{display:grid;grid-template-columns:repeat(auto-fit,minmax(138px,1fr));gap:10px}
.card{border:1px solid var(--rule);border-radius:5px;padding:10px 12px;display:flex;flex-direction:column;gap:4px;min-width:0}
.card.made{background:var(--made-wash);border-color:var(--made)} .card.lost{background:var(--lost-wash);border-color:var(--lost)}
.bname{font:700 .8rem var(--display);letter-spacing:.06em;text-transform:uppercase}
.says{display:flex;justify-content:space-between;gap:8px;font-size:.84rem}
.says b{font-family:var(--mono);font-size:.84rem}
.rn{font:.72rem var(--mono);color:var(--muted)}
.verdict{font:700 .86rem var(--display)} .card.made .verdict{color:var(--made)} .card.lost .verdict{color:var(--lost)}
.card.ok .verdict{color:var(--muted)}
.total{display:flex;flex-wrap:wrap;gap:6px 18px;align-items:baseline;font-family:var(--display)}
.total b{font-size:1.6rem}
.total.lost b{color:var(--lost)} .total.ok b{color:var(--made)}
.two{display:grid;grid-template-columns:1fr 1fr;gap:14px}
.two>div{min-width:0;display:flex;flex-direction:column;gap:6px}
.two h3{font-size:1.02rem}
@media (max-width:620px){.two{grid-template-columns:1fr}.panel{padding:20px 16px 22px}}
.foot{font-size:.85rem;color:var(--muted);padding-top:22px;display:flex;flex-direction:column;gap:6px}
"""


def page():
    gain1, loss4 = GAIN[1], -GAIN[4]
    lost = -TOTAL
    return f"""<title>BRIE's Missing Sand</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Familjen+Grotesk:wght@500;600;700&family=JetBrains+Mono:wght@400;500;700&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<style>{CSS}</style>
<main>

<header class="cover">
  <div class="eyebrow">CASCADE &#183; BRIE alongshore step &#183; a walk-through</div>
  <h1>Where does BRIE's sand go at a sharp bend?</h1>
  <p class="quote">"BRIE looks at both sides of each section, so no sand can be made or lost."</p>
  <p class="lede">That is half right. BRIE uses the <b>positions</b> on both sides, but each section's
  <b>rate</b> comes from its north side only. Where a sharp bend gives two neighbours different rates,
  they disagree about how much sand crossed between them, and the difference is sand made or lost.
  This pamphlet follows one year on six sections of coast to show exactly where.</p>
</header>

<section class="panel">
  <div class="step"><span class="num">Step 1</span><h2>Six sections of coast</h2></div>
  <figure>{coast_svg("Six sections, the middle three sticking out 100 m")}
    <figcaption>Each block is one 500 m section, looking down from above. Taller blocks stick further out to sea.
    Height is exaggerated.</figcaption></figure>
  <p>BRIE splits the coast into sections and keeps one number for each: the shoreline position, <code>x</code>.
  Here sections 2, 3 and 4 form a headland that sticks out 100 m.</p>
  <p class="note">BRIE measures <code>x</code> landward from a fixed line, so the headland has the
  <em>smaller</em> number (900 m). Sections are numbered from 0, as in the code.</p>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">Step 2</span><h2>Each section looks north to find its angle</h2></div>
  <figure>{coast_svg("Each section's angle, measured to its north neighbour", angle_extra())}
    <figcaption>Arrows run from each section to its north neighbour. The label is that section's angle.
    *Section 5's north neighbour is section 0, because BRIE joins the two ends.</figcaption></figure>
  <span class="code"><span class="c"># brie.py 820-828: x[i+1] - x[i], the north neighbour only</span>
theta = 180 * arctan2(x_s[np.r_[1:ny, 0]] - x_s, dy) / pi</span>
  <p>The angle is worked out from one pair only: the section and the one to its north. The south neighbour
  plays no part.</p>
  <p>Only two sections face a slope to the north: section 1 (the coast steps <em>out</em> to the headland,
  {f(THETA[1], 1)}&#176;) and section 4 (the coast steps <em>back</em>, {f(THETA[4], 1)}&#176;). Every other
  section sees flat coast to its north, so its angle is 0&#176;.</p>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">Step 3</span><h2>The angle sets the section's rate</h2></div>
  <span class="code"><span class="c"># brie.py 1294-1309: one rate per section, from its angle</span>
r = max(0, coast_diff[round(90 - theta)] * dt / 2 / dy**2)</span>
  <div class="tw"><table>
    <thead><tr><th>Section</th><th>Angle</th><th>BRIE's table</th><th>Rate r</th><th></th></tr></thead>
    <tbody>{rate_rows()}</tbody></table></div>
  <p>BRIE looks the angle up in a table built from the wave climate and turns it into a rate: how fast this
  section passes sand to its neighbours. Only sections 1 and 4 get one.</p>
  <p>Section 2 is on the headland but gets no rate, because its north side is flat. Section 5 sits right
  next to the headland's edge but also gets no rate, for the same reason.</p>
  <p class="note">Why is a flat side's rate zero? Under our waves (option A, mostly from the north) BRIE's table
  gives {f(COAST_DIFF[90], 0)} at 0&#176;, and the code floors anything negative at zero.</p>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">Step 4</span><h2>Each section moves using both neighbours, but its own rate</h2></div>
  <span class="code"><span class="c"># brie.py 1315-1324</span>
x_new[i] = x[i] + r[i] * (x[i+1] - 2*x[i] + x[i-1])
<span class="c">#        = x[i] + r[i]*(x[i+1] - x[i])  +  r[i]*(x[i-1] - x[i])</span>
<span class="c">#                 with the north side      with the south side</span></span>
  <p>This is the "both sides" line. The positions to the north and to the south both appear.
  But the same <code>r[i]</code>, the one from the north side, multiplies both.</p>
  <figure>{coast_svg("How much beach each section gains or loses", move_extra(), show_x=False)}
    <figcaption>Beach gained (green) or lost (red) by each section in the first half of the year.</figcaption></figure>
  <div class="tw"><table>
    <thead><tr><th>Section</th><th>r</th><th>With north</th><th>With south</th><th>Beach gained</th></tr></thead>
    <tbody>{move_rows()}</tbody></table></div>
  <p>Section 1 gains {f(gain1, sign=False)} m of beach from the headland. Section 4, the headland's north end,
  loses {f(loss4, sign=False)} m. Every other section has r = 0, so it does not move at all, even sections
  2 and 5, which sit right next to a change.</p>
  <p class="note">The table shows beach gained, so the signs are flipped from BRIE's <code>x</code>,
  which counts landward.</p>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">Step 5</span><h2>Check the books at each boundary</h2></div>
  <p>Sand that crosses a boundary should leave one section and arrive in the other. So the two sections
  either side should write down equal and opposite amounts, and their entries should add to zero.</p>
  <div class="ledger">{ledger(shared=False)}</div>
  <p>At boundary 1|2, section 1 writes down {f(gain1, sign=False)} m received. Section 2, with r = 0, writes
  down nothing sent. That sand appears from nowhere.</p>
  <p>At boundary 4|5, section 4 writes down {f(loss4, sign=False)} m sent north. Section 5, with r = 0,
  writes down nothing received. That sand vanishes.</p>
  <div class="total lost"><span>Island total after half a year:</span><b>{f(TOTAL)} m of beach</b>
    <span class="note">about {round(lost * DY * DEPTH, -3):,.0f} m&#179; of sand</span></div>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">Step 6</span><h2>Give each boundary one rate, and the books balance</h2></div>
  <p>A boundary's angle is the angle of the section to its south (step 2), so its rate belongs to that boundary.
  If the section on the far side used the same rate for it, both would book the same amount.</p>
  <div class="ledger">{ledger(shared=True)}</div>
  <div class="total ok"><span>Island total:</span><b>0.00 m</b>
    <span class="note">The headland now loses sand at both ends and sections 2 and 5 take up what crosses.</span></div>
  <p>Section 2 now gives up {f(gain1, sign=False)} m to section 1, and section 5 takes in the
  {f(loss4, sign=False)} m from section 4. Nothing is left over.</p>
</section>
<div class="fold" aria-hidden="true"></div>

<section class="panel">
  <div class="step"><span class="num">So</span><h2>Who is right?</h2></div>
  <div class="two">
    <div><h3>What the committee member is right about</h3>
      <p>BRIE's equation conserves sand, and its step uses both neighbours' positions. On a straight coast it
      loses nothing. On a gently curving coast it loses very little.</p></div>
    <div><h3>What happens at a sharp bend</h3>
      <p>Each section takes its rate from its north side only and uses it on both sides. At a sharp bend, two
      neighbours end up with different rates for the boundary they share, and sand is made or lost there.</p></div>
  </div>
  <p>The worked example above is the first half of one year. BRIE also takes a second half with the same rates.
  The whole year on this coast, from the real <code>Brie.update()</code>, loses {FULL_YEAR_LOSS:.1f} m of beach,
  about {round(FULL_YEAR_LOSS * DY * DEPTH, -3):,.0f} m&#179;. Over 50 years on a 30-section island:</p>
  <div class="tw"><table>
    <thead><tr><th>Coast at the start</th><th>Beach lost in 50 years</th><th>One rate per boundary</th></tr></thead>
    <tbody>
      <tr><td>straight</td><td>0.0 m</td><td>0</td></tr>
      <tr><td>gentle 100 m bump</td><td>10.6 m</td><td>0</td></tr>
      <tr class="lost"><td>sharp 100 m step, 3 sections wide</td><td>117.9 m (about 1.2 million m&#179;)</td><td>0</td></tr>
    </tbody></table></div>
  <p class="key"><b>Why it matters for the groin.</b> A groin makes exactly this kind of sharp bend. The
  fixed-shift groin's results were mostly this error. The wave-driven groin works out the disagreement at the
  boundaries near the groin each year and puts it back, so its total stays within a few metres. BRIE itself is
  not edited.</p>
</section>

<footer class="foot">
  <span>Numbers from unmodified BRIE with the barrier and inlet models off, so a year moves the shoreline only by
  its alongshore step. Option A waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle fraction 0.5), 500 m sections,
  shoreface depth {DEPTH} m.</span>
  <span>Made by <code>hard-structures/groin/2-module-explained/3-brie-sand-balance/6-pamphlet/brie_sand_pamphlet.py</code>.
  Figures for each step are in folders 1 to 5 beside it.</span>
</footer>
</main>
"""


if __name__ == "__main__":
    OUT.write_text(page(), encoding="utf-8")
    print(f"wrote {OUT}  (half-year total {TOTAL:+.3f} m; r1 {RATE[1]:.5f}, r4 {RATE[4]:.5f})")
