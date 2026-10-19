# 3-brie-sand-balance — does BRIE's alongshore step conserve sand?

**The question colleagues ask.** BRIE's shoreline equation conserves sand, so how can anything be "made"
or "lost"? They are right about the equation, and about most coasts. The issue is how the step is written
down on the grid, and it only shows where the coast bends sharply.

**What is right.** On a straight coast BRIE's step conserves sand exactly, and on a gently curving one
almost exactly (a 100 m bump spread over a few kilometres loses 10.6 m of summed shoreline in 50 years,
1.4%).

**What breaks, and where.** Each domain takes one rate, r, from the angle of the boundary to its **north**
neighbour (`brie.py` 820-826, 1294-1309) and uses that rate on **both** of its boundaries (1311-1324). Sand
crossing a boundary should leave one domain and arrive in the next in the same amount. At a bend, though,
the two domains either side of a boundary have different north-boundary angles, so they use different
rates for the boundary they share. Their books disagree, and the difference is sand made or lost. Give each
boundary one rate that both domains use, and every case below conserves to round-off.

**The evidence is unmodified BRIE.** Every run here is a real `brie.Brie` with the barrier and inlet models
off, so `update()` moves the shoreline only by its alongshore step. Option A waves (Hs 2.0 m, Tp 7.5 s,
asymmetry 0.6, high-angle 0.5), 30 domains of 500 m, nothing added or removed. A copy of BRIE's step
reproduces `Brie.update()` to three decimals.

| coast (start) | after 50 years | with one rate per boundary |
|---|---|---|
| straight | 0.000 m | 0 |
| gentle 100 m bump | −10.6 m | 0 |
| sharp 100 m step, 3 domains wide (300 m of sand) | **−117.9 m** (39% of it; about 1.2 million m³) | 0 |

Over 40 years on four shapes (`5-shape-matters/`): a groin-like fillet **gains** 8.8 m, a sharp headland loses
105.8 m, a gentle bump 9.0 m and a bay 51.3 m. The error grows roughly with the square of the turn in the
coast, so it matters at a groin's fillet or a sharp natural bend, and hardly at all elsewhere.

**Why it matters for the groin.** A groin makes exactly the kind of sharp bend where this happens. The
fixed-shift groin's results were mostly this error (+539 m of sand made in 50 years at −20°). The
wave-driven groin's leak correction works out the disagreement at the boundaries near the groin each year
and adds it back, so its total stays within a few metres (`../2-wave-driven/`). BRIE itself is not edited.

## Layout

One folder per question, in the order to explain it. Each holds its script, `figures/`, and `figures/supporting/`
(each figure's PDF and its Tests / How to read / Shows caption in `CAPTIONS.md`). `common.py` holds what they
share: the BRIE runs (alongshore step only), the one-rate-per-boundary comparison and the drawing helpers. GIFs
are git-ignored; rerun a folder's script to make them.

```
common.py                  BRIE runs, the shared-rate solve, drawing helpers (imported, not run)
1-the-idea/                the bookkeeping problem, with no model
    brie_sand_the_idea.png
2-the-code/                "BRIE looks at both sides"; BRIE's three lines worked through on six sections
    brie_sand_both_sides_one_rate.png         positions from both sides, the rate from the north boundary only
    brie_sand_equation_step_by_step.gif       the code line by line, with the coast and the working (~2 min)
    brie_sand_equation_worked_example.png     the same on one page of tables
3-one-year/                one year at a 100 m headland: who books what at each edge
    brie_sand_headland_one_year.png           plain: +3.7 m appears, 22.8 m vanishes, 15 m short by year's end
    brie_sand_one_boundary_books.png          the same with the rule, the code and one rate per boundary
4-over-time/               50 years
    brie_sand_headland_50_years.gif           the headland told step by step, then 50 years beside one rate per boundary
    brie_sand_total_over_time.png             straight, gentle bump, sharp step: the island's total sand over time
5-shape-matters/           how the coast's shape sets the size of the error
    brie_sand_where_it_matters.png            straight 0 m, gentle bump −11 m, sharp headland −118 m
    brie_sand_smooth_medium_sharp.gif         one 100 m feature three ways, 50 years (+ _end.png, the last frame)
    brie_sand_how_sharp.png                   one year against how sharply a step turns the coast
    brie_sand_four_shapes.png                 fillet, headland, bump, bay: each boundary's books and the totals
    brie_sand_four_shapes_step_by_step.gif    the four shapes year by year
6-pamphlet/                the whole story in plain words, one step per panel (the six-section example, first half-year)
    brie_sand_pamphlet.py                     writes the page; static HTML, no script in the page
    brie_sand_pamphlet.html                   open in a browser; also published as a private Claude artifact
6-worksheet/               work BRIE's step by hand: fill-in steps with checks, hints and a sandbox
    brie_sand_worksheet.html                  open in a browser; also published as a private Claude artifact
```

To present it, go 1 → 5 (the worksheet in 6 is for working it through yourself); within a folder, the first figure is the plain one and the rest add detail.

Run (each about a minute; `5-shape-matters` a few minutes, `--no-gif` for its figures only):

    python 1-the-idea/the_idea_figure.py
    python 2-the-code/the_code_figures.py
    python 3-one-year/one_year_figures.py
    python 4-over-time/over_time_figures.py
    python 5-shape-matters/shape_figures.py
