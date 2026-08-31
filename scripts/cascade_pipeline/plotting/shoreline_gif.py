"""Plan-view shoreline animation (GIF) for a completed CASCADE run.

One frame per model year; current shoreline drawn over a year-0 reference
(dashed grey), shaded blue where seaward of the reference, red where
landward. See make_shoreline_gif's docstring for the three y-axis modes.
"""

import dataclasses
import io
import os
import subprocess
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from matplotlib.transforms import blended_transform_factory

from cascade_pipeline.annotations import DEFAULT_ANNOTATIONS, add_geographic_annotations
from cascade_pipeline.domains import DEFAULT_DOMAINS


@dataclasses.dataclass(frozen=True)
class GifConfig:
    """Shared settings for every shoreline-animation job.

    Attributes:
        fps: Frames per second.
        year_stride: 1 = every model year, 2 = every other year, etc.
        annotate: Draw the geographic annotation layer (GIS-axis modes).
        auto_open: Pop the finished GIF open in the OS default viewer.
        keep_frames: Also write individual frame PNGs to disk.
        save_matrix: Save shoreline_m.npy alongside the GIFs in
            make_all_shoreline_gifs, so this run can serve as a future
            baseline for a "difference" job.
        ocean_at_bottom: True -> ocean/seaward at the bottom of the plot,
            sound/landward at the top (Hatteras' real cross-shore layout).
            False -> seaward at the top.
        baseline_label: Shown in a "difference" job's title/y-axis label.
        target_label: Legend label for the optional observed/target line
            drawn by a job that passes target_m.
    """

    fps: int = 3
    year_stride: int = 1
    annotate: bool = True
    auto_open: bool = False
    keep_frames: bool = False
    save_matrix: bool = True
    ocean_at_bottom: bool = True
    baseline_label: str = "baseline"
    target_label: str = "observed target"


DEFAULT_GIF_CONFIG = GifConfig()

# Shared with road_planview.RoadPlanViewStyle.relocated and the plan-view
# animation, so a relocation reads as the same event in all three.
RELOCATION_COLOR = "#FF8C00"


def _open_file(path):
    """Open a file in the OS default application (cross-platform, non-fatal)."""
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)  # type: ignore[attr-defined]  # noqa
        elif sys.platform == "darwin":
            subprocess.run(["open", path], check=False)
        else:
            subprocess.run(["xdg-open", path], check=False)
    except Exception as e:
        print(f"  [GIF] could not auto-open ({e}); open it manually: {path}")


def _slug(text):
    """Filename-safe short label (e.g. 'North Jetty' -> 'NorthJetty')."""
    return "".join(ch for ch in str(text) if ch.isalnum())


def _groin_zoom_window(gdom, pad, domains):
    """Inclusive GIS window centred on one groin, clipped to real domains.

    ANN_GROINS-style values sit on a domain boundary (e.g. 5.5 = the 5/6
    interface), so the window is built around that interface, not a cell
    centre.
    """
    lo = int(max(domains.first_gis_id, np.floor(gdom - pad)))
    hi = int(min(domains.last_gis_id, np.ceil(gdom + pad)))
    return lo, hi


def _select_groins(which, annotations):
    """Resolve a job's "which" key into [(name, domain), ...] from annotations.groins.

    Args:
        which: "all" (default, every groin -- new structures are picked up
            with no config change), a single name, or a list of names.
        annotations: AnnotationConfig.
    """
    if not annotations.groins:
        raise ValueError("range='groin' requires at least one entry in annotations.groins.")

    if isinstance(which, str) and which.strip().lower() == "all":
        return list(annotations.groins.items())

    names = [which] if isinstance(which, str) else list(which)
    missing = [n for n in names if n not in annotations.groins]
    if missing:
        raise ValueError(
            f"groin(s) {missing} not in annotations.groins (have: {list(annotations.groins)})."
        )
    return [(n, annotations.groins[n]) for n in names]


def _resolve_gif_domain_range(spec, pad, groin, domains, annotations):
    """Turn a GIF job's "range" value into concrete plotting coordinates.

    Args:
        spec: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi).
        pad: half-width in domains; only read for "groin"/"groin_span".
        groin: (name, domain) tuple, required when spec == "groin". Supplied
            by make_all_shoreline_gifs, which fans one "groin" job out into
            one call per structure.
        domains: DomainGeometry.
        annotations: AnnotationConfig (for "groin"/"groin_span").

    Returns:
        (pad_lo, pad_hi, axis_kind, x_lo, x_hi, tag):
            pad_lo, pad_hi: padded-array slice bounds (pad_hi exclusive).
            axis_kind: "gis" (x-axis in GIS domain IDs) or "pad" (padded
                index).
            x_lo, x_hi: inclusive x-axis limits in whichever space
                axis_kind names.
            tag: short string for the output filename.
    """
    if isinstance(spec, str):
        key = spec.strip().lower()
        if key == "real":
            return (domains.start_real_index, domains.end_real_index, "gis",
                    domains.first_gis_id, domains.last_gis_id,
                    f"D{domains.first_gis_id}-{domains.last_gis_id}")
        if key == "all":
            return (0, domains.total_domains, "pad",
                    0, domains.total_domains - 1,
                    f"ALL{domains.total_domains}pad")
        if key == "groin":
            if groin is None:
                raise ValueError(
                    "range='groin' must be resolved per-structure; call it via "
                    "make_all_shoreline_gifs, or pass groin=(name, domain)."
                )
            gname, gdom = groin
            gis_lo, gis_hi = _groin_zoom_window(gdom, pad, domains)
            return (domains.gis_to_pad(gis_lo), domains.gis_to_pad(gis_hi) + 1, "gis",
                    gis_lo, gis_hi, f"groinZoom_{_slug(gname)}_D{gis_lo}-{gis_hi}")
        if key == "groin_span":
            # One window enclosing EVERY groin -- useful once structures interact.
            doms = list(annotations.groins.values())
            if not doms:
                raise ValueError("range='groin_span' requires entries in annotations.groins.")
            gis_lo, _ = _groin_zoom_window(min(doms), pad, domains)
            _, gis_hi = _groin_zoom_window(max(doms), pad, domains)
            return (domains.gis_to_pad(gis_lo), domains.gis_to_pad(gis_hi) + 1, "gis",
                    gis_lo, gis_hi, f"groinSpan_D{gis_lo}-{gis_hi}")
        raise ValueError(
            f"GIF range={spec!r} not understood; use 'real', 'all', 'groin', "
            f"'groin_span', or (lo, hi)."
        )

    gis_lo, gis_hi = int(spec[0]), int(spec[1])
    if gis_lo > gis_hi:
        gis_lo, gis_hi = gis_hi, gis_lo
    if gis_lo < domains.first_gis_id or gis_hi > domains.last_gis_id:
        raise ValueError(
            f"GIF range=({gis_lo}, {gis_hi}) falls outside the real domains "
            f"{domains.first_gis_id}-{domains.last_gis_id}. Use 'all' to include buffers."
        )
    return (domains.gis_to_pad(gis_lo), domains.gis_to_pad(gis_hi) + 1, "gis",
            gis_lo, gis_hi, f"D{gis_lo}-{gis_hi}")


def make_shoreline_gif(shoreline_m, run,
                        domain_range="real", mode="displacement", pad=9, groin=None,
                        baseline_m=None, baseline_label=None,
                        target_m=None, target_label=None,
                        relocations=None,
                        fps=None, stride=None, annotate=None, auto_open=None,
                        keep_frames=None,
                        domains=DEFAULT_DOMAINS,
                        annotations=DEFAULT_ANNOTATIONS,
                        gif_config=DEFAULT_GIF_CONFIG):
    """Plan-view animation of the shoreline, one frame per model year.

    Current shoreline drawn over a year-0 reference (dashed grey); shaded
    blue where seaward of the reference, red where landward. Axis
    orientation is set by gif_config.ocean_at_bottom: ocean at bottom /
    sound at top by default, matching Hatteras' real cross-shore layout.

    SIGN CONVENTION: build_shoreline_matrix() returns raw x_s_TS (x10 for
    m), NO sign flip applied -- in BRIE/Barrier3D, x_s_TS INCREASES as the
    shoreline retreats landward. run.flip_sign_model is applied internally
    here so larger = more seaward, matching compute_change_rate's sign
    handling. Do not pre-flip shoreline_m before passing it in.

    MODE:
      "displacement" (default) - each domain's change from its OWN year-0
        position; all start at 0, so a few m/yr stays legible across all 90
        domains. Use for the full island.
      "position" - cross-shore position relative to year-0 alongshore mean.
        Keeps planform shape, but x_s_TS is referenced to each domain's own
        local grid, so the alongshore "shape" is partly an artefact of how
        the initial DEMs were extracted, and ~80 m of 40-yr change can be
        nearly invisible against a planform range of several hundred m.
        Fine for a narrow groin-zone window like (1, 30); use
        "displacement" for the full island.
      "difference" - this run minus baseline_m, in displacement terms:
        (run - run_yr0) - (base - base_yr0). Differencing displacements
        (not raw positions) cancels any per-domain initial-condition offset
        between runs, so year 0 is exactly zero. With a no-groin baseline
        this isolates the groin's effect: blue = seaward of baseline
        (updrift fillet), red = landward (downdrift notch).

    Args:
        shoreline_m: 2-D array [n_years, domains.total_domains] from
            cascade_pipeline.shoreline.build_shoreline_matrix(cascade,
            to_meters=True); raw x_s_TS convention.
        run: RunInfo for this run.
        domain_range: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi).
        pad: half-width in domains; only read when domain_range is
            "groin"/"groin_span".
        groin: (name, domain) tuple, required for domain_range="groin".
        baseline_m: same-shape matrix from a previous run (same raw x_s_TS
            convention); required for mode="difference".
        target_m: 1-D array of length domains.total_domains giving an
            observed/target shoreline position to draw as a static line in
            every frame -- e.g. the surveyed end-year island position, so
            "how close did we get" is readable off the animation. Same raw
            x_s_TS convention as shoreline_m (meters, increasing landward);
            run.flip_sign_model and the mode's own referencing are applied to
            it exactly as they are to the run, so it lands in the same frame
            as the shoreline. Ignored by mode="difference", where the y-axis
            is a run-minus-run quantity the target has no place in.
        target_label: legend label for that line; defaults to
            gif_config.target_label.
        Everything else defaults to gif_config.

    Returns:
        Path to the saved GIF, or None if the job was skipped.
    """
    try:
        from PIL import Image
    except Exception as e:  # Pillow ships with matplotlib, so this is unlikely
        print(f"  [GIF] Pillow unavailable ({e}); skipping animation.")
        return None

    fps = gif_config.fps if fps is None else fps
    stride = gif_config.year_stride if stride is None else stride
    annotate = gif_config.annotate if annotate is None else annotate
    auto_open = gif_config.auto_open if auto_open is None else auto_open
    keep_frames = gif_config.keep_frames if keep_frames is None else keep_frames
    baseline_label = gif_config.baseline_label if baseline_label is None else baseline_label
    target_label = gif_config.target_label if target_label is None else target_label

    mode = str(mode).strip().lower()
    if mode not in ("position", "displacement", "difference"):
        raise ValueError(
            f"GIF mode={mode!r} must be 'position', 'displacement', or 'difference'."
        )
    if mode == "difference" and baseline_m is None:
        print("  [GIF] mode='difference' needs a baseline run "
              "(pass baseline_m); skipping this job.")
        return None

    pad_lo, pad_hi, axis_kind, x_lo, x_hi, range_tag = _resolve_gif_domain_range(
        domain_range, pad, groin, domains, annotations
    )

    raw = np.asarray(shoreline_m, dtype=float)
    if raw.ndim != 2 or raw.shape[1] < pad_hi:
        print(f"  [GIF] shoreline matrix shape {raw.shape} too small for "
              f"padded slice [{pad_lo}:{pad_hi}]; skipping animation.")
        return None

    # x_s_TS increases landward -> flip so that larger = more seaward.
    flip = -1.0 if run.flip_sign_model else 1.0
    pos = (flip * raw)[:, pad_lo:pad_hi]
    n_years = pos.shape[0]
    if n_years < 2:
        print("  [GIF] fewer than 2 model years; skipping animation.")
        return None

    x = np.arange(x_lo, x_hi + 1)
    init = pos[0].copy()

    if mode == "difference":
        base = np.asarray(baseline_m, dtype=float)
        if base.shape != raw.shape:
            print(f"  [GIF] baseline shape {base.shape} != run shape {raw.shape}; "
                  f"the two runs must share domain count and length. Skipping.")
            return None
        base_pos = (flip * base)[:, pad_lo:pad_hi]
        # Difference the DISPLACEMENTS so any initial-condition offset cancels.
        series = (pos - init[None, :]) - (base_pos - base_pos[0][None, :])
        up_word = "landward" if gif_config.ocean_at_bottom else "seaward"
        ylabel = f"\u0394 shoreline vs {baseline_label} (m)\n{up_word} \u25b2"
        lbl_up, lbl_dn = (f"Seaward of {baseline_label}", f"Landward of {baseline_label}")
        lbl_ref = "No difference"
    elif mode == "displacement":
        series = pos - init[None, :]
        up_word = "landward" if gif_config.ocean_at_bottom else "seaward"
        ylabel = f"Shoreline displacement since year 0 (m)\n{up_word} \u25b2"
        lbl_up, lbl_dn = "Accretion (seaward of year 0)", "Erosion (landward of year 0)"
        lbl_ref = f"Year 0 ({run.start_year}) shoreline"
    else:
        series = pos - np.nanmean(init)
        up_word = "landward" if gif_config.ocean_at_bottom else "seaward"
        ylabel = f"Cross-shore position (m, rel. year-0 mean)\n{up_word} \u25b2"
        lbl_up, lbl_dn = "Accretion (seaward of year 0)", "Erosion (landward of year 0)"
        lbl_ref = f"Year 0 ({run.start_year}) shoreline"

    ref = series[0]  # year-0 reference; identically zero except in "position"

    # -- Optional observed/target line ----------------------------------------
    # Referenced exactly as the run is, so it shares the run's axis. "position"
    # subtracts the same year-0 alongshore mean; "displacement" subtracts the
    # same per-domain year-0 position, which turns the target into the
    # displacement the run is trying to reproduce. Both lines therefore sit on
    # the same year-0 base, so the gap between them is the model's misfit even
    # though that base is itself an initial condition.
    target = None
    if target_m is not None:
        tgt = np.asarray(target_m, dtype=float).ravel()
        if mode == "difference":
            print("  [GIF] mode='difference' plots a run-minus-run axis; "
                  "target_m has no place on it and is not drawn.")
        elif tgt.size != raw.shape[1]:
            print(f"  [GIF] target_m has {tgt.size} values but the run has "
                  f"{raw.shape[1]} domains; drawing without the target line.")
        else:
            tgt_pos = (flip * tgt)[pad_lo:pad_hi]
            target = tgt_pos - (init if mode == "displacement" else np.nanmean(init))

    # Fixed axes across all frames so the shoreline doesn't jitter frame-to-frame.
    spread = series if target is None else np.vstack([series, target[None, :]])
    ymin, ymax = float(np.nanmin(spread)), float(np.nanmax(spread))
    ypad = (ymax - ymin) * 0.10 or 1.0
    if gif_config.ocean_at_bottom:
        ylim = (ymax + ypad, ymin - ypad)
    else:
        ylim = (ymin - ypad, ymax + ypad)

    n_dom = len(x)
    figsize = (float(np.clip(6.0 + 0.115 * n_dom, 9.0, 18.0)), 5.0)

    frames_dir = os.path.join(run.run_dir, "gif_frames")
    if keep_frames:
        os.makedirs(frames_dir, exist_ok=True)

    year_idx = list(range(0, n_years, max(int(stride), 1)))
    if year_idx[-1] != n_years - 1:
        year_idx.append(n_years - 1)  # always land on the final year

    be_lbl = "on" if run.background_erosion_on else "off"
    hs_lbl = f"Hs={run.Hs} m  |  " if run.Hs is not None else ""

    frames = []
    for t in year_idx:
        cur = series[t]

        fig, ax = plt.subplots(figsize=figsize, dpi=110)
        # Fixed margins (NOT bbox_inches="tight") so every frame is byte-identical
        # in size -- mismatched frame dimensions break GIF assembly.
        fig.subplots_adjust(left=0.085, right=0.985, top=0.90, bottom=0.235)
        fig.patch.set_facecolor("white")
        ax.set_facecolor("white")

        # -- Geographic context ------------------------------------------------
        if axis_kind == "gis" and annotate:
            add_geographic_annotations(ax, annotations)
        elif axis_kind == "pad":
            ax.axvspan(-0.5, domains.start_real_index - 0.5, color="red", alpha=0.10, zorder=0)
            ax.axvspan(domains.end_real_index - 0.5, domains.total_domains - 0.5,
                       color="red", alpha=0.10, zorder=0)
            ax.axvline(domains.start_real_index - 0.5, color="k", ls="--", lw=1.0, alpha=0.5, zorder=2)
            ax.axvline(domains.end_real_index - 0.5, color="k", ls="--", lw=1.0, alpha=0.5, zorder=2)
            if annotate:
                trans_pad = blended_transform_factory(ax.transData, ax.transAxes)
                for span_label, (d_lo, d_hi) in annotations.town_spans.items():
                    ax.axvspan(domains.gis_to_pad(d_lo) - 0.5, domains.gis_to_pad(d_hi) + 0.5,
                               color=annotations.color_town_span, alpha=0.14, zorder=0)
                    ax.text((domains.gis_to_pad(d_lo) + domains.gis_to_pad(d_hi)) / 2.0, 0.90,
                            span_label, transform=trans_pad, ha="center", va="top",
                            fontsize=8, color="0.25", fontweight="bold",
                            bbox=dict(boxstyle="round,pad=0.2", fc="white",
                                      ec="none", alpha=0.85))
                for gname, dom in annotations.groins.items():
                    ax.axvline(domains.gis_to_pad(dom), color=annotations.color_groin, lw=1.1, ls=":",
                               alpha=0.85, zorder=2)
                    ax.text(domains.gis_to_pad(dom), annotations.groin_label_y, gname,
                            transform=trans_pad, ha="center", va="top", fontsize=7,
                            color=annotations.color_groin, rotation=90,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                      ec="none", alpha=0.80))

        # -- Shoreline + erosion/accretion shading -----------------------------
        ax.fill_between(x, cur, ref, where=(cur >= ref), interpolate=True,
                        color="#1565C0", alpha=0.28, zorder=1, label=lbl_up)
        ax.fill_between(x, cur, ref, where=(cur < ref), interpolate=True,
                        color="#B71C1C", alpha=0.22, zorder=1, label=lbl_dn)
        ax.plot(x, ref, color="0.70", ls="--", lw=1.2, zorder=3, label=lbl_ref)
        if target is not None:
            ax.plot(x, target, color="#00897B", ls="-.", lw=1.6, zorder=3.5,
                    label=target_label)
        ax.plot(x, cur, color="#1a2a3a", lw=2.0, zorder=4, label="Shoreline")

        # -- Roadway relocations, in the year they happen ----------------------
        # An EVENT, not a state: _road_relocated_TS is raised in the year the
        # roadway manager moves the road, so a frame shows that year's moves
        # only. Drawn on the shoreline itself because that is the line whose
        # arrival at the road caused them.
        if relocations is not None and t < len(relocations):
            moved = np.flatnonzero(np.asarray(relocations[t], dtype=bool))
            in_view = [p for p in moved if pad_lo <= p < pad_hi]
            if in_view:
                mx = [x_lo + (p - pad_lo) for p in in_view]
                my = [cur[p - pad_lo] for p in in_view]
                ax.plot(mx, my, linestyle="none", marker="v", markersize=9,
                        color=RELOCATION_COLOR, markeredgecolor="white",
                        markeredgewidth=0.7, zorder=6,
                        label=f"relocated ({len(in_view)})")

        # -- Axes --------------------------------------------------------------
        ax.set_xlim(x_lo - 0.5, x_hi + 0.5)
        ax.set_ylim(*ylim)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(10 if n_dom > 40 else 5))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5 if n_dom > 40 else 1))
        ax.tick_params(axis="both", which="major", labelsize=9, direction="in", length=5)
        ax.tick_params(axis="both", which="minor", direction="in", length=3)
        ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
        ax.spines[["top", "right"]].set_visible(False)

        if axis_kind == "pad":
            ax.set_xlabel(
                f"{run.model_name} domain index (buffers included, 0\u2013{domains.total_domains - 1})"
                f"   |   \u2190 {annotations.low_end_label}      {annotations.high_end_label} \u2192",
                fontsize=11, fontweight="bold")
            top_ax = ax.secondary_xaxis("top")
            tp, tl = [], []
            for gid in range(domains.first_gis_id, domains.last_gis_id + 1, 10):
                tp.append(domains.gis_to_pad(gid))
                tl.append(str(gid))
            top_ax.set_xticks(tp)
            top_ax.set_xticklabels(tl, fontsize=8)
            top_ax.set_xlabel(f"GIS Domain ID ({domains.first_gis_id}\u2013{domains.last_gis_id})",
                              fontsize=9)
        else:
            ax.set_xlabel(
                f"{run.model_name} Model Domain ({int(domains.domain_spacing_m)} m alongshore)"
                f"   |   \u2190 {annotations.low_end_label}      {annotations.high_end_label} \u2192",
                fontsize=11, fontweight="bold")

        ax.set_ylabel(ylabel, fontsize=10.5, fontweight="bold")

        # -- Year stamp + title + legend ---------------------------------------
        ax.text(0.985, 0.045, f"{run.start_year + t}  (year {t})",
                transform=ax.transAxes, ha="right", va="bottom",
                fontsize=13, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7"))
        ax.set_title(
            (f"Shoreline difference \u2014 {run.run_name} minus {baseline_label}  |  "
             f"{run.start_year}\u2013{run.end_year}  |  {hs_lbl}BE={be_lbl}"
             if mode == "difference" else
             f"Shoreline evolution \u2014 {annotations.region_name}  |  "
             f"{run.start_year}\u2013{run.end_year}  |  "
             f"{hs_lbl}BE={be_lbl}  |  {run.run_name}"),
            fontsize=11, fontweight="bold", color="#1a2a3a", pad=(14 if axis_kind == "gis" else 26))
        # Legend outside (below) the axes so it never covers the town / shoal
        # labels, which sit at fixed axes fractions inside the panel.
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.115),
                  bbox_transform=ax.transAxes, fontsize=8.5, framealpha=0.95,
                  edgecolor="#cccccc", frameon=True,
                  ncol=5 if target is not None else 4)

        # -- Render to an in-memory PNG (constant size, no disk churn) ----------
        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor="white")
        buf.seek(0)
        frames.append(Image.open(buf).convert("RGB"))
        if keep_frames:
            fig.savefig(os.path.join(frames_dir, f"frame_{t:03d}.png"),
                        facecolor="white")
        plt.close(fig)

    sizes = {f.size for f in frames}
    if len(sizes) > 1:
        print(f"  [GIF] frame sizes differ ({sizes}); aborting animation.")
        return None

    gif_path = os.path.join(run.run_dir, f"{run.run_name}_shoreline_{mode}_{range_tag}.gif")
    frames[0].save(gif_path, save_all=True, append_images=frames[1:],
                   duration=int(1000.0 / max(float(fps), 0.1)), loop=0, optimize=True)

    print(f"  [GIF] saved shoreline animation ({len(frames)} frames, {mode}, "
          f"{range_tag}): {gif_path}")
    if keep_frames:
        print(f"  [GIF] frames kept in: {frames_dir}")
    if auto_open:
        _open_file(gif_path)
    return gif_path


def make_all_shoreline_gifs(shoreline_m, run, jobs,
                             baseline_npy=None, target_m=None,
                             relocations=None,
                             domains=DEFAULT_DOMAINS,
                             annotations=DEFAULT_ANNOTATIONS,
                             gif_config=DEFAULT_GIF_CONFIG):
    """Run every job in `jobs` for one completed run.

    Also saves the shoreline matrix (when gif_config.save_matrix) so this
    run can serve as the baseline for a later difference GIF. Each job is
    independent: one failing (bad range, missing baseline) logs and is
    skipped rather than taking the others -- or the run -- down with it.

    Args:
        shoreline_m: [n_years, domains.total_domains] array, raw x_s_TS
            convention (from build_shoreline_matrix, to_meters=True).
        run: RunInfo for this run.
        jobs: List of dicts, same shape as the original GIF_JOBS list, e.g.
            [dict(range="real", mode="displacement"),
             dict(range="groin", mode="difference", pad=9)]
        baseline_npy: Path to a previous run's *_shoreline_matrix.npy to
            difference against. None skips "difference" jobs.
        target_m: Observed/target shoreline position (raw x_s_TS convention,
            meters, one value per padded domain) drawn as a static line on
            every "position"/"displacement" job. A job may override it with
            its own "target" key, or opt out with "target": False.
        relocations: Optional (n_years, domains.total_domains) boolean array,
            True in the year a domain's roadway is relocated. Passed to every
            job; each marks only the domains inside its own window. None -- the
            default, and what a run without roadway management should pass --
            draws no markers and leaves the GIFs exactly as they were.

    Returns:
        List of saved GIF paths.
    """
    if gif_config.save_matrix:
        npy_path = os.path.join(run.run_dir, f"{run.run_name}_shoreline_matrix.npy")
        np.save(npy_path, np.asarray(shoreline_m, dtype=float))
        print(f"  Saved shoreline matrix: {npy_path}")

    # Load the baseline once, not per job.
    baseline_m = None
    wants_diff = any(str(j.get("mode", "")).lower() == "difference" for j in jobs)
    if wants_diff:
        if not baseline_npy:
            print("  [GIF] baseline_npy is None; 'difference' jobs will be skipped.")
        elif not os.path.exists(baseline_npy):
            print(f"  [GIF] baseline not found: {baseline_npy}; "
                  f"'difference' jobs will be skipped.")
        else:
            baseline_m = np.load(baseline_npy)
            print(f"  [GIF] baseline loaded ({baseline_m.shape}): {baseline_npy}")

    made = []
    for job in jobs:
        # A "groin" job fans out into one GIF per structure in
        # annotations.groins, so adding a second groin needs no change here
        # or in `jobs`. Use "which" to restrict, or range="groin_span" for
        # one enclosing window.
        if str(job.get("range", "")).strip().lower() == "groin":
            try:
                targets = [{"groin": g} for g in _select_groins(job.get("which", "all"), annotations)]
            except Exception as e:
                print(f"  [GIF] job {job} skipped ({e}).")
                continue
        else:
            targets = [{"groin": None}]

        # A job's own "target" wins over the shared one; "target": False opts
        # out, which is how a job stays a clean model-only figure.
        job_target = job.get("target", target_m)
        if job_target is False:
            job_target = None

        for extra in targets:
            try:
                path = make_shoreline_gif(
                    shoreline_m, run,
                    domain_range=job.get("range", "real"),
                    mode=job.get("mode", "displacement"),
                    pad=job.get("pad", 9),
                    groin=extra["groin"],
                    baseline_m=baseline_m,
                    target_m=job_target, target_label=job.get("target_label"),
                    relocations=relocations,
                    fps=job.get("fps"), stride=job.get("stride"),
                    annotate=job.get("annotate"), auto_open=job.get("auto_open"),
                    keep_frames=job.get("keep_frames"),
                    domains=domains, annotations=annotations, gif_config=gif_config,
                )
                if path:
                    made.append(path)
            except Exception as e:
                label = f"{job} [{extra['groin'][0]}]" if extra["groin"] else str(job)
                print(f"  [GIF] job {label} failed ({e}); continuing.")
    return made
