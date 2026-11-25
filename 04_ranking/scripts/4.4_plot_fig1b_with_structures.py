#!/usr/bin/env python3
"""Step 4.4 script: draw Fig. 1B and Fig. 1C with representative AF3 structures and PAE.

Combines the structural ranking tracks from Step 4.1 with representative
ChimeraX multimer renders (identical 1:1 physical scale across all 8 views)
and AF3 PAE matrices from Step 4.3 for ranks #1 (XP_002449692.2),
#18 (NbNRC2a), #128 (Pi5-3), and #297 (XP_052135576.1).
"""

import argparse
import sys
import tempfile
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MplPath
import numpy as np
import pandas as pd
from PIL import Image

BAR_COLORS = {"D_APEX": "#C45B47", "ipTM_LCB": "#15415C"}
CLADE_LANES = [
    ("Pia-2/MIC1", "Pia-2/MIC1", "#2A8C82"),
    ("NRC-S",      "NRC-S",      "#6B5AA6"),
    ("Pii-2",      "Pii-2",      "#D69F2D"),
]
LANE_BASELINE = "#D2D6DC"

PLDDT_BINS = [
    ("#0053D6", "> 90"),
    ("#65CBF3", "70\u201390"),
    ("#FFDB13", "50\u201370"),
    ("#FF7D45", "< 50"),
]

PAE_CMAP = LinearSegmentedColormap.from_list(
    "af3_reference_green",
    ["#087f3a", "#2fa24f", "#74c476", "#c7e9c0", "#f7fcf5"],
    N=256,
)

SELECTED_MODELS = [
    {
        "rank": 1,
        "expected_protein_id": "XP_002449692.2",
        "title": "#1  XP_002449692.2",
        "subtitle": "Top-ranked",
        "subtitle_color": "#4B5563",
        "accent_color": "#15415C",
        "top_img": "rank1_top.png",
        "side_img": "rank1_side.png",
        "pae_img": "4.3_rank001_XP_002449692_pae.png",
    },
    {
        "rank": 18,
        "expected_protein_id": "NbNRC2a",
        "title": "#18  NbNRC2a",
        "subtitle": "Benchmark",
        "subtitle_color": "#4B5563",
        "accent_color": "#15415C",
        "top_img": "NbNRC2_top.png",
        "side_img": "NbNRC2_side.png",
        "pae_img": "4.3_rank018_NbNRC2a_pae.png",
    },
    {
        "rank": 128,
        "expected_protein_id": "XP_015612169.1",
        "title": "#128  Pi5-3",
        "subtitle": "Pii-2 sensor clade",
        "subtitle_color": "#B37D14",
        "accent_color": "#D69F2D",
        "top_img": "Pi5-3_top.png",
        "side_img": "Pi5-3_side.png",
        "pae_img": "4.3_rank128_Pi5-3_pae.png",
    },
    {
        "rank": 297,
        "expected_protein_id": "XP_052135576.1",
        "title": "#297  XP_052135576.1",
        "subtitle": "Pia-2/MIC1 sensor clade",
        "subtitle_color": "#20736A",
        "accent_color": "#2A8C82",
        "top_img": "rank297_top.png",
        "side_img": "rank297_side.png",
        "pae_img": "4.3_rank297_XP_052135576_pae.png",
    },
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", type=Path,
                   default=Path("04_ranking/data/results/4.1_ranked_scores.tsv"),
                   help="Ranked structural score TSV from Step 4.1.")
    p.add_argument("--chimerax-dir", type=Path,
                   default=Path("04_ranking/data/inputs/4.4_chimerax_images"),
                   help="Directory containing ChimeraX top/side renders.")
    p.add_argument("--pae-dir", type=Path,
                   default=Path("04_ranking/data/results"),
                   help="Directory containing Step 4.3 PAE PNGs.")
    p.add_argument("--output-prefix", type=Path,
                   default=Path("04_ranking/data/results/4.4_fig1b_structures"),
                   help="Output prefix (writes .pdf, .svg, .png).")
    p.add_argument("--dapex-color", default=BAR_COLORS["D_APEX"],
                   help="Bar color for D_APEX track (default: %(default)s).")
    p.add_argument("--iptm-color", default=BAR_COLORS["ipTM_LCB"],
                   help="Bar color for ipTM track (default: %(default)s).")
    p.add_argument("--dpi", type=int, default=300,
                   help="Raster output DPI (default: %(default)s).")
    return p.parse_args()


def _register_macos_helvetica_neue(tmp_font_dir: Path):
    """Optionally extract macOS HelveticaNeue.ttc faces if available; ignore otherwise."""
    ttc_path = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if not ttc_path.exists():
        return
    try:
        from fontTools.ttLib import TTCollection
        ttc = TTCollection(str(ttc_path))
        for idx, fname in (
            (0, "HelveticaNeue-Regular.ttf"),
            (1, "HelveticaNeue-Bold.ttf"),
            (10, "HelveticaNeue-Medium.ttf"),
        ):
            if idx < len(ttc.fonts):
                out_ttf = tmp_font_dir / fname
                ttc.fonts[idx].save(str(out_ttf))
                mpl.font_manager.fontManager.addfont(str(out_ttf))
    except Exception:
        pass


def configure_style(tmp_font_dir: Path):
    _register_macos_helvetica_neue(tmp_font_dir)
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.default": "regular",
        "axes.linewidth": 0.75,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def load_aligned_cropped_images(img_paths, align_y="center", white_thresh=248, pad_px=24):
    """Tightly crop each ChimeraX render individually, then pad onto a uniform white canvas.

    Preserves true 1:1 molecular pixel scale across all models while keeping
    consistent vertical alignment (center for top view, bottom baseline for side view).
    """
    crops = []
    for path in img_paths:
        im = Image.open(path).convert("RGB")
        arr = np.asarray(im)
        mask = (arr < white_thresh).any(axis=2)
        ys, xs = np.where(mask)
        if len(xs) == 0:
            raise RuntimeError(f"no non-white content found in {path}")
        crop = im.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1))
        crops.append(crop)

    max_w = max(c.width for c in crops)
    max_h = max(c.height for c in crops)
    canvas_w = max_w + 2 * pad_px
    canvas_h = max_h + 2 * pad_px

    out = []
    for crop in crops:
        canvas = Image.new("RGB", (canvas_w, canvas_h), (255, 255, 255))
        dx = (canvas_w - crop.width) // 2
        if align_y == "center":
            dy = (canvas_h - crop.height) // 2
        elif align_y == "bottom":
            dy = canvas_h - pad_px - crop.height
        else:
            raise ValueError(f"unsupported align_y={align_y}")
        canvas.paste(crop, (dx, dy))
        out.append(np.asarray(canvas))
    return out


def _group_contiguous(indices, max_gap=2):
    groups = []
    for idx in indices:
        val = int(idx)
        if not groups or val > groups[-1][-1] + max_gap:
            groups.append([val])
        else:
            groups[-1].append(val)
    return [(g[0], g[-1]) for g in groups]


def load_cropped_pae(path: Path):
    """Detect the exact square PAE matrix box and chain border fractions."""
    im = Image.open(path).convert("RGB")
    arr = np.asarray(im)
    dark = (arr[:, :, 0] < 110) & (arr[:, :, 1] < 110) & (arr[:, :, 2] < 110)
    h_idx = np.where(dark.sum(axis=1) > arr.shape[1] * 0.5)[0]
    v_idx = np.where(dark.sum(axis=0) > arr.shape[0] * 0.5)[0]
    h_groups = _group_contiguous(h_idx)
    v_groups = _group_contiguous(v_idx)
    if len(h_groups) != 8 or len(v_groups) != 8:
        raise RuntimeError(
            f"{path}: expected 8 horizontal and 8 vertical chain borders, "
            f"got h={len(h_groups)}, v={len(v_groups)}"
        )
    y0, y1 = h_groups[0][0], h_groups[-1][1]
    x0, x1 = v_groups[0][0], v_groups[-1][1]
    cropped = im.crop((x0, y0, x1 + 1, y1 + 1))

    box_w = float((x1 + 1) - x0)
    v_centers = [0.5 * (g[0] + g[1]) - x0 for g in v_groups]
    chain_spans = []
    for i in range(6):
        c_mid = 0.5 * (v_centers[i] + v_centers[i + 1]) / box_w
        chain_spans.append((str(i + 1), c_mid))
    lig_mid = 0.5 * (v_centers[6] + v_centers[7]) / box_w
    chain_spans.append(("oleic acids", lig_mid))
    return np.asarray(cropped), chain_spans


def draw_s_connector(ax, x0, y0, x1, y1, color="#4B5563", lw=0.90, stem=0.05, zorder=3):
    """Vertical drop from guide line bottom, then smooth cubic Bezier S-curve to column."""
    y_stem = y0 - stem
    y_mid = 0.5 * (y_stem + y1)
    verts = [
        (x0, y0),
        (x0, y_stem),
        (x0, y_mid),
        (x1, y_mid),
        (x1, y1),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
    ]
    ax.add_patch(patches.PathPatch(
        MplPath(verts, codes),
        facecolor="none", edgecolor=color, lw=lw, capstyle="round", zorder=zorder,
    ))


def add_in_axes(fig, W, H, x0, y0, w, h, sharex=None):
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H], sharex=sharex)


def build_figure(args, tmp_font_dir: Path):
    configure_style(tmp_font_dir)

    df = pd.read_csv(args.table, sep="\t", keep_default_na=False)
    for col in ("weighted_rank", "D_APEX", "ipTM_LCB", "sensor_clade", "protein_ID"):
        if col not in df.columns:
            sys.exit(f"error: {args.table} missing column '{col}'")
    df = df.sort_values("weighted_rank").reset_index(drop=True)
    x_max = int(df["weighted_rank"].max())

    model_records = []
    for spec in SELECTED_MODELS:
        row_df = df.loc[df["weighted_rank"].eq(spec["rank"])]
        if len(row_df) != 1:
            sys.exit(f"error: rank {spec['rank']} not uniquely found in {args.table}")
        row = row_df.iloc[0]
        if row["protein_ID"] != spec["expected_protein_id"]:
            sys.exit(
                f"error: rank {spec['rank']} protein_ID={row['protein_ID']} "
                f"!= expected {spec['expected_protein_id']}"
            )
        rec = dict(spec)
        rec["d_apex"] = float(row["D_APEX"])
        rec["iptm"] = float(row["ipTM_LCB"])
        rec["sensor_clade"] = row["sensor_clade"]
        model_records.append(rec)

    top_paths = [args.chimerax_dir / m["top_img"] for m in model_records]
    side_paths = [args.chimerax_dir / m["side_img"] for m in model_records]
    pae_paths = [args.pae_dir / m["pae_img"] for m in model_records]
    for p in top_paths + side_paths + pae_paths:
        if not p.exists():
            sys.exit(f"error: missing image file {p}")

    top_imgs = load_aligned_cropped_images(top_paths, align_y="center", pad_px=24)
    side_imgs = load_aligned_cropped_images(side_paths, align_y="bottom", pad_px=24)
    pae_data = [load_cropped_pae(p) for p in pae_paths]

    W, H = 9.15, 7.46
    x_left, x_right = 0.64, 8.17
    w_main = x_right - x_left  # 7.53 in
    x_lbl_right = 0.56
    x_leg = 8.26

    # Top rank axis header & tracks vertical positions
    y_top_title = H - 0.14  # 7.32 in
    y_apex_top = 6.94
    h_apex = 0.38
    y_apex_bot = y_apex_top - h_apex  # 6.56 in

    gap_tracks = 0.065
    h_iptm = 0.345
    y_iptm_top = y_apex_bot - gap_tracks  # 6.495 in
    y_iptm_bot = y_iptm_top - h_iptm      # 6.15 in

    gap_clade = 0.24
    h_clade = 0.56
    y_clade_top = y_iptm_bot - gap_clade  # 5.91 in
    y_clade_bot = y_clade_top - h_clade   # 5.35 in

    fig = plt.figure(figsize=(W, H), dpi=args.dpi, facecolor="white")

    ax_ov = fig.add_axes([0, 0, 1, 1], facecolor="none")
    ax_ov.set_xlim(0, W)
    ax_ov.set_ylim(0, H)
    ax_ov.axis("off")

    # Track 1: D_APEX (A) with top rank x-axis
    ax_apex = add_in_axes(fig, W, H, x_left, y_apex_bot, w_main, h_apex)
    ax_apex.set_facecolor("none")
    ranks = df["weighted_rank"].to_numpy()
    ax_apex.bar(ranks, df["D_APEX"], width=1.0, color=args.dapex_color, edgecolor="none", zorder=2)
    ax_apex.set_xlim(0.5, x_max + 0.5)
    ax_apex.set_ylim(0, float(df["D_APEX"].max()) * 1.12)
    ax_apex.set_yticks([])
    for s in ("right", "left", "bottom"):
        ax_apex.spines[s].set_visible(False)
    ax_apex.spines["top"].set_visible(True)
    ax_apex.spines["top"].set_linewidth(0.75)
    ax_apex.spines["top"].set_color("#1A1A1A")
    ax_apex.xaxis.set_ticks_position("top")

    selected_ranks = [m["rank"] for m in model_records]
    base_ticks = [50, 100, 150, 200, 250]
    all_ticks = sorted(set(selected_ranks + base_ticks))
    ax_apex.set_xticks(all_ticks)
    ax_apex.set_xticklabels([str(t) for t in all_ticks], fontsize=7.8)
    ax_apex.tick_params(
        axis="x", top=True, bottom=False, labeltop=True, labelbottom=False,
        direction="out", length=3.2, width=0.75, pad=2.2, colors="#4B5563",
    )
    for lbl, tval in zip(ax_apex.get_xticklabels(), all_ticks):
        if tval in selected_ranks:
            lbl.set_fontweight("bold")
            lbl.set_color("#111827")
            lbl.set_fontsize(8.0)

    ax_ov.text(
        0.5 * (x_left + x_right), y_top_title, "Composite structural rank",
        ha="center", va="center", fontsize=8.8, color="#1A1A1A", zorder=10,
    )

    # Track 2: ipTM
    ax_iptm = add_in_axes(fig, W, H, x_left, y_iptm_bot, w_main, h_iptm, sharex=ax_apex)
    ax_iptm.set_facecolor("none")
    ax_iptm.bar(ranks, df["ipTM_LCB"], width=1.0, color=args.iptm_color, edgecolor="none", zorder=2)
    ax_iptm.set_xlim(0.5, x_max + 0.5)
    ax_iptm.set_ylim(0, float(df["ipTM_LCB"].max()) * 1.06)
    ax_iptm.set_yticks([])
    ax_iptm.tick_params(axis="x", which="both", top=False, bottom=False, labeltop=False, labelbottom=False, length=0)
    for s in ("top", "right", "left", "bottom"):
        ax_iptm.spines[s].set_visible(False)

    # Track 3: Sensor clades
    ax_clade = add_in_axes(fig, W, H, x_left, y_clade_bot, w_main, h_clade, sharex=ax_apex)
    ax_clade.set_facecolor("none")
    for y_idx, (key, display, color) in enumerate(CLADE_LANES):
        ax_clade.hlines(y_idx, 0.5, x_max + 0.5, color=LANE_BASELINE, lw=0.5, zorder=1)
        sel = df.loc[df["sensor_clade"].eq(key), "weighted_rank"]
        ax_clade.bar(sel, height=0.68, bottom=y_idx - 0.34, width=1.0, color=color, edgecolor="none", zorder=2)
        ax_clade.text(x_max + 0.5 + x_max * 0.012, y_idx, display, color=color,
                      fontsize=7.6, va="center", ha="left")
    ax_clade.text(x_max + 0.5 + x_max * 0.012, -0.78, "Sensor clades",
                  fontsize=8.0, fontweight="bold", color="#1A1A1A", va="center", ha="left")
    ax_clade.set_ylim(len(CLADE_LANES) - 0.45, -0.82)
    ax_clade.set_yticks([])
    ax_clade.tick_params(axis="x", which="both", top=False, bottom=False, labeltop=False, labelbottom=False, length=0)
    for s in ("top", "right", "left", "bottom"):
        ax_clade.spines[s].set_visible(False)

    # Left margin labels for top tracks
    ax_ov.text(x_lbl_right, 0.5 * (y_apex_bot + y_apex_top), r"$\mathrm{D_{APEX}}$ (Å)",
               ha="right", va="center", fontsize=8.4, color="#1A1A1A")
    ax_ov.text(x_lbl_right, 0.5 * (y_iptm_bot + y_iptm_top), "ipTM",
               ha="right", va="center", fontsize=8.4, color="#1A1A1A")

    for r in selected_ranks:
        for ax in (ax_apex, ax_iptm, ax_clade):
            ax.axvline(r, color="#374151", lw=0.72, ls=(0, (2.6, 2.0)),
                       alpha=0.80, zorder=4, clip_on=False)

    # Four structure + PAE columns
    col_gap = 0.21
    col_w = (w_main - 3 * col_gap) / 4.0  # 1.725 in
    col_centers = [x_left + k * (col_w + col_gap) + col_w / 2.0 for k in range(4)]

    y_col_top = 4.94

    def rank_to_in_x(rank_val):
        return x_left + ((rank_val - 0.5) / float(x_max)) * w_main

    for k, m in enumerate(model_records):
        xr = rank_to_in_x(m["rank"])
        xc = col_centers[k]
        draw_s_connector(
            ax_ov, xr, y_clade_bot, xc, y_col_top,
            color="#4B5563", lw=0.90, stem=0.05, zorder=3,
        )
        cx0 = x_left + k * (col_w + col_gap)
        cx1 = cx0 + col_w
        ax_ov.plot([cx0 + 0.04, cx1 - 0.04], [y_col_top, y_col_top],
                   color=m["accent_color"], lw=1.45, solid_capstyle="butt", zorder=5)
        ax_ov.plot(xc, y_col_top, marker="o", markersize=3.2,
                   color=m["accent_color"], markeredgewidth=0, zorder=6)

    y_hdr_title = 4.81
    y_hdr_sub   = 4.675
    y_hdr_stat  = 4.545

    y_top_hi, y_top_lo = 4.44, 3.10    # h = 1.34 in
    y_side_hi, y_side_lo = 3.02, 1.82  # h = 1.20 in
    pae_side = 1.44
    y_pae_lo = 0.12
    y_pae_hi = y_pae_lo + pae_side     # 1.56 in
    y_pae_chain_lbl = y_pae_hi + 0.025

    # Single shared physical scale across ALL 8 ChimeraX views (Top + Side)
    top_h_px, top_w_px = top_imgs[0].shape[:2]
    side_h_px, side_w_px = side_imgs[0].shape[:2]
    chimerax_px_per_in = max(
        top_w_px / (col_w * 0.94),
        top_h_px / (y_top_hi - y_top_lo),
        side_w_px / (col_w * 0.94),
        side_h_px / (y_side_hi - y_side_lo),
    )
    top_draw_w = top_w_px / chimerax_px_per_in
    top_draw_h = top_h_px / chimerax_px_per_in
    side_draw_w = side_w_px / chimerax_px_per_in
    side_draw_h = side_h_px / chimerax_px_per_in

    for k, m in enumerate(model_records):
        xc = col_centers[k]
        ax_ov.text(xc, y_hdr_title, m["title"], ha="center", va="center",
                   fontsize=8.6, fontweight="bold", color="#111827")
        ax_ov.text(xc, y_hdr_sub, m["subtitle"], ha="center", va="center",
                   fontsize=7.0, fontweight="medium", color=m["subtitle_color"])
        stat_txt = f"ipTM = {m['iptm']:.2f},  " + r"$\mathrm{D_{APEX}}$" + f" = {m['d_apex']:.1f} Å"
        ax_ov.text(xc, y_hdr_stat, stat_txt, ha="center", va="center",
                   fontsize=6.8, color="#4B5563")

        y_top_c = 0.5 * (y_top_hi + y_top_lo)
        ax_ov.imshow(
            top_imgs[k],
            extent=(
                xc - top_draw_w / 2.0, xc + top_draw_w / 2.0,
                y_top_c - top_draw_h / 2.0, y_top_c + top_draw_h / 2.0,
            ),
            aspect="auto",
            interpolation="antialiased",
            zorder=4,
        )

        y_side_c = 0.5 * (y_side_hi + y_side_lo)
        ax_ov.imshow(
            side_imgs[k],
            extent=(
                xc - side_draw_w / 2.0, xc + side_draw_w / 2.0,
                y_side_c - side_draw_h / 2.0, y_side_c + side_draw_h / 2.0,
            ),
            aspect="auto",
            interpolation="antialiased",
            zorder=4,
        )

        pae_img, chain_spans = pae_data[k]
        px0 = xc - pae_side / 2.0
        px1 = xc + pae_side / 2.0
        ax_ov.imshow(
            pae_img,
            extent=(px0, px1, y_pae_lo, y_pae_hi),
            aspect="auto",
            interpolation="nearest",
            zorder=4,
        )
        for lbl_txt, frac_x in chain_spans:
            fs = 5.25 if lbl_txt == "oleic acids" else 5.8
            ax_ov.text(
                px0 + frac_x * pae_side, y_pae_chain_lbl, lbl_txt,
                ha="center", va="bottom", fontsize=fs, color="#222222", zorder=6,
            )

    ax_ov.text(x_lbl_right, 0.5 * (y_top_hi + y_top_lo), "Top view",
               ha="right", va="center", fontsize=8.4, color="#1A1A1A")
    ax_ov.text(x_lbl_right, 0.5 * (y_side_hi + y_side_lo), "Side view",
               ha="right", va="center", fontsize=8.4, color="#1A1A1A")
    ax_ov.text(
        x_lbl_right, 0.5 * (y_pae_hi + y_pae_lo),
        "Predicted\naligned\nerror\n(PAE)",
        ha="right", va="center", fontsize=7.8, linespacing=1.22, color="#1A1A1A",
    )

    # Right margin legends
    y_plddt_top = 3.42
    ax_ov.text(x_leg, y_plddt_top, "pLDDT", ha="left", va="bottom",
               fontsize=8.0, fontweight="bold", color="#1A1A1A")
    sw_size = 0.11
    for i, (sw_col, sw_lbl) in enumerate(PLDDT_BINS):
        sy = y_plddt_top - 0.16 - i * 0.175
        ax_ov.add_patch(patches.FancyBboxPatch(
            (x_leg, sy - sw_size / 2.0), sw_size, sw_size,
            boxstyle="round,pad=0,rounding_size=0.015",
            facecolor=sw_col, edgecolor="#6B7280", lw=0.35, zorder=5,
        ))
        ax_ov.text(x_leg + sw_size + 0.06, sy, sw_lbl, ha="left", va="center",
                   fontsize=7.3, color="#222222")

    cb_x0 = x_leg + 0.01
    cb_w = 0.125
    cb_y0 = 0.23
    cb_h = 1.18
    cb_y1 = cb_y0 + cb_h
    ax_ov.text(x_leg, cb_y1 + 0.11, "PAE (Å)", ha="left", va="bottom",
               fontsize=8.0, fontweight="bold", color="#1A1A1A")

    grad = np.linspace(0.0, 31.0, 256)[:, None]
    ax_ov.imshow(
        grad, cmap=PAE_CMAP, vmin=0.0, vmax=31.0, origin="upper",
        extent=(cb_x0, cb_x0 + cb_w, cb_y0, cb_y1),
        aspect="auto", interpolation="bilinear", zorder=4,
    )
    ax_ov.add_patch(patches.Rectangle(
        (cb_x0, cb_y0), cb_w, cb_h,
        facecolor="none", edgecolor="#222222", lw=0.55, zorder=5,
    ))
    for t_ang in (0, 10, 20, 30):
        ty = cb_y1 - (t_ang / 31.0) * cb_h
        ax_ov.plot([cb_x0 + cb_w, cb_x0 + cb_w + 0.035], [ty, ty],
                   color="#222222", lw=0.55, zorder=5)
        ax_ov.text(cb_x0 + cb_w + 0.065, ty, str(t_ang),
                   ha="left", va="center", fontsize=7.2, color="#222222")

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix, kwargs in ((".pdf", {}), (".svg", {}), (".png", {"dpi": args.dpi})):
        out_p = Path(f"{args.output_prefix}{suffix}")
        fig.savefig(out_p, facecolor="white", **kwargs)
        written.append(out_p)
    plt.close(fig)
    return written


def main():
    args = parse_args()
    with tempfile.TemporaryDirectory() as td:
        written = build_figure(args, Path(td))
    print(
        f"Step 4.4 plot_fig1b_with_structures completed: "
        f"dapex={args.dapex_color}, iptm={args.iptm_color}, outputs={len(written)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
