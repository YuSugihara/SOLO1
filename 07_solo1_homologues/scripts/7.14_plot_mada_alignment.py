#!/usr/bin/env python3
"""Step 7.14 script: align NbNRC4 and AetSOLO1 and plot the MADA-motif region.

Aligns full-length NbNRC4 (RefPlantNLR NbNRC4c = GenBank MK692737.1) and
AetSOLO1 (XP_020170543.1) with MAFFT L-INS-i (--localpair --maxiterate 1000)
and draws the BoxShade alignment of the N-terminal MADA motif with red frames
marking the three residues substituted in the SOLO1^EEE variants.
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.Align import substitution_matrices

COL_IDENT_BG = "#000000"
COL_SIM_BG   = "#888888"
COL_NON_BG   = "#FFFFFF"
COL_TEXT_BLK = "#000000"
COL_TEXT_WHT = "#FFFFFF"
COL_MUT_RED  = "#D90000"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refplantnlr", type=Path,
                   default=Path("06_plant_nlr_nbarc_tree/data/inputs/6.0_refplantnlr/RefPlantNLR_NLR.fasta"),
                   help="RefPlantNLR full-length protein FASTA containing NbNRC4c.")
    p.add_argument("--solo1-fasta", type=Path,
                   default=Path("07_solo1_homologues/data/results/7.1_solo1_clade_full_length.fasta"),
                   help="SOLO1 clade full-length protein FASTA containing XP_020170543.1.")
    p.add_argument("--mafft-bin", default="mafft",
                   help="MAFFT executable name or path (default: %(default)s).")
    p.add_argument("--output-alignment", type=Path,
                   default=Path("07_solo1_homologues/data/results/7.14_nrc4c_aetsolo1.mafft_linsi.fasta"),
                   help="Output path for the full-length MAFFT L-INS-i alignment FASTA.")
    p.add_argument("--output-prefix", type=Path,
                   default=Path("07_solo1_homologues/data/results/7.14_mada_alignment"),
                   help="Output prefix (writes .pdf, .svg, .png).")
    p.add_argument("--matrix", default="BLOSUM62",
                   help="Biopython substitution matrix name (default: %(default)s).")
    p.add_argument("--sim-threshold", type=float, default=0.0,
                   help="Substitution matrix score threshold (> threshold) for Similar shading.")
    p.add_argument("--n-cols", type=int, default=24,
                   help="Alignment columns covering NbNRC4 MADA motif 1-21 (default: %(default)s).")
    p.add_argument("--dpi", type=int, default=600,
                   help="Raster output DPI (default: %(default)s).")
    return p.parse_args()


def _register_macos_helvetica_neue(tmp_font_dir: Path):
    """Optionally extract macOS HelveticaNeue.ttc faces if present; ignore otherwise."""
    ttc_path = Path("/System/Library/Fonts/HelveticaNeue.ttc")
    if not ttc_path.exists():
        return
    try:
        from fontTools.ttLib import TTCollection
        ttc = TTCollection(str(ttc_path))
        for idx, fname in (
            (0, "HelveticaNeue-Regular.ttf"),
            (1, "HelveticaNeue-Bold.ttf"),
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
        "font.monospace": ["Menlo", "DejaVu Sans Mono", "Courier New", "monospace"],
        "mathtext.default": "regular",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def extract_sequences(ref_fasta: Path, solo1_fasta: Path):
    nrc4_seq = None
    for rec in SeqIO.parse(ref_fasta, "fasta"):
        if rec.id == "NbNRC4c":
            nrc4_seq = str(rec.seq)
            break
    if nrc4_seq is None:
        sys.exit(f"error: NbNRC4c not found in {ref_fasta}")

    aet_seq = None
    for rec in SeqIO.parse(solo1_fasta, "fasta"):
        if "XP_020170543.1" in rec.id:
            aet_seq = str(rec.seq)
            break
    if aet_seq is None:
        sys.exit(f"error: XP_020170543.1 not found in {solo1_fasta}")

    return nrc4_seq, aet_seq


def run_mafft_linsi(mafft_bin: str, nrc4_seq: str, aet_seq: str, out_fasta: Path):
    resolved = shutil.which(mafft_bin) or mafft_bin
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", encoding="utf-8") as tmp_in:
        tmp_in.write(f">NbNRC4c\n{nrc4_seq}\n>AetSOLO1_XP_020170543.1\n{aet_seq}\n")
        tmp_in.flush()
        aligned = subprocess.check_output(
            [resolved, "--localpair", "--maxiterate", "1000", "--quiet", tmp_in.name],
            text=True,
        )
    out_fasta.write_text(aligned, encoding="utf-8")
    recs = {r.id: str(r.seq) for r in SeqIO.parse(out_fasta, "fasta")}
    return recs["NbNRC4c"], recs["AetSOLO1_XP_020170543.1"]


def residue_positions(aln_seq: str):
    """Return 1-based ungapped residue index (or None for gap) per alignment column."""
    pos = 0
    out = []
    for ch in aln_seq:
        if ch == "-":
            out.append(None)
        else:
            pos += 1
            out.append(pos)
    return out


def is_similar(a: str, b: str, matrix, threshold: float) -> bool:
    if a == "-" or b == "-" or a == b:
        return False
    return float(matrix[a, b]) > threshold


def plot_alignment(aln_nrc4: str, aln_aet: str, args, tmp_font_dir: Path):
    configure_style(tmp_font_dir)
    matrix = substitution_matrices.load(args.matrix)
    n_cols = args.n_cols

    assert aln_nrc4[:n_cols] == "MADAVV---NFLVENLLQLLTDNV", (
        f"Unexpected NbNRC4c MADA alignment prefix: {aln_nrc4[:n_cols]}"
    )

    pos_nrc4 = residue_positions(aln_nrc4)
    pos_aet = residue_positions(aln_aet)

    col_for_nrc4 = {p: c for c, p in enumerate(pos_nrc4[:n_cols], start=1) if p is not None}
    col_for_aet = {p: c for c, p in enumerate(pos_aet[:n_cols], start=1) if p is not None}

    expected_pairs = [
        (9,  12, "L", "L", "L12E"),
        (13, 16, "L", "L", "L16E"),
        (17, 20, "L", "I", "I20E"),
    ]
    mut_cols = []
    for p_nrc, p_aet, r_nrc, r_aet, mut_lbl in expected_pairs:
        c_nrc = col_for_nrc4[p_nrc]
        c_aet = col_for_aet[p_aet]
        assert c_nrc == c_aet, (
            f"Column mismatch: NbNRC4 {p_nrc} at col {c_nrc} != AetSOLO1 {p_aet} at col {c_aet}"
        )
        assert aln_nrc4[c_nrc - 1] == r_nrc, (
            f"Expected NbNRC4 {r_nrc}{p_nrc}, got {aln_nrc4[c_nrc - 1]}"
        )
        assert aln_aet[c_aet - 1] == r_aet, (
            f"Expected AetSOLO1 {r_aet}{p_aet}, got {aln_aet[c_aet - 1]}"
        )
        mut_cols.append((c_nrc, f"{r_nrc}{p_nrc}", mut_lbl))

    end_res_nrc4 = max(p for p in pos_nrc4[:n_cols] if p is not None)
    end_res_aet = max(p for p in pos_aet[:n_cols] if p is not None)

    W, H = 4.95, 1.54
    fig = plt.figure(figsize=(W, H), dpi=args.dpi, facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1], facecolor="white")
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    x_grid0 = 1.00
    cell_w = 0.142
    x_grid1 = x_grid0 + n_cols * cell_w

    def col_x0(c):
        return x_grid0 + (c - 1) * cell_w

    def col_cx(c):
        return x_grid0 + (c - 0.5) * cell_w

    y_mada_lbl  = 1.42
    y_mada_line = 1.31
    y_ref_lbl   = 1.15
    cell_h      = 0.22
    y_aet_bot   = 0.59
    y_nrc4_bot  = y_aet_bot + cell_h
    y_mut_lbl   = 0.34
    y_key_y     = 0.10

    # Top MADA motif bracket (NbNRC4 residues 1-21 -> cols 1-24)
    ax.text(0.5 * (x_grid0 + x_grid1), y_mada_lbl, "MADA motif (NbNRC4)",
            ha="center", va="center", fontsize=8.2, color=COL_TEXT_BLK)
    tick_h = 0.035
    ax.plot([x_grid0, x_grid0, x_grid1, x_grid1],
            [y_mada_line - tick_h, y_mada_line, y_mada_line, y_mada_line - tick_h],
            color=COL_TEXT_BLK, lw=0.85, solid_capstyle="butt", zorder=3)

    # Sequence labels and ungapped bounds
    y_nrc4_cy = y_nrc4_bot + cell_h / 2.0
    y_aet_cy  = y_aet_bot + cell_h / 2.0

    ax.text(0.75, y_nrc4_cy, "NbNRC4", ha="right", va="center",
            fontsize=8.6, fontweight="bold", color=COL_TEXT_BLK)
    ax.text(0.75, y_aet_cy, "AetSOLO1", ha="right", va="center",
            fontsize=8.6, fontweight="bold", color=COL_TEXT_BLK)

    ax.text(x_grid0 - 0.06, y_nrc4_cy, "1", ha="right", va="center", fontsize=7.6, color=COL_TEXT_BLK)
    ax.text(x_grid1 + 0.06, y_nrc4_cy, str(end_res_nrc4), ha="left", va="center", fontsize=7.6, color=COL_TEXT_BLK)

    ax.text(x_grid0 - 0.06, y_aet_cy, "1", ha="right", va="center", fontsize=7.6, color=COL_TEXT_BLK)
    ax.text(x_grid1 + 0.06, y_aet_cy, str(end_res_aet), ha="left", va="center", fontsize=7.6, color=COL_TEXT_BLK)

    # BoxShade tight monospace blocks (Identical = #000, Similar > threshold = #888)
    for c in range(1, n_cols + 1):
        a = aln_nrc4[c - 1]
        b = aln_aet[c - 1]
        x0 = col_x0(c)
        xc = col_cx(c)
        ident = (a == b and a != "-")
        sim = is_similar(a, b, matrix, args.sim_threshold)

        if ident:
            bg, fg, fw = COL_IDENT_BG, COL_TEXT_WHT, "bold"
        elif sim:
            bg, fg, fw = COL_SIM_BG, COL_TEXT_WHT, "bold"
        else:
            bg, fg, fw = COL_NON_BG, COL_TEXT_BLK, "normal"

        for yb, ch in ((y_nrc4_bot, a), (y_aet_bot, b)):
            if bg != COL_NON_BG:
                ax.add_patch(patches.Rectangle(
                    (x0, yb), cell_w, cell_h,
                    facecolor=bg, edgecolor="none", zorder=2,
                ))
            ax.text(
                xc, yb + cell_h / 2.0, ch,
                ha="center", va="center",
                fontfamily="monospace", fontsize=8.8, fontweight=fw, color=fg, zorder=4,
            )

    # Red frames around EEE columns + mutation callouts
    for c, nrc_lbl, mut_lbl in mut_cols:
        x0 = col_x0(c)
        xc = col_cx(c)

        ax.add_patch(patches.Rectangle(
            (x0, y_aet_bot), cell_w, 2.0 * cell_h,
            facecolor="none", edgecolor=COL_MUT_RED, lw=1.35, zorder=6,
        ))

        ax.text(xc, y_ref_lbl, nrc_lbl, ha="center", va="center",
                fontsize=7.3, color=COL_TEXT_BLK, zorder=5)

        ax.annotate(
            "",
            xy=(xc, y_mut_lbl + 0.068),
            xytext=(xc, y_aet_bot - 0.015),
            arrowprops=dict(arrowstyle="-|>,head_length=0.22,head_width=0.14",
                            color=COL_MUT_RED, lw=0.90),
            zorder=6,
        )
        ax.text(xc, y_mut_lbl, mut_lbl, ha="center", va="center",
                fontsize=7.8, fontweight="bold", color=COL_MUT_RED, zorder=6)

    # Inline BoxShade key
    sw = 0.105
    thresh_str = f"{args.sim_threshold:g}"
    sim_lbl = f"Similar ({args.matrix} > {thresh_str})"
    key_items = [
        (x_grid0,        COL_IDENT_BG, "none",      0.0,  "Identical"),
        (x_grid0 + 0.78, COL_SIM_BG,   "none",      0.0,  sim_lbl),
        (x_grid0 + 2.26, "none",       COL_MUT_RED, 1.25, "EEE substitution"),
    ]
    for kx, bg_c, edge_c, lw_c, lbl in key_items:
        ax.add_patch(patches.Rectangle(
            (kx, y_key_y - sw / 2.0), sw, sw,
            facecolor=bg_c, edgecolor=edge_c, lw=lw_c, zorder=5,
        ))
        ax.text(kx + sw + 0.055, y_key_y, lbl,
                ha="left", va="center", fontsize=7.0, color=COL_TEXT_BLK, zorder=5)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": args.dpi})):
        out_p = Path(f"{args.output_prefix}{suffix}")
        fig.savefig(out_p, facecolor="white", **kw)
        written.append(out_p)
    plt.close(fig)
    return written, mut_cols


def main():
    args = parse_args()
    nrc4_seq, aet_seq = extract_sequences(args.refplantnlr, args.solo1_fasta)
    aln_nrc4, aln_aet = run_mafft_linsi(args.mafft_bin, nrc4_seq, aet_seq, args.output_alignment)
    with tempfile.TemporaryDirectory() as td:
        written, mut_cols = plot_alignment(aln_nrc4, aln_aet, args, Path(td))
    print(
        f"Step 7.14 plot_mada_alignment completed: cols={args.n_cols}, "
        f"matrix={args.matrix}>{args.sim_threshold:g}, "
        f"sites={[m[2] for m in mut_cols]}, outputs={len(written)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
