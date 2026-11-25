#!/usr/bin/env python3
"""Colour bars for Fig. S5B (sequence position), S5C (pLDDT) and S5D (PAE).

Panel B uses the ChimeraX 'rainbow' palette (N terminus blue to C terminus red).
Panel C uses the four pLDDT bands set in chimrax_commands.txt:
    pLDDT >= 90  #2767CE      70-90  #55ABCC
    50-70        #E4C412      < 50   #EE7642
The bands are drawn at equal width, highest confidence on the left, and are
labelled with the band boundaries alone.
Panel D uses the palette, range and ticks of 7.12_plot_af3_pae.py, so that one
bar serves all five matrices.
Each bar is written twice: 'white' for placing on the black panel background,
'dark' for placing on a white margin.
"""
import argparse
from pathlib import Path

RAINBOW = ["#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000"]
PLDDT_COLOURS = ["#2767CE", "#55ABCC", "#E4C412", "#EE7642"]
PLDDT_BOUNDS = ["100", "90", "70", "50", "0"]
PAE = ["#087F3A", "#2FA24F", "#74C476", "#C7E9C0", "#F7FCF5"]
PAE_VMAX, PAE_STEP = 31.0, 5
W, BAR_H, FS = 260.0, 17.0, 20.0
GAP, FONT = 4.0, "Arial, Helvetica, sans-serif"
# Room for descenders, and a margin wide enough that a caption still fits when
# the renderer substitutes a wider font for Arial.
PAD = 0.35 * FS


def text_width(text, font_size):
    """Generous advance width of a label, used only to size the canvas."""
    return 0.62 * font_size * len(text) + font_size


def header(w, h):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.4g}pt" '
            f'height="{h:.4g}pt" viewBox="0 0 {w:.4g} {h:.4g}">')


def canvas(caption):
    """Canvas wide enough for the caption, with the bar centred in it."""
    width = max(W, text_width(caption, FS))
    return width, (width - W) / 2.0


def label(x, y, text, ink, anchor="middle"):
    return (f'<text x="{x:.4g}" y="{y:.4g}" text-anchor="{anchor}" '
            f'font-family="{FONT}" font-size="{FS:.4g}" fill="{ink}">{text}</text>')


def sequence_bar(ink):
    caption = "sequence position"
    width, x0 = canvas(caption)
    y_caption = BAR_H + GAP + 2 * FS + GAP / 2
    h = y_caption + PAD
    stops = "".join(f'<stop offset="{i/(len(RAINBOW)-1)*100:.4g}%" stop-color="{c}"/>'
                    for i, c in enumerate(RAINBOW))
    y_tick = BAR_H + GAP + FS
    return (header(width, h) +
            f'<defs><linearGradient id="seq" x1="0" x2="1">{stops}</linearGradient></defs>'
            f'<rect x="{x0:.4g}" y="0" width="{W:.4g}" height="{BAR_H:.4g}" fill="url(#seq)"/>'
            f'<rect x="{x0:.4g}" y="0" width="{W:.4g}" height="{BAR_H:.4g}" fill="none" '
            f'stroke="{ink}" stroke-width="1"/>' +
            label(x0, y_tick, "N", ink, "start") +
            label(x0 + W, y_tick, "C", ink, "end") +
            label(width / 2, y_caption, caption, ink) +
            '</svg>')


def plddt_bar(ink):
    """Equal-width bands, highest confidence on the left, boundary labels only."""
    caption = "pLDDT"
    width, x0 = canvas(caption)
    y_caption = BAR_H + GAP + 2 * FS + GAP / 2
    h = y_caption + PAD
    seg = W / len(PLDDT_COLOURS)
    y_tick = BAR_H + GAP + FS
    parts = [header(width, h)]
    for i, colour in enumerate(PLDDT_COLOURS):
        parts.append(f'<rect x="{x0 + i*seg:.4g}" y="0" width="{seg:.4g}" '
                     f'height="{BAR_H:.4g}" fill="{colour}"/>')
    last = len(PLDDT_BOUNDS) - 1
    for i, bound in enumerate(PLDDT_BOUNDS):
        anchor = "start" if i == 0 else ("end" if i == last else "middle")
        parts.append(label(x0 + i * seg, y_tick, bound, ink, anchor))
    parts.append(label(width / 2, y_caption, caption, ink))
    return "".join(parts) + "</svg>"


def pae_bar(ink):
    """Matches vmin=0, vmax=31 and ticks every 5 of 7.12_plot_af3_pae.py."""
    caption = "Expected Position Error (Å)"
    width, x0 = canvas(caption)
    tick_h = 4.0
    y_caption = BAR_H + tick_h + GAP + 2 * FS + GAP / 2
    h = y_caption + PAD
    stops = "".join(f'<stop offset="{i/(len(PAE)-1)*100:.4g}%" stop-color="{c}"/>'
                    for i, c in enumerate(PAE))
    y_tick = BAR_H + tick_h + GAP + FS
    parts = [header(width, h),
             f'<defs><linearGradient id="pae" x1="0" x2="1">{stops}</linearGradient></defs>',
             f'<rect x="{x0:.4g}" y="0" width="{W:.4g}" height="{BAR_H:.4g}" fill="url(#pae)"/>',
             f'<rect x="{x0:.4g}" y="0" width="{W:.4g}" height="{BAR_H:.4g}" fill="none" '
             f'stroke="{ink}" stroke-width="1"/>']
    ticks = [t for t in range(0, int(PAE_VMAX) + 1, PAE_STEP) if t <= PAE_VMAX]
    for t in ticks:
        x = x0 + t / PAE_VMAX * W
        parts.append(f'<line x1="{x:.4g}" y1="{BAR_H:.4g}" x2="{x:.4g}" '
                     f'y2="{BAR_H + tick_h:.4g}" stroke="{ink}" stroke-width="1"/>')
        anchor = "start" if t == ticks[0] else ("end" if t == ticks[-1] else "middle")
        parts.append(label(x, y_tick, str(t), ink, anchor))
    parts.append(label(width / 2, y_caption, caption, ink))
    return "".join(parts) + "</svg>"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Where to write the colour bars.")
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for tag, ink in (("white", "#FFFFFF"), ("dark", "#232323")):
        (out / f"colorbar_S5B_sequence_{tag}.svg").write_text(sequence_bar(ink))
        (out / f"colorbar_S5C_plddt_{tag}.svg").write_text(plddt_bar(ink))
        (out / f"colorbar_S5D_pae_{tag}.svg").write_text(pae_bar(ink))
        written += 3
    print(f"Step 7.13 make_colorbars completed: files={written}, output_dir={out}",
          flush=True)


if __name__ == "__main__":
    main()
