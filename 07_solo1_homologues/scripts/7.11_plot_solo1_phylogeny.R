#!/usr/bin/env Rscript
#' Step 7.11 script: plot the SOLO1 clade phylogeny and domain architecture for Fig. S5A.
#'
#' Draws the midpoint-rooted maximum-likelihood phylogeny of the 86 full-length
#' SOLO1 homologues alongside stacked domain-length bars (CC, NB-ARC, LRR, and
#' C-terminal region), printed C-terminal lengths, grass-lineage branch colours,
#' ultrafast bootstrap triangles (> threshold), integrated-domain annotations
#' (HMA, PB1, DUF761, Nudix), and dashed connectors to the five homologues
#' modelled in Fig. S5, B to D.

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
  library(tibble)
  library(svglite)
  library(ragg)
  library(systemfonts)
})

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--tree", type = "character",
    default = "07_solo1_homologues/data/results/7.6_solo1_clade.treefile",
    help = "Newick maximum-likelihood tree file [default %default]"),
  make_option("--domains", type = "character",
    default = "07_solo1_homologues/data/results/7.7_domain_length_summary.tsv",
    help = "Domain-length summary TSV from Step 7.7 [default %default]"),
  make_option("--clades", type = "character",
    default = "07_solo1_homologues/data/results/7.8_seqid_to_bop_pacmad.tsv",
    help = "Grass-clade assignment TSV from Step 7.8 [default %default]"),
  make_option("--labels", type = "character",
    default = "07_solo1_homologues/data/results/7.9_itol_labels.tsv",
    help = "Tip display-label TSV from Step 7.9 [default %default]"),
  make_option("--output-prefix", type = "character",
    default = "07_solo1_homologues/data/results/7.11_solo1_phylogeny",
    help = "Output prefix; writes .pdf, .svg, .png, and _plotted_values.tsv [default %default]"),
  make_option("--lrr-source", type = "character", default = "anchor",
    help = "Source for LRR bar length: 'anchor' (lrr_length_aa) or 'nlrtracker' (lrr_length_domain_aa) [default %default]"),
  make_option("--font-family", type = "character", default = "Helvetica",
    help = "Explicit sans-serif font family name [default %default]"),
  make_option("--bootstrap-threshold", type = "double", default = 80,
    help = "Bootstrap threshold above which internal nodes are marked [default %default]"),
  make_option("--show-model-connectors", type = "character", default = "FALSE",
    help = "Draw dashed connectors for the five homologues in Fig. S5, B to D (TRUE/FALSE) [default %default]"),
  make_option("--width", type = "double", default = 5.4,
    help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 8.4,
    help = "Figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 600,
    help = "PNG resolution in DPI [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_bool <- function(x, name) {
  val <- toupper(trimws(as.character(x)))
  if (val %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (val %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop(sprintf("Invalid boolean for --%s: %s", name, x), call. = FALSE)
}

lrr_source <- tolower(trimws(opt$`lrr-source`))
if (!lrr_source %in% c("anchor", "nlrtracker")) {
  stop("--lrr-source must be 'anchor' or 'nlrtracker'", call. = FALSE)
}
show_connectors <- parse_bool(opt$`show-model-connectors`, "show-model-connectors")
bs_thresh       <- opt$`bootstrap-threshold`

pick_sans_font <- function(preferred) {
  avail <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
  cands <- unique(c(preferred, "Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"))
  for (cand in cands) {
    if (cand %in% avail) {
      if (!identical(cand, preferred)) {
        message(sprintf("Requested font '%s' not installed; using '%s'", preferred, cand))
      }
      return(cand)
    }
  }
  warning(sprintf("None of (%s) are installed; falling back to 'sans'", paste(cands, collapse = ", ")),
          call. = FALSE)
  "sans"
}
sys_font <- pick_sans_font(opt$`font-family`)

# ---------------------------------------------------------------------------
# Midpoint rooting using ape (verified against exact longest tip-to-tip path)
# ---------------------------------------------------------------------------
midpoint_root_ape <- function(phy) {
  phy <- ape::unroot(phy)
  D <- ape::cophenetic.phylo(phy)
  idx <- which(D == max(D), arr.ind = TRUE)[1L, ]
  t1 <- idx[1L]
  t2 <- idx[2L]
  half <- D[t1, t2] / 2.0
  path <- ape::nodepath(phy, t1, t2)
  cum <- 0.0
  for (i in seq_len(length(path) - 1L)) {
    a <- path[i]
    b <- path[i + 1L]
    eidx <- which((phy$edge[, 1L] == a & phy$edge[, 2L] == b) |
                  (phy$edge[, 1L] == b & phy$edge[, 2L] == a))[1L]
    el <- phy$edge.length[eidx]
    if (cum + el >= half - 1e-12) {
      u <- phy$edge[eidx, 1L]
      v <- phy$edge[eidx, 2L]
      pos_from_u <- if (a == u) (half - cum) else (el - (half - cum))
      pos_from_u <- max(0.0, min(el, pos_from_u))

      tr_r <- ape::root(phy, node = v, resolve.root = TRUE, edgelabel = TRUE)
      root_nd <- ape::Ntip(tr_r) + 1L
      r_edges <- which(tr_r$edge[, 1L] == root_nd)
      if (length(r_edges) != 2L) {
        stop("Expected bifurcating root after midpoint rooting", call. = FALSE)
      }
      v_tips <- if (v <= ape::Ntip(phy)) phy$tip.label[v] else ape::extract.clade(phy, v)$tip.label
      c1 <- tr_r$edge[r_edges[1L], 2L]
      c1_tips <- if (c1 <= ape::Ntip(tr_r)) tr_r$tip.label[c1] else ape::extract.clade(tr_r, c1)$tip.label
      if (setequal(c1_tips, v_tips)) {
        tr_r$edge.length[r_edges[1L]] <- el - pos_from_u
        tr_r$edge.length[r_edges[2L]] <- pos_from_u
      } else {
        tr_r$edge.length[r_edges[2L]] <- el - pos_from_u
        tr_r$edge.length[r_edges[1L]] <- pos_from_u
      }
      return(tr_r)
    }
    cum <- cum + el
  }
  stop("Failed to locate midpoint along longest path", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Load & verify inputs (read-only)
# ---------------------------------------------------------------------------
for (p in c(opt$tree, opt$domains, opt$clades, opt$labels)) {
  if (!file.exists(p)) stop(sprintf("Required input file not found: %s", p), call. = FALSE)
}

tr_raw <- ape::read.tree(opt$tree)
if (ape::Ntip(tr_raw) != 86L) {
  stop(sprintf("Expected 86 tips in tree, got %d", ape::Ntip(tr_raw)), call. = FALSE)
}

tr_mid <- midpoint_root_ape(tr_raw)

# Verify pairwise distances and bootstrap support values are strictly preserved
D_orig <- ape::cophenetic.phylo(tr_raw)
D_mid  <- ape::cophenetic.phylo(tr_mid)[rownames(D_orig), colnames(D_orig)]
cophen_diff <- max(abs(D_orig - D_mid))
if (cophen_diff > 1e-11) {
  stop(sprintf("Midpoint rooting altered pairwise distances by %.3e", cophen_diff), call. = FALSE)
}

root_nd <- ape::Ntip(tr_mid) + 1L
root_kids <- tr_mid$edge[tr_mid$edge[, 1L] == root_nd, 2L]
depths_all <- ape::node.depth.edgelength(tr_mid)[seq_len(ape::Ntip(tr_mid))]
names(depths_all) <- tr_mid$tip.label
d_sides <- vapply(root_kids, function(k) {
  tips <- if (k <= ape::Ntip(tr_mid)) tr_mid$tip.label[k] else ape::extract.clade(tr_mid, k)$tip.label
  max(depths_all[tips])
}, numeric(1))
if (abs(diff(d_sides)) > 1e-11) {
  stop(sprintf("Root-to-tip max depths differ across root sides: %.6f vs %.6f", d_sides[1], d_sides[2]), call. = FALSE)
}

bs_orig <- sort(suppressWarnings(as.numeric(tr_raw$node.label)))
bs_mid  <- sort(suppressWarnings(as.numeric(tr_mid$node.label)))
if (!identical(bs_orig, bs_mid)) {
  stop("Bootstrap node labels changed during rooting!", call. = FALSE)
}

# Ladderize with right = TRUE so the 23-tip PACMAD/Oryza clade sits at top (y = 64..86)
tr_lad <- ape::ladderize(tr_mid, right = TRUE)

dom_df <- readr::read_tsv(
  opt$domains,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)
clade_df <- readr::read_tsv(
  opt$clades,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)
lbl_df <- readr::read_tsv(
  opt$labels,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)

for (nm in c("dom_df", "clade_df", "lbl_df")) {
  df_chk <- get(nm)
  if (nrow(df_chk) != 86L || !setequal(tr_lad$tip.label, df_chk$seq_id)) {
    stop(sprintf("%s must match the 86 tree tips 1-to-1", nm), call. = FALSE)
  }
}

# Integrated domain and Fig. S5 B-D modelled accessions
ID_TABLE <- tibble::tribble(
  ~protein_accession, ~integrated_domain,
  "XP_002449692.2",   "HMA",
  "XP_051188689.1",   "HMA",
  "XP_051188856.1",   "HMA",
  "XP_044964102.1",   "PB1",
  "XP_048559952.1",   "DUF761",
  "XP_020170543.1",   "DUF761",
  "XP_044335242.1",   "DUF761",
  "XP_048547975.1",   "Nudix",
  "XP_048554414.1",   "Nudix"
)
MODELED_ACCESSIONS <- c(
  "XP_020170543.1",
  "XP_044964102.1",
  "XP_051188689.1",
  "XP_002449692.2",
  "XP_048547975.1"
)

# ---------------------------------------------------------------------------
# Color palettes (exact Step 7.7 / 7.8 / Fig. S5A harmonies)
# ---------------------------------------------------------------------------
CLADE_COLORS <- c(
  "PACMAD"             = "#CC79A7",
  "BOP:Triticodae"     = "#0072B2",
  "BOP:Poodae"         = "#E69F00",
  "BOP:other"          = "#009E73",
  "Poaceae_unresolved" = "#FF7F00",
  "Non-Poaceae"        = "#33A02C",
  "Unresolved"         = "#BDBDBD",
  "Mixed"              = "#333333"
)

DOMAIN_COLORS <- c(
  "CC"                = "#3B6FB6",
  "NB-ARC"            = "#D68C45",
  "LRR"               = "#4A9A8F",
  "C-terminal region" = "#7E6BA8"
)

ID_COLORS <- c(
  "HMA"    = "#C83E4D",
  "PB1"    = "#C94B8B",
  "DUF761" = "#2E8B57",
  "Nudix"  = "#0097A7"
)

BS_TRIANGLE_COL <- "#E52421"

# ---------------------------------------------------------------------------
# Unanimous descendant consensus grass-lineage branch colouring
# ---------------------------------------------------------------------------
n_tip <- ape::Ntip(tr_lad)
n_tot <- n_tip + ape::Nnode(tr_lad)

clade_idx <- match(tr_lad$tip.label, clade_df$seq_id)
tip_clades <- clade_df$clade_detail[clade_idx]
names(tip_clades) <- tr_lad$tip.label

node_clade <- character(n_tot)
node_clade[seq_len(n_tip)] <- tip_clades
for (nd in (n_tip + 1L):n_tot) {
  d_tips <- ape::extract.clade(tr_lad, nd)$tip.label
  u <- unique(tip_clades[d_tips])
  node_clade[nd] <- if (length(u) == 1L) u else "Mixed"
}

p_tree_base <- ggtree(tr_lad, ladderize = FALSE, linewidth = 0.42, lineend = "square")
p_tree_base$data$branch_col <- unname(CLADE_COLORS[node_clade[p_tree_base$data$node]])
p_tree_base$mapping <- utils::modifyList(p_tree_base$mapping, aes(color = branch_col))

tree_data <- p_tree_base$data

tip_df <- tree_data %>%
  filter(isTip) %>%
  mutate(y_pos = y, tip_depth = x) %>%
  arrange(y_pos) %>%
  left_join(dom_df, by = c("label" = "seq_id")) %>%
  left_join(clade_df %>% select(seq_id, species, clade_coarse, clade_detail, color_hex),
            by = c("label" = "seq_id")) %>%
  left_join(lbl_df, by = c("label" = "seq_id")) %>%
  left_join(ID_TABLE, by = "protein_accession") %>%
  mutate(
    protein_length_aa   = as.integer(protein_length_aa),
    cc_length_aa        = as.integer(cc_length_aa),
    nbarc_length_aa     = as.integer(nbarc_length_aa),
    lrr_anchor_aa       = as.integer(lrr_length_aa),
    lrr_nlrtracker_aa   = as.integer(lrr_length_domain_aa),
    lrr_plotted_aa      = if (lrr_source == "anchor") lrr_anchor_aa else lrr_nlrtracker_aa,
    cterm_length_aa     = as.integer(lrr_end_to_end_length_aa),
    total_bar_aa        = cc_length_aa + nbarc_length_aa + lrr_plotted_aa + cterm_length_aa,
    integrated_domain   = ifelse(is.na(integrated_domain), "", integrated_domain),
    modeled_in_s5_b_d   = protein_accession %in% MODELED_ACCESSIONS,
    sp_part             = sub(" [^ ]+$", "", display_label),
    acc_part            = sub("^.* ", "", display_label),
    expr_label          = ifelse(
      modeled_in_s5_b_d,
      sprintf('paste(bolditalic("%s"), bold(" %s"))', sp_part, acc_part),
      sprintf('paste(italic("%s"), " %s")', sp_part, acc_part)
    )
  )

# Confirm 5 modelled proteins sit monotonically top-to-bottom so connectors never cross
modeled_check <- tip_df %>% filter(modeled_in_s5_b_d) %>% arrange(desc(y_pos))
if (!identical(modeled_check$protein_accession, MODELED_ACCESSIONS[c(4, 3, 2, 1, 5)])) {
  stop("Ladderized tree does not place the 5 modelled proteins in expected non-crossing order!", call. = FALSE)
}

# Check bootstrap internal nodes
int_nodes <- tree_data %>%
  filter(!isTip) %>%
  mutate(bs = suppressWarnings(as.numeric(label))) %>%
  filter(!is.na(bs) & bs > bs_thresh)
n_bs <- nrow(int_nodes)

# ---------------------------------------------------------------------------
# Shared vertical bounds across tree & domain-bar panels
# ---------------------------------------------------------------------------
Y_LIMS    <- c(0.2, 86.8)
MAX_DEPTH <- max(tip_df$tip_depth)
TOP_MARG_PT    <- 38
BOTTOM_MARG_PT <- 22

# ---------------------------------------------------------------------------
# Left Panel: Midpoint-rooted Phylogeny (ggtree)
# ---------------------------------------------------------------------------
make_left_triangles <- function(df, dx, dy) {
  if (nrow(df) == 0L) {
    return(tibble::tibble(id = integer(0), x = numeric(0), y = numeric(0)))
  }
  n <- nrow(df)
  tibble::tibble(
    id = rep(seq_len(n), each = 3L),
    x  = as.numeric(rbind(df$x - 0.55 * dx, df$x + 0.55 * dx, df$x + 0.55 * dx)),
    y  = as.numeric(rbind(df$y,             df$y + dy,        df$y - dy))
  )
}

leader_end <- MAX_DEPTH * 1.035
label_x0   <- MAX_DEPTH * 1.055
tree_xmax  <- MAX_DEPTH * 2.30

if (any(is.na(int_nodes$branch))) {
  stop("Missing branch midpoint coordinates in ggtree data", call. = FALSE)
}
bs_poly_df <- make_left_triangles(int_nodes %>% mutate(x = branch), dx = 0.090, dy = 0.42)

scale_x0 <- 0.04 * MAX_DEPTH
scale_x1 <- scale_x0 + 0.5

p_tree <- p_tree_base +
  # Subtle dotted grey leader lines aligning every tip to the right label column
  geom_segment(
    data = tip_df,
    aes(x = tip_depth, xend = leader_end, y = y_pos, yend = y_pos),
    inherit.aes = FALSE,
    color = "#CBD5E1",
    linewidth = 0.22,
    linetype = "dotted"
  ) +
  # Left-pointing ultrafast bootstrap > 80 red triangles on internal nodes
  geom_polygon(
    data = bs_poly_df,
    aes(x = x, y = y, group = id),
    inherit.aes = FALSE,
    fill = BS_TRIANGLE_COL,
    color = NA
  ) +
  # Tip labels: italic genus species + upright accession (bold for 5 modelled homologues)
  geom_text(
    data = tip_df,
    aes(x = label_x0, y = y_pos, label = expr_label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = 0,
    vjust = 0.5,
    size = 1.95,
    family = sys_font,
    color = "#1F2937"
  ) +
  scale_color_identity() +
  # Branch-length scale bar (0.5 substitutions/site) cleanly above tree top
  annotate("segment", x = scale_x0, xend = scale_x1, y = 87.6, yend = 87.6,
           color = "#1F2937", linewidth = 0.55) +
  annotate("text", x = (scale_x0 + scale_x1) / 2, y = 88.8, label = "0.5",
           size = 3.80, family = sys_font, color = "#1F2937", hjust = 0.5, vjust = 0) +
  coord_cartesian(xlim = c(-0.02 * MAX_DEPTH, tree_xmax), ylim = Y_LIMS, expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = TOP_MARG_PT, r = 0, b = BOTTOM_MARG_PT, l = 6, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Right Panel: Stacked Domain-Length Bars + Numbers + IDs + Model Connectors
# ---------------------------------------------------------------------------
bar_h <- 0.84
bar_df <- bind_rows(
  tip_df %>% transmute(y_pos, xmin = 0, xmax = cc_length_aa,
                       dom = "CC", fill = DOMAIN_COLORS["CC"], val = cc_length_aa),
  tip_df %>% transmute(y_pos, xmin = cc_length_aa, xmax = cc_length_aa + nbarc_length_aa,
                       dom = "NB-ARC", fill = DOMAIN_COLORS["NB-ARC"], val = nbarc_length_aa),
  tip_df %>% transmute(y_pos, xmin = cc_length_aa + nbarc_length_aa,
                       xmax = cc_length_aa + nbarc_length_aa + lrr_plotted_aa,
                       dom = "LRR", fill = DOMAIN_COLORS["LRR"], val = lrr_plotted_aa),
  tip_df %>% transmute(y_pos, xmin = cc_length_aa + nbarc_length_aa + lrr_plotted_aa,
                       xmax = total_bar_aa,
                       dom = "C-terminal region", fill = DOMAIN_COLORS["C-terminal region"], val = cterm_length_aa)
) %>%
  filter(xmax > xmin)

in_bar_df <- bar_df %>%
  filter(dom %in% c("CC", "NB-ARC", "LRR") & (xmax - xmin) >= 90) %>%
  mutate(xmid = 0.5 * (xmin + xmax), lbl = as.character(val))

max_total_aa <- max(tip_df$total_bar_aa)
x_id_bar0  <- max_total_aa + 175
x_id_bar1  <- x_id_bar0 + 48
x_id_lbl   <- x_id_bar1 + 32
x_conn_dot <- x_id_lbl + 320
x_conn_end <- x_conn_dot + 260
bar_xmax   <- if (show_connectors) x_conn_end + 360 else x_id_lbl + 560

# Group contiguous ID runs in y-order (1 = bottom, 86 = top)
tip_ordered <- tip_df %>% arrange(y_pos)
id_runs <- list()
cur_dom <- ""
cur_y0 <- NA_real_
cur_y1 <- NA_real_
for (i in seq_len(nrow(tip_ordered))) {
  d <- tip_ordered$integrated_domain[i]
  y <- tip_ordered$y_pos[i]
  if (nzchar(d)) {
    if (d == cur_dom && abs(y - cur_y1 - 1.0) < 1e-6) {
      cur_y1 <- y
    } else {
      if (nzchar(cur_dom)) {
        id_runs[[length(id_runs) + 1L]] <- tibble::tibble(
          domain = cur_dom, ymin = cur_y0 - 0.40, ymax = cur_y1 + 0.40, ymid = 0.5 * (cur_y0 + cur_y1)
        )
      }
      cur_dom <- d
      cur_y0 <- y
      cur_y1 <- y
    }
  } else if (nzchar(cur_dom)) {
    id_runs[[length(id_runs) + 1L]] <- tibble::tibble(
      domain = cur_dom, ymin = cur_y0 - 0.40, ymax = cur_y1 + 0.40, ymid = 0.5 * (cur_y0 + cur_y1)
    )
    cur_dom <- ""
  }
}
if (nzchar(cur_dom)) {
  id_runs[[length(id_runs) + 1L]] <- tibble::tibble(
    domain = cur_dom, ymin = cur_y0 - 0.40, ymax = cur_y1 + 0.40, ymid = 0.5 * (cur_y0 + cur_y1)
  )
}
id_runs_df <- bind_rows(id_runs) %>%
  mutate(col = unname(ID_COLORS[domain]))

# Domain column header midpoints derived from dataset medians
med_cc    <- stats::median(tip_df$cc_length_aa)
med_nb    <- stats::median(tip_df$nbarc_length_aa)
med_lrr   <- stats::median(tip_df$lrr_plotted_aa)
hx_cc     <- 0.5 * med_cc - 8
hx_nb     <- med_cc + 0.5 * med_nb + 15
hx_lrr    <- med_cc + med_nb + 0.5 * med_lrr + 10
hx_cterm  <- med_cc + med_nb + med_lrr + 60

x_ticks <- seq(0, floor(max_total_aa / 500) * 500, by = 500)

p_bars <- ggplot() +
  geom_rect(
    data = bar_df,
    aes(xmin = xmin, xmax = xmax, ymin = y_pos - bar_h / 2, ymax = y_pos + bar_h / 2, fill = fill),
    color = NA
  ) +
  scale_fill_identity() +
  geom_text(
    data = in_bar_df,
    aes(x = xmid, y = y_pos, label = lbl),
    size = 1.45,
    family = sys_font,
    color = "white",
    hjust = 0.5,
    vjust = 0.5
  ) +
  # Printed C-terminal length beside every bar
  geom_text(
    data = tip_df,
    aes(x = total_bar_aa + 22, y = y_pos, label = as.character(cterm_length_aa)),
    size = 1.45,
    family = sys_font,
    color = "#1F2937",
    hjust = 0,
    vjust = 0.5
  ) +
  # Integrated-domain vertical marker bars & bold labels (strictly aligned column)
  geom_rect(
    data = id_runs_df,
    aes(xmin = x_id_bar0, xmax = x_id_bar1, ymin = ymin, ymax = ymax, fill = col),
    color = NA
  ) +
  geom_text(
    data = id_runs_df,
    aes(x = x_id_lbl, y = ymid, label = domain, color = col),
    size = 3.70,
    fontface = "bold",
    family = sys_font,
    hjust = 0,
    vjust = 0.5
  ) +
  scale_color_identity() +
  # Domain headers at top
  annotate("text", x = 0.5 * max_total_aa, y = 91.2, label = "Domain size (aa)",
           size = 3.80, fontface = "bold", family = sys_font, color = "#111827", hjust = 0.5, vjust = 0) +
  annotate("text", x = c(hx_cc, hx_nb, hx_lrr, hx_cterm), y = 87.8,
           label = c("CC", "NB-ARC", "LRR", "C-terminal region"),
           color = unname(DOMAIN_COLORS[c("CC", "NB-ARC", "LRR", "C-terminal region")]),
           size = 2.30, fontface = "bold", family = sys_font,
           hjust = c(0.5, 0.5, 0.5, 0.0), vjust = 0) +
  # Clean bottom length scale
  annotate("segment", x = 0, xend = max(x_ticks), y = 0.2, yend = 0.2,
           color = "#333333", linewidth = 0.35) +
  annotate("segment", x = x_ticks, xend = x_ticks, y = 0.2, yend = -0.35,
           color = "#333333", linewidth = 0.32) +
  annotate("text", x = x_ticks, y = -0.65, label = as.character(x_ticks),
           size = 3.00, family = sys_font, color = "#333333", hjust = 0.5, vjust = 1) +
  annotate("text", x = max(x_ticks) + 190, y = -0.65, label = "aa",
           size = 3.00, family = sys_font, color = "#333333", hjust = 0, vjust = 1)

if (show_connectors) {
  conn_df <- modeled_check %>%
    arrange(desc(y_pos)) %>%
    mutate(
      y_target = seq(77, 11, length.out = 5),
      x_mid    = x_conn_dot + c(34, 68, 102, 68, 34)
    )
  p_bars <- p_bars +
    # Short dotted trace between printed C-terminal number and ID bar (never crossing ID text)
    geom_segment(
      data = conn_df,
      aes(x = total_bar_aa + 62, xend = x_id_bar0 - 10, y = y_pos, yend = y_pos),
      color = "#CBD5E1", linewidth = 0.28, linetype = "dotted"
    ) +
    geom_point(
      data = conn_df,
      aes(x = x_conn_dot, y = y_pos),
      shape = 16, size = 1.15, color = "#1F2937"
    ) +
    geom_segment(
      data = conn_df,
      aes(x = x_conn_dot, xend = x_mid, y = y_pos, yend = y_pos),
      color = "#333333", linewidth = 0.42, linetype = "dashed"
    ) +
    geom_segment(
      data = conn_df,
      aes(x = x_mid, xend = x_mid, y = y_pos, yend = y_target),
      color = "#333333", linewidth = 0.42, linetype = "dashed"
    ) +
    geom_segment(
      data = conn_df,
      aes(x = x_mid, xend = x_conn_end, y = y_target, yend = y_target),
      color = "#333333", linewidth = 0.42, linetype = "dashed"
    ) +
    annotate("text", x = x_conn_end + 10, y = conn_df$y_target,
             label = conn_df$protein_accession,
             size = 1.75, fontface = "bold", family = sys_font, color = "#333333", hjust = 0, vjust = 0.5)
  bar_xmax <- x_conn_end + 225
}

p_bars <- p_bars +
  coord_cartesian(xlim = c(0, bar_xmax), ylim = Y_LIMS, expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = TOP_MARG_PT, r = 8, b = BOTTOM_MARG_PT, l = 0, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Bottom Legend: Frameless, count-free spacious 2x3 layout
# ---------------------------------------------------------------------------
leg_branches <- tibble::tribble(
  ~x,  ~y,   ~col,                            ~lbl,
  2,   0.70, CLADE_COLORS["PACMAD"],          "PACMAD clade",
  96,  0.70, CLADE_COLORS["BOP:Triticodae"],  "BOP clade (Triticodae)",
  2,   0.26, CLADE_COLORS["BOP:Poodae"],      "BOP clade (Poodae)",
  96,  0.26, CLADE_COLORS["BOP:other"],       "BOP clade (other)"
)

leg_tri_df <- make_left_triangles(tibble::tibble(x = 214, y = 0.70), dx = 7.0, dy = 0.125)

p_leg <- ggplot() +
  geom_segment(data = leg_branches, aes(x = x, xend = x + 13, y = y, yend = y, color = col),
               linewidth = 1.90, lineend = "round") +
  scale_color_identity() +
  geom_text(data = leg_branches, aes(x = x + 18, y = y, label = lbl),
            hjust = 0, vjust = 0.5, size = 3.15, family = sys_font, color = "#1F2937") +
  # Left-pointing bootstrap triangle key (row 1, col 3)
  geom_polygon(data = leg_tri_df, aes(x = x, y = y, group = id),
               fill = BS_TRIANGLE_COL, color = NA) +
  annotate("text", x = 224, y = 0.70, label = sprintf("Bootstrap > %g", bs_thresh),
           hjust = 0, vjust = 0.5, size = 3.15, family = sys_font, color = "#1F2937")

if (show_connectors) {
  p_leg <- p_leg +
    annotate("segment", x = 207, xend = 220, y = 0.26, yend = 0.26,
             color = "#333333", linewidth = 0.52, linetype = "dashed") +
    annotate("text", x = 224, y = 0.26, label = "Modelled in Fig. S5B-D",
             hjust = 0, vjust = 0.5, size = 3.15, family = sys_font, color = "#1F2937")
}

p_leg <- p_leg +
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = 2, r = 8, b = 4, l = 6, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Assemble with patchwork
# ---------------------------------------------------------------------------
top_row <- (p_tree | p_bars) + plot_layout(widths = c(1.42, 1.00))
combined <- top_row / patchwork::wrap_elements(full = p_leg) +
  plot_layout(heights = c(1.0, 0.088))

# ---------------------------------------------------------------------------
# SVG Sanitizer (for Affinity / Inkscape compatibility) & Exports
# ---------------------------------------------------------------------------
sanitize_svg_for_editors <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  out <- character(length(lines))
  k <- 0L
  in_clip <- FALSE

  inject_style <- function(line, prop) {
    m_pos <- regexpr("style='[^']*'", line)
    if (m_pos[1L] == -1L) return(line)
    m_len <- attr(m_pos, "match.length")
    style_attr <- substr(line, m_pos[1L], m_pos[1L] + m_len - 1L)
    inner <- trimws(substr(style_attr, 8L, nchar(style_attr) - 1L))
    if (nzchar(inner) && !grepl(";$", inner)) inner <- paste0(inner, ";")
    new_style <- sprintf("style='%s %s'", inner, prop)
    paste0(substr(line, 1L, m_pos[1L] - 1L), new_style, substr(line, m_pos[1L] + m_len, nchar(line)))
  }

  for (line in lines) {
    if (grepl("<clipPath\\b", line)) in_clip <- TRUE
    if (!in_clip) {
      if (grepl("^\\s*<rect\\b", line) && grepl("stroke:\\s*none", line) && !grepl("fill:", line)) {
        next
      }
      if (grepl("^\\s*<(rect|circle|polygon|path|polyline)\\b", line) && !grepl("fill:", line)) {
        line <- inject_style(line, "fill: none;")
      }
      if (grepl("^\\s*<(line|polyline|path)\\b", line) && !grepl("stroke-linecap:", line)) {
        line <- inject_style(line, "stroke-linecap: round;")
      }
      if (grepl("^\\s*<text\\b", line) && !grepl("fill:", line)) {
        line <- inject_style(line, "fill: #000000;")
      }
    }
    if (grepl("</clipPath>", line)) in_clip <- FALSE
    k <- k + 1L
    out[k] <- line
  }
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(out[seq_len(k)], con, useBytes = TRUE)
}

out_prefix <- opt$`output-prefix`
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pdf_path <- paste0(out_prefix, ".pdf")
svg_path <- paste0(out_prefix, ".svg")
png_path <- paste0(out_prefix, ".png")
tsv_path <- paste0(out_prefix, "_plotted_values.tsv")

if (isTRUE(capabilities()[["aqua"]])) {
  if (identical(sys_font, "Helvetica Neue")) {
    grDevices::quartzFonts(
      "Helvetica Neue" = grDevices::quartzFont(c(
        "HelveticaNeue", "HelveticaNeue-Bold", "HelveticaNeue-Italic", "HelveticaNeue-BoldItalic"
      ))
    )
  }
  grDevices::quartz(type = "pdf", file = pdf_path, width = opt$width, height = opt$height, family = sys_font, bg = "white")
} else {
  grDevices::cairo_pdf(pdf_path, width = opt$width, height = opt$height, family = sys_font, bg = "white")
}
print(combined)
grDevices::dev.off()

svglite::svglite(
  svg_path,
  width = opt$width,
  height = opt$height,
  bg = "white",
  system_fonts = list(sans = sys_font),
  fix_text_size = FALSE
)
print(combined)
grDevices::dev.off()
sanitize_svg_for_editors(svg_path)

ragg::agg_png(png_path, width = opt$width, height = opt$height, units = "in", res = opt$dpi, background = "white")
print(combined)
grDevices::dev.off()

# Deterministic plotted-values table (sorted by seq_id, fixed columns, no timestamps)
plotted_tsv <- tip_df %>%
  arrange(label) %>%
  transmute(
    seq_id                   = label,
    protein_accession        = protein_accession,
    display_label            = display_label,
    tree_row                 = as.integer(round(y_pos)),
    tip_depth                = sprintf("%.6f", tip_depth),
    clade_coarse             = clade_coarse,
    clade_detail             = clade_detail,
    branch_color_hex         = color_hex,
    protein_length_aa        = protein_length_aa,
    cc_length_aa             = cc_length_aa,
    nbarc_length_aa          = nbarc_length_aa,
    lrr_length_plotted_aa    = lrr_plotted_aa,
    lrr_source               = lrr_source,
    lrr_length_anchor_aa     = lrr_anchor_aa,
    lrr_length_nlrtracker_aa = lrr_nlrtracker_aa,
    cterm_length_aa          = cterm_length_aa,
    integrated_domain        = integrated_domain,
    modeled_in_s5_b_to_d     = ifelse(modeled_in_s5_b_d, "TRUE", "FALSE")
  )
readr::write_tsv(plotted_tsv, tsv_path, na = "")

cat(sprintf(
  "Step 7.11 plot_solo1_phylogeny completed: tips=%d, bootstrap_gt%g=%d, lrr_source=%s, font=%s, outputs=4\n",
  n_tip, bs_thresh, n_bs, lrr_source, sys_font
), file = stderr())
