#!/usr/bin/env Rscript
#' Step 5.7 script: plot the candidate NB-ARC phylogeny for Fig. 2A.
#'
#' Draws the upward-facing midpoint-rooted phylogeny with plant-order branch
#' colouring, sensor-clade and SOLO1 background shading, cell-death screening
#' markers, and ipTM_LCB bar track using R, ggtree, and patchwork. Strictly
#' preserves exact Newick branch lengths, topology, bootstrap support values,
#' tip set (296), and structural/assay measurements.

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
    default = "05_candidate_nbarc_tree/data/results/5.5_candidate_nbarc.midpoint_rooted.treefile",
    help = "Midpoint-rooted Newick tree file [default %default]"),
  make_option("--annotations", type = "character",
    default = "05_candidate_nbarc_tree/data/results/5.6_candidate_tip_annotations.tsv",
    help = "Unified tip annotation table from Step 5.6 [default %default]"),
  make_option("--output-prefix", type = "character",
    default = "05_candidate_nbarc_tree/data/results/5.7_candidate_phylogeny",
    help = "Output path prefix without extension; writes .pdf, .svg, .png [default %default]"),
  make_option("--width", type = "double", default = 7.8,
    help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 3.32,
    help = "Figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 600,
    help = "PNG resolution in DPI [default %default]"),
  make_option("--iptm-threshold", type = "double", default = 0.5,
    help = "Threshold splitting ipTM bar colors [default %default]"),
  make_option("--iptm-color-mode", type = "character", default = "split",
    help = "ipTM bar coloring: 'split' (above vs below threshold) or 'single' (uniform navy) [default %default]"),
  make_option("--plant-order-mode", type = "character", default = "branch",
    help = "Plant order display: 'branch' (branch color), 'strip' (bar above tree), or 'none' [default %default]"),
  make_option("--show-clade-shades", type = "character", default = "TRUE",
    help = "Show background washes for sensor clades and SOLO1 (TRUE/FALSE) [default %default]"),
  make_option("--show-callouts", type = "character", default = "TRUE",
    help = "Show species/rank callout labels for Cell death == Yes (TRUE/FALSE) [default %default]"),
  make_option("--bootstrap-mode", type = "character", default = "major",
    help = "Bootstrap labels: 'major' (backbone & clade MRCAs), 'all95' (dots >=95), or 'none' [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

parse_bool <- function(x, name) {
  val <- toupper(trimws(as.character(x)))
  if (val %in% c("TRUE", "T", "1", "YES")) return(TRUE)
  if (val %in% c("FALSE", "F", "0", "NO")) return(FALSE)
  stop(sprintf("Invalid boolean for --%s: %s", name, x), call. = FALSE)
}

show_clade_shades <- parse_bool(opt$`show-clade-shades`, "show-clade-shades")
show_callouts     <- parse_bool(opt$`show-callouts`, "show-callouts")
iptm_color_mode   <- tolower(trimws(opt$`iptm-color-mode`))
order_mode        <- tolower(trimws(opt$`plant-order-mode`))
bs_mode           <- tolower(trimws(opt$`bootstrap-mode`))

if (!iptm_color_mode %in% c("split", "single")) {
  stop("--iptm-color-mode must be 'split' or 'single'", call. = FALSE)
}
if (!order_mode %in% c("branch", "strip", "none")) {
  stop("--plant-order-mode must be 'branch', 'strip', or 'none'", call. = FALSE)
}
if (!bs_mode %in% c("major", "all95", "none")) {
  stop("--bootstrap-mode must be 'major', 'all95', or 'none'", call. = FALSE)
}

# Pick clean sans-serif font available across macOS and Linux systems
pick_sans_font <- function() {
  for (cand in c("Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans")) {
    matched <- tryCatch(systemfonts::match_fonts(cand), error = function(e) NULL)
    if (!is.null(matched$path) && nzchar(matched$path)) {
      return(cand)
    }
  }
  "sans"
}
sys_font <- pick_sans_font()

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
tree_path <- opt$tree
ann_path  <- opt$annotations

for (p in c(tree_path, ann_path)) {
  if (!file.exists(p)) stop(sprintf("Required input file not found: %s", p), call. = FALSE)
}

tr <- ape::read.tree(tree_path)
if (ape::Ntip(tr) != 296L) {
  stop(sprintf("Expected 296 tips in tree, got %d", ape::Ntip(tr)), call. = FALSE)
}

# Record original root-to-tip distances & sorted edge lengths BEFORE ladderize
orig_depths <- setNames(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))], tr$tip.label)
orig_edges_sorted <- sort(tr$edge.length)

# Ladderize with right = FALSE so L->R order matches Fig. 2A (Pii-2 -> NRC-S -> SOLO1 -> Pia-2/MIC1)
tr_lad <- ape::ladderize(tr, right = FALSE)
if (max(abs(sort(tr_lad$edge.length) - orig_edges_sorted)) > 1e-12) {
  stop("Ladderization altered branch lengths!", call. = FALSE)
}

ann_df <- readr::read_tsv(
  ann_path,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)
req_cols <- c(
  "tip_label", "species", "assembly", "protein_ID",
  "plant_order", "plant_order_group", "sensor_clade", "is_solo1_clade",
  "weighted_rank", "ipTM_LCB", "D_APEX", "screening_stage", "Cell_death"
)
missing_cols <- setdiff(req_cols, colnames(ann_df))
if (length(missing_cols) > 0L) {
  stop(sprintf("Missing required columns in annotations TSV: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}
if (nrow(ann_df) != 296L || !setequal(tr_lad$tip.label, ann_df$tip_label)) {
  stop("Annotations TSV must have an exact 1-to-1 match with the 296 tree tip labels!", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Color palette & shared horizontal geometry
# ---------------------------------------------------------------------------
# Colorblind-safe branch palette verified across normal, deuteranopia,
# protanopia, and tritanopia (min pairwise CIEDE2000 >= 25.3 normal, >= 11.5 CVD)
BRANCH_ORDER_COLORS <- c(
  "Poales"         = "#0B5394",
  "Solanales"      = "#D87B00",
  "Lamiales"       = "#00947E",
  "Caryophyllales" = "#B8195A",
  "Other"          = "#969EA8",
  "Mixed"          = "#1E1E1E"
)

STRIP_ORDER_COLORS <- c(
  "Poales"         = "#0B5394",
  "Solanales"      = "#D87B00",
  "Lamiales"       = "#00947E",
  "Caryophyllales" = "#B8195A",
  "Other"          = "#D3D7DC"
)

IPTM_COLORS <- if (iptm_color_mode == "single") {
  c("gt" = "#15415C", "le" = "#15415C")
} else {
  c("gt" = "#15415C", "le" = "#CBD0D6")
}
CD_YES_COLOR  <- "#D0242A"
NOCD_COLOR    <- "#525252"
NOCD_LINE     <- "#94A3B8"

# Compute unanimous descendant consensus plant order per node on tr_lad
top4_orders <- c("Poales", "Solanales", "Lamiales", "Caryophyllales")
n_tip <- ape::Ntip(tr_lad)
n_tot <- n_tip + ape::Nnode(tr_lad)

ann_idx <- match(tr_lad$tip.label, ann_df$tip_label)
expected_group <- ifelse(ann_df$plant_order[ann_idx] %in% top4_orders, ann_df$plant_order[ann_idx], "Other")
if (!identical(as.character(ann_df$plant_order_group[ann_idx]), as.character(expected_group))) {
  stop("plant_order_group column does not match expected top-4 grouping of plant_order!", call. = FALSE)
}

tip_ord_disp <- ann_df$plant_order_group[ann_idx]
names(tip_ord_disp) <- tr_lad$tip.label

node_order <- character(n_tot)
node_order[seq_len(n_tip)] <- tip_ord_disp
for (nd in (n_tip + 1L):n_tot) {
  d_tips <- ape::extract.clade(tr_lad, nd)$tip.label
  u <- unique(tip_ord_disp[d_tips])
  node_order[nd] <- if (length(u) == 1L) u else "Mixed"
}

# ---------------------------------------------------------------------------
# Build base ggtree plot & verify exact coordinate/branch-length fidelity
# ---------------------------------------------------------------------------
if (order_mode == "branch") {
  p_tree_base <- ggtree(tr_lad, ladderize = FALSE, linewidth = 0.34, lineend = "square")
  p_tree_base$data$branch_col <- unname(BRANCH_ORDER_COLORS[node_order[p_tree_base$data$node]])
  p_tree_base$mapping <- utils::modifyList(p_tree_base$mapping, aes(color = branch_col))
} else {
  p_tree_base <- ggtree(tr_lad, ladderize = FALSE, color = "#232323", linewidth = 0.28, lineend = "square")
}

tree_data <- p_tree_base$data

tip_df <- tree_data %>%
  filter(isTip) %>%
  mutate(
    x_pos     = y,
    tip_depth = x
  ) %>%
  arrange(x_pos)

if (nrow(tip_df) != 296L) stop("Tip count mismatch in ggtree data", call. = FALSE)
depth_diff <- max(abs(tip_df$tip_depth - orig_depths[tip_df$label]))
if (depth_diff > 1e-9) {
  stop(sprintf("ggtree tip depths deviate from Newick root-to-tip distances by %.3e", depth_diff), call. = FALSE)
}

# Join tip annotations (autoactivity merged into Cell_death == Yes without distinction)
tip_df <- tip_df %>%
  left_join(ann_df, by = c("label" = "tip_label")) %>%
  mutate(
    species_id    = species,
    weighted_rank = as.integer(weighted_rank),
    ipTM_LCB      = as.numeric(ipTM_LCB),
    is_solo1      = is_solo1_clade == "True",
    order         = plant_order,
    order_disp    = factor(plant_order_group, levels = c(top4_orders, "Other")),
    iptm_group    = ifelse(ipTM_LCB > opt$`iptm-threshold`, "gt", "le"),
    is_tested_sr  = screening_stage == "structural ranking",
    is_rank12     = weighted_rank %in% c(1L, 2L)
  )

if (any(is.na(tip_df$weighted_rank)) || any(is.na(tip_df$ipTM_LCB)) || any(is.na(tip_df$order))) {
  stop("Missing score or plant order annotation for one or more tree tips!", call. = FALSE)
}

# Verify clade monophyly & Fig. 2A left-to-right ordering
check_clade_span <- function(df_sub, cname) {
  xs <- sort(df_sub$x_pos)
  if (length(xs) == 0 || length(xs) != (max(xs) - min(xs) + 1L)) {
    stop(sprintf("Clade '%s' is not contiguous/monophyletic in tip layout!", cname), call. = FALSE)
  }
  c(xmin = min(xs), xmax = max(xs), n = length(xs))
}

sp_pii2 <- check_clade_span(filter(tip_df, sensor_clade == "Pii-2"), "Pii-2")
sp_nrcs <- check_clade_span(filter(tip_df, sensor_clade == "NRC-S"), "NRC-S")
sp_solo <- check_clade_span(filter(tip_df, is_solo1), "SOLO1")
sp_pia2 <- check_clade_span(filter(tip_df, sensor_clade == "Pia-2/MIC1"), "Pia-2/MIC1")

if (!(sp_pii2["xmin"] < sp_nrcs["xmin"] && sp_nrcs["xmin"] < sp_solo["xmin"] && sp_solo["xmin"] < sp_pia2["xmin"])) {
  stop("Tree ladderization does not match Fig. 2A L->R clade order (Pii-2 -> NRC-S -> SOLO1 -> Pia-2/MIC1)!", call. = FALSE)
}

tested_df <- tip_df %>% filter(is_tested_sr)
n_yes <- sum(tested_df$Cell_death == "Yes")
n_no  <- sum(tested_df$Cell_death == "No")
if (nrow(tested_df) != 34L || n_yes != 11L || n_no != 23L) {
  stop(sprintf("Unexpected structural ranking tested counts: total=%d, Yes=%d, No=%d",
               nrow(tested_df), n_yes, n_no), call. = FALSE)
}

n_iptm_gt <- sum(tip_df$iptm_group == "gt")
n_iptm_le <- sum(tip_df$iptm_group == "le")
cat(sprintf("Verified 296 tips | max depth diff=%.2e | ipTM > %.2f: %d, <= %.2f: %d | Tested (SR): 34 (Yes=%d, No=%d)\n",
            depth_diff, opt$`iptm-threshold`, n_iptm_gt, opt$`iptm-threshold`, n_iptm_le, n_yes, n_no))

# Airy pastel sensor clade washes so colored branches remain crisp & readable
CLADE_WASH <- tibble::tribble(
  ~clade,       ~xmin,                 ~xmax,                 ~fill,      ~face,
  "Pii-2",      sp_pii2["xmin"] - 0.5, sp_pii2["xmax"] + 0.5, "#FAF0D7", "plain",
  "NRC-S",      sp_nrcs["xmin"] - 0.5, sp_nrcs["xmax"] + 0.5, "#ECE7F6", "plain",
  "SOLO1",      sp_solo["xmin"] - 0.5, sp_solo["xmax"] + 0.5, "#E3E6EA", "bold",
  "Pia-2/MIC1", sp_pia2["xmin"] - 0.5, sp_pia2["xmax"] + 0.5, "#E1F2EE", "plain"
)

X_LIMITS  <- c(0.5, 296.5)
MAX_DEPTH <- max(tip_df$tip_depth)

# Uniform left & right margins across all stacked panels so [0.5, 296.5] aligns pixel-for-pixel.
# Left side carries both inner ticks/counts (x = -3.4) and outer track titles (x = -26.5).
LEFT_MARG_PT  <- 90
RIGHT_MARG_PT <- 14
TITLE_X       <- -26.5

# ---------------------------------------------------------------------------
# Bottom Panel: Upward Phylogram (ggtree + coord_flip)
# ---------------------------------------------------------------------------
shade_ymin <- -0.045 * MAX_DEPTH
shade_ymax <- MAX_DEPTH * 1.015

p_tree <- p_tree_base

if (show_clade_shades) {
  shade_layer <- geom_rect(
    data = CLADE_WASH,
    aes(xmin = shade_ymin, xmax = shade_ymax, ymin = xmin, ymax = xmax, fill = fill),
    inherit.aes = FALSE,
    color = NA
  )
  p_tree$layers <- c(list(shade_layer), p_tree$layers)
  p_tree <- p_tree +
    scale_fill_identity() +
    geom_text(
      data = CLADE_WASH,
      aes(x = shade_ymin - 0.045 * MAX_DEPTH, y = (xmin + xmax) / 2, label = clade, fontface = face),
      inherit.aes = FALSE,
      size = 3.00,
      family = sys_font,
      color = "#222222",
      vjust = 1.0,
      hjust = 0.5
    )
}

# Direct tip-to-marker vertical connectors & subtle tip-end dots for the 34 tested branches
# All 11 Cell death Yes share unified crimson #D0242A; Rank 1-2 emphasized by size & black rim.
tested_conn_df <- tested_df %>%
  mutate(
    conn_col  = ifelse(Cell_death == "Yes", CD_YES_COLOR, NOCD_LINE),
    conn_lty  = ifelse(is_rank12, "dashed", "dotted"),
    conn_lw   = case_when(is_rank12 ~ 0.36, Cell_death == "Yes" ~ 0.22, TRUE ~ 0.18),
    conn_alp  = case_when(is_rank12 ~ 0.78, Cell_death == "Yes" ~ 0.45, TRUE ~ 0.45),
    dot_col   = case_when(is_rank12 ~ "#111111", Cell_death == "Yes" ~ CD_YES_COLOR, TRUE ~ "#64748B"),
    dot_fill  = ifelse(Cell_death == "Yes", CD_YES_COLOR, "#64748B"),
    dot_shape = ifelse(is_rank12, 21L, 16L),
    dot_sz    = case_when(is_rank12 ~ 1.25, Cell_death == "Yes" ~ 0.72, TRUE ~ 0.55),
    dot_str   = ifelse(is_rank12, 0.42, 0.0)
  )

p_tree <- p_tree +
  geom_segment(
    data = tested_conn_df,
    aes(x = tip_depth, xend = MAX_DEPTH * 1.015, y = x_pos, yend = x_pos,
        color = conn_col, linetype = conn_lty, linewidth = conn_lw, alpha = conn_alp),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = tested_conn_df,
    aes(x = tip_depth, y = x_pos, color = dot_col, fill = dot_fill,
        shape = dot_shape, size = dot_sz, stroke = dot_str),
    inherit.aes = FALSE
  ) +
  scale_color_identity() +
  scale_shape_identity() +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  scale_size_identity() +
  scale_alpha_identity()

if (!show_clade_shades) {
  p_tree <- p_tree + scale_fill_identity()
}

# Bootstrap support markers (strictly verbatim from Newick node.label)
int_nodes <- tree_data %>%
  filter(!isTip) %>%
  mutate(
    bs = suppressWarnings(as.numeric(label)),
    n_desc = vapply(node, function(nd) length(ape::extract.clade(tr_lad, nd)$tip.label), integer(1))
  )

if (bs_mode == "major") {
  mrca_clades <- c(
    ape::getMRCA(tr_lad, tip_df$label[tip_df$sensor_clade == "Pii-2"]),
    ape::getMRCA(tr_lad, tip_df$label[tip_df$sensor_clade == "NRC-S"]),
    ape::getMRCA(tr_lad, tip_df$label[tip_df$sensor_clade == "Pia-2/MIC1"])
  )
  cand_bs <- int_nodes %>%
    filter(!is.na(bs)) %>%
    filter(node %in% mrca_clades | n_desc >= 44L | (n_desc == 15L & y >= 40 & y <= 55)) %>%
    mutate(is_mrca = node %in% mrca_clades) %>%
    arrange(desc(is_mrca), desc(n_desc))

  kept_idx <- integer(0)
  for (i in seq_len(nrow(cand_bs))) {
    xi <- cand_bs$x[i] / MAX_DEPTH
    yi <- cand_bs$y[i] / 296
    too_close <- FALSE
    for (k in kept_idx) {
      xk <- cand_bs$x[k] / MAX_DEPTH
      yk <- cand_bs$y[k] / 296
      if (abs(xi - xk) < 0.115 && abs(yi - yk) < 0.080) {
        too_close <- TRUE
        break
      }
    }
    if (!too_close) kept_idx <- c(kept_idx, i)
  }
  bs_df <- cand_bs[kept_idx, , drop = FALSE] %>%
    mutate(bs_txt = sprintf("%.0f", bs))

  p_tree <- p_tree +
    geom_point(
      data = bs_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 15,
      size = 1.0,
      color = "#111111"
    ) +
    geom_text(
      data = bs_df,
      aes(x = x, y = y, label = bs_txt),
      inherit.aes = FALSE,
      size = 2.20,
      family = sys_font,
      color = "#111111",
      nudge_x = -0.014 * MAX_DEPTH,
      nudge_y = 2.4,
      hjust = 0,
      vjust = 1
    )
} else if (bs_mode == "all95") {
  bs_df <- int_nodes %>% filter(!is.na(bs) & bs >= 95)
  p_tree <- p_tree +
    geom_point(
      data = bs_df,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 16,
      size = 0.85,
      color = "#111111"
    )
}

# Vertical branch-length scale bar of exact length 0.5 substitutions/site on the left
scale_x0 <- 0.05 * MAX_DEPTH
scale_x1 <- scale_x0 + 0.5
p_tree <- p_tree +
  annotate("segment", x = scale_x0, xend = scale_x1, y = -3.4, yend = -3.4,
           color = "#222222", linewidth = 0.45) +
  annotate("text", x = (scale_x0 + scale_x1) / 2, y = -5.8,
           label = "0.5", size = 2.40, family = sys_font, color = "#222222",
           hjust = 1, vjust = 0.5) +
  coord_flip(
    xlim = c(-0.185 * MAX_DEPTH, MAX_DEPTH * 1.015),
    ylim = X_LIMITS,
    expand = FALSE,
    clip = "off"
  ) +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = 0, r = RIGHT_MARG_PT, b = 2, l = LEFT_MARG_PT, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Optional Plant Order Strip (only when --plant-order-mode strip)
# ---------------------------------------------------------------------------
if (order_mode == "strip") {
  p_order <- ggplot(tip_df, aes(x = x_pos, y = 0.5, fill = order_disp)) +
    geom_tile(width = 1.0, height = 0.88, color = NA) +
    scale_fill_manual(values = STRIP_ORDER_COLORS, guide = "none") +
    annotate("text", x = TITLE_X, y = 0.5, label = "Plant order",
             hjust = 1, vjust = 0.5, size = 3.00, family = sys_font, color = "#111111") +
    coord_cartesian(xlim = X_LIMITS, ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void(base_family = sys_font) +
    theme(
      plot.margin = margin(t = 1, r = RIGHT_MARG_PT, b = 0, l = LEFT_MARG_PT, unit = "pt")
    )
}

# ---------------------------------------------------------------------------
# Cell death marker track (unified crimson #D0242A for Yes; Rank 1-2 via size + black rim)
# ---------------------------------------------------------------------------
cd_plot_df <- tested_df %>%
  mutate(
    y_lane     = ifelse(Cell_death == "Yes", 0.68, 0.25),
    conn_col   = ifelse(Cell_death == "Yes", CD_YES_COLOR, NOCD_LINE),
    conn_lty   = ifelse(is_rank12, "dashed", "dotted"),
    conn_lw    = case_when(is_rank12 ~ 0.36, Cell_death == "Yes" ~ 0.22, TRUE ~ 0.18),
    conn_alp   = case_when(is_rank12 ~ 0.78, Cell_death == "Yes" ~ 0.45, TRUE ~ 0.45),
    y_conn_top = y_lane
  ) %>%
  arrange(x_pos)

p_cd <- ggplot() +
  # Subtle horizontal lane guides
  annotate("segment", x = 0.5, xend = 296.5, y = 0.68, yend = 0.68,
           color = "#ECEEF0", linewidth = 0.25) +
  annotate("segment", x = 0.5, xend = 296.5, y = 0.25, yend = 0.25,
           color = "#ECEEF0", linewidth = 0.25) +
  # Vertical connectors from bottom of track (aligned with tree tips) up to each marker
  geom_segment(
    data = cd_plot_df,
    aes(x = x_pos, xend = x_pos, y = 0.0, yend = y_conn_top,
        color = conn_col, linetype = conn_lty, linewidth = conn_lw, alpha = conn_alp)
  ) +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  # No cell death (23): crisp small open circles on exact y = 0.25
  geom_point(
    data = cd_plot_df %>% filter(Cell_death == "No"),
    aes(x = x_pos, y = y_lane),
    shape = 21, size = 0.82, stroke = 0.44, fill = "white", color = NOCD_COLOR
  ) +
  # Regular Cell death Yes (9): unified crimson downward triangles
  geom_point(
    data = cd_plot_df %>% filter(Cell_death == "Yes" & !is_rank12),
    aes(x = x_pos, y = y_lane),
    shape = 25, size = 1.45, stroke = 0.25, fill = CD_YES_COLOR, color = "white"
  ) +
  # Rank 1 & Rank 2 (2): SAME crimson fill, larger size + thick black outline
  geom_point(
    data = cd_plot_df %>% filter(is_rank12),
    aes(x = x_pos, y = y_lane),
    shape = 25, size = 2.40, stroke = 0.75, fill = CD_YES_COLOR, color = "#111111"
  ) +
  # Left-side inner count labels & outer track title on the SAME side
  annotate("text", x = -3.4, y = c(0.68, 0.25),
           label = c(sprintf("Yes (%d)", n_yes), sprintf("No (%d)", n_no)),
           hjust = 1, vjust = 0.5, size = 2.30, family = sys_font, color = "#444444") +
  annotate("text", x = TITLE_X, y = 0.465, label = "Cell death",
           hjust = 1, vjust = 0.5, size = 3.00, family = sys_font, color = "#111111") +
  coord_cartesian(xlim = X_LIMITS, ylim = c(0.0, 0.94), expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = 1, r = RIGHT_MARG_PT, b = 0, l = LEFT_MARG_PT, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Top Panel: Compact ipTM bars (0, 0.5, 1) + unified monotonic callouts
# ---------------------------------------------------------------------------
solve_1d_spacing <- function(x_raw, min_dx = 16.5, xmin = 6, xmax = 248, n_iter = 600) {
  x <- x_raw
  n <- length(x)
  if (n <= 1L) return(x)
  for (it in seq_len(n_iter)) {
    for (i in seq_len(n - 1L)) {
      gap <- x[i + 1L] - x[i]
      if (gap < min_dx) {
        push <- 0.5 * (min_dx - gap)
        x[i]     <- x[i] - push
        x[i + 1L] <- x[i + 1L] + push
      }
    }
    x <- x + 0.04 * (x_raw - x)
    x[1L] <- max(x[1L], xmin)
    for (i in 2:n) x[i] <- max(x[i], x[i - 1L] + min_dx)
    x[n] <- min(x[n], xmax)
    for (i in (n - 1L):1L) x[i] <- min(x[i], x[i + 1L] - min_dx)
  }
  x
}

callout_df <- tip_df %>%
  filter(is_tested_sr & Cell_death == "Yes") %>%
  arrange(x_pos)

y_bar_max <- 1.0
y_top <- if (show_callouts) 2.95 else 1.08

p_iptm <- ggplot(tip_df, aes(x = x_pos, y = ipTM_LCB, fill = iptm_group)) +
  geom_col(width = 0.88, color = NA) +
  scale_fill_manual(values = IPTM_COLORS, guide = "none") +
  # Dashed threshold line at ipTM = 0.5
  annotate("segment", x = 0.5, xend = 296.5,
           y = opt$`iptm-threshold`, yend = opt$`iptm-threshold`,
           color = "#4B5563", linewidth = 0.28, linetype = "dashed") +
  # Left y-axis spanning [0, 1] with ticks at 0, 0.5, 1 + outer track title on LEFT
  annotate("segment", x = -0.8, xend = -0.8, y = 0, yend = y_bar_max,
           color = "#333333", linewidth = 0.32) +
  annotate("segment", x = -2.2, xend = -0.8, y = c(0, 0.5, 1.0), yend = c(0, 0.5, 1.0),
           color = "#333333", linewidth = 0.30) +
  annotate("text", x = -3.4, y = c(0, 0.5, 1.0), label = c("0", "0.5", "1"),
           hjust = 1, vjust = 0.5, size = 2.30, family = sys_font, color = "#333333") +
  annotate("text", x = TITLE_X, y = 0.50, label = "ipTM",
           hjust = 1, vjust = 0.5, size = 3.00, family = sys_font, color = "#111111")

if (show_callouts && nrow(callout_df) > 0L) {
  callout_df <- callout_df %>%
    mutate(
      x_lab  = solve_1d_spacing(x_pos, min_dx = 19.0, xmin = 4, xmax = 242),
      y0     = ipTM_LCB + 0.015,
      yb1    = 1.02,
      yb2    = 1.24,
      y_txt  = 1.30,
      line_col = CD_YES_COLOR,
      line_w   = ifelse(is_rank12, 0.50, 0.28),
      line_alp = ifelse(is_rank12, 1.0, 0.62),
      txt_col  = ifelse(is_rank12, "#000000", "#2D2D2D"),
      txt_sz   = ifelse(is_rank12, 2.55, 2.30),
      sp_clean = gsub("_", " ", species_id),
      # Compact uniform format: Genus species (#rank) without Symbol space glyphs
      expr_lbl = ifelse(
        is_rank12,
        sprintf('paste(bolditalic("%s"), bold(" (#%d)"))', sp_clean, weighted_rank),
        sprintf('paste(italic("%s"), " (#%d)")', sp_clean, weighted_rank)
      )
    )

  p_iptm <- p_iptm +
    geom_segment(
      data = callout_df %>% filter(!is_rank12),
      aes(x = x_pos, xend = x_pos, y = y0, yend = yb1, color = line_col, linewidth = line_w, alpha = line_alp),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.048, "in"), angle = 22, ends = "first", type = "closed")
    ) +
    geom_segment(
      data = callout_df %>% filter(is_rank12),
      aes(x = x_pos, xend = x_pos, y = y0, yend = yb1, color = line_col, linewidth = line_w, alpha = line_alp),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.048, "in"), angle = 20, ends = "first", type = "closed")
    ) +
    geom_segment(
      data = callout_df,
      aes(x = x_pos, xend = x_lab, y = yb1, yend = yb2, color = line_col, linewidth = line_w, alpha = line_alp),
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = callout_df,
      aes(x = x_lab, xend = x_lab, y = yb2, yend = y_txt - 0.02, color = line_col, linewidth = line_w, alpha = line_alp),
      inherit.aes = FALSE
    ) +
    scale_color_identity() +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    geom_text(
      data = callout_df,
      aes(x = x_lab, y = y_txt, label = expr_lbl, color = txt_col, size = txt_sz),
      inherit.aes = FALSE,
      parse = TRUE,
      angle = 30,
      hjust = 0,
      vjust = 0.22,
      family = sys_font
    ) +
    scale_size_identity()
}

p_iptm <- p_iptm +
  coord_cartesian(xlim = X_LIMITS, ylim = c(0, y_top), expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(
    plot.margin = margin(t = 3, r = RIGHT_MARG_PT, b = 1, l = LEFT_MARG_PT, unit = "pt")
  )

# ---------------------------------------------------------------------------
# Bottom Legend: Single crisp framed rectangular box for Plant order only
# (ipTM and Cell death are self-explanatory directly on their tracks)
# ---------------------------------------------------------------------------
p_leg <- NULL
if (order_mode != "none") {
  ord_counts <- table(tip_df$order_disp)
  if (order_mode == "branch") {
    ord_items <- tibble::tribble(
      ~name,            ~x,
      "Poales",         121,
      "Solanales",      186,
      "Lamiales",       257,
      "Caryophyllales", 325,
      "Other",          414,
      "Mixed",          470
    ) %>%
      mutate(
        col = unname(BRANCH_ORDER_COLORS[name]),
        lbl = ifelse(name == "Mixed", "Mixed clades", sprintf("%s (%d)", name, as.integer(ord_counts[name])))
      )

    p_leg <- ggplot() +
      annotate("rect", xmin = 4, xmax = 530, ymin = 0.08, ymax = 0.92,
               fill = "#FAFBFC", color = "#C4CBD4", linewidth = 0.40) +
      annotate("text", x = 12, y = 0.50, label = "Plant order (branch colour):",
               fontface = "bold", hjust = 0, vjust = 0.5, size = 2.55, family = sys_font, color = "#111111") +
      geom_segment(data = ord_items, aes(x = x, xend = x + 9, y = 0.50, yend = 0.50, color = col),
                   inherit.aes = FALSE, linewidth = 1.25, lineend = "round") +
      scale_color_identity() +
      geom_text(data = ord_items, aes(x = x + 12, y = 0.50, label = lbl),
                inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.48, family = sys_font, color = "#222222") +
      coord_cartesian(xlim = c(0, 534), ylim = c(0, 1), expand = FALSE, clip = "off")
  } else {
    ord_items <- tibble::tribble(
      ~name,            ~x,
      "Poales",         122,
      "Solanales",      192,
      "Lamiales",       268,
      "Caryophyllales", 340,
      "Other",          430
    ) %>%
      mutate(
        fill = unname(STRIP_ORDER_COLORS[name]),
        lbl  = sprintf("%s (%d)", name, as.integer(ord_counts[name]))
      )

    p_leg <- ggplot() +
      annotate("rect", xmin = 48, xmax = 486, ymin = 0.08, ymax = 0.92,
               fill = "#FAFBFC", color = "#C4CBD4", linewidth = 0.40) +
      annotate("text", x = 58, y = 0.50, label = "Plant order:",
               fontface = "bold", hjust = 0, vjust = 0.5, size = 2.55, family = sys_font, color = "#111111") +
      geom_rect(data = ord_items, aes(xmin = x, xmax = x + 9, ymin = 0.30, ymax = 0.70, fill = fill),
                inherit.aes = FALSE, color = NA) +
      scale_fill_identity() +
      geom_text(data = ord_items, aes(x = x + 12, y = 0.50, label = lbl),
                inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 2.48, family = sys_font, color = "#222222") +
      coord_cartesian(xlim = c(0, 534), ylim = c(0, 1), expand = FALSE, clip = "off")
  }

  p_leg <- p_leg +
    theme_void(base_family = sys_font) +
    theme(
      plot.margin = margin(t = 2, r = 14, b = 2, l = 14, unit = "pt")
    )
}

# ---------------------------------------------------------------------------
# Assemble with patchwork
# ---------------------------------------------------------------------------
h_iptm <- if (show_callouts) 1.18 else 0.44
if (order_mode == "strip") {
  combined <- p_iptm / p_cd / p_order / p_tree / patchwork::wrap_elements(full = p_leg) +
    plot_layout(heights = c(h_iptm, 0.34, 0.18, 1.48, 0.28))
} else if (order_mode == "branch") {
  combined <- p_iptm / p_cd / p_tree / patchwork::wrap_elements(full = p_leg) +
    plot_layout(heights = c(h_iptm, 0.34, 1.52, 0.28))
} else {
  combined <- p_iptm / p_cd / p_tree +
    plot_layout(heights = c(h_iptm, 0.34, 1.65))
}

# ---------------------------------------------------------------------------
# Export vector PDF, editor-safe SVG, and high-DPI PNG
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
    paste0(
      substr(line, 1L, m_pos[1L] - 1L),
      new_style,
      substr(line, m_pos[1L] + m_len, nchar(line))
    )
  }

  for (line in lines) {
    if (grepl("<clipPath\\b", line)) in_clip <- TRUE
    if (!in_clip) {
      # Drop transparent patchwork sub-panel background <rect>s without fill
      # so editors that ignore .svglite CSS selectors (e.g. Affinity) do not
      # render them with the SVG default black fill.
      if (grepl("^\\s*<rect\\b", line) &&
          grepl("stroke:\\s*none", line) &&
          !grepl("fill:", line)) {
        next
      }
      if (grepl("^\\s*<(rect|circle|polygon|path|polyline)\\b", line) &&
          !grepl("fill:", line)) {
        line <- inject_style(line, "fill: none;")
      }
      if (grepl("^\\s*<(line|polyline|path)\\b", line) &&
          !grepl("stroke-linecap:", line)) {
        line <- inject_style(line, "stroke-linecap: round;")
      }
      if (grepl("^\\s*<text\\b", line) &&
          !grepl("fill:", line)) {
        line <- inject_style(line, "fill: #000000;")
      }
    }
    if (grepl("</clipPath>", line)) in_clip <- FALSE
    k <- k + 1L
    out[k] <- line
  }
  writeLines(out[seq_len(k)], path, useBytes = FALSE)
}

out_prefix <- opt$`output-prefix`
out_dir <- dirname(out_prefix)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pdf_path <- paste0(out_prefix, ".pdf")
svg_path <- paste0(out_prefix, ".svg")
png_path <- paste0(out_prefix, ".png")

if (isTRUE(capabilities()[["aqua"]])) {
  if (identical(sys_font, "Helvetica Neue")) {
    grDevices::quartzFonts(
      "Helvetica Neue" = grDevices::quartzFont(c(
        "HelveticaNeue",
        "HelveticaNeue-Bold",
        "HelveticaNeue-Italic",
        "HelveticaNeue-BoldItalic"
      ))
    )
  }
  grDevices::quartz(
    type = "pdf",
    file = pdf_path,
    width = opt$width,
    height = opt$height,
    family = sys_font,
    bg = "white"
  )
} else {
  grDevices::cairo_pdf(
    pdf_path,
    width = opt$width,
    height = opt$height,
    family = sys_font,
    bg = "white"
  )
}
print(combined)
grDevices::dev.off()

svglite::svglite(svg_path, width = opt$width, height = opt$height, bg = "white", fix_text_size = FALSE)
print(combined)
grDevices::dev.off()
sanitize_svg_for_editors(svg_path)

ragg::agg_png(png_path, width = opt$width, height = opt$height, units = "in", res = opt$dpi, background = "white")
print(combined)
grDevices::dev.off()

cat(sprintf("Wrote:\n  %s\n  %s\n  %s\n", pdf_path, svg_path, png_path))
