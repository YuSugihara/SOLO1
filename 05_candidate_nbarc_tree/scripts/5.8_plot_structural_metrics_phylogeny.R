#!/usr/bin/env Rscript
#' Step 5.8 script: plot the structural metrics across the candidate phylogeny for Fig. S2.
#'
#' Displays all six AlphaFold 3 structural metrics across the 296 NLR-ID candidates
#' plus the NbNRC2a benchmark alongside the upward-facing midpoint-rooted NB-ARC
#' phylogeny, plant-order branch colouring, clade background shading, major
#' bootstrap support values, and the 3-class functional screening cell-death
#' phenotypes (DV-dependent, Autoactive, No cell death).

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

option_list <- list(
  make_option("--tree", type = "character",
    default = "05_candidate_nbarc_tree/data/results/5.5_candidate_nbarc.midpoint_rooted.treefile",
    help = "Midpoint-rooted Newick tree file [default %default]"),
  make_option("--annotations", type = "character",
    default = "05_candidate_nbarc_tree/data/results/5.6_candidate_tip_annotations.tsv",
    help = "Unified tip annotation table from Step 5.6 [default %default]"),
  make_option("--scores", type = "character",
    default = "04_ranking/data/results/4.1_ranked_scores.tsv",
    help = "Full 297-row structural ranking table (with NbNRC2a) [default %default]"),
  make_option("--output-prefix", type = "character",
    default = "05_candidate_nbarc_tree/data/results/5.8_structural_metrics_phylogeny",
    help = "Output path prefix without extension; writes .pdf, .svg, .png [default %default]"),
  make_option("--cd-mode", type = "character", default = "3class",
    help = "Cell death display mode: '3class' (DV-dependent/Autoactive/No cell death) or '2class' (Yes/No)"),
  make_option("--width", type = "double", default = 7.8,
    help = "Figure width in inches"),
  make_option("--height", type = "double", default = 5.10,
    help = "Figure height in inches"),
  make_option("--dpi", type = "integer", default = 600,
    help = "PNG resolution in DPI")
)

opt <- parse_args(OptionParser(option_list = option_list))
cd_mode <- tolower(trimws(opt$`cd-mode`))
if (!cd_mode %in% c("3class", "2class")) {
  stop("--cd-mode must be '3class' or '2class'", call. = FALSE)
}

pick_sans_font <- function() {
  for (cand in c("Helvetica Neue", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans")) {
    matched <- tryCatch(systemfonts::match_fonts(cand), error = function(e) NULL)
    if (!is.null(matched$path) && nzchar(matched$path)) return(cand)
  }
  "sans"
}
sys_font <- pick_sans_font()

# ---------------------------------------------------------------------------
# Load & verify inputs (read-only)
# ---------------------------------------------------------------------------
for (p in c(opt$tree, opt$annotations, opt$scores)) {
  if (!file.exists(p)) stop(sprintf("Required input not found: %s", p), call. = FALSE)
}

tr <- ape::read.tree(opt$tree)
if (ape::Ntip(tr) != 296L) {
  stop(sprintf("Expected 296 tips, got %d", ape::Ntip(tr)), call. = FALSE)
}
orig_depths <- setNames(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))], tr$tip.label)
orig_edges_sorted <- sort(tr$edge.length)
tr_lad <- ape::ladderize(tr, right = FALSE)
if (max(abs(sort(tr_lad$edge.length) - orig_edges_sorted)) > 1e-12) {
  stop("Ladderization altered branch lengths!", call. = FALSE)
}

ann_df <- readr::read_tsv(
  opt$annotations,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)
scores_df <- readr::read_tsv(
  opt$scores,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  na = character()
)

if (nrow(ann_df) != 296L || !setequal(tr_lad$tip.label, ann_df$tip_label)) {
  stop("Annotations TSV must match the 296 tree tips 1-to-1", call. = FALSE)
}
if (nrow(scores_df) != 297L) {
  stop(sprintf("Expected 297 rows in ranked scores TSV, got %d", nrow(scores_df)), call. = FALSE)
}

nrc2a_row <- scores_df %>% filter(protein_ID == "NbNRC2a")
if (nrow(nrc2a_row) != 1L) {
  stop(sprintf("Expected exactly 1 NbNRC2a benchmark row, got %d", nrow(nrc2a_row)), call. = FALSE)
}

# ---------------------------------------------------------------------------
# Colors & shared geometry (identical branch palette to Fig. 2A)
# ---------------------------------------------------------------------------
BRANCH_ORDER_COLORS <- c(
  "Poales"         = "#0B5394",
  "Solanales"      = "#D87B00",
  "Lamiales"       = "#00947E",
  "Caryophyllales" = "#B8195A",
  "Other"          = "#969EA8",
  "Mixed"          = "#1E1E1E"
)

METRIC_RED   <- "#B8252A"
METRIC_NAVY  <- "#15415C"
CD_DV_COL    <- "#D0242A"
CD_AUTO_COL  <- "#009E73"
NOCD_COL     <- "#525252"
NOCD_LINE    <- "#94A3B8"
SEP_X        <- 300.5
NRC_X        <- 306.0
NRC_W        <- 4.8
X_LIMITS     <- c(0.5, 311.5)

LEFT_MARG_PT  <- 126
RIGHT_MARG_PT <- 10
METRIC_TITLE_X <- -29

# Plant-order consensus coloring
top4_orders <- c("Poales", "Solanales", "Lamiales", "Caryophyllales")
n_tip <- ape::Ntip(tr_lad)
n_tot <- n_tip + ape::Nnode(tr_lad)
ann_idx <- match(tr_lad$tip.label, ann_df$tip_label)
tip_ord_disp <- ann_df$plant_order_group[ann_idx]
names(tip_ord_disp) <- tr_lad$tip.label

node_order <- character(n_tot)
node_order[seq_len(n_tip)] <- tip_ord_disp
for (nd in (n_tip + 1L):n_tot) {
  d_tips <- ape::extract.clade(tr_lad, nd)$tip.label
  u <- unique(tip_ord_disp[d_tips])
  node_order[nd] <- if (length(u) == 1L) u else "Mixed"
}

p_tree_base <- ggtree(tr_lad, ladderize = FALSE, linewidth = 0.34, lineend = "square")
p_tree_base$data$branch_col <- unname(BRANCH_ORDER_COLORS[node_order[p_tree_base$data$node]])
p_tree_base$mapping <- utils::modifyList(p_tree_base$mapping, aes(color = branch_col))

tree_data <- p_tree_base$data
tip_df <- tree_data %>%
  filter(isTip) %>%
  mutate(x_pos = y, tip_depth = x) %>%
  arrange(x_pos)

depth_diff <- max(abs(tip_df$tip_depth - orig_depths[tip_df$label]))
if (depth_diff > 1e-9) {
  stop(sprintf("Tip depth mismatch: %.3e", depth_diff), call. = FALSE)
}

scores_join <- scores_df %>%
  filter(protein_ID != "NbNRC2a") %>%
  select(protein_ID, SD_THETA_ROT, S_PROTO, Sum_CONTACTS, BSA_INT_PROTO)

tip_df <- tip_df %>%
  left_join(ann_df, by = c("label" = "tip_label")) %>%
  left_join(scores_join, by = "protein_ID") %>%
  mutate(
    weighted_rank = as.integer(weighted_rank),
    D_APEX        = as.numeric(D_APEX),
    SD_THETA_ROT  = as.numeric(SD_THETA_ROT),
    S_PROTO       = as.numeric(S_PROTO),
    Sum_CONTACTS  = as.numeric(Sum_CONTACTS),
    BSA_INT_PROTO = as.numeric(BSA_INT_PROTO),
    ipTM_LCB      = as.numeric(ipTM_LCB),
    is_solo1      = is_solo1_clade == "True",
    order_disp    = factor(plant_order_group, levels = c(top4_orders, "Other")),
    is_tested_sr  = screening_stage == "structural ranking",
    is_rank12     = weighted_rank %in% c(1L, 2L)
  )

metric_cols <- c("D_APEX", "SD_THETA_ROT", "S_PROTO", "Sum_CONTACTS", "BSA_INT_PROTO", "ipTM_LCB")
for (mc in metric_cols) {
  if (any(is.na(tip_df[[mc]])) || is.na(as.numeric(nrc2a_row[[mc]]))) {
    stop(sprintf("Missing numeric values in metric '%s'", mc), call. = FALSE)
  }
}

# Verify clade spans (Pii-2 -> NRC-S -> SOLO1 -> Pia-2/MIC1)
check_clade_span <- function(df_sub, cname) {
  xs <- sort(df_sub$x_pos)
  if (length(xs) == 0 || length(xs) != (max(xs) - min(xs) + 1L)) {
    stop(sprintf("Clade '%s' is not monophyletic/contiguous", cname), call. = FALSE)
  }
  c(xmin = min(xs), xmax = max(xs), n = length(xs))
}
sp_pii2 <- check_clade_span(filter(tip_df, sensor_clade == "Pii-2"), "Pii-2")
sp_nrcs <- check_clade_span(filter(tip_df, sensor_clade == "NRC-S"), "NRC-S")
sp_solo <- check_clade_span(filter(tip_df, is_solo1), "SOLO1")
sp_pia2 <- check_clade_span(filter(tip_df, sensor_clade == "Pia-2/MIC1"), "Pia-2/MIC1")

tested_df <- tip_df %>% filter(is_tested_sr)
n_dv   <- sum(tested_df$screening_class == "DV-dependent cell death")
n_auto <- sum(tested_df$screening_class == "Autoactive")
n_nocd <- sum(tested_df$screening_class == "No cell death")
n_yes  <- sum(tested_df$Cell_death == "Yes")
if (nrow(tested_df) != 34L || n_dv != 8L || n_auto != 3L || n_nocd != 23L || n_yes != 11L) {
  stop(sprintf("Unexpected screened counts: total=%d DV=%d Auto=%d No=%d Yes=%d",
               nrow(tested_df), n_dv, n_auto, n_nocd, n_yes), call. = FALSE)
}

cat(sprintf("Verified 296 tips + NbNRC2a | depth diff=%.2e | Screened 34: DV=%d, Auto=%d, No=%d\n",
            depth_diff, n_dv, n_auto, n_nocd))

# ---------------------------------------------------------------------------
# Metric Track Builder
# ---------------------------------------------------------------------------
build_metric_panel <- function(col_name, label_expr, fill_col, y_max, ticks, ref_line = NULL) {
  df_m <- tibble::tibble(
    x_pos = tip_df$x_pos,
    val   = tip_df[[col_name]]
  )
  nrc_val <- as.numeric(nrc2a_row[[col_name]])
  if (max(df_m$val, nrc_val) > y_max) {
    stop(sprintf("Data exceeds y_max=%.2f for %s (max=%.2f)", y_max, col_name, max(df_m$val, nrc_val)), call. = FALSE)
  }

  p <- ggplot(df_m, aes(x = x_pos, y = val)) +
    # Baseline rules
    annotate("segment", x = 0.5, xend = 296.5, y = 0, yend = 0,
             color = "#CBD5E1", linewidth = 0.22) +
    annotate("segment", x = NRC_X - NRC_W / 2, xend = NRC_X + NRC_W / 2, y = 0, yend = 0,
             color = "#CBD5E1", linewidth = 0.22) +
    geom_col(width = 0.88, fill = fill_col, color = NA) +
    annotate("rect", xmin = NRC_X - NRC_W / 2, xmax = NRC_X + NRC_W / 2,
             ymin = 0, ymax = nrc_val, fill = fill_col, color = NA) +
    annotate("segment", x = SEP_X, xend = SEP_X, y = 0, yend = y_max,
             color = "#94A3B8", linewidth = 0.32)

  if (!is.null(ref_line)) {
    p <- p + annotate("segment", x = 0.5, xend = 296.5, y = ref_line, yend = ref_line,
                      color = "#4B5563", linewidth = 0.25, linetype = "dashed")
  }

  p +
    # Crisp left y-axis spine & ticks
    annotate("segment", x = -1.0, xend = -1.0, y = 0, yend = y_max,
             color = "#333333", linewidth = 0.30) +
    annotate("segment", x = -2.6, xend = -1.0, y = ticks, yend = ticks,
             color = "#333333", linewidth = 0.28) +
    annotate("text", x = -4.2, y = ticks, label = as.character(ticks),
             hjust = 1, vjust = c(0.15, 0.50, 0.85), size = 2.05, family = sys_font, color = "#333333") +
    annotate("text", x = METRIC_TITLE_X, y = y_max * 0.50, label = label_expr,
             parse = is.language(label_expr),
             hjust = 1, vjust = 0.5, size = 2.65, family = sys_font, color = "#111111") +
    coord_cartesian(xlim = X_LIMITS, ylim = c(0, y_max), expand = FALSE, clip = "off") +
    theme_void(base_family = sys_font) +
    theme(plot.margin = margin(t = 3.2, r = RIGHT_MARG_PT, b = 2.8, l = LEFT_MARG_PT, unit = "pt"))
}

p_d_apex <- build_metric_panel(
  "D_APEX", expression(italic(D)[APEX]~"(" * ring(A) * ")"), METRIC_RED, y_max = 138, ticks = c(0, 50, 100)
)
p_sd_rot <- build_metric_panel(
  "SD_THETA_ROT", expression(sigma[theta~ROT]~"(" * degree * ")"), METRIC_RED, y_max = 63, ticks = c(0, 30, 60)
)
p_s_proto <- build_metric_panel(
  "S_PROTO", expression(italic(S)[PROTO]~"(" * ring(A) * ")"), METRIC_RED, y_max = 6.2, ticks = c(0, 3, 6)
)
p_contacts <- build_metric_panel(
  "Sum_CONTACTS", expression(Sigma[CONTACTS]), METRIC_NAVY, y_max = 400, ticks = c(0, 200, 400)
)
p_bsa <- build_metric_panel(
  "BSA_INT_PROTO", expression(italic(BSA)["INTER-PROTO"]~"(" * ring(A)^2 * ")"), METRIC_NAVY, y_max = 3500, ticks = c(0, 1500, 3000)
)
p_iptm <- build_metric_panel(
  "ipTM_LCB", "ipTM", METRIC_NAVY, y_max = 1.0, ticks = c(0, 0.5, 1), ref_line = 0.5
)

# ---------------------------------------------------------------------------
# Cell Death Phenotype Track (3-class default matching Fig. S2 legend & Fig. 2A grammar)
# ---------------------------------------------------------------------------
if (cd_mode == "3class") {
  cd_plot_df <- tested_df %>%
    mutate(
      y_lane = case_when(
        screening_class == "DV-dependent cell death" ~ 0.82,
        screening_class == "Autoactive"              ~ 0.50,
        TRUE                                         ~ 0.18
      ),
      conn_col = case_when(
        screening_class == "DV-dependent cell death" ~ CD_DV_COL,
        screening_class == "Autoactive"              ~ CD_AUTO_COL,
        TRUE                                         ~ NOCD_LINE
      ),
      conn_lty = ifelse(is_rank12, "dashed", "dotted"),
      conn_lw  = case_when(
        is_rank12                                    ~ 0.36,
        screening_class == "DV-dependent cell death" ~ 0.22,
        screening_class == "Autoactive"              ~ 0.22,
        TRUE                                         ~ 0.18
      ),
      conn_alp = ifelse(is_rank12, 0.78, 0.48),
      dot_col  = case_when(
        is_rank12                                    ~ "#111111",
        screening_class == "DV-dependent cell death" ~ CD_DV_COL,
        screening_class == "Autoactive"              ~ CD_AUTO_COL,
        TRUE                                         ~ "#64748B"
      ),
      dot_fill = case_when(
        screening_class == "DV-dependent cell death" ~ CD_DV_COL,
        screening_class == "Autoactive"              ~ CD_AUTO_COL,
        TRUE                                         ~ "#64748B"
      ),
      dot_shape = ifelse(is_rank12, 21L, 16L),
      dot_sz = case_when(
        is_rank12                                    ~ 1.25,
        screening_class == "DV-dependent cell death" ~ 0.72,
        screening_class == "Autoactive"              ~ 0.72,
        TRUE                                         ~ 0.55
      ),
      dot_str = ifelse(is_rank12, 0.42, 0.0)
    ) %>%
    arrange(x_pos)

  p_cd <- ggplot() +
    annotate("segment", x = 0.5, xend = 296.5, y = c(0.82, 0.50, 0.18), yend = c(0.82, 0.50, 0.18),
             color = "#ECEEF0", linewidth = 0.25) +
    annotate("segment", x = SEP_X, xend = SEP_X, y = 0, yend = 1.0,
             color = "#94A3B8", linewidth = 0.32) +
    geom_segment(
      data = cd_plot_df,
      aes(x = x_pos, xend = x_pos, y = 0.0, yend = y_lane,
          color = conn_col, linetype = conn_lty, linewidth = conn_lw, alpha = conn_alp)
    ) +
    scale_color_identity() +
    scale_linetype_identity() +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    # No cell death (23)
    geom_point(
      data = filter(cd_plot_df, screening_class == "No cell death"),
      aes(x = x_pos, y = y_lane),
      shape = 21, size = 0.82, stroke = 0.44, fill = "white", color = NOCD_COL
    ) +
    # Autoactive (3)
    geom_point(
      data = filter(cd_plot_df, screening_class == "Autoactive"),
      aes(x = x_pos, y = y_lane),
      shape = 21, size = 1.65, stroke = 0.28, fill = CD_AUTO_COL, color = "white"
    ) +
    # Regular DV-dependent (6)
    geom_point(
      data = filter(cd_plot_df, screening_class == "DV-dependent cell death" & !is_rank12),
      aes(x = x_pos, y = y_lane),
      shape = 25, size = 1.45, stroke = 0.25, fill = CD_DV_COL, color = "white"
    ) +
    # Rank #1 & #2 SOLO1 DV-dependent hits (2)
    geom_point(
      data = filter(cd_plot_df, is_rank12),
      aes(x = x_pos, y = y_lane),
      shape = 25, size = 2.25, stroke = 0.70, fill = CD_DV_COL, color = "#111111"
    ) +
    # Clean left-margin lane labels aligned with metric y-axis ticks (no extra marker glyphs)
    annotate("text", x = -4.2, y = c(0.82, 0.50, 0.18),
             label = c(sprintf("DV-dependent (%d)", n_dv),
                       sprintf("Autoactive (%d)", n_auto),
                       sprintf("No cell death (%d)", n_nocd)),
             hjust = 1, vjust = 0.5, size = 2.05, family = sys_font, color = "#333333") +
    annotate("text", x = -57, y = 0.50, label = "Cell death",
             hjust = 1, vjust = 0.5, size = 2.65, family = sys_font, color = "#111111") +
    coord_cartesian(xlim = X_LIMITS, ylim = c(0.0, 1.02), expand = FALSE, clip = "off") +
    theme_void(base_family = sys_font) +
    theme(plot.margin = margin(t = 3.0, r = RIGHT_MARG_PT, b = 0.5, l = LEFT_MARG_PT, unit = "pt"))
} else {
  cd_plot_df <- tested_df %>%
    mutate(
      y_lane   = ifelse(Cell_death == "Yes", 0.68, 0.25),
      conn_col = ifelse(Cell_death == "Yes", CD_DV_COL, NOCD_LINE),
      conn_lty = ifelse(is_rank12, "dashed", "dotted"),
      conn_lw  = case_when(is_rank12 ~ 0.36, Cell_death == "Yes" ~ 0.22, TRUE ~ 0.18),
      conn_alp = ifelse(is_rank12, 0.78, 0.45),
      dot_col  = case_when(is_rank12 ~ "#111111", Cell_death == "Yes" ~ CD_DV_COL, TRUE ~ "#64748B"),
      dot_fill = ifelse(Cell_death == "Yes", CD_DV_COL, "#64748B"),
      dot_shape = ifelse(is_rank12, 21L, 16L),
      dot_sz   = case_when(is_rank12 ~ 1.25, Cell_death == "Yes" ~ 0.72, TRUE ~ 0.55),
      dot_str  = ifelse(is_rank12, 0.42, 0.0)
    ) %>%
    arrange(x_pos)

  p_cd <- ggplot() +
    annotate("segment", x = 0.5, xend = 296.5, y = c(0.68, 0.25), yend = c(0.68, 0.25),
             color = "#ECEEF0", linewidth = 0.25) +
    annotate("segment", x = SEP_X, xend = SEP_X, y = 0, yend = 0.94,
             color = "#94A3B8", linewidth = 0.32) +
    geom_segment(
      data = cd_plot_df,
      aes(x = x_pos, xend = x_pos, y = 0.0, yend = y_lane,
          color = conn_col, linetype = conn_lty, linewidth = conn_lw, alpha = conn_alp)
    ) +
    scale_color_identity() + scale_linetype_identity() + scale_linewidth_identity() + scale_alpha_identity() +
    geom_point(data = filter(cd_plot_df, Cell_death == "No"), aes(x = x_pos, y = y_lane),
               shape = 21, size = 0.82, stroke = 0.44, fill = "white", color = NOCD_COL) +
    geom_point(data = filter(cd_plot_df, Cell_death == "Yes" & !is_rank12), aes(x = x_pos, y = y_lane),
               shape = 25, size = 1.45, stroke = 0.25, fill = CD_DV_COL, color = "white") +
    geom_point(data = filter(cd_plot_df, is_rank12), aes(x = x_pos, y = y_lane),
               shape = 25, size = 2.25, stroke = 0.70, fill = CD_DV_COL, color = "#111111") +
    annotate("text", x = -4.2, y = c(0.68, 0.25),
             label = c(sprintf("Yes (%d)", n_yes), sprintf("No (%d)", n_nocd)),
             hjust = 1, vjust = 0.5, size = 2.15, family = sys_font, color = "#333333") +
    annotate("text", x = METRIC_TITLE_X, y = 0.465, label = "Cell death",
             hjust = 1, vjust = 0.5, size = 2.65, family = sys_font, color = "#111111") +
    coord_cartesian(xlim = X_LIMITS, ylim = c(0.0, 0.94), expand = FALSE, clip = "off") +
    theme_void(base_family = sys_font) +
    theme(plot.margin = margin(t = 2.0, r = RIGHT_MARG_PT, b = 0.5, l = LEFT_MARG_PT, unit = "pt"))
}

# ---------------------------------------------------------------------------
# Bottom Panel: Upward Phylogeny (exact Fig. 2A layout, colors, clades & BS)
# ---------------------------------------------------------------------------
MAX_DEPTH <- max(tip_df$tip_depth)
shade_ymin <- -0.045 * MAX_DEPTH
shade_ymax <- MAX_DEPTH * 1.015

CLADE_WASH <- tibble::tribble(
  ~clade,       ~xmin,                 ~xmax,                 ~fill,      ~face,
  "Pii-2",      sp_pii2["xmin"] - 0.5, sp_pii2["xmax"] + 0.5, "#FAF0D7", "plain",
  "NRC-S",      sp_nrcs["xmin"] - 0.5, sp_nrcs["xmax"] + 0.5, "#ECE7F6", "plain",
  "SOLO1",      sp_solo["xmin"] - 0.5, sp_solo["xmax"] + 0.5, "#E3E6EA", "bold",
  "Pia-2/MIC1", sp_pia2["xmin"] - 0.5, sp_pia2["xmax"] + 0.5, "#E1F2EE", "plain"
)

p_tree <- p_tree_base
shade_layer <- geom_rect(
  data = CLADE_WASH,
  aes(xmin = shade_ymin, xmax = shade_ymax, ymin = xmin, ymax = xmax, fill = fill),
  inherit.aes = FALSE, color = NA
)
p_tree$layers <- c(list(shade_layer), p_tree$layers)

p_tree <- p_tree +
  scale_fill_identity() +
  geom_text(
    data = CLADE_WASH,
    aes(x = shade_ymin - 0.045 * MAX_DEPTH, y = (xmin + xmax) / 2, label = clade, fontface = face),
    inherit.aes = FALSE, size = 2.95, family = sys_font, color = "#222222", vjust = 1.0, hjust = 0.5
  ) +
  geom_segment(
    data = cd_plot_df,
    aes(x = tip_depth, xend = MAX_DEPTH * 1.015, y = x_pos, yend = x_pos,
        color = conn_col, linetype = conn_lty, linewidth = conn_lw, alpha = conn_alp),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = cd_plot_df,
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

# Major bootstrap support squares & labels
int_nodes <- tree_data %>%
  filter(!isTip) %>%
  mutate(
    bs = suppressWarnings(as.numeric(label)),
    n_desc = vapply(node, function(nd) length(ape::extract.clade(tr_lad, nd)$tip.label), integer(1))
  )

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
bs_df <- cand_bs[kept_idx, , drop = FALSE] %>% mutate(bs_txt = sprintf("%.0f", bs))

p_tree <- p_tree +
  geom_point(data = bs_df, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 15, size = 0.95, color = "#111111") +
  geom_text(data = bs_df, aes(x = x, y = y, label = bs_txt), inherit.aes = FALSE,
            size = 2.15, family = sys_font, color = "#111111",
            nudge_x = -0.014 * MAX_DEPTH, nudge_y = 2.4, hjust = 0, vjust = 1) +
  # Vertical divider & NbNRC2a benchmark label on the right
  annotate("segment", x = shade_ymin, xend = shade_ymax, y = SEP_X, yend = SEP_X,
           color = "#94A3B8", linewidth = 0.32) +
  annotate("text", x = MAX_DEPTH * 0.48, y = NRC_X, label = "NbNRC2a",
           angle = 90, size = 2.65, family = sys_font, color = "#222222", hjust = 0.5, vjust = 0.5)

scale_x0 <- 0.05 * MAX_DEPTH
scale_x1 <- scale_x0 + 0.5
p_tree <- p_tree +
  annotate("segment", x = scale_x0, xend = scale_x1, y = -3.4, yend = -3.4,
           color = "#222222", linewidth = 0.45) +
  annotate("text", x = (scale_x0 + scale_x1) / 2, y = -6.0, label = "0.5",
           size = 2.30, family = sys_font, color = "#222222", hjust = 1, vjust = 0.5) +
  coord_flip(xlim = c(-0.20 * MAX_DEPTH, MAX_DEPTH * 1.015), ylim = X_LIMITS, expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(plot.margin = margin(t = 0, r = RIGHT_MARG_PT, b = 2, l = LEFT_MARG_PT, unit = "pt"))

# ---------------------------------------------------------------------------
# Bottom Legend: Single framed Plant order bar matching Fig. 2A
# ---------------------------------------------------------------------------
ord_counts <- table(tip_df$order_disp)
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
  coord_cartesian(xlim = c(0, 534), ylim = c(0, 1), expand = FALSE, clip = "off") +
  theme_void(base_family = sys_font) +
  theme(plot.margin = margin(t = 2, r = 14, b = 2, l = 14, unit = "pt"))

# ---------------------------------------------------------------------------
# Assemble with patchwork
# ---------------------------------------------------------------------------
combined <- p_d_apex / p_sd_rot / p_s_proto / p_contacts / p_bsa / p_iptm / p_cd / p_tree /
  patchwork::wrap_elements(full = p_leg) +
  plot_layout(heights = c(0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.44, 1.38, 0.25))

# ---------------------------------------------------------------------------
# Export PDF, Affinity-safe SVG, and 600 DPI PNG
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

svglite::svglite(svg_path, width = opt$width, height = opt$height, bg = "white", fix_text_size = FALSE)
print(combined)
grDevices::dev.off()
sanitize_svg_for_editors(svg_path)

ragg::agg_png(png_path, width = opt$width, height = opt$height, units = "in", res = opt$dpi, background = "white")
print(combined)
grDevices::dev.off()

cat(sprintf("Wrote:\n  %s\n  %s\n  %s\n", pdf_path, svg_path, png_path))
