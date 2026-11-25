#!/usr/bin/env Rscript
#' Step 6.9 script: plot the 14,364-tip plant NLR NB-ARC phylogeny for Fig. S4.
#'
#' Draws the open-fan circular phylogeny using native R / ggtree with:
#'   1. Inner N-terminal domain ring (CC with ID, CC without ID, TIR, RPW8, Unclassified).
#'   2. Outer Clade / ID highlight track (monophyletic SOLO1 clade n = 136, DDE 176, BED 160).
#'   3. Domain- and support-backed highlights for SOLO1, MIC1, MIC2, and MIC3.
#'   4. Tip-to-ring dotted radial guides and non-overlapping callouts for 21 reference NLRs.
#'   5. Two-section vertical legend, 1.0 scale bar, and Type 1 Helvetica vector PDF output.

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggnewscale)
  library(dplyr)
  library(readr)
  library(grid)
  library(ragg)
})

# ---------------------------------------------------------------------------
# CLI options
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--tree", type = "character",
    default = "06_plant_nlr_nbarc_tree/data/results/6.6_merged_nbarc.vft.tree",
    help = "NB-ARC Newick tree file [default %default]"),
  make_option("--classification", type = "character",
    default = "06_plant_nlr_nbarc_tree/data/results/6.7_n_terminal_domain_classification.tsv",
    help = "N-terminal domain classification table from Step 6.7 [default %default]"),
  make_option("--solo1-ids", type = "character",
    default = "07_solo1_homologues/data/inputs/7.0_solo1_clade_ids.txt",
    help = "SOLO1 clade tip identifier list [default %default]"),
  make_option("--id-gff", type = "character",
    default = "01_nlr_id_dataset/data/results/1.1_nlr_id_candidates.gff",
    help = "Integrated-domain candidate GFF annotation file [default %default]"),
  make_option("--nlrtracker-table", type = "character",
    default = "01_nlr_id_dataset/data/inputs/refseq_nlrtracker_tables.tsv.gz",
    help = "RefSeq NLRtracker classification table [default %default]"),
  make_option("--output-prefix", type = "character",
    default = "06_plant_nlr_nbarc_tree/data/results/6.9_plant_nlr_phylogeny",
    help = "Output path prefix without extension; writes .pdf and .png [default %default]"),
  make_option("--width", type = "double", default = 8.65,
    help = "Figure width in inches [default %default]"),
  make_option("--height", type = "double", default = 6.55,
    help = "Figure height in inches [default %default]"),
  make_option("--dpi", type = "integer", default = 300,
    help = "PNG resolution in DPI [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

tree_path  <- opt$tree
tsv_path   <- opt$classification
solo1_path <- opt$`solo1-ids`
gff_path   <- opt$`id-gff`
nlrt_path  <- opt$`nlrtracker-table`
out_prefix <- opt$`output-prefix`

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
out_png <- paste0(out_prefix, ".png")
out_pdf <- paste0(out_prefix, ".pdf")

# ---- 1. Load tree & metadata ----
tr <- read.tree(tree_path)
df <- read_tsv(tsv_path, na = character(), show_col_types = FALSE)
nlrt <- read_tsv(nlrt_path, show_col_types = FALSE)
solo1_ids <- readLines(solo1_path)
solo1_ids <- solo1_ids[nzchar(solo1_ids)]

cat_vec  <- setNames(df$color_category, df$tip_id)
arch_vec <- setNames(df$domain_architecture, df$tip_id)

# ---- 2. Programmatically find TNL + RNL outgroup clade and root ----
# Largest edge-bipartition containing all RPW8 tips and zero CC/CC_WITH_O tips
rpw8_tips_all <- tr$tip.label[cat_vec[tr$tip.label] == "RPW8"]
cc_tips_all   <- tr$tip.label[cat_vec[tr$tip.label] %in% c("CC", "CC_WITH_O")]

n_tip <- Ntip(tr)
n_tot <- n_tip + Nnode(tr)
is_rpw8 <- integer(n_tot)
is_cc   <- integer(n_tot)
sub_sz  <- integer(n_tot)
is_rpw8[seq_len(n_tip)] <- as.integer(tr$tip.label %in% rpw8_tips_all)
is_cc[seq_len(n_tip)]   <- as.integer(tr$tip.label %in% cc_tips_all)
sub_sz[seq_len(n_tip)]  <- 1L

tr_po <- reorder(tr, order = "postorder")
for (e_i in seq_len(nrow(tr_po$edge))) {
  p <- tr_po$edge[e_i, 1]
  c <- tr_po$edge[e_i, 2]
  is_rpw8[p] <- is_rpw8[p] + is_rpw8[c]
  is_cc[p]   <- is_cc[p]   + is_cc[c]
  sub_sz[p]  <- sub_sz[p]  + sub_sz[c]
}

tot_rpw8 <- length(rpw8_tips_all)
tot_cc   <- length(cc_tips_all)
best_node <- NA_integer_
best_side <- NA_character_
best_sz   <- -1L

for (nd in seq_len(n_tot)) {
  if (nd == n_tip + 1L) next
  if (is_rpw8[nd] == tot_rpw8 && is_cc[nd] == 0L && sub_sz[nd] > best_sz) {
    best_sz <- sub_sz[nd]; best_node <- nd; best_side <- "sub"
  }
  c_rpw8 <- tot_rpw8 - is_rpw8[nd]
  c_cc   <- tot_cc   - is_cc[nd]
  c_sz   <- n_tip    - sub_sz[nd]
  if (c_rpw8 == tot_rpw8 && c_cc == 0L && c_sz > best_sz) {
    best_sz <- c_sz; best_node <- nd; best_side <- "comp"
  }
}
stopifnot(!is.na(best_node))

sub_tips <- if (best_node <= n_tip) tr$tip.label[best_node] else extract.clade(tr, best_node)$tip.label
outgroup_tips <- if (best_side == "sub") sub_tips else setdiff(tr$tip.label, sub_tips)
cat(sprintf("Outgroup TNL+RNL clade: %d tips\n", length(outgroup_tips)))

tr_r <- root(tr, outgroup = outgroup_tips, resolve.root = TRUE)
tr_r <- ladderize(tr_r, right = FALSE)

p_check <- ggtree(tr_r, layout = "fan", open.angle = 20)
d_chk <- p_check$data
rpw8_y_med <- median(d_chk$y[d_chk$isTip & d_chk$label %in% rpw8_tips_all])
if (rpw8_y_med > n_tip / 2) {
  tr_r <- ladderize(tr_r, right = TRUE)
  p_check <- ggtree(tr_r, layout = "fan", open.angle = 20)
  d_chk <- p_check$data
}
cat(sprintf("Outgroup RPW8 median y = %.1f (of %d)\n",
            median(d_chk$y[d_chk$isTip & d_chk$label %in% rpw8_tips_all]), n_tip))

# ---- 3. Compute exact MRCAs via ape::getMRCA ----
tip_df <- d_chk %>% filter(isTip)
tip_x  <- setNames(tip_df$x, tip_df$label)
tip_y  <- setNames(tip_df$y, tip_df$label)

m_solo1 <- getMRCA(tr_r, solo1_ids)
m_rnl   <- getMRCA(tr_r, rpw8_tips_all)

out_tir <- intersect(outgroup_tips, tr_r$tip.label[cat_vec[tr_r$tip.label] == "TIR"])
m_tnl   <- getMRCA(tr_r, out_tir)

tn_tips <- setdiff(tr_r$tip.label[cat_vec[tr_r$tip.label] == "TIR" & arch_vec[tr_r$tip.label] == "TN"], outgroup_tips)
m_tn    <- getMRCA(tr_r, tn_tips)

node_sup <- function(nd) {
  idx <- nd - n_tip
  if (!is.null(tr_r$node.label) && idx >= 1L && idx <= length(tr_r$node.label)) {
    suppressWarnings(as.numeric(tr_r$node.label[idx]))
  } else {
    NA_real_
  }
}
parent_of <- setNames(tr_r$edge[, 1], tr_r$edge[, 2])

tip_acc <- sub("^.*?((?:XP|NP|YP|WP)_[0-9]+\\.[0-9]+).*$", "\\1", tip_df$label)

# MIC2 (Bailey et al. 2018 Inner MIC2 / C23, DDE clade): core MRCA of DDE hits near ZAR1..Pit (169/176 = 96.0%)
gff_lines <- readLines(gff_path)
dde_lines <- gff_lines[!startsWith(gff_lines, "#") & grepl("DDE", gff_lines, fixed = TRUE)]
dde_acc <- unique(sub(".*\\|", "", vapply(strsplit(dde_lines, "\t"), `[`, "", 1)))
all_dde_tips <- tip_df$label[tip_acc %in% dde_acc]
cand_dde <- tip_df$label[tip_df$label %in% all_dde_tips & tip_df$y >= 700 & tip_df$y <= 2000]
m_mic2 <- getMRCA(tr_r, cand_dde)

# MIC3 (Bailey et al. 2018 C24, BED clade): 100% monophyletic zf-BED NLRtracker hits + Xa1/Xo1/Yr5/Yr7/Rph15
nlrt_bed <- nlrt %>%
  filter(grepl("BED", Domain, ignore.case = TRUE)) %>%
  mutate(acc = sub(".*\\|", "", seqname)) %>%
  distinct(acc)
all_bed_tips <- tip_df$label[tip_acc %in% nlrt_bed$acc]
m_mic3 <- getMRCA(tr_r, union(all_bed_tips, c("Xa1_CONL", "Xo1_CONL", "Yr5_CONL", "Yr7_CONL", "Rph15_CONL")))

# MIC1 (Bailey et al. 2018 Inner MIC1 / C16): walk up from core MRCA while CC_WITH_O >= 64%
m_mic1 <- getMRCA(tr_r, c("RPG5_CNLO", "Pia-2_CNLO", "Pi-ta_CNLO", "YrU1_OCNLO"))
repeat {
  p_nd <- unname(parent_of[as.character(m_mic1)])
  if (is.na(p_nd)) break
  desc_p <- extract.clade(tr_r, p_nd)$tip.label
  frac_id <- mean(cat_vec[desc_p] == "CC_WITH_O", na.rm = TRUE)
  if (frac_id >= 0.64) {
    m_mic1 <- p_nd
  } else {
    break
  }
}

clade_summary <- function(nm, nd) {
  desc <- extract.clade(tr_r, nd)$tip.label
  ys <- tip_y[desc]
  sup <- node_sup(nd)
  cat(sprintf("  %-6s node=%5d  sup=%.3f  n=%4d  CC_ID=%3d (%.1f%%)  DDE=%3d/%d  BED=%3d/%d  y=(%.0f..%.0f)\n",
              nm, nd, ifelse(is.na(sup), -1, sup), length(desc), sum(cat_vec[desc] == "CC_WITH_O"),
              100 * mean(cat_vec[desc] == "CC_WITH_O"),
              sum(desc %in% all_dde_tips), length(all_dde_tips),
              sum(desc %in% all_bed_tips), length(all_bed_tips),
              min(ys), max(ys)))
  list(node = nd, desc = desc, y_min = min(ys), y_max = max(ys), y_mid = 0.5 * (min(ys) + max(ys)))
}
cat("Monophyletic MRCA summary:\n")
s_solo1 <- clade_summary("SOLO1", m_solo1)
s_mic1  <- clade_summary("MIC1",  m_mic1)
s_mic2  <- clade_summary("MIC2",  m_mic2)
s_mic3  <- clade_summary("MIC3",  m_mic3)
s_rnl   <- clade_summary("RNL",   m_rnl)
s_tnl   <- clade_summary("TNL",   m_tnl)

non_tn_max_x <- max(tip_df$x[!(tip_df$label %in% tn_tips) & !(tip_df$y >= 2080 & tip_df$y <= 2150)])
tn_max_x     <- max(tip_df$x[tip_df$label %in% tn_tips])
cat(sprintf("Max non-TN radial x = %.3f, Max TN radial x = %.3f\n", non_tn_max_x, tn_max_x))

# ---- 4. Group branches & assign editorial palette ----
tr_grp <- groupClade(tr_r, c(
  SOLO1 = m_solo1,
  MIC1  = m_mic1,
  MIC2  = m_mic2,
  MIC3  = m_mic3,
  RNL   = m_rnl,
  TNL   = m_tnl,
  TN    = m_tn
))

COL_BG    <- "#525E6B"
COL_SOLO1 <- "#C93B2B"
COL_DDE   <- "#4E9A47"
COL_BED   <- "#9E5DB5"
COL_ID    <- "#155774"
COL_CC    <- "#A9C3D1"
COL_TIR   <- "#C88A2E"
COL_RPW8  <- "#7A3E65"
COL_NA    <- "#EAEFF4"
COL_TXT   <- "#22272E"

branch_cols <- c(
  "0"     = COL_BG,
  "SOLO1" = COL_SOLO1,
  "MIC1"  = "#2E6A84",
  "MIC2"  = "#2E6A84",
  "MIC3"  = "#2E6A84",
  "RNL"   = COL_RPW8,
  "TNL"   = COL_TIR,
  "TN"    = COL_TIR
)

# Radial layout geometry
R_RING_IN  <- ceiling(non_tn_max_x * 10) / 10 + 0.65
R_RING_OUT <- R_RING_IN + 0.38
R_OUT_IN   <- R_RING_OUT + 0.06
R_OUT_OUT  <- R_OUT_IN + 0.19
R_ELBOW    <- R_OUT_OUT + 0.36
R_LABEL    <- R_OUT_OUT + 0.62

START_DEG <- -12.0
START_RAD <- START_DEG * pi / 180.0
Y_LIM_HI  <- n_tip * 360.0 / (360.0 - 20.0)

y_to_cw_deg <- function(y) (START_DEG + (y / Y_LIM_HI) * 360.0) %% 360.0

# ---- 5. Build outer ring data + shared outer Clade / ID highlight track ----
dom_levels <- c(
  "CC with ID (1,071)",
  "CC without ID (10,265)",
  "TIR (146)",
  "RPW8 (40)",
  "Unclassified (2,842)"
)
dom_map <- c(
  "CC_WITH_O" = "CC with ID (1,071)",
  "CC"        = "CC without ID (10,265)",
  "TIR"       = "TIR (146)",
  "RPW8"      = "RPW8 (40)",
  "NA"        = "Unclassified (2,842)"
)

ring_df <- tip_df %>%
  mutate(
    dom_cat = factor(dom_map[cat_vec[label]], levels = dom_levels),
    r_in  = R_RING_IN,
    r_out = R_RING_OUT
  )

ring_cols <- c(
  "CC with ID (1,071)"     = COL_ID,
  "CC without ID (10,265)" = COL_CC,
  "TIR (146)"              = COL_TIR,
  "RPW8 (40)"              = COL_RPW8,
  "Unclassified (2,842)"   = COL_NA
)

outer_levels <- c("SOLO1 clade (136)", "DDE (176)", "BED (160)")
outer_cols <- c(
  "SOLO1 clade (136)" = COL_SOLO1,
  "DDE (176)"         = COL_DDE,
  "BED (160)"         = COL_BED
)

outer_track_df <- bind_rows(
  tip_df %>% filter(label %in% s_solo1$desc) %>%
    transmute(y = y, x = R_OUT_IN, xend = R_OUT_OUT, outer_cat = factor("SOLO1 clade (136)", levels = outer_levels)),
  tip_df %>% filter(label %in% all_dde_tips) %>%
    transmute(y = y, x = R_OUT_IN, xend = R_OUT_OUT, outer_cat = factor("DDE (176)", levels = outer_levels)),
  tip_df %>% filter(label %in% all_bed_tips) %>%
    transmute(y = y, x = R_OUT_IN, xend = R_OUT_OUT, outer_cat = factor("BED (160)", levels = outer_levels))
)

# ---- 6. Curated Reference NLR callouts + repulsion solver ----
ref_specs <- tibble::tribble(
  ~lbl,             ~tip_id,         ~col,      ~face,
  "ZAR1",           "ZAR1_CNL",      COL_TXT,   "plain",
  "Pit",            "Pit_CNL",       COL_TXT,   "plain",
  "Pm3b",           "Pm3b_CNL",      COL_TXT,   "plain",
  "Pii-1 (helper)", "Pii-1_CNL",     COL_TXT,   "plain",
  "Pii-2 (sensor)", "Pii-2_CNL",     COL_TXT,   "plain",
  "Tsn1",           "Tsn1_ONL",      COL_TXT,   "plain",
  "Rp1-D",          "Rp1-D_CNL",     COL_TXT,   "plain",
  "Xa1",            "Xa1_CONL",      COL_TXT,   "plain",
  "Yr5",            "Yr5_CONL",      COL_TXT,   "plain",
  "Piz-t",          "Piz-t_CNL",     COL_TXT,   "plain",
  "Pik-2 (helper)", "Pik-2_CNL",     COL_TXT,   "plain",
  "Pia-1 (helper)", "Pia-1_CNL",     COL_TXT,   "plain",
  "Rpg5",           "RPG5_CNLO",     COL_TXT,   "plain",
  "Pia-2 (sensor)", "Pia-2_CNLO",    COL_TXT,   "plain",
  "Pi-ta",          "Pi-ta_CNLO",    COL_TXT,   "plain",
  "YrU1",           "YrU1_OCNLO",    COL_TXT,   "plain",
  "Pik-1 (sensor)", "Pik-1_CONL",    COL_TXT,   "plain",
  "Lr10",           "Lr10_CNL",      COL_TXT,   "plain",
  "Sr35",           "Sr35_CNL",      COL_TXT,   "plain",
  "Sr33",           "Sr33_CNL",      COL_TXT,   "plain",
  "MLA10",          "MLA10_CNL",     COL_TXT,   "plain"
) %>%
  mutate(x0 = tip_x[tip_id], y0 = tip_y[tip_id])

stopifnot(!any(is.na(ref_specs$y0)))

extra_refs <- tibble::tibble(
  lbl  = c("RNL", "TNL", "SOLO1 clade"),
  tip_id = NA_character_,
  col  = c(COL_RPW8, COL_TIR, COL_SOLO1),
  face = c("bold", "bold", "bold"),
  x0   = NA_real_,
  y0   = c(s_rnl$y_mid, s_tnl$y_mid, s_solo1$y_mid)
)

callouts <- bind_rows(ref_specs, extra_refs) %>%
  arrange(y0)

min_sep_y <- 6.2 / 360.0 * Y_LIM_HI
ys <- callouts$y0
n_c <- length(ys)
req_sep <- rep(min_sep_y, n_c - 1)
for (i in seq_len(n_c - 1)) {
  if (callouts$lbl[i] == "RNL" && callouts$lbl[i + 1] == "TNL") req_sep[i] <- min_sep_y * 1.28
  if (callouts$lbl[i] == "MLA10" && grepl("SOLO1", callouts$lbl[i + 1])) req_sep[i] <- min_sep_y * 1.65
}
y_lo <- 20.0
y_hi <- n_tip - 220.0
orig_ys <- ys
for (iter in 1:500) {
  moved <- FALSE
  for (i in seq_len(n_c - 1)) {
    diff <- ys[i + 1] - ys[i]
    if (diff < req_sep[i]) {
      push <- 0.5 * (req_sep[i] - diff) + 0.5
      ys[i]     <- ys[i] - push
      ys[i + 1] <- ys[i + 1] + push
      moved <- TRUE
    }
  }
  ys <- 0.94 * ys + 0.06 * orig_ys
  ys <- pmin(pmax(ys, y_lo), y_hi)
  if (!moved) break
}
ys[1] <- max(ys[1], y_lo)
for (i in 2:n_c) {
  if (ys[i] - ys[i - 1] < req_sep[i - 1]) ys[i] <- ys[i - 1] + req_sep[i - 1]
}
if (ys[n_c] > y_hi) {
  ys[n_c] <- y_hi
  for (i in (n_c - 1):1) {
    if (ys[i + 1] - ys[i] < req_sep[i]) ys[i] <- ys[i + 1] - req_sep[i]
  }
}

outgroup_max_x <- max(tip_df$x[tip_df$label %in% outgroup_tips])
callouts <- callouts %>%
  mutate(
    y1 = ys,
    deg = y_to_cw_deg(y1),
    hjust_val = ifelse(deg > 0 & deg < 180, 0, 1),
    is_solo1 = grepl("SOLO1", lbl),
    is_out   = lbl %in% c("RNL", "TNL"),
    r_anchor = R_OUT_OUT + 0.04,
    r_end    = ifelse(is_out, max(R_LABEL, outgroup_max_x + 0.12), R_LABEL),
    r_txt    = ifelse(is_solo1, R_LABEL + 0.14, r_end + 0.10),
    sz_txt   = ifelse(is_solo1, 4.20, ifelse(face == "bold", 3.70, 3.45))
  )

solo1_arc_df <- data.frame(
  y_start = s_solo1$y_min - 15,
  y_end   = s_solo1$y_max + 15,
  x_arc   = 0.5 * (R_OUT_IN + R_OUT_OUT)
)

tn_df <- tip_df %>% filter(label %in% tn_tips)
tn_lbl_df <- data.frame(
  x = tn_max_x + 0.16,
  y = median(tn_df$y),
  lbl = "TN",
  hjust_val = ifelse(y_to_cw_deg(median(tn_df$y)) < 180, 0, 1)
)

tip_guides_df <- callouts %>% filter(!is.na(x0))

mic_lbl_df <- tibble::tibble(
  lbl = c("MIC1", "MIC2", "MIC3"),
  x   = c(R_RING_IN - 0.95, R_RING_IN - 0.65, R_RING_IN - 0.85),
  y   = c(11490, s_mic2$y_mid, 6225)
)

# ---- 7. Assemble ggtree plot ----
X_MIN <- -0.25
X_MAX <- R_LABEL + 0.26

p <- ggtree(tr_grp, aes(color = group), layout = "fan", open.angle = 20, size = 0.14) +
  geom_hilight(node = s_solo1$node, fill = COL_SOLO1, alpha = 0.14, extendto = R_OUT_OUT) +
  geom_hilight(node = s_mic1$node,  fill = COL_ID,    alpha = 0.085, extendto = R_RING_IN - 0.01) +
  geom_hilight(node = s_mic2$node,  fill = COL_ID,    alpha = 0.085, extendto = R_RING_IN - 0.01) +
  geom_hilight(node = s_mic3$node,  fill = COL_ID,    alpha = 0.085, extendto = R_RING_IN - 0.01) +
  scale_color_manual(values = branch_cols, guide = "none") +
  ggnewscale::new_scale_color() +
  geom_segment(
    data = tip_guides_df,
    aes(x = x0, xend = R_RING_IN, y = y0, yend = y0),
    color = "#5A6572", linewidth = 0.28, linetype = "dotted", alpha = 0.85, inherit.aes = FALSE
  ) +
  geom_point(
    data = tip_guides_df,
    aes(x = x0, y = y0),
    color = COL_TXT, size = 0.80, inherit.aes = FALSE
  ) +
  geom_segment(
    data = ring_df %>% filter(dom_cat == "Unclassified (2,842)"),
    aes(x = r_in, xend = r_out, y = y, yend = y, color = dom_cat),
    linewidth = 0.18, inherit.aes = FALSE
  ) +
  geom_segment(
    data = ring_df %>% filter(dom_cat == "CC without ID (10,265)"),
    aes(x = r_in, xend = r_out, y = y, yend = y, color = dom_cat),
    linewidth = 0.20, inherit.aes = FALSE
  ) +
  geom_segment(
    data = ring_df %>% filter(dom_cat == "CC with ID (1,071)"),
    aes(x = r_in, xend = r_out, y = y, yend = y, color = dom_cat),
    linewidth = 0.28, inherit.aes = FALSE
  ) +
  geom_segment(
    data = ring_df %>% filter(dom_cat == "TIR (146)"),
    aes(x = r_in, xend = r_out, y = y, yend = y, color = dom_cat),
    linewidth = 0.34, inherit.aes = FALSE
  ) +
  geom_segment(
    data = ring_df %>% filter(dom_cat == "RPW8 (40)"),
    aes(x = r_in, xend = r_out, y = y, yend = y, color = dom_cat),
    linewidth = 0.42, inherit.aes = FALSE
  ) +
  scale_color_manual(
    name = "N-terminal domain",
    values = ring_cols,
    breaks = dom_levels,
    drop = FALSE,
    guide = guide_legend(
      ncol = 1,
      order = 1,
      override.aes = list(linewidth = 4.2, alpha = 1)
    )
  ) +
  ggnewscale::new_scale_color() +
  geom_segment(
    data = outer_track_df,
    aes(x = x, xend = xend, y = y, yend = y, color = outer_cat),
    linewidth = 0.18, inherit.aes = FALSE
  ) +
  scale_color_manual(
    name = "Clade / ID",
    values = outer_cols,
    breaks = outer_levels,
    drop = FALSE,
    guide = guide_legend(
      ncol = 1,
      order = 2,
      override.aes = list(linewidth = 4.2, alpha = 1)
    )
  ) +
  geom_segment(
    data = data.frame(y = 1, yend = n_tip, r_in = R_RING_IN, r_out = R_RING_OUT),
    aes(x = r_in, xend = r_in, y = y, yend = yend),
    color = "#C3CBD6", linewidth = 0.28, inherit.aes = FALSE
  ) +
  geom_segment(
    data = data.frame(y = 1, yend = n_tip, r_in = R_RING_IN, r_out = R_RING_OUT),
    aes(x = r_out, xend = r_out, y = y, yend = yend),
    color = "#C3CBD6", linewidth = 0.28, inherit.aes = FALSE
  ) +
  geom_segment(
    data = solo1_arc_df,
    aes(x = x_arc, xend = x_arc, y = y_start, yend = y_end),
    color = COL_SOLO1, linewidth = 2.6, lineend = "butt", inherit.aes = FALSE
  ) +
  geom_text(
    data = mic_lbl_df,
    aes(x = x, y = y, label = lbl),
    color = COL_ID, fontface = "bold", family = "Helvetica", size = 3.9, inherit.aes = FALSE
  ) +
  geom_point(
    data = callouts %>% filter(!is_solo1),
    aes(x = r_anchor, y = y0),
    color = callouts$col[!callouts$is_solo1], size = 0.85, inherit.aes = FALSE
  ) +
  geom_segment(
    data = callouts %>% filter(!is_solo1),
    aes(x = r_anchor, xend = R_ELBOW, y = y0, yend = y0),
    color = callouts$col[!callouts$is_solo1], linewidth = 0.30, inherit.aes = FALSE
  ) +
  geom_segment(
    data = callouts %>% filter(!is_solo1),
    aes(x = R_ELBOW, xend = r_end, y = y0, yend = y1),
    color = callouts$col[!callouts$is_solo1], linewidth = 0.30, inherit.aes = FALSE
  ) +
  geom_segment(
    data = callouts %>% filter(is_solo1),
    aes(x = r_anchor, xend = R_LABEL, y = y0, yend = y1),
    color = COL_SOLO1, linewidth = 0.65, inherit.aes = FALSE
  ) +
  geom_text(
    data = callouts,
    aes(x = r_txt, y = y1, label = lbl, hjust = hjust_val, fontface = face, size = sz_txt),
    color = callouts$col, family = "Helvetica", vjust = 0.5, inherit.aes = FALSE
  ) +
  scale_size_identity() +
  geom_text(
    data = tn_lbl_df,
    aes(x = x, y = y, label = lbl, hjust = hjust_val),
    color = COL_TIR, fontface = "bold", family = "Helvetica", size = 3.6, vjust = 0.5, inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(X_MIN, X_MAX), expand = c(0, 0), oob = scales::oob_keep) +
  coord_polar(theta = "y", start = START_RAD, direction = 1, clip = "off") +
  theme_void(base_family = "Helvetica") +
  theme(
    text                  = element_text(family = "Helvetica"),
    plot.background       = element_rect(fill = "white", color = NA),
    plot.margin           = margin(12, 136, 12, 36),
    legend.position       = c(0.95, 0.01),
    legend.justification  = c(0, 0),
    legend.box            = "vertical",
    legend.box.background = element_rect(fill = "#FFFFFF", color = "#CBD2DC", linewidth = 0.35),
    legend.box.margin     = margin(6, 9, 6, 9),
    legend.background     = element_blank(),
    legend.margin         = margin(2, 0, 2, 0),
    legend.spacing.y      = unit(0.20, "cm"),
    legend.title          = element_text(size = 10.0, face = "bold", color = COL_TXT, margin = margin(b = 4)),
    legend.text           = element_text(size = 9.2, color = COL_TXT),
    legend.key            = element_rect(fill = NA, color = NA),
    legend.key.height     = unit(0.44, "cm"),
    legend.key.width      = unit(0.44, "cm")
  )

draw_plot_with_scalebar <- function(p_obj, x_min, x_max) {
  print(p_obj)
  len_npc <- 0.4 / (x_max - x_min)
  vp_txt <- paste(capture.output(grid::current.vpTree()), collapse = " ")
  vp_name <- regmatches(vp_txt, regexpr("panel\\.[0-9-]+", vp_txt))
  stopifnot(length(vp_name) == 1L && nzchar(vp_name))
  grid::seekViewport(vp_name)
  x_mid <- 1.09
  y_sb  <- 0.430
  x0 <- x_mid - 0.5 * len_npc
  x1 <- x_mid + 0.5 * len_npc
  grid::grid.segments(x0 = unit(x0, "npc"), x1 = unit(x1, "npc"),
                      y0 = unit(y_sb, "npc"), y1 = unit(y_sb, "npc"),
                      gp = gpar(col = COL_TXT, lwd = 1.3))
  grid::grid.segments(x0 = unit(c(x0, x1), "npc"), x1 = unit(c(x0, x1), "npc"),
                      y0 = unit(y_sb - 0.009, "npc"), y1 = unit(y_sb + 0.009, "npc"),
                      gp = gpar(col = COL_TXT, lwd = 1.1))
  grid::grid.text("1.0", x = unit(x_mid, "npc"), y = unit(y_sb + 0.017, "npc"),
                  just = c("center", "bottom"),
                  gp = gpar(col = COL_TXT, fontsize = 9.6, fontfamily = "Helvetica"))
  grid::upViewport(0)
}

cat("Rendering PNG and PDF...\n")
agg_png(out_png, width = opt$width, height = opt$height, units = "in", res = opt$dpi, background = "white")
draw_plot_with_scalebar(p, X_MIN, X_MAX)
invisible(dev.off())

grDevices::pdf(out_pdf, width = opt$width, height = opt$height, bg = "white", family = "Helvetica", useDingbats = FALSE)
draw_plot_with_scalebar(p, X_MIN, X_MAX)
invisible(dev.off())

cat("Saved:\n ", out_png, "\n ", out_pdf, "\n")
