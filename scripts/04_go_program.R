# ===============================================================================
# 04 — GO ENRICHMENT OF THE PROGRAM + PATHWAY PRESERVATION
# ===============================================================================
# Derived from the stage-1 program and the direct pulse-vs-continuous contrast:
#   1. GO enrichment of the continuous IL-4 program, induced and repressed, drawn
#      as per-arm bars, with each term's not-below-continuous fraction alongside
#      (MAIN figure). A dot-plot rendering of the same terms is also written.
#   2. Pathway preservation across more terms, as a supplementary panel and a CSV.
#   3. Per-pathway amplitude, as a supplementary uniformity panel.
#
# Depends only on deg_program.csv and master_expressed_genes.csv. Shared
# enrichment machinery is in scripts/_go_helpers.R.
#
# Input : data/deg_program.csv, data/master_expressed_genes.csv
# Output: results/go_program_main.pdf             (MAIN)
#         results/go_program_dotplot_alt.pdf      (same terms, dot-plot rendering)
#         results/go_pathway_preservation.pdf     (supplement, more terms)
#         results/go_pathway_preservation.csv
#         results/amplitude_by_pathway_supp.pdf   (supplement)
#         results/amplitude_by_pathway.csv
#         results/go_program_results.xlsx
# ===============================================================================

source(file.path("scripts", "_go_helpers.R"))

prog <- read.csv(file.path(data_dir, "deg_program.csv"), stringsAsFactors = FALSE)
psig <- function(p) !is.na(p) & p < 0.05
prog_groups <- list(
  # Sign + padj only. padj already tests H0: |log2FC| <= log2(1.5) (set in 01), so
  # an extra > 1 point-estimate gate would reinstate exactly the mismatch 01 removed.
  "Program: up"   = prog$gene_symbol[(prog$log2fc_c100_4h >  0 & psig(prog$padj_c100_4h)) |
                                     (prog$log2fc_c100_24h >  0 & psig(prog$padj_c100_24h))],
  "Program: down" = prog$gene_symbol[(prog$log2fc_c100_4h <  0 & psig(prog$padj_c100_4h)) |
                                     (prog$log2fc_c100_24h <  0 & psig(prog$padj_c100_24h))]
)
cat("\n[04] GO enrichment — stage-1 continuous program:\n")
prog_ego <- lapply(names(prog_groups), function(nm) run_enrich(prog_groups[[nm]], nm))
names(prog_ego) <- names(prog_groups)

cat("\n[04] term selection:\n")
prog_terms <- lapply(names(prog_ego), function(nm) { cat(sprintf("  %s\n", nm)); select_terms(prog_ego[[nm]]) })
names(prog_terms) <- names(prog_ego)

# The dot plot is an alternative rendering of the same prog_terms as the main bar
# figure below, so the two cannot disagree about which pathways are shown.
dotplot_panel(prog_ego, prog_terms,
              "GO enrichment of the continuous IL-4 program",
              paste(COLLAPSE_NOTE,
                    "Stage-1 program (>1.5-fold, padj < 0.05); induced and repressed separately", sep = "\n"),
              file.path(out_dir, "go_program_dotplot_alt.pdf"), 8.2, 6.0)

# ---- MAIN FIGURE: per-arm bars, preservation printed at the bar end -----------
# Built after pres_all below is available; see the bar_panel call further down.

# ===============================================================================
# PATHWAY-LEVEL SHORTFALL — is any pathway selectively under-reproduced?
# ===============================================================================
# Deliberately NOT a two-condition dot plot. Two reasons:
#  1. SELECTION ASYMMETRY. Picking top-N per condition independently makes both
#     columns look equally full whether or not anything is preserved, because a
#     continuous term that dies under pulse simply never becomes a row. Terms here
#     are anchored on the CONTINUOUS program and the pulse is measured against them.
#  2. p-VALUES AND COUNTS ARE NOT COMPARABLE across separately-run enrichments with
#     different gene-set sizes -- the pulse sets are much smaller, so a smaller padj
#     under pulse can be pure power. Worse, counts can agree while the GENES differ:
#     a term can carry a similar number of continuous and pulse genes while sharing
#     only half of them, and so look preserved when it is not.
# So the metric is the fraction of a term's continuous-program genes whose pulse
# response is NOT SIGNIFICANTLY BELOW continuous at 24 h, judged by the DIRECT
# pulse-vs-continuous contrast. Immune to set-size and to count coincidences.
#
# WHAT THIS IS NOT. It is the complement of a significance test, not an equivalence
# test: a gene counts toward it because the difference was not detected, which
# reflects statistical power as well as similarity. Do not describe it as genes the
# pulse "reproduced" or "preserved" -- those words assert equivalence the test
# cannot establish. The Methods state this limitation explicitly.
#
# The dashed line is the program-wide rate, so a pathway to its right is above the
# program average. HIGHER IS BETTER and terms sort accordingly. Do not invert this
# to a shortfall rate without also changing the axis label and the sort direction.
# The main figure carries the rate as an annotation column; the supplementary panel
# extends it to the top PRES_SUPP_N terms per arm and the full table goes to CSV.
# Ordering is always by continuous significance, never by preservation.
PRES_SUPP_N <- 15
pulse_groups <- list(
  "Program: up"   = prog$gene_symbol[(prog$log2fc_p100_4h >  0 & psig(prog$padj_p100_4h)) |
                                     (prog$log2fc_p100_24h >  0 & psig(prog$padj_p100_24h))],
  "Program: down" = prog$gene_symbol[(prog$log2fc_p100_4h <  0 & psig(prog$padj_p100_4h)) |
                                     (prog$log2fc_p100_24h <  0 & psig(prog$padj_p100_24h))]
)
cat("\n[04] pulse program enrichment (only for the 'still enriched' flag):\n")
pulse_ego <- lapply(names(pulse_groups), function(nm) run_enrich(pulse_groups[[nm]], paste("pulse", nm)))
names(pulse_ego) <- names(pulse_groups)

# Computed for EVERY collapsed continuous term, not just the displayed ones.
#
# CAUGHT-UP GENES COUNT AS PRESERVED. The three bands from scripts/_common.R are
# matched, caught_up and persistent_shortfall; a gene that reaches continuous
# levels by 24 h after a delay has not fallen short, so the first two both belong in
# the numerator. Equivalently: the numerator is simply the genes NOT significantly
# below continuous at 24 h -- the 4 h test only decides which of the two labels a
# gene carries, never the count.
#
# WHY THIS PANEL AND NOT A DIRECT-COMPARISON ORA. Running enrichGO with the union of
# shortfall and matched as the universe is rigorous but returns only BROAD terms
# ("signalling", "developmental process"): in a universe that small a specific term
# retains too few members to reach significance, and lowering minGSSize does not
# help because the limit is power, not filtering. Anchoring on terms enriched in the
# FULL program keeps their membership intact, so the terms stay specific, and the
# preservation rate supplies the rigour.
#
# MIN_CLASSIFIED is a precision floor: suppress rates computed on too few genes.
# NB with T = 0 every program gene falls into one of the three bands, so
# n_classified equals the term's gene count and this acts as a term-size floor.
MIN_CLASSIFIED <- 15

# Per-gene 24h amplitude ratio (pulse / continuous), the threshold-free quantity
# script 03 regresses. Carried into the CSV as med_ratio so the fate-based % can be
# checked against it; the panel itself still plots the fate-based %.
.mg <- read.csv(file.path(data_dir, "master_expressed_genes.csv"), stringsAsFactors = FALSE)
.ratio_of <- function(syms) {
  i <- match(syms, .mg$gene_symbol); i <- i[!is.na(i)]
  cc <- .mg$log2fc_c100_24h[i]; pp <- .mg$log2fc_p100_24h[i]
  k <- !is.na(cc) & !is.na(pp) & abs(cc) > 0.585      # avoid dividing by ~0
  if (!any(k)) return(NA_real_)
  median(pmax(pmin(pp[k] / cc[k], 2), -1), na.rm = TRUE)   # clamp wild ratios
}
fb_pres <- pulse_fate_bands(data_dir)
pres_all <- do.call(rbind, lapply(names(prog_ego), function(nm) {
  cont <- prune_redundant(prog_ego[[nm]], OVERLAP_CUT)
  dir  <- ifelse(grepl("up", nm), "Up", "Down")
  sh   <- fb_pres$gene_symbol[fb_pres$direction == dir & fb_pres$band == "persistent_shortfall"]
  mt   <- fb_pres$gene_symbol[fb_pres$direction == dir & fb_pres$band == "matched"]
  cu   <- fb_pres$gene_symbol[fb_pres$direction == dir & fb_pres$band == "caught_up"]
  r <- do.call(rbind, lapply(seq_len(nrow(cont)), function(i) {
    gs <- strsplit(cont$geneID[i], "/")[[1]]
    a <- length(intersect(gs, sh))                       # persistent shortfall
    b <- length(intersect(gs, mt)) + length(intersect(gs, cu))   # matched + caught up
    data.frame(direction = sub("Program: ", "", nm), Term = cont$Description[i],
               p_cont = cont$p.adjust[i], n_classified = a + b,
               n_preserved = b, n_shortfall = a,
               pct = if (a + b > 0) 100 * b / (a + b) else NA_real_,
               med_ratio = .ratio_of(gs),
               # DATA ONLY -- carried in the CSV and workbook, not drawn. It
               # compares two separately-run enrichments of different-sized gene
               # sets, so it largely tracks power rather than preservation.
               sig_in_pulse = cont$Description[i] %in% pulse_ego[[nm]]$Description,
               stringsAsFactors = FALSE)
  }))
  # report what the precision floor costs, so it is never a silent truncation
  ok <- !is.na(r$pct) & r$n_classified >= MIN_CLASSIFIED
  cat(sprintf("  %-14s %3d collapsed terms -> %3d kept at n_classified >= %d (%d dropped; median n dropped = %.0f)\n",
              nm, nrow(r), sum(ok), MIN_CLASSIFIED, sum(!ok),
              if (any(!ok)) median(r$n_classified[!ok]) else NA_real_))
  r[ok, ]
}))
# reference line = the program-wide shortfall rate, per direction
baseline <- do.call(rbind, lapply(c("Up", "Down"), function(dir) {
  a <- sum(fb_pres$direction == dir & fb_pres$band == "persistent_shortfall")
  b <- sum(fb_pres$direction == dir & fb_pres$band %in% c("matched", "caught_up"))
  # b / (a+b), NOT a / (a+b): `pct` counts genes NOT below continuous, so the
  # reference line must be that same rate or it lands on the mirror-image value.
  data.frame(direction = ifelse(dir == "Up", "up", "down"), baseline = 100 * b / (a + b),
             stringsAsFactors = FALSE)
}))
write.csv(pres_all[order(pres_all$direction, pres_all$p_cont), ],
          file.path(out_dir, "go_pathway_preservation.csv"), row.names = FALSE)
cat(sprintf("\nSaved: results/go_pathway_preservation.csv (%d terms — full supplementary table)\n",
            nrow(pres_all)))
print(baseline, row.names = FALSE)
# THE NUMBERS THE MANUSCRIPT QUOTES ARE PRINTED HERE, so there is one authoritative
# source for them. Take them from this line rather than re-deriving them by hand;
# an earlier draft quoted values no configuration of this pipeline reproduces.
cat("\n  --- summary for the Results text ---\n")
for (dd in c("up", "down")) {
  s <- pres_all$pct[pres_all$direction == dd]
  cat(sprintf("  %-5s %3d enriched terms | MEDIAN %.0f%% not significantly below continuous | range %.0f-%.0f%% | program-wide %.0f%%\n",
              dd, length(s), median(s), min(s), max(s),
              baseline$baseline[baseline$direction == dd]))
}

LAB <- c(up = "Up (activated)", down = "Down (repressed)")

# ---- MAIN FIGURE ---------------------------------------------------------------
# Per-arm bars with the not-below-continuous fraction at the bar end, which does
# the work of a shared-axis dot plot and a separate preservation panel at once:
#  - per-arm faceting avoids the ordering artefact a shared axis forces, where a
#    term must be placed by its best padj across arms and neither column reads in
#    its own order;
#  - the preservation rates are tightly clustered within each arm, so no pathway is
#    selectively preserved; printing them in a column shows this without a second
#    panel.
# The annotation is the preservation RATE, not "still enriched under pulse": that
# flag splits terms into groups whose preservation is statistically
# indistinguishable, so drawing it would imply a difference the genes do not
# support.
bar_df <- do.call(rbind, lapply(names(prog_terms), function(nm) {
  dir <- sub("Program: ", "", nm)
  e   <- prog_ego[[nm]]; tms <- prog_terms[[nm]]
  e   <- e[match(tms, e$Description), ]
  pr  <- pres_all[pres_all$direction == dir, ]
  i   <- match(tms, pr$Term)
  data.frame(direction = LAB[[dir]], Term = tms,
             neglogp = -log10(e$p.adjust),
             label = ifelse(is.na(i), sprintf("%d genes", e$Count),
                            sprintf("%.0f%% (%d/%d)", pr$pct[i], pr$n_preserved[i], pr$n_classified[i])),
             stringsAsFactors = FALSE)
}))
bar_df$direction <- factor(bar_df$direction, levels = LAB)
bar_panel(bar_df, file.path(out_dir, "go_program_main.pdf"),
          width = 7.6, height = 5.6, base_size = 10,
          title = "GO enrichment of the continuous IL-4 program",
          subtitle = paste(COLLAPSE_NOTE,
            "Printed beside each bar: % of that pathway's genes not significantly below continuous at 24 h.",
            sep = "\n"))


preservation_panel <- function(terms_by_dir, file, width, height, base_size, title, subtitle,
                               order_ref = NULL) {
  # TERMS ARE FIXED BY CONTINUOUS ENRICHMENT, NOT BY PRESERVATION. Do not rank this
  # panel by `pct`: ranking by a proportion preferentially surfaces small, noisy
  # terms, and it would show a DIFFERENT pathway list from the enrichment panel
  # beside it, so the two could no longer be read together. The caller passes the
  # term list; this only reports the rate for each, in continuous-significance
  # order.
  d <- do.call(rbind, lapply(names(terms_by_dir), function(dd) {
    s <- pres_all[pres_all$direction == dd &
                  pres_all$Term %in% terms_by_dir[[dd]], ]
    miss <- setdiff(terms_by_dir[[dd]], s$Term)
    if (length(miss))
      cat(sprintf("  [pres] %-5s dropped below the n>=%d floor: %s\n",
                  dd, MIN_CLASSIFIED, paste(miss, collapse = "; ")))
    s[order(s$p_cont), ]
  }))
  b <- baseline
  d$direction <- factor(LAB[d$direction], levels = LAB)
  b$direction <- factor(LAB[b$direction], levels = LAB)
  # terms repeat across facets (chemotaxis is top in both), so key on direction
  # Row height comes from `order_ref` when supplied (the enrichment panel's own
  # level order) so the two figures are row-for-row identical; otherwise fall back
  # to continuous significance within each facet.
  rank_of <- if (is.null(order_ref)) -d$p_cont else match(d$Term, order_ref)
  d$key <- factor(paste0(as.integer(d$direction), "|", d$Term),
                  levels = paste0(as.integer(d$direction), "|", d$Term)[order(d$direction, rank_of)])
  p <- ggplot(d, aes(x = pct, y = key)) +
    geom_vline(data = b, aes(xintercept = baseline), linetype = "dashed", color = "grey45") +
    geom_segment(aes(x = 0, xend = pct, yend = key), color = "grey78", linewidth = 0.5) +
    geom_point(shape = 21, size = 3.2, fill = "#2C7FB8", colour = "grey20", stroke = 0.4) +
    geom_text(aes(label = sprintf("%d/%d", n_preserved, n_classified)), hjust = -0.4,
              size = base_size * 0.24, colour = "grey35") +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^[0-9]+\\|", "", x)) +
    scale_x_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, subtitle = subtitle, x = "% of pathway genes not significantly below continuous", y = NULL) +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          plot.title = element_text(size = base_size * 0.95),
          plot.subtitle = element_text(size = base_size * 0.72, colour = "grey35"),
          panel.grid.major.y = element_blank())
  for (ff in levels(d$direction)) {
    lv <- rev(levels(d$key)); lv <- lv[lv %in% d$key[d$direction == ff]]
    cat(sprintf("  [yorder] %s [%s] : %s\n", basename(file), ff,
                paste(sub("^[0-9]+\\|", "", lv), collapse = " | ")))
  }
  ggsave(file, p, width = width, height = height)
  cat(sprintf("Saved: %s\n", file))
}

FULL_NOTE <- paste("Pathways enriched in the continuous program. Point = % of that pathway's",
                   "genes whose pulse response is NOT significantly below continuous at 24 h,",
                   "by the direct pulse-vs-continuous contrast (n / n classified). This is the",
                   "complement of a significance test, not an equivalence test. Dashed line =",
                   "program-wide rate.", sep = "\n")

pres_supp_terms <- lapply(c(up = "Program: up", down = "Program: down"), function(nm) {
  s <- pres_all[pres_all$direction == sub("Program: ", "", nm), ]
  head(s$Term[order(s$p_cont)], PRES_SUPP_N)
})
prog_order <- term_order(prog_ego, prog_terms)
# supplement — more terms, full detail. A main-figure preservation panel is not
# drawn: its content is the annotation column of the bar figure above.
preservation_panel(pres_supp_terms, file.path(out_dir, "go_pathway_preservation.pdf"),
                   order_ref = prog_order,
                   width = 7.5, height = 8, base_size = 10,
                   title = "Pathway-level reproduction of the continuous IL-4 program",
                   subtitle = FULL_NOTE)

# Pulse-program enrichment is exported but never plotted: adjusted p-values are not
# comparable across enrichments of different-sized gene sets.
names(pulse_ego) <- c("Pulse: up", "Pulse: down")
# ===============================================================================
# SUPPLEMENT — amplitude reproduction is uniform across pathways
# ===============================================================================
# One point per enriched GO term: the median amplitude ratio (pulse/continuous
# log2FC at 24 h) of that term's genes. The result is the ABSENCE of structure --
# the observed SD of term medians equals the SD expected from sampling alone given
# the median term size, so excess between-pathway variance is negligible and the
# tails are what that sampling noise looks like, not pathway biology.
#
# Metric is the RATIO, not the preservation percentage: preservation is a
# thresholded call from the direct contrast and is strongly confounded by expression
# level, whereas the ratio is continuous and unthresholded. Genes with
# |continuous log2FC| <= 0.585 are excluded, since the denominator is then near zero.
#
# NB the amplitude ratio correlates with the size of the continuous response, so a
# gene set sitting low on that axis can look selectively deficient without being so.
# Script 10 tests the EGR2 result against this confound.
.mg2 <- .mg[is_annotatable_gene(.mg$ensembl_gene_id, data_dir), ]
gene_ratio_tab <- function(dir) {
  sgn <- if (dir == "up") 1 else -1
  k <- !is.na(.mg2$padj_c100_24h) & .mg2$padj_c100_24h < 0.05 &
       sign(.mg2$log2fc_c100_24h) == sgn & abs(.mg2$log2fc_c100_24h) > 0.585
  g <- .mg2[k, ]
  data.frame(direction = LAB[[dir]], ratio = g$log2fc_p100_24h / g$log2fc_c100_24h,
             stringsAsFactors = FALSE)
}
gr <- do.call(rbind, lapply(c("up", "down"), gene_ratio_tab))
gr <- gr[is.finite(gr$ratio), ]; gr$direction <- factor(gr$direction, levels = LAB)
base_med <- tapply(gr$ratio, gr$direction, median)

term_tab <- do.call(rbind, lapply(names(prog_ego), function(nm) {
  dir  <- sub("Program: ", "", nm)
  cont <- prune_redundant(prog_ego[[nm]], OVERLAP_CUT)
  do.call(rbind, lapply(seq_len(nrow(cont)), function(i) {
    gs <- strsplit(cont$geneID[i], "/")[[1]]
    if (length(gs) < 15) return(NULL)
    data.frame(direction = LAB[[dir]], Term = cont$Description[i], n = length(gs),
               med_ratio = .ratio_of(gs), stringsAsFactors = FALSE)
  }))
}))
term_tab <- term_tab[!is.na(term_tab$med_ratio), ]
term_tab$direction <- factor(term_tab$direction, levels = LAB)

# ---- Is the term-to-term spread larger than chance? ---------------------------
# PERMUTATION, not a closed-form SE. The statistic is the SD of term medians. The
# null holds each term's gene MEMBERSHIP fixed -- so term sizes and the (heavy)
# overlap between terms are preserved exactly -- and permutes the gene-level
# amplitude ratios. Permutation is done WITHIN quintiles of continuous-response
# magnitude, so the null retains the ratio-vs-magnitude relationship and excess
# spread cannot be manufactured by pathways differing in magnitude composition.
#
# A large p alone would only say "not detected". The reported detection limit is
# what makes this a positive statement: it is the smallest true between-term SD
# that this design would flag with >=80% power, obtained by adding N(0, delta)
# offsets to the term medians under the null.
PERM_B  <- 4000
POW_B   <- 300
set.seed(1)

uniformity <- function(dir) {
  sgn <- if (dir == "up") 1 else -1
  k <- !is.na(.mg2$padj_c100_24h) & .mg2$padj_c100_24h < 0.05 &
       sign(.mg2$log2fc_c100_24h) == sgn & abs(.mg2$log2fc_c100_24h) > 0.585
  g <- .mg2[k, ]
  g$ratio <- pmax(pmin(g$log2fc_p100_24h / g$log2fc_c100_24h, 2), -1)
  g <- g[is.finite(g$ratio), ]

  cont <- prune_redundant(prog_ego[[paste0("Program: ", dir)]], OVERLAP_CUT)
  idx  <- lapply(strsplit(cont$geneID, "/"), function(s) match(intersect(s, g$gene_symbol), g$gene_symbol))
  idx  <- idx[lengths(idx) >= 15]
  med_of <- function(v) vapply(idx, function(i) median(v[i]), 0)

  strat <- cut(abs(g$log2fc_c100_24h), quantile(abs(g$log2fc_c100_24h), 0:5/5), include.lowest = TRUE)
  perm  <- function() { v <- g$ratio
    for (s in levels(strat)) { j <- which(strat == s); v[j] <- sample(v[j]) }; v }

  obs  <- sd(med_of(g$ratio))
  null <- replicate(PERM_B, sd(med_of(perm())))
  q95  <- quantile(null, 0.95)
  det  <- NA_real_
  for (delta in c(0.02, 0.03, 0.04, 0.05, 0.075, 0.10)) {
    if (mean(replicate(POW_B, sd(med_of(perm()) + rnorm(length(idx), 0, delta))) > q95) >= 0.8) {
      det <- delta; break }
  }
  data.frame(direction = LAB[[dir]], n_terms = length(idx), obs_sd = obs,
             null_mean = mean(null), p_perm = mean(null >= obs),
             excess_sd = sqrt(max(0, obs^2 - mean(null^2))), detect_sd = det,
             stringsAsFactors = FALSE)
}
unif <- do.call(rbind, lapply(c("up", "down"), uniformity))
unif$direction <- factor(unif$direction, levels = LAB)
write.csv(unif, file.path(out_dir, "amplitude_uniformity_test.csv"), row.names = FALSE)
cat("\nSaved: results/amplitude_uniformity_test.csv\n")
print(unif, row.names = FALSE, digits = 3)

note <- data.frame(direction = unif$direction,
                   label = sprintf("%d terms | SD %.3f vs %.3f under a permutation null (p = %.2f)\nexcludes between-pathway SD > %.2f",
                                   unif$n_terms, unif$obs_sd, unif$null_mean, unif$p_perm, unif$detect_sd),
                   stringsAsFactors = FALSE)

pa <- ggplot(term_tab, aes(med_ratio, y = 0)) +
  geom_vline(data = data.frame(direction = factor(names(base_med), levels = LAB), b = base_med),
             aes(xintercept = b), linetype = "dashed", colour = "grey45") +
  geom_jitter(height = 0.28, width = 0, alpha = 0.75, size = 1.9, colour = "#34495E") +
  geom_text(data = note, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.03, vjust = 1.8, size = 2.5, colour = "grey30") +
  facet_wrap(~ direction, ncol = 1) +
  scale_y_continuous(limits = c(-0.55, 0.75), breaks = NULL) +
  labs(x = "median amplitude ratio of the pathway's genes (pulse / continuous, 24 h)", y = NULL,
       title = "No pathway is selectively under-reproduced",
       subtitle = paste("One point per enriched GO term (>=15 genes). Dashed line = program-wide median.",
                        "Term-to-term spread does not exceed a permutation null that preserves term size,",
                        "term overlap and the ratio-vs-magnitude relationship.", sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(out_dir, "amplitude_by_pathway_supp.pdf"), pa, width = 6.6, height = 5.2)
cat("Saved: results/amplitude_by_pathway_supp.pdf\n")
write.csv(term_tab[order(term_tab$direction, term_tab$med_ratio), ],
          file.path(out_dir, "amplitude_by_pathway.csv"), row.names = FALSE)
cat("Saved: results/amplitude_by_pathway.csv\n")

all_ego <- c(prog_ego, pulse_ego)
# in_figure must mean "DRAWN IN THIS GROUP'S COLUMN", which is not the same as
# "selected for this group's top-N list": a term selected for one group's top-N is
# also drawn wherever else it is significant. So flag against the UNION of the
# panel's term lists, not the per-group list. Each sheet already contains only terms
# significant in that group, so membership in the union is exactly "drawn here".
panel_of <- list(prog_terms)
panel_union <- function(nm) {
  for (pl in panel_of) if (nm %in% names(pl)) return(unique(unlist(pl)))
  character()
}
wb <- createWorkbook()
for (nm in names(all_ego)) {
  sn <- substr(gsub("[^A-Za-z0-9 _]", "_", nm), 1, 31)   # Excel-safe sheet name
  # ALL significant terms are written out, with a flag marking the ones the figure
  # shows, so the redundancy collapse is auditable.
  out <- all_ego[[nm]]
  if (nrow(out)) out$in_figure <- out$Description %in% panel_union(nm)
  addWorksheet(wb, sn); writeData(wb, sn, out)
}
# Preservation table also goes in the workbook, so the xlsx is self-contained.
addWorksheet(wb, "Preservation")
writeData(wb, "Preservation", pres_all[order(pres_all$direction, pres_all$p_cont), ])
saveWorkbook(wb, file.path(out_dir, "go_program_results.xlsx"), overwrite = TRUE)
cat("Saved: results/go_program_results.xlsx\n")
