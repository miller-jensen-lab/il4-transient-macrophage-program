# ===============================================================================
# 05 — TF REGULONS AND PULSE AMPLITUDE
# ===============================================================================
# Does the pulse under-deliver the targets of any particular transcription factor?
# For every CollecTRI regulon with enough targets in the program, this compares the
# amplitude ratio (pulse/continuous log2FC at 24 h) of its targets with the rest of
# that arm. Individual gene ratios are noise-dominated; a mean over tens of targets
# is not, which is what makes a regulon-level test viable where a per-gene
# classification is not.
#
# The result is a negative, and its power limit must be reported with it. The median
# regulon is small enough that BH significance requires a deviation larger than the
# EGR2 effect reported in 06, so this scan cannot resolve an effect the size of that
# mechanism. It excludes strong network-level selectivity, not modest effects. The
# per-regulon confidence intervals are more informative than the scan-level verdict
# and are written to tf_regulon_amplitude.csv.
#
# Two reasons not to interpret a single top hit: the top-ranked regulon of a large
# scan is inflated by winner's curse, and CollecTRI regulons are not independent, so
# related factors sharing most of their targets are one observation, not two.
#
# Related: pathway-level uniformity in 04; EGR2 dependence, tested against external
# knockout data rather than a network inference, in 06.
#
# Input : data/master_expressed_genes.csv, data/collectri_human_network.csv
# Output: results/tf_regulon_amplitude.pdf, results/tf_regulon_amplitude.csv,
#         results/tf_enrichment_results.xlsx
# ===============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(openxlsx); library(ggrepel) })
source(file.path("scripts", "_common.R"))

MIN_TARGETS <- 15     # below this a regulon mean is too noisy to be informative

# SD of the per-gene amplitude ratio across the program, used only to size the
# resolution band on the figure: the SE of a regulon mean over n targets is
# RATIO_SD/sqrt(n). Measured from the ratio distribution in 03; recompute if the
# gene universe changes. BH_SE_MULT is the multiple of that SE a typical regulon
# must exceed to survive BH across all regulons tested.
RATIO_SD    <- 0.17
BH_SE_MULT  <- 3.3
PADJ        <- 0.05

data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

m <- read.csv(file.path(data_dir, "master_expressed_genes.csv"), stringsAsFactors = FALSE)
m <- m[is_annotatable_gene(m$ensembl_gene_id, data_dir), ]
m$SYM <- toupper(m$gene_symbol)

# CollecTRI is a HUMAN network; symbols are matched case-insensitively, standard for
# conserved TF-target relationships but an assumption worth stating in the Methods.
net <- read.csv(file.path(data_dir, "collectri_human_network.csv"), stringsAsFactors = FALSE)
names(net)[1:2] <- c("source", "target")
net$source <- toupper(net$source); net$target <- toupper(net$target)

arm_table <- function(arm) {
  sgn <- if (arm == "induced") 1 else -1
  k <- !is.na(m$padj_c100_24h) & m$padj_c100_24h < PADJ &
       sign(m$log2fc_c100_24h) == sgn & abs(m$log2fc_c100_24h) > 1
  d <- m[k, ]
  d$ratio <- d$log2fc_p100_24h / d$log2fc_c100_24h
  d[is.finite(d$ratio), ]
}

res <- do.call(rbind, lapply(c("induced", "repressed"), function(arm) {
  d <- arm_table(arm)
  cat(sprintf("[05] %-9s arm: %d genes, program mean ratio %.3f\n", arm, nrow(d), mean(d$ratio)))
  out <- do.call(rbind, lapply(unique(net$source), function(tf) {
    i <- d$SYM %in% net$target[net$source == tf]
    if (sum(i) < MIN_TARGETS) return(NULL)
    tt <- t.test(d$ratio[i], d$ratio[!i])
    # SCALE. t.test estimates mean(targets) - mean(NON-targets), but `delta` below
    # is mean(targets) - mean(ALL genes in the arm). These are not the same
    # quantity; they differ by a constant factor:
    #   delta = (n_non / n_all) * [mean(targets) - mean(non-targets)]
    # so the interval is rescaled onto delta's scale. Left unscaled it is up to
    # ~25% too wide for the largest regulons and describes something other than
    # the effect printed beside it. The p-value is unaffected by the rescaling.
    ci <- tt$conf.int * ((nrow(d) - sum(i)) / nrow(d))
    # CI retained: the plausible range for this regulon's deviation, which is more
    # informative than the BH verdict (see header).
    data.frame(arm = arm, tf = tf, n_targets = sum(i),
               regulon_mean = mean(d$ratio[i]), program_mean = mean(d$ratio),
               delta = mean(d$ratio[i]) - mean(d$ratio),
               ci_lo = ci[1], ci_hi = ci[2],
               excludes_egr2_sized = max(abs(ci)) < 0.055,
               p = tt$p.value, stringsAsFactors = FALSE)
  }))
  out$padj <- p.adjust(out$p, "BH")
  out[order(out$p), ]
}))
write.csv(res, file.path(out_dir, "tf_regulon_amplitude.csv"), row.names = FALSE)
cat("Saved: results/tf_regulon_amplitude.csv\n")
for (a in unique(res$arm)) {
  s <- res[res$arm == a, ]
  cat(sprintf("  %-9s %d regulons (>=%d targets); significant after BH: %d; largest |delta| %.3f (%s)\n",
              a, nrow(s), MIN_TARGETS, sum(s$padj < PADJ),
              max(abs(s$delta)), s$tf[which.max(abs(s$delta))]))
}
key <- res[res$tf %in% c("NFKB1", "RELA", "STAT6", "SPI1", "STAT1", "PPARG"), ]
if (nrow(key)) { cat("\n  regulons of interest:\n")
  print(key[, c("arm","tf","n_targets","delta","ci_lo","ci_hi","excludes_egr2_sized","padj")],
        row.names = FALSE, digits = 2) }

# ---- figure -------------------------------------------------------------------
res$arm <- factor(res$arm, levels = c("induced", "repressed"),
                  labels = c("Induced arm", "Repressed arm"))
# The band is the RESOLUTION LIMIT, not a confidence region: the smallest deviation
# a typical regulon in that arm could have been called significant at, after BH
# across all regulons tested. Anything inside it is beyond this analysis to detect.
band <- do.call(rbind, lapply(levels(res$arm), function(a) {
  s <- res[res$arm == a, ]
  lim <- BH_SE_MULT * RATIO_SD / sqrt(median(s$n_targets))
  data.frame(arm = factor(a, levels = levels(res$arm)),
             lo = -lim, hi = lim, stringsAsFactors = FALSE)
}))
EGR2_EFFECT <- 0.055   # the adjusted effect from script 06, for scale
lab <- do.call(rbind, lapply(levels(res$arm), function(a) {
  s <- res[res$arm == a, ]
  data.frame(arm = factor(a, levels = levels(res$arm)),
             label = sprintf("%d regulons | %d significant after BH | resolves |delta| > %.2f for a typical regulon",
                             nrow(s), sum(s$padj < PADJ),
                             BH_SE_MULT * RATIO_SD / sqrt(median(s$n_targets))),
             stringsAsFactors = FALSE)
}))
p <- ggplot(res, aes(delta, n_targets)) +
  geom_rect(data = band, aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey88") +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_vline(xintercept = c(-EGR2_EFFECT, EGR2_EFFECT), linetype = "dashed",
             colour = "#B9600F", linewidth = 0.4) +
  geom_point(aes(colour = padj < PADJ), size = 1.9, alpha = 0.85) +
  ggrepel::geom_text_repel(data = res[res$tf %in% c("NFKB1", "RELA", "STAT6", "SPI1"), ],
                           aes(label = tf), size = 2.5, colour = "grey25",
                           min.segment.length = 0, seed = 1, max.overlaps = Inf) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.7, size = 2.6, colour = "grey30") +
  facet_wrap(~ arm, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("FALSE" = "#34495E", "TRUE" = "#C0392B"),
                      labels = c("FALSE" = "not significant", "TRUE" = "BH < 0.05"), name = NULL) +
  labs(x = "regulon mean amplitude ratio minus the program mean",
       y = "targets in the program",
       title = "No regulon shows a large selective deficit",
       subtitle = paste("Each point is a CollecTRI regulon with >=15 targets in that arm.",
                        "Shaded band = what this analysis cannot resolve for a typical regulon;",
                        "orange dashed lines = the EGR2 effect size (script 06), for scale.",
                        sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(out_dir, "tf_regulon_amplitude.pdf"), p, width = 7.0, height = 6.0)
cat("Saved: results/tf_regulon_amplitude.pdf\n")

wb <- createWorkbook()
addWorksheet(wb, "regulon amplitude"); writeData(wb, "regulon amplitude", res)
saveWorkbook(wb, file.path(out_dir, "tf_enrichment_results.xlsx"), overwrite = TRUE)
cat("Saved: results/tf_enrichment_results.xlsx\n")
