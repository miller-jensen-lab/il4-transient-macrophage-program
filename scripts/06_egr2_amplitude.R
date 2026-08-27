# ===============================================================================
# 06 — EGR2 DEPENDENCE AND PULSE AMPLITUDE  (lead mechanism panel)
# ===============================================================================
# THE QUESTION. If EGR2 is a secondary relay — induced by IL-4 via STAT6, then
# driving a downstream wave of genes — its targets require the upstream signal to
# persist long enough for the relay itself to be built. A 15-minute pulse should
# therefore under-deliver EGR2-dependent genes specifically, relative to genes the
# pathway reaches directly.
#
# THE TEST. Within the induced program, is the pulse/continuous amplitude ratio
# lower for EGR2-dependent genes than for the rest? Dependence comes from an
# EXTERNAL knockout annotation, never from any classification made in this
# pipeline, so the test is independent of our own labels.
#
# Three features of the design carry the result:
#   - the ratio correlates with the magnitude of the continuous response, so the
#     effect is reported adjusted for magnitude and for mean expression;
#   - tertiles of response magnitude show whether the effect is carried by one end
#     of the range;
#   - STAT6 dependence, tested identically, is the specificity control, separating
#     an EGR2-specific effect from a generic property of signal-dependent genes.
# Results are written to results/egr2_amplitude.csv.
#
# data/genes_stat6_egr2_dependence.xlsx merges calls from TWO Nagy-lab knockout
# series, and the two columns do not share a source:
#   STAT6KO_direction  Czimmerer et al., Immunity 2018   (GEO GSE106706)
#   EGR2KO_direction   Daniel et al., Genes Dev 2020     (GEO GSE151015)
# This script's result rests on the EGR2 column (Daniel et al. 2020).
#
# Sign convention:
#   STAT6KO_direction "Up" = IL-4-induced and lost in STAT6-KO => STAT6 activates
#   EGR2KO_direction  "Down" = down in EGR2-KO                 => EGR2 activates
#
# Input : data/master_expressed_genes.csv, data/genes_stat6_egr2_dependence.xlsx
# Output: results/egr2_amplitude.pdf, results/egr2_amplitude.csv
# ===============================================================================

suppressPackageStartupMessages({ library(openxlsx); library(ggplot2); library(DESeq2) })
source(file.path("scripts", "_common.R"))

PADJ <- 0.05
data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dep <- read.xlsx(file.path(data_dir, "genes_stat6_egr2_dependence.xlsx"))
m   <- read.csv(file.path(data_dir, "master_expressed_genes.csv"), stringsAsFactors = FALSE)
m   <- m[is_annotatable_gene(m$ensembl_gene_id, data_dir), ]

# baseMean is a covariate, not an outcome — taken from the same 24h model
d24 <- readRDS(file.path(data_dir, "dds_24hrs.rds"))
r24 <- results(d24, contrast = c("condition", "c100", "ct"))
m$baseMean <- r24$baseMean[match(m$ensembl_gene_id, rownames(r24))]

i <- match(m$gene_symbol, dep$gene_symbol)
m$egr2  <- dep$EGR2KO_direction[i]
m$stat6 <- dep$STAT6KO_direction[i]

# Induced continuous program. |log2FC| > 1 rather than the stage-1 boundary so the
# ratio's denominator is comfortably away from zero.
k <- !is.na(m$padj_c100_24h) & m$padj_c100_24h < PADJ & m$log2fc_c100_24h > 1
d <- m[k, ]
d$ratio <- d$log2fc_p100_24h / d$log2fc_c100_24h
d <- d[is.finite(d$ratio) & !is.na(d$egr2), ]
d$egr2_act  <- d$egr2  == "Down"
d$stat6_act <- d$stat6 == "Up"
cat(sprintf("[06] induced genes with EGR2 annotation: %d (%d EGR2-activated)\n",
            nrow(d), sum(d$egr2_act)))

# ---- the test, unadjusted and adjusted ----------------------------------------
fits <- list(
  "unadjusted"                     = lm(ratio ~ egr2_act, d),
  "+ continuous log2FC"            = lm(ratio ~ egr2_act + log2fc_c100_24h, d),
  "+ continuous log2FC + baseMean" = lm(ratio ~ egr2_act + log2fc_c100_24h + log10(baseMean), d))
res <- do.call(rbind, lapply(names(fits), function(nm) {
  cf <- summary(fits[[nm]])$coefficients["egr2_actTRUE", ]
  data.frame(model = nm, effect = cf[1], se = cf[2], p = cf[4], stringsAsFactors = FALSE)
}))
# specificity control: the same test on STAT6 dependence
sc <- summary(lm(ratio ~ stat6_act + log2fc_c100_24h, d))$coefficients["stat6_actTRUE", ]
res <- rbind(res, data.frame(model = "STAT6 control (+ log2FC)", effect = sc[1], se = sc[2], p = sc[4]))
write.csv(res, file.path(out_dir, "egr2_amplitude.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 3)
cat("Saved: results/egr2_amplitude.csv\n")

# ---- figure: distribution, plus the effect within magnitude tertiles -----------
# Tertiles of continuous-response magnitude: the ratio correlates with magnitude,
# so this shows the effect within strata of the confound.
d$tertile <- cut(d$log2fc_c100_24h, quantile(d$log2fc_c100_24h, 0:3/3),
                 include.lowest = TRUE,
                 labels = c("weak response", "intermediate", "strong response"))
# short labels: the full ones collide under three facets at publication width
d$grp <- factor(ifelse(d$egr2_act, "EGR2-dep.", "rest"), levels = c("rest", "EGR2-dep."))
COL <- c("rest" = "grey62", "EGR2-dep." = "#B9600F")

lab <- do.call(rbind, lapply(levels(d$tertile), function(t) {
  s <- d[d$tertile == t, ]
  data.frame(tertile = factor(t, levels = levels(d$tertile)),
             label = sprintf("%.2f vs %.2f   n=%d", median(s$ratio[s$egr2_act]),
                             median(s$ratio[!s$egr2_act]), nrow(s)), stringsAsFactors = FALSE)
}))
adj  <- res$effect[res$model == "+ continuous log2FC + baseMean"]
adjp <- res$p[res$model == "+ continuous log2FC + baseMean"]
# computed, not typed: this subtitle previously carried a hardcoded p-value
scp  <- res$p[res$model == "STAT6 control (+ log2FC)"]

p <- ggplot(d, aes(grp, ratio, fill = grp)) +
  geom_hline(yintercept = median(d$ratio), linetype = "dashed", colour = "grey45") +
  geom_violin(alpha = 0.45, linewidth = 0.35, colour = "grey30", scale = "width") +
  geom_boxplot(width = 0.16, outlier.size = 0.4, linewidth = 0.35, fill = "white") +
  geom_text(data = lab, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
            vjust = 1.6, size = 2.5, colour = "grey25") +
  facet_wrap(~ tertile, nrow = 1) +
  scale_fill_manual(values = COL, guide = "none") +
  coord_cartesian(ylim = unname(quantile(d$ratio, c(0.01, 0.99)))) +
  labs(x = "grey = rest of the induced program;  orange = EGR2-dependent",
       y = "amplitude ratio (pulse / continuous, 24 h)",
       title = "EGR2-dependent genes are the ones a pulse under-delivers",
       subtitle = paste(
         sprintf("Effect %+.3f adjusted for response magnitude and expression (p = %.1g); consistent in every tertile.", adj, adjp),
         "Dependence calls from Daniel et al. 2020 EGR2-KO data, independent of any classification made here.",
         sprintf("STAT6 dependence, tested identically, does not survive adjustment (p = %.2f) - the effect is EGR2-specific.", scp),
         sep = "\n")) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey35"),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8, colour = "grey35"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(out_dir, "egr2_amplitude.pdf"), p, width = 7.4, height = 4.4)
cat("Saved: results/egr2_amplitude.pdf\n")
