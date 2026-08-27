# ===============================================================================
# 07 — CURATED TRANSCRIPTION-FACTOR HEATMAP
# ===============================================================================
# log2 fold-change (vs untreated) for a curated set of IL-4-regulated transcription
# factors, across all four contrasts. Ten induced and ten repressed, chosen a priori
# from the alternative-activation literature rather than selected from these data --
# so this panel is descriptive and carries no multiple-testing claim of its own.
# Cells whose adjusted p-value does not clear PADJ are labelled n.s.
#
# The gene list and its two blocks are fixed here. Do not reorder them to match a
# result; the point of an a priori list is that it is independent of the outcome.
#
# Input : data/master_expressed_genes.csv  (from 01)
# Output: results/tf_heatmap.pdf, results/tf_heatmap.csv
# ===============================================================================

suppressPackageStartupMessages({ library(ggplot2) })
source(file.path("scripts", "_common.R"))

PADJ <- 0.05
LIM  <- 4      # colour saturates here; matches the published scale

TF_UP <- c("Batf3", "Atf3", "Klf4", "Hivep3", "Fosl2",
           "Mitf", "Bhlhe40", "Pparg", "Egr2", "Irf4")
TF_DN <- c("Tox2", "Nfatc2", "Relb", "Pou2f2", "Dbp",
           "Irf7", "Maf", "Smad3", "Spic", "Nfkbiz")

# Two-line labels keep the four columns legible at publication width without
# rotating them.
CONDS <- list(
  "Transient\n4h"   = c("log2fc_p100_4h",  "padj_p100_4h"),
  "Transient\n24h"  = c("log2fc_p100_24h", "padj_p100_24h"),
  "Continuous\n4h"  = c("log2fc_c100_4h",  "padj_c100_4h"),
  "Continuous\n24h" = c("log2fc_c100_24h", "padj_c100_24h"))

data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

mexp <- file.path(data_dir, "master_expressed_genes.csv")
if (!file.exists(mexp)) stop("Missing ", mexp, " — run scripts/01 first.")
m <- read.csv(mexp, stringsAsFactors = FALSE)
m <- m[is_annotatable_gene(m$ensembl_gene_id, data_dir), ]

# Fail loudly rather than silently dropping a row from a fixed 20-gene panel.
missing <- setdiff(c(TF_UP, TF_DN), m$gene_symbol)
if (length(missing)) stop("TFs absent from the expressed set: ", paste(missing, collapse = ", "))

rows <- list()
for (blk in c("Increased", "Decreased")) {
  for (gn in if (blk == "Increased") TF_UP else TF_DN) {
    r <- m[m$gene_symbol == gn, ][1, ]
    for (cd in names(CONDS)) rows[[length(rows) + 1]] <- data.frame(
      block = blk, gene = gn, cond = cd,
      lfc = r[[CONDS[[cd]][1]]], padj = r[[CONDS[[cd]][2]]],
      stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
d$block <- factor(d$block, levels = c("Increased", "Decreased"))
d$cond  <- factor(d$cond,  levels = names(CONDS))
# levels reversed so the first gene listed renders at the TOP of its block
d$gene  <- factor(d$gene, levels = rev(c(TF_UP, TF_DN)))
d$sig   <- !is.na(d$padj) & d$padj < PADJ

write.csv(d, file.path(out_dir, "tf_heatmap.csv"), row.names = FALSE)
cat("Saved: results/tf_heatmap.csv\n")
cat(sprintf("[07] %d TFs x %d contrasts; %d cells significant at padj < %g\n",
            length(TF_UP) + length(TF_DN), length(CONDS), sum(d$sig), PADJ))

p <- ggplot(d, aes(cond, gene, fill = pmax(pmin(lfc, LIM), -LIM))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(data = d[!d$sig, ], aes(label = "n.s."), size = 2.4, colour = "grey25") +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-LIM, LIM), breaks = seq(-LIM, LIM, 2),
                       name = expression(log[2]~FC)) +
  labs(x = NULL, y = NULL,
       title = "IL-4-regulated transcription factors",
       subtitle = sprintf("log2 fold-change vs untreated\nn.s. = padj >= %g; colour saturates at +/-%g", PADJ, LIM)) +
  theme_bw(base_size = 10) +
  theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"),
        strip.text.y.left = element_text(angle = 0, size = 8),
        strip.placement = "outside", panel.grid = element_blank(),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8, face = "italic"))
ggsave(file.path(out_dir, "tf_heatmap.pdf"), p, width = 5.4, height = 5.6)
cat("Saved: results/tf_heatmap.pdf\n")
