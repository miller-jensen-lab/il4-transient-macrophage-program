# ===============================================================================
# 02 — VOLCANO PLOTS + PER-TIMEPOINT VENN DIAGRAMS (stage-1 program description)
# ===============================================================================
# Descriptive panels for the stage-1 DEG program:
#   (1) Volcano plots for pulse and continuous IL-4 at 4h and 24h (2x2).
#   (2) Venn diagrams of continuous vs pulse DEGs, separated by direction
#       (Up / Down) WITHIN each timepoint (4h and 24h).
#
# Universe and threshold match 01 and 03; see 01 for the definition.
# NB: these Venns are PER TIMEPOINT (continuous-4h vs pulse-4h, etc.).
#
# Input : data/master_expressed_genes.csv   (all expressed genes; log2FC + padj, from 01)
# Output: results/volcano_4panel.pdf
#         results/venn_by_timepoint.pdf       (per-timepoint 2-way Venns, Up/Down)
#         results/venn_by_timepoint.csv
#         results/venn_4way.pdf               (4-condition Venn, Up & Down; limma)
#         results/venn_4way.csv               (all 4-way intersection counts)
# ===============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(grid); library(VennDiagram); library(ggrepel) })
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

# The padj in master_expressed_genes.csv already tests H0: |log2FC| <= log2(1.5),
# so significance alone defines a DEG. LFC is retained only to draw the reference
# lines at the tested fold-change -- it is NOT an additional filter.
LFC  <- log2(1.5)
PADJ <- 0.05
YCAP <- 50       # clip -log10(padj) at this value for volcano display
YTOP <- YCAP * 1.16   # headroom above the ceiling for the up/down count annotations

# Volcano x-axis, fixed and asymmetric because the induced arm reaches much further
# than the repressed arm; XMAX = 18 covers the observed range. `oob = squish` is a
# guard: if the data ever exceed this, points pile up at the edge instead of
# silently vanishing, and the console prints a count.
XMIN <- -7
XMAX <- 18

# Marker genes come from scripts/_common.R (MARKER_GENES) so the volcano and the
# amplitude scatter in 03 label the identical set.

source(file.path("scripts", "_common.R"))   # canonical gene universe

data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

d <- read.csv(file.path(data_dir, "master_expressed_genes.csv"), stringsAsFactors = FALSE)

# protein-coding restriction, identical to scripts/01
d <- d[is_annotatable_gene(d$ensembl_gene_id, data_dir), ]
cat(sprintf("[02] expressed genes after protein-coding filter: %d\n", nrow(d)))

# the four contrasts, as (column suffix, condition label, timepoint)
contrasts <- list(
  c100_4h  = list(cond = "Continuous", tp = "4h"),  p100_4h  = list(cond = "Pulse", tp = "4h"),
  c100_24h = list(cond = "Continuous", tp = "24h"), p100_24h = list(cond = "Pulse", tp = "24h"))

status <- function(lfc, padj) {
  sig <- !is.na(padj) & padj < PADJ & !is.na(lfc)
  ifelse(sig & lfc > 0, "Up",
  ifelse(sig & lfc < 0, "Down", "NS"))
}
# Display class: significant points are BLACK regardless of direction. Direction is
# still carried by `status` because the per-facet up/down counts depend on it.
sig_class <- function(st) ifelse(st == "NS", "NS", "Significant")

# ===============================================================================
# (1) Volcano plots — 2x2 (timepoint rows x condition cols)
# ===============================================================================
vlist <- lapply(names(contrasts), function(lab) {
  lfc <- d[[paste0("log2fc_", lab)]]; padj <- d[[paste0("padj_", lab)]]
  keep <- !is.na(lfc) & !is.na(padj)
  y <- -log10(padj[keep]); y[is.infinite(y) | y > YCAP] <- YCAP
  data.frame(gene = d$gene_symbol[keep], log2fc = lfc[keep], neglogp = y,
             status = status(lfc[keep], padj[keep]),
             condition = contrasts[[lab]]$cond, timepoint = contrasts[[lab]]$tp,
             stringsAsFactors = FALSE)
})
v <- do.call(rbind, vlist)
v$condition <- factor(v$condition, levels = c("Continuous", "Pulse"))
v$timepoint <- factor(v$timepoint, levels = c("4h", "24h"))
v$status    <- factor(v$status, levels = c("Up", "Down", "NS"))
v$sig       <- factor(sig_class(as.character(v$status)), levels = c("Significant", "NS"))

# per-facet up/down counts (annotation)
cnt <- aggregate(gene ~ condition + timepoint + status, v[v$status != "NS", ], length)
lab_df <- do.call(rbind, lapply(split(cnt, list(cnt$condition, cnt$timepoint), drop = TRUE), function(s) {
  data.frame(condition = s$condition[1], timepoint = s$timepoint[1],
             up   = sum(s$gene[s$status == "Up"]), down = sum(s$gene[s$status == "Down"]))
}))

# marker-gene coordinates, one row per gene per facet
mk <- v[v$gene %in% MARKER_GENES, ]
mk$x <- pmin(pmax(mk$log2fc, XMIN), XMAX)   # keep labels inside the panel
oob <- sum(v$log2fc < XMIN | v$log2fc > XMAX, na.rm = TRUE)
if (oob > 0) cat(sprintf("[02] volcano: %d point(s) outside [%g, %g] squished to the axis edge\n",
                         oob, XMIN, XMAX))

pv <- ggplot(v, aes(log2fc, neglogp, color = sig)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(data = mk, aes(x = x, y = neglogp), inherit.aes = FALSE,
             size = 1.1, shape = 21, fill = "white", colour = "black", stroke = 0.45) +
  ggrepel::geom_text_repel(data = mk, aes(x = x, y = neglogp, label = gene), inherit.aes = FALSE,
                           size = 2.5, colour = "black", segment.colour = "grey45",
                           segment.size = 0.25, min.segment.length = 0,
                           box.padding = 0.35, point.padding = 0.2,
                           max.overlaps = Inf, seed = 1,
                           ylim = c(NA, YCAP)) +
  # Counts sit in a headroom strip ABOVE the clipped points. Genes whose padj is at
  # the YCAP ceiling would otherwise collide with them, which is most of the 24h
  # markers -- hence YTOP rather than placing the counts at YCAP.
  geom_text(data = lab_df, aes(x = XMAX, y = YTOP, label = paste0("up ", up)),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3, color = "black") +
  geom_text(data = lab_df, aes(x = XMIN, y = YTOP, label = paste0("down ", down)),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3, color = "black") +
  facet_grid(timepoint ~ condition) +
  scale_color_manual(values = c(Significant = "black", NS = "grey80"), name = NULL) +
  scale_x_continuous(limits = c(XMIN, XMAX), oob = scales::squish,
                     breaks = seq(-6, 18, by = 4)) +
  scale_y_continuous(limits = c(0, YTOP), breaks = seq(0, YCAP, by = 10)) +
  labs(title = "IL-4 differential expression (vs untreated)",
       subtitle = sprintf("DEG: >1.5-fold & padj < %.2f;  -log10(padj) clipped at %g", PADJ, YCAP),
       x = expression(log[2]~fold~change), y = expression(-log[10](p[adj]))) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(out_dir, "volcano_4panel.pdf"), pv, width = 7.5, height = 7)
cat("Saved: results/volcano_4panel.pdf\n")

# ===============================================================================
# (2) Per-timepoint Venn diagrams: continuous vs pulse, Up & Down, at 4h and 24h
# ===============================================================================
deg <- function(tp, cond, dir) {
  lab  <- paste0(if (cond == "Continuous") "c100_" else "p100_", tp)
  lfc  <- d[[paste0("log2fc_", lab)]]; padj <- d[[paste0("padj_", lab)]]
  sig  <- !is.na(padj) & padj < PADJ & !is.na(lfc)
  if (dir == "Up") sig & lfc > 0 else sig & lfc < 0
}

regions <- do.call(rbind, lapply(c("4h", "24h"), function(tp) do.call(rbind, lapply(c("Up", "Down"), function(dir) {
  C <- deg(tp, "Continuous", dir); P <- deg(tp, "Pulse", dir)
  data.frame(timepoint = tp, direction = dir,
             n_continuous = sum(C), n_pulse = sum(P), n_both = sum(C & P),
             continuous_only = sum(C & !P), pulse_only = sum(P & !C), stringsAsFactors = FALSE)
}))))
write.csv(regions, file.path(out_dir, "venn_by_timepoint.csv"), row.names = FALSE)
cat("Saved: results/venn_by_timepoint.csv\n"); print(regions, row.names = FALSE)

venn_panel <- function(tp, dir, fill) {
  r <- regions[regions$timepoint == tp & regions$direction == dir, ]
  # Continuous first here, deliberately NOT matching the 4-way Venn's pulse-first
  # order: pulse DEGs are very nearly a subset of continuous, so these render as
  # nested circles with no meaningful left-right axis. Swapping the arguments does
  # not move pulse to the left, it only pushes the category labels off their
  # tuned cat.pos angles.
  draw.pairwise.venn(area1 = r$n_continuous, area2 = r$n_pulse, cross.area = r$n_both,
    category = c("Cont.", "Pulse"), fill = fill, alpha = 0.55, lwd = 1,
    cex = 1.1, fontfamily = "sans", cat.cex = 1, cat.fontfamily = "sans",
    cat.pos = c(-25, 25), margin = 0.08, ind = FALSE)
}

pdf(file.path(out_dir, "venn_by_timepoint.pdf"), width = 8, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 3,
  heights = unit(c(2, 1, 1), c("lines", "null", "null")),
  widths  = unit(c(2.4, 1, 1), c("lines", "null", "null")))))
grid.text("Continuous vs pulse DEGs by timepoint  (>1.5-fold, padj < 0.05)",
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3), gp = gpar(fontsize = 12, fontface = "bold"))
grid.text("Up (activated)",   vp = viewport(layout.pos.row = 1, layout.pos.col = 2), y = 0.1, gp = gpar(col = "#C0392B", fontface = "bold"))
grid.text("Down (repressed)", vp = viewport(layout.pos.row = 1, layout.pos.col = 3), y = 0.1, gp = gpar(col = "#2980B9", fontface = "bold"))
cells <- list(list(2, "4h"), list(3, "24h"))
for (cc in cells) {
  grid.text(cc[[2]], vp = viewport(layout.pos.row = cc[[1]], layout.pos.col = 1), gp = gpar(fontface = "bold"), rot = 90)
  pushViewport(viewport(layout.pos.row = cc[[1]], layout.pos.col = 2)); grid.draw(venn_panel(cc[[2]], "Up",   c("#E8A39B", "#C0392B"))); popViewport()
  pushViewport(viewport(layout.pos.row = cc[[1]], layout.pos.col = 3)); grid.draw(venn_panel(cc[[2]], "Down", c("#A3C4DD", "#2980B9"))); popViewport()
}
invisible(dev.off())
cat("Saved: results/venn_by_timepoint.pdf\n")

# ===============================================================================
# (3) 4-condition Venn (one diagram, all combinatorial overlaps) via limma,
#     separately for Up and Down. Pooling directions would miscount genes that
#     are induced in one condition and repressed in another, so we split.
#     Each diagram is restricted to genes that are a DEG in >= 1 of the four
#     conditions (so the "outside" region is 0 and not a distracting huge count).
# ===============================================================================
conds <- c("Continuous", "Pulse")
mat <- function(dir) {
  # Set order is PULSE first, then CONTINUOUS, each 4h before 24h. This puts the
  # two pulse ellipses adjacent and the two continuous ellipses adjacent, so the
  # within-condition timepoint overlap reads directly off the diagram.
  m <- cbind(
    "Pulse 4h"  = as.integer(deg("4h",  "Pulse",      dir)),
    "Pulse 24h" = as.integer(deg("24h", "Pulse",      dir)),
    "Cont 4h"   = as.integer(deg("4h",  "Continuous", dir)),
    "Cont 24h"  = as.integer(deg("24h", "Continuous", dir)))
  m[rowSums(m) > 0, , drop = FALSE]   # keep genes that are a DEG in >=1 condition
}
m_up <- mat("Up"); m_dn <- mat("Down")

# counts table (all 16 patterns) for both directions -> CSV
# unclass() drops the "VennCounts" class so as.data.frame keeps the 16x5 layout
vc_up <- as.data.frame(unclass(limma::vennCounts(m_up))); vc_up$direction <- "Up"
vc_dn <- as.data.frame(unclass(limma::vennCounts(m_dn))); vc_dn$direction <- "Down"
write.csv(rbind(vc_up, vc_dn), file.path(out_dir, "venn_4way.csv"), row.names = FALSE)
cat("Saved: results/venn_4way.csv\n")

pdf(file.path(out_dir, "venn_4way.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
# Hue encodes condition (pulse orange, continuous blue), lightness encodes
# timepoint (light 4h, dark 24h), matching the column order above.
VENN4_COL <- c("#E8A33D", "#B9600F", "#7FB3D5", "#1F618D")
limma::vennDiagram(m_up, circle.col = VENN4_COL,
                   names = colnames(m_up), cex = c(1, 0.8, 0.6), main = "Up (activated)")
limma::vennDiagram(m_dn, circle.col = VENN4_COL,
                   names = colnames(m_dn), cex = c(1, 0.8, 0.6), main = "Down (repressed)")
mtext("4-condition DEG overlap  (>1.5-fold, padj < 0.05; DEG in >=1 condition)",
      outer = TRUE, cex = 1.1, font = 2)
invisible(dev.off())
cat("Saved: results/venn_4way.pdf\n")
