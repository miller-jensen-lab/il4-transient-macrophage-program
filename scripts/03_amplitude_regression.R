# ===============================================================================
# 03 — PULSE-VS-CONTINUOUS AMPLITUDE REGRESSION (threshold-free primary claim)
# ===============================================================================
# Regresses pulse log2FC on continuous log2FC over the stage-1 program, separately
# for INDUCED and REPRESSED genes at each timepoint. The slope is the fraction of
# the continuous amplitude that the pulse reproduces.
#
# A threshold-free alternative to per-gene classification: one slope with a
# confidence interval per arm per timepoint. Unlike a ratio-based measure it has no
# denominator that is systematically smaller for repressed genes, which would
# inflate the apparent difference between arms. The recovery between timepoints is
# what the interpretation rests on: induction climbs toward continuous, repression
# does not.
#
# Authoritative slopes are in results/amplitude_regression.csv, rewritten every run.
#
# SENSITIVITY: the same slopes are recomputed on apeglm-shrunken estimates. apeglm
# fits an independent prior per coefficient, so a weaker coefficient is shrunk
# harder -- at 4 h the pulse (y) axis compresses far more than the continuous (x)
# axis, while at 24 h both are barely touched. Read the 4 h shrunken slope as a
# lower bound, not as independent evidence that the 4 h response is weak; that
# would be circular. NB apeglm, NOT type="normal" — the latter over-shrinks at n=3.
#
# Input : data/master_expressed_genes.csv, data/dds_{4,24}hrs.rds
# Output: results/amplitude_regression.pdf
#         results/amplitude_regression.csv   (slopes, CIs, n, both estimators)
# ===============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(DESeq2); library(ggrepel) })
source(file.path("scripts", "_common.R"))   # canonical gene universe

# Arm membership uses the SAME criterion as the stage-1 program in 01: the padj in
# master_expressed_genes.csv already tests H0: |log2FC| <= log2(1.5), so a gene is
# selected on padj alone and assigned to an arm by the SIGN of its continuous
# response. Adding a |log2FC| cut here would re-introduce the point-estimate
# filtering that 01 removed, and would fit the headline on a different population
# from the program it describes.
PADJ <- 0.05

data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

m <- read.csv(file.path(data_dir, "master_expressed_genes.csv"), stringsAsFactors = FALSE)

# PROTEIN-CODING ONLY, as everywhere else: the slopes must be fitted on the same
# gene universe as the figures they support. The repressed arm is sensitive to this
# -- a single high-leverage pseudogene with a large continuous and small pulse
# response can move the 4 h slope by several points.
m <- m[is_annotatable_gene(m$ensembl_gene_id, data_dir), ]
cat(sprintf("[03] protein-coding expressed genes: %d\n", nrow(m)))

# ---- shrunken estimates for the sensitivity arm -------------------------------
have_apeglm <- requireNamespace("apeglm", quietly = TRUE)
if (!have_apeglm) warning("apeglm not installed - shrunken sensitivity slopes skipped", call. = FALSE)
shrunk <- NULL
if (have_apeglm) {
  d4 <- readRDS(file.path(data_dir, "dds_4hrs.rds")); d24 <- readRDS(file.path(data_dir, "dds_24hrs.rds"))
  s <- function(dds, cf) {
    r <- as.data.frame(lfcShrink(dds, coef = cf, type = "apeglm", quiet = TRUE))
    setNames(r$log2FoldChange, rownames(r))
  }
  shrunk <- list(c100_4h  = s(d4,  "condition_c100_vs_ct"), p100_4h  = s(d4,  "condition_p100_vs_ct"),
                 c100_24h = s(d24, "condition_c100_vs_ct"), p100_24h = s(d24, "condition_p100_vs_ct"))
}

fit_one <- function(tp, direction, estimator) {
  cc_mle <- m[[paste0("log2fc_c100_", tp)]]
  qq     <- m[[paste0("padj_c100_",   tp)]]
  # The ARM is always defined on the MLE + its adjusted p-value, so both estimators
  # are fitted on exactly the same genes and the comparison is like-for-like.
  keep <- !is.na(qq) & qq < PADJ & !is.na(cc_mle) &
          (if (direction == "induced") cc_mle > 0 else cc_mle < 0)
  if (estimator == "MLE") {
    x <- cc_mle; y <- m[[paste0("log2fc_p100_", tp)]]
  } else {
    x <- unname(shrunk[[paste0("c100_", tp)]][m$ensembl_gene_id])
    y <- unname(shrunk[[paste0("p100_", tp)]][m$ensembl_gene_id])
  }
  keep <- keep & !is.na(x) & !is.na(y)
  f  <- lm(y[keep] ~ x[keep]); ci <- confint(f)[2, ]
  data.frame(timepoint = tp, arm = direction, estimator = estimator, n = sum(keep),
             slope = unname(coef(f)[2]), lo = ci[1], hi = ci[2],
             r2 = summary(f)$r.squared, stringsAsFactors = FALSE)
}

ests <- if (have_apeglm) c("MLE", "shrunk") else "MLE"
res <- do.call(rbind, lapply(ests, function(e)
         do.call(rbind, lapply(c("4h", "24h"), function(tp)
           do.call(rbind, lapply(c("induced", "repressed"), fit_one, tp = tp, estimator = e))))))
write.csv(res, file.path(out_dir, "amplitude_regression.csv"), row.names = FALSE)
cat("Saved: results/amplitude_regression.csv\n")
print(res, row.names = FALSE, digits = 3)

# ---- formal test that the two arms differ -------------------------------------
# Interaction of continuous log2FC with arm, on the MLE scale, repression sign-
# flipped so both arms share an "amplitude reproduced" scale.
#
# THIS TEST REQUIRES TWO CORRECTIONS. Plain OLS on y ~ x*arm assumes one common
# error variance and a common linear model over a common x range. Neither holds
# here: the residual SD differs several-fold between arms, the induced arm is
# strongly NONLINEAR, and the two arms span non-overlapping x ranges because the
# induced arm reaches much further. Uncorrected OLS therefore compares slopes
# estimated over different, curvature-relevant stretches of a curved relationship
# using a pooled sigma, and can report an arm difference that is an artifact of the
# variance model rather than a property of the data.
#
# So: report the HC3 heteroscedasticity-consistent p-value, with the matched-range
# refit as the sensitivity analysis. Both are written to
# results/amplitude_arm_test.csv alongside the naive OLS value, so the reader can
# see what each correction does. Do not quote the OLS p-value on its own.
tst <- function(tp) {
  cc <- m[[paste0("log2fc_c100_", tp)]]; pp <- m[[paste0("log2fc_p100_", tp)]]
  qq <- m[[paste0("padj_c100_",   tp)]]
  up <- !is.na(qq) & qq < PADJ & !is.na(cc) & cc > 0
  dn <- !is.na(qq) & qq < PADJ & !is.na(cc) & cc < 0
  d <- rbind(data.frame(x = cc[up], y = pp[up], arm = "induced"),
             data.frame(x = -cc[dn], y = -pp[dn], arm = "repressed"))
  d <- d[!is.na(d$x) & !is.na(d$y), ]; d$arm <- factor(d$arm, levels = c("induced", "repressed"))
  f     <- lm(y ~ x * arm, data = d)
  p_ols <- summary(f)$coefficients["x:armrepressed", 4]
  p_hc3 <- lmtest::coeftest(f, vcov. = sandwich::vcovHC(f, type = "HC3"))["x:armrepressed", 4]
  xm    <- max(d$x[d$arm == "repressed"])          # match the amplitude range
  fm    <- lm(y ~ x * arm, data = d[d$x <= xm, ])
  p_mat <- lmtest::coeftest(fm, vcov. = sandwich::vcovHC(fm, type = "HC3"))["x:armrepressed", 4]
  data.frame(timepoint = tp,
             sd_induced = sd(resid(f)[d$arm == "induced"]),
             sd_repressed = sd(resid(f)[d$arm == "repressed"]),
             p_ols = p_ols, p_hc3 = p_hc3,
             matched_xmax = xm,
             slope_induced_matched   = unname(coef(fm)["x"]),
             slope_repressed_matched = unname(coef(fm)["x"] + coef(fm)["x:armrepressed"]),
             p_hc3_matched = p_mat, stringsAsFactors = FALSE)
}
if (!requireNamespace("sandwich", quietly = TRUE) || !requireNamespace("lmtest", quietly = TRUE))
  stop("packages 'sandwich' and 'lmtest' are required for the heteroscedasticity-robust ",
       "interaction test; install them before running this script")
arm_test <- do.call(rbind, lapply(c("4h", "24h"), tst))
write.csv(arm_test, file.path(out_dir, "amplitude_arm_test.csv"), row.names = FALSE)
cat("\nArm-difference test (HC3-robust; matched range as sensitivity):\n")
print(arm_test, row.names = FALSE, digits = 3)
cat("Saved: results/amplitude_arm_test.csv\n")

# ---- figure -------------------------------------------------------------------
pd <- do.call(rbind, lapply(c("4h", "24h"), function(tp) {
  cc <- m[[paste0("log2fc_c100_", tp)]]; pp <- m[[paste0("log2fc_p100_", tp)]]
  qq <- m[[paste0("padj_c100_", tp)]]
  k  <- !is.na(qq) & qq < PADJ & !is.na(cc) & !is.na(pp)
  data.frame(timepoint = tp, cont = cc[k], pulse = pp[k], gene = m$gene_symbol[k],
             arm = ifelse(cc[k] > 0, "induced", "repressed"))
}))
pd$timepoint <- factor(pd$timepoint, levels = c("4h", "24h"))
lab <- res[res$estimator == "MLE", ]
lab$timepoint <- factor(lab$timepoint, levels = c("4h", "24h"))
lab$x <- ifelse(lab$arm == "induced", 10, -1); lab$y <- ifelse(lab$arm == "induced", -3, 12)

# The identity line is the whole visual reference, so the axes MUST be equally
# scaled. Unequal ranges tilt the 45-degree reference and make every slope look
# shallower than it is.
AX <- c(floor(min(pd$cont, pd$pulse)), ceiling(max(pd$cont, pd$pulse)))
# Same marker genes as the volcano in 02, so a reader can follow Arg1 / Retnla /
# Ifit3b from the volcano into this panel and see them move toward identity.
mk <- pd[pd$gene %in% MARKER_GENES, ]
for (tp in levels(pd$timepoint)) {
  miss <- setdiff(MARKER_GENES, mk$gene[mk$timepoint == tp])
  if (length(miss))
    cat(sprintf("[03] %s: marker not in the plotted set (continuous not significant): %s\n",
                tp, paste(miss, collapse = ", ")))
}

p <- ggplot(pd, aes(cont, pulse, colour = arm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey55") +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_point(data = mk, aes(cont, pulse), inherit.aes = FALSE,
             size = 1.3, shape = 21, fill = "white", colour = "black", stroke = 0.45) +
  ggrepel::geom_text_repel(data = mk, aes(cont, pulse, label = gene), inherit.aes = FALSE,
                           size = 2.4, colour = "black", segment.colour = "grey45",
                           segment.size = 0.25, min.segment.length = 0,
                           box.padding = 0.4, point.padding = 0.2,
                           max.overlaps = Inf, seed = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, formula = y ~ x) +
  geom_text(data = lab, aes(x = x, y = y, colour = arm,
                            label = sprintf("slope %.2f [%.2f-%.2f]", slope, lo, hi)),
            inherit.aes = FALSE, size = 2.7, hjust = 0.5, show.legend = FALSE) +
  facet_wrap(~ timepoint) +
  coord_fixed(ratio = 1, xlim = AX, ylim = AX) +
  scale_colour_manual(values = c(induced = "#C0392B", repressed = "#2980B9"), name = NULL) +
  labs(title = "Amplitude of the pulse response relative to continuous IL-4",
       subtitle = paste("Each point is a gene significantly regulated by continuous IL-4 (formal test, >1.5-fold, padj < 0.05).",
                        "Slope = fraction of the continuous amplitude reproduced; dotted line is identity.",
                        sep = "\n"),
       x = expression(continuous~log[2]~FC), y = expression(pulse~log[2]~FC)) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(out_dir, "amplitude_regression.pdf"), p, width = 7.2, height = 4.4)
cat("Saved: results/amplitude_regression.pdf\n")
