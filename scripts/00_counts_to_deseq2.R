# ===============================================================================
# 00 — FEATURECOUNTS -> DESeq2 OBJECTS
# ===============================================================================
# Builds the 4h and 24h DESeq2 objects from the GEO-deposited featureCounts table
# (GSE322520), restricted to the 100 ng/mL dose used in the manuscript.
#
# Upstream DESeq2/normalization logic by Iyad Sayed Issa ("Feature Counts to
# DDS_{4,24}hrs.R", repo miller-jensen-lab/il4-pulse-macrophage-rnaseq); adapted
# here to (a) read the public GEO counts as the single canonical input, (b) use
# the 100 ng/mL samples only (the 10 ng/mL dose arms were not part of the
# published analysis), and (c) use portable paths.
#
# Design:
#   4h : ~ batch + condition   (sequencing batch: replicate 1 = run 1,
#         replicates 2 & 3 = run 2). Conditions: ct, c100, p100.
#   24h: ~ condition           (single sequencing batch). Conditions: ct, c100, p100.
#
# Input : data/GSE322520_gene_counts_combined_4h_24h.txt   (featureCounts format)
# Output: data/dds_4hrs.rds, data/dds_24hrs.rds            (regenerable; gitignored)
#         data/vst_4hrs_batchcorrected.csv, data/vst_24hrs.csv  (for clustering/heatmaps)
# ===============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

data_dir <- file.path(getwd(), "data")

# ---- Counts input: GEO archive, or the per-timepoint files committed here --------
# The GEO submission (GSE322520) is a single 24 MB combined matrix and is gitignored.
# If it is absent, fall back to the two committed featureCounts_*_cleaned.txt files so
# the pipeline runs from a clean clone with no network access and no GEO credentials.
# The two sources carry identical row order and identical counts in the columns used
# here; the assertions below re-check that on every run rather than trusting it.
# The cleaned 4h file also carries the 10 ng/mL arms (c10_*/p10_*), which this
# analysis excludes -- only the 100 ng/mL columns are mapped across.
geo_file <- file.path(data_dir, "GSE322520_gene_counts_combined_4h_24h.txt")
fc_4h    <- file.path(data_dir, "featureCounts_4hours_cleaned.txt")
fc_24h   <- file.path(data_dir, "featureCounts_24hours_cleaned.txt")

if (file.exists(geo_file)) {
  # featureCounts layout: col 1 = Geneid, cols 2-6 = Chr/Start/End/Strand/Length,
  # remaining columns = per-sample raw counts.
  raw <- read.delim(geo_file, header = TRUE, check.names = FALSE)
  rownames(raw) <- raw$Geneid
  cat("[00] counts source: GEO combined matrix\n")
} else if (file.exists(fc_4h) && file.exists(fc_24h)) {
  a <- read.delim(fc_4h,  header = TRUE, check.names = FALSE)
  b <- read.delim(fc_24h, header = TRUE, check.names = FALSE)
  if (!identical(a$EnsemblID, b$EnsemblID))
    stop("featureCounts_4hours_cleaned.txt and featureCounts_24hours_cleaned.txt ",
         "do not share row order; cannot combine them safely")
  # map the committed column names onto the GEO names the rest of this script uses
  ren <- c(Ct_1 = "ct_4h_1", Ct_2 = "ct_4h_2", Ct_3 = "ct_4h_3",
           c_100_1 = "c100_4h_1", c_100_2 = "c100_4h_2", c_100_3 = "c100_4h_3",
           p_100_1 = "p100_4h_1", p_100_2 = "p100_4h_2", p_100_3 = "p100_4h_3",
           Untreated1 = "ct_24h_1", Untreated2 = "ct_24h_2", Untreated3 = "ct_24h_3",
           C1 = "c100_24h_1", C2 = "c100_24h_2", C3 = "c100_24h_3",
           P1 = "p100_24h_1", P2 = "p100_24h_2", P3 = "p100_24h_3")
  src <- cbind(a[, setdiff(colnames(a), "EnsemblID"), drop = FALSE],
               b[, setdiff(colnames(b), "EnsemblID"), drop = FALSE])
  missing <- setdiff(unname(ren), colnames(src))
  if (length(missing))
    stop("Committed count files lack expected columns: ", paste(missing, collapse = ", "))
  raw <- data.frame(Geneid = a$EnsemblID, src[, unname(ren), drop = FALSE],
                    check.names = FALSE, stringsAsFactors = FALSE)
  colnames(raw)[-1] <- names(ren)
  rownames(raw) <- raw$Geneid
  cat("[00] counts source: committed featureCounts_*_cleaned.txt (GEO file absent)\n")
} else {
  stop("No counts input found. Expected either\n  ", geo_file,
       "\nor both of\n  ", fc_4h, "\n  ", fc_24h)
}

# ---- Sample selection + renaming (100 ng/mL only) ----
# 4h columns and the 24h untreated columns are renamed; 24h C*/P* are kept as-is.
cols_4h  <- c(Ct_1 = "ct_1", Ct_2 = "ct_2", Ct_3 = "ct_3",
              c_100_1 = "c100_1", c_100_2 = "c100_2", c_100_3 = "c100_3",
              p_100_1 = "p100_1", p_100_2 = "p100_2", p_100_3 = "p100_3")
cols_24h <- c(Untreated1 = "untreated1", Untreated2 = "untreated2", Untreated3 = "untreated3",
              C1 = "C1", C2 = "C2", C3 = "C3", P1 = "P1", P2 = "P2", P3 = "P3")

condition_of <- function(x)
  ifelse(grepl("^c100|^C[0-9]", x), "c100",
  ifelse(grepl("^p100|^P[0-9]", x), "p100",
  ifelse(grepl("^ct|^untreated", x), "ct", NA_character_)))

# ---- Minimum-expression filter -------------------------------------------------
# Applied to the UNION across timepoints so both objects keep identical rownames
# (01 merges them by row), and BEFORE DESeq() so the dispersion trend is not fit
# through the tens of thousands of near-zero rows.
MIN_COUNT   <- 10   # a gene must reach this count...
MIN_SAMPLES <- 3    # ...in at least this many samples (= smallest group size)

count_matrix <- function(colmap) {
  missing <- setdiff(names(colmap), colnames(raw))
  if (length(missing)) stop("Missing expected sample columns: ", paste(missing, collapse = ", "))
  m <- as.matrix(raw[, names(colmap)])
  colnames(m) <- unname(colmap)
  mode(m) <- "integer"
  m
}
m4 <- count_matrix(cols_4h); m24 <- count_matrix(cols_24h)
stopifnot(identical(rownames(m4), rownames(m24)))
keep <- (rowSums(m4 >= MIN_COUNT) >= MIN_SAMPLES) | (rowSums(m24 >= MIN_COUNT) >= MIN_SAMPLES)
cat(sprintf("[00] Expression filter (>=%d counts in >=%d samples, either timepoint): %d of %d genes kept\n",
            MIN_COUNT, MIN_SAMPLES, sum(keep), length(keep)))
m4 <- m4[keep, , drop = FALSE]; m24 <- m24[keep, , drop = FALSE]

build_dds <- function(m, with_batch) {
  cond <- factor(condition_of(colnames(m)), levels = c("ct", "c100", "p100"))
  if (anyNA(cond)) stop("Unmapped condition for: ", paste(colnames(m)[is.na(cond)], collapse = ", "))
  cd <- data.frame(row.names = colnames(m), condition = cond)

  if (with_batch) {
    cd$batch <- factor(ifelse(grepl("_1$", colnames(m)), "run1", "run2"))
    dds <- DESeqDataSetFromMatrix(m, cd, design = ~ batch + condition)
  } else {
    dds <- DESeqDataSetFromMatrix(m, cd, design = ~ condition)
  }
  DESeq(dds)
}

cat("[00] Building 4h DESeq2 object (~ batch + condition, 100 ng/mL)...\n")
dds_4hrs  <- build_dds(m4,  with_batch = TRUE)
cat("[00] Building 24h DESeq2 object (~ condition, 100 ng/mL)...\n")
dds_24hrs <- build_dds(m24, with_batch = FALSE)

saveRDS(dds_4hrs,  file.path(data_dir, "dds_4hrs.rds"))
saveRDS(dds_24hrs, file.path(data_dir, "dds_24hrs.rds"))

# ---- VST matrices for downstream clustering / heatmaps (4h batch-corrected) ----
vst4    <- assay(vst(dds_4hrs, blind = FALSE))
vst4_bc <- limma::removeBatchEffect(vst4, batch = colData(dds_4hrs)$batch)
write.csv(vst4_bc, file.path(data_dir, "vst_4hrs_batchcorrected.csv"), quote = FALSE)
write.csv(assay(vst(dds_24hrs, blind = FALSE)), file.path(data_dir, "vst_24hrs.csv"), quote = FALSE)

cat(sprintf("[00] Done. 4h: %d samples, 24h: %d samples, %d genes.\n",
            ncol(dds_4hrs), ncol(dds_24hrs), nrow(dds_4hrs)))
