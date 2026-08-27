# ===============================================================================
# 01 — DEG PROGRAM  (stage 1 of 2: formal test that |log2FC| exceeds 1.5-fold)
# ===============================================================================
# Extracts the four DESeq2 contrasts (pulse / continuous vs control at 4h & 24h),
# calls DEGs by TESTING H0: |log2FC| <= log2(1.5) at either timepoint, and builds
# the UNION "program" of genes activated and/or repressed by continuous and/or
# pulse IL-4. Also extracts the direct pulse-vs-continuous contrast.
#
# The resulting "program" is the basis for the volcano plots and Venn diagrams (02),
# the amplitude regression (03), and the pathway analyses (04).
#
# Differential-expression integration logic by Iyad Sayed Issa ("DDS to Union
# Table"); adapted for the integrated 100 ng/mL pipeline and portable paths.
#
# Input : data/dds_4hrs.rds, data/dds_24hrs.rds   (from 00)
# Output: data/deg_program.csv               (union DEG program)
#         data/master_expressed_genes.csv         (all expressed genes; log2FC + padj)
# ===============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(biomaRt)
})

# ---- Parameters ----
# Stage-1 DEGs are called by a FORMAL TEST of the fold-change claim, not by
# filtering a point estimate. `results(lfcThreshold = x)` tests H0: |log2FC| <= x
# (TREAT; McCarthy & Smyth 2009), so the reported padj and FDR refer to the claim
# the manuscript actually makes. Pairing a p-value that tests H0: log2FC = 0 with
# a separate |log2FC| filter asks two different questions and the FDR covers
# neither -- a point estimate can clear the line while its interval does not.
# 1.5-fold keeps the tested claim modest enough to retain a broad program.
PROGRAM_FC  <- 1.5                  # fold-change the response must exceed
PROGRAM_LFC <- log2(PROGRAM_FC)     # = 0.585, passed to results(lfcThreshold=)
PADJ        <- 0.05

source(file.path("scripts", "_common.R"))   # canonical gene universe

data_dir  <- file.path(getwd(), "data")
dds_4hrs  <- readRDS(file.path(data_dir, "dds_4hrs.rds"))
dds_24hrs <- readRDS(file.path(data_dir, "dds_24hrs.rds"))

# ---- Expression filter lives in 00, BEFORE DESeq() -----------------------------
# 00 applies >=10 counts in >=3 samples across ALL samples INCLUDING untreated, so
# nothing further is needed here. Filtering on stimulated samples alone would drop
# genes expressed in untreated and silenced by IL-4, biasing against repressed genes.
cat(sprintf("[01] Genes carried over from 00 (filtered there): %d\n", nrow(counts(dds_4hrs))))

# ---- Extract the four contrasts (each vs control 'ct') ----
contrasts <- list(p100_4h  = list(dds_4hrs,  "p100"), c100_4h  = list(dds_4hrs,  "c100"),
                  p100_24h = list(dds_24hrs, "p100"), c100_24h = list(dds_24hrs, "c100"))
res <- NULL
for (lab in names(contrasts)) {
  r <- as.data.frame(results(contrasts[[lab]][[1]],
                             contrast = c("condition", contrasts[[lab]][[2]], "ct"),
                             alpha = PADJ, lfcThreshold = PROGRAM_LFC))
  tab <- data.frame(ensembl_gene_id = rownames(r), stringsAsFactors = FALSE)
  tab[[paste0("log2fc_", lab)]] <- r$log2FoldChange
  tab[[paste0("padj_",   lab)]] <- r$padj
  res <- if (is.null(res)) tab else merge(res, tab, by = "ensembl_gene_id", all = TRUE)
}

# ---- Direct pulse-vs-continuous contrast ---------------------------------------
# The four contrasts above are each vs untreated. Any statement about the pulse
# "deviating from" continuous must rest on THIS contrast, never on comparing two
# independently thresholded vs-control lists -- that is the
# difference-between-significant-and-not-significant error (Gelman & Stern 2006).
# Naming: `pvc` = pulse vs continuous. A POSITIVE log2fc_pvc means the pulse is
# HIGHER than continuous; padj_pvc < 0.05 means the two conditions differ.
# NOTE this contrast deliberately keeps H0: log2FC = 0 -- no lfcThreshold. The
# question here is "do the two conditions differ at all", not "do they differ by
# more than 1.5-fold". Do NOT propagate PROGRAM_LFC here.
for (lab in c("pvc_4h", "pvc_24h")) {
  dds <- if (lab == "pvc_4h") dds_4hrs else dds_24hrs
  # alpha MUST be passed: it sets the target FDR that DESeq2's independent-filtering
  # step optimises against. Leaving the default 0.1 while every consumer of
  # padj_pvc_* applies 0.05 silently changes which genes are called significant.
  r <- as.data.frame(results(dds, contrast = c("condition", "p100", "c100"), alpha = PADJ))
  tab <- data.frame(ensembl_gene_id = rownames(r), stringsAsFactors = FALSE)
  tab[[paste0("log2fc_", lab)]] <- r$log2FoldChange
  tab[[paste0("padj_",   lab)]] <- r$padj
  res <- merge(res, tab, by = "ensembl_gene_id", all = TRUE)
}

# ---- Gene symbols: biomaRt pinned to Ensembl 110 (matches the GRCm39.110 ----
#      alignment genome). Cached to CSV so reruns are offline / reproducible.
ENSEMBL_VERSION <- 110
map_file <- file.path(data_dir, "ensembl110_gene_map.csv")
if (file.exists(map_file)) {
  gene_map <- read.csv(map_file, stringsAsFactors = FALSE)
  cat(sprintf("[01] Loaded cached Ensembl-%d gene map (%d genes)\n", ENSEMBL_VERSION, nrow(gene_map)))
} else {
  cat(sprintf("[01] Querying biomaRt (Ensembl %d) for gene symbols...\n", ENSEMBL_VERSION))
  mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = ENSEMBL_VERSION)
  gene_map <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id", values = unique(res$ensembl_gene_id), mart = mart)
  gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
  write.csv(gene_map, map_file, row.names = FALSE)
  cat(sprintf("[01] Cached gene map -> %s (%d genes)\n", map_file, nrow(gene_map)))
}
res$gene_symbol <- gene_map$mgi_symbol[match(res$ensembl_gene_id, gene_map$ensembl_gene_id)]
res <- res[, c("ensembl_gene_id", "gene_symbol",
               "log2fc_p100_4h", "padj_p100_4h", "log2fc_p100_24h", "padj_p100_24h",
               "log2fc_c100_4h", "padj_c100_4h", "log2fc_c100_24h", "padj_c100_24h",
               "log2fc_pvc_4h", "padj_pvc_4h", "log2fc_pvc_24h", "padj_pvc_24h")]
write.csv(res, file.path(data_dir, "master_expressed_genes.csv"), row.names = FALSE)

# ---- OR-gate DEG selection (stage 1) ----
# padj ALONE -- the fold-change requirement is inside the null (see PROGRAM_LFC).
#
# SCOPE OF THE 5%: BH is applied within each contrast, so the FDR is controlled per
# contrast, NOT across the union taken here. The program is the union of four
# individually controlled sets and does not itself carry a 5% guarantee. Nothing
# downstream depends on it doing so, but do not describe it as a 5%-FDR gene set.
# Controlling the union would need a per-gene omnibus p-value across the four
# contrasts (Simes or similar), adjusted across genes -- which would change program
# membership and therefore every downstream figure.
hit <- function(lfc, padj) (padj < PADJ) %in% TRUE
pulse <- hit(res$log2fc_p100_4h, res$padj_p100_4h) | hit(res$log2fc_p100_24h, res$padj_p100_24h)
cont  <- hit(res$log2fc_c100_4h, res$padj_c100_4h) | hit(res$log2fc_c100_24h, res$padj_c100_24h)
program <- res[pulse | cont, ]

# ---- Restrict to protein-coding genes (see scripts/_common.R) ----
program <- program[is_annotatable_gene(program$ensembl_gene_id, data_dir), ]

# NA log2FC means DESeq2 could not estimate the effect in that contrast; it does NOT
# mean "no change". Left as NA: 0 would be indistinguishable from a real zero.
# Every consumer guards with !is.na().

write.csv(program, file.path(data_dir, "deg_program.csv"), row.names = FALSE)
cat(sprintf("[01] Stage-1 DEG program (formal test H0: |log2FC| <= %.3f i.e. %.1f-fold, padj < %g): %d genes\n",
            PROGRAM_LFC, PROGRAM_FC, PADJ, nrow(program)))
cat(sprintf("       pulse-significant: %d   continuous-significant: %d\n",
            sum(pulse, na.rm = TRUE), sum(cont, na.rm = TRUE)))
