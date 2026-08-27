# ===============================================================================
# _common.R — definitions shared across the analysis scripts
# ===============================================================================
# Sourced by the numbered analysis scripts and by _go_helpers.R. Not standalone.
#
# The enrichment background depends on the gene universe as much as on the gene set,
# so the filter is defined here only.
# ===============================================================================

# ---------------------------------------------------------------------------
# GENE UNIVERSE: protein-coding, by Ensembl gene_biotype.
#
# Filter on BIOTYPE, never on gene symbols. A symbol regex (^Gm[0-9]+, *Rik,
# ^LOC and similar) is wrong in both directions: it drops real protein-coding genes
# carrying provisional Riken names, and keeps non-coding entries whose symbols
# happen not to match. Gm2a is the GM2 activator, not "Gm" + "2a".
#
# Biotype comes from data/ensembl110_gene_map.csv, parsed from the Ensembl 110
# GTF -- the same release the reads were aligned to, so there is no version skew
# and no network dependency.
#
# protein_coding ONLY; IG_*/TR_* immunoglobulin and TCR genes are excluded and are
# negligibly expressed in BMDM.
PROTEIN_CODING_ONLY <- TRUE

gene_biotype_map <- function(data_dir) {
  f <- file.path(data_dir, "ensembl110_gene_map.csv")
  if (!file.exists(f)) stop("Missing ", f)
  m <- read.csv(f, stringsAsFactors = FALSE)
  if (!"gene_biotype" %in% colnames(m))
    stop("ensembl110_gene_map.csv has no gene_biotype column - rebuild it from the Ensembl 110 GTF")
  m
}

# TRUE for genes to KEEP.
is_annotatable_gene <- function(ensembl_ids, data_dir) {
  m <- gene_biotype_map(data_dir)
  bt <- m$gene_biotype[match(ensembl_ids, m$ensembl_gene_id)]
  !is.na(bt) & bt == "protein_coding"
}

# The canonical expressed-gene universe: protein-coding expressed genes.
# Every enrichment background derives from this.
expressed_symbols <- function(data_dir) {
  f <- file.path(data_dir, "master_expressed_genes.csv")
  if (!file.exists(f)) stop("Missing ", f, " - run scripts/01 first.")
  d <- read.csv(f, stringsAsFactors = FALSE)
  s <- d$gene_symbol[is_annotatable_gene(d$ensembl_gene_id, data_dir)]
  unique(s[!is.na(s) & s != ""])
}

# ---------------------------------------------------------------------------
# Marker genes labelled on the descriptive panels. Defined once so 02 and 03 label
# the same set.
#   canonical M2 induction   Arg1, Retnla, Chil3
#   negative regulators      Klf4, Egr2
#   repressed                Cx3cr1, and the ISGs Ifit3b, Cxcl10
MARKER_GENES <- c("Arg1", "Retnla", "Chil3", "Klf4", "Cx3cr1", "Ifit3b", "Egr2", "Cxcl10")

# ---------------------------------------------------------------------------
# Partition of the continuous IL-4 program by what the PULSE did to each gene,
# using the DIRECT pulse-vs-continuous contrast (padj_pvc_*, produced by 01).
# Used by 04 to compute the per-pathway preservation percentages.
#
#   matched               indistinguishable from continuous at both timepoints
#   delayed (caught up)   below continuous at 4h but NOT at 24h
#   reduced at 24h        below continuous at 24h, still responds vs untreated
#   absent                below continuous at 24h, no significant response at all
#
# `band` collapses these into the three sets used for enrichment. Use
# persistent_shortfall (reduced + absent) for "which functions are under-reproduced";
# caught_up genes recover by 24h and do not belong in that set.
# "Below continuous" is always evaluated IN THE CONTINUOUS DIRECTION, so a gene the
# pulse drives further than continuous is never counted as a shortfall.
# ---------------------------------------------------------------------------
# T = extra point-estimate threshold. Default 0: padj already comes from the formal
# 1.5-fold test in 01, so significance plus the right sign IS the stage-1 program.
# Callers may raise T to subset on a larger point estimate, but nothing here does.
pulse_fate_bands <- function(data_dir, T = 0, padj = 0.05) {
  f <- file.path(data_dir, "deg_program.csv")
  if (!file.exists(f)) stop("Missing ", f, " — run scripts/01 first.")
  d <- read.csv(f, stringsAsFactors = FALSE)
  sg <- function(p) !is.na(p) & p < padj
  out <- do.call(rbind, lapply(c("Up", "Down"), function(dir) {
    s   <- if (dir == "Up") 1 else -1
    c4  <- d$log2fc_c100_4h * s;  c24 <- d$log2fc_c100_24h * s
    p4  <- d$log2fc_p100_4h * s;  p24 <- d$log2fc_p100_24h * s
    inC <- (c4 >= T & sg(d$padj_c100_4h)) | (c24 >= T & sg(d$padj_c100_24h))
    short4  <- sg(d$padj_pvc_4h)  & (d$log2fc_pvc_4h  * s) < 0
    short24 <- sg(d$padj_pvc_24h) & (d$log2fc_pvc_24h * s) < 0
    presp   <- (p4 > 0 & sg(d$padj_p100_4h)) | (p24 > 0 & sg(d$padj_p100_24h))
    fate <- ifelse(!short4 & !short24, "matched",
            ifelse(short4 & !short24,  "delayed (caught up)",
            ifelse(presp,              "reduced at 24h", "absent")))
    data.frame(gene_symbol = d$gene_symbol[inC], ensembl_gene_id = d$ensembl_gene_id[inC],
               direction = dir, fate = fate[inC],
               ratio = (pmax(p4, p24, na.rm = TRUE) / pmax(c4, c24, na.rm = TRUE))[inC],
               stringsAsFactors = FALSE)
  }))
  out$band <- ifelse(out$fate %in% c("reduced at 24h", "absent"), "persistent_shortfall",
              ifelse(out$fate == "delayed (caught up)",           "caught_up", "matched"))
  out
}
