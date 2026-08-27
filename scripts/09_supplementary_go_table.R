# ===============================================================================
# SUPPLEMENTARY TABLE — GO enrichment of the continuous IL-4 program
# ===============================================================================
# A supplement covering ONLY the stage-1 continuous program, induced and repressed.
# The shortfall-vs-matched comparison is a different question and a separate file;
# adjusted p-values are not comparable across separately-run enrichments.
#
# Packaging step, not an analysis: everything here is read back from the workbook
# that go_program.R writes, so the table cannot disagree with the figures. Run
# go_program.R first.
#
# Each sheet lists EVERY significant term, not only the ones drawn, with two audit
# columns:
#   figure_representative  the term survived the gene-overlap redundancy collapse
#   in_figure              the term is actually drawn in that column of the panel
# and, where the term carries enough classified genes, the preservation statistics
# from the companion panel, so one file answers both "what is the program" and
# "how much of it is not significantly below continuous".
#
# Input : results/go_program_results.xlsx   (from go_program.R)
# Output: results/SupplementaryTable_GO_program.xlsx
# ===============================================================================

suppressPackageStartupMessages(library(openxlsx))

in_file  <- file.path(getwd(), "results", "go_program_results.xlsx")
out_dir  <- file.path(getwd(), "results")
if (!file.exists(in_file))
  stop("Missing ", in_file, " - run scripts/04_go_program.R first.")

sheets <- getSheetNames(in_file)
want   <- c("Program_ up" = "GO induced", "Program_ down" = "GO repressed")
missing <- setdiff(names(want), sheets)
if (length(missing))
  stop("go_program_results.xlsx lacks expected sheet(s): ", paste(missing, collapse = ", "))

pres <- if ("Preservation" %in% sheets) read.xlsx(in_file, "Preservation") else NULL

# Enrichment run directly on the PULSE gene sets. The power caveat and the
# interpretation are stated in full on the Legend sheet.
pulse <- list(up   = if ("Pulse_ up"   %in% sheets) read.xlsx(in_file, "Pulse_ up")   else NULL,
              down = if ("Pulse_ down" %in% sheets) read.xlsx(in_file, "Pulse_ down") else NULL)

# Terms that survived the redundancy collapse, recomputed with the SAME rule the
# figure uses so the audit column cannot drift from the panel.
source(file.path("scripts", "_go_helpers.R"))

tidy_one <- function(sheet, direction) {
  d <- read.xlsx(in_file, sheet)
  if (!nrow(d)) return(NULL)
  kept <- prune_redundant(d, OVERLAP_CUT)$Description
  out <- data.frame(
    GO_ID       = d$ID,
    Term        = d$Description,
    gene_ratio  = d$GeneRatio,      # term genes in the program / program genes tested
    bg_ratio    = d$BgRatio,        # term genes in the background / background size
    n_genes     = d$Count,
    pvalue      = signif(d$pvalue,   3),
    padj        = signif(d$p.adjust, 3),
    qvalue      = signif(d$qvalue,   3),
    figure_representative = d$Description %in% kept,
    in_figure   = if ("in_figure" %in% names(d)) d$in_figure else NA,
    stringsAsFactors = FALSE)
  if (!is.null(pres)) {
    p <- pres[pres$direction == direction, ]
    i <- match(out$Term, p$Term)
    out$pct_genes_not_below_continuous <- round(p$pct[i], 1)
    out$n_genes_not_below_continuous   <- p$n_preserved[i]
    out$n_genes_classified            <- p$n_classified[i]
  }
  pd <- pulse[[direction]]
  if (!is.null(pd)) {
    i <- match(out$Term, pd$Description)
    out$enriched_in_pulse_genes <- !is.na(i)
    out$padj_pulse_gene_set     <- signif(pd$p.adjust[i], 3)
    out$n_genes_pulse           <- pd$Count[i]
  }
  out$genes <- d$geneID
  out[order(out$padj), ]
}

tabs <- Map(tidy_one, names(want), c("up", "down"))
names(tabs) <- unname(want)

legend <- data.frame(Field = c(
  "Contents", "Gene sets", "Test", "Background", "Redundancy collapse",
  "figure_representative", "in_figure", "pct_genes_not_below_continuous",
  "enriched_in_pulse_genes / padj_pulse_gene_set / n_genes_pulse",
  "gene_ratio / bg_ratio", "Software", "Source"),
  Description = c(
  sprintf("GO Biological Process over-representation for the continuous IL-4 program: %d significant terms for induced genes and %d for repressed, listed in full.",
          nrow(tabs[["GO induced"]]), nrow(tabs[["GO repressed"]])),
  "Induced = genes significantly UP under continuous IL-4 at 4h and/or 24h; repressed = significantly DOWN. Significance is padj < 0.05 from a DESeq2 test of H0: |log2FC| <= log2(1.5); the two directions are tested separately.",
  "clusterProfiler::enrichGO, ontology BP, Benjamini-Hochberg adjusted, then simplify(cutoff = 0.7, by = 'p.adjust') to remove semantically redundant terms.",
  "All expressed protein-coding genes (Ensembl 110 biotype) that map to an Entrez ID, NOT the whole genome.",
  sprintf("simplify() compares GO graph topology and does not see shared genes, so terms sharing more than %.0f%% of their genes are additionally collapsed. The most significant member of each redundant cluster is kept as its representative; ties are broken toward the more specific term.", 100 * OVERLAP_CUT),
  "TRUE if the term survived that gene-overlap collapse as the representative of its group. FALSE terms are real enrichments, not errors - they are redundant with a listed term.",
  sprintf("TRUE if the term is drawn in that column of the figure (top %d representatives by adjusted p-value).", TOPN),
  "Of the term's program genes, the percentage whose pulsed response is NOT significantly below continuous at 24h, judged by the direct pulsed-vs-continuous contrast. Note this is the complement of a significance test, not an equivalence test: a gene counts here because a difference was not detected, which reflects statistical power as well as similarity. Blank where the term has too few classified genes for the rate to be precise.",
  "The same GO test run on the PULSE gene sets instead of the continuous ones. TRUE means the term is also significantly enriched when only pulse-regulated genes are tested. Read with care: the pulse sets are smaller than the continuous ones, so FALSE partly reflects reduced power and adjusted p-values are NOT comparable between the two columns. The interpretable signal is the asymmetry between arms.",
  "gene_ratio = term genes found in this gene set / gene set size. bg_ratio = term genes in the background / background size.",
  sprintf("R %s.%s, clusterProfiler, org.Mm.eg.db.", R.version$major, R.version$minor),
  "Generated by scripts/09_supplementary_go_table.R from results/go_program_results.xlsx"),
  stringsAsFactors = FALSE)
# Field and Description are matched positionally; a length mismatch silently
# shifts every description onto the wrong field, so fail loudly instead.
stopifnot(nrow(legend) == 12L, !anyNA(legend$Field), !anyNA(legend$Description))

wb <- createWorkbook()
addWorksheet(wb, "Legend"); writeData(wb, "Legend", legend)
setColWidths(wb, "Legend", cols = 1:2, widths = c(30, 130))
addStyle(wb, "Legend", createStyle(wrapText = TRUE, valign = "top"),
         rows = 2:(nrow(legend) + 1), cols = 2, gridExpand = TRUE)
for (nm in names(tabs)) {
  d <- tabs[[nm]]; if (is.null(d)) next
  addWorksheet(wb, nm); writeData(wb, nm, d)
  setColWidths(wb, nm, cols = 1:ncol(d), widths = "auto")
  setColWidths(wb, nm, cols = which(names(d) == "genes"), widths = 60)
  freezePane(wb, nm, firstActiveRow = 2, firstActiveCol = 3)
  addStyle(wb, nm, createStyle(textDecoration = "bold"), rows = 1,
           cols = 1:ncol(d), gridExpand = TRUE)
}
f <- file.path(out_dir, "SupplementaryTable_GO_program.xlsx")
saveWorkbook(wb, f, overwrite = TRUE)
cat(sprintf("Saved: results/SupplementaryTable_GO_program.xlsx\n"))
for (nm in names(tabs))
  cat(sprintf("  %-14s %4d significant terms; %3d survive the collapse; %2d drawn in the figure\n",
              nm, nrow(tabs[[nm]]), sum(tabs[[nm]]$figure_representative),
              sum(tabs[[nm]]$in_figure %in% TRUE)))
