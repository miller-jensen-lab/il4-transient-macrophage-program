# ===============================================================================
# _go_helpers.R — shared GO over-representation machinery
# ===============================================================================
# Sourced by the two GO scripts. Not a standalone script.
#
# Enrichment machinery for 04: the expressed-gene background, the redundancy
# collapse, and the panel builders.
#
# Kept out of _common.R so scripts that do no enrichment don't load
# clusterProfiler and org.Mm.eg.db.
# ===============================================================================

suppressPackageStartupMessages({
  library(openxlsx); library(ggplot2)
  library(clusterProfiler); library(org.Mm.eg.db); library(AnnotationDbi)
})

TOPN <- 8   # top terms per group (by adjusted p-value) shown in the figure

# ---- Redundancy control for the figure -------------------------------------
# clusterProfiler::simplify() removes SEMANTICALLY similar terms, but semantic
# similarity is computed from GO graph topology and does NOT see shared genes: it
# keeps `chemotaxis` alongside `taxis` even when the two carry identical gene sets
# and identical p-values, and lowering its cutoff only discards other terms.
# OVERLAP_CUT additionally collapses by gene overlap, which is what actually makes
# the figure redundant.
OVERLAP_CUT <- 0.7   # collapse a term sharing > this fraction of its genes with a kept term

source(file.path("scripts", "_common.R"))   # canonical gene universe

data_dir <- file.path(getwd(), "data")
out_dir  <- file.path(getwd(), "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  m <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = symbols,
        keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
  unique(na.omit(unname(m)))
}

# Expressed-gene background, not whole-genome
bg_sym    <- expressed_symbols(data_dir)   # canonical universe; see scripts/_common.R
bg_entrez <- to_entrez(bg_sym)
cat(sprintf("[go] background universe: %d expressed genes (Entrez)\n", length(bg_entrez)))

run_enrich <- function(symbols, label) {
  ez <- to_entrez(symbols)
  cat(sprintf("  %s: %d genes (%d Entrez)\n", label, length(symbols), length(ez)))
  ego <- enrichGO(gene = ez, universe = bg_entrez, OrgDb = org.Mm.eg.db, ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                  minGSSize = 10, maxGSSize = 500, readable = TRUE)
  # Warn LOUDLY on failure. A silent fallback here returns the unsimplified table,
  # which looks identical apart from carrying far more terms -- easy to miss.
  tryCatch(as.data.frame(simplify(ego, cutoff = 0.7, by = "p.adjust")),
           error = function(e) {
             warning(sprintf("simplify() FAILED for '%s' (%s); using unsimplified terms",
                             label, conditionMessage(e)), call. = FALSE)
             as.data.frame(ego)
           })
}

# ---- Selection for the figure: top-N by p.adjust, redundant terms collapsed ----
# Greedy walk in p.adjust order (ties -> more specific term first). A term
# overlapping an already-kept term by more than OVERLAP_CUT is dropped as
# redundant. Applied uniformly to every group; display order stays by p.adjust.
#
# The walk keeps the first term it meets in each redundancy cluster and never
# replaces it, so the sort order alone decides which label survives. Sorting by
# p.adjust means a term is only ever DROPPED as redundant, never promoted over a
# better-supported one.
prune_redundant <- function(d, cut) {
  if (!nrow(d)) return(d)
  d  <- d[order(d$p.adjust, d$Count), ]   # most significant first; ties -> more specific
  gl <- strsplit(d$geneID, "/")
  keep <- integer(0)
  for (i in seq_len(nrow(d))) {
    redundant <- FALSE
    for (k in keep) {
      ov <- length(intersect(gl[[i]], gl[[k]])) / min(length(gl[[i]]), length(gl[[k]]))
      if (ov > cut) { redundant <- TRUE; break }
    }
    if (!redundant) keep <- c(keep, i)
  }
  d[sort(keep), ]
}

select_terms <- function(d) {
  if (!nrow(d)) return(character())
  p <- prune_redundant(d, OVERLAP_CUT)
  cat(sprintf("    (%d significant -> %d after redundancy collapse; showing %d)\n",
              nrow(d), nrow(p), min(TOPN, nrow(p))))
  head(p$Description, TOPN)
}

ratio <- function(x) { p <- as.numeric(strsplit(x, "/")[[1]]); p[1] / p[2] }

# Canonical term order for a set of enrichment groups, returned as ggplot factor
# LEVELS (levels[1] renders at the BOTTOM). A term is placed by the group in which
# it is most significant, then by that significance, most significant at the top.
# Shared by dotplot_panel and preservation_panel (05) so a term appears at the same
# height in both.
term_order <- function(egos, terms) {
  sel <- unique(unlist(terms)); rows <- list()
  for (term in sel) for (nm in names(egos)) {
    d <- egos[[nm]]; row <- d[d$Description == term, ]
    rows[[length(rows) + 1]] <- data.frame(Term = term, Group = nm,
      neglogp = if (nrow(row)) -log10(row$p.adjust[1]) else NA_real_,
      stringsAsFactors = FALSE) }
  pd <- do.call(rbind, rows)
  best <- do.call(rbind, lapply(split(pd, pd$Term), function(s) {
    s <- s[which.max(replace(s$neglogp, is.na(s$neglogp), -1)), ]
    data.frame(Term = s$Term[1], gi = match(s$Group[1], names(egos)),
               neglogp = s$neglogp[1], stringsAsFactors = FALSE) }))
  best[order(-best$gi, best$neglogp), "Term"]
}

dotplot_panel <- function(egos, terms, title, subtitle, file, width, height) {
  sel <- unique(unlist(terms))
  rows <- list()
  for (term in sel) for (nm in names(egos)) {
    d <- egos[[nm]]; row <- d[d$Description == term, ]
    rows[[length(rows) + 1]] <- data.frame(
      Term = term, Group = nm,
      neglogp   = if (nrow(row)) -log10(row$p.adjust[1]) else NA_real_,
      GeneRatio = if (nrow(row)) ratio(row$GeneRatio[1]) else NA_real_,
      stringsAsFactors = FALSE)
  }
  pd <- do.call(rbind, rows)
  pd$Term  <- factor(pd$Term, levels = term_order(egos, terms))
  for (gg in names(egos)) {
    lv <- rev(levels(pd$Term))
    lv <- lv[lv %in% pd$Term[pd$Group == gg & !is.na(pd$neglogp)]]
    cat(sprintf("  [yorder] %s [%s] : %s\n", basename(file), gg, paste(lv, collapse = " | ")))
  }
  pd$Group <- factor(pd$Group, levels = names(egos))
  # `drop = FALSE` on the x scale keeps a group with no enriched terms as an empty
  # column rather than silently dropping it -- but blank space is something the
  # reader has to interpret, so label those columns explicitly.
  empty <- names(egos)[vapply(terms, length, 1L) == 0]
  p <- ggplot(pd[!is.na(pd$neglogp), ], aes(Group, Term, size = GeneRatio, color = neglogp)) +
    geom_point() +
    scale_color_gradient(low = "steelblue", high = "darkred", name = expression(-log[10](p[adj]))) +
    scale_x_discrete(drop = FALSE) +
    scale_size_continuous(name = "gene ratio", range = c(1.5, 7)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  if (length(empty))
    p <- p + annotate("text", x = empty, y = (nlevels(pd$Term) + 1) / 2,
                      label = "no enriched\nterms", size = 2.5, colour = "grey45", lineheight = 0.95)
  ggsave(file, p, width = width, height = height)
  cat(sprintf("Saved: %s\n", file))
}

COLLAPSE_NOTE <- sprintf(paste("Top %d terms per group by adjusted p-value; expressed-gene background.",
                               "Terms sharing >%.0f%% of their genes collapsed to the most significant one.", sep = "\n"),
                         TOPN, 100 * OVERLAP_CUT)


# ---------------------------------------------------------------------------
# BAR PANEL — one facet per arm, bars = -log10(padj), annotation printed at the
# bar end. Used for the program figure in place of the shared-axis dot plot.
#
# One facet per arm so each panel orders terms by its own padj; dotplot_panel's
# shared y-axis forces a single cross-arm ordering, in which a term can sit above a
# stronger one because the other arm set their relative order. The cost is that
# terms enriched in both arms no longer align across facets; they remain in the
# workbook.
#
# NB the -log10(padj) range within an arm is narrow, so bar LENGTH does little work
# and the ordering carries the message -- a property of ranking by significance.
bar_panel <- function(d, file, width, height, base_size, title, subtitle,
                      xlab = expression(-log[10](p[adj]))) {
  d$key <- factor(paste0(as.integer(d$direction), "|", d$Term),
                  levels = paste0(as.integer(d$direction), "|", d$Term)[order(d$direction, d$neglogp)])
  xmax <- max(d$neglogp) * 1.40           # headroom for the printed annotation
  p <- ggplot(d, aes(x = neglogp, y = key, fill = direction)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = label), hjust = -0.12, size = base_size * 0.26, colour = "grey25") +
    facet_wrap(~ direction, ncol = 1, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^[0-9]+\\|", "", x)) +
    scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = setNames(c("#C0392B", "#2980B9"), levels(d$direction))) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_bw(base_size = base_size) +
    theme(plot.title = element_text(size = base_size * 0.95),
          plot.subtitle = element_text(size = base_size * 0.72, colour = "grey35"),
          panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold"))
  for (ff in levels(d$direction)) {
    lv <- rev(levels(d$key)); lv <- lv[lv %in% d$key[d$direction == ff]]
    cat(sprintf("  [yorder] %s [%s] : %s\n", basename(file), ff,
                paste(sub("^[0-9]+\\|", "", lv), collapse = " | ")))
  }
  ggsave(file, p, width = width, height = height)
  cat(sprintf("Saved: %s\n", file))
}
