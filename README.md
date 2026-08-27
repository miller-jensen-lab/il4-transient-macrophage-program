# IL-4 transient macrophage program

Analysis code for the paper

> _Transient IL-4 stimulation durably activates an alternative macrophage program and modulates the inflammatory response_ — currently under review.

Bone marrow-derived macrophages (BMDMs) were stimulated with IL-4 either **continuously** or as a brief **pulse**, and profiled by bulk RNA-seq at **4 h** and **24 h**. The central question the pipeline answers is *how much of the continuous IL-4 program a transient pulse reproduces, and when*.

## Quick start

```r
# from the repository root
for (f in list.files("scripts", "^[0-9]{2}_.*\\.R$", full.names = TRUE)) source(f)
```

Scripts are numbered in the order they run. Each writes to `data/` (intermediates) or `results/` (figures and tables); both are regenerable, and `results/` is not version-controlled.

**Once the environment is installed, the analysis runs offline.** The GEO submission is a single combined counts matrix that is not yet public, so `00` falls back to the per-timepoint `featureCounts_*_cleaned.txt` files committed here — the two sources produce byte-identical DESeq2 objects. No *analysis* step needs the network: the Ensembl annotation is cached in `data/ensembl110_gene_map.csv`, and the biomaRt call is a fallback that fires only if that file is deleted. Installing the packages in the first place does need network access — see Requirements.

## Pipeline

| # | Script | Does | Manuscript item |
|---|--------|------|-----------------|
| 00 | `counts_to_deseq2` | Counts → two DESeq2 objects (4 h, 24 h); expression filter | — (foundation) |
| 01 | `contrasts_and_program` | Six contrasts; defines the DEG **program** | — (foundation) |
| 02 | `volcano_venn` | Volcano panels; per-timepoint and 4-way Venns | **Figure 5A**, **Figure S4A** |
| 03 | `amplitude_regression` | What fraction of the continuous amplitude the pulse reproduces | **Figure 5B** |
| 04 | `go_program` | GO enrichment of the program; preservation printed on the bars | **Figure 5C** |
| 05 | `tf_regulon_amplitude` | Does any TF regulon deviate from the program-wide amplitude? | **Figure S4C** |
| 06 | `egr2_amplitude` | EGR2-dependent genes vs the rest, by amplitude ratio | **Figure S4B** |
| 07 | `tf_heatmap` | Curated a-priori TF panel, log2FC across all four contrasts | **Figure 6A** |
| 08 | `supplementary_deg_table` | Supplementary DEG workbook | **Data File S1** |
| 09 | `supplementary_go_table` | Supplementary GO workbook | **Data File S2** |
| 10 | `session_info` | Records the R and package versions that produced the results | — (provenance) |

`_common.R` (gene universe, fate bands, marker genes) and `_go_helpers.R` (enrichment machinery) are sourced by the numbered scripts and are not run directly.

Every numbered script produces something that appears in the paper, or is a
prerequisite for one that does. Scripts 02–04 describe how much of the program the
pulse reproduces; 05 and 06 ask whether any part of it is reproduced selectively;
07 shows a curated transcription-factor panel chosen a priori.

## How DEGs are defined

This is the decision most likely to be queried, so it is stated plainly.

A gene enters the **program** if it rejects **H₀: |log₂FC| ≤ log₂(1.5)** at BH-adjusted *p* < 0.05 in at least one of the four IL-4-versus-untreated contrasts. That null is tested directly by `DESeq2::results(lfcThreshold = log2(1.5))` — the TREAT approach of McCarthy & Smyth (2009) — rather than by filtering point estimates after a test against zero, so the reported FDR covers the fold-change claim itself.

A consequence: every program gene necessarily has |log₂FC| > 0.585 in whichever contrast made it significant. But that is a *property of the result*, not a filter, and no separate fold-change cut is applied downstream — re-applying one would reinstate the mismatch this design removes.

**What the 5% is, and is not.** BH adjustment is applied *within each contrast*, so the false-discovery rate is controlled per contrast and not across their union. The program is the union of four individually controlled gene sets, and it does not itself carry a 5% guarantee — 45% of it enters on the strength of a single contrast. Nothing downstream depends on the union being exactly 5%, but do not describe the program as a 5%-FDR gene set. Controlling the union would need a per-gene omnibus *p*-value across the four contrasts, adjusted across genes.

Other decisions worth knowing before modifying anything:

- **Two experiments, modelled separately.** 4 h and 24 h were run independently, each with its own untreated control (*n* = 3 per condition, IL-4 at 100 ng/mL). They are never pooled: the 4 h model carries a batch term, the 24 h model does not. Comparisons *between* timepoints are therefore across experiments — see the note in `03`.
- **Expression filter before dispersion estimation**: ≥ 10 counts in ≥ 3 samples (the smallest group size), applied as a union across timepoints so both objects keep identical rownames.
- **Protein-coding only**, by Ensembl 110 `gene_biotype` — the release the reads were aligned to. A gene-symbol regex was used previously and was wrong in both directions; read the comment in `_common.R` before changing this.
- **A direct pulse-versus-continuous contrast** is computed alongside the versus-untreated ones, and every claim about the pulse *deviating from* continuous rests on it rather than on comparing two independently thresholded lists. It keeps the conventional null of no difference, because it asks a different question.
- **Enrichment backgrounds** are the expressed protein-coding universe, never the whole genome, and are shared across the TF and GO analyses.

## Data

| File | Source |
|------|--------|
| `featureCounts_{4,24}hours_cleaned.txt` | This study. Raw counts, 56,941 Ensembl genes. Also deposited at GEO **GSE322520**, which is private until publication — these committed files are what make the pipeline runnable in the meantime. |
| `ensembl110_gene_map.csv` | Ensembl release 110 — gene symbols and biotypes, cached so the pipeline runs offline. |
| `collectri_human_network.csv` | [CollecTRI](https://omnipathdb.org/) via OmniPath: 42,990 signed TF–target interactions across 1,185 TFs. The human network is used because no mouse version exists; symbols are matched case-insensitively (`Stat6` → `STAT6`). |
| `genes_stat6_egr2_dependence.xlsx` | Dependence calls drawn from two Nagy-lab knockout series, one per column. **`STAT6KO_direction`** — Czimmerer Z et al., *Immunity* 48(1):75–90.e6, 2018 ([doi](https://doi.org/10.1016/j.immuni.2017.12.010), GEO **GSE106706**). **`EGR2KO_direction`** — Daniel B et al., *Genes Dev* 34(21–22):1474–1492, 2020 ([doi](https://doi.org/10.1101/gad.343038.120), GEO **GSE151015**). Used by `06`; the EGR2 column carries that script's result. |

`deg_program.csv` is committed so the downstream scripts can be run without regenerating the DESeq2 objects first.

## Requirements

**R 4.4.x — use 4.4.1 if you can.** The lockfile pins R 4.4.1 with Bioconductor 3.19,
which is an R-4.4 release; a later R pulls a different Bioconductor, and as noted below
that can change which GO terms are returned. `renv.lock` records the exact version and
source of all 159 packages in the dependency tree.

`renv::restore()` downloads those packages and **requires network access** (or a
pre-populated renv cache) the first time. Only the analysis itself is offline.

```r
# from the repository root — restores the recorded environment
renv::restore()
```

Opening R in this directory activates the project library automatically (via `.Rprofile`).
`renv/library/` is not version-controlled; `renv.lock` and `renv/activate.R` are.

**Why pinning matters here more than usual.** GO term-to-gene mappings live inside versioned
annotation snapshots (`org.Mm.eg.db`, `GO.db`). Which genes belong to a term changes between
Bioconductor releases, so re-running script 04 under a different release can return a different
list of terms, not merely different *p*-values. The differential expression results are far more
stable; the enrichment figures are the fragile part. Script 10 additionally writes
`results/sessionInfo.txt` on every run, so the outputs document their own provenance even if the
lockfile is not used.

## Outputs

Figures and tables are written to `results/`, which is gitignored — regenerate by running the scripts. `results/sessionInfo.txt` records the environment that produced them. The two workbooks intended as manuscript supplements are `SupplementaryTable_DEG_program.xlsx` (every DEG with log₂FC and adjusted *p* for all six contrasts) and `SupplementaryTable_GO_program.xlsx` (every enriched GO term for the program, with audit columns showing which are drawn in the figure and why). Both carry a `Legend` sheet so they are self-describing away from this repository.

## Citation

> Sayed Issa I., Kelsey I.A., et al. Transient IL-4 stimulation durably activates an alternative macrophage program and modulates the inflammatory response. *Journal TBD*. 2026.

A DOI will be added here once the paper is published.

## License

Released into the public domain under [The Unlicense](LICENSE).
