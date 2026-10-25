# snRNA aging microglia

Reproducible, path-independent R code for the single-nucleus and single-cell transcriptomic panels in the manuscript **Microglial GPR183 drives an aging-amplified state in neuropathic pain**.

The scripts are organized by the V7 artwork. They do not contain local computer or server paths. Each script accepts an input RDS file and an output directory:

```bash
Rscript R/Figure1/figure1_snrna_panels.R data/annotated_spinal_cord.rds results/Figure1
```

Every plotting script writes vector PDF, 600-dpi PNG and CSV source-data files where applicable.

## Scope and evidence boundaries

- Cell-level violin plots reproduce the manuscript display, but nuclei are not independent biological replicates. Sample-level summaries are exported for replicate-aware inference.
- Cell-type fractions are calculated within each biological sample before group summaries.
- Seurat `FindMarkers` output is suitable for exploratory cell-level visualization. Confirmatory differential-expression testing should use pseudobulk counts with biological sample as the unit.
- CellChat networks are descriptive unless independently replicated libraries are available. The scripts do not add artificial error bars or P values.
- Figure 2f is spatial transcriptomics, not snRNA-seq, and is intentionally outside the core code collection.
- Figures 4 and 5, Supplementary Figure 1 and the formal Supplementary Figure 5 contain no snRNA-seq panels.

## Figure to code map

| Figure panels | Analysis | Script |
|---|---|---|
| Fig. 1d-e | Cell-type UMAP and microglial marker FeaturePlot | `R/Figure1/figure1_snrna_panels.R` |
| Fig. 1f | Per-sample microglial abundance | `R/Figure1/figure1_snrna_panels.R` |
| Fig. 1g, i, k | Priming, pathological and homeostatic scores | `R/Figure1/figure1_snrna_panels.R` |
| Fig. 1h, j | Microglial program DotPlots | `R/Figure1/figure1_snrna_panels.R` |
| Fig. 2b-d | Integrated microglial atlas and two-axis state model | `R/Figure2/figure2_snrna_panels.R` |
| Fig. 2e | Gpr183 expression by spinal cell type | `R/Figure2/figure2_snrna_panels.R` |
| Fig. 2i-j | Temporal Gpr183 expression and Gpr183-status priming score | `R/Figure2/figure2_snrna_panels.R` |
| Fig. 3d-e | Gpr183 cKO microglial programs and scores | `R/Figure3/figure3_snrna_panels.R` |
| Fig. 3i-j | Agonist-associated programs and scores | `R/Figure3/figure3_snrna_panels.R` |
| Fig. 6a-d | M3 state-space shift, scores, DEGs and programs | `R/Figure6/figure6_state_and_deg.R` |
| Fig. 6e-f | Microglia-outgoing CellChat networks | `R/Figure6/figure6_cellchat.R` |
| Fig. 6g | Monocle3 trajectory and pseudotime | `R/Figure6/figure6_trajectory.R` |
| Supplementary Fig. 2a-i | QC, annotation, composition, DEGs and scores | `R/Supplementary/supplementary_figure2.R` |
| Supplementary Fig. 3a-g | Integrated atlas validation, time course and priming ranking | `R/Supplementary/supplementary_figure3.R` |
| Supplementary Fig. 4a-e | Independent dataset validation | `R/Supplementary/supplementary_figure4_independent_dataset.R` |
| Supplementary Fig. 6a-k | Gpr183-positive microglial programs | `R/Supplementary/supplementary_figure6_Gpr183_positive.R` |
| Supplementary Fig. 7a-e, g | Perturbation state space, enrichment inputs and network | `R/Supplementary/supplementary_figure7.R` |
| Supplementary Fig. 7f | Young-group trajectories | `R/Figure6/figure6_trajectory.R` with the requested conditions |

## Expected Seurat metadata

The exact columns required by each script are validated at runtime. The common columns are:

| Column | Meaning | Example values |
|---|---|---|
| `sample_id` | Biological library or animal identifier | `Mouse01` |
| `condition` | Experimental group | `Sham_Young`, `SNI_Aged`, `SNI_Aged_M3` |
| `cell_type` | Major spinal-cord cell type | `Microglia`, `Astrocytes` |
| `microglial_state` | Integrated microglial state | `Homeostatic`, `Inflammatory` |
| `priming_score` | Frozen priming-axis score | numeric |
| `pathological_score` | Frozen pathological-axis score | numeric |
| `timepoint` | Time after injury in the temporal dataset | `Naive`, `dpi3`, `dpi14`, `dpi150` |
| `sex` | Sex, when present in a public dataset | `Male`, `Female` |

The analysis assumes normalized RNA expression is stored in the Seurat RNA assay `data` layer and raw counts remain available in `counts`.

## Installation

```bash
conda env create -f environment.yml
conda activate snrna-aging-microglia
```

CellChat is installed separately:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("jinworks/CellChat")
```

## Suggested execution order

1. Prepare and annotate the full spinal-cord object with `R/00_preprocess_and_annotate.R` or provide an equivalent object.
2. Freeze the gene lists and state-model definitions with `R/01_construct_state_model.R` before regenerating manuscript panels.
3. Run Figure 1 and Supplementary Figure 2 from the full annotated object.
4. Run Figure 2 and Supplementary Figures 3 and 6 from the integrated microglial atlas.
5. Run Figure 3 from the perturbation-aware microglial object.
6. Run Figure 6 and Supplementary Figure 7 from the pharmacological projection object, precomputed CellChat list and monocle3 object.

GO and ranked enrichment used in Figure 1g and Supplementary Figures 2g, 3g and 7d-e can be generated with `R/02_enrichment_analysis.R`. The tested-gene universe must be supplied for over-representation analysis; it should not be replaced by the full genome.

## Reproducibility notes

- Set random seeds before stochastic PCA, UMAP, clustering, module scoring and trajectory learning.
- Record package versions with `sessionInfo()` for each release.
- Do not recompute a frozen score model while reproducing a manuscript figure; store the defining genes and cutoffs with the analysis object.
- Public datasets should be documented with accession, library-level metadata and download date.
- Do not commit identifiers, credentials, local absolute paths or large raw matrices.

## Tests

```bash
Rscript tests/test_helpers.R
```
