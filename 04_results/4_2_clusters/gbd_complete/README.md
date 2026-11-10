# §4.2 — GBD Hungary complete DTW clustering: code, audited results and figures

**This folder supersedes the older one-factor-only script one directory above.** It is based on the existing 2026-11-10 full experiment archive, which had seven independently fitted single-population models × two sexes × nine causes (126 fitted models), and 14 sex-model DTW clustering analyses.

## Reproduction order and rationale

1. Fit each of seven Table-1 mean structures separately on actual GBD `deaths_gbd` and `exposure_proxy`. Do **not** substitute `fitted_count`.
2. Extract cause-specific period factors: **factor 2** for two-factor Poisson/Binomial LC, **factor 1** for all others. Canonicalize LC factors with centered temporal indices and SVD, and note that rotation changes individual factor interpretations.
3. z-standardize each factor across the **training** years (1980–2015) within its cause and sex. This ensures scale does not dominate a similarity measure, while preserving temporal shape.
4. Compute pairwise symmetric2 DTW with squared local cost, diagonal weight 2, horizontal/vertical weight 1, and `sqrt(accumulated_cost/(N+M))`. No time window in main analysis; five-year window used for sensitivity.
5. Perform average-linkage hierarchical clustering and exact-enumeration k-medoids on the same DTW matrix, for K=2…8. Choose max average silhouette in this project. Compare seven auxiliary internal CVIs, medoid starting points and rotated LC factors.
6. Make separate Figure 5 (male second LC factor alignment), Figure 6 (female), Figure 7 (male CBD first factor grouped), Figure 8 (male two-factor Binomial LC second-factor grouped), and Tables 3–4 analogues. These are **GBD analogues**, not numerical reproductions of the original HCSO figures.

## Figures — view directly in GitHub

**Figure 5: male Poisson LC2, CVD and DIG.** [SVG](figures/Figure5_GBD_Male_DTW.svg). Distance **1.2283**, rank **35/36** among male cause pairs by smallest distance. **Unlike the published paper, these two GBD causes are *not* a close pair.** The smallest observed male pair is DIG–INJ, distance **0.2381**.

**Figure 6: female Poisson LC2, CVD and DKD.** [SVG](figures/Figure6_GBD_Female_DTW.svg). Distance **0.2044**, rank **1/36**; however DKD in GBD is not the same definition as original ICD chapter IV.

**Figure 7: male CBD first-factor groups.** [SVG](figures/Figure7_GBD_Male_CBD.svg). Two clusters: CVD+NEO+DIG+CRD+INJ+SHV+RIT versus NEU+DKD, silhouette **0.7144**. Original male CBD Figure 7 reported **six** clusters on HCSO ICD groups, not two.

**Figure 8: male two-factor Binomial Lee–Carter second-factor groups.** [SVG](figures/Figure8_GBD_Male_Binomial_LC2.svg). Two groups: CVD+NEU and seven other causes; silhouette **0.5154**. Original male Binomial LC2 Figure 8 reported **five** clusters; the GBD grouping differs from the original and from our CBD Figure 7.\n\nThe SVGs are plotted from verified, standardized factor values rounded to 0.01 for compact text-based GitHub storage; full-resolution source CSVs and PNGs are in the referenced result archive, with all output computed before rounding.

## Audited aggregate results

| Model | Male K / Silhouette | Female K / Silhouette |
| --- | --- | --- |
| Poisson LC1 | 2 / 0.638 | 2 / 0.705 |
| Binomial LC1 | 2 / 0.638 | 2 / 0.705 |
| Poisson LC2 | 2 / 0.513 | 3 / 0.414 |
| Binomial LC2 | 2 / 0.515 | 3 / 0.413 |
| CBD | 2 / 0.714 | 2 / 0.541 |
| Plat2 | 2 / 0.700 | 2 / 0.522 |
| Plat3 | 2 / 0.726 | 2 / 0.578 |

See [cluster_summary_selected.csv](results/cluster_summary_selected.csv) and [cluster_membership.csv](results/cluster_membership.csv) for the actual numbers and each of the 126 memberships. Each row is a different sex-model factor-space clustering; group numbers do not have a global identity across models.

### Important stability limitations

- Rotating LC2 factors by 15° preserves the fitted surface to numerical precision (max error **2.22e−16**) yet changes second-factor clustering (Poisson LC2 ARI **0.182 male**, **−0.050 female**). This is a critical identification limit.
- A five-year Sakoe–Chiba window changes female LC2 from 3 clusters to 2 (ARI **0.353**); female CBD and Plat2 change from 2 to 3.
- Excluding DKD in the female sensitivity check leaves memberships of the remaining causes unchanged in this particular result (ARI=1).
- The original male HCSO results were six clusters for most models, five for LC2; female original clusters were two. We get male K=2 throughout and female K=3 for LC2. Differences cannot be attributed to a single factor without matched data/age/category/DTW and K-selection controlled studies.
- A **silhouette score is not a forecasting accuracy measure**. No §4.3 multi-population forecasts or held-out prediction gains are claimed.

## How to reproduce from the CSV

**RStudio users**: [run_4_2_complete_RStudio_launcher.R](run_4_2_complete_RStudio_launcher.R) lets you use graphical file selectors, preview the exact input/output paths, then invoke the following **Python** pipeline with a single confirmation. It is an R front-end to Python/SciPy, **not** an independently estimated R implementation. Install the listed Python packages first; it writes outputs next to the scripts in `results/`.


```bash
python -m pip install numpy pandas scipy matplotlib scikit-learn
OPENBLAS_NUM_THREADS=1 python run_all_models.py /path/to/fitted_rates.csv results
OPENBLAS_NUM_THREADS=1 python verify_and_plot.py results
OPENBLAS_NUM_THREADS=1 python complete_cvi.py results
OPENBLAS_NUM_THREADS=1 python plot_paper_figures.py
```

Input is the previously uploaded `fitted_rates(1).csv`, 3,888 rows (36 years × 6 ages × 9 causes × 2 sexes). The program verifies row count, years and selected causes. Full original CSV is not stored here, because the user's GBD file is not in the repository. Main outputs include all fitted parameters, 8,424 period-factor values, 4,536 standardized selected-factor values, 504 distances, 196 validity candidates, membership tables, alignment paths, sensitivity files and plots.

**Dataset and inference warning:** GBD's estimated deaths are not integer registration counts; the Poisson/Binomial structures are working pseudo-likelihoods here. Binomial exposure `E0 ≈ Ec + D/2` is estimated, not observed. Numeric optimization checks of 126/126 models are **not** evidence for the models' scientific validity or future accuracy.

**Input provenance**: [manifest.json](results/manifest.json); [validation.json](results/validation.json). The validated pre-existing complete archive is available via the project's [ChatGPT Library](https://chatgpt.com/library?file_id=libfile_9a69c3210b3c81918a86c5d8ac395012&account_id=personal), not duplicated as a binary ZIP in GitHub.

**Paper**: [Varga (2025), §4.2, Figures 5–8, Tables 3–4](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf#page=14). See also the detailed [Notion analysis](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2).
