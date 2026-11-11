# §4.2｜七模型与死因聚类：完整 R 代码、实跑结果与图 5–8

**2026-11-10：已改为纯 R，并从同一份原始 GBD 输入端到端重跑。** 126 组拟合、14 组聚类全部完成；与迁移前存档的14组成员归属全部ARI=1。详见[实跑核验](R_MIGRATION_VALIDATION.md)和[Notion中文代码精读](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)。

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

**Figure 8: male two-factor Binomial Lee–Carter second-factor groups.** [SVG](figures/Figure8_GBD_Male_Binomial_LC2.svg). Two groups: CVD+NEU and seven other causes; silhouette **0.5154**. Original male Binomial LC2 Figure 8 reported **five** clusters; the GBD grouping differs from the original and from our CBD Figure 7. 

图形已由R重绘，使用完整精度参数，不提前舍入。PNG与SVG均由同一个R绘图函数生成。

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

See [cluster_summary.csv](results/cluster_summary.csv) and [cluster_membership.csv](results/cluster_membership.csv) for the actual numbers and each of the 126 memberships. Each row is a different sex-model factor-space clustering; group numbers do not have a global identity across models.

### Important stability limitations

- Rotating LC2 factors by 15° preserves the fitted surface to numerical precision (max error **2.22e−16**) yet changes second-factor clustering (Poisson LC2 ARI **0.182 male**, **−0.050 female**). This is a critical identification limit.
- A five-year Sakoe–Chiba window changes female LC2 from 3 clusters to 2 (ARI **0.353**); female CBD and Plat2 change from 2 to 3.
- Excluding DKD in the female sensitivity check leaves memberships of the remaining causes unchanged in this particular result (ARI=1).
- The original male HCSO results were six clusters for most models, five for LC2; female original clusters were two. We get male K=2 throughout and female K=3 for LC2. Differences cannot be attributed to a single factor without matched data/age/category/DTW and K-selection controlled studies.
- A **silhouette score is not a forecasting accuracy measure**. No §4.3 multi-population forecasts or held-out prediction gains are claimed.

## 在 RStudio 中完整运行

1. 下载整个仓库，打开本目录 [run_4_2_complete.R](run_4_2_complete.R)，点击 **Source**。
2. 选择此前的 `fitted_rates.csv`（3888行，含 `deaths_gbd`、`exposure_proxy`）。
3. 程序会重新拟合、聚类、核验和绘图，自动在输入文件旁建立时间戳结果目录。
4. 直接查看 Console 的14组结果和新目录中的 Figures 5–8。不需要安装第三方包。

```r
options(section42.autorun = FALSE)
source("run_4_2_complete.R", encoding = "UTF-8")
run_section42("D:/mortality/fitted_rates.csv", "D:/mortality/new_R_results")
```

或在命令行：

```bash
Rscript run_4_2_complete.R fitted_rates.csv new_R_results
```

输出目录必须不存在或为空。原始GBD输入未公开上传；`fitted_count`不用于估计。已清理多余的旧RStudio启动器，现在只有 `run_4_2_complete.R` 作为本章正式入口。原四个Python脚本已从当前工作树移除，历史可由Git记录追溯。

## R 文件逐项对应

| 文件 | 内容 |
|---|---|
| [§4.1/model_functions.R](../../4_1_single_population/model_functions.R) | 数据检查、七模型目标函数/解析梯度、三起点估计、SVD（属于§4.1） |
| [§4.1/run_all_models.R](../../4_1_single_population/run_all_models.R) | 126组单人口模型拟合、参数导出（属于§4.1） |
| [cluster_all_models.R](cluster_all_models.R) | 14组聚类、旋转/窗口/删项敏感性（属于§4.2） |
| [clustering_functions.R](clustering_functions.R) | DTW动态规划、hclust、精确medoids、七项CVI、ARI |
| [complete_cvi.R](complete_cvi.R) | 196组CVI核对、294次PAM起点检查 |
| [verify_and_plot.R](verify_and_plot.R) | 有限差分梯度、路径枚举、结果完整性、迁移对照 |
| [plot_paper_figures.R](plot_paper_figures.R) | Figures 5–8、热图、树状图、全模型曲线 |
| [Notion R代码精读](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2) | 详细教学已放在Notion，GitHub只保留源码与验证结果 |

## 当前 R 实跑文件

`results/` 中的全部CSV已更新为R计算结果：拟合诊断、alpha/beta/kappa、标准化因子、504个距离、196个候选CVI、主分组、敏感性、294次PAM检查，以及迁移数值对照。环境与设置见 `sessionInfo.txt`、`manifest.txt`；检查见 `validation.csv`。

[完整核验与误差限度](R_MIGRATION_VALIDATION.md)记录：最大缩放梯度3.143e−9、三起点最大拟合面差4.505e−6、14组归属全部ARI=1。优化器和随机数发生器不同，因此数值尾数并非逐位相等；没有把旧结果改名当作R结果。

[Notion解释与结果](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)｜[原论文](https://doi.org/10.1080/03461238.2025.2480846)

