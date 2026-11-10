# §4.2｜The clusters of causes of death

**当前完整实验版本（2026-11-10）：** [点击查看七模型完整工作目录（方法、可运行 Python、Figure 5–7 单图、Table 3/4 对照与 CSV）](gbd_complete/README.md)。已经从以前真实结果包核验：7 模型 × 2 性别 × 9 死因 = 126 拟合；14 次聚类；九类死因全部保留；Figures 5、6、7、8 已分别输出并同步 Notion。此处不是原 HCSO 数据逐值复刻，且 §4.3 预测尚未执行。

**本目录下的** [run_4_2_dtw.R](run_4_2_dtw.R) **是较早的单因子 Poisson LC RStudio 原型，只覆盖 1/7 模型**，不能用于重现当前 14 组聚类；不要从其早期 K=2..6 等参数推断最新实验已经遵守相同规则。完整七模型 Python 版本使用 K=2..8、SVD 识别和精确 k-medoids，结果和所有改动见上述链接。

## 当前 Figure 与结果索引

- [Figure 5｜男性 Poisson LC2 第二时期因子 CVD–DIG DTW 对齐](gbd_complete/figures/Figure5_GBD_Male_DTW.svg)：距离 1.2283，36 对中的第 35；与原文“近似配对”不一致。
- [Figure 6｜女性 Poisson LC2 第二时期因子 CVD–DKD DTW 对齐](gbd_complete/figures/Figure6_GBD_Female_DTW.svg)：距离 0.2044，36 对中的第 1；但 DKD 与原文 ICD IV 非同类。
- [Figure 7｜男性 CBD 第一时期因子的全部九类分组](gbd_complete/figures/Figure7_GBD_Male_CBD.svg)：2 簇，Silhouette 0.714；原文男性 CBD 为 6 簇。\n- [Figure 8｜男性双因子 Binomial LC 第二时期因子的全部九类分组](gbd_complete/figures/Figure8_GBD_Male_Binomial_LC2.svg)：2 簇，Silhouette 0.515；原文男性 Binomial LC2 为 5 簇。
- [14 组实测 K 与 Silhouette](gbd_complete/results/cluster_summary_selected.csv)｜[126 条实测死因归属](gbd_complete/results/cluster_membership.csv)｜[验证数字](gbd_complete/results/validation.json)｜[完整 Notion 复现报告](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)。

---

## 以下保留早期 RStudio 原型说明（历史文档，非当前项目进度）

**本目录包含 §4.2 的所有内容，不按单独 Figure 建 Part。** 包括原文公式（6）–（8）、Figures 5–8、Tables 3–4 所对应的 DTW 最优对齐、距离矩阵、聚类与可视化。

## 1. 原文方法的正确顺序

1. **先拟合单人口模型**，取每种死因的时期因子 κ，男性和女性分别处理。
2. **每条 κ 在训练期内进行 z 标准化**（减均值、除标准差），比较形状而非原始尺度。
3. **在同一个模型、同一个因子编号内**，计算不同死因的两两 DTW 距离。原文式（6）单变量局部代价 `abs(kappa_i[h]-kappa_j[v])`；式（7）使用 `l_p` 范数，文中采用 `p=2`；式（8）沿最优路径累计并归一化。
4. 两种聚类法：层次聚类和 `k-medoids (PAM)`；比较不同簇数和有效性指标；对男女分别给出死因分组。
5. 原文大多数模型使用 **第一时期因子**；**双因子 Poisson/Binomial LC 则使用第二时期因子**。不能混用模型或因子。
6. 分组仅为后续建立多人口模型提供候选共享结构；不证明相同医学病因，也不证明会提高预测精度。

## 2. 本次可运行的版本

[run_4_2_dtw.R](run_4_2_dtw.R) 当前可对**修正后的单因子 Poisson LC（1/7 模型）**完整执行上述分析。原文所有七种模型的 Tables 3–4 尚不能直接复现，因为另外六种模型的 §4.1 参数仍未拟合。

输入：在 RStudio Source 代码，先 `file.choose()` 选新版 §4.1 的 `kappa.csv`，再选同一次运行的 `fit_diagnostics.csv`。程序要求 9 类 × 2 性别 × 36 年，每模型/性别均收敛，并检查时期因子可标准化。

所用 R 包：`dtw`、`cluster`。若未安装，先运行：

```r
install.packages(c("dtw","cluster"))
```

**具体且可复核的软件选项（不是声称完全等于作者的所有默认设置）**：
- 两两 DTW：局部欧氏距离（单维即绝对差）、`step.pattern=dtw::symmetric2`、使用 `normalizedDistance`（按 N+M 归一化）。
- 敏感性分析：使用 `sakoechiba` 窗口限制最大错位为5年，再计算一套 DTW。
- 层次聚类：`hclust(..., method="average")`；PAM：`cluster::pam`。
- 簇数范围 K=2..6（九类死因）；分别给出平均 Silhouette 和 Dunn 指标；**以平均 Silhouette 最大选择首选方案**，平局时取较少簇。这与作者综合多指标/聚类树/三次起点的人工决策并不完全等价。
- 对 `beta_flag!=OK` 的死因标记并输出另外一份“剔除异常因子后的敏感性聚类”（如果剩余至少四类）。完整聚类结果包含异常原因，因此**仅用于探索**，不应直接作为最终共享因子的依据。

## 3. RStudio 先预览、确认后保存

程序先在 Console 输出：
- 18 条 κ 的标准化核查；
- 男女各自选择的最佳聚类方法、簇数；
- 每种死因的所属簇；
- 距离最小的时期因子对；
- 异常 β 对应的死因及敏感性分析。

Plots 历史依次展示：男女树状图与 Silhouette 对比；男性最佳 DTW 对齐及路径；女性最佳 DTW 对齐及路径；男女性别内各簇的标准化时期因子曲线。这些属于**同一个 §4.2 Part**。

确认结果后输入 `y`，再通过文件夹选择窗口决定保存位置。父文件夹里创建 `Section_4_2_DTW_clusters`，包含：

| 文件 | 内容 |
|---|---|
| standardized_kappa.csv | 男女九条标准化 κ（36 年） |
| dtw_pairwise_distances.csv | 全窗口和五年窗口的两两距离 |
| cluster_validity_scores.csv | K=2..6 的层次/PAM Silhouette、Dunn |
| cluster_membership.csv | 各死因所属簇及异常标记 |
| flagged_beta_factors.csv | 不稳定 β 诊断记录 |
| sensitivity_exclude_flagged.csv | 排除异常 β 对应死因的聚类 |
| closest_factor_pairs.csv | 男女各自选择的最接近时期因子对 |
| clustering_overview.png | 树状图、簇数评价 |
| alignment_male.png / alignment_female.png | 最优对齐路径示例 |
| standardized_clusters_male.png / standardized_clusters_female.png | 按簇展示标准化因子 |
| method_and_scope.txt | 完整的参数设置和适用范围说明 |

## 4. 正式科研解释的边界

**不能直接套用原文 Table 3/4 的组别**：原文是九个 ICD 主类、1970–2021、60+；本项目使用 GBD Level 2 九类、1980–2015、60–89，因此分组可能完全不同，且不以相同为目标。

仅用 Poisson LC 的单因素聚类，目前完成的是 **§4.2 的一条模型分支的方法性复现**，而不是原论文七条模型分支的全部实证复刻。先检查女性异常的来源再决定最终分组，后续再补双因子 LC 的第二时期因子聚类和其他模型。

参考：Varga (2025) §4.2，公式（6）–（8）、Figures 5–8、Tables 3–4：[PDF](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf#page=14)；Giorgino (2009) `dtw`：[CRAN vignette](https://cran.r-project.org/web/packages/dtw/vignettes/dtw.pdf)。
