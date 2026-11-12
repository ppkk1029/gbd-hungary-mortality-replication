# §4.2｜DTW 死因时间序列聚类：唯一新版 R 实现

**代码已更新。原有 Silhouette-only 聚类文件、旧 Figure5–8 和旧迁移文档已删除。** 本目录只保留正在使用的 R 代码与方法说明。

## 正文对应

1. 用 §4.1 七种单人口模型的时期因子，男性/女性分别进行比较。Poisson/Binomial LC2 对第二时期因子聚类，其他模型对第一因子聚类。
2. 每条序列单独 Z 标准化。按公式 (6)–(8) 用 `dtw` 计算 p=2、`symmetric2`、`sqrt(累积平方成本 / (n+m))` 距离，并用 `dtwclust::dtw_basic` 的原始成本逐项验证。
3. 同一 DTW 距离用于 `dtwclust` 的层次聚类（平均联接）和 PAM，候选组数 K=2:8，PAM 记录三个随机起点。
4. 用 `dtwclust::cvi` 计算 **Sil、D、COP、DB、DBstar、CH、SF** 七项原生指标，输出全部 196 个候选方案、树状图、距离矩阵和三起点稳定性。
5. **论文没有公布七指标如何合并为最终 K 的固定公式。** 代码不自动根据轮廓系数、投票或自造评分选择 K；必须结合原文正文与在线附录，填写 14 组审查记录，才能输出正式最终分组。
6. 输出 Figure 5–8 的 **GBD 类比图**、女性配套图、其他模型的标准化组内曲线、七项指标对 K 的变化图及组内对齐路径表。图 5–6 区分未标准化 κ（展示）与标准化 κ（聚类）；女性 DKD 不能冒充作者的 ICD-IV。

## 最简单运行方式：读取已保存的 §4.1 时期因子

先在 RStudio 安装软件包（只需要一次）：

```r
install.packages(c("dtw","dtwclust"))
```

保持项目文件夹结构。在 RStudio 中打开 `run_4_2_complete.R`，然后运行：

```r
options(section42.autorun = FALSE)
source("run_4_2_complete.R", encoding="UTF-8")
out <- run_section42()
```

`run_section42()` 会自动读取：

`04_results/4_1_single_population/results/kappa_all_models.csv`

并在本目录的 `outputs/Section42_日期时间/` 下生成独立结果。**不会重跑 126 个单人口拟合，不会读取任何旧的 §4.2 聚类结果。**

要从原始 GBD 死亡数和暴露量全部重估，可显式传入：

```r
out <- run_section42("D:/data/fitted_rates.csv")
```

## 最终选 K：必须有理由，不能把原文 K 直接写死

第一阶段输出 `cluster_selection_template.csv`（14 行）、`cluster_validity_scores.csv`（196 行）、`candidate_membership.csv`、`pam_three_start_check.csv`、`metric_preferences.csv` 和每种模型的七指标图与树状图。先查看这些数据；将模板另存为 `reviewed_selection.csv`，逐行填写：

- `K`：2–8 的整数；
- `method`：`hierarchical_average` 或 `pam`；
- `reason`：对七项 CVI、树状图、单成员组和 PAM 稳定性的具体判断依据，不能为空。

之后运行：

```r
finalize_section42(out, file.path(out,"reviewed_selection.csv"))
```

此时才生成正式 `cluster_membership.csv`、`cluster_summary.csv`、`cluster_selection_decisions.csv` 和相应图表。后续 §4.3 要传入 **这个新输出目录**，不能使用历史聚类结果。

## 代码清单

| 文件 | 功能 |
| --- | --- |
| [run_4_2_complete.R](run_4_2_complete.R) | 唯一入口；从 §4.1 因子或原始 GBD 输入运行、两阶段审查 |
| [kappa_import.R](kappa_import.R) | 严格验证并导入 126 组已拟合时期因子 |
| [clustering_functions.R](clustering_functions.R) | 双包 DTW 距离一致性、层次/PAM、七项 CVI |
| [cluster_all_models.R](cluster_all_models.R) | 14 组模型×性别，全体候选方案与人工选 K |
| [complete_cvi.R](complete_cvi.R) | 按指标输出最优 K、完整性审查 |
| [plot_paper_figures.R](plot_paper_figures.R) | Figure5–8 类比、在线附录类型图表与对齐 |
| [verify_and_plot.R](verify_and_plot.R) | R 数值自检和计算结果完整性核验 |
| [selftest_section42.R](selftest_section42.R) | 自动化合成数据测试 |

## 仍不可声称与作者完全相同的部分

- 作者是 HCSO 1970–2021、单岁年龄及 ICD 九个主类；这里用 GBD 1980–2015、五岁年龄组和 GBD 二级死因。原始输入不等价，最优组数不应被强制匹配原文。
- 作者未完整披露层次联接选择、DTW 窗口细节、七指标最终汇总规则与全部随机种子；本版将 `average`、无窗口等明确列为可复现的额外假定。
- 现有 §4.1 的双因子 LC 使用 SVD 识别规范，可能改变第二因子的聚类；不能通过 §4.2 代码凭空恢复作者未公开的模型约束。
- 经审查确认的最终 (K) 在真实数据计算并完成理由记录之前，不对表 3、表 4 的成员归属做任何真实性宣称。

论文：https://doi.org/10.1080/03461238.2025.2480846
