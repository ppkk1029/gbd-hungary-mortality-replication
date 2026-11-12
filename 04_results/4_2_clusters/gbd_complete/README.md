# §4.2｜按 Varga (2025) 方法重写的 DTW 聚类（GBD 方法复现）

> **2026-11-12 修订**。此目录的正式 R 程序已改为使用 `dtw` 与 `dtwclust`，不再使用“Silhouette 最大自动选择最终 K”的旧规则。**GBD 与作者的 HCSO/ICD 数据不同，不能承诺复现相同 K 或 Figure 5–8 数值。**
>
> **结果状态**：仓库内现有 `results/*.csv` 与 `figures/*` 是 **2026-11-10 旧版规则的历史运行产物**，并非此版代码的输出。旧结果保留以免影响已有 §4.3 数据依赖；不要与本次新方法混用。新运行会写到独立的时间戳目录。

## 与论文 §4.2 的对应关系

| 环节 | 论文可确认的做法 | 本版处理 | 边界 |
|---|---|---|---|
| 单人口模型 | 七种模型，分别取时期因子 | 复用 §4.1 的全部七模型估计，LC2 用第2因子，其余用第1因子 | GBD 原始数据口径不同 |
| DTW | 标准化，公式(6)–(8)，p=2 | `dtw::dtw()` 在平方局部成本下用 symmetric2，并除以序列长度和开平方；`dtwclust::dtw_basic(normalize=FALSE,sqrt.dist=FALSE)` 的未标准化累计成本逐对核对后进行**同一外部归一化** | 原文没有充分披露所有窗口与参数细节；无窗口是显式假定 |
| 层次聚类 | Hierarchical + dendrogram | `dtwclust::tsclust(type="hierarchical")`，平均联接；通过 proxy 注册与距离矩阵一致的 `paperDTW_p2` | 原文没有披露联接方式，average 是可核验的额外假定 |
| PAM | Partition around medoids, 三个起点 | `dtwclust::tsclust(type="partitional",centroid="pam")`，使用相同的 `paperDTW_p2` 距离，逐 K 运行三组可复现种子 | 并不声称作者使用了相同随机种子 |
| 聚类质量 | Silhouette、Dunn、COP、DB、DBstar、CH、SF | **七项全部由 `dtwclust::cvi()` 计算**，另保留算法×K×成员归属及单成员组数量 | 部分指标在某种分组下可能不可计算；记录 NA 和原因 |
| 最终 K | 树状图、两种算法、7 CVI、重复起点综合判断 | 默认 **researcher review**：输出全部候选和填写表，经人工审核后才能生成正式结果 | **原文没有公布七项指标的唯一合并公式，严禁虚构权重、投票法或硬编码 K=5/6** |
| 画图 | Fig.5–8 + 在线附录 | 原始 κ 展示与标准化 DTW 路径分开；树状图/热图/聚类组图 | 女性 GBD 的 DKD 不等于原文 ICD IV，明确写在 Figure6 类比图上 |

## 如何在 RStudio 使用（两阶段；不重复拟合）

首先安装依赖（一次即可）：

```r
install.packages(c("dtw","dtwclust"))
```

**第一阶段：生成全部证据**

保持仓库的目录结构。在 RStudio 打开 `run_4_2_complete.R`。首次运行可以直接点 **Source** 选择本地原始 `fitted_rates.csv`；或者指定明确输入路径：

```r
options(section42.autorun=FALSE)
source("run_4_2_complete.R",encoding="UTF-8")
out <- run_section42("D:/mortality/fitted_rates.csv",
                     "D:/mortality/section42_revised")
```

程序重新拟合七个单人口模型，输出独立目录的候选诊断：

- `cluster_validity_scores.csv`：**196 行**（14 组×2 算法×7 个候选 K）及 **七项 dtwclust 原生 CVI**。
- `candidate_membership.csv`：全部 196 个候选方案的成员标签（1764 行）。
- `pam_three_start_check.csv`：三起点的目标函数和与最优起点的 ARI（294 行）。
- `metric_preferences.csv`：每项 CVI 各自支持的 K（并不意味着多数投票）。
- `dendrograms_Male/Female.png`、`distance_heatmaps_Male/Female.png`、两性 DTW 对齐图。
- `PREVIEW_NOT_FINAL_*.png`：仅供观察的候选分组，**不是最终 Figure 7/8**。
- `cluster_selection_template.csv`：需要审阅后填写 `K`、`method`、`reason`，共 **14 行**。

**第二阶段：结合完整评价和树状图，填写最终 K**

检查原文 Table 3/4 仅作对照，不要强行复制作者的 K。把 `cluster_selection_template.csv` 另存为 `reviewed_selection.csv`，逐行填写：

- `K`：整数 2–8；
- `method`：`hierarchical_average` 或 `pam`；
- `reason`：写下为什么支持这个 K（建议注明哪些 CVI 支持、是否有大量单成员组）。

然后不需要再拟合 126 组模型：

```r
finalize_section42("D:/mortality/section42_revised",
                   "D:/mortality/reviewed_selection.csv")
```

此时才输出正式的 `cluster_summary.csv`、`cluster_membership.csv`（§4.3 可读取的结构）、`cluster_selection_decisions.csv` 以及 Figure 7/8 的 **GBD 类比图**。

命令行用法：

```bash
Rscript run_4_2_complete.R fitted_rates.csv section42_revised
Rscript run_4_2_complete.R --finalize section42_revised reviewed_selection.csv
```

*可选探索模式*：`run_section42(input,out,selection_mode="exploratory_consensus")` 会按**七 CVI 在各算法内的平均名次**生成明确标记为探索性的分组。这是本仓库的附加规则，**不是论文原文规则，不能把它当作确认后的最终 K**。

## 代码分工

- [run_4_2_complete.R](run_4_2_complete.R)：唯一入口、两阶段选 K、无需重复拟合的最终确认函数。
- [clustering_functions.R](clustering_functions.R)：`dtw`、`dtwclust`、PAM 三起点、七项原生 CVI、候选排名。
- [cluster_all_models.R](cluster_all_models.R)：14 组实验、全部候选、独立的审查/最终决策输出。
- [complete_cvi.R](complete_cvi.R)：候选数和起点结果核验、每个指标的最优 K 对照表。
- [verify_and_plot.R](verify_and_plot.R)：梯度、DTW 路径及输出完整性检查。
- [plot_paper_figures.R](plot_paper_figures.R)：对齐图、树状图、热图、组内曲线。
- [selftest_section42.R](selftest_section42.R)：合成数据的包接口与计算一致性测试（不需要私有 GBD 输入）。
- [历史验收记录](R_MIGRATION_VALIDATION.md)：仅证明上一版纯 R 迁移成功，**不证明本次新代码已重跑完整 GBD**。

## 尚不能通过 §4.2 代码消除的差异

1. 数据：原文 HCSO 1970–2021、单岁年龄与 ICD 死因；当前是 GBD 1980–2015 训练段、六个五岁年龄组、GBD 二级原因。原始 HCSO 数据未在本仓库提供。
2. LC2 的因子识别：§4.1 当前使用 SVD 正交规范化，旋转可不改变拟合面却改变第二因子聚类。需要额外对照作者的完整估计和约束，不能在 §4.2 强行“修好”。
3. 原文没有提供唯一的七指标合成和最终选 K 规则，也没有给出所有代码/随机种子。所以本仓库**只对可验证的方法步骤给出确定性实现，最终 K 要有可追踪理由**。

参考：Varga (2025), *Scandinavian Actuarial Journal*, https://doi.org/10.1080/03461238.2025.2480846；dtwclust CVI 文档 https://asardaes.r-universe.dev/dtwclust/doc/manual.html 。
