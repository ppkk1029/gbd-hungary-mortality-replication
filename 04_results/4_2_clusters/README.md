# §4.2｜DTW 死因聚类：按原文审查式选 K 的 R 修订版

> **本目录正式代码已替换（2026-11-12）。** 新版依赖 CRAN 包 `dtw`、`dtwclust`；使用七项原生 CVI、层次聚类与 PAM 三起点。作者未提供唯一的七项指标合并公式，因此默认生成全部证据，再由研究者确认最终 K；不再使用旧版 Silhouette 最大值直接决定 K。**仓库旧结果 CSV/图片是旧方法历史产物，不是新版复跑结果。**

| 文件 | 用途 |
|---|---|
| **[正式入口 run_4_2_complete.R](gbd_complete/run_4_2_complete.R)** | 从原始 GBD 输入估计七个单人口模型，生成所有候选；再通过 `finalize_section42()` 复用计算结果确认 K |
| [clustering_functions.R](gbd_complete/clustering_functions.R) | dtw 距离、dtwclust 分组、七项原生 CVI、三次 PAM 起点 |
| [cluster_all_models.R](gbd_complete/cluster_all_models.R) | 两性七模型全部候选分组与明确的选 K 决策 |
| [complete_cvi.R](gbd_complete/complete_cvi.R) | 评价指标支持的 K 与审计 |
| [plot_paper_figures.R](gbd_complete/plot_paper_figures.R) | 原始 κ 与标准化距离分开、对齐图、树状图和组内曲线 |
| [selftest_section42.R](gbd_complete/selftest_section42.R) | 合成数据核验；由 GitHub Actions 自动运行 |

[完整使用步骤、方法边界和结果状态](gbd_complete/README.md) · [历史 R 迁移验收](gbd_complete/R_MIGRATION_VALIDATION.md)

**特别提醒：** 原文 HCSO 1970–2021 的单岁年龄和 ICD 死因分类，与当前 GBD 1980–2015、六个五岁年龄组及 GBD 死因类别不同。新方法的 K 不必与原文的男性/女性五组或六组相等。重新审查后的分组可能影响 §4.3 的模型拟合任务数，不能继续按旧版固定的 30 组解释。
