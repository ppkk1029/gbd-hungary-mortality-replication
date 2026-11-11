# §4.2｜死因DTW聚类（Figures 5–8，Tables 3–4）

**只保留本章聚类与诊断代码。** 七模型估计模块已归入[§4.1](../4_1_single_population/README.md)；§4.2入口仍自动从原始D/E重新拟合全部七模型，再运行聚类，避免混用时期因子版本。

| 部分 | 入口 |
| --- | --- |
| **唯一正式入口** | [run_4_2_complete.R](gbd_complete/run_4_2_complete.R) |
| 选取κ、所有模型聚类与敏感性 | [cluster_all_models.R](gbd_complete/cluster_all_models.R) |
| 标准化、DTW、层次聚类、精确medoids | [clustering_functions.R](gbd_complete/clustering_functions.R) |
| CVI与随机起点核验 | [complete_cvi.R](gbd_complete/complete_cvi.R) |
| 梯度/路径及结果核验 | [verify_and_plot.R](gbd_complete/verify_and_plot.R) |
| Figure 5–8与附加诊断图 | [plot_paper_figures.R](gbd_complete/plot_paper_figures.R) |

[完整结果、Figure5–8、14组聚类表](gbd_complete/README.md) · [纯R实跑验收](gbd_complete/R_MIGRATION_VALIDATION.md) · [Notion R代码逐步学习笔记](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)

目前14组模型×性别聚类、504个DTW距离、196种候选CVI、294项多起点比较已保存；仅是GBD的历史趋势聚类，不证明§4.3共享模型可提高预测。
