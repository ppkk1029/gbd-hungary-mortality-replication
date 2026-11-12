# §4.2｜DTW 死因聚类

本章节只保留一套正在使用的 R 代码：

**[打开最新版 run_4_2_complete.R](gbd_complete/run_4_2_complete.R)** · [完整 RStudio 运行说明](gbd_complete/README.md)

原先只按 Silhouette 选 K 的聚类结果 CSV、SVG 和迁移记录已删除，避免误用； §4.1 拟合的时期因子已归到 [§4.1/results](../4_1_single_population/results/README.md)。新版使用 R 的 `dtw`、`dtwclust`；两种聚类方法、三个 PAM 起点与七项 CVI 全部可输出。最终 K 需要查看树状图与全部指标后明确记录选择理由，不能写死作者的 K。

本仓库仍为 **GBD 方法复现**，并非使用原文 HCSO 数据的数值复现。
