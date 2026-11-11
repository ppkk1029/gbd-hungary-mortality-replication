# 05｜多人口死亡率模型

§4.3 的纯 R 联合估计代码、Figure 9 的 GBD 类比绘图及作者附录第 9 目录 52 张源图全量索引，现整理在：

**[§4.3｜multi-population full R implementation](../04_results/4_3_multi_population/README.md)**

本目录保留为项目章节入口，避免破坏 §4.1/4.2 既有路径。

## 本节究竟做什么？

1. 使用 §4.2 的聚类结果，针对每个性别、模型与非单例簇，联合估计一个共同时期因子。
2. 保留 Table 2 不同模型所需的死因特异参数：Carter–Lee 共享第 1 因子；两因子 LC 共享第 2 因子；CBD、Plat2 和 Plat3 共享第 1 因子。
3. 生成 Figure 9 类型的「联合共享 vs 各死因独立时期因子」比较图，并为每个 GBD 非单例簇绘制参数图和对照图。
4. 对照作者附录 9 的 52 张原图，按模型、性别、原始 ICD 簇分配解释每张原图。原图为 26 组×2，GBD 的分组与原文不同，预计为 27 组×2=54 张。
5. 使用拟似然梯度、按死因训练 MAPE、死亡量占比等核验；预测样本外优劣仍待 §4.5。

## 运行状态（2026-11-11）

**R 源码与原图目录已写好，但本环境没有完整的 fitted_rates.csv 输入与 R 可执行文件。** 目前 §4.3 的 27 个联合组并未在这里重跑，也没有发布实际的 GBD Figure 9 或 54 张新图。不应把源代码提交写成实验复现完成，更不应把作者 PNG 误称为本次 GBD 新图。

[§4.3 RStudio 一键入口](../04_results/4_3_multi_population/run_4_3_complete.R)｜[作者全量图索引](../04_results/4_3_multi_population/author_appendix_inventory.R)｜[§4.2 R 实跑结果](../04_results/4_2_clusters/gbd_complete/README.md)｜[Notion 学习笔记](https://app.notion.com/p/415eb49c94ed81e0bd70c9093f5380a0)
