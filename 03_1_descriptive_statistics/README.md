# §3.1｜描述性统计（Figures 1–2）

| 任务 | 正式R文件 | 验证状态 |
| --- | --- | --- |
| Figure 1：死因死亡人数构成（堆积面积图） | [figure1/run.R](figure1/run.R) · [说明](figure1/README.md) | Level1与Level2人数一致性、历史数值已核验；RStudio完整重跑日志尚未入库 |
| Figure 2：年龄标准化分死因死亡率（ASMR） | [figure2/run.R](figure2/run.R) · [说明](figure2/README.md) · [1980人口权重](figure2/verified_age_weights_1980.csv) | 权重已独立核验；尚缺完整R重跑日志及1672行正式输出的仓库验收 |

Figure 1统计占全部死亡人数之比；Figure 2使用固定1980年年龄构成的加权年龄别死亡率。Figure 2的ASMR已汇总年龄，**不能直接用于§4.1拟合年龄载荷**。数据为GBD Hungary（1980–2023、60–89岁六组）而非HCSO原文数值。

**运行目录要求：** Figure1脚本使用相对仓库根路径 `data/raw/`，请从仓库根目录运行；Figure2脚本通过 `file.choose()` 选择Population和Deaths文件。详细原理已移至[Notion Figure1](https://app.notion.com/p/413eb49c94ed81c688bde90150342930)与[Notion Figure2](https://app.notion.com/p/413eb49c94ed81c7a043dd878c2dbaee)。

[数据核验](../00_data/README.md) · [§4.1建模](../04_results/4_1_single_population/README.md)
