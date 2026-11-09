# 匈牙利分死因死亡率预测：GBD 2023 方法复现（R）

本项目面向初学者，按 Varga (2025) 论文的研究顺序，用 GBD 数据复现分死因死亡率预测方法。**只制作原文实际涉及的分析与图表**；每章解释研究问题、计算公式、R 代码、实际结果以及与论文的异同。尚未运行的分析明确标记为“待复现”。

**原文**：Lívia Varga, *Forecasting cause-of-death mortality with single- and multi-population models in Hungary*, *Scandinavian Actuarial Journal*. [DOI](https://doi.org/10.1080/03461238.2025.2480846)

**使用数据**：IHME GBD 2023，Hungary，1980–2023 年，男性和女性，60–89 岁六个五岁年龄组，19 个在该年龄范围有数据的 GBD 死因类别；Deaths 的 Number 与 Rate。属于**方法复现**，不是原始 HCSO 数据的数值复现。

## 按论文顺序学习

| 章节 | 对应论文 | 当前状态 |
|---|---|---|
| [00 数据核对](00_data/README.md) | §3 数据准备 | 已有核对代码 |
| [01 死因构成](01_figure1/README.md) | §3.1 Figure 1 | 已有计算代码；缺 All causes 分母，暂不生成正式占比图 |
| [02 年龄标准化死亡率](02_figure2/README.md) | §3.1 Figure 2 | 已有 ASMR 计算和九面板图代码；待本地原始 CSV 执行与结果验证 |
| [03 单人口死亡率模型](03_single_population/README.md) | 单人口模型与拟合 | 待复现 |
| [04 DTW 聚类](04_dtw_clustering/README.md) | 时期因子聚类 | 待复现 |
| [05 多人口死亡率模型](05_multi_population/README.md) | 聚类后的联合建模 | 待复现 |
| [06 模型比较及预测](06_forecast_comparison/README.md) | 拟合、回测、预测结果 | 待复现 |

## 如何运行前两节

1. 从 [GBD Results Tool](https://vizhub.healthdata.org/gbd-results/) 下载原始 CSV，将其保存在本地 `data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv`；仓库不会公开上传该原始文件。
2. 在 RStudio 中将工作目录设为仓库根目录。
3. 逐段运行 `00_data/check_data.R` 核对数据。
4. 运行 `01_figure1/run.R` 计算逐年各死因死亡人数；如另行下载同口径 All causes 数据，保存为 `data/raw/gbd_all_causes.csv`，才会生成 Figure 1 的百分比堆积面积图。
5. 运行 `02_figure2/run.R`，程序会按 Number/Rate 反推并核验 1980 年人口代理，使用统一年龄权重计算各死因 ASMR，并在 `outputs/02_figure2` 输出图和 CSV；**当前代码尚未使用原始 CSV 运行，因此不宣称已获得实际结果**。
6. 之后各章仅在真实计算完成后加入程序、结果及图表。

## 必须区分的三个量

- **Number**：某类死亡人数的 GBD 估计。
- **Rate**：每 10 万人口死亡率；要换成每人年近似中心死亡率，除以 100,000。
- **死因占全部死亡比例**：指定死因死亡人数 / 同年龄、年份范围的 **All causes** 死亡人数。19 类的内部占比不能冒充全部死亡比例。

## 与原文的主要差异

| 方面 | Varga 原文 | 本项目 |
|---|---|---|
| 来源 | 匈牙利 HCSO 登记及协调数据 | GBD 估计值 |
| 历史年份 | 1970–2021 | 1980–2023 |
| 死因 | 9 类 ICD 章节 | GBD 原生分类；最终筛选规则要预先明确 |
| 年龄 | 单岁年龄，90+ 合并 | 60–89 岁六个五岁组，不含 90+ |
| 人口暴露量 | central / initial exposure | 当前原始文件仅包含死亡人数与率 |
| 模型回测 | 2016–2019 | 训练期截止 2015，测试期同样选 2016–2019 |

### 方法与可复现性说明

GBD 估计死亡人数不完全等同于原始登记计数。对基于 Poisson/Binomial 的模型，必须明确分布假设的局限；不同 GBD 轮次的历史重估也可能影响回顾性预测。选死因、选模型及 DTW 聚类不得利用测试期结果。

本仓库没有公开原始 GBD CSV；使用者应从官方数据源获取并遵守署名和再分发规则。**代码尚未在本次上传中通过 RStudio 实际运行测试**。已有数值核查来自上传 CSV 的独立检查，后续以本地 R 运行结果为准。
