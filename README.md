# 匈牙利分死因死亡率预测：GBD 2023 方法复现（R）

本项目面向初学者，按 Varga (2025) 论文的研究顺序，用 GBD 数据复现分死因死亡率预测方法。**只制作原文实际涉及的分析与图表**；每章解释研究问题、计算公式、R 代码、实际结果以及与论文的异同。尚未运行的分析明确标记为“待复现”。

**原文**：Lívia Varga, *Forecasting cause-of-death mortality with single- and multi-population models in Hungary*, *Scandinavian Actuarial Journal*. [DOI](https://doi.org/10.1080/03461238.2025.2480846)

**使用数据**：IHME GBD 2023，Hungary，1980–2023 年，男性和女性，60–89 岁六个五岁年龄组，19 个在该年龄范围有数据的 GBD 死因类别；Deaths 的 Number 与 Rate。属于**方法复现**，不是原始 HCSO 数据的数值复现。

## 按原文小节组织 Part（不是一图一个 Part）

| 论文小节 | 研究内容和材料 | 当前状态 |
|---|---|---|
| [§3.1 描述性统计](03_1_descriptive_statistics/README.md) | Figure 1 死因构成 + Figure 2 年龄标准化死亡率，**同属一小节** | 用户已在 RStudio 绘制 GBD Figure 1、Figure 2；原文件与质控记录尚未全部上传 |
| [§4.1 单人口模型拟合](04_results/4_1_single_population/README.md) | Figure 3、Figure 4 及七种模型的参数和诊断 | 用户 Figure 3 截图显示女性糖尿病与肾脏疾病的 β 极端值；已新增数值诊断 + 稳健重估脚本，尚需运行验证 |
| [§4.2 死因聚类](04_results/4_2_clusters/README.md) | DTW 公式（6）–（8）、对齐路径、层次/PAM、簇数选择、Figures 5–8、Tables 3–4 | 已实现**单因子 Poisson LC** 的完整聚类流程与异常剔除敏感性；七模型全部复现仍缺其他六种因子 |
| [§4.3 多人口模型](05_multi_population/README.md) | 类内共享时期因子建模 | 待复现 |
| [§4.4–4.7 模型比较与预测](06_forecast_comparison/README.md) | 拟合、ARIMA、回测及长期预测 | 待复现 |

**数据检查**：[00 数据核对](00_data/README.md)。早期按 Figure 分开的 `01_figure1`、`02_figure2`、`03_single_population` 保留用于旧链接兼容；**以新的小节索引和脚本为准**。

## 2026-11-10 新增：Figure 3 实测诊断与异常 beta 分图

- [§4.1 **我们的实际 Figure 3 分析报告**：18 组数值、异常分图理由、模型与数据边界](04_results/4_1_single_population/FIGURE3_VERIFIED_RESULTS.md)
- [女性糖尿病／肾脏疾病 beta 独立图（真实结果）](04_results/4_1_single_population/female_diabetes_beta_standard.svg) 与 [18 组 CSV](04_results/4_1_single_population/figure3_constraint_diagnostics_2026-11-10.csv)；[Figure 3 六面板一键绘图代码](04_results/4_1_single_population/plot_figure3_standard.R)
- 该异常组正负抵消比 **0.00334884**，按标准约束归一化时 `max(abs(beta))=98.849445`；其余 17 组最大不超过 **0.329914**。**九类死因不从估计或分析中剔除**，仅将异常的女性 beta 另绘以免遮蔽其他曲线。
- Figure 4 的二因子 Poisson LC 必须重新拟合；其他五种单人口模型连同一因子 Binomial LC 共六种尚待完成。单因子 Poisson LC **1/7**；§4.2 DTW 算法代码已就绪、真实聚类与七模型比较**尚未运行核验**。

## Figure 4｜双因子 Poisson LC（代码已就绪；真实结果待用户本地运行）

- [Figure 4 R 拟合与导出程序](04_results/4_1_single_population/run_4_1_poisson_lc_2factor.R)：读取此前 Figure 3 输出的 `fitted_rates.csv`，重新拟合实际 GBD Deaths/Exposure，18组模型，输出男女各五条参数类别图和 deviance 比较。
- [详细执行说明、数学识别约束与证据范围](04_results/4_1_single_population/FIGURE4_REPLICATION_PLAN_AND_VERIFICATION.md)。
- [独立模拟自测图](04_results/4_1_single_population/figure4_synthetic_numeric_selftest.svg) 明确标注 **synthetic, not actual GBD**。当前会话无法读入 Library 的 CSV 原始字节；**不能把程序准备完毕或模拟测试误写成 Figure 4 真实拟合已完成**。

## 继续复现的正确顺序（全部交互选文件）

1. 若要检查你现有截图的女性 β≈18，先运行 [§4.1 旧参数诊断](04_results/4_1_single_population/diagnose_previous_beta.R)，在 `file.choose()` 窗口选择**旧版 beta.csv、kappa.csv**；输出异常抵消比和尺度校验，**不会覆盖旧结果**。
2. 再运行 [§4.1 修正的 Poisson LC](04_results/4_1_single_population/run_4_1_poisson_lc.R)，依次交互选择 Population、Deaths，先在 Console/Plots 检查重估参数和异常标记，确认 `y` 后再选择结果保存父文件夹；生成新的 `kappa.csv`、`fit_diagnostics.csv`。
3. 最后运行 [§4.2 DTW](04_results/4_2_clusters/run_4_2_dtw.R)，依次选择**新版** `kappa.csv` 和 `fit_diagnostics.csv`；展示男女 DTW 对齐、树状图、PAM 和聚类表后，输入 `y` 再选择保存文件夹。
4. §4.2 不混用不同模型的时间因子；原文双因子 Poisson/Binomial LC 后续应聚类**第二个**因子。当前仅实现单因子 Poisson LC 的聚类方法，不能声称七种模型均已复现。

## 必须区分的三个量

- **Number**：某类死亡人数的 GBD 估计。
- **Rate**：每 10 万人口死亡率；要换成每人年近似中心死亡率，除以 100,000。
- **死因占全部死亡比例**：指定死因死亡人数 / 同年龄、年份范围的 **All causes** 死亡人数。本项目 Level 1 三大类和 19 个 Level 2 死因已在 528 个年×性别×年龄单元核验一致。Figure 1 使用 Level 1 合计作为全部死亡人数分母，2019 年约 99,182.03 人（60–89 岁男女合计）。

## 与原文的主要差异

| 方面 | Varga 原文 | 本项目 |
|---|---|---|
| 来源 | 匈牙利 HCSO 登记及协调数据 | GBD 估计值 |
| 历史年份 | 1970–2021 | 1980–2023 |
| 死因 | 9 类 ICD 章节 | Level 2，按 1980–2015 训练期死亡人数前九类，**不重组 ICD** |
| 年龄 | 单岁年龄，90+ 合并 | 60–89 岁六个五岁组，不含 90+ |
| 人口暴露量 | 1970 年年中人口用于标准化；模型还涉及各年 exposure | 1980–2023 GBD Population 已收到并核验；可用于标准化权重及年度一致性检查，但不等于 HCSO central exposure |
| 模型回测 | 2016–2019 | 训练期截止 2015，测试期同样选 2016–2019 |

### 方法与可复现性说明

GBD 估计死亡人数不完全等同于原始登记计数。对基于 Poisson/Binomial 的模型，必须明确分布假设的局限；不同 GBD 轮次的历史重估也可能影响回顾性预测。选死因、选模型及 DTW 聚类不得利用测试期结果。

本仓库没有公开原始 GBD CSV；使用者应从官方数据源获取并遵守署名和再分发规则。**§3.1 已由用户在 RStudio 显示图形；§4.1 已由用户显示单因子 Poisson LC 图，但女性一组参数异常需要重估。**新版本 §4.1、§4.2 代码尚未在完整 Deaths CSV 上端到端运行，后续以本地诊断、数值一致性和图表为准。
