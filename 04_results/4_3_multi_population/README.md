# §4.3｜多人口共享因子模型：R 语言联合估计、验证与作者附录对照

> **状态（2026-11-12）**：七类多人口模型的 R 解析梯度和两类合成数据优化实测已通过 [GitHub Actions](https://github.com/ppkk1029/gbd-hungary-mortality-replication/actions/workflows/section43-r-synthetic-tests.yml)。已修正 Lee–Carter 的识别约束以遵循原文 Table 2，而非只锚定一个年龄点。但**尚未使用真实 GBD 3888 行 fitted_rates.csv 完成联合估计**：输入文件不在公开仓库，本环境中没有原始数据。请勿将测试通过当作图 9、原文结果或预留样本预测已成功复现。

[原论文 PDF](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf)｜[作者附录第 9 目录：52 张 PNG](https://github.com/LiviaVarga/Forecasting-cause-of-death-mortality-with-single-and-multi-population-models-in-Hungary/tree/main/9.%20The%20fitted%20parameters%20of%20the%20multi-population%20models)｜[§4.2 已实跑完整 R 结果](../4_2_clusters/gbd_complete/README.md)

## 研究问题（与 §4.2 的衔接）

§4.1 分死因单独估计死亡率；§4.2 按各模型的时期因子用 DTW 进行分簇；§4.3 要检验：同一簇联合使用各死因的年龄×年份观测，是否能用少量共同时期参数解释趋势。关键操作是 **联合优化似然**，而不是把独立 κ 平均。

例：男性 Poisson LC1 中 NEU/DKD 同组。联合模型分别拟合 NEU/DKD 的两套年龄参数，却联合拟合唯一一条 κ(t)；更新 κ 的似然梯度同时包含两种死因。因此死亡量高的死因可能对共同因子产生更大的影响。

## Table 2：七种共享模型的正确公式

以下 j 是一个 DTW 簇内的死因；t 是年份，x 是年龄；不带 j 的时期因子表示共享。

| 模型 | 线性预测子（简写） | 共享 | 死因特有 |
|---|---|---|---|
| Poisson Carter–Lee | α(x,j)+β1(x,j)κ1(t) | κ1 | α、β1 |
| Binomial Carter–Lee | 同上，用 logit 概率 | κ1 | α、β1 |
| Poisson LC2 | α(x,j)+β1(x,j)κ1(t,j)+β2(x,j)κ2(t) | **κ2** | α、β1、β2、κ1(t,j) |
| Binomial LC2 | 与 Poisson LC2 同线性结构 | **κ2** | 同上 |
| Binomial CBD | κ1(t)+(x−x̄)κ2(t,j) | **κ1** | κ2(t,j) |
| reduced Plat2 | α(x,j)+κ1(t)+(x̄−x)κ2(t,j) | **κ1** | α、κ2(t,j) |
| reduced Plat3 | α(x,j)+κ1(t)+(x̄−x)κ2(t,j)+(x̄−x)_+κ3(t,j) | **κ1** | α、κ2(t,j)、κ3(t,j) |

实现中使用已有 §4.2 R 参数初始化，并严格执行原文 **Table 2 的 Lee–Carter 识别条件：每一种死因 j、每一个因子 i 的年龄系数之和为 1（Σₓβₓⱼ⁽ⁱ⁾=1）；所有相应时期因子按年份求和为 0**。这与此前“只固定一个年龄的 β”不等价，尤其在时期因子跨死因共享时，会改变联合模型允许的拟合范围。双因子 LC 初始化先通过代数重参数化将各死因两条 β 的和均调整为 1，再联合估计共享的 κ²；该重参数化在共享之前保持各死因的单人口拟合面。CBD 不含 α。有限差分梯度检验是强制步骤；优化后的最大缩放梯度仍须检查。年龄基函数沿用 §4.2 的 /10 缩放，因此不能直接比较与作者使用不同标度的 κ 数值。

本项目 GBD 2023 模拟死亡数为连续估计，而非 HCSO 登记整数。对于 Binomial，本项目仍沿用 §4.2 的 N=exposure_proxy+D/2 近似：方法性拟似然对照，不应当写成精确的原文 HCSO 极大似然。

## 输入与模型数量

- 1980–2015 训练、2016–2019 预留测试；年龄 60–89，六个五岁年龄组；男女分别拟合；九种 GBD Level-2 死因。
- 7 模型×2 性别，§4.2 的 CSV 确定共 30 个簇。其中 27 簇至少两个死因，可估计共同参数；另外 3 簇是单独死因，直接保留单人口模型。
- 原文使用九个 ICD-10 章节，与 GBD 原分类**不可一一对应**；结果只能对照机制和图型，不能声称各分类数值相等。

## Figure 9：为什么黑线与彩线分开？

原文 Figure 9 以三因子 reduced Plat 为例，**黑线**是同簇各死因联合估计的共同时期因子 κ1(t)；**彩线**是各死因各自单独估计的 κ1(t,j)。彩线不是多人口 Plat3 模型中估计出的 κ1(t,j)：联合模型的第一个因子恰恰不带 j。

本项目 Figure9_GBD 类比图也按此读取，一簇一个面板。观察趋势方向、拐点、黑线接近哪一类死因；用 common_vs_single_diagnostics.csv 的相关系数、RMSE 配合 cause_fit_metrics.csv 的死亡量占比和训练 MAPE 判断共同因子是否较偏向大死因。**这不是预测优劣验证**，2016–2019 回测要另作。

## 作者附录第 9 目录全部 52 张：明确区分两类

作者共有 **26 个 HCSO 簇×每簇 2 张图=52 张 PNG**。我把原始 GitHub 里 52 个真实文件名，与自动生成的文件名做过逐一交叉核对，0 缺失、0 多余。复核脚本 [author_appendix_inventory.R](author_appendix_inventory.R)，运行后生成含全部 52 条来源链接的 author_original_appendix_52_files.csv。

| 图名开头 | 共26张 | 说明 |
|---|---:|---|
| The fitted parameters ... | 26 | 展示共享模型年龄参数、共同时期参数和死因特异时期参数的估计结果 |
| The common and cause-specific mortality indices ... | 26 | 比较联合估计的共同趋势与独立估计的各死因趋势，判断聚类内部差异是否较大 |

26组场景：Poisson / Binomial Carter–Lee、CBD、Plat2、Plat3 五种模型，每种 3 个男性簇+1个女性簇（共20组）；Poisson/Binomial 双因子 LC 两种模型，每种 2 个男性簇+1个女性簇（共6组）。原文男性 ICD 簇包含 I_XIV、II_XI、IX_XX 等；双因子 LC 使用不同的 II_IX_XI_XX 组；女性包括 IX_XI_XIV_XX、IV_IX_XI_XIV_XX。不能直接挪作 GBD 的死因标签。

本项目因为 §4.2 有 **27 个 GBD 非单例簇**，同样对每组生成两张图，因此预计生成 **54 张新的 GBD 同类型图**；这不是作者原 52 张的像素复刻，而是逐类、逐模型重新绘制的同类型结果。每张文件名包含性别、模型、cluster 编号。本地运行之后由 appendix_image_manifest.csv 列出相对路径。

## 在 RStudio 中执行（完整仓库）

打开 [run_4_3_complete.R](run_4_3_complete.R)，点击 Source；选择与 §4.2 相同的原始 fitted_rates.csv。保持本目录、相邻的§4.1七模型目录和§4.2聚类目录完整，确保 §4.2/results 内有 parameters_all_models.csv、kappa_all_models.csv、cluster_membership.csv。脚本将依次检查数据→加载单人口参数及簇→逐簇联合优化→数值梯度验证→输出拟合误差及参数→绘 Figure9 / GBD 附录图。

若需要自己设路径、逐步运行，先关闭自动入口：

~~~r
options(section43.autorun=FALSE)
source("run_4_3_complete.R",encoding="UTF-8")
run_section43(
  input="D:/mortality/fitted_rates.csv",
  out="D:/mortality/Section_4_3_NEW",
  source42="D:/mortality/gbd-hungary-mortality-replication/04_results/4_2_clusters/gbd_complete/results",
  restarts=1L,
  maxit=3500L,
  strict=TRUE,
  forecast_preview=FALSE
)
~~~

新建输出目录必须不存在或为空；命令行入口为 Rscript run_4_3_complete.R /path/fitted_rates.csv /path/new_results。预计主输出包括 fit_summary.csv、cause_fit_metrics.csv、joint_age_parameters.csv、joint_period_indices.csv、finite_difference_gradient_checks.csv、author_original_appendix_52_files.csv、appendix_image_manifest.csv、common_vs_single_diagnostics.csv、Figure9_GBD...png 及 54 张 GBD 附录对照图。RDS 和运行环境/输入哈希也同时保存。检查 fit_summary.csv 中所有联合组状态为 CONVERGED 后才能讨论拟合图。

Figure10_GBD_Male... 和 FigureA4_GBD_Female... 是可选的 **§4.5 预测预览**，不是 §4.3 的已验证预测；**现在默认关闭**（只有显式 `forecast_preview=TRUE` 才绘制）。其简化 ARIMA/RW 漂移尚未完成原文逐项异常值处理与 ARIMA 诊断，不能声称完成 §4.5。

## 自动验证与实际复现的区别

在仓库根目录运行（不需要私有 GBD 原始文件）：

```bash
Rscript 04_results/4_3_multi_population/selftest_4_3.R
```

脚本利用**明确标为合成数据**的年龄×年份×死因张量，针对 7 种 Table 2 模型逐一核对：模型维度、年龄 β 求和约束、时期 κ 中心化、目标函数与解析梯度的有限差分一致性；另对 Poisson Carter–Lee 和三因子 Plat 执行真正的优化程序，并核对作者附录第 9 目录的 52 个文件名。GitHub Actions 每次修改上述 R 源码都会运行同样的检查。[测试脚本](selftest_4_3.R)｜[CI 工作流](../../.github/workflows/section43-r-synthetic-tests.yml)。

这些检查不能替代真实 GBD 拟合（尤其不能证明 27 组联合模型均会收敛），也不能检验原文 HCSO 曲线一致性。原文使用整数登记死亡数、官方期初/期中暴露量和 Newton–Raphson；本仓库使用 GBD 估计数、代理暴露量和 R 的 L-BFGS-B / nlminb 目标函数优化。二者即使拥有相同模型公式也不保证出现相同参数数值。

## 核验限制与下一步

- 已实测：7 种模型的合成数据梯度，2 种模型的合成数据优化，30 个旧分组中 27 个应联合拟合的分组个数，作者附录 52 张图的文件目录核对，以及 R 语法检查。
- **没有核验：本环境运行 §4.3 联合拟合结果、每组是否优化收敛、54 张新图是否成功生成、样本外是否更准。**
- 原始 HCSO 人口/死亡数未获得；目前是与原文机制相同但数据定义不同的 GBD 方法性复现。
- 模型较高维，初次 RStudio 运行若梯度检查不通过，应先排查估计的尺度、锚点和数值稳定性，不应取消检查把失败结果写成论文实证。
