# §4.1｜The fitted parameters of the single-population models

**章节单位：§4.1（一个小节一个 Part）。** 原文 Figure 3（Poisson LC 参数）和 Figure 4（男性双因子 Poisson LC 参数）都属于 §4.1；后续所有单人口模型参数图、诊断与对照也放本目录，不再一张图单独建章节。

## 一、先诊断你截图的女性 β 大于 15

截图中“Diabetes and kidney diseases（糖尿病与肾脏疾病）”女性 β 约从 −11 增至 +18，配套 κ 接近 0。**这是数值异常信号，但不等于原始 Deaths CSV 已证实有错**。原版程序只检查 `sum(beta)=1` 和优化收敛；当正负 β 大量抵消时，可能产生绝对值很大的 β 与振幅很小的 κ，同时保持 `βκ` 不变：

例：`15×0.02=0.3×1`。只有直接查看旧 `beta.csv`、`kappa.csv`、拟合率与数据，才能判断是识别尺度或非稳定最优。

**先运行** [diagnose_previous_beta.R](diagnose_previous_beta.R)：在 RStudio 交互选择旧版 `beta.csv`、`kappa.csv`；查看各死因、性别的 `max_abs_beta_old`、`cancellation_ratio=abs(sum(beta))/sum(abs(beta))`、`sd_kappa_old`、`identical_fit_max_abs_error`。如比值小于 0.08，标记 REVIEW。检测不会覆盖文件，不会改动死亡率拟合；在 R 的 Console 和 View 中查看即可。

## 二、修正后的拟合

[run_4_1_poisson_lc.R](run_4_1_poisson_lc.R) ：RStudio 打开后 Source，依次用 `file.choose()` 选择 Population CSV、Deaths CSV；无需固定文件名。

原模型：`D_xt ≈ E_xt exp(alpha_x + beta_x*kappa_t)`。GBD Deaths Number 是模型输出的连续估计，不是真实整数登记计数，因此本脚本使用 Poisson 均值结构的拟似然式拟合。采用 `1980–2015` 训练，`2016–2019` 留作以后回测，`2020–2023` 不参与拟合。女性/男性分别为9类拟合，共18个独立模型。

### 关键更改：使用**稳定的、但与原文不同的识别约束**

原来 `sum(beta)=1`、`sum(kappa)=0`；某些年龄载荷正负抵消会使第一条约束病态。

新版在优化时让 **SVD 最大绝对载荷对应的一个参考年龄 β=1**，并令 `sum(kappa)=0`；优化后转换为：

```
max(abs(beta)) = 1
mean(kappa) = 0
```

仅改变参数的代表尺度，**数学上不改变 βκ 或同一参数解的拟合死亡率**；脚本逐单元格检查重新归一化前后的 `alpha+beta*kappa` 是否一致。但由于数值优化也重新执行，新版结果不保证与之前拟合在每个单元格完全相同，因此仍需查看偏差与收敛质量。

新版用 3 个起点、BFGS 解析梯度优化，记录收敛状态、梯度、deviance 和 `beta_cancellation_ratio`。异常 `beta_flag=REVIEW_near_zero_beta_sum` 不是论文效果结论，而是后续 §4.2 的风险提示。若有异常，仍可**探索性**分析，并分别做剔除异常的敏感性分析；绝不能直接据此宣布跨死因共享趋势可靠。

## 三、结果在 RStudio 如何出现？

先在 Console 展示参数和诊断，Plots 生成一个六面板参数图（男性与女性各 α、β、κ），View 打开参数和诊断表；确认保存输入 `y` 后，再选择输出父文件夹，自动创建 `Figure3_Poisson_LC_results`（这是**结果文件夹名**，不是单独的论文小节或 Part）。

输出 `alpha.csv`、`beta.csv`、`kappa.csv`、`fit_diagnostics.csv`、`fitted_rates.csv`、`selected_causes.csv`、`figure3_poisson_lc.png`、`run_summary.txt`。

## 三点五、2026-11-10 实际 Figure 3 展示结果（已核验 CSV）

- [**结果分析：18 组模型诊断、图 3 单独分图理由与 Figure 4/DTW 进度**](FIGURE3_VERIFIED_RESULTS.md)
- [**真实数据女性糖尿病／肾脏疾病 beta 独立图**](female_diabetes_beta_standard.svg)；[18 组诊断 CSV](figure3_constraint_diagnostics_2026-11-10.csv)
- [**主图与独立图一键重绘 R 脚本**](plot_figure3_standard.R)：所有九类在 alpha、kappa 均保留；女性 beta 主面板暂不叠加异常曲线，改由独立图展示。

本轮恢复 `sum(beta)=1` 后：女性糖尿病与肾脏疾病抵消比 **0.00334884**；`max(abs(beta))` 为 **98.849445**；其余17组的最大值不超过 **0.329914**。这支持“共轴会压扁其他曲线”的绘图决策，不证明模型拟合失败；新拟合与预测比较仍待运行。

## 三点六、Figure 4 双因子模型的下一步实证流程（2026-11-10）

- [RStudio 实际重新拟合脚本：run_4_1_poisson_lc_2factor.R](run_4_1_poisson_lc_2factor.R)。
- [完整方法、原始输入文件、输出清单、识别约束、数值自测与证据状态](FIGURE4_REPLICATION_PLAN_AND_VERIFICATION.md)。
- 选择**已生成的 fitted_rates.csv**，脚本真正使用其中的 `deaths_gbd` 和 `exposure_proxy` 重新拟合 18 个**双因子**模型；`fitted_count` 仅用于与已有单因子模型比较，不用于重估。
- 输出男性 Figure 4 参数图和女性对应参数图，以及训练期 deviance 比较图、18 组诊断、第二时期因子 CSV（供后续 DTW）。
- 本轮独立 Python 数值自测使用模拟 Poisson 数据（**不是匈牙利真实拟合**），梯度最大误差约 `1.21e-12`，标准化重建误差 `8.88e-16`。[仅模拟的测试图](figure4_synthetic_numeric_selftest.svg)。此处验证方法，不等于已实证完成 Figure 4。
- 当前运行环境无法从用户 Library 取得 **fitted_rates.csv** 原始字节，且没有 R 运行时，故 **Figure 4 的实际 18 组模型尚未计算；不能将自测数字写成预测/实证提升**。

## 2026-11-10｜Figure 4 **真实 GBD 数据已拟合并取得图片**

- [**本项目真实 Figure 4 实验报告（含18组定量对照和真实图）**](FIGURE4_REAL_GBD_RESULTS_2026-11-10.md)
- [男性真实 Figure 4 五面板](figure4_male_actual_data.svg) | [女性真实 Figure 4 五面板](figure4_female_actual_data.svg) | [18 组实测 CSV](figure4_actual_summary_2026-11-10.csv)
- 用户实际上传的 `fitted_rates(1).csv` 共3888行。用其 `deaths_gbd`、`exposure_proxy` 而非 `fitted_count`，独立 Python/SciPy 重新拟合了9死因×2性别共18个双因子模型。**真实数据结果：18/18 training deviance 下降；伪 AIC 9/18 更低；伪 BIC 6/18 更低。**
- 女性 Diabetes and kidney diseases：一因子 deviance 568.8344 → 双因子 18.1007；训练改善550.7337，伪 BIC 差 -346.47。尚不能推断预测更优。
- [**更新后的 RStudio R 程序**](run_4_1_poisson_lc_2factor.R) 已修复用户截图的 `outdir` 未定义错误：不再询问 `y`，自动在输入文件所在目录新建时间戳结果文件夹并保存图表；**男性5参数Figure4作为最后一张图**，不再停在对比柱状图。
- 这里的实际参数由 Python/SciPy 拟合，RStudio 的优化结果需再次运行并交叉核验；下面原先“未完成 Figure 4”的说明是更早阶段记录，以此处及实验报告为准。

## 四、还缺什么？

- 原文 §4.1 涉及七种模型：Poisson LC、Binomial LC、双因子 Poisson LC、双因子 Binomial LC、CBD、2因子/3因子 reduced Plat。目前 Figure 3 的一因子 Poisson LC 已有实测结果（**1/7**），Figure 4 双因子 Poisson LC 已有可运行的拟合程序但尚未在真实 CSV 上完成实证。
- Figure 4 与另外模型的参数尚未估计；第 4.2 节当前只能对已有的 Poisson LC κ 做真实聚类，不能制造其他六套 κ。
- 本轮已有并核验 **先前拟合生成的 beta、kappa 与 18 组约束诊断**，但仍没有原始 Deaths CSV 供重新运行 Figure 4 或六个其他模型；新的重拟合及 §4.2 聚类尚未执行。
- 原文对照须保留 GBD Level 2 / 年龄 60–89 / 1980–2015，与 HCSO ICD 九章 / 60+ / 1970–2021 的区别。

原文第 12–13 页 Figure 3、Figure 4：[PDF](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf#page=13)。
