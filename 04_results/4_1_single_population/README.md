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

## 四、还缺什么？

- 原文 §4.1 涉及七种模型：Poisson LC、Binomial LC、双因子 Poisson LC、双因子 Binomial LC、CBD、2因子/3因子 reduced Plat。目前只实现 Poisson LC（**1/7**）。
- Figure 4 与另外模型的参数尚未估计；第 4.2 节当前只能对已有的 Poisson LC κ 做真实聚类，不能制造其他六套 κ。
- 本轮没有可用的 Deaths CSV 供我执行完整拟合；上面的数值质量检查会在用户本地真实运行后才产生结果。
- 原文对照须保留 GBD Level 2 / 年龄 60–89 / 1980–2015，与 HCSO ICD 九章 / 60+ / 1970–2021 的区别。

原文第 12–13 页 Figure 3、Figure 4：[PDF](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf#page=13)。
