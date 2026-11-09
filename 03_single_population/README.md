# 03｜单人口死亡率模型：§4.1 / Figure 3

**本章已新增程序：[run_poisson_lc.R](run_poisson_lc.R)**。先复现 Poisson Lee–Carter 单因子模型和 Figure 3 的参数图；其他六种单人口模型以及双因子 Figure 4 尚未实现。

**状态**：代码已提交，已经做公式及计算结构的静态审查；没有在用户的完整 Population + Deaths 文件上实际运行，故目前没有真实的参数或 Figure 3 结果。用户运行成功并核验诊断表后才能标记本章第一步完成。

## 一、与原文的关系

Varga（2025）在 §4.1 分别对九个 ICD 死因、男女两性拟合七种单人口模型，Figure 3 展示 Poisson LC 的年龄基准参数、年龄敏感度和时期因子；Figure 4 展示男性的双因子 Poisson LC。随后 §4.2 才对时期因子进行 DTW 聚类。

**本项目选择**：GBD 2023 Hungary，1980—2023；六个五岁年龄组（60–64 … 85–89）；两性分别拟合；19 类 Level 2 死因中，沿用 Figure 1/2 的规则，以 1980—2015 年死亡总数选出九类。不进行 ICD 九章重组，不追求原文的相同参数数值。

本章第一阶段拟合 1980—2015；2016—2019 留作回测；2020—2023 暂时排除疫情影响。原文 Figure 3 的拟合年份更长（1970—2021），所以两图的历史区间不同。

## 二、模型与参数

对一种死因、一个性别，拟合

$$
D_{x,t}\sim \operatorname{Poisson}(E_{x,t}m_{x,t}),\qquad
\log m_{x,t}=\alpha_x+\beta_x\kappa_t.
$$

- D：GBD 的 Deaths Number；
- E：同一年龄、年份、性别的 GBD Population，作为暴露量代理；
- α：该死因在不同年龄的基础对数死亡率（年龄效应）；
- β：各年龄对时间变化的响应系数；
- κ：逐年的时期因子，后续 DTW 使用这条时间序列。

模型的两个识别约束：

$$
\sum_x\beta_x=1,\qquad \sum_t\kappa_t=0.
$$

例子：β 为 0.2，κ 从 5 降到 4，则该年龄的拟合死亡率乘以 exp(−0.2)，即约下降 18.1%。如果 β 为负，就不能按这个方向解释。

GBD Number 是模型估计，可能不是整数，因此**本程序是借助 Poisson 均值与 deviance 的拟似然/准似然式拟合，不是原文使用整数登记死亡数的严格 Poisson 计数似然**。程序使用 BFGS + 解析梯度优化，作者使用 Newton–Raphson；统计均值结构相同，但优化算法不完全相同。不能把这里的结果称作对作者参数的逐值复刻。

## 三、数据

从个人电脑通过两个窗口依次选择 Population CSV 和 Deaths（Number + Rate）CSV。两文件不限路径、不限文件名：

1. Population：44 年×2 性别×6 年龄组=528 行；人口点估计应为正。
2. Deaths：44 年×2 性别×6 年龄组×19 死因×2 指标=20064 行，含 Number 和 Rate。
3. 通过年份、年龄和性别严格匹配人口暴露量；按照 Number/Population×100000 与下载的年龄别 Rate 进行一致性核验。
4. LC 拟合使用 Number 与 Population，**不使用 Figure 2 的固定 1980 年标准人口权重**。固定标准权重只用于描述性 ASMR。

## 四、如何运行（先展示，再保存）

在 RStudio 打开 [run_poisson_lc.R](run_poisson_lc.R)，点击 Source 运行整份脚本。也可以在 Console 运行下列命令，然后第一次选择该 R 脚本：

    source(file.choose())

之后先选择 Population CSV，再选择 Deaths CSV。不要只运行代码中间的一小段。

程序检查两文件后，对**九个死因×两种性别=18 个模型**进行参数估计。计算完成时：

- Console 显示九类死因、质量检查、18 个模型的收敛代码和拟合诊断。
- RStudio 数据查看器展示参数和诊断表。
- **Plots 窗格先展示六面板图**：男性 α、β、κ；女性 α、β、κ；底部统一图例。

仅当你在 Console 输入 y 后，才弹出保存文件夹的选择窗口。选择父文件夹后，在那里新建 Figure3_Poisson_LC_results 子文件夹；已存在时会再次确认是否覆盖。任何输入数据路径与保存位置都不预先强制指定。

## 五、选择保存之后的结果文件

| 文件名 | 内容 |
|---|---|
| alpha.csv | 2 性别×9 死因×6 年龄的 α |
| beta.csv | 2 性别×9 死因×6 年龄的 β |
| kappa.csv | 2 性别×9 死因×36 年的 κ，后续 DTW 输入 |
| fit_diagnostics.csv | 18 个模型的收敛状态、Poisson deviance、梯度指标 |
| fitted_rates.csv | 逐年龄、年、性别、死因的死亡率与拟合率 |
| selected_causes.csv | 按训练期死亡人数确定的九类 |
| figure3_poisson_lc.png | 经过 RStudio 预览后确认保存的 Figure 3 风格图 |
| run_summary.txt | 本次运行日期、年份范围、数据核验和方法限制 |

这些文件目前**尚未由本项目实测生成**；具体估计结果待用户 RStudio 运行后核对。

## 六、先检查什么，再解释结果？

- 每个模型 sum(β)=1、sum(κ)=0；拟合率为正，拟合 deviance 为有限值。
- BFGS 收敛代码均为 0；若任何模型失败，程序停止，不要直接使用其 κ 做 DTW。
- β 可以有正有负；κ 可以下降或上升，不强制死因趋势一致。
- 六个年龄组远少于原文细龄数据，可能增加估计不稳定性；后续还应做多起点或敏感性检查。
- **参数图只说明历史拟合，不说明预测准确率**；预测须在 2016—2019 年回测。

下一讲继续补原文 Figure 4 的双因子 Poisson LC，之后完整实现另外的模型并开展 §4.2 的 DTW 聚类。

论文：[Varga（2025）§4.1 / Figure 3](https://unipub.lib.uni-corvinus.hu/11046/1/Forecastingcause-of-deathmortalitywithsingle-andmulti-populationmodelsinHungary.pdf)。
