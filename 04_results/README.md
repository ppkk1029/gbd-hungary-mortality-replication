# 第 4 节｜Results and discussion（按原论文小节组织）

**重要规则：按论文的小节组织一个独立 Part，不再把 Figure 3、Figure 4、Figure 5、Figure 6、Tables 3–4 分割成不同 Part。** 每个小节可包含多张图、多个表和多个模型；保留以前的脚本作为历史版本，但新代码以本目录为准。

- [§4.1 The fitted parameters of the single-population models](4_1_single_population/README.md)：单人口模型参数。当前实现了 **单因子 Poisson Lee–Carter 的修正版**、旧版异常 β 诊断与 Figure 3 风格参数图；其余六个模型尚待实现。
- [§4.2 The clusters of causes of death](4_2_clusters/README.md)：DTW 标准化、距离、最优对齐路径、层次聚类、PAM、簇数选择、Tables 3–4 风格分组汇总，全部属于本小节。当前计算适用于已经实现的 **单因子 Poisson LC**；其余模型暂未得到可信时期因子，不虚构其分组。

## 已知数据和结果边界

GBD 2023：Hungary，1980–2023，60–89 岁六个五岁组，两性分开，19 个 Level 2 死因中按训练期 1980–2015 死亡人数选九类。

Population CSV 已检查：528 行且真实有效；本轮附件尚未含 Deaths 原始 CSV 或 Figure 3 的 `beta.csv`/`kappa.csv`/`fit_diagnostics.csv`，故**不能确认女性异常的具体数值原因，也不能给出真实 §4.2 分组结果**。用户提供的 Figure 3 截图是 RStudio 的实际绘图结果，但异常模型需要复核。

## 运行顺序

先运行 §4.1 的 `diagnose_previous_beta.R` 来检测旧参数是否发生 β 极端缩放（只需旧 `beta.csv` 与 `kappa.csv`），然后运行新的 `run_4_1_poisson_lc.R` 用 Population 与 Deaths 重新拟合、查看拟合诊断并保存。

当且仅当 §4.1 成功产生**新版** `kappa.csv` 与 `fit_diagnostics.csv` 后，运行 `4_2_clusters/run_4_2_dtw.R`。所有脚本都先在 RStudio 预览图表、确认后通过对话框选择输出文件夹，无固定输入路径。

论文：[Varga (2025)](https://doi.org/10.1080/03461238.2025.2480846)。
