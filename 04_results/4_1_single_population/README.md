# §4.1｜七种单人口死亡率模型（Figures 3–4）

所有**模型拟合函数和126组纯R重估驱动**现在集中于本章。图3/图4仍有各自独立R脚本和实际GBD图、CSV。图形专用脚本与七模型共同估计器的参数识别方式可能不同，**不能直接拼接两类脚本输出的κ参与DTW**。

## 源代码

| 文件 | 用途 |
| --- | --- |
| [model_functions.R](model_functions.R) | 七模型公式、数据验证、工作似然、解析梯度、SVD规范化和优化 |
| [run_all_models.R](run_all_models.R) | `fit_all_models(input,out)`：只拟合126组单人口模型、导出参数/诊断 |
| [run_4_1_poisson_lc.R](run_4_1_poisson_lc.R) | Figure3单因子Poisson LC；交互选择Population和Deaths |
| [run_4_1_poisson_lc_2factor.R](run_4_1_poisson_lc_2factor.R) | Figure4双因子Poisson LC；从原D/E重新拟合，不用旧fitted_count冒充观测 |
| [diagnose_previous_beta.R](diagnose_previous_beta.R) | 异常β的尺度/抵消比检查（不重新拟合） |
| [plot_figure3_standard.R](plot_figure3_standard.R) | Figure3只改变绘图与参数规范化，不删除女性DKD |

## 结果证据

[Figure3已核验异常报告](FIGURE3_VERIFIED_RESULTS.md) · [女性DKD独立beta图](female_diabetes_beta_standard.svg) · [18组诊断CSV](figure3_constraint_diagnostics_2026-11-10.csv)

[Figure4真实GBD报告（历史独立Python拟合）](FIGURE4_REAL_GBD_RESULTS_2026-11-10.md) · [男性参数图](figure4_male_actual_data.svg) · [女性参数图](figure4_female_actual_data.svg) · [18组汇总CSV](figure4_actual_summary_2026-11-10.csv)

最新完整126组**纯R**运行证据：[R迁移验收](../4_2_clusters/gbd_complete/R_MIGRATION_VALIDATION.md)；七模型参数及诊断CSV仍保存在[§4.2运行产物](../4_2_clusters/gbd_complete/results/)（因为该次运行同时执行了§4.1与§4.2）。

## 核心结果及限制

单因子Poisson LC训练deviance总计约4268.354，双因子约875.941；女性DKD约568.834降至18.101。Plat2→Plat3训练deviance也降低。**上述改善仅说明训练期拟合，更复杂的模型未必提高2016–2019样本外预测。** GBD死亡人数不是登记整数计数，模型应称为工作性伪似然。

**推荐完整入口：** [§4.2/run_4_2_complete.R](../4_2_clusters/gbd_complete/run_4_2_complete.R)，会自动加载本目录两个模块重新估计并执行聚类。长篇教学笔记见[Notion｜§4.1](https://app.notion.com/p/413eb49c94ed811f8f3ad354001239f7)。
