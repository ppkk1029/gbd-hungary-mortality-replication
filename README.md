# 匈牙利分死因死亡率预测｜GBD 2023 方法复现（R）

> **2026-11-12 更新**：已按原论文章节清理目录；§4.3 的七类 R 联合模型通过自动合成数据测试，但还没有真实 GBD 联合拟合图。以 R 代码、实际结果及验收记录为准；完整学习笔记放 Notion。使用 GBD 而非 HCSO，只能称为**方法复现**，不是原文数值重现。

## 研究进度与目录

| 章节 | 研究任务 | 当前状态 |
| --- | --- | --- |
| [00 数据核验](00_data/README.md) | GBD Level 1/2 人数一致性 | 1,056个年龄×性别×年份×指标单元核验通过 |
| [§3.1 描述统计](03_1_descriptive_statistics/README.md) | Figure 1死因构成、Figure 2年龄标准化率 | Figure 1有实际计算结果；Figure 2标准人口权重已核验，完整R重跑产物尚未入库 |
| [§4.1 七模型拟合](04_results/4_1_single_population/README.md) | Figures 3–4及七种单人口模型 | **126组**纯R模型拟合，参数和梯度检查有验收记录 |
| [§4.2 DTW聚类](04_results/4_2_clusters/README.md) | Figures 5–8、Tables 3–4、7项原生CVI与聚类数审查 | **新版代码已改为 dtw/dtwclust，需重新运行真实GBD数据并人工审查最终K**；旧结果是Silhouette-only历史结果 |
| [§4.3 多人口模型](04_results/4_3_multi_population/README.md) | Table 2 类内共享因子联合拟合、Figure 9与附录 | **七类模型的数值梯度及合成数据测试通过；真实联合拟合结果和图像仍待生成** |
| [§4.4–4.7 预测](04_results/4_4_to_4_7_forecast/README.md) | 预测模型选择、2016–2019测试、长期预测 | **未完成正式样本外比较** |

## 唯一推荐入口

- §3.1：[Figure1/run.R](03_1_descriptive_statistics/figure1/run.R)（从仓库根目录运行）、[Figure2/run.R](03_1_descriptive_statistics/figure2/run.R)（交互选原始数据）。
- §4.1：[model_functions.R](04_results/4_1_single_population/model_functions.R)（七模型数学与优化）+ [run_all_models.R](04_results/4_1_single_population/run_all_models.R)（126组拟合）；[Figure 3](04_results/4_1_single_population/run_4_1_poisson_lc.R)与[Figure 4](04_results/4_1_single_population/run_4_1_poisson_lc_2factor.R)另有专用绘图/拟合脚本。
- §4.2：[run_4_2_complete.R](04_results/4_2_clusters/gbd_complete/run_4_2_complete.R)：先安装 `dtw`、`dtwclust`；Source后选择本地 `fitted_rates.csv`，生成14组全部聚类候选、7项CVI及树状图。**须填写审查后的K和方法并调用 `finalize_section42()`，才输出正式分组与Figure 7–8。** 详见[使用说明](04_results/4_2_clusters/gbd_complete/README.md)。
- §4.3：[run_4_3_complete.R](04_results/4_3_multi_population/run_4_3_complete.R)：依赖本地 `fitted_rates.csv` 与§4.2的参数和分组CSV；[selftest_4_3.R](04_results/4_3_multi_population/selftest_4_3.R) 可在不依赖私有数据的情况下用 R 验证七种模型；[GitHub Actions](https://github.com/ppkk1029/gbd-hungary-mortality-replication/actions/workflows/section43-r-synthetic-tests.yml) 自动运行，真实 GBD 联合估计仍待验收。

## 数据与限制

GBD 2023 Hungary、1980–2023、两性、60–89岁六个年龄组，19种Level2死因中按1980–2015训练期死亡人数选择9种；2016–2019留作测试，2020–2023未进入当前训练。GBD死亡人数是连续估计，不等于HCSO整数登记计数；Poisson/Binomial采用工作性伪似然。原始Population/Deaths及完整`fitted_rates.csv`未公开提交仓库，须自行获得同口径输入；查看[数据说明](data/raw/README.md)。**特别注意：§4.2 的 checked-in `results/` 和 `figures/` 均是旧的 Silhouette-only 方法产物，不能视为新版七项 CVI 方法的结果。**

Notion：[Figure 1](https://app.notion.com/p/413eb49c94ed81c688bde90150342930) · [Figure 2](https://app.notion.com/p/413eb49c94ed81c7a043dd878c2dbaee) · [§4.1代码精读](https://app.notion.com/p/413eb49c94ed811f8f3ad354001239f7) · [§4.2代码精读](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)。

旧占位README、无效旧入口、重复聚类表和重复学习文档已清理；正式实验图及数值文件保留。清理前版本见[备份分支](https://github.com/ppkk1029/gbd-hungary-mortality-replication/tree/archive/pre-chapter-cleanup-2026-11-11)及Git提交历史。
