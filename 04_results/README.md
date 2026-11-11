# 第4节｜Results and discussion：一节一目录

| 论文小节 | 内容 | 状态 |
| --- | --- | --- |
| [§4.1](4_1_single_population/README.md) | 126组单人口模型、Figures3–4 | 纯R已完成拟合与数值验收 |
| [§4.2](4_2_clusters/README.md) | DTW、聚类、Figures5–8及CVI | 纯R已实跑14组并保存结果 |
| [§4.3](4_3_multi_population/README.md) | 多人口共享因子联合似然、Figure9及附录图 | R代码完成；尚无完整实跑产物 |
| [§4.4–4.7](4_4_to_4_7_forecast/README.md) | 单/多人口模型预测及测试比较 | 尚未完成 |

依赖关系：§4.1 `model_functions.R` 提供七模型和输入验证，`run_all_models.R` 负责126组拟合；§4.2入口自动加载§4.1函数及本节聚类模块；§4.3从§4.1加载模型数据接口，从§4.2加载聚类与参数CSV。移动仓库时应保留相对目录结构。

本仓库缺原始GBD CSV，保存的验证结果不能取代复跑；数值收敛也不意味着预测正确。
