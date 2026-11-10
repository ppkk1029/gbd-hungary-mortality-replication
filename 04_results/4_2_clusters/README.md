# §4.2｜死因聚类：完整 R 版本

**当前版本：七种模型、全部聚类、Figures 5–8 和核验均用 R 执行。**

- [完整代码与结果说明](gbd_complete/README.md)
- [RStudio 一键入口](gbd_complete/run_4_2_complete.R)：打开后点 Source，选择含 `deaths_gbd` 与 `exposure_proxy` 的 `fitted_rates.csv`。
- [逐步代码讲解：问题→方法→为什么→结果→原文对照](gbd_complete/CODE_WALKTHROUGH_R.md)
- [126组拟合及14组聚类的R实跑核验](gbd_complete/R_MIGRATION_VALIDATION.md)
- [Notion完整报告](https://app.notion.com/p/413eb49c94ed81149a9bfdca526c0bc2)

本目录旧文件名 `run_4_2_dtw.R` 也转向同一个完整R流程，避免误运行仅支持单因子的早期原型。旧实现和早期结果可通过Git历史追溯，当前工作树以此版本为准。

数据仍为GBD Hungary、1980–2015、60–89岁、九类死因。双因子LC取κ2，其余模型取κ1；逐条样本z标准化，平方局部差的symmetric2 DTW，average-linkage与精确枚举medoids，K=2…8，最大平均Silhouette选择。七项CVI与旋转、五年窗口、删除女性DKD的敏感性均保留。

这是方法复现，不是HCSO数值复刻。分组属于后续共享模型的候选结构，§4.3及样本外预测仍未完成。
