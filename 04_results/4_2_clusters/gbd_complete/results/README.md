# 历史结果快照（旧版§4.2）

这些 CSV 是 2026-11-10 在**自写 CVI + 仅按最大 Silhouette 选 K** 的旧版算法下实跑的数据。

**不得称为 2026-11-12 的 dtw / dtwclust 七指标审查版结果。** 保留它们是为了可追溯和防止旧 §4.3 默认输入瞬间失效。

新版运行写入一个独立的输出目录，其中 `cluster_selection_template.csv`、`cluster_validity_scores.csv` 由官方包计算；只有人工确认 K 并调用 `finalize_section42()` 后，才生成新的 `cluster_membership.csv` 和 `cluster_summary.csv`。不要把此目录的旧 CSV 与新版计算输出混合。

参见 [新版使用说明](../README.md)。
