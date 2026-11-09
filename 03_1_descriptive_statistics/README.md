# §3.1｜Descriptive statistics（描述性统计）

**按原论文的 3.1 小节整理为一个 Part**，包含两个不同的描述性统计内容：
- **Figure 1**：九类死因人数占全部死亡人数的构成（堆积面积图）；
- **Figure 2**：各死因 ASMR（年龄标准化死亡率）的历史趋势，男女分别计算。

这里暂时保留以前的代码和图形文件以避免已有 GitHub 链接失效，不新增“Figure 1 Part”和“Figure 2 Part”：
- [Figure 1 原有 R 与详细文档](../01_figure1/README.md)；
- [Figure 2 原有 R 与详细文档](../02_figure2/README.md)。

**本项目当前数据**：GBD Hungary，1980–2023，60–89 岁六组，九类按 1980–2015 年人数选出；原文为 HCSO、1970–2021、60+、九个 ICD 章节。这是方法性复现，不是数值复刻。

Figure 1 的“死因人数占比”与 Figure 2 的“年龄标准化死亡率”不是同一个统计量；§4.1 模型拟合应使用每年人口与 Deaths Number，不应使用 Figure 2 固定的 1980 年标准权重。

后续以 [§4.1 + §4.2 研究结果的分节目录](../04_results/README.md) 为主。
