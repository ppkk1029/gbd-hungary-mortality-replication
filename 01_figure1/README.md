# 01｜Figure 1：主要死因构成（原文 §3.1）

## 原文做了什么？
原文 Figure 1 是按全部死亡人数计算占比的堆积面积图：九类 ICD 死因加 `All others`。它**不是**男女各七条折线图；本项目只做原文这一类图。

## 按 GBD 官方层级理解分母
GBD Protocol（2024）规定：死因分类在同一层级互斥且穷尽；分死因死亡数应汇总成全死因死亡数。因此“完整 Level 2 的总和”在定义上就是 All causes。参见 [官方 Protocol](https://www.healthdata.org/sites/default/files/2024-06/GBD%20Protocol%20060424.pdf) §2、§7。

但当前下载的 `IHME-GBD_2023_DATA-b40ffc46-1.csv` 仅有 **19 个 Level 2 类别的分年龄数据**（44 年 × 2 性别 × 6 年龄组 × 19 类 × Number/Rate）。之前的宽表另外出现 `Maternal and neonatal disorders`、`Mental disorders` 的 All ages 汇总记录，没有对应的六个年龄组；`Sense organ diseases` 没有出现在这份死亡人数表。因此**不能不经核对就认定 19 类合计与 All causes 精确一致**。

## 核心计算
先计算 19 类内部的年度死亡人数：
`D19[t] = sum_c∈19 sum_x sum_s D[x,t,s,c]`

再取训练期前九类作为**GBD 方法演示**，得到：
`p19[t,c] = D[t,c] / D19[t] × 100`。

该比例以 19 类合计为分母，**不是已核验的全死因占比**。正式全死因构成用：
`p_all[t,c] = D[t,c] / D_all[t] × 100`。

## 代码、输出、核验
[run.R](run.R) 的前半段用 `subset()` 提取 Number，用 `aggregate()` 按年和死因汇总。
- 没有单独的 All causes 文件时，代码生成 `figure1_level2_internal_unverified.png`，图题和输出表明确写明分母是 19 类，仅用于初步对照。
- 若有 `data/raw/gbd_all_causes.csv`，代码先逐年对比 `D19` 与 `D_all`，再用官方全死因数计算正式占比，残余用 `All other causes` 表示。
- 不把父类与子类、汇总年龄与五岁年龄组重复相加。

## 现有数据可核对的结果
2019 年六个五岁年龄组、男女合计的 19 类死亡人数合计为 **99,182.03**（GBD 点估计）。这不应直接被标注为“官方 All causes”。如新增官方分母，先比较差值后再报告全部死因覆盖率。

## 与原文对照的限制
原文：1970–2021、HCSO、九类 ICD、包含 90+。本项目：1980–2023、GBD、19 类原生 Level 2、60–89 岁。九类候选的选取原则也与原文不同，不能直接要求图形数值一致。
