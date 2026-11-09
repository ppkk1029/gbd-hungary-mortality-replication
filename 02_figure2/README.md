# 02｜Figure 2：年龄标准化死亡率（ASMR）｜方法复现

**2026-11-09 核查修订：原版用 Number/Rate 反推人口，未使用独立 Population 文件；现已修订为需要独立人口文件。**

## 一、原文与本项目对照

| 项目 | Varga (2025) 原文 | GBD 方法复现 |
|---|---|---|
| 来源 | 匈牙利 HCSO 数据 | GBD 2023 估计 |
| 年份 | 1970–2021 | 1980–2023 |
| 年龄 | 60 岁及以上（含合并的 90+） | 60–89 岁六个五岁年龄组 |
| 死因 | 固定九类 ICD 章节 I/II/IV/VI/IX/X/XI/XIV/XX | 19 类 GBD 死因；九图选训练期人数前九，不能声称 ICD 类别一致 |
| 标准人口 | 1970 年年中人口 | 独立的 1980 年 GBD Population 文件 |
| 男女权重 | 原文未明确两性标准人口是否合并 | 本项目使用男女合并的固定标准年龄权重 |

因此，这是公式和绘图流程的**方法复现**，不是原文的数值复现。

原文链接：https://doi.org/10.1080/03461238.2025.2480846

## 二、必须准备的两份数据

**1. 死亡数据**（已有）：
路径：data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv

GBD 2023，Hungary，1980–2023，Male/Female，60–89 岁六个五岁组；Deaths 的 Number 与年龄别 Rate。若本地文件名变更，请修改 run.R 中的 deaths_path。

**2. 人口数据**（新增，仍未提供）：
从 https://vizhub.healthdata.org/gbd-results/ 选择 GBD 2023 Population、Hungary、**1980 年**、Male/Female、60–64/65–69/70–74/75–79/80–84/85–89 years、**Number**。保存为：
data/raw/gbd_population_1980.csv

本次标准化**只需要基准年 1980 年的六年龄组×两性=12 条人口记录**，不需要另外下载 44 年全部人口用于标准年龄权重。请优先使用相同 GBD 2023 轮次。若下载工具没有 1980 年选项，不可填造人口；应改用官方可获得的基准年份，并统一修改代码和结果说明。

Population CSV 必须有列：year、sex_name、age_name、val；若还有 location_name、measure_name、metric_name，代码会核对地区、Population 和 Number。

## 三、计算公式和关键区别

公式（5）的实现：

ASMR(t,c,s) = sum_age [ GBDrate_per_100k(t,c,s,age) / 100000 × w(age) ]

w(age) = [ Population_1980_Male(age) + Population_1980_Female(age) ] / sum_all_ages [ Population_1980_Male(age) + Population_1980_Female(age) ]

年龄别 Rate 已经使用人口分母计算，不需要再次除人口；人口文件用来计算**固定年龄权重**。

**修订前的问题**：旧程序使用 Deaths Number×100000/Rate 反推一个人口代理，再取跨死因中位数作为年龄权重来源。该代数关系在同一人口分母条件下成立，但反推数值不是已核实的 GBD Population 数据。GBD Number 和 Rate 均属估计输出，不能预先认定反推结果恰好等于官方人口。更不能把这个代理称作作者使用的 HCSO 年中暴露量。

**修订后的处理**：run.R 明确要求独立 Population CSV；缺文件就停止，既不默认反推人口，也不生成未经核查的结果。

## 四、在 RStudio 运行

在仓库根目录打开 RStudio，执行：

    source("02_figure2/run.R")

程序逐步检验：
1. 死亡数据的年份、性别、年龄、19 类死因与 Rate/Number 完整性；
2. Population 中 1980 年 12 条年龄×性别人口是否完整；
3. 六个固定权重加总为 1；
4. 19 死因×男女×44 年共 1672 条 ASMR；
5. 各 ASMR 是否处于年龄别率最小值与最大值之间；
6. 画与原文类似的 3×3 男女折线图（注意九类死因**并非**原文那九类）。

输出目录：outputs/02_figure2/

文件包括 standard_age_weights_1980.csv、asmr_19_causes_1980_2023.csv、figure2_selected_9_causes.csv 和 figure2_gbd_asmr_9_panels.png。

## 五、实际状态

**已修订并提交方法及 R 代码；用户尚未提供 GBD Population 数据，且尚未在真实 CSV 上运行 R，因此没有可靠的 ASMR 实测结果或绘制图片。** 后续应核查六个标准年龄权重、人口量级、九类死因实际名称及与原文趋势的异同。

研究限制：GBD 估计与 HCSO 登记数据不同；只用到 60–89 岁、欠缺 90+；人口标准年用 1980 而非 1970；GBD 分类与原文 ICD 章节不同；原文男女权重口径未在 §3.1 中充分说明。