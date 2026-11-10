# Figure 2｜年龄标准化分死因死亡率：从论文公式到 R 代码（逐步学习笔记）

> **定位：学习与复现说明，不替代运行日志。** 该文档逐块解释同目录 [run.R](run.R) 的实际做法，方便向导师解释“问题—理论—数据—计算—核验—作图—与原文比较”。编写说明不等于已经在本地重新运行验证 R 脚本。原文：[Varga (2025), §3.1, Figure 2](https://doi.org/10.1080/03461238.2025.2480846)；对应 [Figure 1 目录](../01_figure1/README.md)。

## 0. 一句话说清楚这个图的研究问题

- Figure 1：某死因死亡人数 / 全部死因死亡人数，是死亡**构成**，不能直接判断风险。
- Figure 2：固定年龄结构后比较某死因在每年、男性或女性中的**年龄标准化死亡率**（Age-Standardized Mortality Rate，ASMR）。
- 为什么要标准化？如果高龄人口占比增加，即使每个年龄组的死亡风险完全不变，总体粗死亡率也会增加；固定年龄权重，才能更清晰地观察年龄别死亡率随时间的变化。

**两年龄教学例子（不是匈牙利实际数据）**：60–69岁死亡率1%，80–89岁死亡率10%；第一年人口比例80%/20%，第二年30%/70%。两个年龄组死亡率完全未变，但粗死亡率由2.8%变为7.3%。若两年都按80%/20%权重计算，ASMR均为2.8%。

## 1. 论文公式、口径差异和符号

对年龄组 $x$、年份 $t$、性别 $s$、死因 $c$，年龄别中心死亡率约为

$$
m_{x,t,s,c}=\frac{D_{x,t,s,c}}{E^c_{x,t,s}}.
$$

固定标准人口权重 $w_x$（加总等于1）后，

$$
ASMR_{t,s,c}=\sum_{x\in\{60-64,\dots,85-89\}}w_x m_{x,t,s,c}.
$$

原文使用1970年的标准人口、HCSO登记数据、60岁及以上（含90+）、九个ICD章。**本项目自行选择**GBD 2023 Hungary、1980—2023年、六个60—89岁组、1980年男女合并人口的固定权重、19个GBD Level 2死因（作图选其中九个）。男性和女性死亡率**分别计算**，但共享同一套 $w_x$。所以这是**方法复现**，不是原文数值原样复刻。

GBD的下载 Rate 单位是“每10万人”，R程序转换为每人单位：

$$
m_{x,t,s,c}=\frac{\mathrm{Rate}_{x,t,s,c}}{100000}.
$$

**Rate本身已经除过人口，不能再除一次 Population。** 独立 Population 主要用于构造权重和交叉核对 Number/Rate 的一致性。

## 2. 全部代码执行路线（记住先后依赖）

| 步骤 | run.R 位置（以本次读取的版本为准） | 输入→输出 | 要回答的问题 |
| --- | --- | --- | --- |
| 0 | 13–39行 | 选择 Population/Deaths 文件、设年份和年龄 | 数据从哪里来？ |
| 1 | 41–75行 | Population 44年×2性别×6组＝528行 → pop | 人口范围、缺失、重复、异常是否正确？ |
| 2 | 77–98行 | 1980年的12行人口 → 六个权重 standard | 固定年龄结构是什么？ |
| 3 | 100–141行 | Deaths Number/Rate 20,064行 → paired | 19类死因数据是否匹配？ |
| 4 | 143–184行 | Number、Rate、独立 Population → qc | 两种数据来源有没有口径/单位错配？ |
| 5 | 185–214行 | Rate × 权重 → asmr，共1672行 | 每个死因、年份、性别的ASMR怎么算？ |
| 6 | 215–220行 | 1980–2015年死亡人数 → 九类 top9 | 为什么只画这九类？ |
| 7 | 222–302行 | 权重、QC、ASMR → 预览与3×3九面板 | 男女趋势与质量检查结果是否合理？ |
| 8–9 | 303–380行 | 输入 y + 选择目录 → CSV、PNG、摘要 | 用户确认后保存了什么？ |

## 3. 步骤0：选择两个不同的文件

~~~r
population_path <- file.choose()
deaths_path <- file.choose()
years <- 1980:2023
standard_year <- 1980L
sexes <- c("Male", "Female")
ages <- c("60-64 years","65-69 years","70-74 years",
          "75-79 years","80-84 years","85-89 years")
~~~

- 第一次选**Population人口CSV**，第二次选**Deaths CSV（必须包含Number与Rate）**，不是两次都选死亡数据。
- 文件名可更改，但其字段内容、国家、年份、年龄和性别必须符合核验规则。
- Population共 $44\times2\times6=528$ 行；Deaths共 $44\times2\times6\times19\times2=20064$ 行（最后乘2是Number和Rate两种指标）。
- 使用 file.choose() 是方便 RStudio 手动选文件；与 Figure 1 使用固定 data/raw 路径不同。

## 4. 步骤1：读人口并做范围、重复和数值核验

~~~r
pop <- read.csv(population_path, check.names=FALSE,
                stringsAsFactors=FALSE)
stopifnot(all(c("year","location_name","sex_name",
                "age_name","measure_name","metric_name",
                "val","lower","upper") %in% names(pop)))
~~~

随后 run.R 验证：Hungary；measure=Population；metric=Number；年份正好1980–2023；2个性别、6个年龄组；528行；(year,sex,age)组合不能重复；val必须正且有限；lower ≤ val ≤ upper。可选人口组字段必须为 All Population。

**R语法**：read.csv读取表；stopifnot条件不满足就停止；setequal验证两个值集合的成员是否相同；anyDuplicated查重复。这里的目的是排除“不相关数据也能画出一张看似正常的图”。

## 5. 步骤2：固定1980年年龄结构，得到六个权重

~~~r
base <- pop[pop$year == standard_year, ]
standard <- aggregate(population ~ age_name,
                      data=base,FUN=sum)
standard <- standard[match(ages,standard$age_name),]
standard$weight <- standard$population_1980_both_sexes /
                   sum(standard$population_1980_both_sexes)
stopifnot(abs(sum(standard$weight)-1)<1e-12)
~~~

上面最后两句省略了**原脚本第83行**将 aggregate 输出列改名为 population_1980_both_sexes 的过程；完整能直接运行的版本以 [run.R](run.R) 为准。aggregate按年龄把男女人口合并，match固定年龄排序；将每个年龄组人口除以1980年60–89岁人口总和，就是标准权重。

**独立核验的真实Population结果**：

| 年龄 | 1980男女合计人口 | 权重 |
| --- | ---: | ---: |
| 60–64 | 379402.238 | 0.2120072 |
| 65–69 | 536031.002 | 0.2995302 |
| 70–74 | 411662.203 | 0.2300338 |
| 75–79 | 270187.230 | 0.1509786 |
| 80–84 | 142485.453 | 0.0796198 |
| 85–89 | 49804.550 | 0.0278304 |
| 合计 | 1789572.676 | 1 |

举例：$w_{60-64}=379402.238/1789572.676\approx0.2120072$。**以后计算任何年份、任何性别、任何死因，这套年龄权重都不变**。但使用的年龄别死亡率来自各自的年份、性别和死因。

> 注意原脚本还带有一个对人口总和1789572.6757876906的硬编码核验（误差1以内）。它检查是否正在使用这批特定数据，不具有跨数据版本通用性；换数据库或版本需要有依据地更新此项，不能随意删除全部质量核验。

## 6. 步骤3：把Deaths的Number和Rate正确配对

~~~r
keys <- c("year","sex_name","age_name","cause_name")
rate <- gbd[gbd$metric_name=="Rate",c(keys,"val")]
num  <- gbd[gbd$metric_name=="Number",c(keys,"val")]
names(rate)[names(rate)=="val"] <- "rate_per_100k"
names(num)[names(num)=="val"] <- "death_number"
paired <- merge(rate,num,by=keys,all=FALSE,sort=FALSE)
~~~

- gbd 是另一个CSV；先核验年份、国家、年龄、性别、19类死因、Number/Rate、重复/缺失、数值及区间上下界。
- merge **必须同时按年份、性别、年龄和死因四个键**连接，不能只按年份！否则一个年龄组的死亡率可能匹配到另一个年龄组的死亡人数。
- 共10032对同单元Number/Rate，连接后应该仍为10032行。

## 7. 步骤4：独立人口交叉核验（这是质量控制，不是ASMR公式）

~~~r
qc <- merge(paired,
            pop[c("year","sex_name","age_name","population")],
            by=c("year","sex_name","age_name"),
            all.x=TRUE,sort=FALSE)
qc$rate_recomputed_per_100k <-
   100000*qc$death_number/qc$population
qc$relative_error <- ifelse(qc$rate_per_100k>0,
     abs(qc$rate_recomputed_per_100k-qc$rate_per_100k) /
          qc$rate_per_100k,
     ifelse(qc$rate_recomputed_per_100k==0,0,Inf))
~~~

- 用独立人口重算每10万人死亡率，与Deaths文件中的Rate比较；中位相对差>2%或第95百分位>10%时停止。这是**当前程序自设**的质量警戒阈值，不是论文统计检验阈值。
- 当GBD分别报告死亡人数、人口和死亡率的估计均值时，均值的比值未必严格等于比值的均值，所以**不能要求所有单元误差为0**。
- 若异常，优先检查国家、GDB版本、性别、年龄段、死亡指标和Rate单位，不要直接调松阈值掩盖错误。
- 输出 number_rate_population_check.csv 是给导师检查的数据溯源表。

## 8. 步骤5：真正的ASMR核心代码（务必会口头解释）

~~~r
rates_w <- merge(rate,standard[c("age_name","weight")],
                 by="age_name",all.x=TRUE,sort=FALSE)
rates_w$per_person_year <- rates_w$rate_per_100k / 100000
rates_w$weighted <- rates_w$per_person_year * rates_w$weight
asmr <- aggregate(weighted ~ year+sex_name+cause_name,
                  data=rates_w,FUN=sum)
names(asmr)[names(asmr)=="weighted"] <- "asmr_per_person_year"
asmr$asmr_per_100k <- 100000*asmr$asmr_per_person_year
~~~

依次解释：

1. merge：给所有年份的60–64岁数据配上同一个1980年权重0.2120072，其他五组同理；这是关键，因为要保持权重不随年变。
2. Rate/100000：从每10万人死亡率转成每人单位。**不能再次除人口**。
3. weighted = per_person_year × weight：算每个年龄组对标准化率的贡献。
4. aggregate(...year+sex_name+cause_name,FUN=sum)：将六组贡献相加，得到“某一年×某性别×某死因”的一个ASMR。
5. *100000：仅用于同时输出一个更便于阅读的每10万人版本，和原来的每人单位是一回事。

**36? 不是。** Figure 2覆盖44年（1980—2023），19类，男女各一套，因此 $44\times2\times19=1672$ 个ASMR；§4.1/§4.2训练期才是36年（1980—2015）。这两个年份范围不能混。

**实际权重代入演示**（六个率只是符号，不是编造观察值）：

$$
ASMR_{2019,\mathrm{Male},\mathrm{CVD}}=
0.2120072m_{60-64,2019,M,CVD}
+0.2995302m_{65-69,2019,M,CVD}
+\cdots+0.0278304m_{85-89,2019,M,CVD}.
$$

因为它是正权重之和为1的加权平均，所以每个ASMR应处于对应六个年龄别率的最小值与最大值之间；run.R第203–214行明确做了这个检查。

## 9. 步骤6：为什么全部19类算完只画九个面板？

~~~r
rankdata <- num[num$year<=2015,c("cause_name","death_number")]
ranking <- aggregate(death_number~cause_name,data=rankdata,FUN=sum)
ranking <- ranking[order(-ranking$death_number,ranking$cause_name),]
top9 <- head(ranking$cause_name,9)
~~~

采用**仅训练期1980–2015年的累计死亡人数**排序选九类，以与 Figure 1 保持一致并避免事后利用未来信息挑选类别。这是我们的展示规则，**不是原文从ICD章节中挑选九类的完全相同口径**；所有19类的ASMR仍计算并保存，不必只保留九类。

## 10. 步骤7：为什么是3×3九面板，每图男女两条线？

~~~r
par(mfrow=c(3,3))
for (i in seq_along(top9)) {
  z <- asmr[asmr$cause_name==top9[i],]
  male <- z[z$sex_name=="Male",]
  female <- z[z$sex_name=="Female",]
  plot(male$year,male$asmr_per_person_year,
       type="l",lty=1)
  lines(female$year,female$asmr_per_person_year,lty=2)
}
~~~

这段是解释原代码画法的简化示例：原脚本另外还有年份排序、颜色、图例、自动坐标范围、检查1980—2023年完整性等。par(mfrow=c(3,3))分九格；每格固定同一死因，男性为实线、女性为虚线。**比较时优先看同一格内男女与同一死因跨年的变化，不能误以为九格共享同一纵轴范围。**

预览包括：六个1980年权重、Number/Rate/Population差异最大单元、九死因清单、1980和2023 ASMR表，以及九幅图。仅预览时不保存任何文件。

## 11. 步骤8–9：如何本机运行与保存（不要把代码存在当作运行成功）

1. 在RStudio打开当前文件 [02_figure2/run.R](run.R)，点击 Source；或在Console运行 source(file.choose()) 后选择这份R程序。
2. 依次选 Population CSV、Deaths CSV（含Number和Rate）。二者不是同一个文件。
3. 检查 Console：528条 Population、六权重和≈1、10032条Number/Rate配对、质控差异、九个死因、1980和2023 ASMR，检查Plots中的九面板图。
4. **自行判断合理后**在Console输入 y；再选一个保存父文件夹。程序新建 Figure2_ASMR_results 子文件夹。如果已经存在，会再次确认是否覆盖同名文件。
5. 程序最后显示 RESULTS SUCCESSFULLY SAVED 和保存路径，才算把这次运行的文件写出完成。

| 输出文件 | 用途 |
| --- | --- |
| population_1980_2023_clean.csv | 清洗并核验的528条人口数据 |
| standard_age_weights_1980.csv | 六个1980年标准权重 |
| number_rate_population_check.csv | 10032条Number/Rate/Population逐行核验 |
| asmr_19_causes_1980_2023.csv | 1672条按性别/死因/年份的ASMR |
| figure2_selected_9_causes.csv | 选入九面板的死因及其顺序 |
| figure2_gbd_asmr_9_panels.png | 九面板折线图 |
| run_summary.txt | 输入文件、运行日期、权重基准、QC简报 |

## 12. 与原文比较，哪些结论能说，哪些不能说？

**能说**：标准化分析避免年龄构成直接干扰；心血管等死因的长期改善与男女差异可从对应曲线描述；2021年的感染相关突变值得单独说明。

**不能说**：这是原文Figure 2的相同数值；不同分类定义的两组曲线完全可比；ASMR变化就是疾病因果机制变化；Figure 2直接提供后续Poisson LC的输入。

后续 §4.1 的Poisson LC需要**按年龄 × 年份**保留死亡人数 D 和暴露量 E，不能直接把Figure 2已汇总成一个ASMR的结果交给LC拟合；§4.2比较的是模型重新估计出的时期因子κ，不是ASMR线本身。

**现有证据状态**：GBD人口CSV和权重在先前工作中独立核验过，Notion中保留有历史Figure 2图；但截至编写本文件，尚未得到本地RStudio对当前run.R版本完整重跑及全部输出逐值核验的证据。请不要将历史图或上传脚本替代本机运行日志。

## 13. 向导师汇报：1分钟讲清楚计算原理与代码

“Figure 1只反映死亡人数的构成，不等于死亡风险。Figure 2为了比较不同年份、男女不同死因的死亡率趋势，采用直接年龄标准化。我们下载匈牙利GBD Population和Deaths文件，先核对数据的年份、年龄、性别、死因和Number/Rate指标，再把1980年60—89岁男女合并人口按六年龄组计算固定权重。对每一年、每性别、每死因，把GBD各年龄别Rate转成每人单位，逐组乘固定年龄权重，再相加得到ASMR；最后用训练期死亡人数选九个主要死因，分别绘制男女趋势。代码另外检查Number/Population与GBD Rate的一致性，并核验每个ASMR位于六年龄别率的范围内。由于原文使用1970年标准人口、HCSO以及不同死因分类，这里是方法复现，不是原文数值复刻。”

**自查题（建议先独立口述）**：

1. 为什么不能逐年重新计算年龄权重？（会把人口结构变化重新混入年度比较）
2. 男女性别的标准权重为什么相同，但ASMR通常不同？（使用同样权重，却乘各自的年龄别死亡率）
3. GBD原始Rate为什么不能再除Population？（Rate已经是以人口为分母的率）
4. 为什么不能用Figure 2的ASMR直接拟合多年龄的Poisson LC？（年龄维度已被聚合）
5. 为什么这个脚本保存前要先预览？（避免未经质控的图与数值被误当成正式复现）
