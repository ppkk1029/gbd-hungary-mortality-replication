# §4.2 全部 R 代码：问题、方法与逐步讲解

## 1. 本节到底做什么

§4.1分别拟合每一种死因，得到每年的时期因子。§4.2问：同一性别、同一模型下，哪些死因的时期因子走势相似？先将它们分组，才能在§4.3检验是否适合共享时期因子。本节输出的是历史聚类，不是未来预测，也不证明共同医学病因。

本次用GBD Hungary的1980—2015年、60—89岁六个年龄组、男女各九类死因。训练单元是36×6×9×2=3888。测试期不参与标准化、选K或拟合。GBD分类不是原文ICD九章，不能把本项目结果称作HCSO数值复刻。

## 2. 运行与文件分工

下载整个仓库，进入 `04_results/4_2_clusters/gbd_complete`。在RStudio打开 `run_4_2_complete.R`，点击Source，选择现有 `fitted_rates.csv`。运行过程中自动创建新的时间戳结果目录，不需要安装第三方R包。旧的 `run_4_2_complete_RStudio_launcher.R` 也已改成直接运行同一套R代码。

```r
# 若希望自己指定路径、逐阶段研究：
options(section42.autorun = FALSE)
source('run_4_2_complete.R', encoding = 'UTF-8')
run_section42(input = 'D:/mortality/fitted_rates.csv',
              out = 'D:/mortality/section42_R_new')
```

`out`应是新目录；若其中已有文件，程序停止以防覆盖。命令行也可运行：

```bash
Rscript run_4_2_complete.R fitted_rates.csv new_results
```

| 文件 | 负责的工作 | 重点函数 |
|---|---|---|
| run_4_2_complete.R | 顺序调用全部阶段、输出清单与环境记录 | run_section42 |
| model_functions.R | 输入检查、七模型似然、解析梯度、SVD识别 | read_input、make_problem、fit_model、canonical |
| run_all_models.R | 126组拟合、14次聚类与敏感性、CSV导出 | fit_all_models、cluster_all_models |
| clustering_functions.R | 标准化、DTW、层次聚类、精确medoids、七CVI、ARI | z_rows、dtw_distance、choose_clusters |
| complete_cvi.R | 196组CVI重算核对、294次随机起点PAM swap | complete_cvi |
| verify_and_plot.R | 数值自测、输出检查、可选迁移对照 | algorithm_selftest、verify_results、compare_reference |
| plot_paper_figures.R | Figures 5–8、热图、树状图、全部时期曲线 | plot_paper_figures |

所有估计、距离、聚类、绘图和核验都在R中执行。常规运行不需要迁移前的结果文件。

## 3. 先把“输入数据”与“拟合结果”分清

虽然文件叫 `fitted_rates.csv`，真正进入本次估计的列只有：`sex_name, cause_name, age_name, year, deaths_gbd, exposure_proxy`。已有的 `fitted_count`、`fitted_rate_per_100k` 不参与拟合。

```r
raw <- read_input(input)
d <- input_matrices(raw, 'Male', CAUSES[1])
D <- d$D  # 6行年龄 × 36列年份：GBD估计死亡人数
E <- d$E  # 同样排列：暴露量代理
```

`read_input()`检查3888行、主键不重复、完整的年龄/年份/性别/死因集合以及有限正数。`input_matrices()`按明确的年龄顺序和年份排序，避免CSV行顺序变化造成错位。D/E是中心死亡率近似量；E是代理暴露量，不宣称等于HCSO观测人年。

## 4. 七模型怎样在 R 中估计

`fit_all_models()`的三层循环是：性别→死因→模型，合计126组。每一组都从D、E重新估计，不读取既有参数当成新结果。

| 模型 | 线性预测面 | 本节聚类用的因子 |
|---|---|---|
| Poisson LC1 | log(m)=α+β1κ1 | κ1 |
| Binomial LC1 | logit(q)=α+β1κ1 | κ1 |
| Poisson LC2 | log(m)=α+β1κ1+β2κ2 | κ2 |
| Binomial LC2 | logit(q)=α+β1κ1+β2κ2 | κ2 |
| CBD | logit(q)=κ1+(x−均龄)κ2/10 | κ1 |
| Plat2 | log(m)=α+κ1+(均龄−x)κ2/10 | κ1 |
| Plat3 | 上式加max(均龄−x,0)κ3/10 | κ1 |

年龄载荷除以10是内部数值缩放；还原到未缩放的论文载荷时，对应κ2、κ3要除以10。代码中的factor编号不变。GBD死亡数是连续估计量，因此Poisson/Binomial只是工作似然，不能直接作严格整数计数推断。

对Poisson，均值是 `mu = E * exp(eta)`；对Binomial/CBD，近似期初暴露量 `N = E + D/2`，均值是 `mu = N * plogis(eta)`。`plogis()`是R的逆logit函数。输出Binomial的中心死亡率时，用 `q/(1-q/2)` 作一致的近似逆变换。

`make_problem()`把参数整理为优化器所需的向量，并提供 `calc()` 和 `unpack()`。解析梯度来自残差 `res=(mu-D)/sum(D)`，例如时期因子的梯度是 `t(res) %*% B`。LC和Plat最后一期κ由前35期的负和给出，所以其梯度还要减掉最后一期的梯度。这一步不是直接删去一个未知量的梯度。

优化采用三个起点，`stats::optim(method='BFGS')` 后必要时以 `stats::nlminb()`细化。目标函数使用“去掉与参数无关常数后的负对数似然”，即半deviance/总死亡数。这样与旧目标有同一极小点，但减少一个很大常数掩盖微小优化进步的问题。Poisson用 `expm1()` 计算小量以降低相消误差。最终依靠缩放梯度<1e−6判断通过，不能只看优化器返回码。

### 为什么还要 SVD？

同一个LC拟合面可以有多套B、K。`canonical()`先把K的列均值移进α，再对 `B %*% t(K)` 做SVD，令B列正交、奇异值降序，并令每列最大绝对载荷为正。这样不同死因与不同语言的结果可以在同一约定下比较。但这个约定不证明κ2是唯一真实风险因素；旋转检验正是检查该局限。

## 5. 从时期因子到 DTW 距离

```r
factor <- if (endsWith(model, 'LC2')) 2L else 1L
X <- do.call(rbind, lapply(CAUSES, function(c)
  fits[[key(sex, c, model)]]$K[, factor]))
Z <- z_rows(X)
d <- distance_matrix(Z)
```

X为9×36：每行一种死因、每列一年。`z_rows()`对每一行减均值、除样本标准差（分母35）。这样主要比较形状，不能跨模型直接比较原κ大小；也不能靠标准化解决符号或旋转问题。

DTW并不是把某一年的死亡率改成年份不同的观测值，而是为两条历史时期因子找一条有序的匹配路径。本项目明确采用局部平方差：

$$
c_{ij}=(z_i-z_j)^2.
$$

累计代价的递推为：

$$
A_{ij}=\min\{A_{i-1,j-1}+2c_{ij},\;A_{i-1,j}+c_{ij},\;A_{i,j-1}+c_{ij}\}.
$$

初始 `A[1,1]=c[1,1]`，其它不可达位置为Inf。端点闭合，主分析无窗口。最终距离 `sqrt(A[n,m]/(n+m))`。`prev`矩阵记录最小值来自哪一格，从右下角反向追踪就得到对齐路径。图上为可读性显示绝对差 `abs(z1-z2)`，但路径和定量距离用平方差；必须区分这两件事。

`distance_matrix()`每次计算九类的36个不同配对，形成对称9×9矩阵。男女×七模型共14组，输出504条距离。窗口5年与L1版本属于敏感性，不混入主距离。

## 6. 如何分组、如何选 K

```r
hc <- hclust(as.dist(d), method = 'average')
lab <- cutree(hc, k = 2)
exact <- exact_medoids(d, k = 2)
a <- choose_clusters(d, Z)
```

层次聚类的average指合并两组时使用所有跨组配对的平均距离。另一条路线选择k个真实序列作medoids，每个死因分给最近的medoid。九条序列允许枚举所有组合，从而得到该距离矩阵上的最小组内目标；这不是用一个随机起点的PAM局部解代替全局最优。

每个性别/模型比较 K=2,…,8 和两种方法，共14个候选；总共196个。主规则是平均Silhouette最大，平局取较小K，再按方法名字排序。单元素簇的Silhouette记0。代码也输出Dunn、COP、DB、DBstar、CH、SF；CH/SF等使用medoid版本定义，另附平方距离CH和欧氏CH，不能把不同几何的指标混成一张选票。

Silhouette是组内与组间分离程度，不是预测误差。原文综合多项指标与树状图做判断；本项目显式使用最大Silhouette，因此选K规则存在差异。任何组号G1、G2只在该行模型/性别内有效。

`complete_cvi()`另从三个随机medoid集合做PAM swap，2×7×7×3=294次，比较其目标与枚举最优的差距。固定seed只保证本R实现复跑可重复；不同语言即使seed相同也不保证随机起点相同。

## 7. Figures 5–8 各自读取什么

| 图 | 人口与模型 | 画什么 | 从哪里来 |
|---|---|---|---|
| 5 | 男性 Poisson LC2 κ2 | CVD–DIG曲线、局部差异、路径 | analyses中的Z和dtw_distance(path=TRUE) |
| 6 | 女性 Poisson LC2 κ2 | CVD–DKD曲线、局部差异、路径 | 同上，改性别与配对 |
| 7 | 男性 CBD κ1 | 每簇一个面板的时期曲线 | analyses中的Z和labels |
| 8 | 男性 Binomial LC2 κ2 | 每簇一个面板的时期曲线 | 同上，改模型与因子 |

所有图同时保存PNG和SVG，绘图使用完整精度的R结果，没有把小数提前截到两位。组数由实际labels决定，没有硬编码必须两簇。CVD–DIG与CVD–DKD作为论文对照配对保留，即使其中一对在GBD并不接近，也不换成“最好看”的配对。

## 8. 怎样检查程序是否真的完成

`validation.csv`记录算法自测、126组梯度检查及旋转前后拟合面差。核心输出数量检查包括8424个完整κ值、4536个用于聚类的标准化值、504个距离和126条死因分组。

独立自测用有限差分验证七模型解析梯度，用所有单调路径的完整枚举验证小型DTW例子，而非只让同一个函数与自己比较。三起点只检查本次初值稳定性，不是似然全局最优证明。

`manifest.txt`记录输入MD5、年份、年龄、seed、DTW和选K约定，`sessionInfo.txt`记录R环境。迁移验收额外读取只读旧CSV，生成 `migration_numeric_comparison.csv` 与 `migration_membership_comparison.csv`；它们不会影响拟合或分组。实际结果与精度见[运行核验](R_MIGRATION_VALIDATION.md)。

## 9. 怎么向导师解释结果与差异

依次讲：为何需要候选共享组→如何拟合/提取κ→为何标准化和DTW→怎样选K→展示图和分组→与原文比较→哪些结论仍未检验。男性Figure 5的指定配对在当前GBD下并不接近；女性Figure 6相对接近但DKD不等于原文ICD IV；Figures 7和8的男性簇数与原文不同。不能把这些压缩成“与原文大致一致”。

双因子同时旋转B、K，拟合面保持不变，但κ2的分组会改变。五年窗口和删除女性DKD的检查也都要报告。历史相似性尚未经过§4.3共享模型和样本外回测验证；当前只能把分组作为候选共享结构。
