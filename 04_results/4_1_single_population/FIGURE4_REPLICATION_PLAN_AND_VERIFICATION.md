# §4.1 | Figure 4 | 双因子 Poisson Lee–Carter：复现方法、执行证据与状态

## 一、原文对应位置
- Varga (2025), Table 1 与 §4.1 Figure 4: **Poisson Lee–Carter with two factors for males**。
- 模型：log(m[x,t,c]) = alpha[x,c] + beta1[x,c]*kappa1[t,c] + beta2[x,c]*kappa2[t,c]。
- 数据分布：D[x,t,c] 的条件均值为 E[x,t]*m[x,t,c]，使用 Poisson 均值结构估计；GBD Deaths Number 属于平滑后的连续估计，故本项目不能把概率解释成严格的真实整数观测 Poisson likelihood。
- 原文使用 HCSO 数据与更多年龄，GBD 复现只能验证方法，不能期待图像数值相同。

## 二、Figure 4 本项目的数据入口

**输入：此前已经在 RStudio 得到的 `fitted_rates.csv`**，其列包括：
`sex_name,cause_name,age_name,year,deaths_gbd,exposure_proxy,fitted_count`。

该文件中 `deaths_gbd` 与 `exposure_proxy` 是用于重新估计的观测输入；`fitted_count` **仅用于和先前单因子拟合比较**，不得将 `fitted_count` 作为双因子的观测死亡人数。

核验：6 年龄组 × 36 年(1980–2015) × 9 死因 × 2 性别 = **3888 行**；每个模型216个单元。2016–2019 留作测试，2020–2023 未用来拟合。

### 运行
1. 下载 [Figure 4 双因子 R 脚本](run_4_1_poisson_lc_2factor.R)；在 RStudio 中打开后点击 Source，选择之前已保存的 **fitted_rates.csv**。不需要重复下载原始 20064 行 GBD 数据。
2. Console 打印 18 组实际模型对比与约束检查；Plots 显示男性五种参数、女性五种参数，以及一因子/二因子的训练期 deviance 差。
3. 检查完成后，在 Console 输入 `y`、选择保存父目录；生成 `Figure4_Poisson_LC_2factor/` 文件夹。
4. 或在已安装 R 的系统终端执行：`Rscript run_4_1_poisson_lc_2factor.R /path/fitted_rates.csv /path/results`。

### 输出及作用

| 文件 | 说明 |
|---|---|
| `figure4_male.png` | 对应原文 Figure 4 的男性 α、β¹、β²、κ¹、κ² 参数曲线；九死因 |
| `figure4_female.png` | 女性同样结构，作为额外复现 |
| `figure4_deviance_comparison.png` | 两模型训练期拟合损失对照（不是预测误差） |
| `fit_diagnostics_2factor.csv` | 18 组 deviance、pseudo-AIC/BIC差、收敛、二因子强度、参数异常 |
| `alpha_2factor.csv, beta_2factor.csv, kappa_2factor.csv` | 18 组模型的新估计参数 |
| `fitted_rates_2factor.csv` | 每个年龄×年单元预测死亡人数、Poisson 形式残差 |
| `kappa_factor2_for_DTW.csv` | §4.2 后续可考虑聚类的第二时期因子，但须先检查稳定性 |

## 三、为可识别性增加的数学约束

单因子 `sum(beta)=1` 在近乎正负抵消时会膨胀。双因子还具有更深的旋转不唯一：任意可逆 2×2 矩阵 H 可以使 `B*K^T=(B H)*(K (H^{-1})^T)^T` 保持不变。因此本项目优化时固定两处参考年龄的 `beta` 子矩阵为单位阵，并令 `mean(kappa1)=mean(kappa2)=0`；优化后对拟合的因子表面作 SVD 正交分解、奇异值降序及符号固定。**与原文 Table 1 的单纯 sum(beta)=1 比例不同**，这是为数值稳定和可重复性采取的额外约定。程序检查转化前后每个年龄×年份的对数死亡率绝对偏差 <1e-7。

如果两奇异值相近，两个因素的次序不稳定；如果第二奇异值非常小，则第二因素可能主要是噪声，须标记。即使拟合更好，样本外不一定更好。

## 四、模型参数与有效自由度

6 个年龄 × 36 年共 216 个观测单元：

- 一因子可识别参数数目：6 + (6−1) + (36−1) = **46**。
- 双因子可识别参数数目：6 + 2×(6+36−1−2) = **84**。

双因子比一因子多38个参数，训练期 Poisson deviance 下降不代表一定具有泛化能力。我们另外比较伪 AIC/BIC，测试期误差仍须另行验证。

## 五、已进行的独立数值自检（不是 GBD 实证）

在用已知双因子结构生成的 Poisson 模拟死亡人数上，用相同的 anchor+center 思路以 Python/SciPy 作独立目标函数/梯度检查（**本次运行环境没有 R，不能说 RStudio 的实际脚本已运行**）：

| 检查 | 实际数值 |
|---|---:|
| 6×36 模拟数据，一因子 Poisson deviance | 3545.795811 |
| 同数据，双因子 Poisson deviance | 152.445482 |
| 双因子最大有限差分/解析梯度误差（抽样参数） | 1.21×10⁻¹² |
| 双因子重参数化后对数死亡率最大误差 | 8.88×10⁻¹⁶ |
| Python 优化是否报告成功 | 是（55次迭代） |

注意：以上只说明**原理自检通过**，不能将它们写作匈牙利女性/男性九死因的改进结果，更不能说明预测优于一因子模型。R 的逐组优化仍需通过用户已有 CSV 实际执行后验证。

[数值自检示意图](figure4_synthetic_numeric_selftest.svg)（注明 synthetic，不是 Figure 4 真实结果）。

## 六、真实数据当前状态

- **已完成**：Figure 3 的 18 组实测参数及异常诊断；Figure 4 全套重估 R 程序提交；双因子实现原理的独立模拟数值自检。
- **未完成**：Figure 4 18 组真实 GBD 数据上的模型重估、男性/女性真实 Figure 4 参数图、和七模型比较。
- **具体受阻**：之前保存的 `fitted_rates.csv` 可在个人 Library 找到，但本会话计算环境无法取得它的原始字节；无法在这里执行双因子实证拟合。不能仅凭一因子参数画出双因子真实结果。
- **后续 §4.2**：二因子 Poisson LC 的第二时期因子可以进入 DTW，但必须确认 SVD 旋转与识别的稳健性，不将未检验的因素曲线用于正式分组。

论文：[Varga (2025) full text](https://www.tandfonline.com/doi/full/10.1080/03461238.2025.2480846)。
