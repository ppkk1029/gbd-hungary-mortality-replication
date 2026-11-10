# §4.2 | RStudio graphical launcher for the COMPLETE seven-model Python experiment.
# This file does NOT reimplement Poisson/Binomial/CBD/Plat estimation in R.
# Instead it invokes the same Python/SciPy pipeline used for the archived
# 126 actual model fits and 14 DTW cluster analyses.
#
# Requirements: Python 3.10+ on PATH with numpy pandas scipy matplotlib scikit-learn
# In RStudio: source this R file; choose the existing run_all_models.py; choose
# the ORIGINAL fitted_rates.csv (columns deaths_gbd/exposure_proxy).
# Output is always a "results" directory NEXT TO the selected Python scripts,
# so that plot_paper_figures.py can load it consistently.
# BEWARE: the output directory is overwritten on confirmation; back it up first.

choose_python <- function() {
  cands <- c(Sys.which("python"), Sys.which("python3"), Sys.which("py"))
  cands <- cands[nzchar(cands)]
  if (!length(cands)) stop("找不到 Python：先安装 Python 3.10+ 和所需包。")
  cands[[1]]
}
cat("步骤 1：选择仓库中的 gbd_complete/run_all_models.py\n")
script <- normalizePath(file.choose(),winslash="/",mustWork=TRUE)
if (basename(script)!="run_all_models.py") stop("请选择 run_all_models.py")
srcdir <- dirname(script)
required <- c("run_all_models.py","verify_and_plot.py",
              "complete_cvi.py","plot_paper_figures.py")
if (!all(file.exists(file.path(srcdir,required)))) stop(
  "该文件夹缺少完整 Python 脚本：",paste(setdiff(required,list.files(srcdir)),collapse=", "))
cat("步骤 2：选择原始 GBD fitted_rates.csv（不得选择只含 fitted_count 的预测输出）\n")
raw <- normalizePath(file.choose(),winslash="/",mustWork=TRUE)
# Only a quick header check here; all 3,888 data cells are checked in Python.
hdr <- names(read.csv(raw,nrows=1,check.names=FALSE))
fields <- c("sex_name","cause_name","age_name","year","deaths_gbd","exposure_proxy")
if (!all(fields %in% hdr)) stop("输入 CSV 缺少必需列：",paste(setdiff(fields,hdr),collapse=", "))
out <- file.path(srcdir,"results")
cat("\n源程序：",srcdir,"\n输入文件：",raw,
    "\n结果目录：",out,
    "\n本操作将估计 126 组模型，执行 14 组聚类并生成 Figure 5/6/7/8。\n")
ans <- tolower(trimws(readline("输入 y 确认执行及允许覆盖同名结果；其他键取消：")))
if (ans!="y") stop("已取消，没有执行模型或写文件。")
py <- choose_python()
run_step <- function(name,argv) {
  cat("\n运行：",name,"\n")
  # system2 supports quoting paths containing spaces.
  lines <- suppressWarnings(system2(py,c(shQuote(file.path(srcdir,name)),
                                     vapply(argv,shQuote,character(1))),
                      stdout=TRUE,stderr=TRUE))
  cat(paste(lines,collapse="\n"),"\n")
  code <- attr(lines,"status")
  if (!is.null(code) && code!=0) stop(name," 返回错误码：",code)
}
run_step("run_all_models.py",c(raw,out))
run_step("verify_and_plot.py",c(out))
run_step("complete_cvi.py",c(out))
run_step("plot_paper_figures.py",character(0))
needed <- c("fit_diagnostics_all_models.csv","cluster_summary.csv",
            "cluster_membership.csv","Figure_5_GBD_Male_Poisson_LC2.png",
            "Figure_6_GBD_Female_Poisson_LC2.png",
            "Figure_7_GBD_Male_CBD_Clusters.png",
            "Figure_8_GBD_Male_Binomial_LC2.png","validation.json")
if(!all(file.exists(file.path(out,needed)))) stop(
  "程序执行完毕但部分结果缺失：",
  paste(setdiff(needed,list.files(out)),collapse=", "))
cat("\n复现完毕，确认结果：\n")
print(read.csv(file.path(out,"cluster_summary.csv"))[,c("sex","model","K","silhouette")])
cat("\n所有文件保存在：",normalizePath(out,winslash="/"),"\n")
