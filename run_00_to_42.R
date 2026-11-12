# ================================================================
# Varga (2025) Hungary GBD method replication | Chapters 00 to 4.2
# RStudio unified launcher (UTF-8). Open this file -> click Source.
# Download the ENTIRE repository and keep its folder structure.
#
# Stage 00: validate GBD L1/L2 when ORIGINAL files are present.
# Stage 3.1: Figure 1 when ORIGINAL files are present.
# Stage 3.1: Figure 2 optional (select Population and Deaths CSV).
# Stage 4.1: optional Figure 3/4 interactive standalone plots.
# Stage 4.1/4.2: 7-model period factors -> DTW distance -> 7 CVIs,
#                hierarchical/PAM, 3 seeds and candidate partitions.
# No auto-K: paper does not disclose a rule combining seven indices.
# ================================================================

# --------------------- USER SETTINGS -----------------------------
RUN_DATA_CHECKS <- TRUE
RUN_FIGURE1 <- TRUE
RUN_FIGURE2 <- FALSE          # TRUE prompts to select Population + Deaths
RUN_FIGURE3 <- FALSE          # TRUE prompts for 2 files; separate visualization
RUN_FIGURE4 <- FALSE          # TRUE requires Figure 3's fitted_rates.csv
REFIT_ALL_41_MODELS <- FALSE  # TRUE refits 126 models from original fitted_rates.csv

# If REFIT_ALL_41_MODELS is TRUE, provide full local path or leave blank
# to choose fitted_rates.csv interactively. Otherwise §4.1 saved factors
# are read from "04_results/4_1_single_population/results/".
FITTED_RATES_CSV <- ""

# Optional: after reviewing all CVIs, fill the 14 rows in
# cluster_selection_template.csv and save as reviewed_selection.csv.
# Only then enable the option below to create formal clusters.
FINALIZE_REVIEWED_K <- FALSE
REVIEWED_SELECTION_CSV <- ""  # e.g. "D:/research/Section42.../reviewed_selection.csv"
PREVIOUS_SECTION42_OUTPUT <- "" # for finalizing an existing run only

# ---------------------- RESOLVE FOLDER ---------------------------
this_file <- local({
  of <- lapply(sys.frames(),function(f) f$ofile)
  of <- Filter(function(x) is.character(x) && length(x)==1L && nzchar(x),of)
  if(length(of)) normalizePath(tail(of,1L)[[1L]],mustWork=TRUE)
  else {
    z <- grep("^--file=",commandArgs(),value=TRUE)
    if(length(z))normalizePath(sub("^--file=","",z[[1L]]),mustWork=TRUE)
    else stop("请打开run_00_to_42.R并点击Source；不能只复制其中几行到Console")
  }
})
repo_root <- dirname(this_file)
stopifnot(file.exists(file.path(repo_root,"04_results","4_1_single_population",
                              "model_functions.R")))
setwd(repo_root)
cat("\nRepository: ",repo_root,"\n",sep="")

# ---------------------- PACKAGES ---------------------------------
for(pkg in c("dtw","dtwclust")) {
  if(!requireNamespace(pkg,quietly=TRUE)) {
    message("正在安装R包：",pkg)
    install.packages(pkg,repos="https://cloud.r-project.org")
  }
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("未成功安装R包 ",pkg,"。请检查网络，手动安装后重新Source。")
}
cat("dtw = ",as.character(packageVersion("dtw")),
    " ; dtwclust = ",as.character(packageVersion("dtwclust")),"\n",sep="")

# Optional finalization of a previously computed §4.2 directory.
old_autorun <- getOption("section42.autorun",TRUE)
options(section42.autorun=FALSE)
source(file.path(repo_root,"04_results","4_2_clusters","gbd_complete",
                 "run_4_2_complete.R"),local=.GlobalEnv,encoding="UTF-8")
options(section42.autorun=old_autorun)
if(FINALIZE_REVIEWED_K && nzchar(PREVIOUS_SECTION42_OUTPUT)) {
  if(!nzchar(REVIEWED_SELECTION_CSV))
    stop("请设置REVIEWED_SELECTION_CSV为填写完的14行决策表")
  finalize_section42(PREVIOUS_SECTION42_OUTPUT,REVIEWED_SELECTION_CSV)
} else {
  # ----------------------- STAGE 00 -------------------------------
  level1 <- file.path(repo_root,"data","raw",
             "IHME-GBD_2023_DATA-97ee48e8-1.csv")
  level2 <- file.path(repo_root,"data","raw",
             "IHME-GBD_2023_DATA-b40ffc46-1.csv")
  raw_available <- file.exists(level1) && file.exists(level2)
  if(RUN_DATA_CHECKS) {
    if(raw_available) {
      cat("\n[00] 运行原始GBD数据一致性核验...\n")
      source(file.path(repo_root,"00_data","check_data.R"),
             local=new.env(parent=globalenv()),encoding="UTF-8")
      source(file.path(repo_root,"00_data","check_level1_level2.R"),
             local=new.env(parent=globalenv()),encoding="UTF-8")
    } else {
      cat("\n[00] 跳过：缺少GBD Level1/Level2原始CSV；没有假装通过核验。\n")
    }
  }

  # ----------------------- STAGE 3.1 ------------------------------
  if(RUN_FIGURE1) {
    if(raw_available) {
      cat("\n[3.1 Figure 1] 绘制死因构成面积图...\n")
      source(file.path(repo_root,"03_1_descriptive_statistics","figure1","run.R"),
             local=new.env(parent=globalenv()),encoding="UTF-8")
    } else cat("\n[3.1 Figure 1] 缺少原始CSV，跳过。\n")
  }
  if(RUN_FIGURE2) {
    if(!interactive()) stop("Figure 2需要桌面RStudio选择文件与确认保存")
    cat("\n[3.1 Figure 2] 即将选择GBD人口CSV和死亡CSV...\n")
    source(file.path(repo_root,"03_1_descriptive_statistics","figure2","run.R"),
           local=new.env(parent=globalenv()),encoding="UTF-8")
  }

  # ----------------------- STAGE 4.1 FIGURES ----------------------
  # These independent visualization scripts may have a distinct factor
  # identification convention. Their kappa MUST NOT replace the standardized
  # selected factor used in §4.2 automatically.
  if(RUN_FIGURE3) {
    if(!interactive()) stop("Figure 3需要交互选择Population和Deaths原始CSV")
    source(file.path(repo_root,"04_results","4_1_single_population",
                    "run_4_1_poisson_lc.R"),
           local=new.env(parent=globalenv()),encoding="UTF-8")
  }
  if(RUN_FIGURE4) {
    if(!interactive()) stop("Figure 4需要交互选择Figure3输出的fitted_rates.csv")
    source(file.path(repo_root,"04_results","4_1_single_population",
                    "run_4_1_poisson_lc_2factor.R"),
           local=new.env(parent=globalenv()),encoding="UTF-8")
  }

  # ----------------------- STAGE 4.1 -> 4.2 ----------------------
  if(REFIT_ALL_41_MODELS) {
    if(!nzchar(FITTED_RATES_CSV)) {
      if(!interactive())stop("请设置FITTED_RATES_CSV为本地fitted_rates.csv路径")
      cat("\n[4.1] 选择原始fitted_rates.csv（内有死亡数与暴露量）...\n")
      FITTED_RATES_CSV <- file.choose()
    }
    section42_input <- normalizePath(FITTED_RATES_CSV,mustWork=TRUE)
    cat("\n[4.1] 即将从原始D/E重估所有126组单人口模型。\n")
  } else {
    section42_input <- file.path(repo_root,"04_results",
        "4_1_single_population","results","kappa_all_models.csv")
    if(!file.exists(section42_input))
      stop("找不到§4.1缓存时期因子：",section42_input,
           "\n请重新下载完整仓库，或设置REFIT_ALL_41_MODELS=TRUE。")
    cat("\n[4.1] 使用已核验的126组单人口模型时期因子（不重新拟合）。\n")
  }
  output_base <- file.path(repo_root,"04_results","4_2_clusters",
                           "gbd_complete","outputs")
  dir.create(output_base,recursive=TRUE,showWarnings=FALSE)
  out_42 <- file.path(output_base,paste0("Section42_",
                format(Sys.time(),"%Y%m%d_%H%M%S")))
  cat("\n[4.2] 计算DTW、层次聚类、PAM、七项CVI和所有候选K...\n")
  run_section42(input=section42_input,out=out_42)
  cat("\n=====================================================\n")
  cat("§4.2 候选分析已运行完成\n")
  cat("输出文件夹：",normalizePath(out_42,winslash="/"),"\n")
  cat("CVI表：cluster_validity_scores.csv（196行）\n")
  cat("PAM起点：pam_three_start_check.csv（294行）\n")
  cat("树状图和CVI图：figures/\n")
  cat("选择K：填写 cluster_selection_template.csv（14行），另存 reviewed_selection.csv\n")
  cat("然后在RStudio Console运行：\n")
  cat("finalize_section42(out_42, file.path(out_42,'reviewed_selection.csv'))\n")
  cat("原文没有公开把7项CVI合成唯一K的公式；不要强制填原文5/6组。\n")
  cat("=====================================================\n")
}
