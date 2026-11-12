# Section 4.2 unique official entry point. Default: reuse §4.1's verified
# period indices, without rerunning all 126 fits. Alternatively pass raw GBD
# fitted_rates.csv and refit. DO NOT use historical Silhouette-only results.
#
# Interactive:
#   options(section42.autorun=FALSE)
#   source("run_4_2_complete.R",encoding="UTF-8")
#   out <- run_section42()    # uses §4.1/results/kappa_all_models.csv
#   # Review all CVIs, dendrograms and membership candidates; then fill CSV
#   finalize_section42(out,file.path(out,"reviewed_selection.csv"))
#
# CLI:
#   Rscript run_4_2_complete.R /path/kappa_all_models.csv /path/output
#   Rscript run_4_2_complete.R /path/fitted_rates.csv /path/output
#   Rscript run_4_2_complete.R --finalize /path/output /path/reviewed_selection.csv
.section42_script_dir <- local({
  files <- lapply(sys.frames(),function(f) f$ofile)
  files <- Filter(function(x)!is.null(x)&&length(x)==1L,files)
  if(length(files)) dirname(normalizePath(tail(files,1L)[[1L]],mustWork=TRUE)) else {
    a <- grep("^--file=",commandArgs(),value=TRUE)
    if(length(a)) dirname(normalizePath(sub("^--file=","",a[1L]),mustWork=TRUE))
    else getwd()
  }
})
.section41_script_dir <- file.path(dirname(dirname(.section42_script_dir)),
                                   "4_1_single_population")
for(.f in c("model_functions.R","run_all_models.R"))
  source(file.path(.section41_script_dir,.f),local=TRUE,encoding="UTF-8")
for(.f in c("kappa_import.R","clustering_functions.R","cluster_all_models.R",
             "complete_cvi.R","verify_and_plot.R","plot_paper_figures.R"))
  source(file.path(.section42_script_dir,.f),local=TRUE,encoding="UTF-8")

finalize_section42 <- function(output_dir,selection_file) {
  check_cluster_packages()
  output_dir <- normalizePath(output_dir,mustWork=TRUE)
  selection_file <- normalizePath(selection_file,mustWork=TRUE)
  path <- file.path(output_dir,"cluster_analyses.rds")
  if(!file.exists(path)) stop("必须先运行run_section42()生成全部候选")
  if(!file.exists(file.path(output_dir,"cluster_validity_scores.csv")))
    stop("七项聚类评价指标未生成，不能确认K")
  analyses <- readRDS(path)
  analyses <- finalize_cluster_selections(analyses,output_dir,selection_file)
  plot_paper_figures(analyses,output_dir)
  summary <- read.csv(file.path(output_dir,"cluster_summary.csv"))
  membership <- read.csv(file.path(output_dir,"cluster_membership.csv"))
  if(nrow(summary)!=14L || nrow(membership)!=126L ||
     any(summary$selection_status!="researcher_reviewed"))
    stop("最终K及成员归属不完整")
  writeLines(c("Section 4.2 cluster decisions manually reviewed",
      paste("selection_file:",selection_file),
      paste("selection_md5:",unname(tools::md5sum(selection_file))),
      "No inferred Varga (2025) seven-CVI selection rule",
      paste("R_version:",as.character(getRversion()))),
      file.path(output_dir,"selection_manifest.txt"))
  message("正式分组已生成：",output_dir)
  invisible(output_dir)
}

run_section42 <- function(input=NULL,out=NULL,selection_file=NULL) {
  check_cluster_packages()
  if(is.null(input)) {
    cache <- file.path(.section41_script_dir,"results","kappa_all_models.csv")
    if(file.exists(cache)) input <- cache
    else if(interactive()) {
      message("选择 §4.1 kappa_all_models.csv 或原始 fitted_rates.csv")
      input <- file.choose()
    } else stop("缺少输入：请传入 §4.1 时期因子或原始GBD数据")
  }
  input <- normalizePath(input,mustWork=TRUE)
  headers <- names(read.csv(input,nrows=0,check.names=FALSE))
  from_kappa <- all(c("sex","cause","model","factor","year","kappa") %in% headers)
  from_raw <- all(c("sex_name","cause_name","age_name","year",
                    "deaths_gbd","exposure_proxy") %in% headers)
  if(from_kappa==from_raw)
    stop("CSV格式不正确：必须是§4.1时期因子，或完整的GBD死亡/暴露量表")
  if(is.null(out))out <- file.path(.section42_script_dir,"outputs",
      paste0("Section42_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
    stop("输出目录非空，拒绝覆盖：",out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  out <- normalizePath(out,mustWork=TRUE)
  message("§4.2输入：",input)
  message("§4.2输出：",out)
  if(from_kappa) {
    fits <- read_single_population_kappa(input)
    copy_single_population_outputs(input,out)
    source_kind <- "Existing independently fitted §4.1 kappa; NO refitting"
    if(!file.exists(file.path(out,"parameters_all_models.csv")))
      message("警告：缺少§4.1参数CSV，此输出不可直接用于§4.3")
  } else {
    fits <- fit_all_models(input,out)
    source_kind <- "Refitted all 126 single-population models from raw GBD D/E"
    checks <- algorithm_selftest()
    write_csv(checks,file.path(out,"algorithm_selftest.csv"))
  }
  analyses <- cluster_all_models(fits,out,selection_file)
  complete_cvi(analyses,out)
  verify_cluster_stage(analyses,out)
  if(!from_kappa)verify_results(fits,analyses,out)
  plot_paper_figures(analyses,out)
  writeLines(c(
    "Section 4.2 / Varga (2025) methodological analogue, not HCSO exact data",
    paste("source:",source_kind),
    paste("input:",input),
    paste("input_md5:",unname(tools::md5sum(input))),
    paste("dtw_version:",as.character(utils::packageVersion("dtw"))),
    paste("dtwclust_version:",as.character(utils::packageVersion("dtwclust"))),
    "Nine GBD level-2 causes, 1980-2015, ages 60-89 in 5-year groups",
    "LC2: second factor; other models: first factor",
    "DTW p=2 squared local cost, symmetric2, sqrt(raw/(n+m)), no window",
    "clustering: hierarchical average + PAM with 3 reproducible starts",
    "candidate K: 2..8; original seven dtwclust CVIs",
    "K decision: human review of dendrogram/CVIs/starts; not auto-voting",
    "No original HCSO/ICD input in repository; cannot claim equal Table 3/4",
    "LC2 factor orientation follows the separately fitted §4.1 SVD convention"
  ),file.path(out,"manifest.txt"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  message("§4.2候选已计算，需审查 cluster_selection_template.csv 后 finalize_section42()：",out)
  invisible(out)
}

if(isTRUE(getOption("section42.autorun",TRUE))) {
  a <- commandArgs(trailingOnly=TRUE)
  if(length(a)>=1L&&identical(a[1L],"--finalize")) {
    if(length(a)<3L)stop("请传 --finalize 输出目录 reviewed_selection.csv")
    finalize_section42(a[2L],a[3L])
  } else run_section42(input=if(length(a)>=1L)a[1L] else NULL,
            out=if(length(a)>=2L)a[2L] else NULL,
            selection_file=if(length(a)>=3L)a[3L] else NULL)
}
