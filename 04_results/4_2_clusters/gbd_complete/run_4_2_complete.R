# §4.2完整纯R入口。RStudio打开本文件点Source，选择 fitted_rates.csv 即可。
# 命令行：Rscript run_4_2_complete.R /path/fitted_rates.csv /path/new_results
# 输入CSV只读取 deaths_gbd / exposure_proxy，不使用 fitted_count。
.section42_script_dir <- local({
 files <- lapply(sys.frames(),function(f) f$ofile)
 files <- Filter(function(x)!is.null(x)&&length(x)==1L,files)
 if(length(files)) dirname(normalizePath(tail(files,1)[[1]],mustWork=TRUE)) else {
  a <- grep('^--file=',commandArgs(),value=TRUE)
  if(length(a)) dirname(normalizePath(sub('^--file=','',a[1]),mustWork=TRUE)) else getwd()
 }
})
.section41_script_dir <- file.path(dirname(dirname(.section42_script_dir)), "4_1_single_population")
for(.f in c('model_functions.R','run_all_models.R')) {
 source(file.path(.section41_script_dir,.f),local=TRUE,encoding='UTF-8')
}
for(.f in c('clustering_functions.R','cluster_all_models.R',
             'complete_cvi.R','verify_and_plot.R','plot_paper_figures.R')) {
 source(file.path(.section42_script_dir,.f),local=TRUE,encoding='UTF-8')
}
run_section42 <- function(input=NULL,out=NULL,reference=NULL) {
 if(is.null(input)) {
  if(!interactive()) stop('命令行必须提供输入CSV路径')
  message('选择原始 fitted_rates.csv，需含 deaths_gbd 与 exposure_proxy。')
  input<-file.choose()
 }
 input<-normalizePath(input,mustWork=TRUE)
 if(is.null(out)) out<-file.path(dirname(input),paste0('Section_4_2_R_',format(Sys.time(),'%Y%m%d_%H%M%S')))
 if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop('输出目录非空，请指定一个新的目录，避免覆盖之前结果：',out)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);out<-normalizePath(out,mustWork=TRUE)
 message('R版开始。输入：',input,'；输出：',out)
 checks<-algorithm_selftest();write_csv(checks,file.path(out,'algorithm_selftest.csv'))
 fits<-fit_all_models(input,out)
 analyses<-cluster_all_models(fits,out)
 complete_cvi(analyses,out)
 verify_results(fits,analyses,out,reference)
 plot_paper_figures(analyses,out,fits)
 writeLines(c('engine: R (no external-language dependency)',paste('R:',getRversion()),
  paste('input_md5:',unname(tools::md5sum(input))),paste('seed:',SEED),
  'years: 1980-2015; ages: 62,67,72,77,82,87; 9 causes; Male/Female',
  'DTW: squared local cost, symmetric2, sqrt(acc/(n+m)), closed endpoints; no main window',
  'LC: centered K, SVD orthogonal B, descending singular values, largest loading positive',
  'Binomial: N=exposure_proxy+deaths_gbd/2 (approximation)',
  'selection: maximum silhouette; tie lower K then alphabetical method; K=2..8',
  'GBD estimated counts: working likelihood, not HCSO registered counts'),file.path(out,'manifest.txt'))
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
 needed<-c('Figure_5_GBD_Male_Poisson_LC2.png','Figure_6_GBD_Female_Poisson_LC2.png',
           'Figure_7_GBD_Male_CBD_Clusters.png','Figure_8_GBD_Male_Binomial_LC2.png')
 stopifnot(all(file.exists(file.path(out,needed))))
 print(read.csv(file.path(out,'cluster_summary.csv'))[,c('sex','model','K','silhouette')],row.names=FALSE)
 message('全部完成：',out)
 if(interactive()) plot_cluster_panels(analyses[['Male Binomial LC2']],'Figure 8 | Male Binomial LC2 | GBD')
 invisible(out)
}
# source() 默认立即运行；开发/分步学习可先 options(section42.autorun=FALSE)。
if(isTRUE(getOption('section42.autorun',TRUE))) {
 .args<-commandArgs(trailingOnly=TRUE)
 run_section42(if(length(.args)>=1L) .args[1] else NULL,
               if(length(.args)>=2L) .args[2] else NULL,
               if(length(.args)>=3L) .args[3] else NULL)
}
