# 兼容旧入口名；已升级为七模型完整R流程。
# 输入由旧版kappa.csv改为含真实deaths_gbd/exposure_proxy的fitted_rates.csv。
local({
 files<-Filter(function(x)!is.null(x)&&length(x)==1L,lapply(sys.frames(),function(f)f$ofile))
 if(length(files)) here<-dirname(normalizePath(tail(files,1)[[1]])) else {
  a<-grep('^--file=',commandArgs(),value=TRUE)
  here<-if(length(a)) dirname(normalizePath(sub('^--file=','',a[1]))) else getwd()
 }
 source(file.path(here,'gbd_complete','run_4_2_complete.R'),encoding='UTF-8')
})
