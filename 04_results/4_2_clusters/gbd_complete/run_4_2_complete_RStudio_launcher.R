# 旧文件名兼容入口：现在完全运行R代码。推荐直接打开 run_4_2_complete.R 点Source。
local({
 files<-Filter(function(x)!is.null(x)&&length(x)==1L,lapply(sys.frames(),function(f)f$ofile))
 if(length(files)) here<-dirname(normalizePath(tail(files,1)[[1]])) else {
  a<-grep('^--file=',commandArgs(),value=TRUE)
  here<-if(length(a)) dirname(normalizePath(sub('^--file=','',a[1]))) else getwd()
 }
 source(file.path(here,'run_4_2_complete.R'),encoding='UTF-8')
})
