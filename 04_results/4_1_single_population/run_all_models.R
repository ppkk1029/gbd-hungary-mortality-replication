# §4.1 pure-R seven-model single-population fit driver; source after model_functions.R.
# 由 run_4_2_complete.R 调用；按阶段阅读可从 fit_all_models() 开始。
fit_all_models <- function(input,out) {
 raw<-read_input(input);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 fits<-list();diags<-list();kp<-list();bp<-list();n<-0L
 for(sex in c('Male','Female')) for(cause in CAUSES) {
  d<-input_matrices(raw,sex,cause)
  for(model in MODELS) {
   f<-fit_model(d$D,d$E,model);nm<-key(sex,cause,model);fits[[nm]]<-f;n<-n+1L
   diags[[n]]<-cbind(data.frame(sex=sex,cause=cause,model=model),f$diag)
   for(j in seq_len(ncol(f$K))) {
    kp[[length(kp)+1L]]<-data.frame(sex=sex,cause=cause,model=model,factor=j,year=YEARS,kappa=f$K[,j])
    bp[[length(bp)+1L]]<-data.frame(sex=sex,cause=cause,model=model,factor=j,age=AGES,alpha=f$a,beta=f$B[,j])
   }
  }
  message('已完成拟合 ',sex,' / ',cause,' (',n,'/126)')
 }
 write_csv(do.call(rbind,diags),file.path(out,'fit_diagnostics_all_models.csv'))
 write_csv(do.call(rbind,kp),file.path(out,'kappa_all_models.csv'))
 write_csv(do.call(rbind,bp),file.path(out,'parameters_all_models.csv'))
 saveRDS(fits,file.path(out,'fitted_models.rds'))
 fits
}
