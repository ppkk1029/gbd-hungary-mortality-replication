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
cluster_all_models <- function(fits,out) {
 summary<-membership<-scores<-sens<-pairs<-zrows<-list();analyses<-list()
 for(sex in c('Male','Female')) for(model in MODELS) {
  factor<-if(endsWith(model,'LC2')) 2L else 1L
  X<-do.call(rbind,lapply(CAUSES,function(c) fits[[key(sex,c,model)]]$K[,factor]))
  Z<-z_rows(X);d<-distance_matrix(Z);a<-choose_clusters(d,Z);lab<-a$labels
  nm<-paste(sex,model);analyses[[nm]]<-list(Z=Z,d=d,hc=a$hc,labels=lab)
  summary[[nm]]<-cbind(data.frame(sex=sex,model=model,factor=factor),a$scores[1,])
  membership[[nm]]<-data.frame(sex=sex,model=model,cause=CAUSES,code=SHORT,cluster=lab)
  scores[[nm]]<-cbind(data.frame(sex=sex,model=model),a$scores)
  for(i in seq_along(CAUSES)) {
   zrows[[length(zrows)+1L]]<-data.frame(sex=sex,model=model,cause=CAUSES[i],year=YEARS,z=Z[i,])
   if(i>1) for(j in seq_len(i-1)) pairs[[length(pairs)+1L]]<-data.frame(sex=sex,model=model,cause1=CAUSES[i],cause2=CAUSES[j],dtw_p2=d[i,j])
  }
  record_sens<-function(tag,ds,Zs,reference) {
   b<-choose_clusters(ds,Zs)
   data.frame(sex=sex,model=model,variant=tag,K=b$scores$K[1],silhouette=b$scores$silhouette[1],ARI=ari(reference,b$labels))
  }
  sens[[length(sens)+1L]]<-record_sens('window5',distance_matrix(Z,window=5),Z,lab)
  sens[[length(sens)+1L]]<-record_sens('L1_legacy',distance_matrix(Z,p=1),Z,lab)
  if(sex=='Female') {keep<-setdiff(seq_along(CAUSES),7L)
   sens[[length(sens)+1L]]<-record_sens('without_DKD',d[keep,keep],Z[keep,],lab[keep])}
  if(factor==2L) for(angle in c(15,30,45)) {
   th<-angle*pi/180;R<-matrix(c(cos(th),sin(th),-sin(th),cos(th)),2,2)
   XR<-do.call(rbind,lapply(CAUSES,function(c) (fits[[key(sex,c,model)]]$K%*%R)[,2]))
   ZR<-z_rows(XR);sens[[length(sens)+1L]]<-record_sens(paste0('rotation_',angle),distance_matrix(ZR),ZR,lab)
  }
  if(model=='Poisson LC2') {
   j<-if(sex=='Male') 3L else 7L;p<-dtw_distance(Z[1,],Z[j,],path=TRUE)$path
   write_csv(data.frame(year_1=YEARS[p[,1]],year_2=YEARS[p[,2]]),file.path(out,paste0('alignment_path_',sex,'.csv')))
  }
  message('已完成聚类 ',sex,' / ',model)
 }
 outputs<-list(cluster_summary=summary,cluster_membership=membership,cluster_validity_scores=scores,
               sensitivity=sens,pairwise_distances=pairs,standardized_kappa=zrows)
 for(nm in names(outputs)) write_csv(do.call(rbind,outputs[[nm]]),file.path(out,paste0(nm,'.csv')))
 saveRDS(analyses,file.path(out,'cluster_analyses.rds'))
 analyses
}
