# §4.2 cluster driver; source after §4.1 fit functions and clustering_functions.R.
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
