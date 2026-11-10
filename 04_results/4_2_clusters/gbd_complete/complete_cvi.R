# 七项CVI已经由choose_clusters计算。这里重算并核对，另记录294个PAM随机起点。
complete_cvi <- function(analyses,out) {
 rows<-checks<-list();set.seed(SEED)
 for(sex in c('Male','Female')) for(model in MODELS) {
  a<-analyses[[paste(sex,model)]];b<-choose_clusters(a$d,a$Z)
  stopifnot(ari(a$labels,b$labels)==1)
  rows[[paste(sex,model)]]<-cbind(data.frame(sex=sex,model=model),b$scores)
  for(k in 2:8) checks[[length(checks)+1L]]<-cbind(data.frame(sex=sex,model=model),pam_check(a$d,k))
 }
 tab<-do.call(rbind,rows);check<-do.call(rbind,checks)
 stopifnot(nrow(tab)==196L,nrow(check)==294L,min(check$gap)> -1e-10)
 write_csv(tab,file.path(out,'cluster_validity_scores.csv'))
 write_csv(check,file.path(out,'pam_three_start_check.csv'))
 invisible(check)
}
