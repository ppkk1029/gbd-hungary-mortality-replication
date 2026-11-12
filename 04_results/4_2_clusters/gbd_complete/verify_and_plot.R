# 独立检查：梯度有限差分、枚举DTW路径、SVD/旋转保持拟合面、输出数量。
algorithm_selftest <- function() {
 x<-c(0,1,3);y<-c(0,2,3);cost<-outer(x,y,'-')^2;vals<-numeric()
 rec<-function(i,j,v) {
  if(i==3L && j==3L) {vals<<-c(vals,v);return(invisible(NULL))}
  for(step in list(c(1,1,2),c(1,0,1),c(0,1,1))) {
   a<-i+step[1];b<-j+step[2];if(a<=3 && b<=3) rec(a,b,v+step[3]*cost[a,b])
  }
 }
 rec(1,1,cost[1,1]);err<-abs(dtw_distance(x,y)-sqrt(min(vals)/6))
 stopifnot(err<1e-12,abs(dtw_distance(x,x))<1e-12,
           abs(dtw_distance(x,y,window=0)-sqrt(2*sum((x-y)^2)/6))<1e-12)
 # 使用构造的计数/暴露量检查七个目标函数的解析梯度，不作实证结果。
 E<-matrix(10000,6,36);D<-E*exp(outer(seq(-4,-2,length.out=6),seq(-.2,.2,length.out=36),'+'))
 ge<-sapply(MODELS,function(model) {
  p<-make_problem(D,E,model);v<-p$start+sin(seq_along(p$start))*.001;h<-1e-5
  numeric_g<-vapply(seq_along(v),function(j) {hi<-lo<-v;hi[j]<-hi[j]+h;lo[j]<-lo[j]-h;(p$calc(hi)$value-p$calc(lo)$value)/(2*h)},0.0)
  max(abs(numeric_g-p$calc(v)$gradient))
 })
 stopifnot(max(ge)<1e-6)
 data.frame(check=c('DTW_bruteforce_error',paste0('gradient_',MODELS)),value=c(err,ge))
}
verify_results <- function(fits,analyses,out,reference=NULL) {
 f<-read.csv(file.path(out,'fit_diagnostics_all_models.csv'))
 stopifnot(nrow(f)==126L,all(f$accepted),length(analyses)==14L,
  nrow(read.csv(file.path(out,'kappa_all_models.csv')))==8424L,
  nrow(read.csv(file.path(out,'standardized_kappa.csv')))==4536L,
  nrow(read.csv(file.path(out,'pairwise_distances.csv')))==504L,
  nrow(read.csv(file.path(out,'cluster_validity_scores.csv')))==196L,
  nrow(read.csv(file.path(out,'candidate_membership.csv')))==1764L,
  nrow(read.csv(file.path(out,'pam_three_start_check.csv')))==294L)
 if(file.exists(file.path(out,'cluster_membership.csv'))) {
  m <- read.csv(file.path(out,'cluster_membership.csv'))
  stopifnot(nrow(m)==126L)
 } else {
  t <- read.csv(file.path(out,'cluster_selection_template.csv'))
  stopifnot(nrow(t)==14L)
 }
 rotation<-0;th<-pi/12;R<-matrix(c(cos(th),sin(th),-sin(th),cos(th)),2,2)
 for(nm in names(fits)) if(endsWith(nm,'LC2')) {a<-fits[[nm]]
  rotation<-max(rotation,abs(a$B%*%t(a$K)-(a$B%*%R)%*%t(a$K%*%R)))}
 stopifnot(rotation<1e-10)
 v<-data.frame(check=c('fits_passed','max_gradient','max_multistart_eta_difference','rotation_surface_max_error'),
  value=c(nrow(f),max(f$max_abs_scaled_gradient),max(f$max_eta_difference_3_starts),rotation))
 write_csv(rbind(algorithm_selftest(),v),file.path(out,'validation.csv'))
 if(!is.null(reference)) compare_reference(out,reference)
 invisible(v)
}
compare_reference <- function(out,reference) {
 # reference是迁移前的只读数值基准；只用于核验，绝不参与参数估计。
 comparisons<-list()
 for(item in list(c('fit_diagnostics_all_models','sex','cause','model','deviance'),
                  c('kappa_all_models','sex','cause','model','factor','year','kappa'),
                  c('standardized_kappa','sex','model','cause','year','z'),
                  c('pairwise_distances','sex','model','cause1','cause2','dtw_p2'))) {
  file<-paste0(item[1],'.csv');value<-tail(item,1);by<-item[2:(length(item)-1)]
  a<-read.csv(file.path(out,file));b<-read.csv(file.path(reference,file));m<-merge(a,b,by=by,suffixes=c('.R','.reference'))
  stopifnot(nrow(m)==nrow(a),nrow(m)==nrow(b))
  delta<-abs(m[[paste0(value,'.R')]]-m[[paste0(value,'.reference')]])
  comparisons[[item[1]]]<-data.frame(file=file,quantity=value,max_abs_difference=max(delta))
 }
 a<-read.csv(file.path(out,'cluster_membership.csv'));b<-read.csv(file.path(reference,'cluster_membership.csv'))
 m<-merge(a,b,by=c('sex','model','cause'))
 # 直接按模型性别切片计算ARI，避免簇号交换被误判为不同。
 aris<-do.call(rbind,lapply(split(m,paste(m$sex,m$model)),function(d)
  data.frame(sex=d$sex[1],model=d$model[1],ARI=ari(d$cluster.x,d$cluster.y))))
 write_csv(do.call(rbind,comparisons),file.path(out,'migration_numeric_comparison.csv'))
 write_csv(aris,file.path(out,'migration_membership_comparison.csv'))
 # Changes in K are expected after replacing the Silhouette-only decision rule.
 invisible(aris)
}
