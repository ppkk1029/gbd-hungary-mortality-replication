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
verify_results <- function(fits,analyses,out) {
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
 invisible(v)
}

# Works for BOTH raw-GBD refitting and already-fitted §4.1 period indices.
verify_cluster_stage <- function(analyses,out) {
  if(length(analyses)!=14L)stop("应完成2性别×7模型共14组")
  files <- c("cluster_validity_scores.csv","candidate_membership.csv",
             "pam_three_start_check.csv","pairwise_distances.csv",
             "standardized_kappa.csv","cluster_selection_template.csv")
  absent <- files[!file.exists(file.path(out,files))]
  if(length(absent))stop("缺失聚类文件：",paste(absent,collapse=","))
  scored <- read.csv(file.path(out,files[1]),stringsAsFactors=FALSE)
  members <- read.csv(file.path(out,files[2]),stringsAsFactors=FALSE)
  starts <- read.csv(file.path(out,files[3]),stringsAsFactors=FALSE)
  distances <- read.csv(file.path(out,files[4]),stringsAsFactors=FALSE)
  std <- read.csv(file.path(out,files[5]),stringsAsFactors=FALSE)
  template <- read.csv(file.path(out,files[6]),stringsAsFactors=FALSE)
  stopifnot(nrow(scored)==196L,nrow(members)==1764L,
    nrow(starts)==294L,nrow(distances)==504L,
    nrow(std)==4536L,nrow(template)==14L)
  if(anyDuplicated(scored[c("sex","model","method","K")]))
    stop("CVI出现重复模型/性别/方法/K")
  if(any(!is.finite(distances$dtw_p2))||
     any(distances$dtw_p2 < 0))stop("DTW距离包含无效数值")
  if(any(!is.finite(starts$objective)) ||
     any(starts$ARI_vs_best < -1 | starts$ARI_vs_best > 1))
    stop("PAM三起点记录错误")
  if(!setequal(names(CVI_DIRECTION),CVI_NAMES))
    stop("七指标方向与名称不一致")
  if(any(!CVI_NAMES %in% names(scored)))
    stop("CVI指标缺失")
  for(key in unique(paste(scored$sex,scored$model,sep="::"))) {
    v <- scored[paste(scored$sex,scored$model,sep="::")==key,]
    for(metric in CVI_NAMES) {
      if(!any(is.finite(v[[metric]])))
        stop("以下模型七项指标之一完全无法计算：",key," ",metric,
             "；请检查dtwclust参数和原始数据")
    }
  }
  if(file.exists(file.path(out,"cluster_membership.csv"))) {
    m <- read.csv(file.path(out,"cluster_membership.csv"))
    dec <- read.csv(file.path(out,"cluster_selection_decisions.csv"))
    if(nrow(m)!=126L||nrow(dec)!=14L||
       any(dec$selection_status!="researcher_reviewed"))
      stop("最终分组缺少审查记录")
  }
  invisible(TRUE)
}
