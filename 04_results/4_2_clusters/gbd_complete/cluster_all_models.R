# Section 4.2: seven model families x two sexes, all candidate partitions.
# NO algorithm is claimed to infer the author's final K: the paper does not
# publish a deterministic rule combining the seven CVIs and dendrograms.

finalize_cluster_selections <- function(analyses,out,selection_file) {
  if(is.null(selection_file)) stop("选K文件缺失；先查看七项指标与树状图")
  manual <- read_selection_file(selection_file,MODELS)
  if(length(analyses)!=14L) stop("应当提供14个性别×模型组")
  sum_rows <- member_rows <- decision_rows <- list()
  for (nm in names(analyses)) {
    a <- analyses[[nm]]
    sex <- sub(" .*","",nm)
    model <- substring(nm,nchar(sex)+2L)
    x <- manual[manual$sex==sex & manual$model==model,]
    if(nrow(x)!=1L) stop("缺少审查记录: ",nm)
    method <- as.character(x$method[1L])
    k <- as.integer(x$K[1L])
    reason <- as.character(x$reason[1L])
    labs <- a$candidate_labels[[paste(method,k,sep="::")]]
    if(is.null(labs)) stop("候选分组不存在: ",nm,", method=",method,", K=",k)
    selected <- a$scores[a$scores$method==method & a$scores$K==k,,drop=FALSE]
    if(nrow(selected)!=1L) stop("所选K不唯一：",nm)
    if(length(labs)!=length(CAUSES)||length(unique(labs))!=k)
      stop("标签数量与K不一致：",nm)
    a$labels <- as.integer(labs)
    a$selection <- data.frame(method=method,K=k,
       selection_status="researcher_reviewed",reason=reason)
    analyses[[nm]] <- a
    sum_rows[[nm]] <- cbind(
       data.frame(sex=sex,model=model,factor=a$factor),
       selected[,c("method","K","Sil","D","COP","DB","DBstar","CH","SF",
                   "singleton_count","objective")],
       data.frame(selection_status="researcher_reviewed",reason=reason))
    member_rows[[nm]] <- data.frame(
       sex=sex,model=model,cause=CAUSES,code=SHORT,cluster=as.integer(labs))
    decision_rows[[nm]] <- data.frame(
       sex=sex,model=model,factor=a$factor,method=method,K=k,
       selection_status="researcher_reviewed",reason=reason)
  }
  write_csv(do.call(rbind,sum_rows),file.path(out,"cluster_summary.csv"))
  write_csv(do.call(rbind,member_rows),file.path(out,"cluster_membership.csv"))
  write_csv(do.call(rbind,decision_rows),file.path(out,"cluster_selection_decisions.csv"))
  saveRDS(analyses,file.path(out,"cluster_analyses.rds"))
  analyses
}

cluster_all_models <- function(fits,out,selection_file=NULL) {
  memberships <- scores <- starts <- pairs <- zrows <- analyses <- templates <- list()
  for (sex in c("Male","Female")) for (model in MODELS) {
    factor <- if(endsWith(model,"LC2")) 2L else 1L
    X <- do.call(rbind,lapply(CAUSES,function(cause) {
      fit <- fits[[key(sex,cause,model)]]
      if(is.null(fit) || is.null(fit$K) || ncol(fit$K)<factor)
        stop("缺少单人口时期因子: ",sex," / ",model," / ",cause)
      fit$K[,factor]
    }))
    rownames(X) <- SHORT
    Z <- z_rows(X)
    d <- distance_matrix(Z,p=2,check_impl=TRUE)
    eva <- evaluate_clusters(Z,d,seed=SEED)
    nm <- paste(sex,model)
    analyses[[nm]] <- list(
       Z=Z,raw_kappa=X,d=d,hc=eva$hc,labels=NULL,
       factor=factor,scores=eva$scores,
       candidate_labels=eva$labels)
    scores[[nm]] <- cbind(
       data.frame(sex=sex,model=model,factor=factor),eva$scores)
    if(nrow(eva$starts))
      starts[[nm]] <- cbind(data.frame(sex=sex,model=model),eva$starts)
    for(k in 2:8) for(method in c("hierarchical_average","pam")) {
      labs <- eva$labels[[paste(method,k,sep="::")]]
      memberships[[length(memberships)+1L]] <-
        data.frame(sex=sex,model=model,method=method,K=k,
                   cause=CAUSES,code=SHORT,cluster=labs)
    }
    templates[[nm]] <- data.frame(
       sex=sex,model=model,factor=factor,
       K=NA_integer_,method="",reason="")
    for(i in seq_along(CAUSES)) {
      zrows[[length(zrows)+1L]] <- data.frame(
        sex=sex,model=model,cause=CAUSES[i],year=YEARS,z=Z[i,])
      if(i>1L) for(j in seq_len(i-1L))
        pairs[[length(pairs)+1L]] <- data.frame(
           sex=sex,model=model,cause1=CAUSES[i],
           cause2=CAUSES[j],dtw_p2=d[i,j])
    }
    message("已计算所有候选: ",sex," / ",model,
            "; 七项CVI、两种方法、K=2..8")
  }
  write_csv(do.call(rbind,scores),file.path(out,"cluster_validity_scores.csv"))
  write_csv(do.call(rbind,memberships),file.path(out,"candidate_membership.csv"))
  write_csv(do.call(rbind,templates),file.path(out,"cluster_selection_template.csv"))
  write_csv(do.call(rbind,pairs),file.path(out,"pairwise_distances.csv"))
  write_csv(do.call(rbind,zrows),file.path(out,"standardized_kappa.csv"))
  write_csv(do.call(rbind,starts),file.path(out,"pam_three_start_check.csv"))
  saveRDS(analyses,file.path(out,"cluster_analyses.rds"))
  if(!is.null(selection_file))
    analyses <- finalize_cluster_selections(analyses,out,selection_file)
  analyses
}
