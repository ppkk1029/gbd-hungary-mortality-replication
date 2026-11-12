# §4.2 driver: 14 sex/model analyses, seven official CVIs, dendrogram, 3 PAM starts.
# Final cluster choice requires documented researcher review by default.
finalize_cluster_selections <- function(analyses, out, selection_file=NULL,
                                      selection_mode="review_required") {
  if(!selection_mode %in% c("review_required","exploratory_consensus"))
    stop("selection_mode必须是 review_required 或 exploratory_consensus")
  manual <- read_selection_file(selection_file,MODELS)
  if (is.null(manual) && selection_mode=="review_required") {
    message("已保存全部聚类候选及建议。论文未公开唯一选K公式，等待核对树状图/CVI后填写 cluster_selection_template.csv。")
    return(analyses)
  }
  sum_rows <- member_rows <- decision_rows <- list()
  for (nm in names(analyses)) {
    a <- analyses[[nm]]
    bits <- strsplit(nm," ",fixed=TRUE)[[1L]]
    sex <- bits[1L]; model <- substring(nm,nchar(sex)+2L)
    if (!is.null(manual)) {
      x <- manual[manual$sex==sex & manual$model==model,]
      method <- as.character(x$method[1L]); k <- as.integer(x$K[1L])
      status <- "researcher_reviewed"; reason <- if("reason" %in% names(x)) as.character(x$reason[1L]) else ""
    } else {
      method <- a$proposal$method[1L]; k <- a$proposal$K[1L]
      status <- "exploratory_equal_rank_NOT_paper_rule"
      reason <- "Two-method mean ranks of seven CVIs; NOT reported as author's decision rule"
    }
    label <- a$candidate_labels[[paste(method,k,sep="::")]]
    if(is.null(label)) stop("所选分组不存在：",nm,", ",method,", K=",k)
    a$labels <- as.integer(label)
    a$selection <- data.frame(method=method,K=k,selection_status=status,reason=reason)
    analyses[[nm]] <- a
    selected <- a$scores[a$scores$method==method & a$scores$K==k,,drop=FALSE]
    if(nrow(selected)!=1L) stop("不能唯一匹配K")
    sum_rows[[nm]] <- cbind(data.frame(sex=sex,model=model,factor=a$factor),
        selected[,c("method","K","Sil","D","COP","DB","DBstar","CH","SF",
                    "singleton_count","objective","mean_indicator_rank")],
        data.frame(selection_status=status,reason=reason))
    member_rows[[nm]] <- data.frame(sex=sex,model=model,
       cause=CAUSES,code=SHORT,cluster=as.integer(label))
    decision_rows[[nm]] <- data.frame(sex=sex,model=model,factor=a$factor,
       method=method,K=k,selection_status=status,reason=reason,
       suggested_K=a$proposal$K[1L],suggested_method=a$proposal$method[1L])
  }
  write_csv(do.call(rbind,sum_rows),file.path(out,"cluster_summary.csv"))
  write_csv(do.call(rbind,member_rows),file.path(out,"cluster_membership.csv"))
  write_csv(do.call(rbind,decision_rows),file.path(out,"cluster_selection_decisions.csv"))
  saveRDS(analyses,file.path(out,"cluster_analyses.rds"))
  analyses
}
cluster_all_models <- function(fits,out,selection_file=NULL,
                               selection_mode="review_required") {
  summaries <- list(); memberships <- list(); scores <- list()
  starts <- list(); pairs <- list(); zrows <- list(); analyses <- list()
  templates <- list(); comparisons <- list()
  for (sex in c("Male","Female")) for (model in MODELS) {
    factor <- if(endsWith(model,"LC2")) 2L else 1L
    X <- do.call(rbind,lapply(CAUSES,function(c)fits[[key(sex,c,model)]]$K[,factor]))
    rownames(X) <- SHORT
    Z <- z_rows(X)
    d <- distance_matrix(Z,p=2,check_impl=TRUE)
    eva <- evaluate_clusters(Z,d,seed=SEED)
    nm <- paste(sex,model)
    analyses[[nm]] <- list(Z=Z,raw_kappa=X,d=d,hc=eva$hc,labels=NULL,
       factor=factor,scores=eva$scores,proposal=eva$proposal,
       candidate_labels=eva$labels)
    sc <- eva$scores
    scores[[nm]] <- cbind(data.frame(sex=sex,model=model,factor=factor),sc)
    if(nrow(eva$starts))starts[[nm]] <-
      cbind(data.frame(sex=sex,model=model),eva$starts)
    for(k in 2:8) for(method in c("hierarchical_average","pam")) {
      lab <- eva$labels[[paste(method,k,sep="::")]]
      memberships[[length(memberships)+1L]] <-
        data.frame(sex=sex,model=model,method=method,K=k,
                   cause=CAUSES,code=SHORT,cluster=lab)
    }
    templates[[nm]] <- data.frame(sex=sex,model=model,
       factor=factor,suggested_K=eva$proposal$K[1L],
       suggested_method=eva$proposal$method[1L],
       K=NA_integer_,method="",reason="")
    for(i in seq_along(CAUSES)) {
      zrows[[length(zrows)+1L]] <-
        data.frame(sex=sex,model=model,cause=CAUSES[i],year=YEARS,z=Z[i,])
      if(i>1L) for(j in seq_len(i-1L))
        pairs[[length(pairs)+1L]] <- data.frame(
          sex=sex,model=model,cause1=CAUSES[i],cause2=CAUSES[j],dtw_p2=d[i,j])
    }
    message("§4.2候选完成 ",sex," / ",model,": ",nrow(sc),
            "个方案；建议仅供参考 K=",eva$proposal$K[1L])
  }
  write_csv(do.call(rbind,scores),file.path(out,"cluster_validity_scores.csv"))
  write_csv(do.call(rbind,memberships),file.path(out,"candidate_membership.csv"))
  write_csv(do.call(rbind,templates),file.path(out,"cluster_selection_template.csv"))
  write_csv(do.call(rbind,pairs),file.path(out,"pairwise_distances.csv"))
  write_csv(do.call(rbind,zrows),file.path(out,"standardized_kappa.csv"))
  write_csv(do.call(rbind,starts),file.path(out,"pam_three_start_check.csv"))
  saveRDS(analyses,file.path(out,"cluster_analyses.rds"))
  # Do not generate an unqualified cluster_membership.csv until a decision is made.
  if(is.null(selection_file) && selection_mode=="review_required") {
    message("已生成候选、七项CVI、三起点和待填写的决策表；正式分组尚未确定。")
    return(analyses)
  }
  finalize_cluster_selections(analyses,out,selection_file,selection_mode)
}
