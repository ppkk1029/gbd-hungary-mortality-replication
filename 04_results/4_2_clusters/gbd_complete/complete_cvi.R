# §4.2 checks; the seven CVIs are computed in clustering_functions.R with
# dtwclust::cvi(), NOT approximated by an unrelated home-made formula.
complete_cvi <- function(analyses,out) {
  f <- file.path(out,"cluster_validity_scores.csv")
  tab <- read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
  if(nrow(tab)!=196L) stop("应有14组×2种算法×7个候选K = 196条CVI")
  if(!all(CVI_NAMES %in% names(tab))) stop("CVI指标列缺失")
  if(!identical(sort(unique(tab$K)),2:8)) stop("候选K不是2:8")
  starts <- read.csv(file.path(out,"pam_three_start_check.csv"))
  if(nrow(starts)!=294L) stop("应有14×7×3=294条PAM起点记录")
  if(any(!is.finite(starts$objective)) || any(starts$ARI_vs_best < -1 | starts$ARI_vs_best > 1))
    stop("PAM三起点检查异常")
  preferences <- list()
  directions <- CVI_DIRECTION
  for(nm in unique(paste(tab$sex,tab$model,tab$method,sep="::"))) {
    x <- tab[paste(tab$sex,tab$model,tab$method,sep="::")==nm,,drop=FALSE]
    for(metric in CVI_NAMES) {
      values <- as.numeric(x[[metric]])
      ok <- which(is.finite(values))
      selected <- if(length(ok)) x$K[ok[which.max(directions[[metric]]*values[ok])]] else NA_integer_
      preferences[[length(preferences)+1L]] <- data.frame(
         sex=x$sex[1L],model=x$model[1L],method=x$method[1L],
         metric=metric,preferred_K=selected,
         direction=if(directions[[metric]]==1) "maximize" else "minimize",
         status=if(length(ok)) "computed" else "not_available")
    }
  }
  write_csv(do.call(rbind,preferences),file.path(out,"metric_preferences.csv"))
  problem <- tab[!is.na(tab$cvi_issues) & nzchar(tab$cvi_issues),
                 c("sex","model","method","K","cvi_issues"),drop=FALSE]
  if(nrow(problem)) write_csv(problem,file.path(out,"cvi_warnings.csv"))
  message("完成CVI审查：196候选，294次PAM起点。计算失败的指标保留NA，不伪造结果。")
  invisible(list(candidates=tab,starts=starts,warnings=problem))
}
