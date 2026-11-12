# §4.2: DTW (dtw) + time-series clustering / official CVIs (dtwclust).
# Replaces earlier hand-coded approximations of seven CVIs and automatic Silhouette-only K.
# Reproduction is METHOD-BASED on GBD; the paper does not publish an automatic K rule.
check_cluster_packages <- function() {
  for (pkg in c("dtw", "dtwclust")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("缺少R包 ", pkg, "：请运行 install.packages('", pkg, "')")
  }
  invisible(TRUE)
}
z_rows <- function(X) {
  s <- apply(X, 1L, sd)
  if (any(!is.finite(s) | s < 1e-12)) stop("时期因子方差过小，不能标准化")
  sweep(sweep(X, 1L, rowMeans(X)), 1L, s, "/")
}
# The local scalar L2 cost is (x_i-y_j)^2. symmetric2 step weights:
# diagonal=2, vertical/horizontal=1; divide accumulated cost by (n+m), take sqrt.
# p=1 is retained ONLY for a sensitivity analysis (not the main paper setting).
dtw_distance <- function(x, y, p = 2, window = Inf, path = FALSE) {
  check_cluster_packages()
  stopifnot(p %in% c(1,2), length(x) > 1L, length(y) > 1L,
            all(is.finite(x)), all(is.finite(y)))
  lcm <- abs(outer(as.numeric(x), as.numeric(y), "-"))^p
  opts <- list(lcm, step.pattern = dtw::symmetric2,
               keep.internals = path, distance.only = !path)
  if (is.finite(window)) {
    opts$window.type <- "sakoechiba"
    opts$window.size <- as.integer(window)
  }
  fit <- do.call(dtw::dtw, opts)
  d <- (fit$distance/(length(x)+length(y)))^(1/p)
  if (!is.finite(d)) stop("DTW无可行路径")
  if (!path) return(as.numeric(d))
  list(distance = as.numeric(d), cost = lcm, acc = fit$costMatrix,
       path = cbind(fit$index1, fit$index2), raw = fit)
}
# Cross-check every pair against dtwclust's official dtw_basic implementation.
distance_matrix <- function(Z, p = 2, window = Inf, check_impl = FALSE) {
  n <- nrow(Z)
  if (n < 2L) stop("至少需要两条序列")
  d <- matrix(0, n, n, dimnames = list(rownames(Z), rownames(Z)))
  for (i in 2:n) for (j in seq_len(i-1L)) {
    x <- Z[i,]; y <- Z[j,]
    d[i,j] <- d[j,i] <- dtw_distance(x,y,p,window)
    if (check_impl && p == 2L) {
      b <- dtwclust::dtw_basic(x,y,norm="L2",
         step.pattern=dtw::symmetric2, normalize=TRUE, sqrt.dist=TRUE,
         window.size=if(is.finite(window)) as.integer(window) else NULL)
      if (!isTRUE(all.equal(d[i,j],as.numeric(b),tolerance=1e-7)))
        stop("dtw和dtwclust距离不一致: ",i," / ",j,": ",d[i,j]," vs ",b)
    }
  }
  d
}
ari <- function(a,b) {
  tab <- table(a,b); choose2 <- function(x) x*(x-1)/2
  total <- choose2(sum(tab)); if(total==0) return(1)
  x <- sum(choose2(tab)); aa <- sum(choose2(rowSums(tab)))
  bb <- sum(choose2(colSums(tab))); expected <- aa*bb/total
  den <- (aa+bb)/2-expected
  if(abs(den)<1e-14) return(as.numeric(identical(as.integer(a),as.integer(b))))
  (x-expected)/den
}
CVI_NAMES <- c("Sil","D","COP","DB","DBstar","CH","SF")
CVI_DIRECTION <- c(Sil=1,D=1,COP=-1,DB=-1,DBstar=-1,CH=1,SF=1)
# tsclust receives the exact same pre-computed distance matrix as hclust.
# PAM is run from three independently reproducible random seeds for each K.
cluster_fit <- function(Z, d, method, k, seed = 20261110L) {
  check_cluster_packages()
  params <- dtwclust::tsclust_args(
    dist=list(norm="L2",step.pattern=dtw::symmetric2,
              normalize=TRUE,sqrt.dist=TRUE))
  if (method=="hierarchical_average") {
    ans <- dtwclust::tsclust(
       series=Z,type="hierarchical",k=as.integer(k),
       distance="dtw_basic",args=params,
       control=dtwclust::hierarchical_control(
         method="average",symmetric=TRUE,distmat=as.dist(d)))
  } else if (method=="pam") {
    ans <- dtwclust::tsclust(
       series=Z,type="partitional",k=as.integer(k),centroid="pam",
       distance="dtw_basic",seed=as.integer(seed),args=params,
       control=dtwclust::partitional_control(
         nrep=1L,symmetric=TRUE,pam.precompute=TRUE,
         distmat=as.dist(d)))
  } else stop("不支持的聚类方法：",method)
  if (is.list(ans) && !methods::is(ans,"TSClusters")) ans <- ans[[1L]]
  ans
}
medoid_objective <- function(d, labs) {
  groups <- split(seq_along(labs),labs)
  sum(vapply(groups,function(g) min(colSums(d[g,g,drop=FALSE])),0.0))
}
cluster_validity <- function(fit) {
  results <- setNames(rep(NA_real_,length(CVI_NAMES)),CVI_NAMES)
  problems <- character()
  for (nm in CVI_NAMES) {
    v <- tryCatch(as.numeric(dtwclust::cvi(fit,type=nm))[1L],
                  error=function(e) {problems <<- c(problems,paste0(nm,": ",conditionMessage(e))); NA_real_})
    if(is.finite(v)) results[nm] <- v
    else if(!any(startsWith(problems,paste0(nm,":"))))
      problems <- c(problems,paste0(nm,": nonfinite"))
  }
  list(values=results,issues=paste(problems,collapse=" | "))
}
# Candidates are not decisions. Use all seven *official* CVIs, two methods,
# three PAM initializations. No hard-coded K=5/6 and no 'Sil only' decision.
evaluate_clusters <- function(Z, d, k_values = 2:(nrow(d)-1L), seed = 20261110L) {
  k_values <- sort(unique(as.integer(k_values)))
  stopifnot(all(k_values >= 2L & k_values <= nrow(d)-1L))
  rows <- list(); labels <- list(); starts <- list(); hc <- hclust(as.dist(d),method="average")
  for (k in k_values) {
    for (method in c("hierarchical_average","pam")) {
      if (method=="pam") {
        candidates <- lapply(1:3,function(s)
          cluster_fit(Z,d,method,k,seed=seed+100L*k+as.integer(s)))
        objectives <- vapply(candidates,function(f) medoid_objective(d,f@cluster),0.0)
        best <- which.min(objectives); fit <- candidates[[best]]
        for (s in 1:3) starts[[length(starts)+1L]] <- data.frame(
          method=method,K=k,start=s,seed=as.integer(seed+100L*k+s),
          objective=objectives[s],best_objective=min(objectives),
          ARI_vs_best=ari(candidates[[s]]@cluster,fit@cluster))
      } else {
        fit <- cluster_fit(Z,d,method,k,seed)
      }
      lab <- as.integer(fit@cluster)
      if(length(lab)!=nrow(Z) || length(unique(lab))!=k)
        stop("聚类标签数量异常：",method," K=",k)
      score <- cluster_validity(fit)
      nm <- paste(method,k,sep="::")
      labels[[nm]] <- lab
      rows[[nm]] <- data.frame(method=method,K=k,
        singleton_count=sum(tabulate(lab)==1L),
        objective=medoid_objective(d,lab),
        t(score$values),cvi_issues=score$issues,check.names=FALSE,
        stringsAsFactors=FALSE)
    }
  }
  tab <- do.call(rbind,rows); rownames(tab) <- NULL
  # Each metric is ranked within each method, never silently replaced by a proxy.
  rank_columns <- paste0("rank_",CVI_NAMES)
  for (nm in CVI_NAMES) tab[[paste0("rank_",nm)]] <- NA_real_
  for(method in unique(tab$method)) {
    sel <- which(tab$method==method)
    for(nm in CVI_NAMES) {
      v <- tab[[nm]][sel]; ok <- which(is.finite(v))
      if(length(ok)) tab[[paste0("rank_",nm)]][sel[ok]] <-
        rank(-CVI_DIRECTION[[nm]]*v[ok],ties.method="average")
    }
  }
  rr <- as.matrix(tab[,rank_columns,drop=FALSE])
  tab$mean_indicator_rank <- rowMeans(rr,na.rm=TRUE)
  tab$mean_indicator_rank[!is.finite(tab$mean_indicator_rank)] <- NA_real_
  # Exploratory suggestion only; NOT a rule published by Varga (2025).
  ks <- sort(unique(tab$K))
  suggest <- data.frame(K=ks,mean_rank=vapply(ks,function(k)
      mean(tab$mean_indicator_rank[tab$K==k],na.rm=TRUE),0.0),
    singleton_min=vapply(ks,function(k) min(tab$singleton_count[tab$K==k]),0L))
  suggest <- suggest[order(suggest$mean_rank,suggest$singleton_min,suggest$K),]
  proposed_k <- suggest$K[1L]
  candidates <- tab[tab$K==proposed_k,]
  candidates <- candidates[order(candidates$mean_indicator_rank,
                                candidates$method!="hierarchical_average"),]
  proposal <- candidates[1L,c("method","K","mean_indicator_rank","singleton_count")]
  list(scores=tab,labels=labels,starts=if(length(starts))do.call(rbind,starts) else data.frame(),
       hc=hc,proposal=proposal,k_summary=suggest)
}
# Run-only feature to ensure all 14 decisions are explicit and auditable.
read_selection_file <- function(path, models, sexes=c("Male","Female")) {
  if (is.null(path)) return(NULL)
  x <- read.csv(path,stringsAsFactors=FALSE)
  if(!all(c("sex","model","method","K") %in% names(x)))
    stop("选K文件须含 sex,model,method,K 四列")
  wanted <- expand.grid(sex=sexes,model=models,stringsAsFactors=FALSE)
  tag <- function(z) paste(z$sex,z$model,sep="::")
  if(nrow(x)!=nrow(wanted) || anyDuplicated(tag(x)) ||
     !setequal(tag(x),tag(wanted))) stop("选K文件须包含14个不重复的性别×模型组合")
  if(any(!x$method %in% c("hierarchical_average","pam")) ||
     any(is.na(x$K)) || any(x$K<2 | x$K>8) || any(x$K!=as.integer(x$K)))
    stop("手工选K须有合法的 method 和整数K(2:8)")
  x
}
