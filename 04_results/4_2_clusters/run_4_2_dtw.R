# §4.2 | The clusters of causes of death (Varga, 2025)
# Current data-supported model: one-factor Poisson LC from our corrected §4.1
# The author clustered 7 model variants; only 1 model is available at this stage.
#
# RStudio: open file and click Source; or use source(file.choose()).
# Dialog 1: select the corrected kappa.csv from §4.1.
# Dialog 2: select corresponding fit_diagnostics.csv.
# Nothing is saved until ALL RStudio figures/tables are previewed and user confirms.
# R packages: dtw, cluster. If missing: install.packages(c("dtw","cluster"))
#
# Note: We use normalized symmetric2 DTW (Euclidean local cost, p=2,
# normalization by N+M), hierarchical average linkage and PAM.
# Original paper did not publish all software options/step constraints;
# our exact settings are documented, not claimed as identical.
# Any series flagged REVIEW is INCLUDED in exploratory full clustering,
# then we ALSO recompute sensitivity clusters with flags excluded.

if (!requireNamespace("dtw", quietly = TRUE) ||
    !requireNamespace("cluster", quietly = TRUE)) {
  stop("请先运行 install.packages(c('dtw','cluster'))，安装后重试。")
}
cat("\n[1/2] 选择更正后 §4.1 的 kappa.csv ...\n")
kp <- read.csv(file.choose(), stringsAsFactors = FALSE)
cat("[2/2] 选择同次 §4.1 的 fit_diagnostics.csv ...\n")
dg <- read.csv(file.choose(), stringsAsFactors = FALSE)

kfields <- c("sex_name", "cause_name", "year", "kappa")
dfields <- c("sex_name", "cause_name", "convergence_code",
             "beta_flag", "beta_cancellation_ratio",
             "max_abs_gradient")
if (!all(kfields %in% names(kp))) stop("kappa.csv 缺少列：", paste(setdiff(kfields,names(kp)),collapse=", "))
if (!all(dfields %in% names(dg))) stop("diagnostics.csv 缺少列；请先运行新版 §4.1 程序。")
stopifnot(!anyNA(kp[kfields]), !anyNA(dg[dfields]),
          !anyDuplicated(kp[c("sex_name","cause_name","year")]),
          !anyDuplicated(dg[c("sex_name","cause_name")]),
          all(is.finite(kp$kappa)), all(is.finite(dg$max_abs_gradient)))
if (!identical(sort(unique(kp$year)), 1980:2015)) {
  stop("时期因子年份必须与 §4.1 训练期 1980-2015 完全一致，不能混入回测年份。")
}
sexes <- c("Male", "Female")
if (!setequal(unique(kp$sex_name), sexes)) stop("男女数据不完整。")
causes <- unique(kp$cause_name[kp$sex_name == "Male"])
if (length(causes) != 9L ||
    !setequal(causes,unique(kp$cause_name[kp$sex_name == "Female"]))) {
  stop("两种性别必须具有相同的九个 GBD Level-2 死因。")
}
causes <- sort(causes)
if (!setequal(paste(dg$sex_name,dg$cause_name,sep="__"),
              paste(rep(sexes,each=9L),rep(causes,times=2L),sep="__"))) {
  stop("kappa 与 diagnostics 的死因/性别不是同一批模型。")
}
if (any(dg$convergence_code != 0L) ||
    any(dg$max_abs_gradient >= 0.02)) {
  stop("存在未收敛 §4.1 模型；暂停聚类，先修复拟合。")
}

# ---- 1. z-standardize EACH cause-specific factor over TRAINING years ----
Z <- list()
z_table <- data.frame()
for (sx in sexes) {
  m <- matrix(NA_real_, nrow=36L, ncol=9L,
              dimnames=list(as.character(1980:2015),causes))
  for (cs in causes) {
    z <- kp[kp$sex_name == sx & kp$cause_name == cs, ]
    z <- z[order(z$year), ]
    if (nrow(z) != 36L || !identical(z$year,1980:2015)) {
      stop("该时期因子年份不完整：",sx," ",cs)
    }
    sigma <- sd(z$kappa)
    if (!is.finite(sigma) || sigma < 1e-8) {
      stop("无法标准化几乎不变的时期因子：",sx," ",cs)
    }
    m[,cs] <- (z$kappa - mean(z$kappa)) / sigma
    z_table <- rbind(z_table,data.frame(
      sex_name=sx,cause_name=cs,year=z$year,
      kappa_raw=z$kappa, kappa_z=m[,cs]))
  }
  Z[[sx]] <- m
}
stopifnot(nrow(z_table)==2L*9L*36L)
cat("\nκ 标准化完成：18 条时间序列，各36年。所有标准化序列均均值0、标准差1。\n")

# ---- 2. Pairwise normalized symmetric2 DTW; also window-size-5 sensitivity ----
dtw_distance <- function(x,y,window=FALSE,return_alignment=FALSE) {
  args <- list(x=as.numeric(x),y=as.numeric(y),step.pattern=dtw::symmetric2,
               dist.method="Euclidean",
               distance.only=!return_alignment)
  if (window) {
    args$window.type <- "sakoechiba"
    args$window.size <- 5L
  }
  if (return_alignment) args$keep.internals <- TRUE
  ans <- do.call(dtw::dtw,args)
  if (return_alignment) return(ans)
  if (!is.finite(ans$normalizedDistance)) stop("DTW 归一化距离未定义。")
  as.numeric(ans$normalizedDistance)
}
build_distances <- function(M,window=FALSE) {
  n <- ncol(M)
  D <- matrix(0,nrow=n,ncol=n,dimnames=list(colnames(M),colnames(M)))
  for (i in seq_len(n-1L)) for (j in (i+1L):n) {
    D[i,j] <- D[j,i] <- dtw_distance(M[,i],M[,j],window=window)
  }
  stopifnot(max(abs(D-t(D)))<1e-10,all(diag(D)==0),all(D>=0))
  D
}
get_pairs <- function(M, D, sx, variant) {
  xy <- which(upper.tri(D),arr.ind=TRUE)
  data.frame(sex_name=sx,dtw_variant=variant,
             cause_1=colnames(M)[xy[,1]],
             cause_2=colnames(M)[xy[,2]],
             distance=D[xy],stringsAsFactors=FALSE)
}
D_full <- lapply(Z,build_distances,window=FALSE)
D_window <- lapply(Z,build_distances,window=TRUE)
pairwise <- do.call(rbind,lapply(sexes,function(sx)
  rbind(get_pairs(Z[[sx]],D_full[[sx]],sx,"full"),
        get_pairs(Z[[sx]],D_window[[sx]],sx,"SakoeChiba_5"))))
pairwise <- pairwise[order(pairwise$sex_name,pairwise$dtw_variant,pairwise$distance),]

# ---- 3. Hierarchical + PAM for K=2..6. Select highest average silhouette. ----
dunn_index <- function(cl,D) {
  groups <- split(seq_along(cl),cl)
  within <- vapply(groups,function(idx) {
    if(length(idx)<2L) return(0)
    max(D[idx,idx,drop=FALSE])
  },numeric(1))
  diam <- max(within)
  cross <- c()
  for(i in 1:(length(groups)-1L)) for(j in (i+1L):length(groups)) {
    cross <- c(cross,min(D[groups[[i]],groups[[j]],drop=FALSE]))
  }
  if (diam==0) return(NA_real_)
  min(cross)/diam
}
cluster_one <- function(D, Kmax=6L) {
  n <- nrow(D)
  hc <- hclust(as.dist(D),method="average")
  if (n < 4L) stop("有效死因少于4类，无法可靠筛选簇数。")
  scores <- data.frame()
  solutions <- list()
  for (kk in 2:min(Kmax,n-1L)) {
    for (meth in c("hierarchical_average","PAM")) {
      if (meth == "hierarchical_average") {
        label <- cutree(hc,k=kk)
      } else {
        set.seed(20261109)
        label <- cluster::pam(as.dist(D),k=kk,diss=TRUE)$clustering
      }
      sil <- cluster::silhouette(label,as.dist(D))
      ss <- mean(sil[,"sil_width"])
      di <- dunn_index(label,D)
      scores <- rbind(scores,data.frame(method=meth,k=kk,
                     mean_silhouette=ss,dunn_index=di))
      solutions[[paste0(meth,"__",kk)]] <- label
    }
  }
  # Only one explicit selection criterion: maximum mean silhouette.
  # Prefer fewer clusters if values tie.
  scores <- scores[order(-scores$mean_silhouette,scores$k,scores$method), ]
  best <- scores[1L,]
  label <- solutions[[paste0(best$method,"__",best$k)]]
  names(label) <- rownames(D)
  list(hc=hc,scores=scores,best=best,label=label)
}
fits <- lapply(D_full, cluster_one)

# ---- 4. A full-data clustering AND an exclusion sensitivity for flagged beta ----
flagged <- dg[dg$beta_flag != "OK",c("sex_name","cause_name",
                                    "beta_flag","beta_cancellation_ratio")]
labels <- data.frame()
for (sx in sexes) {
  sol <- fits[[sx]]
  labels <- rbind(labels,data.frame(
    sex_name=sx,cause_name=names(sol$label),
    cluster=as.integer(sol$label),
    method=sol$best$method,k=sol$best$k,
    flagged_beta=names(sol$label) %in% flagged$cause_name[flagged$sex_name==sx]))
}
sensitivity <- data.frame()
for (sx in sexes) {
  excluded <- flagged$cause_name[flagged$sex_name==sx]
  keep <- setdiff(causes, excluded)
  if (length(keep) >= 4L && length(excluded) > 0L) {
    reduced <- D_full[[sx]][keep,keep,drop=FALSE]
    alt <- cluster_one(reduced,Kmax=min(6L,length(keep)-1L))
    sensitivity <- rbind(sensitivity,data.frame(
      sex_name=sx,cause_name=names(alt$label),
      cluster_after_excluding_flagged=as.integer(alt$label),
      method=alt$best$method,k=alt$best$k,
      excluded_causes=paste(excluded,collapse="; ")))
  }
}

# ---- 5. Select closest unflagged pair for an example alignment ----
example_alignments <- list()
example_pairs <- data.frame()
for (sx in sexes) {
  eligible <- setdiff(causes,
                      flagged$cause_name[flagged$sex_name==sx])
  if(length(eligible)<2L) eligible <- causes
  D <- D_full[[sx]][eligible,eligible,drop=FALSE]
  idx <- which(upper.tri(D),arr.ind=TRUE)
  ix <- idx[which.min(D[idx]),,drop=FALSE]
  cs1 <- rownames(D)[ix[1,1]]
  cs2 <- colnames(D)[ix[1,2]]
  al <- dtw_distance(Z[[sx]][,cs1],Z[[sx]][,cs2],
                     return_alignment=TRUE)
  example_alignments[[sx]] <- list(al=al,cs1=cs1,cs2=cs2)
  example_pairs <- rbind(example_pairs,data.frame(
    sex_name=sx,cause_1=cs1,cause_2=cs2,
    normalized_dtw=al$normalizedDistance))
}

# ---- 6. Preview all results and MULTIPLE plots, not isolated figure chapters ----
score_table <- do.call(rbind,lapply(sexes,function(sx)
  cbind(sex_name=sx,fits[[sx]]$scores)))
cat("\n============= §4.2 PREVIEW — NO FILES SAVED =============\n")
cat("\n各性别簇数与方法选择：\n")
print(do.call(rbind,lapply(sexes,function(sx)
  cbind(sex_name=sx,fits[[sx]]$best))),row.names=FALSE)
cat("\n逐死因分组（有异常 beta 的死因已标识）：\n")
print(labels,row.names=FALSE)
cat("\n最近的可用时期因子对与 DTW 距离：\n")
print(example_pairs,row.names=FALSE)
if (nrow(flagged)>0) {
  cat("\n警告：以下时期因子有严重的 beta 正负抵消，不能直接给聚类赋予强解释：\n")
  print(flagged,row.names=FALSE)
  cat("\n已经另外计算排除异常死因的敏感性聚类（若剩余至少4类）。\n")
}
if (interactive()) {
  View(labels,"§4.2 Group membership")
  View(score_table,"§4.2 Silhouette and Dunn scores")
  View(pairwise,"§4.2 Pairwise DTW distances")
}

plot_overview <- function() {
  old <- par(no.readonly=TRUE)
  on.exit(par(old))
  par(mfrow=c(2,2),mar=c(4.4,4.3,3.5,1.2),oma=c(0,0,2.0,0))
  for(sx in sexes) {
    sol <- fits[[sx]]
    plot(sol$hc,main=paste(sx,"hierarchical average / normalized DTW"),
         xlab="",sub="",cex=0.65)
    rect.hclust(sol$hc,k=sol$best$k,border="gray55")
    s <- sol$scores
    yy <- range(c(-0.25,1,s$mean_silhouette))
    plot(NULL,xlim=c(1.75,6.25),ylim=yy,
         xlab="Number of clusters (K)",
         ylab="Average silhouette", main=paste(sx,"cluster validity"))
    for(meth in unique(s$method)) {
      xx <- s[s$method==meth,]
      points(xx$k,xx$mean_silhouette,pch=if(meth=="PAM")16 else 1,
             type="b",lty=if(meth=="PAM")1 else 2)
    }
    abline(h=0,col="gray75")
    legend("bottomright",legend=c("PAM","Hierarchical"),
           pch=c(16,1),lty=c(1,2),bty="n",cex=0.75)
  }
  mtext("§4.2 | DTW-based cause clustering (1980–2015)",
        outer=TRUE,side=3,line=0.5,font=2,cex=1.15)
}
plot_alignment <- function(sx) {
  old <- par(no.readonly=TRUE); on.exit(par(old))
  e <- example_alignments[[sx]]
  m <- Z[[sx]]
  al <- e$al
  par(mfrow=c(1,2),mar=c(4,4,3,1),oma=c(0,0,2,0))
  yrs <- 1980:2015
  plot(yrs,m[,e$cs1],type="l",lwd=2,col="#C05A8E",ylim=
       range(m[,c(e$cs1,e$cs2)]),xlab="Year",
       ylab="Standardized kappa",main=paste(sx,"closest pair"))
  lines(yrs,m[,e$cs2],lwd=2,col="#008F94",lty=2)
  legend("topleft",legend=c(e$cs1,e$cs2),lty=1:2,
         col=c("#C05A8E","#008F94"),bty="n",cex=0.65)
  image(seq_along(yrs),seq_along(yrs),al$localCostMatrix,
        xlab=paste("Year index:",e$cs1),
        ylab=paste("Year index:",e$cs2),
        main="Absolute local cost & optimal path",
        col=grDevices::hcl.colors(50,"YlOrRd",rev=TRUE))
  lines(al$index1,al$index2,lwd=2)
  mtext(paste("§4.2 DTW alignment:",sx),
        outer=TRUE,side=3,line=0.5,font=2)
}
plot_clusters <- function(sx) {
  old <- par(no.readonly=TRUE); on.exit(par(old))
  M <- Z[[sx]]
  cl <- fits[[sx]]$label
  groups <- split(names(cl),cl)
  nr <- ceiling(length(groups)/2)
  par(mfrow=c(nr,2),mar=c(3.8,3.9,3.0,1.0),oma=c(0,0,2,0))
  for(i in seq_along(groups)) {
    cs <- groups[[i]]
    plot(1980:2015,M[,cs[1]],type="n",ylim=range(M[,cs,drop=FALSE]),
         xlab="Year",ylab="z(kappa)",
         main=paste(sx,"Cluster",names(groups)),cex.main=0.9)
    for(j in seq_along(cs)) lines(1980:2015,M[,cs[j]],
                                  lwd=1.7,col=j)
    legend("topright",legend=cs,col=seq_along(cs),
           lty=1,lwd=1.5,cex=0.53,bty="n")
  }
  mtext(paste("§4.2 Standardized factors grouped by cause —",sx),
        outer=TRUE,side=3,line=0.5,font=2)
}
cat("\nRStudio Plots 历史中将依次展示：簇数选择、男女最佳对齐、男女聚类后的时期因子。\n")
plot_overview()
for(sx in sexes) plot_alignment(sx)
for(sx in sexes) plot_clusters(sx)

answer <- tolower(trimws(readline(
  "\n检查 §4.2 表格与 Plots 后输入 y 保存；其他输入取消："
)))
if (identical(answer,"y")) {
  if (.Platform$OS.type=="windows") {
    parent <- utils::choose.dir(default=getwd(),
              caption="请选择 §4.2 结果保存的父文件夹")
  } else if (requireNamespace("tcltk",quietly=TRUE)) {
    parent <- tcltk::tk_choose.dir(default=getwd(),
              caption="请选择 §4.2 结果保存的父文件夹")
  } else stop("无法打开文件夹选择器；请使用桌面版 RStudio。")
  if (is.na(parent) || !nzchar(parent)) {
    cat("已取消保存。结果仍在当前 R 会话内。\n")
  } else {
    out <- file.path(parent,"Section_4_2_DTW_clusters")
    proceed <- TRUE
    if (dir.exists(out)) {
      proceed <- identical(tolower(trimws(readline(
       "保存文件夹已存在。输入 y 允许覆盖同名文件，其他输入取消："))),"y")
    }
    if (proceed) {
      dir.create(out,recursive=TRUE,showWarnings=FALSE)
      write.csv(z_table,file.path(out,"standardized_kappa.csv"),row.names=FALSE)
      write.csv(pairwise,file.path(out,"dtw_pairwise_distances.csv"),row.names=FALSE)
      write.csv(score_table,file.path(out,"cluster_validity_scores.csv"),row.names=FALSE)
      write.csv(labels,file.path(out,"cluster_membership.csv"),row.names=FALSE)
      write.csv(flagged,file.path(out,"flagged_beta_factors.csv"),row.names=FALSE)
      write.csv(sensitivity,file.path(out,"sensitivity_exclude_flagged.csv"),row.names=FALSE)
      write.csv(example_pairs,file.path(out,"closest_factor_pairs.csv"),row.names=FALSE)
      save_png <- function(filename,fn) {
        grDevices::png(file.path(out,filename),width=2200,height=1650,res=180)
        tryCatch(fn(),finally=grDevices::dev.off())
      }
      save_png("clustering_overview.png",plot_overview)
      for(sx in sexes) {
        save_png(paste0("alignment_",tolower(sx),".png"),
                 function()plot_alignment(sx))
        save_png(paste0("standardized_clusters_",tolower(sx),".png"),
                 function()plot_clusters(sx))
      }
      writeLines(c(
        "Varga (2025) §4.2 GBD Level-2 METHOD replication",
        "CURRENT SCOPE: one-factor Poisson Lee-Carter ONLY; NOT all seven original model variants.",
        "Input: corrected §4.1 kappa.csv + fit_diagnostics.csv.",
        "Time: 1980-2015; GBD Hungary, ages 60-89, female/male separate.",
        "All factors z-standardized within each cause and sex.",
        "DTW Euclidean local cost (p=2), symmetric2, normalized by N+M.",
        "Sensitivity: Sakoe-Chiba window size=5.",
        "Clustering: hclust average linkage and deterministic cluster::pam.",
        "k=2..6 chosen by maximum average silhouette (also report Dunn).",
        "All flagged beta series included in FULL exploratory clustering.",
        "Flag-exclusion sensitivity computed where feasible.",
        "DTW similarity does not imply shared medical causal mechanism.",
        paste("Created:",Sys.time()),
        paste("Flagged series:",nrow(flagged))),
        file.path(out,"method_and_scope.txt"))
      cat("\n§4.2 结果已保存：",normalizePath(out,winslash="/",mustWork=TRUE),
          "\n",sep="")
      print(list.files(out))
    } else cat("取消覆盖，未保存文件。\n")
  }
} else {
  cat("未保存文件；聚类结果保留在 RStudio 对象 labels、fits 等中。\n")
}
