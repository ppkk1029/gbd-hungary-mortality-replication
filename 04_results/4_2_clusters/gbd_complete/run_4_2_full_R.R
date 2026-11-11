# Varga (2025) §4.2 | COMPLETE R-ONLY GBD method replication.
# All seven single-population fits, 14 DTW cluster analyses, Figures 5–8.
# No Python is invoked. Base R only (stats, utils, graphics, grDevices).
# Input: ORIGINAL fitted_rates.csv with deaths_gbd and exposure_proxy.
# Study this together with the Notion §4.2 R code lesson.
#
# RStudio: Source this file -> select original CSV -> choose output folder.
# Console: Rscript run_4_2_full_R.R input.csv result_directory
# Caution: GBD deaths are estimated/noninteger; binomial is pseudo-likelihood.
# HCSO vs GBD and K selection are different. This is a METHOD replication.
# This implementation has NOT been independently executed in an R runtime here.

SEED <- 20261110L
YEARS <- 1980:2015
AGE_LABELS <- c("60-64 years", "65-69 years", "70-74 years",
                "75-79 years", "80-84 years", "85-89 years")
AGE_MID <- seq(62, 87, 5)
SEXES <- c("Male", "Female")
CAUSES <- c("Cardiovascular diseases", "Neoplasms", "Digestive diseases",
            "Chronic respiratory diseases", "Neurological disorders",
            "Unintentional injuries", "Diabetes and kidney diseases",
            "Self-harm and interpersonal violence",
            "Respiratory infections and tuberculosis")
SHORT <- c("CVD", "NEO", "DIG", "CRD", "NEU", "INJ", "DKD", "SHV", "RIT")
MODELS <- c("Poisson LC1", "Binomial LC1", "Poisson LC2", "Binomial LC2",
            "CBD", "Plat2", "Plat3")
set.seed(SEED)

# ----- 1. Identify fitted factors (LC models must be anchored) -----
canonical_lc <- function(alpha, B, K) {
  center <- colMeans(K)
  alpha <- as.numeric(alpha + B %*% center)
  K <- sweep(K, 2L, center, "-")
  S <- svd(B %*% t(K))
  r <- ncol(B)
  B <- S$u[, seq_len(r), drop = FALSE]
  K <- sweep(S$v[, seq_len(r), drop = FALSE], 2L, S$d[seq_len(r)], "*")
  for (j in seq_len(r)) {
    ix <- which.max(abs(B[, j]))
    if (B[ix, j] < 0) {
      B[, j] <- -B[, j]
      K[, j] <- -K[, j]
    }
  }
  list(alpha = alpha, B = B, K = K, singular = S$d)
}

# ----- 2. Estimate one single-sex / single-cause model -----
fit_one <- function(D, E, model, starts = 3L) {
  stopifnot(is.matrix(D), is.matrix(E), all(dim(D) == c(6L, 36L)),
            all(dim(E) == dim(D)), all(is.finite(D)), all(D > 0),
            all(is.finite(E)), all(E > 0))
  A <- nrow(D); Tn <- ncol(D)
  binom <- grepl("^Binomial", model) || identical(model, "CBD")
  lc <- grepl("LC[12]$", model)
  N <- if (binom) E + D / 2 else E
  if (binom && any(N <= D)) stop("Binomial needs N>D; verify deaths/exposure.")
  z <- if (binom) log(D / (N - D)) else log(D / E)
  total <- sum(D)
  a0 <- rowMeans(z)

  if (lc) {
    r <- as.integer(substr(model, nchar(model), nchar(model)))
    sv <- svd(sweep(z, 1L, a0, "-"))
    U <- sv$u[, seq_len(r), drop = FALSE]
    V <- sv$v[, seq_len(r), drop = FALSE]
    if (r == 1L) {
      anchor <- which.max(abs(U[, 1L]))
    } else {
      candidates <- combn(seq_len(A), r)
      detabs <- apply(candidates, 2L, function(ix)
        abs(det(U[ix, , drop = FALSE])))
      anchor <- candidates[, which.max(detabs)]
    }
    free <- setdiff(seq_len(A), anchor)
    M <- U[anchor, , drop = FALSE]
    B0 <- U %*% solve(M)
    K0 <- sweep(V, 2L, sv$d[seq_len(r)], "*") %*% t(M)
    cent <- colMeans(K0)
    a0 <- as.numeric(a0 + B0 %*% cent)
    K0 <- sweep(K0, 2L, cent, "-")
    pack <- function(a, B, K)
      c(a, as.vector(B[free, , drop = FALSE]),
        as.vector(K[-Tn, , drop = FALSE]))
    unpack <- function(v) {
      a <- v[seq_len(A)]
      B <- matrix(0, nrow = A, ncol = r)
      B[anchor, ] <- diag(r)
      off <- A + length(free) * r
      B[free, ] <- matrix(v[(A+1L):off], nrow = length(free), ncol = r)
      K <- matrix(0, nrow = Tn, ncol = r)
      K[-Tn, ] <- matrix(v[(off+1L):length(v)], nrow = Tn-1L, ncol = r)
      K[Tn, ] <- -colSums(K[-Tn, , drop = FALSE])
      list(a = a, B = B, K = K)
    }
    init <- pack(a0, B0, K0)
  } else {
    r <- if (identical(model, "Plat3")) 3L else 2L
    B0 <- cbind(rep(1, A), (mean(AGE_MID)-AGE_MID)/10)
    if (identical(model, "CBD")) B0[, 2L] <- -B0[, 2L]
    if (r == 3L)
      B0 <- cbind(B0, pmax((mean(AGE_MID)-AGE_MID)/10, 0))
    if (identical(model, "CBD")) a0 <- rep(0, A)
    K0 <- t(qr.solve(B0, sweep(z, 1L, a0, "-")))
    if (!identical(model, "CBD")) {
      cent <- colMeans(K0)
      a0 <- as.numeric(a0 + B0 %*% cent)
      K0 <- sweep(K0, 2L, cent, "-")
    }
    pack <- function(a, B, K) {
      if (identical(model, "CBD")) return(as.vector(K))
      c(a, as.vector(K[-Tn, , drop = FALSE]))
    }
    unpack <- function(v) {
      if (identical(model, "CBD")) {
        return(list(a = rep(0, A), B = B0,
                    K = matrix(v, nrow = Tn, ncol = r)))
      }
      a <- v[seq_len(A)]
      K <- matrix(0, nrow = Tn, ncol = r)
      K[-Tn, ] <- matrix(v[-seq_len(A)], nrow = Tn-1L, ncol = r)
      K[Tn, ] <- -colSums(K[-Tn, , drop = FALSE])
      list(a = a, B = B0, K = K)
    }
    init <- pack(a0, B0, K0)
  }

  log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  calc <- function(v) {
    q <- unpack(v)
    eta <- sweep(q$B %*% t(q$K), 1L, q$a, "+")
    if (binom) {
      mu <- N * plogis(eta)
      obj <- sum(N * log1pexp(eta) - D * eta) / total
    } else {
      mu <- E * exp(pmin(pmax(eta, -700), 50))
      obj <- sum(mu - D * eta) / total
    }
    residual <- (mu - D) / total
    ga <- rowSums(residual)
    gk <- t(residual) %*% q$B
    if (lc) {
      gb <- residual %*% q$K
      g <- c(ga, as.vector(gb[free, , drop = FALSE]),
             as.vector(sweep(gk[-Tn, , drop = FALSE], 2L,
                             gk[Tn, ], "-")))
    } else if (identical(model, "CBD")) {
      g <- as.vector(gk)
    } else {
      g <- c(ga, as.vector(sweep(gk[-Tn, , drop = FALSE],
                                2L, gk[Tn, ], "-")))
    }
    list(value = obj, gradient = g, eta = eta, mu = mu)
  }

  initlist <- list(init)
  for (i in 2:starts) initlist[[i]] <- init + rnorm(length(init), sd = 0.015)
  fits <- lapply(initlist, function(v) {
    tryCatch(optim(v, fn = function(p) calc(p)$value,
                   gr = function(p) calc(p)$gradient, method = "BFGS",
                   control = list(maxit = 4000, reltol = 1e-12)),
             error = function(e) list(value=Inf, convergence=999,
                                      par=v, message=conditionMessage(e)))
  })
  ib <- which.min(vapply(fits, function(x) x$value, numeric(1)))
  if (!is.finite(fits[[ib]]$value)) stop("All starts failed: ", model)
  best <- unpack(fits[[ib]]$par)
  cb <- calc(fits[[ib]]$par)
  besteta <- cb$eta
  grad <- max(abs(cb$gradient))
  deviations <- vapply(fits, function(o) {
    if (!is.finite(o$value)) return(NA_real_)
    max(abs(calc(o$par)$eta - besteta))
  }, numeric(1))
  if (lc) {
    canon <- canonical_lc(best$a, best$B, best$K)
    best$a <- canon$alpha; best$B <- canon$B; best$K <- canon$K
  }
  reconstructed <- sweep(best$B %*% t(best$K), 1L, best$a, "+")
  err <- max(abs(reconstructed - besteta))
  accepted <- is.finite(grad) && grad < 1e-6 && err < 1e-6
  if (!accepted) warning(sprintf("Numerical check failed for %s: gradient %.4g, reconstruction %.4g",
                                 model, grad, err))
  list(a=best$a, B=best$B, K=best$K,
       diag=data.frame(model=model, accepted=accepted,
                       max_abs_scaled_gradient=grad,
                       convergence=fits[[ib]]$convergence,
                       objective=fits[[ib]]$value,
                       objective_range_3_starts=diff(range(vapply(fits, function(o) o$value, numeric(1)))),
                       max_eta_difference_3_starts=max(deviations,na.rm=TRUE),
                       reconstruction_error=err))
}

# ----- 3. Dynamic time warping; squared local cost + symmetric2 -----
dtw_r <- function(x, y, p=2, window=NULL, keep_path=FALSE) {
  n <- length(x); m <- length(y)
  stopifnot(all(is.finite(x)), all(is.finite(y)), n>0, m>0)
  cost <- abs(outer(x,y,"-"))^p
  F <- matrix(Inf,n,m)
  prev <- matrix(0L,n,m) # 1 diagonal, 2 up, 3 left
  F[1,1] <- cost[1,1]
  for (i in seq_len(n)) for (j in seq_len(m)) {
    if (i==1L && j==1L) next
    if (!is.null(window) && abs(i-j)>window) next
    opts <- c(
      if (i>1 && j>1) F[i-1,j-1]+2*cost[i,j] else Inf,
      if (i>1) F[i-1,j]+cost[i,j] else Inf,
      if (j>1) F[i,j-1]+cost[i,j] else Inf
    )
    z <- which.min(opts)
    F[i,j] <- opts[z];prev[i,j] <- z
  }
  if (!is.finite(F[n,m])) stop("DTW alignment path not found.")
  d <- (F[n,m]/(n+m))^(1/p)
  if (!keep_path) return(d)
  path <- matrix(integer(0),ncol=2)
  i <- n;j <- m
  repeat {
    path <- rbind(c(i,j),path)
    if (i==1L && j==1L) break
    direction <- prev[i,j]
    if (direction==1L) {i<-i-1L;j<-j-1L}
    else if (direction==2L) {i<-i-1L}
    else if (direction==3L) {j<-j-1L}
    else stop("Invalid DTW predecessor.")
  }
  colnames(path) <- c("i","j")
  list(distance=d,path=path,local_cost=cost,accumulated=F)
}
distance_matrix <- function(Z,p=2,window=NULL) {
  n <- nrow(Z);D <- matrix(0,n,n)
  for (i in seq_len(n)) if (i>1L)
    for (j in seq_len(i-1L))
      D[i,j] <- D[j,i] <- dtw_r(Z[i,],Z[j,],p,window)
  dimnames(D) <- list(rownames(Z),rownames(Z))
  D
}

# ----- 4. Cluster against a precomputed DTW distance matrix -----
medoids_exact <- function(d,k) {
  cand <- combn(seq_len(nrow(d)),k)
  objective <- apply(cand,2L,function(m) sum(apply(d[,m,drop=FALSE],1L,min)))
  med <- cand[,which.min(objective)]
  lab <- max.col(-d[,med,drop=FALSE],ties.method="first")
  lab[med] <- seq_len(k)
  list(lab=lab,med=med,objective=min(objective))
}
adjusted_rand <- function(a,b) {
  tab <- table(a,b)
  C2 <- function(z) z*(z-1)/2
  nij <- sum(C2(tab)); ai <- sum(C2(rowSums(tab))); bj <- sum(C2(colSums(tab)))
  total <- C2(length(a))
  expval <- if (total==0) 0 else ai*bj/total
  den <- (ai+bj)/2-expval
  if (den==0) return(if (identical(outer(a,a,"=="),outer(b,b,"=="))) 1 else 0)
  (nij-expval)/den
}
silhouette_mean <- function(d,lab) {
  n <- nrow(d);out <- numeric(n)
  for (i in seq_len(n)) {
    same <- which(lab==lab[i] & seq_len(n)!=i)
    if (length(same)==0L) {out[i]<-0;next}
    a <- mean(d[i,same])
    others <- setdiff(unique(lab),lab[i])
    b <- min(vapply(others,function(k) mean(d[i,lab==k]),numeric(1)))
    out[i] <- if (max(a,b)==0) 0 else (b-a)/max(a,b)
  }
  mean(out)
}
validity <- function(d,lab,Z) {
  groups <- split(seq_len(nrow(d)),lab)
  k <- length(groups);n <- nrow(d)
  meds <- vapply(groups,function(g) g[which.min(rowSums(d[g,g,drop=FALSE]))],integer(1))
  scat <- mapply(function(g,m) mean(d[g,m]),groups,as.list(meds))
  diam <- max(vapply(groups,function(g) max(d[g,g,drop=FALSE]),numeric(1)))
  sep <- min(vapply(combn(seq_len(k),2,simplify=FALSE),
                    function(ij) min(d[groups[[ij[1]]],groups[[ij[2]]],drop=FALSE]),numeric(1)))
  db <- mean(vapply(seq_len(k),function(i)
    max((scat[i]+scat[-i])/pmax(d[meds[i],meds[-i]],1e-15)),numeric(1)))
  global <- which.min(rowSums(d))
  within <- sum(mapply(function(g,m) sum(d[g,m]^2),groups,as.list(meds)))
  between <- sum(mapply(function(g,m) length(g)*d[m,global]^2,groups,as.list(meds)))
  centers <- lapply(groups,function(g) colMeans(Z[g,,drop=FALSE]))
  center <- colMeans(Z)
  we <- sum(mapply(function(g,c) sum(sweep(Z[g,,drop=FALSE],2,c,"-")^2),
                   groups,centers))
  be <- sum(mapply(function(g,c) length(g)*sum((c-center)^2),groups,centers))
  obj <- sum(mapply(function(g,m) sum(d[g,m]),groups,as.list(meds)))
  data.frame(silhouette=silhouette_mean(d,lab),
             dunn=sep/max(diam,1e-15),
             DB_medoid=db,
             CH_squared_medoid=(between/(k-1))/(within/(n-k)),
             CH_euclidean=(be/(k-1))/(we/(n-k)),
             objective=obj)
}
cluster_choose <- function(d,Z) {
  hc <- hclust(as.dist(d),method="average")
  rows <- list();labs <- list()
  for (k in 2:(nrow(d)-1L)) for (method in c("average","kmedoids_exact")) {
    lab <- if (method=="average") cutree(hc,k) else medoids_exact(d,k)$lab
    score <- validity(d,lab,Z)
    rows[[length(rows)+1L]] <- data.frame(method=method,K=k,score)
    labs[[paste(method,k)]] <- lab
  }
  tab <- do.call(rbind,rows)
  ix <- order(-tab$silhouette,tab$K,tab$method)[1L]
  list(tab=tab,best=tab[ix,,drop=FALSE],
       lab=labs[[paste(tab$method[ix],tab$K[ix])]],hc=hc)
}
zscore <- function(X) {
  s <- apply(X,1L,sd)
  if (any(!is.finite(s) | s < 1e-12)) stop("Invalid/constant period factor.")
  sweep(sweep(X,1L,rowMeans(X),"-"),1L,s,"/")
}

# ----- 5. Figures 5–8 in base R, from actual results -----
plot_alignment <- function(Z,sex,other,number,out) {
  a <- Z[1,];b <- Z[other,]
  path <- dtw_r(a,b,keep_path=TRUE)
  png(file.path(out,sprintf("Figure_%d_%s_DTW_R.png",number,sex)),
      width=1900,height=620,res=150)
  op<-par(mfrow=c(1,3),mar=c(4,4,3,1))
  plot(YEARS,a,type="l",lwd=2,col="#1b6f90",xlab="Year",ylab="z(kappa)",
       main=sprintf("Figure %d | %s | Poisson LC2",number,sex))
  lines(YEARS,b,lwd=2,col="#cb713c")
  legend("topright",legend=c("CVD",SHORT[other]),col=c("#1b6f90","#cb713c"),
         lty=1,bty="n")
  image(seq_len(36),seq_len(36),path$local_cost,xlab="CVD year index",
        ylab="Other cause year index",main=paste("DTW =",round(path$distance,4)),
        col=hcl.colors(55,"YlOrRd"))
  lines(path$path[,1],path$path[,2],col="white",lwd=3)
  plot(YEARS,a,type="l",lwd=2,col="#1b6f90",xlab="Year",ylab="z(kappa)",
       ylim=range(c(a,b+4)),main="Optimal correspondence")
  lines(YEARS,b+4,col="#cb713c",lwd=2)
  for (ii in seq(1,nrow(path$path),by=2)) {
    ij<-path$path[ii,]
    segments(YEARS[ij[1]],a[ij[1]],YEARS[ij[2]],b[ij[2]]+4,
             col=adjustcolor("gray50",alpha.f=.55))
  }
  par(op);dev.off()
  write.csv(data.frame(year_1=YEARS[path$path[,1]],year_2=YEARS[path$path[,2]]),
            file.path(out,paste0("alignment_path_",sex,"_R.csv")),row.names=FALSE)
}
plot_groups <- function(Z,lab,sex,model,fig_no,out) {
  ng <- length(unique(lab))
  png(file.path(out,sprintf("Figure_%d_%s_%s_R.png",fig_no,sex,
                           gsub(" ","_",model))),
      width=850*ng,height=670,res=150)
  op<-par(mfrow=c(1,ng),mar=c(4,4,3,1))
  yy <- range(Z)
  pal <- grDevices::hcl.colors(nrow(Z),"Dark 3")
  for (g in sort(unique(lab))) {
    ix <- which(lab==g)
    plot(YEARS,Z[ix[1],],type="l",col=pal[ix[1]],lwd=2,
         ylim=yy,xlab="Year",ylab="Standardized kappa",
         main=paste0("Cluster ",g))
    if (length(ix)>1L) for (j in ix[-1L])
      lines(YEARS,Z[j,],col=pal[j],lwd=2)
    legend("topleft",legend=SHORT[ix],col=pal[ix],lty=1,lwd=2,
           cex=.76,bty="n")
    abline(h=0,lty=3,col="gray75")
  }
  mtext(sprintf("GBD Hungary | %s | %s",sex,model),outer=TRUE,line=-1,cex=1)
  par(op);dev.off()
}

# ----- 6. Driver and fully archived CSVs -----
main <- function(input, out) {
  raw <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
  fields <- c("sex_name","cause_name","age_name","year",
              "deaths_gbd","exposure_proxy")
  if (!all(fields %in% names(raw)))
    stop("Missing columns: ",paste(setdiff(fields,names(raw)),collapse=", "))
  if (nrow(raw)!=3888L || anyDuplicated(raw[c("sex_name","cause_name","age_name","year")]))
    stop("Wrong row count (need 3888), or duplicate sex/cause/age/year.")
  if (!setequal(raw$year,YEARS) || !setequal(raw$age_name,AGE_LABELS) ||
      !setequal(raw$sex_name,SEXES) || !setequal(raw$cause_name,CAUSES))
    stop("Scope mismatch: require exactly Hungary GBD 9 causes, sexes, ages, 1980:2015.")
  if (!all(is.finite(raw$deaths_gbd)) || any(raw$deaths_gbd<=0) ||
      !all(is.finite(raw$exposure_proxy)) || any(raw$exposure_proxy<=0))
    stop("Nonpositive or missing deaths/exposure.")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  kappas<-list();params<-list();diagnostics<-list();fits<-list()
  for (sex in SEXES) for (cause in CAUSES) {
    f <- raw[raw$sex_name==sex & raw$cause_name==cause,]
    f <- f[order(match(f$age_name,AGE_LABELS),f$year),]
    if (nrow(f)!=216L) stop("Each sex/cause needs 6x36 cells.")
    D <- matrix(f$deaths_gbd,nrow=6L,ncol=36L,byrow=TRUE)
    E <- matrix(f$exposure_proxy,nrow=6L,ncol=36L,byrow=TRUE)
    for (model in MODELS) {
      fit <- fit_one(D,E,model)
      fits[[paste(sex,cause,model,sep="||")]] <- fit
      diagnostics[[length(diagnostics)+1L]] <- cbind(
        data.frame(sex=sex,cause=cause),fit$diag)
      for (j in seq_len(ncol(fit$K))) {
        kappas[[length(kappas)+1L]] <- data.frame(
          sex=sex,cause=cause,model=model,factor=j,year=YEARS,
          kappa=fit$K[,j])
        params[[length(params)+1L]] <- data.frame(
          sex=sex,cause=cause,model=model,factor=j,age=AGE_MID,
          alpha=fit$a,beta=fit$B[,j])
      }
    }
    message("Fitted ",sex," / ",cause)
  }
  write.csv(do.call(rbind,diagnostics),
            file.path(out,"fit_diagnostics_all_models.csv"),row.names=FALSE)
  write.csv(do.call(rbind,kappas),file.path(out,"kappa_all_models.csv"),row.names=FALSE)
  write.csv(do.call(rbind,params),file.path(out,"parameters_all_models.csv"),row.names=FALSE)

  summary<-list();members<-list();scores<-list();pairs<-list();zs<-list();sens<-list()
  for (sex in SEXES) for (model in MODELS) {
    factor <- if (grepl("LC2$",model)) 2L else 1L
    X <- t(vapply(CAUSES,function(ca)
      fits[[paste(sex,ca,model,sep="||")]]$K[,factor],numeric(36)))
    rownames(X) <- CAUSES
    Z <- zscore(X)
    d <- distance_matrix(Z)
    chosen <- cluster_choose(d,Z)
    best <- chosen$best
    lab <- chosen$lab
    summary[[length(summary)+1L]] <- data.frame(
      sex=sex,model=model,factor=factor,best)
    scores[[length(scores)+1L]] <- data.frame(
      sex=sex,model=model,chosen$tab)
    members[[length(members)+1L]] <- data.frame(
      sex=sex,model=model,cause=CAUSES,code=SHORT,
      cluster=as.integer(lab))
    for (i in seq_along(CAUSES)) {
      zs[[length(zs)+1L]] <- data.frame(
        sex=sex,model=model,cause=CAUSES[i],year=YEARS,z=Z[i,])
      if (i>1L) for (j in seq_len(i-1L)) {
        pairs[[length(pairs)+1L]] <- data.frame(
          sex=sex,model=model,cause1=CAUSES[i],
          cause2=CAUSES[j],dtw_p2=d[i,j])
      }
    }
    for (tag in c("window5","L1")) {
      ds <- if (tag=="window5") distance_matrix(Z,p=2,window=5) else distance_matrix(Z,p=1)
      alt <- cluster_choose(ds,Z)
      sens[[length(sens)+1L]] <- data.frame(
        sex=sex,model=model,variant=tag,K=alt$best$K,
        silhouette=alt$best$silhouette,ARI=adjusted_rand(lab,alt$lab))
    }
    if (sex=="Female") {
      keep <- setdiff(seq_len(9),7L)
      alt <- cluster_choose(d[keep,keep,drop=FALSE],Z[keep,,drop=FALSE])
      sens[[length(sens)+1L]] <- data.frame(
        sex=sex,model=model,variant="without_DKD",K=alt$best$K,
        silhouette=alt$best$silhouette,ARI=adjusted_rand(lab[keep],alt$lab))
    }
    if (grepl("LC2$",model)) for (angle in c(15,30,45)) {
      theta <- angle*pi/180
      R <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
      XR <- t(vapply(CAUSES,function(ca) {
        orig <- fits[[paste(sex,ca,model,sep="||")]]$K
        (orig %*% R)[,2L]
      },numeric(36)))
      ZR <- zscore(XR)
      alt <- cluster_choose(distance_matrix(ZR),ZR)
      sens[[length(sens)+1L]] <- data.frame(
        sex=sex,model=model,variant=paste0("rotation_",angle),
        K=alt$best$K,silhouette=alt$best$silhouette,
        ARI=adjusted_rand(lab,alt$lab))
    }
    png(file.path(out,paste0("dendrogram_",sex,"_",gsub(" ","_",model),".png")),
        width=1050,height=650,res=140)
    plot(chosen$hc,labels=SHORT,main=paste(sex,model,"K=",best$K),
         xlab="",sub="Average linkage / DTW distance")
    dev.off()
    if (model=="Poisson LC2")
      plot_alignment(Z,sex,if (sex=="Male") 3L else 7L,
                     if (sex=="Male") 5L else 6L,out)
    if (sex=="Male" && model=="CBD") plot_groups(Z,lab,sex,model,7L,out)
    if (sex=="Male" && model=="Binomial LC2")
      plot_groups(Z,lab,sex,model,8L,out)
    message("Clustered ",sex," / ",model," : K=",best$K)
  }
  tables <- list(cluster_summary=summary,cluster_membership=members,
                 cluster_validity_scores=scores,pairwise_distances=pairs,
                 standardized_kappa=zs,sensitivity=sens)
  for (name in names(tables))
    write.csv(do.call(rbind,tables[[name]]),
              file.path(out,paste0(name,".csv")),row.names=FALSE)
  meta <- c(sprintf("Input: %s",normalizePath(input,winslash="/")),
            paste("Run time:",Sys.time()),
            paste("R version:",R.version.string),
            "Data: GBD 2023 Hungary; 1980-2015; ages 60-89",
            "Models: 7 x sexes 2 x causes 9 = 126 fits",
            "DTW: zscore within cause, symmetric2 squared local p=2, sqrt(cost/72)",
            "Clustering: average linkage and exhaustive k-medoids, K=2..8",
            "Selection: maximum average silhouette; tie lower K",
            "Binomial E0: Ec + D/2 approximate pseudo-likelihood",
            "Not HCSO paper numeric replication. No multi-population forecasts here.")
  writeLines(meta,file.path(out,"run_manifest_R.txt"))
  tab <- do.call(rbind,summary)
  cat("\n=== §4.2 R FULL WORKFLOW SUMMARY ===\n")
  print(tab[,c("sex","model","K","silhouette")],row.names=FALSE)
  cat("\nSaved to: ",normalizePath(out,winslash="/"),"\n",sep="")
  invisible(tab)
}

# ---- Self-test: DTW dynamic programming against the hand-computed example ----
stopifnot(abs(dtw_r(c(0,1,3),c(0,2,3))-sqrt(2/6))<1e-12)
stopifnot(adjusted_rand(c(1,1,2,2),c(2,2,1,1))==1)

args <- commandArgs(trailingOnly=TRUE)
if (length(args)>=2L) {
  input <- normalizePath(args[1L],winslash="/",mustWork=TRUE)
  out <- args[2L]
} else if (interactive()) {
  cat("请选择含 ORIGINAL deaths_gbd 和 exposure_proxy 的 fitted_rates.csv\n")
  input <- normalizePath(file.choose(),winslash="/",mustWork=TRUE)
  out <- file.path(dirname(input),"Section4_2_full_R_results")
  cat("预计输出目录：",out,"\n")
  ans <- tolower(trimws(readline("输入 y 执行完整七模型R复现并允许覆盖同名输出：")))
  if (ans!="y") stop("已取消；没有运行模型。")
} else {
  stop("Usage: Rscript run_4_2_full_R.R input.csv output_directory")
}
main(input,out)
