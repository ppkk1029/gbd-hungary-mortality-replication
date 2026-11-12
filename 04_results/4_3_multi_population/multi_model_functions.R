# §4.3 - Joint maximum-(pseudo)likelihood for all seven multi-population models.
# Pure base R. Requires the §4.2 model_functions.R constants and read_input().
# Data: GBD estimated deaths / population proxy, NOT integer registry counts.
# Table 2: LC1 shares k1; LC2 shares k2; CBD / Plat2 / Plat3 share k1.
# No averaging of separate fitted time series is reported as a joint fit.

get_single_fit <- function(pp, kk, sex, cause, model) {
  p <- pp[pp$sex == sex & pp$cause == cause & pp$model == model, ]
  k <- kk[kk$sex == sex & kk$cause == cause & kk$model == model, ]
  if (!nrow(p) || !nrow(k)) stop("§4.2单人口参数缺失：",sex,"/",cause,"/",model)
  nf <- if (model == "Plat3") 3L else if (model %in% c("Plat2","CBD") ||
                                            endsWith(model,"LC2")) 2L else 1L
  a <- numeric(length(AGES)); B <- matrix(0,length(AGES),nf)
  K <- matrix(0,length(YEARS),nf)
  for (h in seq_len(nf)) {
    ph <- p[p$factor == h, ]; kh <- k[k$factor == h, ]
    ph <- ph[match(AGES,ph$age), ]; kh <- kh[match(YEARS,kh$year), ]
    if (nrow(ph)!=length(AGES) || nrow(kh)!=length(YEARS) ||
        anyNA(ph$beta) || anyNA(kh$kappa)) stop("单人口参数年龄/年份不齐：",model)
    B[,h] <- ph$beta; K[,h] <- kh$kappa
    if (h==1L) a <- ph$alpha
  }
  list(a=a,B=B,K=K)
}

# All cause-specific and common temporal vectors have their explicit constraints.
# LC2 reference age/loadings additionally eliminate within-cause rotation:
# B1[anchor1_j,j]=1, B2[anchor1_j,j]=0, B2[anchor2_ref,ref]=1.
make_multi_problem <- function(D,E,model,singles) {
  dd <- dim(D); A<-dd[1]; TT<-dd[2]; J<-dd[3]
  stopifnot(identical(dim(E),dd),length(singles)==J,J>=2L,
            A==length(AGES),TT==length(YEARS),all(D>0),all(E>0))
  lc1 <- endsWith(model,"LC1")
  lc2 <- endsWith(model,"LC2")
  cbd <- model=="CBD"
  plat <- model %in% c("Plat2","Plat3")
  binom <- startsWith(model,"Binomial") || cbd
  N <- if(binom) E+D/2 else E
  if(binom && any(N<=D)) stop("Binomial approximated exposure <= deaths")
  n_total <- sum(D)
  w <- vapply(seq_len(J),function(j)sum(D[,,j]),0.0)
  ref <- which.max(w)
  w <- w/sum(w)

  alpha <- sapply(singles,function(x)x$a)
  if (is.null(dim(alpha))) alpha <- matrix(alpha,A,J)
  kshare_index <- if(lc2) 2L else 1L
  ks0 <- vapply(singles,function(x)x$K[,kshare_index],numeric(TT))
  # Align LC signs: single-population SVD signs are not globally comparable.
  if (lc1 || lc2) for(j in seq_len(J)) {
    if(j!=ref && cor(ks0[,j],ks0[,ref])<0) {
      ks0[,j] <- -ks0[,j]
      singles[[j]]$B[,kshare_index] <- -singles[[j]]$B[,kshare_index]
    }
  }
  ks <- drop(ks0 %*% w)
  if (sd(ks)<1e-10) ks <- ks0[,ref]

  anchor1 <- integer(J); anchor2 <- NA_integer_
  mask1 <- mask2 <- NULL
  B1 <- B2 <- NULL
  Kc <- NULL
  if (lc1) {
    B1 <- sapply(singles,function(x)x$B[,1])
    if(is.null(dim(B1))) B1 <- matrix(B1,A,J)
    anchor1[ref] <- which.max(abs(B1[,ref]))
    scale <- B1[anchor1[ref],ref]
    if(abs(scale)<1e-9) stop("LC1 reference loading close to zero")
    B1 <- B1/scale; ks <- ks*scale
    mask1 <- matrix(TRUE,A,J); mask1[anchor1[ref],ref]<-FALSE
  } else if (lc2) {
    B1 <- sapply(singles,function(x)x$B[,1])
    B2 <- sapply(singles,function(x)x$B[,2])
    if(is.null(dim(B1))) B1 <- matrix(B1,A,J)
    if(is.null(dim(B2))) B2 <- matrix(B2,A,J)
    Kc <- array(0,c(TT,J,1L))
    for(j in seq_len(J)) {
      anchor1[j] <- which.max(abs(B1[,j]))
      sc <- B1[anchor1[j],j]
      if(abs(sc)<1e-10) stop("LC2 primary reference loading near zero")
      B1[,j] <- B1[,j]/sc
      Kc[,j,1] <- singles[[j]]$K[,1]*sc
      # Add shared k2 into cause-specific k1, then zero one loading in B2.
      h <- B2[anchor1[j],j]
      Kc[,j,1] <- Kc[,j,1]+h*ks
      B2[,j] <- B2[,j]-h*B1[,j]
    }
    opts <- setdiff(seq_len(A),anchor1[ref])
    anchor2 <- opts[which.max(abs(B2[opts,ref]))]
    sc <- B2[anchor2,ref]
    if(abs(sc)<1e-10) stop("LC2 second loading is numerically unidentified")
    B2 <- B2/sc; ks <- ks*sc
    mask1 <- mask2 <- matrix(TRUE,A,J)
    for(j in seq_len(J)) {
      mask1[anchor1[j],j] <- FALSE
      mask2[anchor1[j],j] <- FALSE
    }
    mask2[anchor2,ref] <- FALSE
  } else {
    S <- if(model=="Plat3") 2L else 1L
    Kc <- array(0,c(TT,J,S))
    for(j in seq_len(J)) for(h in seq_len(S))
      Kc[,j,h] <- singles[[j]]$K[,h+1L]
  }

  if (!cbd) {
    # Center each time effect and absorb its level into cause-specific alpha.
    level <- mean(ks); ks <- ks-level
    for(j in seq_len(J)) {
      b <- if(lc1) B1[,j] else if(lc2) B2[,j] else rep(1,A)
      alpha[,j] <- alpha[,j]+b*level
    }
    if(lc2) for(j in seq_len(J)) {
      level <- mean(Kc[,j,1]);Kc[,j,1] <- Kc[,j,1]-level
      alpha[,j] <- alpha[,j]+B1[,j]*level
    }
    if(plat) {
      V <- cbind((mean(AGES)-AGES)/10)
      if(model=="Plat3") V <- cbind(V,pmax(V[,1],0))
      for(j in seq_len(J)) for(h in seq_len(ncol(V))) {
        level <- mean(Kc[,j,h]);Kc[,j,h]<-Kc[,j,h]-level
        alpha[,j] <- alpha[,j]+V[,h]*level
      }
    }
  }
  V <- if(cbd) cbind((AGES-mean(AGES))/10) else if(plat) {
    v <- (mean(AGES)-AGES)/10
    if(model=="Plat3") cbind(v,pmax(v,0)) else cbind(v)
  } else NULL
  S <- if(lc2) 1L else if(cbd||plat) ncol(V) else 0L
  cent <- !cbd
  # Fixed coordinates are omitted from the optimization vector.
  sizes <- integer(0)
  pack <- function(u) {
    if(lc1) return(c(as.vector(u$a),u$B1[mask1],u$ks[-TT]))
    if(lc2) return(c(as.vector(u$a),u$B1[mask1],u$B2[mask2],
                     as.vector(u$Kc[-TT,,1,drop=FALSE]),u$ks[-TT]))
    out <- if(cbd) u$ks else c(as.vector(u$a),u$ks[-TT])
    for(h in seq_len(S)) {
      tmp <- u$Kc[,,h]
      out <- c(out,if(cbd)as.vector(tmp) else as.vector(tmp[-TT,,drop=FALSE]))
    }
    out
  }
  unpack <- function(v) {
    pos <- 0L
    take <- function(n) {i<-pos+seq_len(n);pos<<-pos+n;v[i]}
    if(lc1) {
      a<-matrix(take(A*J),A,J)
      b<-matrix(0,A,J);b[mask1]<-take(sum(mask1));b[anchor1[ref],ref]<-1
      h<-take(TT-1L);k<-c(h,-sum(h))
      u<-list(a=a,B1=b,ks=k)
    } else if(lc2) {
      a<-matrix(take(A*J),A,J)
      b1<-matrix(0,A,J);b1[mask1]<-take(sum(mask1))
      b2<-matrix(0,A,J);b2[mask2]<-take(sum(mask2))
      for(j in seq_len(J))b1[anchor1[j],j]<-1
      b2[anchor2,ref]<-1
      q<-matrix(take((TT-1L)*J),TT-1L,J)
      kc<-array(0,c(TT,J,1L));kc[-TT,,1]<-q;kc[TT,,1]<- -colSums(q)
      h<-take(TT-1L);k<-c(h,-sum(h))
      u<-list(a=a,B1=b1,B2=b2,Kc=kc,ks=k)
    } else {
      a<-if(cbd)matrix(0,A,J) else matrix(take(A*J),A,J)
      if(cbd) {k<-take(TT)} else {h<-take(TT-1L);k<-c(h,-sum(h))}
      kc<-array(0,c(TT,J,S))
      for(h in seq_len(S)) {
        if(cbd) kc[,,h]<-matrix(take(TT*J),TT,J) else {
          q<-matrix(take((TT-1L)*J),TT-1L,J)
          kc[-TT,,h]<-q;kc[TT,,h]<- -colSums(q)
        }
      }
      u<-list(a=a,Kc=kc,ks=k)
    }
    if(pos!=length(v))stop("Incorrect joint parameter vector length")
    u
  }
  u0 <- list(a=alpha,B1=B1,B2=B2,Kc=Kc,ks=ks)
  # Guard the pointwise beta anchors in initial values.
  if(lc1) u0$B1[anchor1[ref],ref]<-1
  if(lc2) {
    for(j in seq_len(J)) {u0$B1[anchor1[j],j]<-1;u0$B2[anchor1[j],j]<-0}
    u0$B2[anchor2,ref]<-1
  }
  initial <- pack(u0)

  calc <- function(v,return_details=FALSE) {
    u<-unpack(v)
    ga<-matrix(0,A,J);gb1<-gb2<-matrix(0,A,J)
    gks<-numeric(TT);gkc<-array(0,c(TT,J,max(1L,S)))
    etas<-mus<-array(0,c(A,TT,J));value<-0
    for(j in seq_len(J)) {
      et <- matrix(u$a[,j],A,TT)
      if(lc1) et<-et+tcrossprod(u$B1[,j],u$ks)
      else if(lc2) et<-et+tcrossprod(u$B1[,j],u$Kc[,j,1])+
                         tcrossprod(u$B2[,j],u$ks)
      else {
        et<-et+matrix(rep(u$ks,each=A),A,TT)
        for(h in seq_len(S)) et<-et+tcrossprod(V[,h],u$Kc[,j,h])
      }
      if(any(!is.finite(et)))
        stop("Non-finite linear predictor in multi-population likelihood")
      if(binom) {
        mu<-N[,,j]*plogis(et)
        value<-value-sum(D[,,j]*plogis(et,log.p=TRUE)+
          (N[,,j]-D[,,j])*plogis(et,lower.tail=FALSE,log.p=TRUE))/n_total
      } else {
        # At eta<=30 this is EXACTLY the Poisson log-likelihood and score.
        # Above 30, continue exp(eta) by its quadratic Taylor upper branch,
        # matching function and first derivative at 30. This prevents
        # numerical overflow without introducing a false zero-gradient plateau.
        excess<-pmax(et-30,0)
        ecap<-exp(pmin(et,30))
        mu<-E[,,j]*ecap*(1+excess)
        value<-value+sum(E[,,j]*ecap*(1+excess+0.5*excess^2)-
                         D[,,j]*et)/n_total
      }
      R<-(mu-D[,,j])/n_total
      ga[,j]<-rowSums(R)
      if(lc1) {
        gb1[,j]<-drop(R%*%u$ks)
        gks<-gks+drop(crossprod(u$B1[,j],R))
      } else if(lc2) {
        gb1[,j]<-drop(R%*%u$Kc[,j,1])
        gb2[,j]<-drop(R%*%u$ks)
        gkc[,j,1]<-drop(crossprod(u$B1[,j],R))
        gks<-gks+drop(crossprod(u$B2[,j],R))
      } else {
        gks<-gks+colSums(R)
        for(h in seq_len(S))gkc[,j,h]<-drop(crossprod(V[,h],R))
      }
      etas[,,j]<-et;mus[,,j]<-mu
    }
    grad<-if(lc1) c(as.vector(ga),gb1[mask1],gks[-TT]-gks[TT])
      else if(lc2) c(as.vector(ga),gb1[mask1],gb2[mask2],
        as.vector(sweep(gkc[-TT,,1,drop=FALSE],2,gkc[TT,,1],'-')),
        gks[-TT]-gks[TT])
      else {
        g<-if(cbd) gks else c(as.vector(ga),gks[-TT]-gks[TT])
        for(h in seq_len(S)) {
          tmp<-gkc[,,h]
          g<-c(g,if(cbd) as.vector(tmp) else
                   as.vector(sweep(tmp[-TT,,drop=FALSE],2,tmp[TT,],'-')))
        }
        g
      }
    if(!return_details) return(list(value=value,gradient=grad))
    list(value=value,gradient=grad,u=u,eta=etas,mu=mus)
  }
  list(start=initial,unpack=unpack,calc=calc,reference=ref,
       anchor1=anchor1,anchor2=anchor2,binom=binom,N=N,
       causes=names(singles),J=J,age_basis=V,model=model,
       sum_deaths=sum(D),death_weights=w)
}

fit_multi_group <- function(problem,D,E,maxit=3500L,restarts=1L,tol=2e-5) {
  starts<-list(problem$start)
  if(restarts>1L) {
    set.seed(SEED+round(problem$sum_deaths)%%10000L)
    for(i in 2:restarts)starts[[i]]<-problem$start+
      rnorm(length(problem$start),sd=0.003)
  }
  fit_one<-function(start) {
    f<-optim(start,function(z)problem$calc(z)$value,
       function(z)problem$calc(z)$gradient,method="L-BFGS-B",
       control=list(maxit=maxit,factr=1e6,pgtol=tol/3))
    if(max(abs(problem$calc(f$par)$gradient))>tol) {
      g<-try(nlminb(f$par,function(z)problem$calc(z)$value,
          function(z)problem$calc(z)$gradient,
          control=list(iter.max=maxit,eval.max=maxit*2L,
                       rel.tol=1e-12,x.tol=1e-10)),silent=TRUE)
      if(!inherits(g,"try-error") && is.finite(g$objective) &&
         g$objective <= f$value+1e-9)
        f<-list(par=g$par,value=g$objective,convergence=g$convergence,
                message="nlminb-refined")
    }
    f
  }
  fs<-lapply(starts,fit_one)
  f<-fs[[which.min(vapply(fs,function(x)x$value,0.0))]]
  out<-problem$calc(f$par,return_details=TRUE)
  gmax<-max(abs(out$gradient))
  u<-out$u
  rr<-array(0,dim(D))
  if(problem$binom) {
    qq<-out$mu/problem$N
    rr<-qq/(1-qq/2)
    dev<-2*sum(D*log(D/out$mu)+(problem$N-D)*
                log((problem$N-D)/(problem$N-out$mu)))
  } else {
    rr<-out$mu/E
    dev<-2*sum(D*log(D/out$mu)-D+out$mu)
  }
  metrics<-data.frame(model=problem$model,J=problem$J,
    n_parameters=length(f$par),objective=out$value,deviance=dev,
    max_scaled_gradient=gmax,optim_code=f$convergence,
    converged=(gmax<=tol),mape_percent=100*mean(abs(rr-D/E)/(D/E)),
    largest_cause_death_share=max(problem$death_weights),
    largest_cause_index=problem$reference,
    stringsAsFactors=FALSE)
  list(u=u,eta=out$eta,mu=out$mu,rate=rr,metrics=metrics,
       problem=problem,optimizer=f,gradient=out$gradient)
}

# Independent finite-difference check for the *joint* analytic score.
check_multi_gradient <- function(problem,eps=1e-6,ncheck=12L) {
  theta<-problem$start
  ind<-unique(round(seq(1,length(theta),length.out=min(ncheck,length(theta)))))
  score<-problem$calc(theta)$gradient
  err<-vapply(ind,function(i) {
    plus<-minus<-theta
    plus[i]<-plus[i]+eps;minus[i]<-minus[i]-eps
    fd<-(problem$calc(plus)$value-problem$calc(minus)$value)/(2*eps)
    abs(fd-score[i])/max(1,abs(fd),abs(score[i]))
  },0.0)
  data.frame(model=problem$model,parameters_checked=length(ind),
             maximum_relative_gradient_error=max(err),
             gradient_check_ok=max(err)<2e-4)
}
