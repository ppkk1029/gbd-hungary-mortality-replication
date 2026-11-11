# 七种死亡率模型：只使用 R 自带 stats；不调用其它语言。
SEED <- 20261110L
MODELS <- c('Poisson LC1','Binomial LC1','Poisson LC2','Binomial LC2','CBD','Plat2','Plat3')
CAUSES <- c('Cardiovascular diseases','Neoplasms','Digestive diseases',
 'Chronic respiratory diseases','Neurological disorders','Unintentional injuries',
 'Diabetes and kidney diseases','Self-harm and interpersonal violence',
 'Respiratory infections and tuberculosis')
SHORT <- c('CVD','NEO','DIG','CRD','NEU','INJ','DKD','SHV','RIT')
AGES <- seq(62,87,5); YEARS <- 1980:2015
AGE_NAMES <- paste0(seq(60,85,5),'-',seq(64,89,5),' years')
write_csv <- function(x,path) write.csv(x,path,row.names=FALSE,na='',fileEncoding='UTF-8')
key <- function(sex,cause,model) paste(sex,cause,model,sep='::')
read_input <- function(path) {
 d <- read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
 req <- c('sex_name','cause_name','age_name','year','deaths_gbd','exposure_proxy')
 if(!all(req %in% names(d))) stop('输入缺列：',paste(setdiff(req,names(d)),collapse=', '))
 # 明确逐个验证单元，不能依赖输入行顺序；fitted_count 从不进入估计。
 if(nrow(d)!=3888L || anyDuplicated(d[req[1:4]])) stop('需3888个不重复观测单元')
 stopifnot(setequal(d$sex_name,c('Male','Female')),setequal(d$cause_name,CAUSES),
           setequal(d$age_name,AGE_NAMES),setequal(d$year,YEARS))
 if(any(!is.finite(d$deaths_gbd)|d$deaths_gbd<=0) ||
    any(!is.finite(d$exposure_proxy)|d$exposure_proxy<=0)) stop('死亡数/暴露量必须为有限正数')
 d
}
input_matrices <- function(d,sex,cause) {
 f <- d[d$sex_name==sex & d$cause_name==cause,]
 f <- f[order(match(f$age_name,AGE_NAMES),f$year),]
 stopifnot(nrow(f)==216L)
 list(D=matrix(f$deaths_gbd,6,36,byrow=TRUE),E=matrix(f$exposure_proxy,6,36,byrow=TRUE))
}
canonical <- function(a,B,K) {
 cen <- colMeans(K); a <- a + drop(B %*% cen); K <- sweep(K,2,cen)
 s <- svd(B %*% t(K)); r <- ncol(B)
 B <- s$u[,seq_len(r),drop=FALSE]
 K <- sweep(s$v[,seq_len(r),drop=FALSE],2,s$d[seq_len(r)],'*')
 for(j in seq_len(r)) { sg <- if(B[which.max(abs(B[,j])),j]>=0) 1 else -1
  B[,j] <- B[,j]*sg; K[,j] <- K[,j]*sg }
 list(a=a,B=B,K=K,s=s$d)
}
# 参数化：LC固定r个年龄载荷为单位阵，消除优化时旋转冗余；最终统一为SVD约定。
make_problem <- function(D,E,model) {
 A <- nrow(D); TT <- ncol(D); lc <- grepl('LC',model)
 binom <- startsWith(model,'Binomial') || model=='CBD'
 N <- if(binom) E+D/2 else E
 if(binom && any(N<=D)) stop('Binomial期初暴露量必须大于死亡人数')
 z <- if(binom) log(D/(N-D)) else log(D/E)
 if(lc) {
  r <- as.integer(substr(model,nchar(model),nchar(model)))
  a0 <- rowMeans(z); s <- svd(z-a0)
  U <- s$u[,seq_len(r),drop=FALSE]
  choices <- combn(A,r,simplify=FALSE)
  anchor <- choices[[which.max(vapply(choices,function(ix) abs(det(U[ix,,drop=FALSE])),0.0))]]
  M <- U[anchor,,drop=FALSE]; B0 <- U %*% solve(M)
  K0 <- sweep(s$v[,seq_len(r),drop=FALSE],2,s$d[seq_len(r)],'*') %*% t(M)
  free <- setdiff(seq_len(A),anchor); nb <- length(free)*r
  pack <- function(a,B,K) c(a,as.vector(B[free,,drop=FALSE]),as.vector(K[-TT,,drop=FALSE]))
  unpack <- function(v) {
   B <- matrix(0,A,r); B[anchor,] <- diag(r); B[free,] <- matrix(v[A+seq_len(nb)],length(free),r)
   K <- matrix(0,TT,r); K[-TT,] <- matrix(v[(A+nb+1):length(v)],TT-1,r); K[TT,] <- -colSums(K[-TT,,drop=FALSE])
   list(a=v[seq_len(A)],B=B,K=K)
  }
 } else {
  r <- if(model=='Plat3') 3L else 2L
  B0 <- cbind(1,(mean(AGES)-AGES)/10)
  if(model=='CBD') B0[,2] <- -B0[,2]
  if(r==3L) B0 <- cbind(B0,pmax((mean(AGES)-AGES)/10,0))
  a0 <- if(model=='CBD') rep(0,A) else rowMeans(z)
  K0 <- t(qr.solve(B0,z-a0))
  pack <- function(a,B,K) if(model=='CBD') as.vector(K) else c(a,as.vector(K[-TT,,drop=FALSE]))
  unpack <- function(v) {
   if(model=='CBD') return(list(a=rep(0,A),B=B0,K=matrix(v,TT,r)))
   K <- matrix(0,TT,r); K[-TT,] <- matrix(v[-seq_len(A)],TT-1,r); K[TT,] <- -colSums(K[-TT,,drop=FALSE])
   list(a=v[seq_len(A)],B=B0,K=K)
  }
 }
 total <- sum(D)
 calc <- function(v) {
  u <- unpack(v); eta <- u$a+u$B %*% t(u$K)
  if(binom) { mu <- N*plogis(eta)
   val <- sum(D*(log(D/N)-plogis(eta,log.p=TRUE)) +
              (N-D)*(log1p(-D/N)-plogis(eta,lower.tail=FALSE,log.p=TRUE)))/total
  } else { mu <- E*exp(eta); ld <- log(E/D)+eta; val <- sum(D*(expm1(ld)-ld))/total }
  res <- (mu-D)/total; ga <- rowSums(res); gk <- t(res) %*% u$B
  g <- if(lc) c(ga,as.vector((res %*% u$K)[free,,drop=FALSE]),
                 as.vector(sweep(gk[-TT,,drop=FALSE],2,gk[TT,]))) else if(model=='CBD') as.vector(gk) else
                  c(ga,as.vector(sweep(gk[-TT,,drop=FALSE],2,gk[TT,])))
  list(value=val,gradient=g)
 }
 list(start=pack(a0,B0,K0),unpack=unpack,calc=calc,lc=lc,binom=binom,N=N)
}
fit_model <- function(D,E,model) {
 p <- make_problem(D,E,model); set.seed(SEED)
 starts <- c(list(p$start),lapply(1:2,function(i) p$start+rnorm(length(p$start),0,.015)))
 fits <- lapply(starts,function(v) {
  f <- optim(v,function(x)p$calc(x)$value,function(x)p$calc(x)$gradient,
             method='BFGS',control=list(reltol=1e-16,maxit=5000))
  # 对近似平坦的似然进一步收紧梯度，而非仅信任optim的返回码。
  for(pass in 1:3) {
   if(max(abs(p$calc(f$par)$gradient))<1e-10) break
   g <- nlminb(f$par,function(x)p$calc(x)$value,function(x)p$calc(x)$gradient,
               control=list(rel.tol=1e-15,x.tol=1e-12,iter.max=3000,eval.max=6000))
   if(g$objective<=f$value+1e-13) f <- list(par=g$par,value=g$objective,convergence=g$convergence)
  }; f
 })
 best <- fits[[which.min(vapply(fits,`[[`,0.0,'value'))]]
 u <- p$unpack(best$par); eta <- u$a+u$B %*% t(u$K)
 grad <- max(abs(p$calc(best$par)$gradient))
 if(p$lc) u <- canonical(u$a,u$B,u$K) else u$s <- svd(u$B %*% t(u$K))$d
 mu <- if(p$binom) p$N*plogis(eta) else E*exp(eta)
 dev <- if(p$binom) 2*sum(D*log(D/mu)+(p$N-D)*log((p$N-D)/(p$N-mu))) else 2*sum(D*log(D/mu)-D+mu)
 rate <- if(p$binom) plogis(eta)/(1-plogis(eta)/2) else exp(eta)
 diag <- data.frame(deviance=dev,max_abs_scaled_gradient=grad,accepted=grad<1e-6,
  objective_range_3_starts=diff(range(vapply(fits,`[[`,0.0,'value'))),
  max_eta_difference_3_starts=max(vapply(fits,function(f) {v<-p$unpack(f$par);max(abs(v$a+v$B%*%t(v$K)-eta))},0.0)),
  mape_percent=100*mean(abs(rate-D/E)/(D/E)),factor2_factor1_ratio=u$s[2]/u$s[1],
  reconstruction_error=max(abs(u$a+u$B%*%t(u$K)-eta)))
 if(!diag$accepted) stop('未通过梯度核验：',model,' gradient=',grad)
 c(u,list(diag=diag))
}
