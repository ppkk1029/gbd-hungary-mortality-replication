# DTW、聚类、七项CVI与敏感性；均为R实现，局部代价和归一化显式写出。
z_rows <- function(X) {
 s <- apply(X,1,sd)
 if(any(!is.finite(s)|s<1e-12)) stop('时期因子方差过小，不能标准化')
 sweep(sweep(X,1,rowMeans(X)),1,s,'/')
}
dtw_distance <- function(x,y,p=2,window=Inf,path=FALSE) {
 n <- length(x); m <- length(y); cost <- abs(outer(x,y,'-'))^p
 acc <- matrix(Inf,n,m); prev <- matrix(0L,n,m); acc[1,1] <- cost[1,1]
 for(i in seq_len(n)) for(j in seq_len(m)) {
  if((i==1 && j==1) || abs(i-j)>window) next
  opts <- c(if(i>1 && j>1) acc[i-1,j-1]+2*cost[i,j] else Inf,
            if(i>1) acc[i-1,j]+cost[i,j] else Inf,
            if(j>1) acc[i,j-1]+cost[i,j] else Inf)
  k <- which.min(opts); acc[i,j] <- opts[k]; prev[i,j] <- k
 }
 d <- (acc[n,m]/(n+m))^(1/p)
 if(!is.finite(d)) stop('无可行DTW路径')
 if(!path) return(d)
 pts <- matrix(c(n,m),1,2); i<-n;j<-m
 while(i>1 || j>1) {k<-prev[i,j];if(k %in% c(1,2)) i<-i-1L;if(k %in% c(1,3)) j<-j-1L;pts<-rbind(c(i,j),pts)}
 list(distance=d,cost=cost,acc=acc,path=pts)
}
distance_matrix <- function(Z,p=2,window=Inf) {
 n<-nrow(Z); d<-matrix(0,n,n)
 for(i in 2:n) for(j in seq_len(i-1)) d[i,j]<-d[j,i]<-dtw_distance(Z[i,],Z[j,],p,window)
 d
}
label_medoids <- function(d,med) {
 lab<-max.col(-d[,med,drop=FALSE],ties.method='first');lab[med]<-seq_along(med);lab
}
exact_medoids <- function(d,k) {
 opts<-combn(nrow(d),k,simplify=FALSE)
 vals<-vapply(opts,function(m) sum(apply(d[,m,drop=FALSE],1,min)),0.0)
 med<-opts[[which.min(vals)]];list(labels=label_medoids(d,med),medoids=med,objective=min(vals))
}
silhouette_mean <- function(d,lab) {
 vals<-vapply(seq_len(nrow(d)),function(i) {
  own<-which(lab==lab[i]);if(length(own)==1L) return(0)
  a<-mean(d[i,setdiff(own,i)]);other<-setdiff(unique(lab),lab[i])
  b<-min(vapply(other,function(k) mean(d[i,lab==k]),0.0))
  if(max(a,b)==0) 0 else (b-a)/max(a,b)
 },0.0);mean(vals)
}
ari <- function(a,b) {
 tab<-table(a,b);choose2<-function(x)x*(x-1)/2
 total<-choose2(sum(tab));x<-sum(choose2(tab));aa<-sum(choose2(rowSums(tab)));bb<-sum(choose2(colSums(tab)))
 expected<-aa*bb/total;den<-(aa+bb)/2-expected
 if(abs(den)<1e-14) return(1);(x-expected)/den
}
cluster_metrics <- function(d,lab,Z) {
 groups<-lapply(sort(unique(lab)),function(k) which(lab==k));n<-nrow(d);k<-length(groups)
 meds<-vapply(groups,function(g) g[which.min(rowSums(d[g,g,drop=FALSE]))],1L)
 diam<-max(vapply(groups,function(g) max(d[g,g,drop=FALSE]),0.0))
 pairs<-combn(k,2);sep<-min(apply(pairs,2,function(ab) min(d[groups[[ab[1]]],groups[[ab[2]]],drop=FALSE])))
 scatter<-vapply(seq_len(k),function(j) mean(d[groups[[j]],meds[j]]),0.0)
 mm<-d[meds,meds,drop=FALSE];rat<-outer(scatter,scatter,'+')/pmax(mm,1e-15);diag(rat)<- -Inf
 gm<-which.min(rowSums(d));within<-sum(vapply(seq_len(k),function(j) sum(d[groups[[j]],meds[j]]^2),0.0))
 between<-sum(lengths(groups)*d[meds,gm]^2)
 centers<-do.call(rbind,lapply(groups,function(g) colMeans(Z[g,,drop=FALSE])));grand<-colMeans(Z)
 we<-sum(vapply(seq_len(k),function(j) sum(sweep(Z[groups[[j]],,drop=FALSE],2,centers[j,])^2),0.0))
 be<-sum(lengths(groups)*rowSums(sweep(centers,2,grand)^2))
 cop<-sum(vapply(seq_len(k),function(j) {
  g<-groups[[j]];den<-min(apply(d[setdiff(seq_len(n),g),g,drop=FALSE],2,max))
  sum(d[g,meds[j]])/max(den,1e-15)
 },0.0))/n
 dbstar<-mean(vapply(seq_len(k),function(j) max(scatter[j]+scatter[-j])/max(min(mm[j,-j]),1e-15),0.0))
 lb<-sum(lengths(groups)*d[meds,gm]);lw<-sum(lengths(groups)*scatter)
 sf<- -expm1(-exp(min(6,max(-700,lb/(n*k)-sum(scatter)))))
 data.frame(silhouette=silhouette_mean(d,lab),dunn=sep/max(diam,1e-15),
  DB_medoid=mean(apply(rat,1,max)),DBstar_medoid=dbstar,CH_dtwclust_medoid=(n-k)/(k-1)*lb/lw,
  SF_medoid=sf,CH_squared_medoid=(between/(k-1))/(within/(n-k)),
  CH_euclidean=(be/(k-1))/(we/(n-k)),COP_medoid=cop,objective=lw)
}
choose_clusters <- function(d,Z) {
 hc<-hclust(as.dist(d),method='average');rows<-list();labs<-list()
 for(k in 2:(nrow(d)-1)) for(method in c('average','kmedoids_exact')) {
  lab<-if(method=='average') cutree(hc,k=k) else exact_medoids(d,k)$labels
  nm<-paste(method,k);labs[[nm]]<-lab
  rows[[nm]]<-cbind(data.frame(method=method,K=k),cluster_metrics(d,lab,Z))
 }
 tab<-do.call(rbind,rows);tab<-tab[order(-tab$silhouette,tab$K,tab$method),];rownames(tab)<-NULL
 best<-tab[1,];list(scores=tab,labels=labs[[paste(best$method,best$K)]],hc=hc)
}
pam_check <- function(d,k) {
 global<-exact_medoids(d,k)$objective
 do.call(rbind,lapply(1:3,function(start) {
  med<-sort(sample.int(nrow(d),k));obj<-sum(apply(d[,med,drop=FALSE],1,min))
  repeat {
   best<-obj;bm<-med
   for(a in med) for(b in setdiff(seq_len(nrow(d)),med)) {
    cand<-sort(c(setdiff(med,a),b));v<-sum(apply(d[,cand,drop=FALSE],1,min))
    if(v<best-1e-12) {best<-v;bm<-cand}
   }
   if(best>=obj-1e-12) break
   obj<-best;med<-bm
  }
  data.frame(K=k,start=start,objective=obj,global_objective=global,gap=obj-global)
 }))
}
