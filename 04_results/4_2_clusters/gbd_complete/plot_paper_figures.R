# §4.2 figures and online-appendix analogues (GBD only; NOT the author's HCSO figures).
# Only draw final grouped panels after all 14 K/method choices are reviewed.
PLOT_COLORS <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00",
                 "#56B4E9","#804000","#555555","#A020F0")
plot_device <- function(out,name,draw,width=12,height=6) {
  dir.create(out,showWarnings=FALSE,recursive=TRUE)
  for(ext in c("png","svg")) {
    f <- file.path(out,paste0(name,".",ext))
    if(ext=="png") grDevices::png(f,width=width,height=height,units="in",res=150)
    else grDevices::svg(f,width=width,height=height)
    tryCatch(draw(),finally=grDevices::dev.off())
  }
}
plot_cluster_panels <- function(a,title) {
  labs <- a$labels
  if(is.null(labs)) stop("分组尚未审查，不能绘制最终图")
  groups <- sort(unique(labs))
  n <- length(groups); nr <- if(n<=3L)1L else 2L; nc <- ceiling(n/nr)
  op <- par(mfrow=c(nr,nc),mar=c(4,4,3,1),oma=c(0,0,3,0))
  on.exit(par(op),add=TRUE)
  for(g in groups) {
    ix <- which(labs==g)
    Y <- t(a$Z[ix,,drop=FALSE])
    matplot(YEARS,Y,type="l",lty=1,col=PLOT_COLORS[ix],lwd=2,
       xlab="Year",ylab="Z-scored period factor",
       ylim=range(a$Z),main=paste0("Cluster ",g," (n=",length(ix),")"))
    grid()
    legend("topleft",SHORT[ix],col=PLOT_COLORS[ix],lty=1,
           lwd=2,bty="n",cex=.75)
  }
  mtext(title,outer=TRUE,cex=1,font=2)
}
plot_alignment <- function(a,sex,number) {
  # Male IX/XI analogue: CVD/DIG. Female IV/IX analogue: DKD/CVD;
  # DKD != endocrine/nutritional/metabolic ICD-IV; explicitly flag.
  i <- match("CVD",SHORT)
  j <- match(if(sex=="Male")"DIG" else "DKD",SHORT)
  stopifnot(!anyNA(c(i,j)))
  x <- a$Z[i,]; y <- a$Z[j,]
  fit <- dtw_distance(x,y,p=2,path=TRUE)
  k1 <- a$raw_kappa[i,]; k2 <- a$raw_kappa[j,]
  op <- par(mar=c(4,4,3,1),oma=c(0,0,3,0))
  on.exit({layout(1);par(op)},add=TRUE)
  layout(matrix(c(1,2,1,3),2,2,byrow=TRUE),widths=c(1,1.4))
  matplot(YEARS,cbind(k1,k2),type="l",lty=1,lwd=2,
          col=PLOT_COLORS[c(i,j)],xlab="Year",ylab="Fitted kappa(2)",
          main="Original period indices, not z-scores")
  legend("topright",SHORT[c(i,j)],col=PLOT_COLORS[c(i,j)],lty=1,bty="n")
  image(YEARS,YEARS,fit$cost,col=hcl.colors(60,"Viridis"),
        xlab=paste(SHORT[i],"year"),ylab=paste(SHORT[j],"year"),
        main=sprintf("Z-score DTW, p=2; distance %.3f",fit$distance))
  lines(YEARS[fit$path[,1L]],YEARS[fit$path[,2L]],col="white",lwd=2)
  matplot(YEARS,cbind(x,y),type="l",lty=1,lwd=2,
          col=PLOT_COLORS[c(i,j)],
          xlab="Year",ylab="Z-scored period factor",
          main="Series actually used for DTW")
  disclaimer <- if(sex=="Female")" | DKD is NOT ICD-IV" else ""
  mtext(paste0("Figure ",number," GBD analogue | ",sex," Poisson LC2",
               disclaimer),outer=TRUE,font=2,cex=.95)
}
plot_cvi_profiles <- function(a,sex,model) {
  op <- par(mfrow=c(4,2),mar=c(4,4,3,1),oma=c(0,0,3,0))
  on.exit(par(op),add=TRUE)
  for(metric in CVI_NAMES) {
    dat <- a$scores
    y <- dat[[metric]]; good <- is.finite(y)
    dir <- if(CVI_DIRECTION[[metric]]==1) "higher" else "lower"
    if(!any(good)) {
      plot.new();title(main=paste(metric,"no finite values"));next
    }
    plot(range(dat$K),range(y[good]),type="n",xlab="K",ylab=metric,
         main=paste(metric,"(",dir,"better)"),xaxt="n")
    axis(1,at=sort(unique(dat$K)))
    grid()
    for(idx in seq_along(c("hierarchical_average","pam"))) {
      method <- c("hierarchical_average","pam")[idx]
      z <- dat[dat$method==method,]
      lines(z$K,z[[metric]],type="b",lty=idx,pch=idx,lwd=1.5,
            col=PLOT_COLORS[idx])
    }
    legend("topright",c("Hierarchical average","PAM"),
           col=PLOT_COLORS[1:2],lty=1:2,pch=1:2,bty="n",cex=.65)
  }
  mtext(paste("Seven original dtwclust CVIs |",sex,model,
              "(not a combined score)"),outer=TRUE,cex=1.05,font=2)
}
export_aligned_pairs <- function(analyses,out) {
  rows <- list()
  for(nm in names(analyses)) {
    a <- analyses[[nm]]
    if(is.null(a$labels)) next
    sex <- sub(" .*","",nm)
    model <- substring(nm,nchar(sex)+2L)
    for(i in seq_len(nrow(a$Z)-1L)) for(j in seq.int(i+1L,nrow(a$Z))) {
      if(a$labels[i]!=a$labels[j]) next
      f <- dtw_distance(a$Z[i,],a$Z[j,],p=2,path=TRUE)
      rows[[length(rows)+1L]] <- data.frame(
        sex=sex,model=model,cluster=a$labels[i],
        cause1=SHORT[i],cause2=SHORT[j],
        year1=YEARS[f$path[,1L]],year2=YEARS[f$path[,2L]],
        local_squared_cost=f$cost[f$path],
        dtw_p2=f$distance)
    }
  }
  if(length(rows)) write_csv(do.call(rbind,rows),
         file.path(out,"within_cluster_dtw_alignments.csv"))
}
plot_paper_figures <- function(analyses,out,fits=NULL) {
  fig <- file.path(out,"figures")
  dir.create(fig,recursive=TRUE,showWarnings=FALSE)
  for(sex in c("Male","Female")) {
    number <- if(sex=="Male")5L else 6L
    a <- analyses[[paste(sex,"Poisson LC2")]]
    plot_device(fig,paste0("Figure_",number,"_GBD_",sex,"_Poisson_LC2"),
                function()plot_alignment(a,sex,number),13,7)
    plot_device(fig,paste0("dendrograms_",sex),function(){
      op <- par(mfrow=c(4,2),mar=c(3,4,3,1))
      on.exit(par(op),add=TRUE)
      for(model in MODELS) {
        z <- analyses[[paste(sex,model)]]
        plot(z$hc,labels=SHORT,main=paste(sex,model),
             xlab="",sub="",cex=.75)
        if(!is.null(z$selection) &&
           identical(z$selection$method[[1L]],"hierarchical_average"))
          rect.hclust(z$hc,k=z$selection$K[[1L]],border=2:8)
      }
    },13,15)
    plot_device(fig,paste0("distance_heatmaps_",sex),function(){
      op <- par(mfrow=c(4,2),mar=c(6,5,3,1))
      on.exit(par(op),add=TRUE)
      for(model in MODELS) {
        d <- analyses[[paste(sex,model)]]$d
        image(1:9,1:9,d,col=hcl.colors(60,"Viridis"),
              zlim=range(d),axes=FALSE,xlab="",ylab="",
              main=paste(sex,model))
        axis(1,at=1:9,labels=SHORT,las=2,cex.axis=.65)
        axis(2,at=1:9,labels=SHORT,las=2,cex.axis=.65)
      }
    },13,15)
    for(model in MODELS) {
      a <- analyses[[paste(sex,model)]]
      name <- paste0("CVI_",sex,"_",gsub(" ","_",model))
      plot_device(fig,name,function()plot_cvi_profiles(a,sex,model),12,13)
      if(!is.null(a$labels)) {
        core <- if(sex=="Male" && model=="CBD")"Figure_7_GBD_Male_CBD"
          else if(sex=="Male" && model=="Binomial LC2")"Figure_8_GBD_Male_Binomial_LC2"
          else paste0("Appendix_standardized_",sex,"_",gsub(" ","_",model))
        plot_device(fig,core,function()plot_cluster_panels(a,
            paste0("GBD | ",sex," | ",model," | reviewed K=",
                   length(unique(a$labels)))),12,7)
      }
    }
  }
  if(all(vapply(analyses,function(a)!is.null(a$labels),TRUE)))
    export_aligned_pairs(analyses,out)
  invisible(NULL)
}
