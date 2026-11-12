# §4.2 plots: author-style diagnostics, always labelled GBD analogues (NOT HCSO).
plot_device <- function(out,name,draw,width=12,height=6) {
  for(ext in c("png","svg")) {
    path <- file.path(out,paste0(name,".",ext))
    if(ext=="png") grDevices::png(path,width=width,height=height,units="in",res=150)
    else grDevices::svg(path,width=width,height=height)
    tryCatch(draw(),finally=grDevices::dev.off())
  }
}
PALETTE <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00",
             "#56B4E9","#804000","#555555","#A020F0")
plot_cluster_panels <- function(a,title,labs=NULL) {
  if(is.null(labs)) labs <- a$labels
  if(is.null(labs)) stop("未确认聚类K，不能绘制正式组内图")
  groups <- sort(unique(labs)); n <- length(groups)
  nr <- if(n<=3L) 1L else 2L; nc <- ceiling(n/nr)
  op <- par(mfrow=c(nr,nc),mar=c(4,4,3,1),oma=c(0,0,3,0))
  on.exit(par(op),add=TRUE)
  for(g in groups) {
    ix <- which(labs==g)
    matplot(YEARS,t(a$Z[ix,,drop=FALSE]),type="l",lty=1,
       col=PALETTE[ix],lwd=2,xlab="Year",ylab="Standardized kappa",
       ylim=range(a$Z),main=paste("Cluster",g))
    grid()
    legend("topleft",SHORT[ix],col=PALETTE[ix],lty=1,
           lwd=2,bty="n",cex=.75)
  }
  mtext(title,outer=TRUE,cex=1.05,font=2)
}
plot_alignment <- function(a,sex,number) {
  # Male: CVD and DIG, matches ICD IX vs XI at category level.
  # Female: CVD and GBD DKD, DOES NOT equal original ICD IX and IV.
  j <- if(sex=="Male") 3L else 7L
  x <- a$Z[1L,]; y <- a$Z[j,]
  fit <- dtw_distance(x,y,p=2,path=TRUE)
  k1 <- a$raw_kappa[1L,];k2 <- a$raw_kappa[j,]
  suffix <- if(sex=="Female") " (DKD is NOT ICD-IV)" else ""
  op <- par(mar=c(4,4,3,1),oma=c(0,0,3,0))
  on.exit({layout(1);par(op)},add=TRUE)
  layout(matrix(c(1,2,1,3),2,2,byrow=TRUE),widths=c(1,1.4))
  matplot(YEARS,cbind(k1,k2),type="l",lty=1,lwd=2,
    col=PALETTE[c(1L,j)],xlab="Year",ylab="Original fitted kappa(2)",
    main="Unstandardized period indices")
  legend("topright",SHORT[c(1L,j)],col=PALETTE[c(1L,j)],lty=1,bty="n")
  # The path was computed on standardized data. The heatmap is squared local cost.
  image(YEARS,YEARS,fit$cost,col=hcl.colors(60,"Viridis"),
    xlab="CVD year",ylab=paste(SHORT[j],"year"),
    main=sprintf("Squared-cost DTW p=2; d=%.3f",fit$distance))
  lines(YEARS[fit$path[,1L]],YEARS[fit$path[,2L]],col="white",lwd=2)
  plot(YEARS,x,type="l",lwd=2,col=PALETTE[1L],
       xlab="Year",ylab="Standardized kappa(2)",ylim=range(c(x,y)),
       main="Alignment computed on z-scores")
  lines(YEARS,y,col=PALETTE[j],lwd=2)
  mtext(paste0("Figure ",number," GBD ANALOGUE | ",sex," | Poisson LC2",
       suffix),outer=TRUE,font=2,cex=1.0)
}
plot_paper_figures <- function(analyses,out,fits=NULL) {
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  for(sex in c("Male","Female")) {
    number <- if(sex=="Male") 5L else 6L
    a <- analyses[[paste(sex,"Poisson LC2")]]
    plot_device(out,paste0("Figure_",number,"_GBD_",sex,"_Poisson_LC2"),
                function() plot_alignment(a,sex,number),13,7)
    plot_device(out,paste0("dendrograms_",sex),function(){
      op<-par(mfrow=c(4,2),mar=c(3,4,3,1))
      on.exit(par(op),add=TRUE)
      for(model in MODELS) {
        a <- analyses[[paste(sex,model)]]
        plot(a$hc,labels=SHORT,main=paste(sex,model),
             xlab="",sub="",cex=.8)
        if(!is.null(a$selection)) {
          rect.hclust(a$hc,k=a$selection$K[1L],border=2:8)
        }
      }
    },13,15)
    plot_device(out,paste0("distance_heatmaps_",sex),function(){
      op<-par(mfrow=c(4,2),mar=c(6,5,3,1))
      on.exit(par(op),add=TRUE)
      for(model in MODELS) {
        d <- analyses[[paste(sex,model)]]$d
        image(1:9,1:9,d,col=hcl.colors(60,"Viridis"),
              zlim=range(d),axes=FALSE,xlab="",ylab="",
              main=paste(sex,model))
        axis(1,1:9,SHORT,las=2,cex.axis=.7)
        axis(2,1:9,SHORT,las=2,cex.axis=.7)
      }
    },13,15)
  }
  # Before the author-style manual review, make only *provisional* preview panels.
  for(sex in c("Male","Female")) for(model in c("CBD","Binomial LC2")) {
    a <- analyses[[paste(sex,model)]]
    if(is.null(a$labels)) {
      s <- a$proposal
      labs <- a$candidate_labels[[paste(s$method[1L],s$K[1L],sep="::")]]
      name <- paste0("PREVIEW_NOT_FINAL_",sex,"_",gsub(" ","_",model))
      plot_device(out,name,function() plot_cluster_panels(a,
          paste("CVI-rank preview; NOT paper K |",sex,model),labs))
    } else {
      idx <- if(sex=="Male") {
        if(model=="CBD")"Figure_7_GBD_Male_CBD_Clusters"
        else "Figure_8_GBD_Male_Binomial_LC2"
      } else paste0("Appendix_GBD_Female_",gsub(" ","_",model))
      plot_device(out,idx,function() plot_cluster_panels(a,
        paste("GBD analogue; reviewed/explicit selection |",sex,model)))
    }
  }
  invisible(NULL)
}
