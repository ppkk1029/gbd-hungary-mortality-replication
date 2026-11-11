# §4.3 figure builder, base R graphics. No Python, no copying of author PNGs.
# Each non-singleton GBD cluster gets two images in the style of appendix folder 9.
# Figure 9 is the GBD analogue of the Plat3 common-versus-single-period comparison.

pretty_model <- function(m) switch(m,
  "Poisson LC1"="Poisson Carter-Lee",
  "Binomial LC1"="Binomial Carter-Lee",
  "Poisson LC2"="Poisson 2-factor LC (shared second index)",
  "Binomial LC2"="Binomial 2-factor LC (shared second index)",
  "CBD"="Binomial CBD (shared first index)",
  "Plat2"="Reduced Plat 2 (shared first index)",
  "Plat3"="Reduced Plat 3 (shared first index)")
pic_colors <- function(n) grDevices::hcl.colors(n,"Dark 3")

period_series <- function(fit,model,j) {
  f<-fit$u
  if(endsWith(model,"LC1")) return(list(common=f$ks,specific=NULL))
  if(endsWith(model,"LC2"))return(list(common=f$ks,specific=f$Kc[,j,1]))
  list(common=f$ks,specific=f$Kc[,j,,drop=FALSE])
}
z_compare <- function(v) (v-mean(v))/max(sd(v),1e-12)
period_overlay <- function(fit,singles,model,causes,main="",with_legend=TRUE) {
  # LC factor units are not identified across independent and joint fits.
  # Shape (z-score) comparison is meaningful; raw level comparison is NOT.
  factor <- if(endsWith(model,"LC2"))2L else 1L
  shared <- fit$u$ks
  local <- sapply(singles,function(x)x$K[,factor])
  if(is.null(dim(local)))local<-matrix(local,ncol=length(causes))
  normalized <- startsWith(model,"Poisson LC") || startsWith(model,"Binomial LC")
  if(normalized) {
    shared<-z_compare(shared)
    local<-apply(local,2,z_compare)
    ylabel<-"z-standardized period index (shape only)"
  } else ylabel<-"period index (fixed age basis)"
  yy<-range(c(shared,local),finite=TRUE)
  plot(YEARS,shared,type="l",lwd=2.9,col="black",ylim=yy,
       xlab="Calendar year",ylab=ylabel,main=main)
  cc<-pic_colors(length(causes))
  for(j in seq_along(causes))lines(YEARS,local[,j],col=cc[j],lwd=1.4)
  if(with_legend)legend("topright",legend=c("Joint common",causes),
    col=c("black",cc),lty=1,lwd=c(3,rep(1.5,length(cc))),
    cex=0.59,bg="white",ncol=if(length(causes)>5)2 else 1)
  data.frame(cause=causes,shape_corr=vapply(seq_along(causes),
      function(j)cor(shared,local[,j]),0.0),
    comparison_rmse=vapply(seq_along(causes),
      function(j)sqrt(mean((shared-local[,j])^2)),0.0),
    comparison_scale=if(normalized)"z-standardized" else "original fixed age scale")
}
draw_parameter_panels <- function(fit,model,causes,sex,group) {
  u<-fit$u;J<-length(causes);cc<-pic_colors(J)
  old<-par(no.readonly=TRUE);on.exit(par(old))
  par(mfrow=c(2,3),mar=c(3.2,3.9,3.5,1.1),oma=c(1,0,2,0),
      cex.axis=0.83,cex.lab=0.86)
  drawlines <- function(x,z,ttl,ylab,colors=cc) {
    if(is.vector(z)) z<-matrix(z,ncol=1)
    matplot(x,z,type="l",lty=1,lwd=1.4,col=colors,
            xlab=if(length(x)==length(AGES))"age midpoint" else "year",
            ylab=ylab,main=ttl)
  }
  if(model!="CBD") drawlines(AGES,u$a,"Cause-specific alpha","alpha")
  if(endsWith(model,"LC1")||endsWith(model,"LC2"))
    drawlines(AGES,u$B1,"Cause-specific beta1","beta1")
  if(endsWith(model,"LC2"))
    drawlines(AGES,u$B2,"Cause-specific beta2","beta2")
  drawlines(YEARS,u$ks,"Common period index","kappa common",colors="black")
  if(endsWith(model,"LC2"))
    drawlines(YEARS,u$Kc[,,1],"Specific first index","kappa1(j)")
  else if(model=="CBD"||model=="Plat2"||model=="Plat3")
    drawlines(YEARS,u$Kc[,,1],"Specific second index","kappa2(j)")
  if(model=="Plat3")
    drawlines(YEARS,u$Kc[,,2],"Specific third index","kappa3(j)")
  plot.new()
  legend("topleft",legend=causes,col=cc,lty=1,lwd=1.6,cex=.73,bty="n")
  mtext(paste0("GBD | ",sex," | ",pretty_model(model),
               " | cluster ",group),outer=TRUE,side=3,line=0.8,cex=1.12)
}

save_plot_png <- function(path,width=2100,height=1400,fun) {
  grDevices::png(path,width=width,height=height,res=170)
  tryCatch(fun(),finally=grDevices::dev.off())
}

make_appendix_figures <- function(fits,single_lookup,out) {
  d<-file.path(out,"figures","author_appendix_09_GBD_analogues")
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
  inventory<-list();similarity<-list()
  for(nm in names(fits)) {
    fit<-fits[[nm]];model<-fit$model;sex<-fit$sex;grp<-fit$group
    causes<-fit$causes;singles<-single_lookup[[nm]]
    stub<-paste0(sex,"_",gsub("[^A-Za-z0-9]","",model),
                 "_cluster",grp)
    comparison<-file.path(d,paste0(stub,"_common_vs_single.png"))
    params<-file.path(d,paste0(stub,"_joint_parameters.png"))
    sim<-NULL
    save_plot_png(comparison,fun=function() {
      sim<<-period_overlay(fit,singles,model,causes,
        main=paste("Common vs separate single-population indices:",nm))
    })
    save_plot_png(params,fun=function()
      draw_parameter_panels(fit,model,causes,sex,grp))
    sim$sex<-sex;sim$model<-model;sim$cluster<-grp
    similarity[[nm]]<-sim
    inventory[[nm]]<-data.frame(sex=sex,model=model,cluster=grp,
      causes=paste(causes,collapse="; "),original_class="appendix 9 image pair",
      comparison_png=file.path("figures","author_appendix_09_GBD_analogues",
                               basename(comparison)),
      parameters_png=file.path("figures","author_appendix_09_GBD_analogues",
                               basename(params)),stringsAsFactors=FALSE)
    message("§4.3 appendix pair: ",nm)
  }
  if(length(inventory)) {
    write_csv(do.call(rbind,inventory),file.path(out,"appendix_image_manifest.csv"))
    write_csv(do.call(rbind,similarity),file.path(out,"common_vs_single_diagnostics.csv"))
  }
  invisible(d)
}
draw_figure9_analogue <- function(fits,single_lookup,out) {
  chosen<-names(fits)[grepl("^Male Plat3 cluster|^Female Plat3 cluster",names(fits))]
  if(!length(chosen))return(invisible(NULL))
  chosen<-head(chosen,4)
  path<-file.path(out,"figures","Figure9_GBD_common_vs_single_Plat3.png")
  save_plot_png(path,width=2300,height=1700,fun=function() {
    old<-par(no.readonly=TRUE);on.exit(par(old))
    par(mfrow=c(2,2),mar=c(3.7,4.2,4.0,1.1),oma=c(1,1,2,0))
    for(nm in chosen) {
      f<-fits[[nm]]
      period_overlay(f,single_lookup[[nm]],f$model,f$causes,
         main=paste(f$sex,"cluster",f$group))
    }
    if(length(chosen)<4)for(i in seq_len(4-length(chosen)))plot.new()
    mtext("Figure 9 | GBD analogue | Plat3 joint common vs independent",
      side=3,outer=TRUE,line=0.6,cex=1.2)
  })
  invisible(path)
}

# Figure 10 (males) and Appendix Figure A4 (females) belong to
# §4.5 forecasting; these are generated as a separately labelled extension.
forecast_rw_drift <- function(z,h=2050L-max(YEARS),order=c(0L,1L,0L)) {
  tt<-seq_along(z)
  fit<-try(stats::arima(z,order=order,xreg=tt,method="ML"),silent=TRUE)
  if(!inherits(fit,"try-error")) {
    pr<-try(stats::predict(fit,n.ahead=h,newxreg=length(z)+seq_len(h)),silent=TRUE)
    if(!inherits(pr,"try-error") && all(is.finite(pr$pred)) &&
       all(is.finite(pr$se)))return(list(mean=as.numeric(pr$pred),
          se=as.numeric(pr$se),method=paste0("ARIMA(",paste(order,collapse=","),")+drift")))
  }
  dz<-diff(z);drift<-mean(dz);sdstep<-sd(dz)
  list(mean=tail(z,1)+seq_len(h)*drift,
       se=sqrt(seq_len(h))*sdstep,method="RW drift fallback")
}
plot_forecast_panel <- function(fit,main,yearend=2050L) {
  v<-fit$u$ks;f<-forecast_rw_drift(v,yearend-max(YEARS))
  yrs<-seq.int(max(YEARS)+1L,yearend)
  z80<-qnorm(.90);z95<-qnorm(.975)
  lo95<-f$mean-z95*f$se;hi95<-f$mean+z95*f$se
  lo80<-f$mean-z80*f$se;hi80<-f$mean+z80*f$se
  yrange<-range(c(v,lo95,hi95))
  plot(YEARS,v,type="n",xlim=range(c(YEARS,yrs)),ylim=yrange,
    main=main,xlab="Year",ylab="Common kappa")
  polygon(c(yrs,rev(yrs)),c(lo95,rev(hi95)),col="gray88",border=NA)
  polygon(c(yrs,rev(yrs)),c(lo80,rev(hi80)),col="gray72",border=NA)
  lines(YEARS,v,col="black",lwd=1.7)
  lines(c(tail(YEARS,1),yrs),c(tail(v,1),f$mean),col="firebrick",lwd=1.6)
  abline(v=max(YEARS),col="gray35",lty=3)
  f
}
draw_fig10_A4_extension <- function(fits,out) {
  for(sex in c("Male","Female")) {
    nm<-names(fits)[startsWith(names(fits),paste(sex,"Poisson LC1 cluster"))]
    if(!length(nm))next
    path<-file.path(out,"figures",
      if(sex=="Male")"Figure10_GBD_Male_Poisson_CarterLee_forecast_PREVIEW.png"
      else "FigureA4_GBD_Female_Poisson_CarterLee_forecast_PREVIEW.png")
    save_plot_png(path,width=1700,height=1250,fun=function() {
      old<-par(no.readonly=TRUE);on.exit(par(old))
      par(mfrow=c(1,length(nm)),mar=c(4,4.5,3.7,1.1),oma=c(1,0,2,0))
      for(k in nm)plot_forecast_panel(fits[[k]],main=k)
      mtext(paste(sex,"| common period index forecasts | 95% / 80% intervals"),
        outer=TRUE,side=3,line=0.5,cex=1.1)
    })
  }
}
