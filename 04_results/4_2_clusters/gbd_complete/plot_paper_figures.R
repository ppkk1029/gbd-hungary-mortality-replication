# 同时导出PNG和SVG；使用未四舍五入的R结果。图中数值由拟合和聚类生成。
plot_device <- function(out,name,draw,width=12,height=6) {
 for(ext in c('png','svg')) {
  path<-file.path(out,paste0(name,'.',ext))
  if(ext=='png') png(path,width=width,height=height,units='in',res=150) else svg(path,width=width,height=height)
  tryCatch(draw(),finally=dev.off())
 }
}
PALETTE <- c('#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#804000','#555555','#A020F0')
plot_cluster_panels <- function(a,title) {
 groups<-sort(unique(a$labels));par(mfrow=c(1,length(groups)),mar=c(4,4,3,1),oma=c(0,0,3,0))
 for(g in groups) {
  ix<-which(a$labels==g)
  matplot(YEARS,t(a$Z[ix,,drop=FALSE]),type='l',lty=1,col=PALETTE[ix],lwd=2,
          xlab='Year',ylab='Standardized period factor',ylim=range(a$Z),main=paste('Group',g))
  grid();legend('topleft',SHORT[ix],col=PALETTE[ix],lty=1,lwd=2,bty='n',cex=.8)
 }
 mtext(title,outer=TRUE,cex=1.1,font=2)
}
plot_alignment <- function(a,sex,number) {
 j<-if(sex=='Male') 3L else 7L;x<-a$Z[1,];y<-a$Z[j,];d<-dtw_distance(x,y,path=TRUE)
 layout(matrix(c(1,2,1,3),2,2,byrow=TRUE),widths=c(1,1.5));par(mar=c(4,4,3,1),oma=c(0,0,3,0))
 matplot(YEARS,cbind(x,y),type='l',lty=1,lwd=2,col=PALETTE[c(1,j)],xlab='Year',ylab='Standardized kappa(2)',main='Two period factors')
 legend('topright',SHORT[c(1,j)],col=PALETTE[c(1,j)],lty=1,bty='n')
 # x轴为第二死因、y轴为CVD；image的z必须转置以匹配轴方向。
 image(YEARS,YEARS,t(abs(outer(x,y,'-'))),col=hcl.colors(60,'Viridis'),
       xlab=paste(SHORT[j],'year'),ylab='CVD year',main=sprintf('Local |z1-z2|; DTW p2 = %.4f',d$distance))
 lines(YEARS[d$path[,2]],YEARS[d$path[,1]],col='white',lwd=2)
 matplot(YEARS,cbind(x,y+4),type='l',lty=1,col=PALETTE[c(1,j)],lwd=2,
         xlab='Year',ylab='Shifted z',main='Correspondences (second series +4)')
 for(k in seq(1,nrow(d$path),2)) segments(YEARS[d$path[k,1]],x[d$path[k,1]],YEARS[d$path[k,2]],y[d$path[k,2]]+4,col=adjustcolor('grey40',.4))
 mtext(paste('Figure',number,'| GBD Hungary |',sex,'| Poisson LC2 | 1980-2015'),outer=TRUE,font=2)
}
plot_paper_figures <- function(analyses,out,fits=NULL) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 for(sex in c('Male','Female')) {
  number<-if(sex=='Male') 5L else 6L;a<-analyses[[paste(sex,'Poisson LC2')]]
  plot_device(out,paste0('Figure_',number,'_GBD_',sex,'_Poisson_LC2'),function()plot_alignment(a,sex,number),13,7)
  plot_device(out,paste0('dendrograms_',sex),function() {
   par(mfrow=c(4,2),mar=c(3,4,3,1))
   for(model in MODELS) {a<-analyses[[paste(sex,model)]];plot(a$hc,labels=SHORT,main=paste(sex,model),xlab='',sub='',cex=.8)}
  },13,15)
  plot_device(out,paste0('heatmaps_',sex),function() {
   par(mfrow=c(2,4),mar=c(6,5,3,1))
   for(model in MODELS) {d<-analyses[[paste(sex,model)]]$d
    image(1:9,1:9,d,col=hcl.colors(60,'Viridis'),zlim=c(0,2),axes=FALSE,xlab='',ylab='',main=paste(sex,model))
    axis(1,1:9,SHORT,las=2,cex.axis=.7);axis(2,1:9,SHORT,las=2,cex.axis=.7)}
  },14,8)
  plot_device(out,paste0('clusters_',sex),function() {
   par(mfrow=c(4,2),mar=c(4,4,3,1))
   for(model in MODELS) {a<-analyses[[paste(sex,model)]]
    matplot(YEARS,t(a$Z),type='l',lty=1,col=PALETTE,xlab='Year',ylab='z',main=paste(sex,model))
    legend('topleft',paste0(SHORT,' G',a$labels),col=PALETTE,lty=1,ncol=3,cex=.6,bty='n')}
  },13,16)
 }
 for(item in list(c('CBD','7','Figure_7_GBD_Male_CBD_Clusters'),c('Binomial LC2','8','Figure_8_GBD_Male_Binomial_LC2'))) {
  a<-analyses[[paste('Male',item[1])]]
  plot_device(out,item[3],function()plot_cluster_panels(a,paste('Figure',item[2],'| Male',item[1],'| GBD 1980-2015')),13,5.5)
 }
 # 补充男女CBD/LC2各簇曲线，每个模型单独输出，避免硬编码最多3簇。
 for(sex in c('Male','Female')) for(model in c('CBD','Binomial LC2')) {
  a<-analyses[[paste(sex,model)]]
  plot_device(out,paste0('grouped_',sex,'_',gsub(' ','_',model)),function()plot_cluster_panels(a,paste(sex,model,'| standardized factors')))
 }
 invisible(NULL)
}
