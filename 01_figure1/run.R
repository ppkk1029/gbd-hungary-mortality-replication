# Figure 1：必须有 All causes 才能计算占全部死因的比例
path <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
stopifnot(file.exists(path))
gbd <- read.csv(path,stringsAsFactors=FALSE)
stopifnot(all(c("year","sex_name","age_name","cause_name","metric_name","val") %in% names(gbd)))
ages <- c("60-64 years","65-69 years","70-74 years",
          "75-79 years","80-84 years","85-89 years")
stopifnot(setequal(unique(gbd$age_name),ages))
stopifnot(!anyDuplicated(gbd[c("year","sex_name","age_name","cause_name","metric_name")]))

deaths <- subset(gbd,metric_name=="Number")
by_year <- aggregate(val ~ year + cause_name,data=deaths,FUN=sum)
train <- subset(by_year,year<=2015)
ranking <- aggregate(val ~ cause_name,data=train,FUN=sum)
ranking <- ranking[order(-ranking$val),]
top9 <- head(ranking$cause_name,9)
print(top9)

outdir <- "outputs/01_figure1"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(ranking,file.path(outdir,"cause_ranking_training_only.csv"),row.names=FALSE)
write.csv(by_year,file.path(outdir,"yearly_number_by_cause.csv"),row.names=FALSE)

allpath <- "data/raw/gbd_all_causes.csv"
if (!file.exists(allpath)) {
  message("尚缺 All causes CSV：只输出人数汇总，不生成 Figure 1 全死因占比图。")
} else {
  total <- read.csv(allpath,stringsAsFactors=FALSE)
  stopifnot(all(c("year","sex_name","age_name","cause_name","metric_name","val") %in% names(total)))
  total <- subset(total,metric_name=="Number" & cause_name=="All causes" &
                        age_name %in% ages & sex_name %in% c("Male","Female"))
  stopifnot(nrow(total)==44*6*2,
            setequal(unique(total$year),1980:2023),
            !anyDuplicated(total[c("year","sex_name","age_name")]))
  mat <- xtabs(val ~ year + cause_name,data=by_year)
  all_by_year <- aggregate(val ~ year,data=total,FUN=sum)
  all_by_year <- all_by_year[order(all_by_year$year),]
  stopifnot(identical(as.integer(rownames(mat)),as.integer(all_by_year$year)))
  nine <- mat[,top9,drop=FALSE]
  rest <- as.numeric(all_by_year$val)-rowSums(nine)
  stopifnot(all(rest >= -1e-5))
  plotdata <- cbind(nine,"All other causes"=pmax(rest,0))
  shares <- 100*plotdata/as.numeric(all_by_year$val)
  stopifnot(max(abs(rowSums(shares)-100))<1e-6)
  write.csv(data.frame(year=as.integer(rownames(shares)),shares,check.names=FALSE),
            file.path(outdir,"figure1_shares_all_causes.csv"),row.names=FALSE)
  years <- as.integer(rownames(shares))
  colors <- grDevices::hcl.colors(ncol(shares),"Set 3")
  png(file.path(outdir,"figure1_gbd.png"),width=1800,height=1000,res=160)
  par(mar=c(5,5,4,14),xpd=TRUE)
  plot(NA,xlim=range(years),ylim=c(0,100),xlab="Year",
       ylab="Share of all deaths (%)",main="GBD Hungary, ages 60-89")
  low <- rep(0,length(years))
  for (j in seq_len(ncol(shares))) {
    upper <- low+shares[,j]
    polygon(c(years,rev(years)),c(low,rev(upper)),
            col=colors[j],border=NA)
    low <- upper
  }
  legend("topright",inset=c(-0.39,0),legend=colnames(shares),
         fill=colors,cex=0.65,bty="n")
  dev.off()
  message("已生成 GBD Figure 1 对照图及占比 CSV。")
}
