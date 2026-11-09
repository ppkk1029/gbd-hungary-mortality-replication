# Figure 1：GBD Level 2 内部构成；如有官方 All causes 则核验并生成正式全死因构成
path <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
stopifnot(file.exists(path))
gbd <- read.csv(path,stringsAsFactors=FALSE)
required <- c("year","sex_name","age_name","cause_name","metric_name","val")
stopifnot(all(required %in% names(gbd)))
ages <- c("60-64 years","65-69 years","70-74 years",
          "75-79 years","80-84 years","85-89 years")
stopifnot(setequal(unique(gbd$age_name),ages))
stopifnot(!anyDuplicated(gbd[c("year","sex_name","age_name","cause_name","metric_name")]))

# 1. 原始 Number 与按年、死因加总
deaths <- subset(gbd,metric_name=="Number")
stopifnot(length(unique(deaths$cause_name))==19)
by_year <- aggregate(val ~ year + cause_name,data=deaths,FUN=sum)
annual_level2 <- aggregate(val ~ year,data=by_year,FUN=sum)
annual_level2 <- annual_level2[order(annual_level2$year),]
names(annual_level2)[2] <- "D19"
print(subset(annual_level2,year==2019)) # 99182.03

# 2. 只使用训练期来固定图中展示的九个候选 GBD 类别
train <- subset(by_year,year<=2015)
ranking <- aggregate(val ~ cause_name,data=train,FUN=sum)
ranking <- ranking[order(-ranking$val),]
top9 <- head(ranking$cause_name,9)
print(top9)
outdir <- "outputs/01_figure1"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(ranking,file.path(outdir,"cause_ranking_training_only.csv"),row.names=FALSE)
write.csv(by_year,file.path(outdir,"yearly_number_by_cause.csv"),row.names=FALSE)

# 3. 获取年份×死因的死亡人数矩阵，年龄和男女已经加总
mat <- xtabs(val ~ year + cause_name,data=by_year)
years <- as.integer(rownames(mat))
stopifnot(identical(years,as.integer(annual_level2$year)))
nine <- mat[,top9,drop=FALSE]
allpath <- "data/raw/gbd_all_causes.csv"
verified <- file.exists(allpath)

# 4. 确认分母：优先用官方 All causes；否则只计算19类内部构成
if (verified) {
  total <- read.csv(allpath,stringsAsFactors=FALSE)
  stopifnot(all(required %in% names(total)))
  total <- subset(total,metric_name=="Number" & cause_name=="All causes" &
                        age_name %in% ages & sex_name %in% c("Male","Female"))
  stopifnot(nrow(total)==44*6*2,
            setequal(unique(total$year),1980:2023),
            !anyDuplicated(total[c("year","sex_name","age_name")]))
  official <- aggregate(val ~ year,data=total,FUN=sum)
  official <- official[order(official$year),]
  stopifnot(identical(as.integer(official$year),years))
  denominator <- official$val
  audit <- data.frame(year=years,D19=annual_level2$D19,
                      All_causes=denominator,
                      difference=denominator-annual_level2$D19)
  write.csv(audit,file.path(outdir,"level2_vs_all_causes_check.csv"),row.names=FALSE)
  print(head(audit))
  stopifnot(all(denominator>0),all(denominator-rowSums(nine)>=-1e-5))
  residual_label <- "All other causes"
  fig_name <- "figure1_gbd_all_causes.png"
  share_name <- "figure1_shares_all_causes.csv"
  plot_title <- "Hungary: share of all deaths (GBD All causes denominator)"
} else {
  denominator <- annual_level2$D19
  stopifnot(all(denominator>0))
  residual_label <- "Other selected Level 2 causes"
  fig_name <- "figure1_level2_internal_unverified.png"
  share_name <- "figure1_level2_internal_shares.csv"
  plot_title <- "Hungary: composition within 19 GBD Level 2 causes (UNVERIFIED)"
  message("没有 All causes：只能画19类内部构成，不能称为全部死因的比例。")
}
rest <- denominator - rowSums(nine)
stopifnot(all(rest >= -1e-5))
display <- cbind(nine,pmax(rest,0))
colnames(display)[ncol(display)] <- residual_label
shares <- 100 * display/denominator
stopifnot(max(abs(rowSums(shares)-100))<1e-5)
write.csv(data.frame(year=years,shares,check.names=FALSE),
          file.path(outdir,share_name),row.names=FALSE)

# 5. 原文同类的堆积面积图，未核验时明确标注分母差异
colours <- grDevices::hcl.colors(ncol(shares),"Set 3")
png(file.path(outdir,fig_name),width=1800,height=1000,res=160)
par(mar=c(5,5,4,14),xpd=TRUE)
plot(NA,xlim=range(years),ylim=c(0,100),
     xlab="Year",ylab="Share (%)",main=plot_title)
low <- rep(0,length(years))
for (j in seq_len(ncol(shares))) {
  upper <- low+shares[,j]
  polygon(c(years,rev(years)),c(low,rev(upper)),
          col=colours[j],border=NA)
  low <- upper
}
legend("topright",inset=c(-0.39,0),legend=colnames(shares),
       fill=colours,cex=0.65,bty="n")
dev.off()
cat("已输出:",file.path(outdir,fig_name),"\n")
