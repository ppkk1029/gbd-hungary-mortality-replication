# 第01章｜Figure 1：全部死亡人数中的主要死因构成（GBD方法复现）
# 从同一GBD2023的Level 1三大类，得到官方层级下的全部死亡人数
# 从Level 2得到各死因死亡人数；首次分析必须先核对L1与L2相等
p1 <- "data/raw/IHME-GBD_2023_DATA-97ee48e8-1.csv"
p2 <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
stopifnot(file.exists(p1),file.exists(p2))
l1 <- read.csv(p1,stringsAsFactors=FALSE)
l2 <- read.csv(p2,stringsAsFactors=FALSE)

# 1. 同一年龄、年份、性别、指标下，两层级的和应一致
keys <- c("year","sex_name","age_name","metric_name")
a <- aggregate(val ~ year+sex_name+age_name+metric_name,data=l1,FUN=sum)
b <- aggregate(val ~ year+sex_name+age_name+metric_name,data=l2,FUN=sum)
names(a)[5] <- "L1"
names(b)[5] <- "L2"
verified <- merge(a,b,by=keys)
stopifnot(nrow(verified)==1056,max(abs(verified$L1-verified$L2))<1e-6)

# 2. Figure 1只用Number（死亡人数），不使用Rate
n1 <- subset(l1,metric_name=="Number")
n2 <- subset(l2,metric_name=="Number")
total <- aggregate(val ~ year,data=n1,FUN=sum)
total <- total[order(total$year),]
bycause <- aggregate(val ~ year+cause_name,data=n2,FUN=sum)

# 3. 仅用训练期(1980-2015)选出9类死亡人数最多的GBD死因
train <- subset(bycause,year<=2015)
rank <- aggregate(val ~ cause_name,data=train,FUN=sum)
rank <- rank[order(-rank$val),]
chosen <- head(rank$cause_name,9)
print(chosen)

# 4. 每年9类死亡数加总；其余所有死因统一归为All other causes
wide <- xtabs(val ~ year+cause_name,data=bycause)
years <- as.integer(rownames(wide))
stopifnot(identical(years,as.integer(total$year)))
nine <- wide[,chosen,drop=FALSE]
others <- total$val-rowSums(nine)
stopifnot(all(others>=-1e-6))
stack <- cbind(nine,"All other causes"=pmax(others,0))
shares <- sweep(stack,1,total$val,"/")*100
stopifnot(max(abs(rowSums(shares)-100))<1e-6)

dir.create("outputs/01_figure1",recursive=TRUE,showWarnings=FALSE)
write.csv(rank,"outputs/01_figure1/train_cause_ranking.csv",row.names=FALSE)
write.csv(data.frame(year=years,shares,check.names=FALSE),
          "outputs/01_figure1/figure1_shares.csv",row.names=FALSE)

# 5. 画与原文Figure 1形式相同的堆积面积图
png("outputs/01_figure1/figure1_gbd.png",width=1800,height=1000,res=160)
par(mar=c(5,5,4,14),xpd=TRUE)
plot(NA,xlim=range(years),ylim=c(0,100),xlab="Year",ylab="Share of all deaths (%)",
     main="GBD 2023 Hungary, ages 60-89 (different causes from Varga)")
cols <- grDevices::hcl.colors(ncol(shares),"Set 3")
low <- rep(0,length(years))
for (j in seq_len(ncol(shares))) {
  up <- low+shares[,j]
  polygon(c(years,rev(years)),c(low,rev(up)),col=cols[j],border=NA)
  low <- up
}
legend("topright",inset=c(-0.39,0),legend=colnames(shares),
       fill=cols,cex=0.65,bty="n")
dev.off()

# 6. 直接输出与原文可对照的前9类覆盖率（不是原文那9类ICD）
coverage <- data.frame(year=years,all_deaths=total$val,
                       top9_deaths=rowSums(nine),
                       top9_share=100*rowSums(nine)/total$val)
write.csv(coverage,"outputs/01_figure1/top9_coverage.csv",row.names=FALSE)
print(subset(coverage,year %in% c(2019,2020,2021)))
cat("已生成 Figure 1 数据表、面积图和9类覆盖率。\n")
