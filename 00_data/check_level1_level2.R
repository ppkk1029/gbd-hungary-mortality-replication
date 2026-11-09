# 第 00 章｜验证 Level 1 与 Level 2 的死亡人数、死亡率一致性
p1 <- "data/raw/IHME-GBD_2023_DATA-97ee48e8-1.csv"
p2 <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
stopifnot(file.exists(p1),file.exists(p2))
l1 <- read.csv(p1,stringsAsFactors=FALSE)
l2 <- read.csv(p2,stringsAsFactors=FALSE)
stopifnot(nrow(l1)==3168,nrow(l2)==20064)

# 一个键表示相同的年份、性别、年龄组和计量方式
keys <- c("year","sex_name","age_name","metric_name")
stopifnot(!anyDuplicated(l1[c(keys,"cause_name")]))
stopifnot(!anyDuplicated(l2[c(keys,"cause_name")]))

# 同一个键下，把不同死因死亡人数（或同年龄组的 Rate）相加
s1 <- aggregate(val ~ year+sex_name+age_name+metric_name,data=l1,FUN=sum)
s2 <- aggregate(val ~ year+sex_name+age_name+metric_name,data=l2,FUN=sum)
names(s1)[5] <- "level1"
names(s2)[5] <- "level2"
audit <- merge(s1,s2,by=keys)
audit$gap <- audit$level1-audit$level2

# 44 年×2 性别×6 年龄组×2 Metrics = 1056 条对照
stopifnot(nrow(audit)==1056,max(abs(audit$gap))<1e-6)
cat("所有组别相等，最大绝对差额 =",max(abs(audit$gap)),"\n")

# 只有 Number 可以直接跨年龄与性别相加；Rate 不可这样合计
counts <- subset(audit,metric_name=="Number")
annual <- aggregate(cbind(level1,level2) ~ year,data=counts,FUN=sum)
annual$gap <- annual$level1-annual$level2
print(subset(annual,year %in% c(1980,2019,2021,2023)))

dir.create("outputs/00_data",recursive=TRUE,showWarnings=FALSE)
write.csv(audit,"outputs/00_data/level1_vs_level2_by_cell.csv",row.names=FALSE)
write.csv(annual,"outputs/00_data/level1_vs_level2_by_year.csv",row.names=FALSE)
