# 请在项目根目录运行
path <- "data/raw/IHME-GBD_2023_DATA-b40ffc46-1.csv"
stopifnot(file.exists(path))
gbd <- read.csv(path, stringsAsFactors = FALSE)

required <- c("year","measure_name","location_name","sex_name","age_name",
              "cause_name","metric_name","val","lower","upper")
stopifnot(all(required %in% names(gbd)))

print(dim(gbd))
print(range(gbd$year))
print(setdiff(1980:2023, unique(gbd$year)))
print(sort(unique(gbd$age_name)))
print(table(gbd$metric_name))
print(length(unique(gbd$cause_name)))

keys <- c("year","location_name","sex_name","age_name","cause_name","metric_name")
n_duplicate <- sum(duplicated(gbd[keys]))
n_missing <- sum(is.na(gbd[required]))
n_negative <- sum(gbd$val < 0)
n_invalid_ui <- sum(gbd$lower > gbd$val | gbd$upper < gbd$val)
print(c(duplicate=n_duplicate,missing=n_missing,
        negative=n_negative,invalid_ui=n_invalid_ui))

ages <- c("60-64 years","65-69 years","70-74 years",
          "75-79 years","80-84 years","85-89 years")
stopifnot(nrow(gbd)==20064, ncol(gbd)==18,
          setequal(unique(gbd$year), 1980:2023),
          length(unique(gbd$year))==44,
          setequal(unique(gbd$sex_name), c("Male","Female")),
          setequal(unique(gbd$age_name), ages),
          length(unique(gbd$cause_name))==19,
          setequal(unique(gbd$metric_name), c("Number","Rate")),
          all(gbd$location_name=="Hungary"),
          all(gbd$measure_name=="Deaths"),
          n_duplicate==0, n_missing==0, n_negative==0, n_invalid_ui==0)
cat("数据基础检查通过：44年 × 2性别 × 6年龄组 × 19死因 × 2指标。\n")
