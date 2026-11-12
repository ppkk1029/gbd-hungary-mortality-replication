# Section 4.2 input adapter: accept previously fitted SINGLE-population
# period indices from §4.1, avoiding redundant refits of 126 models.
# This file does NOT estimate or modify period factors.
read_single_population_kappa <- function(path) {
  df <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
  required <- c("sex","cause","model","factor","year","kappa")
  if(!all(required %in% names(df)))
    stop("时期因子CSV缺少列: ",paste(setdiff(required,names(df)),collapse=","))
  if(anyNA(df[required]) || any(!is.finite(df$kappa)) ||
     any(!is.finite(df$factor)) || any(!is.finite(df$year)))
    stop("时期因子包含NA或无穷值")
  if(anyDuplicated(df[required[1:5]]))
    stop("时期因子sex/cause/model/factor/year键存在重复行")
  if(!setequal(unique(df$sex),c("Male","Female")) ||
     !setequal(unique(df$cause),CAUSES) ||
     !setequal(unique(df$model),MODELS) ||
     !setequal(unique(df$year),YEARS))
    stop("缓存的性别/死因/模型/年份与当前GBD定义不一致")
  factors <- function(model) if(model=="Plat3")3L else
    if(endsWith(model,"LC1"))1L else 2L
  required_n <- length(YEARS)*length(CAUSES)*2L*
    sum(vapply(MODELS,factors,0L))
  if(nrow(df)!=required_n) stop("时期因子记录数错误：应为",
                                 required_n,"，实为",nrow(df))
  fits <- list()
  for(sex in c("Male","Female")) for(cause in CAUSES) for(model in MODELS) {
    z <- df[df$sex==sex & df$cause==cause & df$model==model,]
    nf <- factors(model)
    if(nrow(z)!=length(YEARS)*nf ||
       !setequal(unique(z$factor),seq_len(nf)))
      stop("时期因子缺失或因子数错误：",sex,"/",cause,"/",model)
    K <- matrix(NA_real_,length(YEARS),nf)
    for(j in seq_len(nf)) {
      u <- z[z$factor==j,]
      if(nrow(u)!=length(YEARS) || !setequal(u$year,YEARS))
        stop("缺少年份：",sex,"/",cause,"/",model,"/",j)
      K[,j] <- u$kappa[match(YEARS,u$year)]
    }
    fits[[key(sex,cause,model)]] <- list(K=K)
  }
  if(length(fits)!=126L)stop("导入拟合结果不完整")
  fits
}
copy_single_population_outputs <- function(kappa_path,out) {
  parent <- dirname(kappa_path)
  for(nm in c("kappa_all_models.csv",
              "parameters_all_models.csv","fit_diagnostics_all_models.csv")) {
    src <- file.path(parent,nm)
    if(file.exists(src)) {
      ok <- file.copy(src,file.path(out,nm),overwrite=FALSE)
      if(!ok) stop("无法复制§4.1模型参数: ",nm)
    }
  }
}
