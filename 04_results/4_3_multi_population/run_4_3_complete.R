# §4.3 pure R: fit each cluster jointly (NEVER average kappa as a "joint fit").
# RStudio: open this file and click Source, then choose the SAME fitted_rates.csv
# that was used in §4.2. Run default: 30 clusters = 27 joint + 3 singleton.
# Command: Rscript run_4_3_complete.R /path/fitted_rates.csv /path/new_output

.section43_dir <- local({
  z<-lapply(sys.frames(),function(e)e$ofile)
  z<-Filter(function(e)!is.null(e)&&length(e)==1L,z)
  arg<-grep("^--file=",commandArgs(),value=TRUE)
  if(length(z)) dirname(normalizePath(tail(z,1)[[1]],mustWork=TRUE))
  else if(length(arg)) dirname(normalizePath(sub("^--file=","",arg[1]),mustWork=TRUE))
  else getwd()
})
.section42_dir<-file.path(dirname(.section43_dir),"4_2_clusters","gbd_complete")
.section41_dir <- file.path(dirname(.section43_dir),"4_1_single_population")
if(!file.exists(file.path(.section41_dir,"model_functions.R")))
  stop("找不到§4.1的R模型函数；请下载整个仓库并保持目录结构。")
source(file.path(.section41_dir,"model_functions.R"),local=TRUE,encoding="UTF-8")
source(file.path(.section43_dir,"multi_model_functions.R"),local=TRUE,encoding="UTF-8")
source(file.path(.section43_dir,"plot_appendix.R"),local=TRUE,encoding="UTF-8")
source(file.path(.section43_dir,"author_appendix_inventory.R"),local=TRUE,encoding="UTF-8")

run_section43 <- function(input=NULL,out=NULL,source42=NULL,
                          restarts=1L,strict=TRUE,forecast_preview=TRUE) {
  if(is.null(input)) {
    if(!interactive())stop("命令行需提供 fitted_rates.csv")
    message("选择与§4.2完全相同的原始GBD拟合输入 CSV")
    input<-file.choose()
  }
  input<-normalizePath(input,mustWork=TRUE)
  if(is.null(source42))source42<-file.path(.section42_dir,"results")
  source42<-normalizePath(source42,mustWork=TRUE)
  if(is.null(out))out<-file.path(dirname(input),
    paste0("Section_4_3_R_",format(Sys.time(),"%Y%m%d_%H%M%S")))
  if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
    stop("输出目录非空，禁止覆盖：",out)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(out,"figures"),showWarnings=FALSE)
  dir.create(file.path(out,"rds"),showWarnings=FALSE)
  out<-normalizePath(out,mustWork=TRUE)
  must<-c("parameters_all_models.csv","kappa_all_models.csv","cluster_membership.csv")
  if(!all(file.exists(file.path(source42,must))))
    stop("§4.2单人口拟合与聚类结果CSV缺失: ",
         paste(must[!file.exists(file.path(source42,must))],collapse=", "))
  raw<-read_input(input)
  pp<-read.csv(file.path(source42,must[1]),stringsAsFactors=FALSE)
  kk<-read.csv(file.path(source42,must[2]),stringsAsFactors=FALSE)
  cl<-read.csv(file.path(source42,must[3]),stringsAsFactors=FALSE)
  expect<-expand.grid(sex=c("Male","Female"),model=MODELS,
                      cause=CAUSES,stringsAsFactors=FALSE)
  tag<-function(d)paste(d$sex,d$model,d$cause,sep="::")
  if(nrow(cl)!=nrow(expect) ||
     !setequal(tag(cl),tag(expect)) || anyDuplicated(tag(cl)) ||
     anyNA(cl$cluster))stop("§4.2 cluster membership不是14组×9死因的完整结果")
  if(any(!cl$model %in% MODELS))stop("意外模型名")
  if(!all(pp$sex %in% c("Male","Female")))stop("§4.2 source inconsistent")
  message("Input verified: ",nrow(raw)," data cells; clusters: ",
          length(unique(paste(cl$sex,cl$model,cl$cluster))))
  fits<-list();single_lookup<-list();summ<-list();causemetrics<-list()
  plist<-list();indices<-list();checks<-list();nm_i<-0L
  for(sex in c("Male","Female"))for(model in MODELS) {
    subset_cl<-cl[cl$sex==sex & cl$model==model, ]
    for(g in sort(unique(subset_cl$cluster))) {
      causes<-intersect(CAUSES,subset_cl$cause[subset_cl$cluster==g])
      nm<-paste(sex,model,"cluster",g)
      J<-length(causes)
      singles<-setNames(lapply(causes,function(c)
        get_single_fit(pp,kk,sex,c,model)),causes)
      Ds<-Es<-array(0,c(length(AGES),length(YEARS),J))
      for(j in seq_len(J)) {
        z<-input_matrices(raw,sex,causes[j]);Ds[,,j]<-z$D;Es[,,j]<-z$E
      }
      weights<-vapply(seq_len(J),function(j)sum(Ds[,,j]),0.0)
      if(J==1L) {
        summ[[nm]]<-data.frame(sex=sex,model=model,cluster=g,J=1L,
          status="SINGLETON_NO_SHARING",converged=NA,gradient=NA_real_,
          n_parameters=NA_integer_,mape_percent=NA_real_,
          largest_death_share=1,causes=causes,stringsAsFactors=FALSE)
        message("Skip singleton; retain independent §4.2 fit: ",nm)
        next
      }
      q<-make_multi_problem(Ds,Es,model,singles)
      ck<-check_multi_gradient(q)
      checks[[nm]]<-cbind(data.frame(sex=sex,model=model,cluster=g),ck)
      if(!ck$gradient_check_ok)stop("联合解析梯度未通过检查：",nm)
      f<-fit_multi_group(q,Ds,Es,restarts=restarts)
      f$sex<-sex;f$model<-model;f$group<-g;f$causes<-causes
      fitgood<-isTRUE(f$metrics$converged)
      status<-if(fitgood)"CONVERGED" else "NOT_CONVERGED_CHECK_FIT"
      summ[[nm]]<-data.frame(sex=sex,model=model,cluster=g,J=J,
        status=status,converged=fitgood,
        gradient=f$metrics$max_scaled_gradient,
        n_parameters=f$metrics$n_parameters,
        mape_percent=f$metrics$mape_percent,
        largest_death_share=max(weights)/sum(weights),
        causes=paste(causes,collapse="; "),stringsAsFactors=FALSE)
      message(nm, " | ",status," | gradient=",format(f$metrics$max_scaled_gradient,digits=3))
      if(!fitgood && strict) {
        write_csv(do.call(rbind,summ),file.path(out,"fit_summary_PARTIAL.csv"))
        stop("该组未达阈值，停止发布未核验结果: ",nm,
             ". 可调大 maxit 或为调试显式传 strict=FALSE")
      }
      fits[[nm]]<-f;single_lookup[[nm]]<-singles
      saveRDS(f,file.path(out,"rds",paste0(gsub("[^A-Za-z0-9]","",nm),".rds")))
      for(j in seq_len(J)) {
        s<-singles[[j]]
        eta_s<-s$a+s$B%*%t(s$K)
        mu_s<-if(q$binom)q$N[,,j]*plogis(eta_s) else Es[,,j]*exp(eta_s)
        rate_s<-if(q$binom) (mu_s/q$N[,,j])/(1-mu_s/q$N[,,j]/2)
                   else mu_s/Es[,,j]
        obs<-Ds[,,j]/Es[,,j]
        nm_i<-nm_i+1L
        causemetrics[[nm_i]]<-data.frame(
          sex=sex,model=model,cluster=g,cause=causes[j],
          deaths=sum(Ds[,,j]),death_share_in_cluster=weights[j]/sum(weights),
          MAPE_single_percent=100*mean(abs(rate_s-obs)/obs),
          MAPE_multi_percent=100*mean(abs(f$rate[,,j]-obs)/obs),
          MAPE_delta_multi_minus_single=100*mean(abs(f$rate[,,j]-obs)/obs)-
                                        100*mean(abs(rate_s-obs)/obs),
          status=status,stringsAsFactors=FALSE)
        u<-f$u;beta1<-if(endsWith(model,"LC1")||endsWith(model,"LC2"))
          u$B1[,j] else rep(1,length(AGES))
        beta2<-if(endsWith(model,"LC2"))u$B2[,j] else
          if(model=="CBD"||model=="Plat2"||model=="Plat3")
            q$age_basis[,1] else rep(NA_real_,length(AGES))
        beta3<-if(model=="Plat3")q$age_basis[,2] else
          rep(NA_real_,length(AGES))
        plist[[length(plist)+1L]]<-data.frame(sex=sex,model=model,
          cluster=g,cause=causes[j],age=AGES,alpha=u$a[,j],
          beta1=beta1,beta2=beta2,beta3=beta3)
        if(endsWith(model,"LC2")||model=="CBD"||model=="Plat2"||model=="Plat3"){
          nc<-dim(u$Kc)[3]
          for(h in seq_len(nc))indices[[length(indices)+1L]]<-
            data.frame(sex=sex,model=model,cluster=g,
              index_type="cause_specific",cause=causes[j],
              factor=if(endsWith(model,"LC2"))1L else h+1L,
              year=YEARS,kappa=u$Kc[,j,h])
        }
      }
      indices[[length(indices)+1L]]<-data.frame(
        sex=sex,model=model,cluster=g,index_type="shared",
        cause="",factor=if(endsWith(model,"LC2"))2L else 1L,
        year=YEARS,kappa=f$u$ks)
    }
  }
  stopifnot(length(summ)==30L)
  write_csv(do.call(rbind,summ),file.path(out,"fit_summary.csv"))
  write_csv(do.call(rbind,causemetrics),file.path(out,"cause_fit_metrics.csv"))
  write_csv(do.call(rbind,plist),file.path(out,"joint_age_parameters.csv"))
  write_csv(do.call(rbind,indices),file.path(out,"joint_period_indices.csv"))
  write_csv(do.call(rbind,checks),file.path(out,"finite_difference_gradient_checks.csv"))
  author_appendix_inventory(file.path(out,"author_original_appendix_52_files.csv"))
  make_appendix_figures(fits,single_lookup,out)
  draw_figure9_analogue(fits,single_lookup,out)
  if(forecast_preview)draw_fig10_A4_extension(fits,out)
  lines<-c("§4.3: joint maximum pseudo-likelihood, R base functions only",
    paste("R version:",getRversion()),
    paste("GBD input MD5:",unname(tools::md5sum(input))),
    paste("§4.2 parameters MD5:",unname(tools::md5sum(file.path(source42,must[1])))),
    paste("§4.2 cluster CSV MD5:",unname(tools::md5sum(file.path(source42,must[3])))),
    "1980-2015, 60-89 y (six ages), nine GBD causes, male/female separately",
    "Figure9 & appendix9: actual joint estimator output, not averaged single kappa",
    "LC1 shares first index; LC2 shares second; CBD/Plat2/Plat3 share first",
    "30 model-sex-clusters; 27 joint groups, 3 unshared singletons",
    "Binomial n=GBD_population_proxy + 0.5*GBD_deaths, not official raw E0",
    "Criterion: scaled score gradient, first-order necessary condition only",
    "2016-2019 held out; no predictive accuracy claim in this section",
    "Figure 10 / A4 previews are §4.5 extensions, not §4.3 fitting results",
    "All author images are external references, no source pixels re-used")
  writeLines(lines,file.path(out,"MANIFEST.txt"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  message("§4.3 done: ",out)
  print(read.csv(file.path(out,"fit_summary.csv"))[,c("sex","model","cluster","J","status")],row.names=FALSE)
  invisible(out)
}
if(isTRUE(getOption("section43.autorun",TRUE))) {
  aa<-commandArgs(trailingOnly=TRUE)
  run_section43(input=if(length(aa)>=1)aa[1] else NULL,
                out=if(length(aa)>=2)aa[2] else NULL)
}
