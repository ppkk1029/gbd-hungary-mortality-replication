# Standalone synthetic verification; DOES NOT imply actual GBD refitting.
# Run from repository root: Rscript 04_results/4_2_clusters/gbd_complete/selftest_section42.R
p <- "04_results/4_2_clusters/gbd_complete"
source("04_results/4_1_single_population/model_functions.R",local=TRUE,encoding="UTF-8")
source(file.path(p,"clustering_functions.R"),local=TRUE,encoding="UTF-8")
source(file.path(p,"kappa_import.R"),local=TRUE,encoding="UTF-8")
source(file.path(p,"plot_paper_figures.R"),local=TRUE,encoding="UTF-8")
source(file.path(p,"cluster_all_models.R"),local=TRUE,encoding="UTF-8")
check_cluster_packages()
set.seed(3104)
tt <- seq_len(36)
Z <- rbind(
 sin(tt/7)+tt/100,sin(tt/7)+tt/100+.05,
 cos(tt/8),cos(tt/8)+.05,
 tt/36,tt/36+.01,rev(tt)/36,
 sin(tt/3),cos(tt/3))
Z <- z_rows(Z)
rownames(Z) <- paste0("S",seq_len(nrow(Z)))
x <- c(1,2,3,4);y <- c(1,1,3,4)
d1 <- dtw_distance(x,y,p=2)
raw_d2 <- dtwclust::dtw_basic(x,y,norm="L2",normalize=FALSE,
  sqrt.dist=FALSE,step.pattern=dtw::symmetric2)
d2 <- sqrt(as.numeric(raw_d2)/(length(x)+length(y)))
cat("DTW cross-package squared-cost comparison: dtw=",d1,
    " dtwclust(raw->paper)=",d2,"\n")
stopifnot(isTRUE(all.equal(as.numeric(d1),d2,tolerance=1e-7)))
D <- distance_matrix(Z,p=2,check_impl=TRUE)
stopifnot(max(abs(D-t(D)))<1e-12,all(diag(D)==0))
e <- evaluate_clusters(Z,D,k_values=2:3,seed=42L)
stopifnot(nrow(e$scores)==4L,nrow(e$starts)==6L,
          length(e$labels)==4L,
          all(CVI_NAMES %in% colnames(e$scores)))
stopifnot(all(vapply(e$labels,length,0L)==9L))
# All seven advertised CVIs must actually compute at least one finite value on
# nondegenerate synthetic data; names alone do not establish working formulas.
finite_cvi <- vapply(CVI_NAMES,function(nm) sum(is.finite(e$scores[[nm]])),0L)
print(finite_cvi)
if(any(finite_cvi==0L))
  stop("原生CVI运行未产生有效值: ",paste(names(finite_cvi)[finite_cvi==0L],collapse=","),
       "\n",paste(e$scores$cvi_issues,collapse=" | "))
for (method in c("hierarchical_average","pam")) {
  fit <- cluster_fit(Z,D,method,k=2L,seed=42L)
  stopifnot(length(fit@cluster)==9L,length(unique(fit@cluster))==2L)
}
# Test the human K decision pathway and exact §4.3 membership CSV schema.
out <- tempfile(pattern="section42-test-");dir.create(out)
a <- list()
for(sex in c("Male","Female"))for(model in MODELS) {
  a[[paste(sex,model)]] <- list(
    Z=Z,raw_kappa=Z,d=D,hc=e$hc,
    labels=NULL,factor=if(endsWith(model,"LC2"))2L else 1L,
    scores=e$scores,candidate_labels=e$labels)
}
manual <- expand.grid(sex=c("Male","Female"),model=MODELS,
                      stringsAsFactors=FALSE)
manual$K <- 2L
manual$method <- "hierarchical_average"
manual$reason <- "synthetic smoke test, not a research conclusion"
manual_path <- file.path(out,"reviewed_selection.csv")
write.csv(manual,manual_path,row.names=FALSE)
a2 <- finalize_cluster_selections(a,out,manual_path)
# Plots include original and standardized kappa views and all 14 appendix panels.
plot_paper_figures(a2,out)
stopifnot(file.exists(file.path(out,"figures","Figure_7_GBD_Male_CBD.png")),
          file.exists(file.path(out,"figures","Appendix_standardized_Female_Plat3.svg")),
          file.exists(file.path(out,"within_cluster_dtw_alignments.csv")))
m <- read.csv(file.path(out,"cluster_membership.csv"))
s <- read.csv(file.path(out,"cluster_summary.csv"))
stopifnot(nrow(m)==126L,nrow(s)==14L,all(s$K==2L),
          all(s$selection_status=="researcher_reviewed"))
cat("§4.2 package, DTW equality, 7 CVIs, PAM starts and review-mode smoke test: PASS\n")

# Strictly round-trip all 126 fits, 13 factors per sex/cause, 1980-2015.
fake <- list()
for(sex in c("Male","Female"))for(cause in CAUSES)for(model in MODELS) {
  nf <- if(model=="Plat3")3L else if(endsWith(model,"LC1"))1L else 2L
  for(j in seq_len(nf)) fake[[length(fake)+1L]] <- data.frame(
     sex=sex,cause=cause,model=model,factor=j,year=YEARS,
     kappa=seq_along(YEARS)*j/36 + if(sex=="Male")0.1 else -0.1)
}
fake_path <- file.path(out,"test_single_period_indices.csv")
write.csv(do.call(rbind,fake),fake_path,row.names=FALSE)
loaded <- read_single_population_kappa(fake_path)
stopifnot(length(loaded)==126L,
          nrow(loaded[[key("Male",CAUSES[1],"Plat3")]]$K)==36L,
          ncol(loaded[[key("Male",CAUSES[1],"Plat3")]]$K)==3L)
cat("§4.2 synthetic verification: real DTW/CVI, 3 PAM starts, grouped figures and 126 cached fits: PASS\\n")
