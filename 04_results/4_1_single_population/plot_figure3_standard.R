# §4.1 / Figure 3 — visual-only split of one extreme female beta curve
# Inputs: alpha.csv, beta.csv, kappa.csv, selected_causes.csv,
# from THE SAME single-factor Poisson LC fit (1980–2015).
# This script NEVER refits a model or removes a cause from the data.
# Run: source("04_results/4_1_single_population/plot_figure3_standard.R")
#
# The standard Lee–Carter constraints are sum_x beta[x]=1, mean_t kappa[t]=0.
# If a beta vector has positive/negative near-cancellation, its numerical
# magnitude can become extreme under the sum-to-one normalization.
# We keep all nine causes in alpha and kappa, move only the female
# Diabetes and kidney diseases beta to a separate companion chart.

cat("Select alpha.csv\n"); A <- read.csv(file.choose(), stringsAsFactors=FALSE)
cat("Select beta.csv\n"); B <- read.csv(file.choose(), stringsAsFactors=FALSE)
cat("Select kappa.csv\n"); K <- read.csv(file.choose(), stringsAsFactors=FALSE)
cat("Select selected_causes.csv\n"); C <- read.csv(file.choose(), stringsAsFactors=FALSE)
stopifnot(all(c("sex_name","cause_name","age_midpoint","alpha") %in% names(A)),
          all(c("sex_name","cause_name","age_midpoint","beta") %in% names(B)),
          all(c("sex_name","cause_name","year","kappa") %in% names(K)),
          "cause_name" %in% names(C))
causes <- C$cause_name
sexes <- c("Male","Female")
stopifnot(length(causes)==9L, nrow(A)==108L, nrow(B)==108L,
          nrow(K)==648L, setequal(unique(A$cause_name),causes),
          setequal(unique(B$cause_name),causes),
          setequal(unique(K$cause_name),causes))
flagged <- "Diabetes and kidney diseases"
stopifnot(flagged %in% causes)
alpha0 <- A; beta0 <- B; kappa0 <- K
diagnostics <- data.frame()
for (sx in sexes) for (cs in causes) {
  ib <- which(B$sex_name==sx & B$cause_name==cs)
  ik <- which(K$sex_name==sx & K$cause_name==cs)
  ia <- which(A$sex_name==sx & A$cause_name==cs)
  stopifnot(length(ib)==6L,length(ik)==36L,length(ia)==6L)
  b <- B$beta[ib]; k <- K$kappa[ik]; a <- A$alpha[ia]
  # center the index; move the old mean contribution into alpha
  km <- mean(k); A$alpha[ia] <- a + b*km
  K$kappa[ik] <- k-km
  sc <- sum(b)
  if (!is.finite(sc) || abs(sc)<1e-9) stop("beta sum too close to zero: ",sx," / ",cs)
  B$beta[ib] <- b/sc; K$kappa[ik] <- K$kappa[ik]*sc
  err <- max(abs(outer(b, k, "*")+a -
                 (outer(B$beta[ib],K$kappa[ik],"*")+A$alpha[ia])))
  if(err>1e-8)stop("Reparameterization changed the fitted log-rates")
  diagnostics <- rbind(diagnostics, data.frame(sex_name=sx,cause_name=cs,
       cancellation_ratio=abs(sc)/sum(abs(b)), max_abs_beta_standard=max(abs(B$beta[ib])),
       fitted_log_rate_error=err))
}
out <- "outputs/04_results/4_1_single_population"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(A,file.path(out,"figure3_alpha_standard.csv"),row.names=FALSE)
write.csv(B,file.path(out,"figure3_beta_standard.csv"),row.names=FALSE)
write.csv(K,file.path(out,"figure3_kappa_standard.csv"),row.names=FALSE)
write.csv(diagnostics,file.path(out,"figure3_constraint_diagnostics.csv"),row.names=FALSE)

pal <- grDevices::hcl.colors(length(causes),"Dark 3")
names(pal) <- causes
plot_six <- function() {
 op <- par(no.readonly=TRUE);on.exit(par(op))
 layout(matrix(c(1,2,3,4,5,6,7,7,7),3,byrow=TRUE),heights=c(1,1,.25))
 par(mar=c(3.6,4.2,3.6,1),oma=c(1,1,3.6,.5))
 for(sx in sexes) for(var in c("alpha","beta","kappa")) {
    dat <- switch(var,alpha=A,beta=B,kappa=K)
    if(sx=="Female" && var=="beta") dat <- subset(dat,cause_name!=flagged)
    dat <- subset(dat,sex_name==sx)
    X <- if(var=="kappa")"year" else "age_midpoint"
    vals <- range(dat[[var]])
    plot(NA,xlim=range(dat[[X]]),ylim=vals,xlab=if(var=="kappa")"Year" else "Age midpoint",
         ylab=var,main=paste(sx,var))
    for(cs in causes) {
      sub <- subset(dat,cause_name==cs); if(!nrow(sub))next
      sub <- sub[order(sub[[X]]),]
      lines(sub[[X]],sub[[var]],col=pal[cs],lwd=1.7)
    }
    if(sx=="Female"&&var=="beta")
      mtext("Female diabetes/kidney beta shown separately",side=3,line=.3,cex=.65)
 }
 par(mar=c(.5,1,.5,1))
 plot.new()
 legend("center",legend=causes,col=pal,lty=1,lwd=2,cex=.68,ncol=3,bty="n")
 mtext("Figure 3 — 9 causes / beta companion / GBD 1980–2015 (not HCSO)",outer=TRUE,side=3,line=1,cex=1.0)
}
png(file.path(out,"figure3_all_causes.png"),width=2000,height=1600,res=170)
plot_six();dev.off()
ib <- B$sex_name=="Female" & B$cause_name==flagged
dat <- B[ib,][order(B$age_midpoint[ib]),]
png(file.path(out,"figure3_female_diabetes_beta_companion.png"),width=1100,height=750,res=160)
plot(dat$age_midpoint,dat$beta,type="o",lwd=2,pch=16,
  xlab="Age midpoint",ylab="beta (sum beta=1)",
  main="Female | Diabetes and kidney diseases | separate beta")
abline(h=0,col="grey70");grid()
dev.off()
print(diagnostics)
cat("Figure 3 main and companion chart saved to ",out,"\n")
