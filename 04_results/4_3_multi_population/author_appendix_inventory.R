# Author appendix 9 – enumerate the exact 52 original PNGs as a reference catalogue.
# This is an index of source images, NOT a claim that HCSO data are available.
# Original: https://github.com/LiviaVarga/Forecasting-cause-of-death-mortality-with-single-and-multi-population-models-in-Hungary/tree/main/9.%20The%20fitted%20parameters%20of%20the%20multi-population%20models

author_appendix_inventory <- function(path) {
  rr <- "https://github.com/LiviaVarga/Forecasting-cause-of-death-mortality-with-single-and-multi-population-models-in-Hungary/blob/main/"
  base <- "9. The fitted parameters of the multi-population models/"
  male_carter <- c("I_XIV","II_XI","IX_XX")
  female_carter <- "IX_XI_XIV_XX"
  male_lc2 <- c("I_XIV","II_IX_XI_XX")
  female_lc2 <- "IV_IX_XI_XIV_XX"
  spec <- list(
    list(model="Poisson LC1",single="Poisson Lee–Carter",joint="Poisson Carter–Lee",n=1,
         M=male_carter,F=female_carter),
    list(model="Binomial LC1",single="Binomial Lee–Carter",joint="Binomial Carter–Lee",n=1,
         M=male_carter,F=female_carter),
    list(model="Poisson LC2",single="Poisson Lee–Carter",joint="Poisson Lee–Carter",n=2,
         M=male_lc2,F=female_lc2),
    list(model="Binomial LC2",single="Binomial Lee–Carter",joint="Binomial Lee–Carter",n=2,
         M=male_lc2,F=female_lc2),
    list(model="CBD",single="CBD",joint="CBD",n=2,
         M=male_carter,F=female_carter),
    list(model="Plat2",single="reduced Plat",joint="reduced Plat",n=2,
         M=male_carter,F=female_carter),
    list(model="Plat3",single="reduced Plat",joint="reduced Plat",n=3,
         M=male_carter,F=female_carter)
  )
  rows<-list()
  for(s in spec) for(sex in c("males","females")) {
    groups<-if(sex=="males")s$M else s$F
    for(codes in groups) {
      if(grepl("LC1$",s$model)) {
        par<-paste0("The fitted parameters of the ",s$joint,
                    " model for ",sex,"_",codes,".png")
        comp<-paste0("The common and cause-specific mortality indices of the ",
                     s$single," and Carter–Lee models for ",sex,"_",codes,".png")
      } else if(grepl("LC2$",s$model)) {
        par<-paste0("The fitted parameters of the ",s$joint,
                    " model with two factors for ",sex,"_",codes,"_multi.png")
        comp<-paste0("The common and cause-specific mortality indices of the ",
                      s$single," models with two factors for ",sex,"_",codes,".png")
      } else if(s$model=="CBD") {
        par<-paste0("The fitted parameters of the CBD model for ",sex,"_",codes,"_multi.png")
        comp<-paste0("The common and cause-specific mortality indices of the CBD models for ",
                     sex,"_",codes,".png")
      } else {
        k<-if(s$n==2L)"two" else "three"
        par<-paste0("The fitted parameters of the reduced Plat model with ",k,
                    " factors for ",sex,"_",codes,"_multi.png")
        comp<-paste0("The common and cause-specific mortality indices of the reduced Plat models with ",
                      k," factors for ",sex,"_",codes,".png")
      }
      for(typ in c("joint_parameters","indices_comparison")) {
        f<-if(typ=="joint_parameters")par else comp
        rows[[length(rows)+1L]]<-data.frame(
          sex=if(sex=="males")"Male" else "Female",
          model=s$model,ICD10_chapters=codes,image_type=typ,
          original_basename=f,
          source_url=paste0(rr,utils::URLencode(paste0(base,f),reserved=TRUE)),
          interpretation=if(typ=="joint_parameters")
            "Age/period joint fitted parameters; compare shared and cause-specific components"
            else "A common jointly estimated period index compared with independent cause-specific period indices; different model units need identification checks",
          stringsAsFactors=FALSE)
      }
    }
  }
  ans<-do.call(rbind,rows)
  stopifnot(nrow(ans)==52L,
    sum(ans$image_type=="joint_parameters")==26L,
    sum(ans$image_type=="indices_comparison")==26L)
  write.csv(ans,path,row.names=FALSE,fileEncoding="UTF-8")
  invisible(ans)
}
