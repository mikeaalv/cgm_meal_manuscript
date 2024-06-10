rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(Hmisc)
require(readxl)
require(dplyr)
require(data.table)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
load("heatmap_plot.RData")
# data preparation 
temptab=summtab2[match(subjects,summtab2[,"subject"]),c("subject","ricespiker","potato_vs_grape_val")]
temptab$subject=subjects
stattab=temptab
temptab=mat_mitigator_effect
colnames(temptab)=c("Fat","Fiber","Protein")
stattab=cbind(stattab,temptab[subjects,])
metadatashow[,"subject"]<-metadatashow[,"study_id"]%>%str_remove(string=.,pattern="STUDYID-")%>%as.numeric()
temptab=metadatashow[,c("subject","age_today_avg_all","bmi_avg_all","sspg_status_heyjun","DI_status_heyjun","ethnicity","a1c_avg_all","fbg_avg_all","ffa_avg_heyjun","sspg_avg_all","systolic_avg_all","diastolic_avg_all","ethnicity")]
stattab=cbind(stattab,temptab[match(subjects,temptab[,"subject"]),-1])
temptab=as.data.frame(mat_worstfood[,"Bread",drop=FALSE])
colnames(temptab)="breadspiker"
stattab=cbind(stattab,temptab[subjects,,drop=FALSE])
temptab=mat_quantile[,c("Potatoes","Rice")]
stattab=cbind(stattab,temptab[subjects,])
# 
stattab$"IsAsian"=ifelse(str_detect(stattab[,"ethnicity"],pattern="Asian"),0.5,-0.5)
stattab$"IsBreadspiker"=ifelse(stattab[,"breadspiker"]==1,0.5,-0.5)
stattab$"IsRicespiker"=ifelse(stattab[,"ricespiker"]=="ricespiker",0.5,-0.5)
stattab$"IsMuscleIR"=ifelse(stattab[,"sspg_status_heyjun"]=="IR",0.5,-0.5)
stattab$"IsBetaDys"=ifelse(stattab[,"DI_status_heyjun"]=="BC_dys",0.5,-0.5)
# model list
modelist=list("potato_IR"="Potatoes~IsMuscleIR","potato_beta"="Potatoes~IsBetaDys","Fat_IR"="Fat~IsMuscleIR","Fiber_IR"="Fiber~IsMuscleIR","Protein_IR"="Protein~IsMuscleIR","Fat_Bcell"="Fat~IsBetaDys","Fiber_Bcell"="Fiber~IsBetaDys","Protein_Bcell"="Protein~IsBetaDys","Rice_Asian"="Rice~IsAsian","Fat_Rice"="Fat~IsRicespiker","Fiber_Rice"="Fiber~IsRicespiker","Protein_Rice"="Protein~IsRicespiker","SSPG_PG"="sspg_avg_all~potato_vs_grape_val","FFA_PG"="ffa_avg_heyjun~potato_vs_grape_val","A1C_PG"="a1c_avg_all~potato_vs_grape_val","FBG_PG"="fbg_avg_all~potato_vs_grape_val","Sys_Bread"="systolic_avg_all~IsBreadspiker","Dia_Bread"="diastolic_avg_all~IsBreadspiker")
addterms="+bmi_avg_all+age_today_avg_all"
statcoll=c()
for(modelname in names(modelist)){
    formula_mlm_str=paste0(modelist[[modelname]],addterms)
    formula_mlm=as.formula(formula_mlm_str)
    termsall=attr(terms(lme4:::nobars(formula_mlm)),"term.labels")
    targetterm=termsall[1]
    lmdol=lm(formula_mlm,stattab)
    coeftab=summary(lmdol)$coefficients
    statcoll=rbind(statcoll,data.frame(Y=all.vars(formula_mlm)[1],X=targetterm,coef=coeftab[targetterm,"Estimate"],pval=coeftab[targetterm,"Pr(>|t|)"],formula=formula_mlm_str))
}
save(statcoll,file="comb_model_assocheck.RData")
write.table(statcoll,file="comb_model_assocheck.txt",row.names=FALSE)