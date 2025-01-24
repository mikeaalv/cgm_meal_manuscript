rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(Hmisc)
require(npreg)
require(readxl)
require(plyr)
require(dplyr)
require(data.table)
require(ComplexHeatmap)
require(ggpubr)
require(circlize)
require(performance)
require(lme4)
require(lmerTest)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
feattab=read.table("cgm_foods_manual_features.csv",header=TRUE,sep=",")
feattab=feattab[,c("subject","foods","rep","mitigator","peak_value","baseline_glucose","AUC_above_baseline","time_to_peak","return_to_baseline_time")]
feattab$peak_relative=feattab$peak_value-feattab$baseline_glucose
foodslist=c("Beans","Berries","Bread","Grapes","Pasta","Potatoes","Rice")
feattab=feattab[feattab[,"foods"]%in%foodslist,]
# 
statchelist=c("peak_relative","AUC_above_baseline","time_to_peak","return_to_baseline_time")
# listres=list()
iccdf=c()
for(food in foodslist){
    for(thestat in statchelist){
        set.seed(1)
        subtab=feattab[feattab[,"foods"]==food,c("subject","foods","rep",thestat)]
        colnames(subtab)[4]="feature"
        # model=lmer(feature ~ foods + (1|subject),data=subtab)
        model=lmer(feature ~ (1|subject),data=subtab)
        # listres[[thestat]]=summary(model)
        valvec=as.data.frame(icc(model,ci=TRUE))[,"ICC_unadjusted"]
        dftmp=data.frame(ICC=valvec[1],CI=paste0(formatC(valvec[-1],format="f",digits=2),collapse=" "),feature=thestat,foods=food)
        iccdf=rbind(iccdf,dftmp)
    }
}
save(iccdf,file="icc_tab.RData")
write.table(iccdf,file="icc_tab.txt",row.names=FALSE)