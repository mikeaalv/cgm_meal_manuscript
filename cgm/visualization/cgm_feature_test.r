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
require(dplyr)
require(data.table)
require(ComplexHeatmap)
require(ggpubr)
require(circlize)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
feattab=read.table("cgm_foods_manual_features.csv",header=TRUE,sep=",")
feattab=feattab[,c("subject","foods","rep","mitigator","peak_value","baseline_glucose","AUC_above_baseline","time_to_peak")]
feattab$peak_relative=feattab$peak_value-feattab$baseline_glucose
# foods spikers
food_list_carbs=c("Beans","Berries","Bread","Grapes","Pasta","Potatoes","Rice+Fiber","Rice+Fat","Rice+Protein")
feat_to_compares=c("peak_relative","AUC_above_baseline","time_to_peak")
# treat rice as the baseline/default
stat_coll=c()
for(foods in food_list_carbs){
    for(feature in feat_to_compares){
        xvec=feattab[feattab[,"foods"]=="Rice",feature]
        yvec=feattab[feattab[,"foods"]==foods,feature]
        stat_coll=rbind(stat_coll,data.frame(foods=foods,feature=feature,pval=t.test(xvec,yvec)$p.value))
    }
}
stat_coll$padj=p.adjust(stat_coll$pval,method="fdr")
write.table(stat_coll,file="foodvs_ttest.txt",row.names=FALSE)