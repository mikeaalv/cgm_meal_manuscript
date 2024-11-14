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
require(lsr)
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
        teststat=t.test(xvec,yvec)
        stat_coll=rbind(stat_coll,data.frame(foods=foods,feature=feature,pval=teststat$p.value,statistics=teststat$statistic,ci=paste0(formatC(teststat$conf.int,format="E",digits=2),collapse=" "),Df=teststat$parameter,effectsz=cohensD(xvec,yvec),delta=mean(yvec)-mean(xvec)))
    }
}
stat_coll_padj=c()
for(foodthe in food_list_carbs){
    subtab=stat_coll[stat_coll[,"foods"]==foodthe,]
    subtab$padj=p.adjust(subtab$pval,method="fdr")
    stat_coll_padj=rbind(stat_coll_padj,subtab)
}
write.table(stat_coll_padj,file="foodvs_ttest.txt",row.names=FALSE)
# paired test for mitigators
stat_coll_pair=c()
mitig_comb=c("Rice+Fiber","Rice+Fat","Rice+Protein")
for(foodthe in mitig_comb){
    for(feature in feat_to_compares){
        xtab <- feattab %>% filter(foods %in% c("Rice")) %>% group_by(subject) %>% summarize(meanval=mean(.data[[feature]]),.groups="keep") %>% group_by(subject) %>% as.data.frame()
        ytab <- feattab %>% filter(foods %in% c(foodthe)) %>% group_by(subject) %>% summarize(meanval=mean(.data[[feature]]),.groups="keep") %>% group_by(subject) %>% as.data.frame()
        subjs=ytab[,"subject"]
        xorder=match(subjs,xtab[,"subject"])
        teststat=t.test(xtab[xorder,"meanval"],ytab[,"meanval"],paired=TRUE)
        stat_coll_pair=rbind(stat_coll_pair,data.frame(foods=foodthe,feature=feature,pval=teststat$p.value,statistics=teststat$statistic,ci=paste0(formatC(teststat$conf.int,format="E",digits=2),collapse=" "),Df=teststat$parameter,effectsz=cohensD(xtab[,"meanval"],ytab[,"meanval"]),delta=mean(ytab[,"meanval"])-mean(xtab[xorder,"meanval"])))
    }
}
stat_coll_pair_padj=c()
for(foodthe in mitig_comb){
    subtab=stat_coll_pair[stat_coll_pair[,"foods"]==foodthe,]
    subtab$padj=p.adjust(subtab$pval,method="fdr")
    stat_coll_pair_padj=rbind(stat_coll_pair_padj,subtab)
}
write.table(stat_coll_pair_padj,file="foodvs_ttest.mitig_pair.txt",row.names=FALSE)
