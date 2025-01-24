rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(ggpubr)
require(tidyr)
require(zoo)
# 
projdir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/"
datadir=paste0(projdir,"result/ressistant_starch/result/Results of Resistant starch K-DTRS kit/")
resdir=paste0(projdir,"result/ressistant_starch/result/")
setwd(resdir)
set.seed(101)
# 
timep=c("before","after")
comparisons=list("1"=c("Pasta","Bread"),"2"=c("Pasta","Rice"),"3"=c("Potatoe","Bread"),"4"=c("Potatoe","Rice"))
for(tp in timep){
    tab=read.table(paste0(datadir,"resistant_starch_",tp,".txt"),sep="\t",header=FALSE)
    tab=tab[,c("V1","V2","V9")]
    colnames(tab)=c("food_rep","starchtype","amount")
    tab[tab[,"food_rep"]=="","food_rep"]=NA
    tab[,"food_rep"]=na.locf(tab[,"food_rep"])#fill by last value
    tab[,"food_rep"]<-tab[,"food_rep"] %>% str_replace_all(string=.,pattern="Control Starch",replacement="Control_Starch") %>% str_trim(string=.,side="both")
    sidtab=str_split(string=tab[,"food_rep"],pattern=fixed(" "),simplify=TRUE)
    colnames(sidtab)=c("food","rep")
    tab=cbind(sidtab,tab)
    tab=tab[tab[,"food"]!="Control_Starch",]
    tab[,"food"]=factor(tab[,"food"],levels=c("Pasta","Potatoe","Bread","Rice"))
    # all features
    p<-ggplot(tab,aes(x=food,y=amount,fill=starchtype))+geom_boxplot()+xlab("foods")+ylab("amount (100g food)")
    ggsave(paste0("starch_dist_",tp,".pdf"),plot=p)
    mask2=tab[,"starchtype"]=="Total Starch"
    for(fr in unique(tab[,"food_rep"])){
        mask1=tab[,"food_rep"]==fr
        tab[mask1,"amount"]=tab[mask1,"amount"]/tab[mask1&mask2,"amount"]
    }
    p<-ggboxplot(tab[tab[,"starchtype"]=="RS",],x="food",y="amount",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE,method="t.test")+xlab("foods")+ylab("RS ratio")#,p.adjust.method="none"
    ggsave(paste0("RS_ratio_",tp,".pdf"),plot=p)
    p<-ggboxplot(tab[tab[,"starchtype"]=="SDS",],x="food",y="amount",add="jitter")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE,method="t.test")+xlab("foods")+ylab("SDS ratio")#,p.adjust.method="none"
    ggsave(paste0("SDS_ratio_",tp,".pdf"),plot=p)
    if(tp=="after"){
        savdf=tab[tab[,"starchtype"]=="RS"|tab[,"starchtype"]=="SDS",c("food","rep","starchtype","amount")]
        colnames(savdf)[dim(savdf)[2]]="ratio"
        write.table(savdf,file=paste0("figext5.txt"))
    }
}
