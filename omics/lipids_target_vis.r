rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(igraph)
require(Hmisc)
require(npreg)
require(ppcor)
require(readxl)
require(lipidr)
# 
require(clusterProfiler)
require(org.Hs.eg.db)
require(ReactomePA)
require(metpath)
require(ComplexHeatmap)
require(pathview)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/pathway/")
setwd(resdir)
load("pathway_heatmap_plot.RData")
tabmatch=read.table(paste0(pardir,"data/metabolomics_lipidomics/Lipomatdb.txt"),header=TRUE,sep="\t")
tabmatch[,"Lipid_Name"]=str_replace_all(string=tabmatch[,"Lipid_Name"],pattern="\\_",replacement="/")
# mean lipids length associated with each food 
food_dirc=names(list_sig)
lipid_interest_list=c("TAG","PE","PC","FFA","DAG")
plotdf=c()
for(fd_name in food_dirc){
    parts=str_split(string=fd_name,pattern="corr")[[1]]
    featlist=list_sig[[fd_name]]
    lipidlist=featlist[str_which(string=featlist,pattern="lipidomics")]
    lipidlist=str_remove(string=lipidlist,pattern="lipidomics")
    summtab=tabmatch[match(lipidlist,tabmatch[,"Lipid_Name"]),]
    summtab$raw_name=lipidlist
    for(lipidele in lipid_interest_list){
        ind=str_which(string=summtab[,"Lipid_Name"],pattern=paste0("^",lipidele))
        if(length(ind)==0){
            len=0
        }else{
            len=mean(summtab[ind,"Total_Carb"])
        }
        even_mask=(summtab[ind,"Total_Carb"]%%2==0)
        plotdf=rbind(plotdf,data.frame(food=parts[2],direction=parts[1],class=lipidele,len=len,noddf=sum(!even_mask),neven=sum(even_mask)))
    }
}
p<-ggplot(plotdf,aes(x=food,y=len,fill=class))+geom_bar(stat="identity",position="dodge2")+facet_wrap(~direction,ncol=1)+theme(strip.background=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))#
ggsave("lipids_len.pdf",plot=p)
p<-ggplot(plotdf,aes(x=food,y=noddf,fill=class))+geom_bar(stat="identity",position="dodge2")+facet_wrap(~direction,ncol=1)+theme(strip.background=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))#
ggsave("lipids_odd.pdf",plot=p)
p<-ggplot(plotdf,aes(x=food,y=neven,fill=class))+geom_bar(stat="identity",position="dodge2")+facet_wrap(~direction,ncol=1)+theme(strip.background=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))#
ggsave("lipids_even.pdf",plot=p)
# general numbers statistics
lipidtab=read.table(paste0(pardir,"data/metabolomics_lipidomics/CGM_Plasma_Lipid_Species(1).csv"),header=TRUE)
lipids=lipidtab[,1]
print(table(tabmatch[match(lipids,tabmatch[,"Lipid_Name"]),"Lipid_Class"]))
lipid_sig=unlist(list_sig)
lipid_sig=lipid_sig[str_which(string=lipid_sig,pattern="lipidomics")]
lipid_sig=unique(str_remove(string=lipid_sig,pattern="lipidomics"))
print(table(tabmatch[match(lipid_sig,tabmatch[,"Lipid_Name"]),"Lipid_Class"]))
# test on TAG unsaturated level
statdf=c()
for(fd_name in food_dirc){
    parts=str_split(string=fd_name,pattern="corr")[[1]]
    featlist=list_sig[[fd_name]]
    lipidlist=featlist[str_which(string=featlist,pattern="lipidomics")]
    c_n_tab=lipidlist%>%str_extract(string=.,pattern="TAG\\d+:\\d+")%>%str_remove(string=.,pattern="TAG")%>%str_split(string=.,pattern=":",simplify=TRUE)
    if(dim(c_n_tab)[1]<1){
        next
    }
    c_n_tab=c_n_tab[!is.na(c_n_tab[,1]),,drop=FALSE]
    class(c_n_tab)="numeric"
    if(dim(c_n_tab)[1]<1){
        next
    }
    unsat_ratio_vec=c_n_tab[,2]/c_n_tab[,1]
    nlip=length(unsat_ratio_vec)
    statdf=rbind(statdf,data.frame(food=rep(parts[2],times=nlip),direction=rep(parts[1],times=nlip),prop=unsat_ratio_vec))
}
# reformat
statdf2=c()
grouptab=list("simpsug"=c("Berries","Grapes"),"starch"=c("Beans","Bread","Pasta","Potatoes","Rice"))
for(dirc in c("posi","nega")){
    for(groups in names(grouptab)){
        inds=which(statdf[,"food"]%in%grouptab[[groups]] & statdf[,"direction"]==dirc)
        statdf2=rbind(statdf2,data.frame(comparpart=rep(paste0(groups,"_",dirc),times=length(inds)),prop=statdf[inds,"prop"]))
    }
}
comparisons=list("1"=c("simpsug_posi","starch_posi"),"2"=c("simpsug_posi","simpsug_nega"),"3"=c("starch_posi","starch_nega"))
p<-ggboxplot(statdf2,x="comparpart",y="prop")+stat_compare_means(comparisons=comparisons,label="p.signif",hide.ns=TRUE)
ggsave(plot=p,paste0("lipids_group_wise_compare.pdf"))