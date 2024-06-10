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
require(ComplexHeatmap)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/pathway/")
setwd(resdir)
# 
load("pathway_heatmap_plot.RData")
# select repeated omics features (>=3 in carb spike or >=2 in mitigation effect)
mitiga_mask=str_detect(string=names(list_sig),pattern=fixed("+"))
for(dirc in c("posicorr","negacorr")){
    mask1=str_detect(string=names(list_sig),pattern=fixed(dirc))
    for(type in c("carb","mitigation")){
        if(type=="mitigation"){
            mask2=mitiga_mask
            thres=2
        }else{
            mask2=!mitiga_mask
            thres=3
        }
        loclist=table(unlist(list_sig[mask1&mask2]))
        loclist=loclist[loclist>=thres]
        featlist=c(featlist,names(loclist))
    }
}
featlist=unique(featlist)
for(omic in names(annolist)){
    featplot=featlist[str_detect(string=featlist,pattern=fixed(omic))]
    featplot=str_remove(string=featplot,pattern=fixed(omic))
    omicsind=sapply(featplot,function(x) which(featvec==x)[1])
    locmat=omicsmat[,omicsind]
    colnames(locmat)=featplot
    showmat=cbind(food_sig_mat,locmat)
    colnames(showmat)=c(foodtypes,featplot)
    # 
    cor_mat=matrix(NA,nrow=length(featplot),ncol=length(foodtypes))
    colnames(cor_mat)=foodtypes
    rownames(cor_mat)=featplot
    p_mat=matrix(1,nrow=length(featplot),ncol=length(foodtypes))
    colnames(p_mat)=foodtypes
    rownames(p_mat)=featplot
    for(feat in featplot){
        for(food in foodtypes){
            xvec=showmat[,feat]
            yvec=showmat[,food]
            nonaind=(!is.na(xvec))&(!is.na(yvec))
            if(length(unique(yvec[nonaind]))<=1){
                next
            }
            corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
            cor_mat[feat,food]=corres$estimate
            p_mat[feat,food]=corres$p.value
        }
    }
    # 
    h=Heatmap(cor_mat,name="corrmat",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,
            cell_fun=function(j,i,x,y,w,h,fill){
                if(p_mat[i,j]<0.05){
                    grid.text("*",x,y)
                }
            },
            col=circlize::colorRamp2(c(-1,0,1),c("blue","white","red")))#
    pdf(paste0("omics_corr_heatmap_highsig.",omic,".pdf"))
    draw(h,main_heatmap="corrmat")
    dev.off()
}