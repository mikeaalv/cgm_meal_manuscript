# PCA and association wiht metadata 
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(readxl)
require(dplyr)
require(data.table)
require(reshape2)
require(ggfortify)
require(circlize)
require(ComplexHeatmap)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/omics/")
dirdata=paste0(pardir,"data/metadata/")
setwd(resdir)
# 
load(paste0(resdir,"metabolomics_lipidomics_processed.RData"))
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
grouplist=c("sex.factor","ethnicity","a1c_t2d_status_bl","sspg_status_heyjun","di_3classes_heyjun","ie_3_classes_heyjun","FFA_3classes_heyjun","hepatic_ir_3classes_heyjun")
numlist=c("age_today_avg_all","bmi_avg_all","a1c_avg_all","fbg_avg_all","sspg_heyjun","modified_DI_heyjun","ie_heyjun","ffa_avg_heyjun","hepatic_IR_heyjun")
metadata=metadata[,c("study_id",grouplist,numlist)]
metadata$transformedid<-metadata$study_id%>%str_replace_all(string=.,pattern="STUDYID-",replacement="")%>%as.numeric()
for(groupfeat in grouplist){
    metadata[which(metadata[,groupfeat]=="Unknown"),groupfeat]=NA
}
omics=names(omicslist_arra)
stat_tab_list=list()
list_pca=list()
for(omic in omics){
    mat=as.data.frame(t(omicslist_arra[[omic]]))
    pca_res=prcomp(mat,center=FALSE)
    list_pca[[omic]]=pca_res
    scores=pca_res$x
    # metadata reformat
    samples=rownames(mat)
    subjectids<-samples%>%str_extract(string=.,pattern="^X\\d+")%>%str_replace_all(string=.,pattern="X",replacement="")
    meta_tranform_df=as.data.frame(metadata[match(subjectids,metadata[,"transformedid"]),])
    # PCA plot with metabolic subtypes and sex, 
    for(grp in grouplist){
        locdf=as.data.frame(cbind(scores[,c(1,2)],meta_tranform_df[,grp]))
        colnames(locdf)=c("PC1","PC2",grp)
        locdf[,grp]=as.factor(locdf[,grp])
        # 
        plotdf=mat
        plotdf$class=meta_tranform_df[,grp]
        p<-autoplot(pca_res,data=plotdf,colour="class")
        ggsave(paste0(omic,"_",grp,".pdf"),plot=p)
        for(pci in c(1,2)){
            locdf_in=locdf[,c(pci,3)]
            colnames(locdf_in)=c("y","type")
            lmodel=lm(y~type,data=locdf_in)
            corr=sqrt(summary(lmodel)$r.squared)
            pval=summary(lmodel)$coefficients[2,"Pr(>|t|)"]
            stat_tab_list[[length(stat_tab_list)+1]]=data.frame(omics=omic,feature=grp,pc=pci,test="lm_level",value=corr,pvalue=pval)
        }     
    }
    # PCA plot with numerical values
    for(grp in numlist){
        locdf=as.data.frame(cbind(scores[,c(1,2)],meta_tranform_df[,grp]))
        colnames(locdf)=c("PC1","PC2",grp)
        # 
        plotdf=mat
        plotdf$val=meta_tranform_df[,grp]
        p<-autoplot(pca_res,data=plotdf,colour="val")
        ggsave(paste0(omic,"_",grp,"_num.pdf"),plot=p) 
    }
    # correlation PC1-2 with metabolic testing result, age, bmi
    for(numfeat in numlist){
        locdf=as.data.frame(cbind(scores[,c(1,2)],meta_tranform_df[,numfeat]))
        for(pci in c(1,2)){
            corres=cor.test(x=locdf[,pci],y=locdf[,3],method="spearman")
            corr=corres$estimate
            pval=corres$p.value
            stat_tab_list[[length(stat_tab_list)+1]]=data.frame(omics=omic,feature=numfeat,pc=pci,test="correlation",value=corr,pvalue=pval)
        }  
    }
    if(omic!="proteomics"){
        savdf=cbind(scores[,c(1,2)],meta_tranform_df[,c("sex.factor","hepatic_ir_3classes_heyjun")])
        colnames(savdf)[3:4]=c("sex","hepatic_ir_classes")
        write.table(savdf,file=paste0("figext6",omic,".txt"),row.names=FALSE)
    }
}
stat_tab=Reduce(rbind,stat_tab_list)
stat_tab$padj=p.adjust(stat_tab[,"pvalue"],method="fdr")
save(stat_tab,file=paste0(resdir,"assoc_stat.RData"))
# heatmap of correlations with PC 1 and 2 in metabolomics and lipidomics
tab_cor=stat_tab[stat_tab[,"test"]=="correlation"&stat_tab[,"omics"]%in%c("metabolomics","lipidomics"),]
clinicfeats=unique(tab_cor[,"feature"])
cormat=matrix(NA,nrow=length(clinicfeats),ncol=4)
rownames(cormat)=clinicfeats
colnames(cormat)=c("metabolomicsPC1","metabolomicsPC2","lipidomicsPC1","lipidomicsPC2")
pmat=cormat
for(rowi in seq(nrow(tab_cor))){
    record=tab_cor[rowi,]
    colid=paste0(record[,"omics"],"PC",record[,"pc"])
    pmat[record[,"feature"],colid]=record[,"pvalue"]
    cormat[record[,"feature"],colid]=record[,"value"]
}
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("cor_pc_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# loading plot
featlist=list("metabolomics"=annolist[["metabolomics"]][,"Compound.name"],"lipidomics"=annolist[["lipidomics"]][,"Compound.name"],"proteomics"=annolist[["proteomics"]][,"UniProt"])
for(omic in names(list_pca)){
    pca_res=list_pca[[omic]]
    featname=featlist[[omic]]
    loadingmat=pca_res$rotation
    rownames(loadingmat)=featname
    ind=order(abs(loadingmat[,1]),decreasing=TRUE)
    barplot(loadingmat[ind[1:50],1])
    # fviz_contrib(pca_res, choice="var", axes=1,top=10)
}
