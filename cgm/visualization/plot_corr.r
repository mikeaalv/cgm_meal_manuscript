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
require(circlize)
# modify to you local path for it to work
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/cgm_meal/")
setwd(resdir)
# 
load("heatmap_plot.RData")
col_fun=colorRamp2(c(-0.5,0,0.5),c("blue","white","red"))
# whole correlation matirx between features and spikes
foodmat=merge(mat_quantile,mat_mitigator_effect,by="row.names")
fullmat=merge(metadatashow_plot,foodmat,by.x="row.names",by.y="Row.names")
numcollist=lapply(metadatashow_plot,function(x) {class(x)=="numeric"|class(x)=="integer"})
numcol=names(numcollist)[unlist(numcollist)]
foodcol=colnames(foodmat)[-1]
p_mat=matrix(NA,nrow=length(numcol),ncol=length(foodcol))
colnames(p_mat)=foodcol
rownames(p_mat)=numcol
for(foodele in foodcol){
    for(numele in numcol){
        xvec=fullmat[,foodele]
        yvec=fullmat[,numele]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        p_mat[numele,foodele]=corres$p.value
    }
}
# correlation heatmap between mitigator effect and metadata  
featlist=c("modified_DI_heyjun","systolic_avg_all","diastolic_avg_all")
n_mitigat=dim(mat_mitigator_effect)[2]
cor_mat=matrix(NA,nrow=length(featlist),ncol=n_mitigat)
colnames(cor_mat)=colnames(mat_mitigator_effect)
rownames(cor_mat)=featlist
p_mat=cor_mat
for(feat in featlist){
    for(mitigi in seq(n_mitigat)){
        xvec=mat_mitigator_effect[,mitigi]
        yvec=metadatashow_plot[,feat]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[feat,mitigi]=corres$estimate
        p_mat[feat,mitigi]=corres$p.value
    }
}
pdf(paste0("heatmap_mitigatoreffect_metadata.pdf"))
print(Heatmap(cor_mat,name="correlation with metadata",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
stattab=c()
p_mat_adj=p.adjust(c(p_mat),method="fdr")
dim(p_mat_adj)=dim(p_mat)
rownames(p_mat_adj)=rownames(p_mat)
colnames(p_mat_adj)=colnames(p_mat)
for(feat in rownames(p_mat_adj)){
    for(foods in colnames(p_mat_adj)){
        stattab=rbind(stattab,data.frame(foods=foods,features=feat,correlation=cor_mat[feat,foods],pvalue=p_mat[feat,foods],p.adjust=p_mat_adj[feat,foods]))
    }
}
# correlation heatmap between carb quantile and metadata  
featlist=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","insulin_fasting_avg_all","C_peptide_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean","HOMA_IR_heyjun","HOMA_S_heyjun","ldl_avg_all","ldl_hdl_ratio_avg_all","triglyceride_avg_all","cholesterol_total_avg_all","systolic_avg_all","ffa_avg_heyjun","hepatic_IR_heyjun","modified_DI_heyjun")
n_mitigat=dim(mat_quantile)[2]
cor_mat=matrix(NA,nrow=length(featlist),ncol=n_mitigat)
colnames(cor_mat)=colnames(mat_quantile)
rownames(cor_mat)=featlist
p_mat=cor_mat
for(feat in featlist){
    for(mitigi in seq(n_mitigat)){
        xvec=mat_quantile[,mitigi]
        yvec=metadatashow_plot[,feat]
        nonaind=(!is.na(xvec))&(!is.na(yvec))
        if(length(unique(yvec[nonaind]))<=1){
            next
        }
        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
        cor_mat[feat,mitigi]=corres$estimate
        p_mat[feat,mitigi]=corres$p.value
    }
}
pdf(paste0("heatmap_spikequantile_metadata.pdf"))
print(Heatmap(cor_mat,name="correlation with metadata",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
p_mat_adj=p.adjust(c(p_mat),method="fdr")
dim(p_mat_adj)=dim(p_mat)
rownames(p_mat_adj)=rownames(p_mat)
colnames(p_mat_adj)=colnames(p_mat)
for(feat in rownames(p_mat_adj)){
    for(foods in colnames(p_mat_adj)){
        stattab=rbind(stattab,data.frame(foods=foods,features=feat,correlation=cor_mat[feat,foods],pvalue=p_mat[feat,foods],p.adjust=p_mat_adj[feat,foods]))
    }
}
nameexch=c("DI","systolic bp","diastolic bp","A1C","FBG","OGTT @120 mins","SSPG","Fasting insulin","C peptide","Glucagon","HOMA IR","HOMA S","LDL","LDL/HDL","triglyceride","Total cholesterol","FFA","Hepatic IR")
names(nameexch)=c("modified_DI_heyjun","systolic_avg_all","diastolic_avg_all","a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","insulin_fasting_avg_all","C_peptide_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean","HOMA_IR_heyjun","HOMA_S_heyjun","ldl_avg_all","ldl_hdl_ratio_avg_all","triglyceride_avg_all","cholesterol_total_avg_all","ffa_avg_heyjun","hepatic_IR_heyjun")
for(repl in names(nameexch)){
    stattab[stattab[,"features"]==repl,"features"]=nameexch[repl]
}
write.table(stattab,file="spike_meta_corr_tab_peak.txt",row.names=FALSE)
# (summarized) stack spike correlated with metadata
glucfeat=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","insulin_fasting_avg_all","C_peptide_metabolic_testing_mean","Glucagon_pmol_L_metabolic_testing_mean","HOMA_IR_heyjun","HOMA_S_heyjun","systolic_avg_all","hepatic_IR_heyjun","modified_DI_heyjun")
lipidsfeat=c("ldl_avg_all","ldl_hdl_ratio_avg_all","triglyceride_avg_all","cholesterol_total_avg_all")
featlist=c(glucfeat,lipidsfeat)
n_mitigat=1
cor_mat=matrix(NA,nrow=length(featlist),ncol=n_mitigat)
colnames(cor_mat)="summarized"
rownames(cor_mat)=featlist
p_mat=cor_mat
# 
mat_quantile_vec=c(mat_quantile)
for(feat in featlist){
    xvec=mat_quantile_vec
    yvec=rep(metadatashow_plot[,feat],times=dim(mat_quantile)[2])
    nonaind=(!is.na(xvec))&(!is.na(yvec))
    if(length(unique(yvec[nonaind]))<=1){
        next
    }
    corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
    cor_mat[feat,]=corres$estimate
    p_mat[feat,]=corres$p.value
}
pdf(paste0("heatmap_spikequantile_metadata_summar_glucosemetab.pdf"))
cor_mat_tmp=cor_mat[glucfeat,,drop=FALSE]
p_mat_tmp=p_mat[glucfeat,,drop=FALSE]
print(Heatmap(cor_mat_tmp,name="correlation with metadata",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat_tmp[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
pdf(paste0("heatmap_spikequantile_metadata_summar_lipidmetab.pdf"))
cor_mat_tmp=cor_mat[lipidsfeat,,drop=FALSE]
p_mat_tmp=p_mat[lipidsfeat,,drop=FALSE]
print(Heatmap(cor_mat_tmp,name="correlation with metadata",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat_tmp[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# (summarized) correlation heatmap between stacked mitigator effect and metadata  
featlist=c("modified_DI_heyjun","systolic_avg_all","diastolic_avg_all")
n_mitigat=1
cor_mat=matrix(NA,nrow=length(featlist),ncol=n_mitigat)
colnames(cor_mat)="summarized"
rownames(cor_mat)=featlist
p_mat=cor_mat
mat_mitigator_effect_vec=c(mat_mitigator_effect)
for(feat in featlist){
    xvec=mat_mitigator_effect_vec
    yvec=rep(metadatashow_plot[,feat],times=dim(mat_mitigator_effect)[2])
    nonaind=(!is.na(xvec))&(!is.na(yvec))
    if(length(unique(yvec[nonaind]))<=1){
        next
    }
    corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
    cor_mat[feat,]=corres$estimate
    p_mat[feat,]=corres$p.value
}
pdf(paste0("heatmap_mitigatoreffect_metadata_summar.pdf"))
print(Heatmap(cor_mat,name="correlation with metadata",cluster_columns=FALSE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=FALSE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()