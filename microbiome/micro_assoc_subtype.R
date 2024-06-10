# partial correlation between t2d subtype features and microbiome
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
require(circlize)
require(ComplexHeatmap)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
resdir=paste0(pardir,"result/microbiome/metabolic_subtypes/")
datadir=paste0(pardir,"result/microbiome/eda/")
setwd(resdir)
# omics data 
load(paste0(datadir,"microbiome_processed.RData"))
subjects=rownames(micromat)
nsubjects=length(subjects)
# add alpha diveristy data
div_tab=read.table(paste0(pardir,"result/microbiome/diversity/diversity_table.tsv"),header=TRUE)
colnames(div_tab)[1]="sample"
div_tab$sample=str_remove(string=div_tab$sample,pattern="\\_trimmed\\.fq\\.gz\\_metaphlan\\_bugs\\_list$")
sampfilematchmeta$sample=str_remove(string=sampfilematchmeta$files,pattern="\\_trimmed\\.fq\\.gz$")
alphafeat=c("chao1","gini_index","observed_otus","shannon","simpson")
div_tab=div_tab[,c("sample",alphafeat)]
addmat=matrix(NA,nrow=length(subjects),ncol=length(alphafeat))
rownames(addmat)=subjects
colnames(addmat)=alphafeat
for(subj in subjects){
    filelist=sampfilematchmeta[sampfilematchmeta[,"subject_id"]==subj,"sample"]
    rowind=div_tab[,"sample"]%in%filelist
    tempmat=apply(div_tab[rowind,alphafeat,drop=FALSE],2,mean)
    addmat[subj,]=tempmat
}
summfeat_list[["alpha_div"]]=alphafeat
micromat=cbind(micromat,addmat)
# metadata
metadata=as.data.frame(read.table(file=paste0(pardir,"data/metadata/metadata_clean_all_ver2.tsv"),header=TRUE))
class_num_cols=c("a1c_avg_all","fbg_avg_all","ogtt_t_120_avg_all","sspg_avg_all","modified_DI_heyjun","ie_heyjun","hepatic_IR_heyjun","HOMA_IR_heyjun","ffa_avg_heyjun","HOMA_S_heyjun","HOMA_B_heyjun","albumin_avg_all","creatinine_avg_all","hscrp_avg_all","ldl_avg_all","hdl_avg_all","triglyceride_avg_all","fructosamine_avg_all","insulin_fasting_avg_all","alt_sgpt_avg_all","systolic_avg_all","diastolic_avg_all","bmi_avg_all")
subphenogrp=class_num_cols[1:9]
obslist=c("go","kegg","pathway","pfam","rxn","species_abs_g","alpha_div")
metadatashow=as.data.frame(matrix(NA,nrow=nsubjects,ncol=length(class_num_cols)+1))
rownames(metadatashow)=subjects
colnames(metadatashow)=c("study_id",class_num_cols)
for(subj in subjects){
    match_id=paste0("STUDYID-",str_pad(subj,3,pad="0"))
    metadatashow[subj,]=metadata[metadata[,"study_id"]==match_id,c("study_id",class_num_cols)]
}
metadatashow=metadatashow[,-1]
# correlation result table 
nclass=length(class_num_cols)
tabstat_cor=c()
wvec=weighttab[match(rownames(metadatashow),weighttab[,"subject_id"]),"meanw"]
for(classele in class_num_cols){
    for(typ in obslist){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        features=summfeat_list[[typ]]
        # select features
        if(!(typ%in%c("species_abs_g","alpha_div"))){
            misvec=missing_ratio_list[[typ]]
            pertvec=avg_perc_list[[typ]]
            features=intersect(names(misvec)[misvec<0.5],features)
            # features=intersect(names(pertvec)[pertvec>0.5],features)
        }
        corvec=rep(NA,times=length(features))
        names(corvec)=features
        pval=corvec
        micromatloc=micromat[nonaind,features]
        ind_dochange=apply(micromatloc,2,var)!=0
        for(feat in features[ind_dochange]){
            res=pcor.test(xvec[nonaind],micromatloc[,feat],wvec[nonaind],"spearman")
            pval[feat]=res$p.value
            corvec[feat]=res$estimate
        }
        p.adj=p.adjust(pval,"fdr")
        seleind=which(c(pval<0.05))
        seleind=seleind[!is.na(seleind)]
        len=length(seleind)
        tabstat_cor=rbind(tabstat_cor,data.frame(type=rep(typ,times=len),clinic=rep(classele,times=len),feature=features[seleind],correlation=corvec[seleind],pvalue=pval[seleind],padj=p.adj[seleind]))
    }
}
save(tabstat_cor,file="association_subtype_microbiome.RData")
# heatmap of all species correlation (>1 match)
tabstat_cor_sub=tabstat_cor[tabstat_cor[,"clinic"]%in%subphenogrp,]
counttab=table(tabstat_cor_sub[tabstat_cor_sub[,"type"]=="species_abs_g","feature"])
feature_spec=names(counttab)[counttab>1]
cormat=matrix(NA,nrow=nclass,ncol=length(feature_spec))
rownames(cormat)=class_num_cols
colnames(cormat)=feature_spec
pmat=cormat
wvec=weighttab[match(rownames(metadatashow),weighttab[,"subject_id"]),"meanw"]
for(featg in feature_spec){
    for(classele in class_num_cols){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[classele,featg]=res$estimate
        pmat[classele,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
colnames(cormat)=str_remove(string=colnames(cormat),pattern="^.+g\\_\\_")
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("spec_morethan1.pdf"))
print(Heatmap(cormat,name="spec",cluster_columns=TRUE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
pdf(paste0("spec_morethan1_subpheno.pdf"))
cormat_subp=cormat[subphenogrp,]
pmat_subp=pmat[subphenogrp,]
print(Heatmap(cormat_subp,name="spec",cluster_columns=TRUE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat_subp[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# Lachnospiraceae family corrleaiton plot
allfeat=colnames(micromat)
f_mask=str_detect(string=allfeat,pattern="f\\_\\_Lachnospiraceae")
g_mask=str_detect(string=allfeat,pattern="\\|g\\_\\_")
s_mask=str_detect(string=allfeat,pattern="\\|s\\_\\_")
# for each genera
feature_spec=unique(allfeat[f_mask&g_mask&(!s_mask)])
cormat=matrix(NA,nrow=nclass,ncol=length(feature_spec))
rownames(cormat)=class_num_cols
colnames(cormat)=feature_spec
pmat=cormat
for(featg in feature_spec){
    for(classele in class_num_cols){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[classele,featg]=res$estimate
        pmat[classele,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
colnames(cormat)=str_remove(string=colnames(cormat),pattern="^.+g\\_\\_")
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("Lachnospiraceae_each_genera.pdf"))
print(Heatmap(cormat,name="genera",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# selected species
select_g_mask=str_detect(string=allfeat,pattern="\\|g\\_\\_((Anaerostipes)|(Butyrivibrio)|(Fusicatenibacter)|(Tyzzerella)|(Anaerotignum)|(Lachnospira))\\|")
feature_spec=unique(allfeat[f_mask&s_mask&select_g_mask])
cormat=matrix(NA,nrow=nclass,ncol=length(feature_spec))
rownames(cormat)=class_num_cols
colnames(cormat)=feature_spec
pmat=cormat
for(featg in feature_spec){
    for(classele in class_num_cols){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[classele,featg]=res$estimate
        pmat[classele,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
colnames(cormat)=str_remove(string=colnames(cormat),pattern="^.+s\\_\\_")
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("Lachnospiraceae_select_species.pdf"))
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# plot pathway associations
tabstat_cor_sub=tabstat_cor[tabstat_cor[,"clinic"]%in%subphenogrp&tabstat_cor[,"padj"]<0.1,]
feature_spec=unique(tabstat_cor_sub[!(tabstat_cor_sub[,"type"]%in%c("species_abs_g","alpha_div")),"feature"])
cormat=matrix(NA,nrow=nclass,ncol=length(feature_spec))
rownames(cormat)=class_num_cols
colnames(cormat)=feature_spec
pmat=cormat
wvec=weighttab[match(rownames(metadatashow),weighttab[,"subject_id"]),"meanw"]
for(featg in feature_spec){
    for(classele in class_num_cols){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[classele,featg]=res$estimate
        pmat[classele,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
pmat=t(pmat)
cormat=t(cormat)
# 
idmatchtab=Reduce(rbind,idname_match_list)
rownamevec=rownames(cormat)
matchind=match(rownamevec,idmatchtab[,1])
notnaind=which(!is.na(matchind))
idmattab=idmatchtab[matchind[notnaind],]
rownamevec[notnaind]=idmattab[,2]
rownames(cormat)=rownamevec
rownames(pmat)=rownamevec
# repeat terms keep one
clists=c("Glycine oxidase","3.1.3.89")
for(crep in clists){
    repmask=str_detect(string=rownames(cormat),pattern=fixed(crep))
    repind=c(which(!repmask),which(repmask)[1])
    cormat=cormat[repind,]
    pmat=pmat[repind,]
    idmattab=idmattab[repind,]
}
# remove terms AMBIGUOUS
remmask=str_detect(string=rownames(cormat),pattern=fixed("AMBIGUOUS"))
cormat=cormat[!remmask,]
pmat=pmat[!remmask,]
idmattab=idmattab[!remmask,]
# 
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("pathway_asso_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# plot pathways among species 
thres_sele=0.1
shownames=str_extract(string=rownamevec,pattern="^.*PWY.*$")
shownames=shownames[!is.na(shownames)]
shownames=c(idmattab[,1],shownames)
feature_spec_species=c()
allfeat=colnames(micromat)
for(feat in shownames){
    featselecspec=allfeat[str_detect(string=allfeat,pattern=fixed(paste0(feat,"|")))]
    quanvec=apply(micromat[,featselecspec,drop=FALSE],2,mean)
    nsele=ceiling(length(quanvec)*thres_sele)
    featdom=names(sort(quanvec,decreasing=TRUE)[seq(nsele)])
    feature_spec_species=c(feature_spec_species,featdom)
}
cormat=matrix(NA,nrow=nclass,ncol=length(feature_spec_species))
rownames(cormat)=class_num_cols
colnames(cormat)=feature_spec_species
pmat=cormat
wvec=weighttab[match(rownames(metadatashow),weighttab[,"subject_id"]),"meanw"]
for(featg in feature_spec_species){
    for(classele in class_num_cols){
        xvec=metadatashow[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[classele,featg]=res$estimate
        pmat[classele,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
pmat=t(pmat)
cormat=t(cormat)
# 
rownamevec=rownames(cormat)
for(rowi in seq(nrow(idmattab))){
    rownamevec=str_replace_all(string=rownamevec,pattern=fixed(idmattab[rowi,1]),replacement=idmattab[rowi,2])
}
rownames(cormat)=rownamevec
rownames(pmat)=rownamevec
# 
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("pathway_asso_spec_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()