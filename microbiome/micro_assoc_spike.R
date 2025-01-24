# partial correlation between food spikes and microbiome
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
resdir=paste0(pardir,"result/microbiome/spikes/")
datadir=paste0(pardir,"result/microbiome/eda/")
setwd(resdir)
# spike data 
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
# microbiome data
load(paste0(datadir,"microbiome_processed.RData"))
# 
food_sig_mat=as.data.frame(cbind(mat_quantile,mat_mitigator_effect))
# potato vs grapes
# food_sig_mat$"potato_vs_grape"=unlist(food_sig_mat[,"Potatoes"]/food_sig_mat[,"Grapes"])
# 
nfoodsig=dim(food_sig_mat)[2]
foodtypes=colnames(food_sig_mat)
# subject information intersect
subjects=intersect(rownames(micromat),rownames(food_sig_mat))
food_sig_mat=food_sig_mat[subjects,]
micromat=micromat[subjects,]
nsubjects=length(subjects)
# 
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
# correlation result table 
obslist=c("go","kegg","pathway","pfam","rxn","species_abs_g","alpha_div","species_abs_s")
tabstat_cor=c()
wvec=weighttab[match(rownames(micromat),weighttab[,"subject_id"]),"meanw"]
# 
featvec=unlist(summfeat_list)
cor_mat=matrix(NA,nrow=length(featvec),ncol=length(foodtypes))
colnames(cor_mat)=foodtypes
rownames(cor_mat)=featvec
p_mat=cor_mat
for(classele in foodtypes){
    for(typ in obslist){
        xvec=food_sig_mat[,classele]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        features=summfeat_list[[typ]]
        # select features
        if(!(typ%in%c("species_abs_g","species_abs_s","alpha_div"))){
            misvec=missing_ratio_list[[typ]]
            pertvec=avg_perc_list[[typ]]
            features=intersect(names(misvec)[misvec<0.5],features)
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
        cor_mat[features,classele]=corvec
        p_mat[features,classele]=pval
        p.adj=p.adjust(pval,"fdr")
        seleind=which(c(pval<0.05))
        seleind=seleind[!is.na(seleind)]
        len=length(seleind)
        tabstat_cor=rbind(tabstat_cor,data.frame(type=rep(typ,times=len),spikes=rep(classele,times=len),feature=features[seleind],correlation=corvec[seleind],pvalue=pval[seleind],padj=p.adj[seleind]))
    }
}
save(tabstat_cor,file="association_spike_microbiome.RData")
# heatmap cutoff by padj
siglist=unique(tabstat_cor[tabstat_cor[,"padj"]<0.1,"feature"])
featselect=siglist
p_mat_sele=p_mat[featselect,colnames(mat_quantile)]
cor_mat_sele=cor_mat[featselect,colnames(mat_quantile)]
# id conversion
idmatchtab=Reduce(rbind,idname_match_list)
rownamevec=rownames(cor_mat_sele)
matchind=match(rownamevec,idmatchtab[,1])
notnaind=which(!is.na(matchind))
idmattab=idmatchtab[matchind[notnaind],]
rownamevec[notnaind]=idmattab[,2]
rownames(cor_mat_sele)=rownamevec
rownames(p_mat_sele)=rownamevec
# 
p_mat_sele_fdr=matrix(1,nrow=length(rownamevec),ncol=dim(mat_quantile)[2])
rownames(p_mat_sele_fdr)=rownamevec
colnames(p_mat_sele_fdr)=colnames(cor_mat_sele)
for(foods in colnames(cor_mat_sele)){
    for(feat in rownamevec){
        padj=tabstat_cor[tabstat_cor[,"feature"]==feat&tabstat_cor[,"spikes"]==foods,"padj"]
        if(length(padj)>0){
            p_mat_sele_fdr[feat,foods]=padj
        }
    }
}
#
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("heatmap_spike_micro_padjcut.pdf"))#,width=10,height=20
print(Heatmap(cor_mat_sele,name="correlation with micro",cluster_columns=TRUE,cluster_rows=TRUE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat_sele[i,j]<0.05){
                if(p_mat_sele_fdr[i,j]<0.2){
                    grid.text("**",x,y)
                }else{
                    grid.text("*",x,y)
                }
            }
        }
))
dev.off()
savdf=cor_mat_sele
write.table(savdf,file=paste0("figext7.txt"))
# plot pathways among species 
thres_sele=0.05
feature_spec_species=c()
allfeat=colnames(micromat)
shownames=str_extract(string=rownamevec,pattern="^.*PWY.*$")
shownames=shownames[!is.na(shownames)]
shownames=c(idmattab[,1],shownames)
for(feat in shownames){
    featselecspec=allfeat[str_detect(string=allfeat,pattern=fixed(paste0(feat,"|")))]
    quanvec=apply(micromat[,featselecspec,drop=FALSE],2,mean)
    nsele=ceiling(length(quanvec)*thres_sele)
    featdom=names(sort(quanvec,decreasing=TRUE)[seq(nsele)])
    feature_spec_species=c(feature_spec_species,featdom)
}
cormat=matrix(NA,nrow=nfoodsig,ncol=length(feature_spec_species))
rownames(cormat)=foodtypes
colnames(cormat)=feature_spec_species
pmat=cormat
wvec=weighttab[match(rownames(food_sig_mat),weighttab[,"subject_id"]),"meanw"]
for(featg in feature_spec_species){
    for(classele in foodtypes){
        xvec=food_sig_mat[,classele]
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
pdf(paste0("heatmap_spike_micro_spec_padjcut.pdf"),width=15,height=10)#,
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()