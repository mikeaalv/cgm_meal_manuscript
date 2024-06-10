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
require(foreach)
require(doSNOW)
require(Rodin)
cl<-makeSOCKcluster(8)
registerDoSNOW(cl)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/");
resdir=paste0(pardir,"result/microbiome/asso_other_omics/")
datadir=paste0(pardir,"result/microbiome/eda/")
setwd(resdir)
# other omics data 
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
load(paste0(pardir,"result/omics/metabolomics_lipidomics_processed.RData"))
load(paste0(pardir,"result/pathway/omics_data_clean.RData"))
# microbiome data
load(paste0(datadir,"microbiome_processed.RData"))
# 
metab_ind=which(omicsvec=="metabolomics"|omicsvec=="lipidomics")
metab_mat=omicsmat[,metab_ind]
featvec_metab=featvec[metab_ind]
# add alpha diveristy data
subjects=rownames(micromat)
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
# subject information intersect
subjects=intersect(rownames(micromat),rownames(omicsmat))
metab_mat=metab_mat[subjects,]
micromat=micromat[subjects,]
nsubjects=length(subjects)
obslist=c("go","kegg","pathway","pfam","rxn","species_abs_g","alpha_div")
wvec=weighttab[match(rownames(micromat),weighttab[,"subject_id"]),"meanw"]
# 
featvec_micro=unlist(summfeat_list)
# 
progress<-function(n) cat(sprintf("task %d is complete\n", n))
opts<-list(progress=progress)
tabstat_cor_list<-foreach(ind=seq(length(featvec_metab)),.packages=c("stringr","ppcor"),.options.snow=opts)%dopar%{
    metafeat=featvec_metab[ind]
    tabstat_cor_tmp=c()
    for(typ in obslist){
        colind=which(featvec_metab==metafeat)
        xvec=metab_mat[,colind]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        microfeatures=summfeat_list[[typ]]
        # select features
        if(!(typ%in%c("species_abs_g","alpha_div"))){
            misvec=missing_ratio_list[[typ]]
            pertvec=avg_perc_list[[typ]]
            microfeatures=intersect(names(misvec)[misvec<0.5],microfeatures)
        }
        corvec=rep(NA,times=length(microfeatures))
        names(corvec)=microfeatures
        pval=corvec
        micromatloc=micromat[nonaind,microfeatures]
        ind_dochange=apply(micromatloc,2,var)!=0
        for(microfeat in microfeatures[ind_dochange]){
            res=pcor.test(xvec[nonaind],micromatloc[,microfeat],wvec[nonaind],"spearman")
            pval[microfeat]=res$p.value
            corvec[microfeat]=res$estimate
        }
        p.adj=p.adjust(pval,"fdr")
        seleind=which(c(pval<0.05))
        seleind=seleind[!is.na(seleind)]
        len=length(seleind)
        tabstat_cor_tmp=rbind(tabstat_cor_tmp,data.frame(type=rep(typ,times=len),metabolites=rep(metafeat,times=len),feature=microfeatures[seleind],correlation=corvec[seleind],pvalue=pval[seleind],padj=p.adj[seleind]))
    }
    tabstat_cor_tmp
}
tabstat_cor=Reduce("rbind",tabstat_cor_list)
# id conversion
idmatchtab=Reduce(rbind,idname_match_list)
rownamevec=tabstat_cor[,"feature"]
matchind=match(rownamevec,idmatchtab[,1])
notnaind=which(!is.na(matchind))
rownamevec[notnaind]=idmatchtab[matchind[notnaind],2]
tabstat_cor$feat_description=rownamevec
save(tabstat_cor,file="association_metabolomics_microbiome.RData")
# correlation matrix heatmap
metabo_list=list(carnitine=c("C4:0 AC"),bileacids=c("Glycochenodeoxycholic acid 3-glucuronide","Hydroxy-cholenoic acid","Ketodeoxycholic acid(2)","Chenodeoxycholic Acid(1)"),other=c("Glucose","Hippuric acid","Citramalic acid","Citric acid","Uric acid"))
microb_list=list(lipids=c("3-OXOACYL-ACP-REDUCT-RXN","3-OXOSTEROID-1-DEHYDROGENASE-RXN","KETOACYLCOATHIOL-RXN"),glucogenesis=c("4.1.1.32-RXN","PWY66-399: gluconeogenesis III"),other=c("K07096","3.2.1.23-RXN","RXN-12398","PF00723"))
metabo_vec=unlist(metabo_list)
microb_vec=unlist(microb_list)
# 
cor_mat=matrix(NA,nrow=length(microb_vec),ncol=length(metabo_vec))
colnames(cor_mat)=metabo_vec
rownames(cor_mat)=microb_vec
p_mat=cor_mat
wvec=weighttab[match(rownames(micromat),weighttab[,"subject_id"]),"meanw"]
for(featg in microb_vec){
    for(featmeta in metabo_vec){
        colind=which(featvec_metab==featmeta)
        xvec=metab_mat[,colind]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cor_mat[featg,featmeta]=res$estimate
        p_mat[featg,featmeta]=res$p.value
    }
}
cor_mat_sele=cor_mat[microb_vec,metabo_vec]
p_mat_sele=p_mat[microb_vec,metabo_vec]
micro_feat_names=tabstat_cor[match(microb_vec,tabstat_cor[,"feature"]),"feat_description"]
rownames(cor_mat_sele)=micro_feat_names
rownames(p_mat_sele)=micro_feat_names
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("micro_metab_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cor_mat_sele,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(p_mat_sele[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
# plot pathways among species 
thres_sele=0.1
shownames=microb_vec
feature_spec_species=c()
allfeat=colnames(micromat)
for(feat in shownames){
    featselecspec=allfeat[str_detect(string=allfeat,pattern=fixed(paste0(feat,"|")))]
    quanvec=apply(micromat[,featselecspec,drop=FALSE],2,mean)
    nsele=ceiling(length(quanvec)*thres_sele)
    featdom=names(sort(quanvec,decreasing=TRUE)[seq(nsele)])
    feature_spec_species=c(feature_spec_species,featdom)
}
cormat=matrix(NA,nrow=length(metabo_vec),ncol=length(feature_spec_species))
rownames(cormat)=metabo_vec
colnames(cormat)=feature_spec_species
pmat=cormat
for(featg in feature_spec_species){
    for(featmeta in metabo_vec){
        colind=which(featvec_metab==featmeta)
        xvec=metab_mat[,colind]
        nonaind=(!is.na(xvec))
        if(length(unique(xvec[nonaind]))<=1){
            next
        }
        res=pcor.test(xvec[nonaind],micromat[nonaind,featg],wvec[nonaind],"spearman")
        cormat[featmeta,featg]=res$estimate
        pmat[featmeta,featg]=res$p.value
    }
}
pmat[is.na(pmat)]=1
pmat=t(pmat)
cormat=t(cormat)
# 
rownamevec=rownames(cormat)
micro_feat_names=tabstat_cor[match(microb_vec,tabstat_cor[,"feature"]),"feat_description"]
names(micro_feat_names)=microb_vec
for(prename in microb_vec){
    rownamevec=str_replace_all(string=rownamevec,pattern=fixed(prename),replacement=micro_feat_names[prename])
}
rownames(cormat)=rownamevec
rownames(pmat)=rownamevec
# 
col_fun=colorRamp2(c(-1,0,1),c("blue","white","red"))
pdf(paste0("micro_metab_spec_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()
save(list=ls(all.names=TRUE),file="meta_micro_ass_env.RData")