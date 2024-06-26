# associate general spikes and mitigator effect with microbiomes
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
require(ppcor)
require(data.table)
# 
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
resdir=paste0(pardir,"result/pathway/")
setwd(resdir)
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
# microbiomes
load(paste0(pardir,"result/microbiome/eda/microbiome_processed.RData"))
subjects=intersect(rownames(micromat),rownames(mat_quantile))
mat_quantile=mat_quantile[subjects,]
mat_mitigator_effect=mat_mitigator_effect[subjects,]
micromat=micromat[subjects,]
nsubjects=length(subjects)
# 
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
featvec=unlist(summfeat_list)
# 
obslist=c("go","kegg","pathway","pfam","rxn","species_abs_g","alpha_div")
listvec=list(spike=c(mat_quantile),mitigation=c(mat_mitigator_effect))
listrep=list(spike=dim(mat_quantile)[2],mitigation=dim(mat_mitigator_effect)[2])
wvec=weighttab[match(rownames(micromat),weighttab[,"subject_id"]),"meanw"]
cor_mat=matrix(NA,nrow=length(featvec),ncol=2)
colnames(cor_mat)=names(listvec)
rownames(cor_mat)=featvec
p_mat=cor_mat
tabstat_cor=c()
for(typ in names(listvec)){
    xvec=listvec[[typ]]
    nonaind=(!is.na(xvec))
    if(length(unique(xvec[nonaind]))<=1){
        next
    }
    for(microtyp in obslist){
        features=summfeat_list[[microtyp]]
        # select features
        if(!(microtyp%in%c("species_abs_g","alpha_div"))){
            misvec=missing_ratio_list[[microtyp]]
            pertvec=avg_perc_list[[microtyp]]
            features=intersect(names(misvec)[misvec<0.5],features)
        }
        corvec=rep(NA,times=length(features))
        names(corvec)=features
        pval=corvec
        micromatloc=micromat[,features]
        ind_dochange=apply(micromatloc,2,var)!=0
        for(feat in features[ind_dochange]){
            yvec=rep(micromatloc[,feat],times=listrep[[typ]])
            wvec_rep=rep(wvec,times=listrep[[typ]])
            res=pcor.test(xvec[nonaind],yvec[nonaind],wvec_rep[nonaind],"spearman")
            pval[feat]=res$p.value
            corvec[feat]=res$estimate
        }
        cor_mat[features,typ]=corvec
        p_mat[features,typ]=pval
        p.adj=p.adjust(pval,"fdr")
        seleind=which(c(pval<0.05))
        seleind=seleind[!is.na(seleind)]
        len=length(seleind)
        tabstat_cor=rbind(tabstat_cor,data.frame(microtype=rep(microtyp,times=len),spikes=rep(typ,times=len),feature=features[seleind],correlation=corvec[seleind],pvalue=pval[seleind],padj=p.adj[seleind]))
    }
}
idmatchtab=Reduce(rbind,idname_match_list)
rownamevec=tabstat_cor[,"feature"]
matchind=match(rownamevec,idmatchtab[,1])
notnaind=which(!is.na(matchind))
idmattab=idmatchtab[matchind[notnaind],]
rownamevec[notnaind]=idmattab[,2]
tabstat_cor$description=rownamevec
save(tabstat_cor,file="association_agonsit_spike_microbiome.RData")
tabsub_mitig=tabstat_cor[tabstat_cor[,"spikes"]=="mitigation",]
write.table(tabsub_mitig,file="microb_ass_mitig.tsv",row.names=FALSE)
# species specific quantification 
feat_select=c("k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Klebsiella","k__Bacteria|p__Proteobacteria|c__Betaproteobacteria|o__Burkholderiales|f__Sutterellaceae|g__Sutterella","PWY-6731: starch degradation III","PF05504","GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)","PF14620","GO:0009311","PF09560","GO:0009409","PF12571","GO:0009408","PF16760","PF01758","k__Bacteria|p__Actinobacteria|c__Coriobacteriia|o__Eggerthellales|f__Eggerthellaceae|g__Gordonibacter","GO:0030259","GO:0015941","GO:0047761")
thres_sele=0.05
feature_spec_species=c()
allfeat=colnames(micromat)
for(feat in feat_select){
    featselecspec=allfeat[str_detect(string=allfeat,pattern=fixed(paste0(feat,"|")))]
    quanvec=apply(micromat[,featselecspec,drop=FALSE],2,mean)
    nsele=ceiling(length(quanvec)*thres_sele)
    featdom=names(sort(quanvec,decreasing=TRUE)[seq(nsele)])
    feature_spec_species=c(feature_spec_species,featdom)
}
cormat=matrix(NA,nrow=1,ncol=length(feature_spec_species))
typ="mitigation"
rownames(cormat)=typ
colnames(cormat)=feature_spec_species
pmat=cormat
wvec=weighttab[match(rownames(micromat),weighttab[,"subject_id"]),"meanw"]
xvec=listvec[[typ]]
for(featg in feature_spec_species){
    nonaind=(!is.na(xvec))
    if(length(unique(xvec[nonaind]))<=1){
        next
    }
    yvec=rep(micromat[,featg],times=listrep[[typ]])
    wvec_rep=rep(wvec,times=listrep[[typ]])
    res=pcor.test(xvec[nonaind],yvec[nonaind],wvec_rep[nonaind],"spearman")
    cormat[typ,featg]=res$estimate
    pmat[typ,featg]=res$p.value
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
pdf(paste0("microb_asso_mitigat_spec_heatmap.pdf"))#,width=15,height=20
print(Heatmap(cormat,name="spec",cluster_columns=FALSE,cluster_rows=FALSE,show_row_names=TRUE,show_column_names=TRUE,col=col_fun,
    cell_fun=function(j,i,x,y,w,h,fill){
            if(pmat[i,j]<0.05){
                grid.text("*",x,y)
            }
        }
))
dev.off()