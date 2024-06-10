# food spikes group wise comparison of microbiome
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
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_meal_project/")
resdir=paste0(pardir,"result/microbiome/spikes/")
datadir=paste0(pardir,"result/microbiome/eda/")
setwd(resdir)
# spike data 
load(paste0(pardir,"result/cgm_meal/heatmap_plot.RData"))
# microbiome data
load(paste0(datadir,"microbiome_processed.RData"))
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
# select sig change features
colstore=colnames(micromat)
colnames(micromat)=NULL
groupvec=c("Bread","Grapes","Potatoes","Rice","Mitigator")
testtab=merge(micromat,mat_worstfood[,groupvec[seq(4)]],by="row.names")
rownames(testtab)=testtab[,1]
testtab=testtab[,-1]
mitig_tab=as.data.frame(ifelse(rowMeans(mat_mitigator_effect)>0,1,0))
colnames(mitig_tab)="Mitigator"
testtab=merge(testtab,mitig_tab,by="row.names")
rownames(testtab)=testtab[,1]
testtab=testtab[,-1]
subjects=rownames(testtab)
colnames(testtab)[seq(length(colstore))]=colstore
nfoodsig=length(groupvec)
foodtypes=groupvec
nsubjects=length(subjects)
# lm test result table 
obslist=c("go","kegg","pathway","pfam","rxn","species_abs_g","alpha_div")
wvec=weighttab[match(rownames(testtab),weighttab[,"subject_id"]),"meanw"]
testtab$weight=wvec
progress<-function(n) cat(sprintf("task %d is complete\n", n))
opts<-list(progress=progress)
stattab_list_comb<-foreach(foodsigi=seq(nfoodsig),.options.snow=opts)%dopar%{
    temptab=c()
    for(typ in obslist){
        foodt=foodtypes[foodsigi]
        xvec=testtab[,foodt]
        nonaind=(!is.na(xvec))
        testtab_temp=testtab[nonaind,]
        if(min(table(testtab_temp[,foodt]))<=3){
            next
        }
        features=summfeat_list[[typ]]
        # select features
        if(!(typ%in%c("species_abs_g","alpha_div"))){
            misvec=missing_ratio_list[[typ]]
            pertvec=avg_perc_list[[typ]]
            features=intersect(names(misvec)[misvec<0.5],features)
        }
        micromatloc=testtab_temp[,features]
        ind_dochange=apply(micromatloc,2,var)!=0
        slopvec=c()
        pvec=c()
        featvec=c()
        for(feat in features[ind_dochange]){
            locdf=testtab[,c(feat,foodt,"weight")]
            colnames(locdf)=c("value","group","weight")
            locmat=as.matrix(locdf)
            y=locmat[,1]
            X=locmat[,2:3]
            lmodel=.lm.fit(x=X,y=y)
            # 
            rss=sum(lmodel$residuals^2)
            rdf=length(y)-ncol(X)
            resvar=rss/rdf
            R=chol2inv(lmodel$qr)
            se=sqrt(diag(R)*resvar)
            pval=2*pt(abs(lmodel$coef/se),rdf,lower.tail=FALSE)
            # pval=wilcox.test(value~group,data=locdf)$p.value
            pvec=c(pvec,pval[1])
            slopvec=c(slopvec,lmodel$coefficients[1])
            featvec=c(featvec,feat)
        }
        # padj=rep(NA,times=length(pvec))
        # namsk=is.na(pvec)
        padj=p.adjust(pvec,method="fdr")
        sigmask=which(pvec<0.05)
        len=length(sigmask)
        temptab=rbind(temptab,data.frame(food=rep(foodtypes[foodsigi],times=len),feature=featvec[sigmask],omics=rep(typ,times=len),slope=slopvec[sigmask],pval=pvec[sigmask],padj=padj[sigmask]))
    }
    temptab
}
stattab_lm=Reduce(rbind,stattab_list_comb)
save(stattab_lm,file="test_spike_microbiome_grpcomp.RData")