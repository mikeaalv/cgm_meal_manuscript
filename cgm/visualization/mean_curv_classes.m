% cgm mean curve plot between metabolic subtypes
close all;
clear all;
% ADD PATH
comp='/Users/yuewu/';%the computer user location
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox=[comp 'Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/'];
addpath(genpath(localPaths.public_toolbox));
% 
pardir=[comp '/Library/CloudStorage/Box-Box/Yue Wu''s Files/cgm_meal_project/'];
workdir=[pardir 'result/cgm_meal/'];%%the working folder
cd(workdir);
% read metadata
metadata=readtable([pardir,'data/metadata/metadata_clean_all_ver2.tsv'],'FileType','text','Delimiter',' ','Format','auto');%,'VariableNamingRule','preserve'
grouplist={'a1c_t2d_status_bl','sspg_status_heyjun','di_3classes_heyjun','ie_3_classes_heyjun','FFA_3classes_heyjun','hepatic_ir_3classes_heyjun'};
metadata=metadata(:,['study_id',grouplist]);
metadata{:,'study_id'}=strip(strrep(metadata{:,'study_id'},'STUDYID-',''),'left','0');
for groupfeat=grouplist
    metadata(strcmp(metadata{:,groupfeat},'Unknown'),groupfeat)={'NA'};
end
% merge groups to improve sample size
metadata{strcmp(metadata{:,'a1c_t2d_status_bl'},'T2D'),'a1c_t2d_status_bl'}={'preDM'};
metadata{strcmp(metadata{:,'FFA_3classes_heyjun'},'Adipocyte_IR'),'FFA_3classes_heyjun'}={'Adipocyte_Intermediate'};
% 
load('fpca_input.mat');
metadata2=metadata;
selecind=cellfun(@(x) any(strcmp(x,subject_vec)),metadata2{:,'study_id'});
metadata2=metadata2(selecind,:);
colorlabels={'Rice','Potatoes','Pasta','Grapes','Bread','Berries','Beans','Rice+Fat','Rice+Fiber','Rice+Protein'};
colors=hsv(length(colorlabels)-3);
colors=[colors; repmat(colors(1,:),[3,1])];
linelabels=[repmat({'-'},[1,7]),':','-.','--'];
mask2=ismember(food_and_metigator_vec,colorlabels);
foodsele=food_and_metigator_vec(mask2);
matsele=ymat_shift(:,mask2);
subject_vec_sele=subject_vec(mask2);
for grpname=grouplist
    groupvec=metadata2{:,grpname};
    groupeles=unique(groupvec(~strcmp(groupvec,'NA')));
    h=figure();
    for i=1:length(groupeles)
        groupsep=groupeles(i);
        matchind=find(strcmp(groupvec,groupsep));
        groupids=metadata2{matchind,"study_id"};
        subind=cellfun(@(x) any(strcmp(x,groupids)),subject_vec_sele);
        subplot(1,3,i);
        hold on;
        hlist={};
        i=1;
        for foodc=colorlabels
            foodc=foodc{1};
            food_ind=strcmp(foodsele,foodc);
            locmat=matsele(:,food_ind&subind);
            if strcmp(foodc,"Rice")
                [grpname, groupsep]
                unique(subject_vec_sele(food_ind&subind))
            end
            meanvec=mean(locmat,2);
            styleid=strcmp(colorlabels,foodc);
            p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
            ylim([80 200]);
            hlist{i}=p(1);
            i=i+1;
        end
        legend([hlist{:}],colorlabels);
        title(groupsep{1});
    end
    sgtitle(grpname);
    saveas(h,[workdir,grpname{1},'curv.mean.fig']);
    close all;
end
% target plot di beta cell funiton wiht normal combine with intermediate
grpname='di_3classes_heyjun';
groupvec=metadata2{:,grpname};
groupvec(strcmp(groupvec,'BC_intermediate'))={'BC_normal'};
metadata2.("di_2classes_heyjun")=groupvec;
groupeles=unique(groupvec(~strcmp(groupvec,'NA')));
h=figure();
for i=1:length(groupeles)
    groupsep=groupeles(i);
    matchind=find(strcmp(groupvec,groupsep));
    groupids=metadata2{matchind,"study_id"};
    subind=cellfun(@(x) any(strcmp(x,groupids)),subject_vec_sele);
    subplot(1,3,i);
    hold on;
    hlist={};
    i=1;
    for foodc=colorlabels
        foodc=foodc{1};
        food_ind=strcmp(foodsele,foodc);
        locmat=matsele(:,food_ind&subind);
        meanvec=mean(locmat,2);
        styleid=strcmp(colorlabels,foodc);
        p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
        ylim([80 200]);
        hlist{i}=p(1);
        i=i+1;
    end
    legend([hlist{:}],colorlabels);
    title(groupsep{1});
end
sgtitle(grpname);
saveas(h,[workdir,grpname,'curv.mean_2group.fig']);
close all;
% plot intersection group predm+is vs normal+ir
plotgroup={{"Normal","IR"},{"preDM","IS"}};
compt={"a1c_t2d_status_bl","sspg_status_heyjun"};
h=figure();
i=1;
for pg=plotgroup
    pg=pg{1};
    grpname=[pg{1},pg{2}];
    matchind=find(strcmp(metadata2{:,compt{1}},pg{1})&strcmp(metadata2{:,compt{2}},pg{2}));
    groupids=metadata2{matchind,"study_id"};
    subind=cellfun(@(x) any(strcmp(x,groupids)),subject_vec_sele);
    subplot(1,2,i);
    hold on;
    hlist={};
    j=1;
    for foodc=colorlabels
        foodc=foodc{1};
        food_ind=strcmp(foodsele,foodc);
        locmat=matsele(:,food_ind&subind);
        if strcmp(foodc,"Rice")
            grpname
            unique(subject_vec_sele(food_ind&subind))
        end
        meanvec=mean(locmat,2);
        styleid=strcmp(colorlabels,foodc);
        p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
        ylim([80 200]);
        hlist{j}=p(1);
        j=j+1;
    end
    legend([hlist{:}],colorlabels);
    title(grpname);
    i=i+1;
end
saveas(h,[workdir,'mix_curv.mean.fig']);
close all;
% taregtted mitigator plot 
mitilist={'Rice+Fat','Rice+Fiber','Rice+Protein'};
plotlist=[{'Rice'},mitilist];
withmitig=cellfun(@(x) any(strcmp(x,mitilist)),foodsele);
subwithmitig=unique(subject_vec_sele(withmitig));
withind=cellfun(@(x) any(strcmp(x,subwithmitig)),metadata2{:,"study_id"});
metadata3=metadata2(withind,:);
grouplist={'a1c_t2d_status_bl','sspg_status_heyjun','di_2classes_heyjun','ie_3_classes_heyjun','FFA_3classes_heyjun','hepatic_ir_3classes_heyjun'};
for grpname=grouplist
    groupvec=metadata3{:,grpname};
    groupeles=unique(groupvec(~strcmp(groupvec,'NA')));
    h=figure();
    for i=1:length(groupeles)
        groupsep=groupeles(i);
        matchind=find(strcmp(groupvec,groupsep));
        groupids=metadata3{matchind,"study_id"};
        subind=cellfun(@(x) any(strcmp(x,groupids)),subject_vec_sele);
        subplot(1,3,i);
        hold on;
        hlist={};
        i=1;
        for foodc=plotlist
            foodc=foodc{1};
            food_ind=strcmp(foodsele,foodc);
            locmat=matsele(:,food_ind&subind);
            if strcmp(foodc,"Rice+Fiber")
                [grpname, groupsep]
                unique(subject_vec_sele(food_ind&subind))
            end
            meanvec=mean(locmat,2);
            styleid=strcmp(plotlist,foodc);
            p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
            ylim([80 200]);
            hlist{i}=p(1);
            i=i+1;
        end
        legend([hlist{:}],plotlist);
        title(groupsep{1});
    end
    sgtitle(grpname);
    saveas(h,[workdir,grpname{1},'curv.mean.mititarg.fig']);
    close all;
end
% 
h=figure();
i=1;
for pg=plotgroup
    pg=pg{1};
    grpname=[pg{1},pg{2}];
    matchind=find(strcmp(metadata3{:,compt{1}},pg{1})&strcmp(metadata3{:,compt{2}},pg{2}));
    groupids=metadata3{matchind,"study_id"};
    subind=cellfun(@(x) any(strcmp(x,groupids)),subject_vec_sele);
    subplot(1,2,i);
    hold on;
    hlist={};
    j=1;
    for foodc=plotlist
        foodc=foodc{1};
        food_ind=strcmp(foodsele,foodc);
        locmat=matsele(:,food_ind&subind);
        if strcmp(foodc,"Rice+Fiber")
            grpname
            unique(subject_vec_sele(food_ind&subind))
        end
        meanvec=mean(locmat,2);
        styleid=strcmp(plotlist,foodc);
        p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
        ylim([80 200]);
        hlist{j}=p(1);
        j=j+1;
    end
    legend([hlist{:}],plotlist);
    title(grpname);
    i=i+1;
end
saveas(h,[workdir,'mix_curv.mean.mitig.fig']);
close all;