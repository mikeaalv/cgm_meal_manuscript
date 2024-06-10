% qc check of the processed CGM data
close all;
clear all;
% ADD PATH
comp='/Users/yuewu/';%the computer user location
% Metabolic toolbox toolbox found @  https://github.com/artedison/Edison_Lab_Shared_Metabolomics_UGA
localPaths.public_toolbox=[comp 'Documents/GitHub/Edison_Lab_Shared_Metabolomics_UGA/'];
% some wrappers for fdaM @ https:/fd/github.com/mikeaalv/fda_learn
localPaths.fdalearn=[comp 'Documents/GitHub/fda_learn/'];
% functiona data analysis matlab package fdaM @ https://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/
localPaths.fdam=[comp 'Documents/MATLAB/fdaM/'];
%
addpath(genpath(localPaths.public_toolbox));
addpath(genpath(localPaths.fdalearn));
addpath(genpath(localPaths.fdam));
% 
pardir=[comp '/Library/CloudStorage/Box-Box/Yue Wu''s Files/cgm_meal_project/'];
workdir=[pardir 'result/cgm_meal/'];%%the working folder
cd(workdir);
% load the dataset
datatab=readtable([workdir,'cgm_foods_smooth.txt'],'Format','auto');
% 
subjects=unique(datatab{:,'subject'})';
food_and_metigator=unique(datatab{:,'foods'})';
food=unique(datatab{:,'food'})';
mitigator=unique(datatab{:,'mitigator'})';
reps=unique(datatab{:,'rep'})';
subject_vec=[];
food_vec={};
mitigator_vec={};
food_and_metigator_vec={};
rep_vec=[];
ymat=[];
for subele=subjects
    for foodcomb=food_and_metigator
        for repele=reps
            matchind=find(datatab{:,'subject'}==subele & strcmp(datatab{:,'foods'},foodcomb) & datatab{:,'rep'}==repele);
            if length(matchind)==0
                continue;
            end
            loctab=datatab(matchind,:);
            timvec=loctab{:,'mins_since_start'};
            [timvecsort,sortind]=sort(timvec);
            ymat=[ymat loctab{sortind,'glucose'}];
            subject_vec=[subject_vec subele];
            food_vec=[food_vec unique(loctab{:,'food'})];
            mitigator_vec=[mitigator_vec unique(loctab{:,'mitigator'})];
            food_and_metigator_vec=[food_and_metigator_vec unique(loctab{:,'foods'})];
            rep_vec=[rep_vec repele];
        end
    end
end
% curve vis of rice mitigators
mask2=contains(food_and_metigator_vec,'Rice');
colors={'b','r','g','y'};
colorlabels={'Rice','Rice+Fat','Rice+Fiber','Rice+Protein'};
subject_vec=arrayfun(@num2str,subject_vec,'UniformOutput',false);
for subj=unique(subject_vec)
    subj=subj{1};
    mask1=strcmp(subject_vec,subj);
    foodsloca=food_and_metigator_vec(mask1&mask2);
    matloc=ymat(:,mask1&mask2);
    h=figure();
    hold on;
    hlist={};
    i=1;
    for comb=unique(foodsloca)
        comb=comb{1};
        submask=strcmp(foodsloca,comb);
        p=plot(timvecsort,matloc(:,submask),colors{strcmp(colorlabels,comb)},'LineWidth',2);
        hlist{i}=p(1);
        i=i+1;
    end
    legend([hlist{:}],unique(foodsloca));
    title(subj);
    xlabel("time (mins)");
    ylabel("glucose (mg/dL)");
    saveas(h,[workdir,'curv.' subj '.fig']);
    saveas(h,[workdir,'curv.' subj '.dircsave.pdf']);
    close all;
end
% 
close all;
subj_uniq_vec=unique(subject_vec);
[~,ind]=sort(str2double(subj_uniq_vec));
subj_uniq_vec=subj_uniq_vec(ind);
axlist={};
for subi=1:length(subj_uniq_vec)
    h1=openfig([workdir,'curv.' subj_uniq_vec{subi} '.fig'],'reuse');
    figure(subi)
    axlist{subi}=gca;
end
htot=figure();
for subi=1:length(subj_uniq_vec)
    s1=subplot(11,5,subi);
    fig1=get(axlist{subi},'children');
    copyobj(fig1,s1);
    title(subj_uniq_vec{subi});
end
saveas(htot,[workdir,'curv.tot.fig']);
close all;
% all basic food
colorlabels={'Rice','Potatoes','Pasta','Grapes','Bread','Berries','Beans'};
colors=hsv(length(colorlabels));
mask2=ismember(food_and_metigator_vec,colorlabels);
for subj=unique(subject_vec)
    subj=subj{1};
    mask1=strcmp(subject_vec,subj);
    foodsloca=food_and_metigator_vec(mask1&mask2);
    matloc=ymat(:,mask1&mask2);
    h=figure();
    hold on;
    hlist={};
    i=1;
    for comb=unique(foodsloca)
        comb=comb{1};
        submask=strcmp(foodsloca,comb);
        p=plot(timvecsort,matloc(:,submask),'Color',colors(strcmp(colorlabels,comb),:),'LineWidth',2);
        hlist{i}=p(1);
        i=i+1;
    end
    legend([hlist{:}],unique(foodsloca));
    title(subj);
    saveas(h,[workdir,'curv.' subj '.basicfoods.fig']);
    saveas(h,[workdir,'curv.' subj '.basicfoods.dircsave.pdf']);
    close all;
end
% 
close all;
subj_uniq_vec=unique(subject_vec);
[~,ind]=sort(str2double(subj_uniq_vec));
subj_uniq_vec=subj_uniq_vec(ind);
axlist={};
for subi=1:length(subj_uniq_vec)
    h1=openfig([workdir,'curv.' subj_uniq_vec{subi} '.basicfoods.fig'],'reuse');
    figure(subi)
    axlist{subi}=gca;
end
htot=figure();
for subi=1:length(subj_uniq_vec)
    s1=subplot(11,5,subi);
    fig1=get(axlist{subi},'children');
    copyobj(fig1,s1);
    title(subj_uniq_vec{subi});
end
saveas(htot,[workdir,'curv.tot.basicfoods.fig']);
close all;
% curve vis of other mitigators
otherfoodsmit={'Bread','Grapes','Potatoes','Pasta'};
colors={'b','r','g','y'};
for otherf=otherfoodsmit
    otherf=otherf{1};
    mask2=contains(food_and_metigator_vec,otherf);
    colorlabels={otherf,[otherf '+Fat'],[otherf '+Fiber'],[otherf '+Protein']};
    for subj=unique(subject_vec)
        subj=subj{1};
        mask1=strcmp(subject_vec,subj);
        foodsloca=food_and_metigator_vec(mask1&mask2);
        if ~any(contains(foodsloca,'+'))
            continue
        end
        matloc=ymat(:,mask1&mask2);
        h=figure();
        hold on;
        hlist={};
        i=1;
        for comb=unique(foodsloca)
            comb=comb{1};
            submask=strcmp(foodsloca,comb);
            p=plot(timvecsort,matloc(:,submask),colors{strcmp(colorlabels,comb)},'LineWidth',2);
            hlist{i}=p(1);
            i=i+1;
        end
        legend([hlist{:}],unique(foodsloca));
        title(subj);
        xlabel("time (mins)");
        ylabel("glucose (mg/dL)");
        saveas(h,[workdir,'curv.' subj '.' otherf '.mitigator.fig']);
        close all;
    end
end
% shift 0
modind=find(timvecsort==0);
subj_uniq=unique(subject_vec);
ymat_shift=ymat;
for subj=subj_uniq
    subj=subj{1};
    inds=find(strcmp(subject_vec,subj));
    % inds=find(subject_vec==subj);
    locmat=ymat_shift(:,inds);
    vec0=locmat(modind,:);
    delt0=mean(vec0)-vec0;
    locmat=locmat+delt0;
    ymat_shift(:,inds)=locmat;
end
% plot after shift 
% curve vis of rice mitigators
mask2=contains(food_and_metigator_vec,'Rice');
colors={'b','r','g','y'};
colorlabels={'Rice','Rice+Fat','Rice+Fiber','Rice+Protein'};
for subj=unique(subject_vec)
    subj=subj{1};
    mask1=strcmp(subject_vec,subj);
    foodsloca=food_and_metigator_vec(mask1&mask2);
    matloc=ymat_shift(:,mask1&mask2);
    h=figure();
    hold on;
    hlist={};
    i=1;
    for comb=unique(foodsloca)
        comb=comb{1};
        submask=strcmp(foodsloca,comb);
        p=plot(timvecsort,matloc(:,submask),colors{strcmp(colorlabels,comb)},'LineWidth',2);
        hlist{i}=p(1);
        i=i+1;
    end
    legend([hlist{:}],unique(foodsloca));
    title(subj);
    xlabel("time (mins)");
    ylabel("glucose (mg/dL)");
    saveas(h,[workdir,'curv.' subj 'shift.fig']);
    close all;
end
% all basic food
colorlabels={'Rice','Potatoes','Pasta','Grapes','Bread','Berries','Beans'};
colors=hsv(length(colorlabels));
mask2=ismember(food_and_metigator_vec,colorlabels);
for subj=unique(subject_vec)
    subj=subj{1};
    mask1=strcmp(subject_vec,subj);
    foodsloca=food_and_metigator_vec(mask1&mask2);
    matloc=ymat_shift(:,mask1&mask2);
    h=figure();
    hold on;
    hlist={};
    i=1;
    for comb=unique(foodsloca)
        comb=comb{1};
        submask=strcmp(foodsloca,comb);
        p=plot(timvecsort,matloc(:,submask),'Color',colors(strcmp(colorlabels,comb),:),'LineWidth',2);
        hlist{i}=p(1);
        i=i+1;
    end
    legend([hlist{:}],unique(foodsloca));
    title(subj);
    saveas(h,[workdir,'curv.' subj '.basicfoods.shift.fig']);
    close all;
end
save('fpca_input.mat','ymat','subject_vec','food_vec','mitigator_vec','food_and_metigator_vec','rep_vec','timvecsort','ymat_shift');
% sample quality 
vec_var_mat=[];
foodchecks={'Rice','Bread','Potatoes','Pasta'}
for foodc=foodchecks
    foodc=foodc{1};
    vec_var=[];
    for subj=subj_uniq
        rice_ind=strcmp(food_and_metigator_vec,foodc);
        indi_ind=strcmp(subject_vec,subj);
        % indi_ind=subject_vec==subj;
        locmat=ymat_shift(:,rice_ind&indi_ind);
        meanvec=mean(locmat,2);
        distmat=abs((locmat-meanvec)/mean(meanvec));
        [nx,ny]=size(distmat);
        fillval=sum(distmat(:))/nx/ny;
        if size(locmat,2)<=1
            fillval=NaN;
        end
        vec_var=[vec_var fillval];
    end
    vec_var_mat=[vec_var_mat; vec_var];
end
h=heatmap(subj_uniq,foodchecks,[vec_var_mat]);
saveas(h,[workdir,'qcheatmap.fig']);
close all;
vec_var_tab=array2table(vec_var_mat,'VariableNames',subj_uniq,'RowNames',foodchecks);
writetable(vec_var_tab,'qcheatmap.txt','WriteRowNames',true);
% food variances
vec_var_food_mean=[];
vec_var_food_se=[];
foodvec=unique(food_and_metigator_vec);
foodvec_sele={};
for foodc=foodvec
    foodc=foodc{1};
    vec_var_food=[];
    for subj=subj_uniq
        indi_ind=strcmp(subject_vec,subj);
        food_ind=strcmp(food_and_metigator_vec,foodc);
        locmat=ymat_shift(:,food_ind&indi_ind);
        meanvec=mean(locmat,2);
        distmat=abs((locmat-meanvec)/mean(meanvec));
        [nx,ny]=size(distmat);
        fillval=sum(distmat(:))/nx/ny;
        if size(locmat,2)<=1
            fillval=NaN;
        end
        vec_var_food=[vec_var_food fillval];
    end
    vec_var_food=vec_var_food(~isnan(vec_var_food));
    if length(vec_var_food)>1
        vec_var_food_mean=[vec_var_food_mean mean(vec_var_food)];
        vec_var_food_se=[vec_var_food_se std(vec_var_food)/sqrt(length(subj_uniq))];
        foodvec_sele=[foodvec_sele foodc];
    end
end
foodvec_cate=categorical(foodvec_sele);
[~,indsort]=sort(vec_var_food_mean);
foodvec_cate=reordercats(foodvec_cate,foodvec_sele(indsort));
h=figure()
bar(foodvec_cate,vec_var_food_mean);
hold on
er=errorbar(foodvec_cate,vec_var_food_mean,2*vec_var_food_se)
er.Color=[0 0 0];
er.LineStyle='none';
saveas(h,[workdir,'bar_qc.fig']);
hold off
% mean food curves
colorlabels={'Rice','Potatoes','Pasta','Grapes','Bread','Berries','Beans','Rice+Fat','Rice+Fiber','Rice+Protein'};
colors=hsv(length(colorlabels)-3);
colors=[colors; repmat(colors(1,:),[3,1])];
linelabels=[repmat({'-'},[1,7]),':','-.','--'];
mask2=ismember(food_and_metigator_vec,colorlabels);
foodsele=food_and_metigator_vec(mask2);
matsele=ymat_shift(:,mask2);
h=figure();
hold on;
hlist={};
i=1;
for foodc=colorlabels
    foodc=foodc{1};
    food_ind=strcmp(foodsele,foodc);
    locmat=matsele(:,food_ind);
    meanvec=mean(locmat,2);
    styleid=strcmp(colorlabels,foodc);
    p=plot(timvecsort,meanvec,'Color',colors(styleid,:),'LineWidth',2,'LineStyle',linelabels(styleid));
    hlist{i}=p(1);
    i=i+1;
end
legend([hlist{:}],colorlabels);
saveas(h,[workdir,'curv.mean.fig']);
close all;
% individual food curve with confidence interval
t=tiledlayout(3,4);
for foodc=colorlabels
    foodc=foodc{1};
    food_ind=strcmp(foodsele,foodc);
    locmat=matsele(:,food_ind);
    meanvec=mean(locmat,2);
    % calculate CI
    ste=std(locmat,0,2)/sqrt(size(locmat,2));
    xconf=[timvecsort; timvecsort(end:-1:1)];
    yconf=[meanvec+ste*2; meanvec(end:-1:1)-ste*2];
    nexttile;
    p=fill(xconf,yconf,'red');
    p.FaceColor=[1 0.8 0.8];
    p.EdgeColor='none';
    hold on;
    plot(timvecsort,meanvec,'LineWidth',2);
    title(foodc);
    ylim([80 170]);
end
saveas(t,[workdir,'curv.allfood.fig']);
close all;