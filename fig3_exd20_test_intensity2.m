clear

clear
load('mats/tmp_tau_fmri_day183_correct2_subj_replace_structure.mat','atp_roi3','C_tmp2','date_tmp2','*update*','tau_roi_name','atp_roi_name');

clear subj_id_tau_no_pvc
check_regi=readtable('/home/jagust/fenghan/adni/MR_Image_Analysis/check_rest_pp_mni_lowres2.xlsx');

for lii=1:length(C2_update_full4)
    tmp1(lii)=(find(ismember(C2_update_full2,C2_update_full4(lii))));
end


check1=double(string(table2cell(check_regi(tmp1,2:3))));
idx_test_test=find(abs(check1(:,1))>0.5);


dc4=dc3_update4(idx_test_test);
age4=age3_update4(idx_test_test);
gender4=gender3_update4(idx_test_test);
couple_update4=couple_update4(idx_test_test);
C4_update=C_tmp2(idx_test_test);
id_full4=C2_update_full4(idx_test_test);
tau_roi4=tau_roi2_169_update4(idx_test_test,:);
atp_roi4=atp_roi3(idx_test_test,:);
pos4=pos_update4(idx_test_test);

%%
load('/home/jagust/fenghan/adni/mats/ppg_thr010_nogsr_ft010_thr0385_50_correct_v2.mat')


%%
moca=load('mats/moca_sort.mat','moca_update4_sort','id_moca_tmp4_sort','date_moca_tmp4_sort','C2_update_full2')

for li=1:length(id_full4)
    moca_tmp1(li,1)=(find(ismember(moca.C2_update_full2,id_full4(li))));
end
moca_ppg=moca.moca_update4_sort(moca_tmp1);






%%
idx_exd=find(dc4<3& pos4<1);


dc4(idx_exd)=[];
age4(idx_exd)=[];
gender4(idx_exd)=[];
couple_update4(idx_exd)=[];
C4_update(idx_exd)=[];
id_full4(idx_exd)=[];
tau_roi4(idx_exd,:)=[];
atp_roi4(idx_exd,:)=[];
pos4(idx_exd)=[];


moca_ppg(idx_exd)=[];

locs_all(idx_exd)=[];
tmp1_prin_tmp(idx_exd,:)=[];
tmp2_prin_tmp(idx_exd,:)=[];
tmp1_prin_tmp_subj(idx_exd)=[];
tmp2_prin_tmp_subj(idx_exd)=[];

num_pos(idx_exd)=[];
num_neg(idx_exd)=[];

clear idx_pos idx_neg
idx_pos=find( pos4>0);
idx_neg=find( pos4<1);


%%
load('mats/label_high_low_83.mat')%_2
tau_cross=nanmean(tau_roi4(:,label_high_num),2);%
tau_cross=nanmean(tau_roi4(:,label_low_num),2);%

tau_cross_pos=tau_cross(idx_pos);
tau_cross_neg=tau_cross(idx_neg);
[~,p,~,t]=ttest2(tau_cross_pos,tau_cross_neg)

%% subj wise

load('mats/bin_num_high_low.mat','bin_num_low','bin_num_high')



for li=1:length(locs_all)

    low_l_h_int(li,1)=nanmean(nanmean(tmp1_prin_tmp_subj{li}(bin_num_low,:)));
    high_l_h_int(li,1)=nanmean(nanmean(tmp1_prin_tmp_subj{li}(bin_num_high,:)));
    diff_low_high_l_h_int(li,1)=low_l_h_int(li,1)-high_l_h_int(li,1);

    low_h_l_int(li,1)=nanmean(nanmean(tmp2_prin_tmp_subj{li}(bin_num_low,:)));
    high_h_l_int(li,1)=nanmean(nanmean(tmp2_prin_tmp_subj{li}(bin_num_high,:)));
    diff_low_high_h_l_int(li,1)=low_h_l_int(li,1)-high_h_l_int(li,1);
    
    li
end
figure,imagesc(tmp1_prin_tmp_subj{li}/num_pos(li),[-.25 .25])

%%
clear *int
for li=1:length(locs_all)
    if num_pos(li)==0
        diff_low_high_l_h_int(li,1)=0;
    else
        clear tmp1 tmp_med1
        tmp1=tmp1_prin_tmp_subj{li}/num_pos(li);
        tmp_med1=prctile(tmp1(:),0);
        tmp1(tmp1<tmp_med1)=nan;
        
        low_l_h_int(li,1)=nanmean(nanmean(tmp1(bin_num_low,:)));
        high_l_h_int(li,1)=nanmean(nanmean(tmp1(bin_num_high,:)));
        diff_low_high_l_h_int(li,1)=high_l_h_int(li,1)-low_l_h_int(li,1);
    end
    
    if num_neg(li)==0
        diff_low_high_h_l_int(li,1)=0;
    else
        clear tmp2 tmp_med2
        tmp2=tmp2_prin_tmp_subj{li}/num_neg(li);
        tmp_med2=prctile(tmp2(:),0);
        tmp2(tmp2<tmp_med2)=nan;

        low_h_l_int(li,1)=nanmean(nanmean(tmp2(bin_num_low,:)));
        high_h_l_int(li,1)=nanmean(nanmean(tmp2(bin_num_high,:)));
        diff_low_high_h_l_int(li,1)=high_h_l_int(li,1)-low_h_l_int(li,1);
    end
    li
end
figure,imagesc(tmp1,[-.25 .25])
figure,imagesc(tmp2,[-.25 .25])

%%
%figure,imagesc(tmp1_prin_tmp{1,6},[-.25 .25])
%nanmean(nanmean(tmp1_prin_tmp{li,lm}(1:num1,:)))
%nanmean(nanmean(tmp1_prin_tmp{li,lm}(num2:70,:)))
%% segment wise
%{
clear *_int

for li=1:length(locs_all)

    if num_pos(li)==0
        diff_low_high_l_h_int(li,1)=0;
    else
     for lm=1:length(num_pos(li))
        clear tmp1 tmp2 
        tmp1(lm,1)=nanmean(nanmean(tmp1_prin_tmp{li,lm}(bin_num_low,:)));
        low_l_h_int(li,1)=nanmean(tmp1);
        tmp2(lm,1)=nanmean(nanmean(tmp1_prin_tmp{li,lm}(bin_num_high,:)));
        high_l_h_int(li,1)=nanmean(tmp2);
        diff_low_high_l_h_int(li,1)=low_l_h_int(li,1)-high_l_h_int(li,1);
     end
    end

    if num_neg(li)==0
        diff_low_high_h_l_int(li,1)=0;
    else
    for lm=1:length(num_neg(li))
        clear tmp3 tmp4
        tmp3(lm,1)=nanmean(nanmean(tmp2_prin_tmp{li,lm}(bin_num_low,:)));
        low_h_l_int(li,1)=nanmean(tmp3);
        tmp4(lm,1)=nanmean(nanmean(tmp2_prin_tmp{li,lm}(bin_num_high,:)));
        high_h_l_int(li,1)=nanmean(tmp4);

        diff_low_high_h_l_int(li,1)=low_h_l_int(li,1)-high_h_l_int(li,1);
    end
    end
    li
end
%}
%%
clear tmp_intensity*
tmp_intensity=(diff_low_high_l_h_int+diff_low_high_h_l_int)/2;%
% tmp_intensity=high_h_l_int+high_l_h_int%;
 %tmp_intensity=low_h_l_int+low_l_h_int;%
 

%%
% diff_low_high_l_h_int(abs(diff_low_high_l_h_int)<0.001)=nan;
% diff_low_high_h_l_int(abs(diff_low_high_h_l_int)<0.001)=nan;

tmp_intensity(isnan(tmp_intensity))=0;
idx_nan=find(abs(tmp_intensity)<0.001);


age4(idx_nan)=[];
gender4(idx_nan)=[];

tau_cross(idx_nan)=[];
couple_update4(idx_nan)=[];
moca_ppg(idx_nan)=[];

tmp_intensity(idx_nan)=[];

pos4(idx_nan)=[];
idx_pos=find( pos4>0);
idx_neg=find( pos4<1);

%%
conf3=[age4 gender4 ];% 
conf4=[ones(size(conf3,1),1),conf3];

% diff_low_high_l_h_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*diff_low_high_l_h_int+nanmean(diff_low_high_l_h_int);
% diff_low_high_h_l_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*diff_low_high_h_l_int+nanmean(diff_low_high_h_l_int);
% 
% high_l_h_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*high_l_h_int+nanmean(high_l_h_int);
% high_h_l_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*high_h_l_int+nanmean(high_h_l_int);
% 
% low_l_h_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*low_l_h_int+nanmean(low_l_h_int);
% low_h_l_int_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*low_h_l_int+nanmean(low_h_l_int);
 tmp_intensity2=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tmp_intensity+nanmean(tmp_intensity);
%%
 tmp_int_h=high_h_l_int+high_l_h_int%;
 tmp_int_l=low_h_l_int+low_l_h_int;%

 tmp_int_h=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tmp_int_h+nanmean(tmp_int_h);
 tmp_int_l=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tmp_int_l+nanmean(tmp_int_l);

 SEM_tmp1 = std(tmp_int_h)/sqrt(length(tmp_int_h));
 SEM_tmp2 = std(tmp_int_l)/sqrt(length(tmp_int_l));

[~,p1,~,t1]=ttest2(tmp_int_h,tmp_int_l)
[~,p1,~,t1]=ttest(tmp_int_h,tmp_int_l)

color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];

figure,
h1=bar([1],[ mean(tmp_int_h)])
hold on,
h2=bar([ 2 ],[mean(tmp_int_l)])

er=errorbar([1 2] ,[mean(tmp_int_h) mean(tmp_int_l)]...
    ,[SEM_tmp1 SEM_tmp2]),
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.MarkerSize=13;
er.LineWidth = 2.2; 
title([num2str(p1)])
set(h1, 'FaceColor', color1, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h2, 'FaceColor', color2, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
xticks([1 2 ])
%ylim([-.45 -.1])
%%

% tmp_intensity=(diff_low_high_l_h_int_regr+diff_low_high_h_l_int_regr)/2;
% tmp_intensity=high_h_l_int_regr+high_l_h_int_regr;%
% tmp_intensity=low_h_l_int_regr+low_l_h_int_regr;%
%%
color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];
color3=[200/255 36/255 35/255];

figure('position',[0 0 1200 800],'color','w'),

[r,p]=corr(tau_cross,tmp_intensity2,'rows','pairwise')%,'type','spearman'
ax=subplot(131),
plot(tau_cross,tmp_intensity2,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(tau_cross(idx_pos),tmp_intensity2(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),%ylim([0 25])
set(gca,'linewidth',3),%ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%

[r,p]=corr(couple_update4,tmp_intensity2,'rows','pairwise')%,'type','spearman'
ax=subplot(132),
plot(couple_update4,tmp_intensity2,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(couple_update4(idx_pos),tmp_intensity2(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%

[r,p]=corr(moca_ppg,tmp_intensity2,'rows','pairwise')%,'type','spearman'
ax=subplot(133),
plot(moca_ppg,tmp_intensity2,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(moca_ppg(idx_pos),tmp_intensity2(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%%
load('mats/tau_high_cp_ppg_2_dir_exd20.mat','num_pg_both_regr')
num_pg_both_regr(idx_nan)=[];

figure('position',[0 0 1200 800],'color','w'),

[r,p]=corr(tau_cross,tmp_intensity2,'rows','pairwise')%,'type','spearman'
ax=subplot(121),
plot(tau_cross,tmp_intensity2,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(tau_cross(idx_pos),tmp_intensity2(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),%ylim([0 25])
set(gca,'linewidth',3),%ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%

[r,p]=corr(num_pg_both_regr,tmp_intensity2,'rows','pairwise')%,'type','spearman'
ax=subplot(122),
plot(num_pg_both_regr,tmp_intensity2,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(num_pg_both_regr(idx_pos),tmp_intensity2(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%
