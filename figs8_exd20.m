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
load('mats/label_high_low_83.mat')%_2
tau_cross=nanmean(tau_roi4(:,label_high_num),2);%
%couple_regr4=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tau_cross+nanmean(tau_cross);
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


tau_cross(idx_exd)=[];


[nanmean(age4) std(age4)]
length(gender4)-sum(gender4)
%%
moca=load('mats/moca_sort.mat','moca_update4_sort','id_moca_tmp4_sort','date_moca_tmp4_sort','C2_update_full2')

for li=1:length(id_full4)
    moca_tmp1(li,1)=(find(ismember(moca.C2_update_full2,id_full4(li))));
end
moca_ppg=moca.moca_update4_sort(moca_tmp1);


%%
idx_s1=find(dc4>2& pos4<1);
idx_s2=find(dc4>2& pos4>0);
idx_s3=find(dc4<3& pos4>0);

idx_pos=find( pos4>0);
idx_neg=find( pos4<1);
%%
conf3=[age4 gender4 ];% 
conf4=[ones(size(conf3,1),1),conf3];
tau_cross_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tau_cross+nanmean(tau_cross);


tau_cross_pos=tau_cross_regr(idx_pos);
tau_cross_neg=tau_cross_regr(idx_neg);
tau_cross_s1=tau_cross_regr(idx_s1);
tau_cross_s2=tau_cross_regr(idx_s2);
tau_cross_s3=tau_cross_regr(idx_s3);

moca_ppg_pos=moca_ppg(idx_pos);
moca_ppg_neg=moca_ppg(idx_neg);
moca_ppg_s1=moca_ppg(idx_s1);
moca_ppg_s2=moca_ppg(idx_s2);
moca_ppg_s3=moca_ppg(idx_s3);

%%
color_pos=[255/255 136/255 132/255];
color_neg=[40/255 120/255 181/255];

color_s1=[154/255 201/255 219/255];
color_s2=[248/255 172/255 140/255];
color_s3=[200/255 36/255 35/255];


%% ppg nad tau

figure('position',[0 0 1300 900],'color','w'),

[r,p]=corr(moca_ppg,tau_cross_regr,'rows','pairwise')%,'type','spearman'
ax=subplot(231),
plot(moca_ppg,tau_cross_regr,'.','MarkerSize',20,'MarkerEdgeColor',color_neg),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(moca_ppg_pos,tau_cross_pos,'.','MarkerSize',20,'MarkerEdgeColor',color_pos)
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,

[r,p]=corr(moca_ppg_neg,tau_cross_neg,'rows','pairwise')%,'type','spearman'
ax=subplot(232),
plot(moca_ppg_neg,tau_cross_neg,'.','MarkerSize',20,'MarkerEdgeColor',color_neg),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color_neg;
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,

[r,p]=corr(moca_ppg_pos,tau_cross_pos,'rows','pairwise')%,'type','spearman'
ax=subplot(233),
plot(moca_ppg_pos,tau_cross_pos,'.','MarkerSize',20,'MarkerEdgeColor',color_pos),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color_pos;
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,
%

[r,p]=corr(moca_ppg_s2,tau_cross_s2,'rows','pairwise')%,'type','spearman'
ax=subplot(234),
plot(moca_ppg_s2,tau_cross_s2,'.','MarkerSize',20,'MarkerEdgeColor',color_s2)
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color_s2,
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,

[r,p]=corr(moca_ppg_s3,tau_cross_s3,'rows','pairwise')%,'type','spearman'
ax=subplot(235),
plot(moca_ppg_s3,tau_cross_s3,'.','MarkerSize',20,'MarkerEdgeColor',color_s3)
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color_s3,
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,
