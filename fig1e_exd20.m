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


%%
idx_s1=find(dc4>2& pos4<1);
idx_s2=find(dc4>2& pos4>0);
idx_s3=find(dc4<3& pos4>0);%idx_s4=find(dc3_update<3& pos_update<1);


idx_pos=find( pos4>0);
idx_neg=find( pos4<1);

tau_cross_pos=tau_cross(idx_pos);
tau_cross_neg=tau_cross(idx_neg);

tau_cross_s1=tau_cross(idx_s1);
tau_cross_s2=tau_cross(idx_s2);
tau_cross_s3=tau_cross(idx_s3);


%%

conf3=[age4 gender4];% 
conf4=[ones(size(conf3,1),1),conf3];
couple_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*couple_update4+nanmean(couple_update4);
atp_roi4_regr4=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*atp_roi4+repmat(nanmean(atp_roi4),size(atp_roi4,1),1);


 setenv('OS','Unix')

 %save('mats/tau_high_cp_exd20.mat','id_full4','tau_cross','couple_regr','idx_s*')


%%
color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];
color3=[200/255 36/255 35/255];

figure('position',[0 0 1200 400],'color','w'),

[r,p]=corr(tau_cross,couple_regr,'rows','pairwise')%,'type','spearman'
ax=subplot(131),
plot(tau_cross,couple_regr,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(tau_cross_pos,couple_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
set(gca,'linewidth',3),
box off,


[r,p]=corr(tau_cross_pos,couple_regr(idx_pos),'rows','pairwise')%,'type','spearman'
ax=subplot(132),
plot(tau_cross_pos,couple_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color1;
set(gca,'fontsize', 20),
set(gca,'linewidth',3),ylim([-.7 .1]),yticks([-.9 -.6 -.3 0 .3])
%xlim([-2 20]),xticks([0 5 10 15])
box off,



[r,p]=corr(tau_cross_s3,couple_regr(idx_s3),'rows','pairwise')%,'type','spearman'
ax=subplot(133),
plot(tau_cross_s3,couple_regr(idx_s3),'.','MarkerSize',20,'MarkerEdgeColor',color3),
title([num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color3;
set(gca,'fontsize', 20),
set(gca,'linewidth',3),ylim([-.7 .1]),yticks([-.9 -.6 -.3 0 .3])
%xlim([-2 20]),xticks([0 5 10 15])
box off,
