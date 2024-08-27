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


%%
idx_s1=find(dc4>2& pos4<1);
idx_s2=find(dc4>2& pos4>0);
idx_s3=find(dc4<3& pos4>0);%idx_s4=find(dc3_update<3& pos_update<1);


idx_pos=find( pos4>0);
idx_neg=find( pos4<1);

%%


conf3=[age4 gender4];% 
conf4=[ones(size(conf3,1),1),conf3];
couple_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*couple_update4+nanmean(couple_update4);
atp_roi4_regr4=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*atp_roi4+repmat(nanmean(atp_roi4),size(atp_roi4,1),1);
tau_roi4_regr4=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tau_roi4+repmat(nanmean(tau_roi4),size(tau_roi4,1),1);
%%
 setenv('OS','Unix')
 %% tau_roi4_regr4=atp_roi4_regr4;
 %%
[~,p_t1,~,tstat1]=ttest2(tau_roi4_regr4(idx_pos,:),tau_roi4_regr4(idx_neg,:))

figure,plot(-log(p_t1),'-o'),hline(3)

[h1, crit_p1, adj_p1]=fdr_bh(p_t1,0.05,'pdep','yes')% no

roi1=tau_roi_name(adj_p1<0.05);%adj_p1
%%
r=tstat1.tstat.*(adj_p1<0.001);%adj_p1

mean_diff_1=mean(tau_roi4_regr4(idx_pos,:))-mean(tau_roi4_regr4(idx_neg,:));
%mean_diff_1=mean(tau_roi4(idx_pos,:))-mean(tau_roi4(idx_neg,:));

r=mean_diff_1.*(p_t1<0.05);%adj_p1%

%%
r=nanmean(tau_roi4_regr4);

%%

pg_surf_test=ft_read_cifti('/home/jagust/fenghan/templates_fh/surface/hcp.embed.dscalar.nii');
pg1 = pg_surf_test.x100307_tfmri_motor_level2_cue_hp200_s2;
pg1 = pg1(1:32492*2,:);

roi_fs_test=ft_read_cifti('/home/jagust/fenghan/templates_fh/surface/100307.aparc.32k_fs_LR.dlabel.nii');
test_roi_fs = roi_fs_test.x100307_aparc;
tmp2=unique(test_roi_fs);

parcel_epimsk68=zeros(length(test_roi_fs),68);

for li=1:68
    parcel_epimsk68(test_roi_fs==tmp2(li+1),li)=1/sum(test_roi_fs==tmp2(li+1));
   li 
end
parcel_epimsk68(isnan(pg1),:)=[];
pg1(isnan(pg1),:)=[];
pg_fs=parcel_epimsk68'*pg1;

%%
idx_tmp_l=[32 33 11 34:38 40:41 12 42 12 43 45 44 13 13 13 46:51 11 52:55 12 56:57 39]-1;
idx_tmp_r=[58 59 21 60:64 66:67 22 68 22 69 71 70 23 23 23 72:77 21 78:81 22 82:83 65]-1;
mean_diff_fs=mean_diff_1([idx_tmp_l idx_tmp_r]);
%mean_diff_fs=r([idx_tmp_l idx_tmp_r]);

% 
% mean_diff_fs([31 65])=[];
% pg_fs([31 65])=[];

[r,p]=corr(pg_fs,mean_diff_fs','rows','pairwise')%,'type','spearman'
figure,
ax=subplot(121),
plot(mean_diff_fs,pg_fs,'.','MarkerSize',20,'MarkerEdgeColor',[1/255 1/255 1/255]),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
set(gca,'fontsize', 30),
set(gca,'linewidth',3),ylim([-7 7]),yticks([-6  0 6]),xticks([  0 .4 .8])
box off,

mean_fs=nanmean(tau_roi4_regr4(:,[idx_tmp_l idx_tmp_r]));
[r,p]=corr(pg_fs,mean_fs','rows','pairwise')%,'type','spearman'

ax=subplot(122),
plot(mean_fs,pg_fs,'.','MarkerSize',20,'MarkerEdgeColor',[1/255 1/255 1/255]),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
set(gca,'fontsize', 30),
set(gca,'linewidth',3),ylim([-7 7]),yticks([-6  0 6])
box off,
%% thk and tau, mean
mean_atp=nanmean(atp_roi4(:,[1:34 36:69]))';
[r,p]=corr(mean_atp,mean_diff_fs','rows','pairwise')%,'type','spearman'
figure,
ax=subplot(211),
plot(mean_atp,mean_diff_fs,'.','MarkerSize',20,'MarkerEdgeColor',[1/255 1/255 1/255]),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = [142/255 207/255 201/255];
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,
%% thk and tau, mean diff
% mean_atp_diff_1=mean(atp_roi4(idx_pos,[1:34 36:69]))-mean(atp_roi4(idx_neg,[1:34 36:69]));
% mean_atp_diff_1=mean(atp_roi4(idx_pos,:))-mean(atp_roi4(idx_neg,:));
% r=mean_atp_diff_1;
% save('mats/thk_mean_diff_nofdr_pos_neg.mat','r')%_t_

[r,p]=corr(mean_atp_diff_1',mean_diff_fs','rows','pairwise')%,'type','spearman'
figure,
ax=subplot(211),
plot(mean_atp_diff_1,mean_diff_fs,'.','MarkerSize',20,'MarkerEdgeColor',[1/255 1/255 1/255]),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = [142/255 207/255 201/255];
set(gca,'fontsize', 20),
set(gca,'linewidth',3),
box off,


%% figure 4 d

%%
load('mats/label_high_low_83.mat')
tau_cross=sum(tau_roi4(:,label_high_num),2)/(length(label_high_num));%-sum(tau_roi4(:,label_low_num),2) +length(label_low_num)

tau_cross_pos=tau_cross(idx_pos);
tau_cross_neg=tau_cross(idx_neg);

tau_cross_s1=tau_cross(idx_s1);
tau_cross_s2=tau_cross(idx_s2);
tau_cross_s3=tau_cross(idx_s3);

couple_regr4=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*tau_cross+nanmean(tau_cross);


%%
gs_amp_pos_tmp1=couple_regr4(idx_pos);
gs_amp_pos_tmp2=couple_regr4(idx_neg);


SEM_pg_pos_tmp1 = std(gs_amp_pos_tmp1)/sqrt(length(gs_amp_pos_tmp1));
SEM_pg_pos_tmp2 = std(gs_amp_pos_tmp2)/sqrt(length(gs_amp_pos_tmp2));

[~,p1,~,t1]=ttest2(gs_amp_pos_tmp1,gs_amp_pos_tmp2);

color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];
color3=[200/255 36/255 35/255];

figure,
h1=bar([1],[ mean(gs_amp_pos_tmp2)])
hold on,
h2=bar([ 2 ],[mean(gs_amp_pos_tmp1)])

er=errorbar([2  1] ,[mean(gs_amp_pos_tmp1) mean(gs_amp_pos_tmp2)]...
    ,[SEM_pg_pos_tmp1 SEM_pg_pos_tmp2]),
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.MarkerSize=13;
er.LineWidth = 2.2; 
title([num2str(p1)])
set(h1, 'FaceColor', color2, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h2, 'FaceColor', color1, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
xticks([1 2 ])
box off,%ylim([-.45 -.1])

%%
color3=[40/255 120/255 181/255];
%color3=[154/255 201/255 219/255];
color4=[248/255 172/255 140/255];
color5=[200/255 36/255 35/255];

gs_amp_pos_tmp1=couple_regr4(idx_s1);
gs_amp_pos_tmp2=couple_regr4(idx_s2);
gs_amp_pos_tmp3=couple_regr4(idx_s3);


SEM_pg_pos_tmp1 = std(gs_amp_pos_tmp1)/sqrt(length(gs_amp_pos_tmp1));
SEM_pg_pos_tmp2 = std(gs_amp_pos_tmp2)/sqrt(length(gs_amp_pos_tmp2));
SEM_pg_pos_tmp3 = std(gs_amp_pos_tmp3)/sqrt(length(gs_amp_pos_tmp3));

[~,p1,~,t1]=ttest2(gs_amp_pos_tmp1,gs_amp_pos_tmp2);
[~,p2,~,t2]=ttest2(gs_amp_pos_tmp1,gs_amp_pos_tmp3);
[~,p3,~,t3]=ttest2(gs_amp_pos_tmp3,gs_amp_pos_tmp2);
% ordinary regression
clear test_123
test_123= [ones(length(idx_s1),1); 2*ones(length(idx_s2),1);3*ones(length(idx_s3),1)];

[B,dev,stats] = mnrfit([couple_regr4(idx_s1,1);couple_regr4(idx_s2,1);couple_regr4(idx_s3,1)],test_123,'model','ordinal');
 stats.p(end)
    
figure,
h1=bar([1],[mean(gs_amp_pos_tmp1) ])
hold on,
h2=bar([2],[mean(gs_amp_pos_tmp2) ])
h3=bar([3],[mean(gs_amp_pos_tmp3)])
er=errorbar([1 2 3] ,[mean(gs_amp_pos_tmp1) mean(gs_amp_pos_tmp2) mean(gs_amp_pos_tmp3)]...
    ,[SEM_pg_pos_tmp1 SEM_pg_pos_tmp2 SEM_pg_pos_tmp3]),
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.MarkerSize=13;
er.LineWidth = 2.2; 
set(h1, 'FaceColor', color3, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h2, 'FaceColor', color4, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(h3, 'FaceColor', color5, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
title([num2str(p1),' ',num2str(p2),' ',num2str(p3)])
xticks([1 2 3]),ylim([0 2]),yticks([0 .5 1 1.5 2])
box off