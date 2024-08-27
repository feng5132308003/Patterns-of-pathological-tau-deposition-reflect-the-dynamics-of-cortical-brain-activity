
clear

load('mats/tau_high_cp_ppg_2_dir_exd20.mat','id_full4','tau_cross','couple_update4','num_pg_both_regr','moca_ppg','idx_s*','idx_neg','idx_pos','num_neg_pg','num_pg')
load('mats/fig2_exd20.mat','idx_s*')

load('mats/age_sex_95.mat','age4','gender4')

conf3=[age4 gender4 ];% atp_roi4
conf4=[ones(size(conf3,1),1),conf3];

%% figure 6A
%%
num_pg_both_regr_P=num_pg_both_regr(idx_pos);
num_pg_both_regr_N=num_pg_both_regr(idx_neg);


SEM_num_pg_both_regr_P = std(num_pg_both_regr_P)/sqrt(length(num_pg_both_regr_P));
SEM_num_pg_both_regr_N = std(num_pg_both_regr_N)/sqrt(length(num_pg_both_regr_N));

[~,p1,~,t1]=ttest2(num_pg_both_regr_P,num_pg_both_regr_N);

color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];
color3=[200/255 36/255 35/255];

figure,
h1=bar([1],[ mean(num_pg_both_regr_N)])
hold on,
h2=bar([ 2 ],[mean(num_pg_both_regr_P)])

er=errorbar([1 2] ,[mean(num_pg_both_regr_N) mean(num_pg_both_regr_P)]...
    ,[SEM_num_pg_both_regr_N SEM_num_pg_both_regr_P]),
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
color4=[248/255 172/255 140/255];
color5=[200/255 36/255 35/255];

num_pg_both_regr_S1=num_pg_both_regr(idx_s1);
num_pg_both_regr_S2=num_pg_both_regr(idx_s2);
num_pg_both_regr_S3=num_pg_both_regr(idx_s3);


SEM_num_pg_both_regr_S1 = std(num_pg_both_regr_S1)/sqrt(length(num_pg_both_regr_S1));
SEM_num_pg_both_regr_S2 = std(num_pg_both_regr_S2)/sqrt(length(num_pg_both_regr_S2));
SEM_num_pg_both_regr_S3 = std(num_pg_both_regr_S3)/sqrt(length(num_pg_both_regr_S3));

[~,p1,~,t1]=ttest2(num_pg_both_regr_S1,num_pg_both_regr_S2);
[~,p2,~,t2]=ttest2(num_pg_both_regr_S1,num_pg_both_regr_S3);
[~,p3,~,t3]=ttest2(num_pg_both_regr_S3,num_pg_both_regr_S2);
% ordinary regression
clear test_123
test_123= [ones(length(idx_s1),1); 2*ones(length(idx_s2),1);3*ones(length(idx_s3),1)];

[B,dev,stats] = mnrfit([num_pg_both_regr(idx_s1,1);num_pg_both_regr(idx_s2,1);num_pg_both_regr(idx_s3,1)],test_123,'model','ordinal');
 stats.p(end)
    
figure,
h1=bar([1],[mean(num_pg_both_regr_S1) ])
hold on,
h2=bar([2],[mean(num_pg_both_regr_S2) ])
h3=bar([3],[mean(num_pg_both_regr_S3)])
er=errorbar([1 2 3] ,[mean(num_pg_both_regr_S1) mean(num_pg_both_regr_S2) mean(num_pg_both_regr_S3)]...
    ,[SEM_num_pg_both_regr_S1 SEM_num_pg_both_regr_S2 SEM_num_pg_both_regr_S3]),
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
xticks([1 2 3]),%ylim([0 2]),yticks([0 .5 1 1.5 2])
box off



set(gca,'linewidth',3),
%xlim([-2 20]),xticks([0 5 10 15])
box off,
%% figure 6B
%%
color1=[255/255 136/255 132/255];
color2=[40/255 120/255 181/255];
color3=[200/255 36/255 35/255];

num_pg_pos_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*num_pg'+nanmean(num_pg);
num_pg_neg_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*num_neg_pg'+nanmean(num_neg_pg);


num_pg_pos_regr_P=num_pg_pos_regr(idx_pos);
num_pg_pos_regr_N=num_pg_pos_regr(idx_neg);

num_pg_neg_regr_P=num_pg_neg_regr(idx_pos);
num_pg_neg_regr_N=num_pg_neg_regr(idx_neg);



SEM_num_pg_pos_regr_P = std(num_pg_pos_regr_P)/sqrt(length(num_pg_pos_regr_P));
SEM_num_pg_pos_regr_N = std(num_pg_pos_regr_N)/sqrt(length(num_pg_pos_regr_N));

SEM_num_pg_neg_regr_P = std(num_pg_neg_regr_P)/sqrt(length(num_pg_neg_regr_P));
SEM_num_pg_neg_regr_N = std(num_pg_neg_regr_N)/sqrt(length(num_pg_neg_regr_N));

[~,p1,~,t1]=ttest2(num_pg_pos_regr_P,num_pg_pos_regr_N)
[~,p2,~,t2]=ttest2(num_pg_neg_regr_P,num_pg_neg_regr_N)

figure,
h1=bar([1],[mean(num_pg_pos_regr_N) ])
hold on,
h2=bar([4],[ mean(num_pg_neg_regr_N) ])
h3=bar([2],[mean(num_pg_pos_regr_P)  ])
h4=bar([5],[mean(num_pg_neg_regr_P) ])

er=errorbar([1 4 2 5] ,[mean(num_pg_pos_regr_N) mean(num_pg_neg_regr_N) mean(num_pg_pos_regr_P) mean(num_pg_neg_regr_P)]...
    ,[SEM_num_pg_pos_regr_N SEM_num_pg_neg_regr_N SEM_num_pg_pos_regr_P SEM_num_pg_neg_regr_P]),
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.MarkerSize=13;
er.LineWidth = 2.2; 
title([num2str(p1)])
set(h1, 'FaceColor', color2, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h2, 'FaceColor', color2, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(h3, 'FaceColor', color1, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h4, 'FaceColor', color1, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
xticks([1 2 4 5])
box off,%ylim([-.45 -.1])

%%
color3=[40/255 120/255 181/255];
color4=[248/255 172/255 140/255];
color5=[200/255 36/255 35/255];

num_pg_pos_regr_S1=num_pg_pos_regr(idx_s1);
num_pg_pos_regr_S2=num_pg_pos_regr(idx_s2);
num_pg_pos_regr_S3=num_pg_pos_regr(idx_s3);

num_pg_neg_regr_S1=num_pg_neg_regr(idx_s1);
num_pg_neg_regr_S2=num_pg_neg_regr(idx_s2);
num_pg_neg_regr_S3=num_pg_neg_regr(idx_s3);

SEM_num_pg_pos_regr_S1 = std(num_pg_pos_regr_S1)/sqrt(length(num_pg_pos_regr_S1));
SEM_num_pg_pos_regr_S2 = std(num_pg_pos_regr_S2)/sqrt(length(num_pg_pos_regr_S2));
SEM_num_pg_pos_regr_S3 = std(num_pg_pos_regr_S3)/sqrt(length(num_pg_pos_regr_S3));

SEM_num_pg_neg_regr_S1 = std(num_pg_neg_regr_S1)/sqrt(length(num_pg_neg_regr_S1));
SEM_num_pg_neg_regr_S2 = std(num_pg_neg_regr_S2)/sqrt(length(num_pg_neg_regr_S2));
SEM_num_pg_neg_regr_S3 = std(num_pg_neg_regr_S3)/sqrt(length(num_pg_neg_regr_S3));


[~,p1,~,t1]=ttest2(num_pg_pos_regr_S1,num_pg_pos_regr_S2)
[~,p2,~,t2]=ttest2(num_pg_pos_regr_S1,num_pg_pos_regr_S3)
[~,p3,~,t3]=ttest2(num_pg_pos_regr_S3,num_pg_pos_regr_S2)

[~,p4,~,t4]=ttest2(num_pg_neg_regr_S1,num_pg_neg_regr_S2)
[~,p5,~,t5]=ttest2(num_pg_neg_regr_S1,num_pg_neg_regr_S3)
[~,p6,~,t6]=ttest2(num_pg_neg_regr_S3,num_pg_neg_regr_S2)

% ordinary regression
clear test_123
test_123= [ones(length(idx_s1),1); 2*ones(length(idx_s2),1);3*ones(length(idx_s3),1)];

[B,dev,stats1] = mnrfit([num_pg_pos_regr(idx_s1,1);num_pg_pos_regr(idx_s2,1);num_pg_pos_regr(idx_s3,1)],test_123,'model','ordinal');
[B,dev,stats2] = mnrfit([num_pg_neg_regr(idx_s1,1);num_pg_neg_regr(idx_s2,1);num_pg_neg_regr(idx_s3,1)],test_123,'model','ordinal');
[stats1.p(end) stats2.p(end)]
 
 
 
figure,
h1=bar([1],[mean(num_pg_pos_regr_S1) ])
hold on,
h2=bar([2],[mean(num_pg_pos_regr_S2) ])
h3=bar([3],[mean(num_pg_pos_regr_S3)])

h4=bar([5],[mean(num_pg_neg_regr_S1) ])
h5=bar([6],[mean(num_pg_neg_regr_S2) ])
h6=bar([7],[mean(num_pg_neg_regr_S3)])


er1=errorbar([1 2 3] ,[mean(num_pg_pos_regr_S1) mean(num_pg_pos_regr_S2) mean(num_pg_pos_regr_S3)]...
    ,[SEM_num_pg_pos_regr_S1 SEM_num_pg_pos_regr_S2 SEM_num_pg_pos_regr_S3]),
er2=errorbar([ 5 6 7] ,[mean(num_pg_neg_regr_S1) mean(num_pg_neg_regr_S2) mean(num_pg_neg_regr_S3)]...
    ,[SEM_num_pg_neg_regr_S1 SEM_num_pg_neg_regr_S2 SEM_num_pg_neg_regr_S3]),
er1.Color = [0 0 0];
er1.LineStyle = 'none'; 
er1.MarkerSize=13;
er1.LineWidth = 2.2; 

er2.Color = [0 0 0];
er2.LineStyle = 'none'; 
er2.MarkerSize=13;
er2.LineWidth = 2.2; 
set(h1, 'FaceColor', color3, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h2, 'FaceColor', color4, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(h3, 'FaceColor', color5, 'EdgeColor',  [0 0 0],'LineWidth',2.2),

set(h4, 'FaceColor', color3, 'EdgeColor', [0 0 0],'LineWidth',2.2),
set(h5, 'FaceColor', color4, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(h6, 'FaceColor', color5, 'EdgeColor',  [0 0 0],'LineWidth',2.2),
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
title([num2str(p1),' ',num2str(p2),' ',num2str(p3)])
xticks([1 2 3 5 6 7]),%ylim([0 2]),yticks([0 .5 1 1.5 2])
box off


%% figure 6e
couple_update4_regr=(diag(ones(size(conf4,1),1))-conf4*inv(conf4'*conf4)*(conf4)')*couple_update4+nanmean(couple_update4);



figure('position',[0 0 1600 900],'color','w'),

[r,p]=corr(couple_update4_regr,num_pg','rows','pairwise')%,'type','spearman'
ax=subplot(241),
plot(num_pg,couple_update4_regr,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([': r=',num2str(r),' p=',num2str(p)])%tau_roi_name{b56_roi}([1:7]),
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(num_pg(idx_pos),couple_update4_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,


[r,p]=corr(couple_update4_regr(idx_pos),num_pg(idx_pos)','rows','pairwise')%,'type','spearman'
ax=subplot(242),
plot(num_pg(idx_pos),couple_update4_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color1;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,


[r,p]=corr(couple_update4_regr(idx_neg),num_pg(idx_neg)','rows','pairwise')%,'type','spearman'
ax=subplot(243),
plot(num_pg(idx_neg),couple_update4_regr(idx_neg),'.','MarkerSize',20,'MarkerEdgeColor',color2),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color2;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,



%num_neg_pg(num_neg_pg>12)=nan;
[r,p]=corr(couple_update4_regr,num_neg_pg','rows','pairwise')%,'type','spearman'
ax=subplot(245),
plot(num_neg_pg,couple_update4_regr,'.','MarkerSize',20,'MarkerEdgeColor',color2),
title([': r=',num2str(r),' p=',num2str(p)])%tau_roi_name{b56_roi}([1:7]),
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = 'k';
hold on,
plot(num_neg_pg(idx_pos),couple_update4_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1)
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,
%save('mats/tmp_data_fig_3d_exd20.mat','num_neg_pg','couple_update4_regr','idx_pos')


[r,p]=corr(couple_update4_regr(idx_pos),num_neg_pg(idx_pos)','rows','pairwise')%,'type','spearman'
ax=subplot(246),
plot(num_neg_pg(idx_pos),couple_update4_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color1;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,


[r,p]=corr(couple_update4_regr(idx_neg),num_neg_pg(idx_neg)','rows','pairwise')%,'type','spearman'
ax=subplot(247),
plot(num_neg_pg(idx_neg),couple_update4_regr(idx_neg),'.','MarkerSize',20,'MarkerEdgeColor',color2),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color2;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,

[r,p]=corr(couple_update4_regr(idx_s3),num_neg_pg(idx_s3)','rows','pairwise')%,'type','spearman'
ax=subplot(248),
plot(num_neg_pg(idx_s3),couple_update4_regr(idx_s3),'.','MarkerSize',20,'MarkerEdgeColor',color5),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color5;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,
%%
num_neg_pg(num_neg_pg>12)=nan;


figure,
[r,p]=corr(couple_update4_regr(idx_pos),num_neg_pg(idx_pos)','rows','pairwise')%,'type','spearman'
ax=subplot(121),
plot(num_neg_pg(idx_pos),couple_update4_regr(idx_pos),'.','MarkerSize',20,'MarkerEdgeColor',color1),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color1;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,



[r,p]=corr(couple_update4_regr(idx_s3),num_neg_pg(idx_s3)','rows','pairwise')%,'type','spearman'
ax=subplot(122),
plot(num_neg_pg(idx_s3),couple_update4_regr(idx_s3),'.','MarkerSize',20,'MarkerEdgeColor',color5),
title(['r=',num2str(r),' p=',num2str(p)])
h2 = lsline(ax);
h2.LineWidth = 3;
h2.Color = color5;
set(gca,'fontsize', 20),xlim([-1 16])
set(gca,'linewidth',3),ylim([-1 .2]),yticks([-.9 -.6 -.3 0 .3])
box off,
