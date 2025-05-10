%% test test bold-csf coupling
clear
 
    load('mats/gs_test_GM_v2.mat','gsl_all','id_test1')
    idx_0=find(nanmean(gsl_all)==0);
    id_test1(idx_0)=[];
    gsl_all(:,idx_0)=[];

    load('mats/csf_v3_auto.mat');
 
  %%
clear pct_tmp_FT
 
gsl=gsl_all;
    for li=1:length(idx_test)%[6:48 50:167]%1:length(list_name)
    close all;
   
    pct_tmp_FT{li}=ft_epi_2D(li,:)'; 
    
    [r_csf_gsl1_FT(li),p_csf_gsl1_FT(li)]=corr(pct_tmp_FT{li},gsl(:,li));   
    
    [cro_r_csf_gsl2_FT(:,li),cro_p_csf_gsl2_FT(:,li)]=xcorr(pct_tmp_FT{li}-nanmean(pct_tmp_FT{li}),gsl(:,li),20,'coeff');  %figure,plot(cro_r_csf_gsl2_FT(:,li))    
 
     % d(gsl)/dt
    tmp1_FT{li}=pct_tmp_FT{li}-nanmean(pct_tmp_FT{li});    
    tmp2_FT{li}=interp1(1:length(pct_tmp_FT{li}),tmp1_FT{li},1.5:1:(length(pct_tmp_FT{li})-0.5),'spline');%figure,plot(tmp1_FT(:,li),'-o'),hold on,plot(1.5:1:139.5,tmp2_FT(:,li),'*')
    [cro_r_csf_gsl3_FT(:,li),cro_p_csf_gsl3_FT(:,li)]=xcorr(tmp2_FT{li},diff(gsl(:,li)),20,'coeff'); % figure,plot(cro_r_csf_gsl3_FT(:,li)),%  figure,plot(diff(gsl_prep4_base(:,li)))  
    [r_csf_gsl3_FT(li),lag(li)]=min(cro_r_csf_gsl3_FT(:,li));
    
    
    li
    
    
    end
 
 
%%
load('mats/id155.mat','id_full4')
 
for lii=1:length(id_full4)
    tmp115_331(lii)=(find(ismember(id_test1,id_full4(lii))));
end
 
cro_r_csf_gsl2_FT115=cro_r_csf_gsl2_FT(:,tmp115_331);
 
 SEM_di_cn_12m = std(cro_r_csf_gsl2_FT115')/sqrt(size(cro_r_csf_gsl2_FT115,2));
 
num=20;
tr=0.607;
%
figure,plot(tr*num:-tr:-tr*num ,mean(cro_r_csf_gsl2_FT115'),'k'),
hold on,
hline(0)
[ph,msg]=jbfill(tr*num:-tr:-tr*num,(mean(cro_r_csf_gsl2_FT115')+SEM_di_cn_12m),(mean(cro_r_csf_gsl2_FT115')-SEM_di_cn_12m),[0.7 .7 .7],[1 1 1],0,0.4)
er.Color = [0 0 0];er.LineStyle = 'none';
yticks([-.4 -.2 0 .2 .4]),
xlim([-15 15])
xticks([-15 -10 -5 0 5 10 15])
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
box off,
 
min(mean(cro_r_csf_gsl2_FT115'))
 
 
 
%% /home/jagust/fenghan/adni/test_coupling_v0.m
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
idx_s1=find(dc4>2& pos4<1);
idx_s2=find(dc4>2& pos4>0);
idx_s3=find(dc4<3& pos4>0);
 
 
idx_pos=find( pos4>0);
idx_neg=find( pos4<1);
 
%% 
cro_r_csf_gsl2_FT115(:,dc4<3& pos4<1)=[];
 
SEM_di_cn_12m = std(cro_r_csf_gsl2_FT115')/sqrt(size(cro_r_csf_gsl2_FT115,2));
 
num=20;
tr=0.607;
%
figure,plot(tr*num:-tr:-tr*num ,mean(cro_r_csf_gsl2_FT115'),'k'),
hold on,
%er=errorbar(tr*num:-tr:-tr*num ,mean(cro_r_csf_gsl2_FT'),SEM_di_cn_12m);
hline(0)
[ph,msg]=jbfill(tr*num:-tr:-tr*num,(mean(cro_r_csf_gsl2_FT115')+SEM_di_cn_12m),(mean(cro_r_csf_gsl2_FT115')-SEM_di_cn_12m),[0.7 .7 .7],[1 1 1],0,0.4)
er.Color = [0 0 0];er.LineStyle = 'none';
yticks([-.4 -.2 0 .2 .4]),
xlim([-15 15])
xticks([-15 -10 -5 0 5 10 15])
set(gca,'fontsize', 30),
set(gca,'linewidth',2.2),
box off,
 
min(mean(cro_r_csf_gsl2_FT115'))
