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

tmp1_prin_tmp_subj(idx_exd)=[];% 56 59 377
tmp2_prin_tmp_subj(idx_exd)=[];
locs_all(idx_exd)=[];
num_pos(idx_exd)=[];% sum(num_pos)
num_neg(idx_exd)=[];% sum(num_neg)
%%

liii=0;
sum_pos_1=0;
sum_neg_1=0;


for li=1:length(locs_all)
    liii=liii+1;
    sum_pos_1=sum_pos_1+tmp1_prin_tmp_subj{(liii)};%figure,imagesc(tmp1_prin_tmp_subj{idx_s1(liii)},[-.5 .5])
    
    sum_neg_1=sum_neg_1+tmp2_prin_tmp_subj{(liii)};
    li
end

mean_seg_pos_1=sum_pos_1/sum(num_pos);
mean_seg_neg_1=sum_neg_1/sum(num_neg);


thr=.25;
figure('color','w'),
subplot(121),imagesc(mean_seg_pos_1,[-thr thr]),colorbar,colormap jet,vline(18)
subplot(122),imagesc(mean_seg_neg_1,[-thr thr]),colorbar,colormap jet,vline(18)
%% v0
%%%%%%%%%
pg_surf_test=ft_read_cifti('/home/jagust/fenghan/templates_fh/surface/hcp.embed.dscalar.nii');
pg1 = pg_surf_test.x100307_tfmri_motor_level2_cue_hp200_s2;
pg1 = pg1(1:32492*2,:);

idx_nan_pg1=find(isnan(pg1));
pg1(idx_nan_pg1)=[];
%pg1(isnan(pg1))=0;

pgd1 = discretize(pg1,prctile(pg1,0:100/70:100));
pgd1 = pgd1-min(pgd1)+1;
pg1_uq=unique(pgd1);


for li=1:size(mean_seg_pos_1,2)
    for lii=1:size(mean_seg_pos_1,1)
        mean1(pgd1==pg1_uq(lii),li)=mean_seg_pos_1(lii,li);
        mean2(pgd1==pg1_uq(lii),li)=mean_seg_neg_1(lii,li);
       
    end
end

%%
%brain=ciftiopen('/home/jagust/fenghan/templates_fh/surface/hcp.embed.dscalar.nii','/srv/local/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command',1);
brain=ciftiopen('/home/jagust/fenghan/templates_fh/surface/S900.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii','/srv/local/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command',1);

new_brain=brain;
new_brain.cdata=mean2;%z
ciftisavereset(new_brain,'/home/jagust/fenghan/adni/figs/ppg_correct/test2_exd20.dscalar.nii','/srv/local/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command',1);
unix('wb_command -cifti-smoothing figs/ppg_correct/test2_exd20.dscalar.nii 4 2 COLUMN figs/ppg_correct/test2_exd20_surfSM.dscalar.nii -left-surface /home/jagust/fenghan/templates_fh/surface/100408.L.inflated.32k_fs_LR.surf.gii -right-surface /home/jagust/fenghan/templates_fh/surface/100408.R.inflated.32k_fs_LR.surf.gii')
brain_test=ciftiopen('figs/ppg_correct/test2_exd20_surfSM.dscalar.nii','/srv/local/workbench/bin_rh_linux64/../exe_rh_linux64/wb_command',1);

%% v1
epi1=map1toN(epi_pos_1/sum(num_pos),mask2);
epi2=map1toN(epi_neg_1/sum(num_neg),mask2);
%% v2
pg1=unique(pgd1);

mean1=zeros(length(pgd1),size(epi_pos_1,2));
mean2=zeros(length(pgd1),size(epi_pos_1,2));


for li=1:size(mean_seg_pos_1,2)
    for lii=1:size(mean_seg_pos_1,1)
        mean1(pgd1==pg1(lii),li)=mean_seg_pos_1(lii,li);
        mean2(pgd1==pg1(lii),li)=mean_seg_neg_1(lii,li);
       
    end
end

epi1=map1toN(mean1,mask2);
epi2=map1toN(mean2,mask2);
%%
nii2 = amri_file_loadnii('/home/jagust/fenghan/templates_fh/volume/MNI152_T1_3mm_brain.nii.gz');
anat = nii2.img;
mask = single(anat~=0);


    nii2.hdr.datatype=16;
    nii2.hdr.dim=[4,61,73,61,35,1,1,1];   
    nii2.img=epi1;
    amri_file_savenii(nii2,['figs/ppg_correct/epi1','.nii']);
    
    
nii2 = amri_file_loadnii('/home/jagust/fenghan/templates_fh/volume/MNI152_T1_3mm_brain.nii.gz');
anat = nii2.img;
mask = single(anat~=0);


    nii2.hdr.datatype=16;
    nii2.hdr.dim=[4,61,73,61,35,1,1,1];   
    nii2.img=epi2;
    amri_file_savenii(nii2,['figs/ppg_correct/epi2','.nii']);
%%
unix('wb_command -volume-to-surface-mapping figs/ppg_correct/epi1.nii /home/jagust/fenghan/templates_fh/surface/100408.L.inflated.32k_fs_LR.surf.gii figs/ppg_correct/epi1_L.func.gii -enclosing')
unix('wb_command -volume-to-surface-mapping figs/ppg_correct/epi1.nii /home/jagust/fenghan/templates_fh/surface/100408.R.inflated.32k_fs_LR.surf.gii figs/ppg_correct/epi1_R.func.gii -enclosing')

unix('wb_command -volume-to-surface-mapping figs/ppg_correct/epi2.nii /home/jagust/fenghan/templates_fh/surface/100408.L.inflated.32k_fs_LR.surf.gii figs/ppg_correct/epi2_L.func.gii -enclosing')
unix('wb_command -volume-to-surface-mapping figs/ppg_correct/epi2.nii /home/jagust/fenghan/templates_fh/surface/100408.R.inflated.32k_fs_LR.surf.gii figs/ppg_correct/epi2_R.func.gii -enclosing')

%
unix('wb_command -cifti-create-dense-scalar figs/ppg_correct/epi1.dscalar.nii -left-metric figs/ppg_correct/epi1_L.func.gii -right-metric figs/ppg_correct/epi1_R.func.gii ')
unix('wb_command -cifti-create-dense-scalar figs/ppg_correct/epi2.dscalar.nii -left-metric figs/ppg_correct/epi2_L.func.gii -right-metric figs/ppg_correct/epi2_R.func.gii ')


unix('wb_command -cifti-smoothing figs/ppg_correct/epi1.dscalar.nii 4 2 COLUMN figs/ppg_correct/epi1_surfSM.dscalar.nii -left-surface /home/jagust/fenghan/templates_fh/surface/100408.L.inflated.32k_fs_LR.surf.gii -right-surface /home/jagust/fenghan/templates_fh/surface/100408.R.inflated.32k_fs_LR.surf.gii')
unix('wb_command -cifti-smoothing figs/ppg_correct/epi2.dscalar.nii 4 2 COLUMN figs/ppg_correct/epi2_surfSM.dscalar.nii -left-surface /home/jagust/fenghan/templates_fh/surface/100408.L.inflated.32k_fs_LR.surf.gii -right-surface /home/jagust/fenghan/templates_fh/surface/100408.R.inflated.32k_fs_LR.surf.gii')
%}
