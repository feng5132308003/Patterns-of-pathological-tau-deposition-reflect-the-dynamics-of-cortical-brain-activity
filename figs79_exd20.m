clear

load('mats/tau_high_cp_ppg_2_dir_exd20.mat','id_full4','tau_cross','couple_update4','num_pg_both_regr','moca_ppg','idx_s*','idx_neg','idx_pos','num_neg_pg','num_pg')
load('mats/tau_high_cp_exd20.mat','id_full4','tau_cross','couple_regr','idx_s*')
%%
% [paths, stats1, stats2] = mediation(double(couple_regr),tau_cross, num_pg_both_regr,  'boottop', 'stats', 'plots')
% [paths, stats1, stats2] = mediation(double(couple_regr(idx_pos)),tau_cross(idx_pos), num_pg_both_regr(idx_pos),  'boottop', 'stats', 'plots')
% [paths, stats1, stats2] = mediation(double(couple_regr(idx_s3)),tau_cross(idx_s3), num_pg_both_regr(idx_s3),  'boottop', 'stats', 'plots')

%%

% [paths, stats1, stats2] = mediation(double(couple_update4),tau_cross, num_pg_both_regr,  'boottop', 'stats', 'plots')
% [paths, stats1, stats2] = mediation(double(couple_update4(idx_pos)),tau_cross(idx_pos), num_pg_both_regr(idx_pos),  'boottop', 'stats', 'plots')
% [paths, stats1, stats2] = mediation(double(couple_update4(idx_s3)),tau_cross(idx_s3), num_pg_both_regr(idx_s3),  'boottop', 'stats', 'plots')

[paths, stats1, stats2] = mediation(double(couple_update4),moca_ppg, num_pg_both_regr,  'boottop', 'stats', 'plots')
[paths, stats1, stats2] = mediation(double(couple_update4(idx_pos)),moca_ppg(idx_pos), num_pg_both_regr(idx_pos),  'boottop', 'stats', 'plots')
[paths, stats1, stats2] = mediation(double(couple_update4(idx_s3)),moca_ppg(idx_s3), num_pg_both_regr(idx_s3),  'boottop', 'stats', 'plots')

%%
[paths, stats1, stats2] = mediation(num_pg_both_regr,tau_cross, double(couple_update4),  'boottop', 'stats', 'plots')
[paths, stats1, stats2] = mediation(num_pg_both_regr(idx_pos), tau_cross(idx_pos),double( couple_update4(idx_pos)), 'boottop', 'stats', 'plots')
[paths, stats1, stats2] = mediation(num_pg_both_regr(idx_s3), tau_cross(idx_s3),double( couple_update4(idx_s3)), 'boottop', 'stats', 'plots')

%% mediation(X, Y, M,...)
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual) KEY
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)