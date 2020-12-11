%% differnt GCM and cGCM methods

%% initialize model and configuration parameters
addpath ../lib
% A = var9_test;
A = var5_test;% var5_test and var9_test are from the mvgc toolbox by Lionel Barnett and Anil K. Seth

n = size(A,1);
SIG = eye(n) + ones(n);
N = 2000;
freq = [0:N-1]/N*2*pi;
%% compute GCM, fGCM_Ent and fGCM_Geweke
gcm = var_to_gcm(A,SIG);%check
fgcm_ent = var_to_fgcm_ent(A,SIG,freq);
fgcm_Geweke = var_to_fgcm_Geweke(A,SIG,freq);
% the mean of fgcm is equal to gcm
mean(fgcm_ent,3)-gcm
mean(fgcm_Geweke,3)-gcm
% 
%% compute cGCM_Std, fcGCM_Std_Ent and fcGCM_Std_Geweke
cgcm_std = var_to_cgcm_std(A,SIG);
fcgcm_std_ent = var_to_fcgcm_std_ent(A,SIG,freq);
fcgcm_std_Geweke = var_to_fcgcm_std_Geweke(A,SIG,freq);
mean(fcgcm_std_ent,3)-cgcm_std
mean(fcgcm_std_Geweke,3)-cgcm_std
%% compute cGCM_SEnt and fcGCM_SEnt
cgcm_sent = var_to_cgcm_sent(A,SIG);
fcgcm_sent = var_to_fcgcm_sent(A,SIG,freq);
mean(fcgcm_sent,3)-cgcm_sent
%% compute cGCM_JEnt and fcGCM_JEnt
cgcm_jent = var_to_cgcm_jent(A,SIG);
fcgcm_jent = var_to_fcgcm_jent(A,SIG,freq);
mean(fcgcm_jent,3)-cgcm_jent
%% compare the three cgcm 
cgcm_jent - cgcm_std
cgcm_std - cgcm_sent



