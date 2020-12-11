function F = var_to_cgcm_jent(A,Sig,U)
%% F = var_to_cgcm_ient(A,Sig,U)
%
% Calculate conditional time-domain cGCM based on VAR parameters
% using joint entropy minimization 
%
%% Arguments
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%
% _output_
%
%     F          Granger causality matrix
%
%
% (C) Lipeng Ning, 2020. 
%% 


[n,~,L] = size(A);
N = n*L;
F = nan(n);


% full state-space model
Sys = var_to_ss(A);
iSys = sys_to_iss(Sys); 

%%
if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);
end

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

F_IJ = nan(1,L);
F_JI = nan(1,L);


for l = 1:L % change for to parfor for parallel computation

    i = I(l);
    j = J(l);
    
    % the (i j) block of iSys
    iSys_ij = [];
    iSys_ij.A = iSys.A;
    iSys_ij.B = iSys.B(:,[i j]);
    iSys_ij.C = iSys.C([i,j],:);
    iSys_ij.D = iSys.D([i,j],[i,j]);
    % the inverse of iSys_ij
    iiSys_ij = sys_to_iss(iSys_ij);% inverse of G_ij
 %    
    % the PSD of the new (i j) is iiSys_ij*Sig_ij*iiSys_ij';
    % the reduced covariance     
    Sig_ij = Sig([i j], [i j]);
    
    %% computer pwGCM from i to j
    
    % compute spectral factoraization of Sys_ij
    % compute pw GCM from i to j
    [~,Sigj] = ss_to_subsys(iiSys_ij,Sig_ij,2);
    Sigj_i = Sig_ij(2,2);
    F_JI(l) = log(Sigj) - log(Sigj_i);
    
    % compute pw GCM from j to i, i.e. F(i,j)
    [~,Sigi] = ss_to_subsys(iiSys_ij,Sig_ij,1);
    Sigi_j = Sig_ij(1,1);
    F_IJ(l) = log(Sigi) - log(Sigi_j);

end

F1 = zeros(n,n);
F1(U>0) = F_IJ;
F2 = zeros(n);
F2(U>0) = F_JI;
F = (F1 + F2')./(U+U');
F = F + diag(nan(1,n));
