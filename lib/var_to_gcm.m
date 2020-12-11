function F = var_to_gcm(A,Sig,U)
%% var_to_gcm
%
% Calculate pairwise time-domain Granger causality based on VAR parameters
%
%
%% Arguments
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     U          indicators for specific connections (optional)
% _output_
%
%     F          Granger causality matrix
%
% (C) Lipeng Ning, 2020. 

%% 


n = size(A,1);

% full regression

% full state-space model
Sys = var_to_ss(A);% information about input variance Sig is not included in Sys


% U: indicator of the connections to be estiamted
if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);% upper triangular matrix, the lower part is also estimated 
end
    

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

F_IJ = nan(1,L);
F_JI = nan(1,L);


for l = 1:L
    
    i = I(l);
    j = J(l);

    [Sys_ij,Sig_ij] = ss_to_subsys(Sys,Sig,[i j]);
    
    [~,Sig_ij_i] = ss_to_subsys(Sys_ij,Sig_ij,1);
    [~,Sig_ij_j] = ss_to_subsys(Sys_ij,Sig_ij,2);
    
    F_IJ(l) = log(Sig_ij_i) - log(Sig_ij(1,1));
    F_JI(l) = log(Sig_ij_j) - log(Sig_ij(2,2));%

end


F1 = zeros(n,n); F1(U>0) = F_IJ; 
F2 = zeros(n,n); F2(U>0) = F_JI;
F = F1+F2';  F = F./(U+U');

F = F + diag(ones(1,n)*nan);

