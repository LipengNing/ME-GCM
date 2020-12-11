function F = var_to_cgcm_sent(A,Sig,U)
%% F = var_to_cgcm_sent(A,Sig,U)
% cGCM using separated entropy-based regression
%
%
%% Arguments
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     U          indicator for connections (optional)
% _output_
%
%     F          Granger causality matrix
%
% (C) Lipeng Ning, 2020. 
%% 

n = size(A,1);

if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);
end

% full state-space model
Sys = var_to_ss(A);

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

F_IJ = nan(1,L);
F_JI = nan(1,L);


for l = 1:L

    i = I(l);
    j = J(l);
    ijo = setdiff(1:n,[i j]);
    
    % computer two sub system 
    % tf from nu to y
    [Sys_i_ijo,Sig_i_ijo] = ss_to_subsys(Sys,Sig,[i ijo]);
    [Sys_j_ijo,Sig_j_ijo] = ss_to_subsys(Sys,Sig,[j ijo]);
    sig_i_ijo = Sig_i_ijo(1,1); % residual of nu_i and nu_j
    sig_j_ijo = Sig_j_ijo(1,1);
    % build tf from eps to nu
    Sys_eps_to_nu = ss_eps_to_nu(Sys,Sys_i_ijo,Sys_j_ijo,i,j);
    
    %% spect-fact for nu_i,nu_j
    [~,Sig_ij] = ss_to_subsys(Sys_eps_to_nu,Sig);
    sig_i_ijo_j = Sig_ij(1,1);
    sig_j_ijo_i = Sig_ij(2,2);
    
    F_IJ(l) = log(sig_i_ijo) - log(sig_i_ijo_j);
    
    F_JI(l) = log(sig_j_ijo) - log(sig_j_ijo_i);
end

F1 = zeros(n,n);
F1(U>0) = F_IJ;
F2 = zeros(n);
F2(U>0) = F_JI;
F = (F1 + F2')./(U+U');
F = F + diag(nan(1,n));

