function f = var_to_fcgcm_sent(A,Sig,freq,U)
%% f = var_to_fcgcm_sent(A,Sig,freq,U)
% frquency-domain cGCM using separated entropy minimization
%
%
%% Arguments
%
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     freq       frequeency
%     U          indicator for connections (optional)
% _output_
%
%     f          frequncy-domain cGCM
%
%
% (C) Lipeng Ning, 2020. 
%% 


n = size(A,1);
N = length(freq);

if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);
end
% full state-space model
Sys = var_to_ss(A);

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);


f_IJ = nan(L,N);
f_JI = nan(L,N);

for l = 1:L

    i = I(l);
    j = J(l);
    ijo = setdiff(1:n,[i j]);
    
    % computer two sub system 
    % tf from nu to y
    [Sys_i_ijo,~] = ss_to_subsys(Sys,Sig,[i ijo]);
    [Sys_j_ijo,~] = ss_to_subsys(Sys,Sig,[j ijo]);

    % build tf from eps to nu
    Sys_eps_to_nu = ss_eps_to_nu(Sys,Sys_i_ijo,Sys_j_ijo,i,j);
    
    %% spect-fact for nu_i,nu_j
    [Sys_ij,Sig_ij] = ss_to_subsys(Sys_eps_to_nu,Sig);


    %% frequency-domain cgcm
    H_ij = ss_to_H(Sys_ij,freq);
    iSys_ij = sys_to_iss(Sys_ij);
    G_ij = ss_to_H(iSys_ij,freq);
    %%
    f_ij_freq = zeros(1,N);
    f_ji_freq = zeros(1,N);
    
    for k = 1:N
        Sij = H_ij(:,:,k)*Sig_ij*H_ij(:,:,k)';
        Si_j = G_ij(1,1,k)\Sig_ij(1,1)/conj(G_ij(1,1,k));
        Sj_i = G_ij(2,2,k)\Sig_ij(2,2)/conj(G_ij(2,2,k));
        f_ij_freq(k) = log(real(Sij(1,1))) - log(real(Si_j));
        f_ji_freq(k) = log(real(Sij(2,2))) - log(real(Sj_i));
    end
     f_IJ(l,:) = f_ij_freq;% both (i,j) and (j,i) are compued within one loop to reduce computational time 
     f_JI(l,:) = f_ji_freq;
end
f1 = zeros(n*n,N);
f1(U>0,:) = f_IJ;
f1 = reshape(f1,n,n,N);
f2 = zeros(n*n,N);
f2(U>0,:) = f_JI;
f2 = reshape(f2,n,n,N);
f2 = permute(f2,[2,1,3]);
L = U';
f = (f1 + f2)./(U(:,:,ones(1,N))+L(:,:,ones(1,N)));


