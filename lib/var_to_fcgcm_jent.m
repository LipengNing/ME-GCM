function f = var_to_fcgcm_jent(A,Sig,freq,U)
%% f = var_to_fcgcm_jent(A,Sig,freq,U)
%
% Calculate conditional frequency-domain GCM (cGCM) based on VAR parameters
% using joint entropy minimization
%
%% Arguments
%
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     freq       frequency resolution
%     U          indicator for connctions (optional)
% _output_
%
%     f          frequency-domain cGCM
%
% (C) Lipeng Ning, 2020. 
%% 


n = size(A,1);
N = length(freq);
%%
if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);
end
% full state-space model
Sys = var_to_ss(A);
iSys = sys_to_iss(Sys); % the G matrix


[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

f_IJ = nan(L,N);
f_JI = nan(L,N);


for l = 1:L

    i = I(l);
    j = J(l);
    
    G = ss_to_H(iSys,freq); 
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
    H_ij = ss_to_H(iiSys_ij,freq); % inverse of G_ij(theta)
    
    Sig_ij = Sig([i j], [i j]);
    
%     
    f_ij_freq = zeros(1,N);
    f_ji_freq = zeros(1,N);
    
    for k = 1:N
        Sij = H_ij(:,:,k)*Sig_ij*H_ij(:,:,k)';
        Si_j = G(i,i,k)\Sig(i,i)/conj(G(i,i,k));
        Sj_i = G(j,j,k)\Sig(j,j)/conj(G(j,j,k));
        f_ij_freq(k) = log(real(Sij(1,1))) - log(real(Si_j));
        f_ji_freq(k) = log(real(Sij(2,2))) - log(real(Sj_i));
    end
     f_IJ(l,:) = f_ij_freq;
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
end

    
    
    
    
