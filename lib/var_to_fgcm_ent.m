function f = var_to_fgcm_ent(A,Sig,freq, U)
%% var_to_fgcm_ent
%
% Calculate pairwise frequency-domain GCM based on VAR parameters
% using the minimum-entropy method
%
%% Arguments
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     freq       frequency 
% _output_
%
%     f          frequency-band (0.01 0.1)Hz averaged valueds
%% Description
%
% (C) Lipeng Ning, 2020. 
%% 


n = size(A,1);
N = length(freq);
% U: indicator of the connections to be estiamted
if(~exist('U','var'))
    U = toeplitz(zeros(1,n),[0 ones(1,n-1)]);% upper triangular matrix, the lower part is also estimated 
end

% full state-space model
Sys = var_to_ss(A);% information about input variance Sig is not included in Sys

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

f_IJ = zeros(L,N);
f_JI = zeros(L,N);

for l = 1:L
    
    i = I(l);
    j = J(l);

    % reduced regression
    
    [Sys_ij,Sig_ij] = ss_to_subsys(Sys,Sig,[i j]);

    H_ij = ss_to_H(Sys_ij,freq);
    
    iSys_ij = sys_to_iss(Sys_ij);
    
    
    %G(i,i)
    iSys_ij_i.A = iSys_ij.A; 
    iSys_ij_i.B = iSys_ij.B(:,1); 
    iSys_ij_i.C = iSys_ij.C(1,:); 
    iSys_ij_i.D = iSys_ij.D(1,1);
    iiSys_ij_i = sys_to_iss(iSys_ij_i);
    
    if(max(abs(eig(iiSys_ij_i.A)))>1)
        iSys_ij_i_mp = ss_to_subsys(iSys_ij_i,1);
        iiSys_ij_i = sys_to_iss(iSys_ij_i_mp);
    end
    iG_ij_i = ss_to_H(iiSys_ij_i,freq);
    
    %G(j,j)
    iSys_ij_j.A = iSys_ij.A; 
    iSys_ij_j.B = iSys_ij.B(:,2); 
    iSys_ij_j.C = iSys_ij.C(2,:); 
    iSys_ij_j.D = iSys_ij.D(2,2);
    iiSys_ij_j = sys_to_iss(iSys_ij_j);
    if(max(abs(eig(iiSys_ij_j.A)))>1)
        iSys_ij_j_mp = ss_to_subsys(iSys_ij_j,1);
        iiSys_ij_j = sys_to_iss(iSys_ij_j_mp);
    end
    iG_ij_j = ss_to_H(iiSys_ij_j,freq); 
    
    f_ij_freq = zeros(1,N);
    f_ji_freq = zeros(1,N);
    
%     for k = 1:N
%         S_ii = H_ij(1,:,k)*Sig_ij*H_ij(1,:,k)';
%         S_jj = H_ij(2,:,k)*Sig_ij*H_ij(2,:,k)';
%         
%         G_ij = inv(H_ij(:,:,k));
%         S_i_j = G_ij(1,1)\Sig_ij(1,1)/conj(G_ij(1,1));
%         S_j_i = G_ij(2,2)\Sig_ij(2,2)/conj(G_ij(2,2));
%         
%         f_ij_freq(k) = real(log(S_ii) - log(S_i_j));
%         f_ji_freq(k) = real(log(S_jj) - log(S_j_i));
%     end
    S_ii = zeros(1,N);
    S_jj = zeros(1,N);
    S_i_j = zeros(1,N);
    S_j_i = zeros(1,N);
    
    for k = 1:N
        S_ii(k) = H_ij(1,:,k)*Sig_ij*H_ij(1,:,k)';
        S_jj(k) = H_ij(2,:,k)*Sig_ij*H_ij(2,:,k)';
        
        S_i_j(k) = iG_ij_i(1,1,k)*Sig_ij(1,1)*conj(iG_ij_i(1,1,k));
        S_j_i(k) = iG_ij_j(1,1,k)*Sig_ij(2,2)*conj(iG_ij_j(1,1,k));
        
        f_ij_freq(k) = real(log(S_ii(k)) - log(S_i_j(k)));
        f_ji_freq(k) = real(log(S_jj(k)) - log(S_j_i(k)));
    end
%     
    
    
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
