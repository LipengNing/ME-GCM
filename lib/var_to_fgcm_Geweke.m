function f = var_to_fgcm_Geweke(A,Sig,freq,U)
%% var_to_fgcm_ent
%
% Calculate pairwise frequency-domain GCM based on VAR parameters
% using the method by Geweke
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

% if(~exist('U','var'))
%     U = toeplitz([0 ones(1,n-1)],zeros(1,n));% upper triangular matrix, the lower part is also estimated 
% end


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

    [Hi,Hj] = stable_fGCM_factor(Sys_ij,Sig_ij,freq);
     H_ij = ss_to_H(Sys_ij,freq);
    
    f_ij_freq = zeros(1,N);
    f_ji_freq = zeros(1,N);
    
    for k = 1:N
        
        S_ii = H_ij(1,:,k)*Sig_ij*H_ij(1,:,k)';
        S_jj = H_ij(2,:,k)*Sig_ij*H_ij(2,:,k)';
        
        % S_i_j = (H_ij(1,1)+H_ij(1,2)*Sig_ij(2,1)/Sig_ij(1,1))*Sig_ij(1,1)*conj((H_ij(1,1)+H_ij(1,2)*Sig_ij(2,1)/Sig_ij(1,1)));
        % the following step ensures that
        % (H_ij(1,1)+H_ij(1,2)*Sig_ij(2,1)/Sig_ij(1,1)) is stable
        
        S_i_j = Hi(1,1,k)*Sig_ij(1,1)*conj(Hi(1,1,k));   
        S_j_i = Hj(1,1,k)*Sig_ij(2,2)*conj(Hj(1,1,k));
        
        f_ij_freq(k) = real(log(real(S_ii)) - log(real(S_i_j)));
        f_ji_freq(k) = real(log(real(S_jj)) - log(real(S_j_i)));
        
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


function [Hi,Hj] = stable_fGCM_factor(Sys_ij,Sig_ij,freq)
% stale verison of (H_ij(1,1)+H_ij(1,2)*Sig_ij(2,1)/Sig_ij(1,1)) for
% computing fGCM_Gewk
% buid the SS for (H_ij(1,1)+H_ij(1,2)*Sig_ij(2,1)/Sig_ij(1,1))
Sys_ij_i.A = blkdiag(Sys_ij.A,Sys_ij.A);
Sys_ij_i.B = [Sys_ij.B(:,1);Sys_ij.B(:,2)*Sig_ij(2,1)/Sig_ij(1,1)];
Sys_ij_i.C = [Sys_ij.C(1,:) Sys_ij.C(1,:)];
Sys_ij_i.D = Sys_ij.D(1,1)+Sys_ij.D(1,2)*Sig_ij(2,1)/Sig_ij(1,1);

z=tzero(Sys_ij_i.A,Sys_ij_i.B,Sys_ij_i.C,Sys_ij_i.D);
if(max(abs(z))>1)
    Sys_ij_i = ss_to_subsys(Sys_ij_i,1);% minimum phase formulation
end
Sys_ij_i = sys_to_iss(Sys_ij_i);
Hi = ss_to_H(Sys_ij_i,freq);


Sys_ij_j.A = blkdiag(Sys_ij.A,Sys_ij.A);
Sys_ij_j.B = [Sys_ij.B(:,2);Sys_ij.B(:,1)*Sig_ij(1,2)/Sig_ij(2,2)];
Sys_ij_j.C = [Sys_ij.C(2,:) Sys_ij.C(2,:)];
Sys_ij_j.D = Sys_ij.D(2,2)+Sys_ij.D(2,1)*Sig_ij(1,2)/Sig_ij(2,2);

z=tzero(Sys_ij_j.A,Sys_ij_j.B,Sys_ij_j.C,Sys_ij_j.D);
if(max(abs(z))>1)
    Sys_ij_j = ss_to_subsys(Sys_ij_j,1);% minimum phase formulation
end
Hj = ss_to_H(Sys_ij_j,freq);

end
