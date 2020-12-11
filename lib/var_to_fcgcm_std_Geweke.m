function f = var_to_fcgcm_std_Geweke(A,Sig,freq,U)
%% f = var_to_fcgcm_std_Geweke(A,Sig,freq,U)
%
% Calculate frequency-domain conditionoal GCM based on VAR parameters using
% the ME method
%
%% Arguments
%
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     freq       frequency 
%     U          indicator matrix of the connections (optional)
% _output_
%
%     f         frequency_domain cGCM
%
% (C) Lipeng Ning, 2020. 
%% 

n = size(A,1);
N = length(freq);

% full state-space model
Sys = var_to_ss(A);% information about input variance Sig is not included in Sys

%%
if(~exist('U','var'))
    U = ones(n) - eye(n);
end

[I,J] = ind2sub([n,n],find(U>0));
L = length(I);

f_sub = nan(L,N);
H = ss_to_H(Sys,freq);

for l = 1:L
    i = I(l);
    j = J(l);
    
    jo = [1:j-1 j+1:n]; % omit j
    
    [Sysj,Sigj] = ss_to_subsys(Sys,Sig,jo);% system removing j 
    LSIGj = log(diag(Sigj));
    
    %Hj = SS2H(Sysj,freq);
    iSysj = sys_to_iss(Sysj);
    Gj = ss_to_H(iSysj,freq);
    
    io = [1:i-1 i+1:n];
    ii = find(jo==i);
    f_ij_freq = zeros(1,N);
        
    for k = 1:N
        Gj_k_aug = zeros(1,n);
        Gj_k_aug(jo) = Gj(ii,:,k);
        H_tilde = Gj_k_aug*H(:,:,k);
        Sig_jk_i = Sig(io,io) - Sig(io,i)/Sig(i,i)*Sig(i,io);
        P_new = Sigj(ii,ii) - H_tilde(io)*Sig_jk_i*H_tilde(io)';
        f_ij_freq(k) = LSIGj(ii) - log(real(P_new));
    end
    f_sub(l,:) = f_ij_freq;
end
f = nan(n*n,N);
f(U>0,:) = f_sub;
f = reshape(f,n,n,N);

end

%%


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
