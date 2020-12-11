function f = var_to_fcgcm_std_ent(A,Sig,freq,U)
%% f = var_to_fcgcm_std_ent(A,Sig,freq,U)
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

f_IJ = zeros(L,N);

iSys = sys_to_iss(Sys);

G = ss_to_H(iSys,freq);

for l = 1:L
    i = I(l);
    j = J(l);
    % reduced regression
    
    jo = [1:j-1 j+1:n]; % omit j
    
    [Sysj,Sigj] = ss_to_subsys(Sys,Sig,jo);% system removing j 
    
    %Hj = SS2H(Sysj,freq);
    iSysj = sys_to_iss(Sysj);
    
    Gj = ss_to_H(iSysj,freq);
    
    ii = find(jo==i);

    f_ij_freq = zeros(1,N);
        
    for k = 1:N
        P_i_io = G(i,i,k)\Sig(i,i)/(conj(G(i,i,k)));
        P_i_ijo = Gj(ii,ii,k)\Sigj(ii,ii)/(conj(Gj(ii,ii,k)));
        f_ij_freq(k) = log(real(P_i_ijo)) - log(real(P_i_io));
    end
    f_IJ(l,:) = f_ij_freq;
    
end

f = nan(n*n,N);
f(U>0,:) = f_IJ;
f = reshape(f,n,n,N);
end

