function F = var_to_cgcm_std(A,Sig,U)
%% var_to_pwcgc_std
%
% Calculate the standard frequency-domain conditional GCM based on VAR parameters using the method from Geweke84
%
%
%% Arguments
% _input_
%
%     A          multivariate AR coefficients
%     Sig        residual covariance matrix  
%     U          indicators for selected connections (optional)
% _output_
%
%     f         frequency_domain cGCM
%
% (C) Lipeng Ning, 2020. 

%% 


n = size(A,1);

% full state-space model
Sys = var_to_ss(A);% information about input variance Sig is not included in Sys

% U: indicator of the connections to be estiamted
if(~exist('U','var'))
    U = ones(n,n) - eye(n);% upper triangular matrix, the lower part is also estimated 
end

[I,J] = ind2sub([n,n],find(U>0));

L = length(I);

F_IJ = zeros(1,L);

for l = 1:L
    LSIG = log(diag(Sig));

    i = I(l);
    j = J(l);
    % reduced regression
    
    jo = [1:j-1 j+1:n]; % omit j
    
    [~,Sigj] = ss_to_subsys(Sys,Sig,jo);% system removing j 
    LSIGj = log(diag(Sigj));
    
    ii = jo==i;
    F_IJ(l) = LSIGj(ii)-LSIG(i);
end

F = nan(n);

F(U>0) = F_IJ;
end
