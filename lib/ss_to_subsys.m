function [Sys_sub,Sigma] = ss_to_subsys(Sys,Sig,subid)
%compute the state-space model of the spectral factorization of a sub
%system of Sys
%input: Sys: the state-space model of the full system
% Sig: covariance of input signal
% subid: the indices of the output corresponding to the subsystem
% output: Sys_sub, the ss model of the spectra-factorization of the
% subsystem
% Sigma: the prediction error variance of Sys_sub
% if 
% (C) Lipeng Ning, 2020

n = size(Sys.C,1);% the dimension of time-series

if nargin<3
    subid = 1:n; % 
end


A = Sys.A';% transpose of system
B = Sys.B*real(sqrtm(Sig));
C = Sys.  C(subid,:);
D = Sys.D(subid,:)*real(sqrtm(Sig));

Q = B*B';
S = B*D';
R = D*D';

%[P,K,~] = idare(A,C',Q,R,S);% for latest Matlab versions
[P,~,K] = dare(A,C',Q,R,S);

Sigma = C*P*C'+R;

Sys_sub.A = A';
Sys_sub.B = K';
Sys_sub.C = C;
Sys_sub.D = eye(size(R));

