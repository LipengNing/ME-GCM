function [Sys_sub,Sigma] = SubsysSpectfact(Sys,subid)
%compute the state-space model of the spectral factorization of a sub
%system of Sys
%input: Sys: the state-space model of the full system
% subid: the indices of the output corresponding to the subsystem
% output: Sys_sub, the ss model of the spectra-factorization of the
% subsystem
% Sigma: the prediction error variance of Sys_sub
% (C) Lipeng Ning, 2020
A = Sys.A';% transpose of system
B = Sys.B;
C = Sys.C(subid,:);
D = Sys.D(subid,:);

Q = B*B';
S = B*D';
R = D*D';

[P,~,K] = dare(A,C',Q,R,S);
Sigma = C*P*C'+R;

Sys_sub.A = A';
Sys_sub.B = K'*real(sqrtm(Sigma));
Sys_sub.C = C;
Sys_sub.D = real(sqrtm(Sigma));

