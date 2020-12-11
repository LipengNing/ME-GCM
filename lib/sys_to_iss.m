function iSys = sys_to_iss(Sys)
% iSys = sys_to_iss(Sys)
% compute the state-space representation of inverse system
% Assume D is invertible
A = Sys.A;
B = Sys.B;
C = Sys.C;
D = Sys.D;

iSys.A = A-B/D*C;
iSys.B = B/D;
iSys.C = -D\C;
iSys.D = inv(D);
%iSys = ss(A-B/D*C,B/D,-D\C,inv(D),-1);
