function Sys = ss_eps_to_nu(Sys_all,Sys_i,Sys_j,i,j)
% build the ss model from input innovation eps to the
% inovation processs of the sub
A = Sys_all.A;
B = Sys_all.B;
C = Sys_all.C;
D = Sys_all.D;
% the inverse of ss model
iSys_i = sys_to_iss(Sys_i);
iSys_j = sys_to_iss(Sys_j);
n = size(A,1);
m = size(B,2);

Ai = iSys_i.A;
Bi = iSys_i.B;
Ci = iSys_i.C;
Di = iSys_i.D;

Aj = iSys_j.A;
Bj = iSys_j.B;
Cj = iSys_j.C;
Dj = iSys_j.D;


ijo = setdiff([1:m],[i j]);

B_ij_aug = zeros(2*n,m);
B_ij_aug(1:n,[i ijo]) = Bi;
B_ij_aug(n+1:end,[j ijo]) = Bj;
        
        
C_ij_aug = [Ci(1,:) zeros(1,n);
            zeros(1,n) Cj(1,:)];
            
Dij_aug = zeros(2,m);
Dij_aug(1,[i ijo]) = Di(1,:);
Dij_aug(2,[j ijo]) = Dj(1,:);


A_aug = [ A zeros(n,2*n);
        B_ij_aug*C blkdiag(Ai,Aj)];
% check stability
lambda = max(abs(eig(A_aug)));
if(lambda>=1)
    A_aug = A_aug/lambda*0.99;
end

B_aug = [B;B_ij_aug*D];


C_aug = [Dij_aug*C C_ij_aug];
D_aug = Dij_aug*D;

Sys.A = A_aug;
Sys.B = B_aug;
Sys.C = C_aug;
Sys.D = D_aug;
      
