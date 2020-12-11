function Sys = var_to_ss(A)
% convert a MVAR model parameterized by AR and Omega to a SS model
[n,~,L] = size(A);
order = L;
N = n*L;
a = compan([1 zeros(1,order)]);
AA = kron(a,eye(n));
A = reshape(A,n,n*L);
AA(1:n,:) = A;

B = [eye(n);zeros(N-n,n)];
C = A;
D = eye(n);

Sys.A = AA;
Sys.B = B;
Sys.C = C;
Sys.D = D;

%Sys = ss(AA,B,C,D,-1);