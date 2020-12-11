function Sys = AR2SS(AR,Omega)
% convert a MVAR model parameterized by AR and Omega to a SS model
% (C) Lipeng Ning (2020)
[n,N] = size(AR);
order = round(N/n-1);
a = compan([1 zeros(1,order)]);
A = kron(a,eye(n));
A(1:n,:) = - AR(:,n+1:end);

B = [eye(n);zeros(N-2*n,n)]*real(sqrtm(Omega));
C = [eye(n) zeros(n,N-2*n)];
D = zeros(n);
Sys = ss(A,B,C,D,-1);