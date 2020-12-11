function A = tf_to_itf(H)
% A = tf_to_itf(H)
%compute frequency-wise inverse of tf


[n,~,N] = size(H);
A = zeros(n,n,N);
for k = 1:N
    A(:,:,k) = inv(H(:,:,k));
end