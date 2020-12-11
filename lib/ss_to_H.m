function H = ss_to_H(Sys,freq)
% transform discrete-time sate-space system Sys to frequency-domain tf
% input: Sys   state-space model
%        freq   evaluation frequency
% output: H   freq-resp function


n = size(Sys.D,1);
N = length(freq);
H = nan(n,n,N);
m = size(Sys.A,1);
i = sqrt(-1);
for k = 1:N
    f = exp(i*freq(k));
    H(:,:,k) = Sys.D + Sys.C/(f*eye(m)-Sys.A)*Sys.B;
end