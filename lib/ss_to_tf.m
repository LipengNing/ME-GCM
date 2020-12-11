function H = ss_to_tf(Sys,fres)
% H = ss_to_tf(Sys,Sig,fres)
% input: Sys: the state space parameter with scaled B matrix
%        fres:the number of frequency points between 0 and pi
% output: H: the frequency response function


freq = linspace(0,pi,fres+1);
% compute the discrete-time state-space model
s = ss(Sys.A, Sys.B, Sys.C, Sys.D,-1);
H = freqresp(s,freq);
