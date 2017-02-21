function [y_noi,noise]=addnoise_gral(y,RMS_noise,seed_w)
% INPUT: 
%   y           : original signal without noise
%   RMS_noise   : RMS of the noise (in percent) (units compatible with y)
%   seed_w      : seed numbers to generate the measurement noise
% OUTPUT:
%   y_noise     : noisy signal 

m=size(y,2);                        % number of outputs
N=size(y,1);                        % number of data samples
y_noi=zeros(N,m);
for jj=1:m
    % Noisy measurement
    randn('state',seed_w(jj,1));                        % Seed to generate noise
    y_noi(:,jj)=y(:,jj)'+RMS_noise(jj)*randn(1,N);      % Measurement + noise
end
noise=y_noi-y;                      % noise vector