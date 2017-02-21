function [rms_error]=rms_error(sig1,sig2)
N=length(sig1);
rms_error=roundn(sqrt(1/N*sum((sig1-sig2).^2))/sqrt(1/N*sum(sig1.^2))*100,-2);