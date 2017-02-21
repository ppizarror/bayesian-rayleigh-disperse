function [y]=motion(iGMfile,iGMdirection,iGMfact,N,exact_p,outputs,col_output,Fs,folder,NstepG)
% INPUT:
%       GMfile          : name of the input file
%       ug              : input motion data
%       exact_p         : vector with the exact values of the parameters
%       outputs         : name of the measured outputs
%       col_output      : ID with the column # of the measured output
%       Fs              : sampling frequency
% OUTPUT:
%       y               : measured output (without noise)
%       record_exact    : recorded outputs for the exact model
% October 9, 2014 (recorders were eliminated)

n=length(exact_p);                      % number of parameters to be identified
m=length(outputs);                      % number of outputs
%%
%% RUN OPENSEES TO GET THE SIMULATED RESPONSE
exactID=fopen('Parameters.txt','w');
for jj=1:n
    fprintf(exactID,'%s\n',['set p',num2str(jj),' ',num2str(exact_p(jj))]);
end
fprintf(exactID,'%s %i \n','set Nsteps',N);             fprintf(exactID,'%s %f \n','set dt',1/Fs);
fprintf(exactID,'%s \n',['set iGMfile ', iGMfile]);     fprintf(exactID,'%s \n',['set iGMdirection ', iGMdirection]);
fprintf(exactID,'%s \n',['set iGMfact ', iGMfact]);     fprintf(exactID,'%s \n',['set WorkingDir "', folder,'"']);        
fprintf(exactID,'%s %i \n','set NstepG', NstepG);       fprintf(exactID,'%s %f \n','set TmaxAnalysis', N/Fs);
fclose(exactID);
system('OpenSees THAnalysis.tcl');
%%
%% MEASUREMENT
y=zeros(N,m);
for jj=1:m
    meas=load(['.\',folder,'\output',num2str(jj),'.out'],'-ascii');
    y(:,jj)=meas(NstepG+1:end,col_output(jj));             % load simulated data and other parameters
end
%%