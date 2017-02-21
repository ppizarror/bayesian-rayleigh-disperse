function [xhat_k,Pk,yhat_k,rhat_k,Pk_rr]=dukf_fem_normalized(iGMfile,iGMdirection,iGMfact,f,xhat_k_1,Pk_1,yk,Q,rhat_k_1,Pk_1_rr,U,T,alpha,kappa,beta,min_lim,max_lim,new_sp,kk,workdir,col_output,exact_p,Fs,NstepG,adaptive,cumulative,step_update)
%% This function estimate unknown time-invariant parameters of a FE model
% residing in OpenSees. It includes the capability of estimate the diagonal
% entries of the measurement noise covariance (adaptive) and has the option
% to do non-sequential updating, both non-cumulative and cumulative
%% Inputs:   
%           GMfile      : name of the input motion
%           GMfact      : factor for units of GMfile
%           f           : function handle for f(x)
%           xhat_k_1    : "a priori" state estimate
%           Pk_1        : "a priori" estimated state covariance
%           yk          : current measurement
%           Q           : process noise covariance 
%           r_hat_k_1   : noise level to construct measurement noise covariance R (initial for adaptive)
%           Pk_1_rr     : initial covariance matrix of the estimate of diagonal entries of R
%           U           : measurement noise covariance matrix of slave filter (for adaptive)
%           T           : process noise covariance of slave filter (for adaptive)
%           u1          : input to state equation
%           alpha       : parameter for the unscented transformation
%           kappa       : parameter for the unscented transformation
%           beta        : parameter for the unscented transformation
%           min_lim     : lower bounds for parameter updating (if the value of
%                         parameter < min_lim --> keep previous time step
%           max_lim     : upper bounds for parameter updating (if the value of
%                         parameter > max_lim --> keep previous time step
%           new_sp      : 'yes' to generate new sigma points in observation prediction
%                         'no' to use sigma points generated in state prediction
%           kk          : time step
%           workdir     : name of working folder
%           recorder    : names of variables to be saved
%           col_output  : column IDs of the output (mx1)
%           mod_par     :
%           exact_p     : exact values for the parameters
%           Fs          : sampling frequency
%           NstepG      : number of time steps to apply gravity in FE model
%           adaptive    : 'yes' for adaptive filtering (estimate measurement noise covariance as well)
%           cumulative  : 'yes' to use the augmented innovation when using non-sequential updating
%           step_update : number defining step size for filtering (update parameters) (1: update at every time step)

%% Output:   xhat_k      : "a posteriori" state estimate
%           Pk          : "a posteriori" state covariance
%           yhat_k      : output estimate
%           rhat_k      : "a posteriori" estimate of diagonal entries of measurement noise covariance
%           Pk_rr       : "a posteriori" state covariance of diagonal entries of measurement noise covariance

%% FEM   :   h           : the nonlinear function h(x) used to evaluate the
%                         sigma points in the observation prediction step,
%                         is represented by the output of the FEM for each
%                         set of sigma points Xhat_k(i)
%%
n=numel(xhat_k_1);                              % Number of states
if strcmp(cumulative,'yes')
    k=step_update;                              % Number defining step size for filtering (updating)
    m=numel(yk)/k;                              % Number of measurements
else
    m=numel(yk);                                % Number of measurements
end
lambda=alpha^2*(n+kappa)-n;
gamma=sqrt(n+lambda);
%% WEIGHTS OF MEAN AND COVARIANCE
Wm=[lambda/(n+lambda) 0.5/(n+lambda)+zeros(1,2*n)];
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);
%%
%% SIGMA POINTS OF THE STATE
Xhat_k_1=sp(xhat_k_1,Pk_1,gamma);
L=size(Xhat_k_1,2);                             % Number of sigma points
%% 
%% TIME UPDATE (PREDICT)
% 1) State prediction
% Project sigma points through the nonlinear fucntion f
Xhat_k=zeros(n,L);
for i=1:L
    Xhat_k(:,i)=f(Xhat_k_1(:,i));
end
% Unscented transformation
[xhat_k_,Pk_,dif_xhat,~]=ut(Xhat_k,Wm,Wc,Q);    % Unscented transformation state
% 2) Observation prediction
% Sigma points
if strcmp(new_sp,'yes')
    Xhat_k=sp(xhat_k_,Pk_,gamma);               % Generate new sigma points
end
%%
%% FINITE ELEMENT MODEL INTERACTION
% Transform the sigma points through the nonlinear function h (FEM)
Yhat_k=zeros(length(yk),L);
parfor i=1:L
    % Generate a TXT file, read it, and assign the parameters in OpenSees
    % Each sigma point (Xhat_k(:,i)) is used as input for the FEM and the
    % corresponding response Yhat_k(:,i) is obtained
    nombre=['Parameters',num2str(i),'.txt'];
    parfileID=fopen(nombre,'w');
    for jj=1:n
        fprintf(parfileID,'%s\n',['set p',num2str(jj),' ',num2str(Xhat_k(jj,i)*exact_p(jj))]);  %% times p_exact
    end
    fprintf(parfileID,'%s %i \n','set Nsteps',kk);            fprintf(parfileID,'%s %f \n','set dt',1/Fs);
    fprintf(parfileID,'%s \n',['set iGMfile ', iGMfile]);     fprintf(parfileID,'%s \n',['set iGMdirection ', iGMdirection]);
    fprintf(parfileID,'%s \n',['set iGMfact ', iGMfact]);     fprintf(parfileID,'%s\n',['set WorkingDir "', workdir, num2str(i),'"']);
    fprintf(parfileID,'%s %i \n','set NstepG', NstepG);       fprintf(parfileID,'%s %f \n','set TmaxAnalysis', kk/Fs);
    fclose(parfileID);
    % RUN OpenSees
    system(['OpenSees THAnalysis',num2str(i),'.tcl']);
end
% Save measured responses (y(k)) and assign it to Yhat_k
for i=1:L
    for mm=1:m
        resp=load(['./',workdir,num2str(i),'/output',num2str(mm),'.out'],'-ascii');
        if strcmp(cumulative,'yes')
            Yhat_k(k*mm-k+1:mm*k,i)=resp(NstepG+kk-k+1:NstepG+kk,col_output(mm));
        else
            Yhat_k(mm,i)=resp(NstepG+kk,col_output(mm));
        end
    end
end
%%
%% Unscented transformation
if strcmp(cumulative,'yes')
    [yhat_k,Pk_y,dif_yhat,Pk_y2]=ut(Yhat_k,Wm,Wc,rhat_k_1);                 % Unscented transformation measurement
else
    [yhat_k,Pk_y,dif_yhat,Pk_y2]=ut(Yhat_k,Wm,Wc,diag(rhat_k_1));           % Unscented transformation measurement    
end
%%

%% NOISE ESTIMATION
if ~strcmp(adaptive,'yes')
    rhat_k=rhat_k_1;
    Pk_rr=[];
else
    e_k=yk-yhat_k;                              % Innovation vector
    z_k=diag(e_k*e_k');                         % measurement of the slave filter
    vec_k=diag(Pk_y2);
    rhat_k_=rhat_k_1;
    Pk_rr_=Pk_1_rr+T;
    zhat_k=rhat_k_+vec_k;
    Pk_zz=Pk_rr_+U;
    Pk_rz=Pk_rr_;
    Kk_r=Pk_rz/Pk_zz;
    rhat_temp=rhat_k_+Kk_r*(z_k-zhat_k);
    if sum(rhat_temp>0)==m
        rhat_k=rhat_temp;
        Pk_rr=Pk_rr_-Kk_r*Pk_zz*Kk_r';
    else
        rhat_k=rhat_k_;
        Pk_rr=Pk_rr_;
    end
    Pk_y=Pk_y2+diag(rhat_k);                    % Update measurement covariance based on measurement noise estimation
end
%%

%% MEASUREMENT UPDATE (CORRECT)
% 1) Cross-covariance
Pk_xy=dif_xhat*diag(Wc)*dif_yhat';
% 2) Kalman gain
Kk=Pk_xy/Pk_y;
% 3) Update state estimate
xhat_k_temp=xhat_k_+Kk*(yk-yhat_k);
xhat_k=zeros(size(xhat_k_temp));
% 4) Update error covariance
Pk_temp=Pk_-Kk*Pk_y*Kk';
Pk=zeros(size(Pk_temp));
% Check difference
if sum(xhat_k_temp>min_lim)==n && sum(xhat_k_temp<max_lim)==n
    xhat_k=xhat_k_temp;
    Pk=Pk_temp;
else
    display('Out of range in update process')
    xhat_k_temp
    xhat_k=xhat_k_1;
    Pk=Pk_1;
end
%%


%% ------------------------------------------------------------------------
function X=sp(x,P,gamma)
% Sigma points around reference point
% Inputs:
%       x       : reference point
%       P       : covariance
%       gamma   : coefficient
% Output:
%       X       : Sigma points
A = gamma*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------
function [xk_,Pk_,dif_x,Pk_2]=ut(Xk,Wm,Wc,Q)
% Unscented Transformation
% Input:
%        X      : sigma points
%       Wm      : weights for mean
%       Wc      : weights for covariance
%        Q      : additive covariance
% Output:
%        xk_    : transformed mean
%        Pk_    : transformed covariance
n=size(Xk,1);   % Lenght of state vector
L=size(Xk,2);   % Number of sigma points
xk_=zeros(n,1);
for i=1:L
    xk_=xk_+Wm(i)*Xk(:,i);
end
dif_x=Xk-xk_(:,ones(1,L));
Pk_=dif_x*diag(Wc)*dif_x'+Q;
Pk_2=dif_x*diag(Wc)*dif_x';
%% ------------------------------------------------------------------------