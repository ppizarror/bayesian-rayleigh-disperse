clear all; clc;  close all;
% matlabpool open 10;

% profile on;
%% -------------------------------------------------------------------
%% ---------------------------- VARIABLES ----------------------------
%% -------------------------------------------------------------------
%% GENERAL
tic;
counter=1;
min_lim=0.2;                        % Lower bound in updating the parameters wrt initial estimate (in %)
max_lim=10;                          % Upper bound in updating the parameters wrt initial estimate (in %)

name='ID_1013';                         % name of file to be saved

% reg1='Gatos_EW_mod';   reg2='Gatos_NS_mod';

% earth1=importdata([reg1,'.txt']); earth2=importdata([reg2,'.txt']);          % load earthquake input

% iGMfile=['"',reg1,'.txt ',reg2,'.txt"'];
% iGMdirection='"1 3"';
% iGMfact='"1 1"';                     % Input-motion scaling factor (with units)
% Fs=50;                              % Sampling frequency
% NstepG=4;                           % Number of time steps to apply gravity
% workdir='Results';                  % Name of folder for results
% step_update=[3];                   % number defining step size to filtering (update parameters) (1: update at every time step)
% cumulative='yes';                  % 'yes' to consider cumulative case
% adaptive='no';                     % 'yes' if the measurement noise covariance want ot be estimated

%% PARAMETERS TO BE IDENTIFIED AND PROCESS NOISE COVARIANCE (Q)
params={'Es','fs','R0','b', 'Ec','fpc','eps0', 'alpha_m','beta_K'};       % Material parameters to be estimated (nx1)
exact_p=[2.0e8,4.14e5,18,0.05, 2.76e7,-4.00e4,-0.0035, 0.1137,0.0026];%0.0931,0.004];%       % Exact material parameters (nx1)
cov_q=[1e-6*ones(1,length(params))];             % Normalized coefficient of variation of process noise (N_q x n)

%% INITIALIZATION OF THE FILTER
x_ini=1.4*ones(length(params),1)'.*exact_p;                  % Initial state (N_x x n) (e.g. [0.9*60,0.9*29000,0.9*0.01,0.9*18; 0.8*60,0.8*29000,0.9*0.01,0.8*18]
cov=0.1*ones(1,length(params));              	% Coefficient of variation (N_x x n) (e.g. [0.3,0.3,0.3,0.3; 0.15,0.15,0.15,0.15]; 
x_ini=bsxfun(@rdivide,x_ini,exact_p);     % Normalize the initial state

%% MEASURED OUTPUTS AND MEASUREMENT NOISE COVARIANCE (R)
outputs={'A2_L_a','A4_L_a','A5_L_a', 'A2_T_a','A4_T_a','A5_T_a'};  % Measured output (mx1)
% col_output=[2,2,2, 4,4,4];                     % IDs of the measured outputs (mx1)                                                                         e.g. [2;2] (2 outputs)
RMS_noise=0.1*ones(1,length(outputs));                    % RMS of measurement noise (in m/s^2) (N_r x m)                                                               e.g. [1 2] (2 outputs 1 noise combination)
RMS_0=0.05*ones(1,length(outputs));                          % RMS to construct measurement noise covariance R (N_r x m) (initial for adaptive)                          e.g. [5 8] (2 outputs 1 noise combination)
coef_r=0.2;                         % CoV to compute initial covariance of r (initial for adaptive) (1 x 1)
% Parameter for adaptive filtering
U=diag(1e-6*ones(1,length(outputs)));                 % Measurement noise covariance of slave filter (mxm) (for adaptive)                             e.g. diag(1e-6*[1 1]) (2 outputs 1 noise combination)
T=diag(1e-8*ones(1,length(outputs)));                 % Process noise covariance of slave filter (mxm) (for adaptive)                                 e.g. diag(1e-8*[1 1]) (2 outputs 1 noise combination)
seed_w=[100:1:100+length(outputs)-1]';                       % Seed to generate the mesurement noise (mx1)                                                            e.g. [100;101]        (2 outputs 1 noise combination)

%% UKF ALGORITHM
alpha=0.01;  % 0.1
kappa=0;    % 3
beta=2;     % 0
new_sp='no';

%% LENGTHS OF VARIABLES
N_noise_level=size(RMS_noise,1);  % Number of different noise level (actual) for measurement
N_q=size(cov_q,1);                  % Number of different process noises (matrix Q)
N_r=size(RMS_0,1);                     % Number of different measurement noises (matrix R)
N_x_ini=size(x_ini,1);              % Number of different initial parameter estimates
N_cov_ini=size(cov,1);              % Number of different initial parameter covariance
n_step_update=length(step_update);  % Number of different step updates

%% LENGTHS OF VECTORS
n=length(params);                   % # of parameters
m=length(outputs);                  % # of measurements
N=length(earth1);                       % # of data samples

%% FIXED VARIABLES
% Process equation
f=@(x)[eye(n)*x];
% Pre-allocation matrices
x_est = zeros(n,N);                     % Parameter estimate 
P_est = zeros(n,n,N);                   % State covariance
y_est = zeros(m,N);                     % Output estimate
r_est = zeros(m,N);                     % Measurement noise covariance estimate
if strcmp(adaptive,'yes')
    Pr_est = zeros(m,m,N);                  % Covariance estimate of r_est
else
    Pr_est=[];
end
time=[0:N-1]/Fs;

%% SIMULATED DATA (without noise)
[y]=motion(iGMfile,iGMdirection,iGMfact,N,exact_p,outputs,col_output,Fs,'Exact',NstepG);

%% ANALYSIS
for c_steps=1:n_step_update
    t_ini=tic;
    ti_cpu=cputime;
    N_max_update=floor(N/step_update(c_steps))*step_update(c_steps);
    if strcmp(cumulative,'yes')
        y_est = zeros(m*step_update(c_steps),N);
    end
    for c_q=1:N_q
        q=cov_q(c_q,:);                      % Normalized standar deviation of process noise (nx1)
        Q=diag(q.^2)                         % Covariance of process noise
        for c_noise_level=1:N_noise_level
            %% SIMULATE DATA (add noise)
            [y_noi,noise]=addnoise_gral(y,RMS_noise(c_noise_level,:),seed_w);
            figure;
            for jj=1:m
                subplot(ceil(m/2),2,jj); plot(time,y(:,jj),'r','linewidth',1); hold all; plot(time, y_noi(:,jj),'k--','linewidth',0.5); plot(time,noise(:,jj),'b'); legend('y','y + noise','noise');
                title(['Noise = ',num2str(RMS_noise(c_noise_level,jj)*10),'% g RMS'],'fontsize',16);
            end
            goal_noise=RMS_noise(c_noise_level,:).^2;            % true covariance of the measurement noise
            display(['True covariance of the measurement noise : ', num2str(goal_noise)]);            
            %%
            for c_r=1:N_r
                r_initial=RMS_0(c_r,:).^2;  % diagonal entries of measurement noise covariance
                r=r_initial'; % column vector
                if strcmp(cumulative,'yes') && strcmp(adaptive,'yes') && step_update(c_steps)>1
                    display('Cummulative and adaptive not supported');
                    break;
                end
                if strcmp(cumulative,'yes')
                    r=diag(kron(r_initial,ones(1,step_update(c_steps))));
                end
                P0_r=coef_r*r_initial;                          % Initial std dev of r
                Pr_initial=diag(P0_r.^2);                       % Initial covariance matrix of the estimate of diagonal entries of R
                Pr=Pr_initial;
                for c_x_ini=1:N_x_ini
                    for c_cov_ini=1:N_cov_ini
                        %% UKF ALGORITHM
                            x_initial=x_ini(c_x_ini,:)';    % To be saved (normalized)
                            x=x_initial;
                            P_initial=diag(cov(c_cov_ini,:).^2);           % Normalized initial state covariance (P0)
                            P=P_initial;
                            % Initial estimate of response
                            [y_initial]=motion(iGMfile,iGMdirection,iGMfact,N,x_initial'.*exact_p,outputs,col_output,Fs,'Initial',NstepG);
                            for kk=step_update(c_steps):step_update(c_steps):N_max_update
                               if strcmp(cumulative,'yes')
                                   YY=reshape(y_noi(kk-step_update(c_steps)+1:kk,:),step_update(c_steps)*m,1);
                               else
                                   YY=y_noi(kk,:)';
                               end
                            % TO RUN MANUAL: xhat_k_1=x; Pk_1=P; yk=YY; rhat_k_1=r; Pk_1_rr=Pr;
                            [x,P,yhat_k,r,Pr]=dukf_fem_normalized(iGMfile,iGMdirection,iGMfact,f,x,P,YY,Q,r,Pr,U,T,alpha,kappa,beta,min_lim,max_lim,new_sp,kk,workdir,col_output,exact_p,Fs,NstepG,adaptive,cumulative,step_update(c_steps));
                            % x=xhat_k; P=Pk ; r=rhat_k; Pr=Pk_rr;
                            display(['Progress: ',num2str(kk/N_max_update*100),'%']);
                              % Save variables
                              x_est(:,kk) = x;                     % Posterior state estimate
                              P_est(:,:,kk) = P;                   % Posterior state covariance
                              if strcmp(cumulative,'yes')
                                y_est(:,kk)=yhat_k;
                              else
                                y_est(1:m,kk)=yhat_k;
                                Pr_est(:,:,kk)=Pr;
                                r_est(:,kk)=r;
                              end
                              x
                            end
                            toc;
                            % Final estimate of response
                            x_final=x; % modify if another criteria (not last estimate) is assumed
                            [y_final]=motion(iGMfile,iGMdirection,iGMfact,N,x_final'.*exact_p,outputs,col_output,Fs,'Final',NstepG);
                            %% FIGURES
                            figure; % Parameters
                            for oo=1:n
                                  subplot(ceil(n/2),2,oo); plot([0 max(time)],[1 1],'k','linewidth',1.5); hold all; plot(time(1:step_update(c_steps):N_max_update), nonzeros(x_est(oo,:)),'r--','linewidth',1.5); ylabel(params{oo});
                            end
                            figure;
                            for oo=1:n
                                subplot(ceil(n/2),2,oo);   plot(time(1:step_update(c_steps):N_max_update), nonzeros(squeeze(P_est(oo,oo,:)))); ylabel(params{oo});
                            end
                            figure; 
                            % RMS errors
                            for pp=1:m
                                RMS_error_ini(pp)=rms_error(y(:,pp),y_initial(:,pp));
                            end
                            for pp=1:m
                                RMS_error_fin(pp)=rms_error(y(:,pp),y_final(:,pp));
                            end                            
                            %% SAVE RESULTS
                            tag_mat=num2str(counter);                        % Tag for name of output mat file
                            nombre_file=[name,'.mat'];
                            RMS_noi=RMS_noise(c_noise_level,:);
                            step_upd=step_update(c_steps);
                            r_ini_for_R=RMS_0(c_r,:);
                            time_process=toc(t_ini);
                            cpu_time=cputime-ti_cpu;
                            save(nombre_file,'min_lim','max_lim','iGMfile','iGMfact','iGMdirection','Fs','NstepG','params','exact_p','outputs','col_output','seed_w','alpha','kappa','beta','new_sp','n','m','N','time','y','N_max_update','y_noi','noise','y_initial','x_est','P_est','y_est','r_est','Pr_est','time_process','cpu_time','step_upd','q','Q','RMS_noi','goal_noise','adaptive','cumulative','r_ini_for_R','r_initial','coef_r','U','T','x_initial','P_initial','x_final','y_final','RMS_error_ini','RMS_error_fin','-mat');
                            counter=counter+1;
%                             x_est = zeros(n,N);                     % Parameter estimate 
%                             P_est = zeros(n,n,N);                   % State covariance
%                             y_est = zeros(m,N);                     % Output estimate
                    end
                end
            end
        end
    end
end



% profile viewer