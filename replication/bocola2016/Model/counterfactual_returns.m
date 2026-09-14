%==========================================================================
%          Counterfactual Experiment: Excess Returns
%                         
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 12/07/2015
%==========================================================================


%=========================================================================
%                              Housekeeping
%=========================================================================

clear all
close all
clc
delete *.asv

path('Smolyak Files',path);
path('Simulation Files',path);
path('Solution Files',path);
path('Matfiles',path);

%=========================================================================
%                         Load model solution
%=========================================================================

load model_solution_posterior

%=========================================================================
%                             Load data
%=========================================================================

load data                                   % Load data for particle filter

quarters        = data(:,1);
data            = data(:,2:4);
data            = [data(:,2),data(:,1),data(:,3)];

%=========================================================================
%                         Define Obbjects
%=========================================================================

Npara          = size(gamma_est,2);
M              = 10;
Time           = size(data(end-7:end,1),1);
Nparticle      = 500000;
Nsim           = 5000;

gdp_gs         = zeros(Time,M,Npara);
liqui_s        = zeros(Time,M,Npara);

kret_s         = zeros(Time,M,Nsim,Npara);
sdf_s          = zeros(Time,M,Nsim,Npara);
R_s            = zeros(Time,M,Nsim,Npara);
gdp_for_gs     = zeros(5,Nsim,Npara);

gdp_gns        = zeros(Time,M,Npara);
liqui_ns       = zeros(Time,M,Npara);
kret_ns        = zeros(Time,M,Nsim,Npara);
sdf_ns         = zeros(Time,M,Nsim,Npara);
R_ns           = zeros(Time,M,Nsim,Npara);
gdp_for_gns    = zeros(5,Nsim,Npara);

filtered_states= zeros(6,Time,M,Npara);

%=========================================================================
%                          Start Iteration
%=========================================================================

disp('                                                                  ');
disp('                   Counterfactual Experiment'                      );
disp('                                                                  ');

l=1;

while l<=10

tic
disp('                                                                  ');
fprintf('Draw: %4.0f of %4.0f\n', l, Npara)
disp('                                                                  ');

param  = param_est(:,l);
ss     = ss_est(:,l);
gamma  = gamma_est(:,l)';
EXP    = EXP_est(:,:,l);

[policies,expectations,residuals]    = model_policies(param,gamma,bounds,...
                                       EXP,TT,coll_points,ss,s,V,haircut(l));

%=========================================================================
%                 Estimate path for model state variables
%=========================================================================

disp('                                                                  ');
disp('                   Starting Particle Filter'                       );
disp('                                                                  ');

gam = gamma';

[filt_state,observation,distr,N_eff] = particle(data,param,gam,bounds,s,V,...
                                       policies,Nparticle,quarters);

%=========================================================================
%            Generate Actual and Counterfactual Trajectories
%=========================================================================

gamz       = param(5);
rhoz       = param(6);   
sigmaz     = param(7);   
rhog       = param(16);
sigmag     = param(17);
rhos       = param(22);
sigmas     = param(23);
s_star     = param(24);

disp('                                                                  ');
disp('        Generate Actual and Counterfactual Trajectories'           );
disp('                                                                  ');

steps          = floor(Nparticle/M);

randn('state',10)
e1             = randn(3,Nsim);
l1             = logistic(1,Nsim,0,1);

% Generate actual trajectories 

disp('                                                                  ');
disp('                     Actual Trajectories'                          );
disp('                                                                  ');

for k=1:M

    fprintf('Iteration: %4.0f of %4.0f\n', k, M)

    state    = squeeze(filt_state(:,k*steps,end-7:end));

    for j=1:Time
       
    x           = state(:,j);    
    
    x(2)        = x(2)/rhoz;
    x(4)        = x(4)/rhog;
    x(6)        = x(6)/rhos;

    parfor m=1:Nsim      
            e2       = [zeros(2,3);e1(:,m)'];       
            l2       = [0;0;l1(m)];            
    [st,obs]         = simul(param,gamma',bounds,policies,TT,ss,s,V,e2,l2,x);
    kret_s(j,k,m,l)  = obs.ret_cap(2:end)';
    R_s(j,k,m,l)     = obs.R(1);
    sdf_s(j,k,m,l)   = obs.sdf(2:end)';
    end
     
    end
 
    
end

% Generate counterfactual trajectories 
counterfactual_states = zeros(6,M);

disp('                                                                  ');
disp('                    Counterfactual Trajectories'                   );
disp('                                                                  ');

for k=1:M
    
    fprintf('Iteration: %4.0f of %4.0f\n', k, M)

    state    = squeeze(filt_state(:,k*steps,end-7:end));
    
    for j=1:Time
    
    state(6,j) = 0;    
    y        =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    y        =  max(y,-1);
    y        =  min(y,1);
    TTT      =  T(y,s); 
    
    kap         = policies.kap(:,1)'*TTT;
    prom        = policies.prom(:,1)'*TTT;
    debt        = policies.debt(:,1)'*TTT;
    XX          = [kap',prom'];
    X_tilde     = XX*V;
    state(1,j+1)= X_tilde(:,1)';
    state(3,j+1)= X_tilde(:,2)';
    state(5,j+1)= debt;
    
    x           = state(:,j);
    x(2)        = x(2)/rhoz;
    x(4)        = x(4)/rhog;
    x(6)        = x(6)/rhos;

    parfor m=1:Nsim      
            e2       = [zeros(2,3);e1(:,m)'];       
            l2       = [0;0;l1(m)];            
    [st,obs]         = simul(param,gamma',bounds,policies,TT,ss,s,V,e2,l2,x);
    kret_ns(j,k,m,l) = obs.ret_cap(2:end)';
    R_ns(j,k,m,l)    = obs.R(1);
    sdf_ns(j,k,m,l)  = obs.sdf(2:end)';
    end
     
    end
   
    counterfactual_states(:,k) = state(:,j);

end

l=l+1;
time=toc;

disp('                                                                  ');
fprintf('Time for iteration: %4.3f minutes\n', (time/60))
disp('                                                                  ');

%=========================================================================
%                            Save Results 
%=========================================================================

save Matfiles/counterfactual kret_s kret_ns R_s R_ns sdf_s sdf_ns 
end

