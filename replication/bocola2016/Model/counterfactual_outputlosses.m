%==========================================================================
%           Counterfactual Experiment: Output Losses
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
M              = 100;
Time           = size(data(end-7:end,1),1);
Nparticle      = 500000;
Nsim           = 2000;

gdp_gs         = zeros(Time,M,Npara);
gdp_for_gs     = zeros(5,Nsim,Npara);

gdp_gns        = zeros(Time,M,Npara);
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

    if l==10
        N_particle=1000000;
    end
    
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

% Generate actual trajectories 

disp('                                                                  ');
disp('                     Actual Trajectories'                          );
disp('                                                                  ');

for k=1:M

    state    = squeeze(filt_state(:,k*steps,end-7:end));
    gdp_last = 0;

    for j=1:Time
    
    state(6,j) = state(6,j);    
    y        =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    y        =  max(y,-1);
    y        =  min(y,1);
    TTT      =  T(y,s); 

    gdp           = policies.gdp(:,1)'*TTT;
    gdp_gs(j,k,l) = gamz + state(2,j) + (gdp-gdp_last);    
    gdp_last      = gdp;
    
    x           = state(:,j);
    x(2)        = x(2)/rhoz;
    x(4)        = x(4)/rhog;
    x(6)        = x(6)/rhos;

    end
 
    
end

% Generate counterfactual trajectories 
counterfactual_states = zeros(6,M);

disp('                                                                  ');
disp('                    Counterfactual Trajectories'                   );
disp('                                                                  ');

for k=1:M
    
    state    = squeeze(filt_state(:,k*steps,end-7:end));
    gdp_last = 0;
    
    for j=1:Time
    
    state(6,j) = 0;    
    y          =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    y          =  max(y,-1);
    y          =  min(y,1);
    TTT        =  T(y,s); 

    gdp            = policies.gdp(:,1)'*TTT;
    gdp_gns(j,k,l) = gamz + state(2,j) + (gdp-gdp_last);    
    gdp_last       = gdp;
    
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

   end
   
    counterfactual_states(:,k) = state(:,j);

end

%=========================================================================
%         Compute Forecasts for output losses for 2012:Q1-2012:Q4 
%=========================================================================

disp('                                                                  ');
disp('                   Forecasts for GDP growth'                       );
disp('                                                                  ');

x_s   = mean(filt_state(:,:,end),2);
x_nos = mean(counterfactual_states,2);

randn('state',10)
e1        = randn(5,3,Nsim);
l1        = logistic(5,Nsim,0,1);

x_s(2)    = x_s(2)/rhoz;
x_s(4)    = x_s(4)/rhog;
x_s(6)    = x_s(6)/rhos;

parfor m=1:Nsim
   
  e2              = [zeros(1,3);e1(:,:,m)];  
  l2              = [0;l1(:,m)];
  [sta,obs]       = simul(param,gamma',bounds,policies,TT,ss,s,V,e2,l2,x_s);
 gdp_for_gs(:,m,l)= obs.gdp_growth(1:end);

end

x_nos(2)    = x_nos(2)/rhoz;
x_nos(4)    = x_nos(4)/rhog;
x_nos(6)    = x_nos(6)/rhos;

parfor m=1:Nsim
 
  e2                = [zeros(1,3);e1(:,:,m)];  
  l2                = [0;l1(:,m)];
  [sta,obs]         = simul(param,gamma',bounds,policies,TT,ss,s,V,e2,l2,x_nos);
 gdp_for_gns(:,m,l) = obs.gdp_growth(1:end);
  
end

l=l+1;
time=toc;

disp('                                                                  ');
fprintf('Time for iteration: %4.3f minutes\n', (time/60))
disp('                                                                  ');

%=========================================================================
%                            Save Results 
%=========================================================================

save Matfiles/output_losses gdp_gs gdp_gns gdp_for_gs gdp_for_gns 

end

