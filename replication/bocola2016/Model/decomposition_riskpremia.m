%==========================================================================
%                        Decomposing risk premia
%                         
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
%                           Load model solution
%=========================================================================

load model_solution_mean

%=========================================================================
%                          Get policy functions
%=========================================================================

[policies,expectations,residuals]    = model_policies(param,gamma,bounds,...
                                       EXP,TT,coll_points,ss,s,V,rec);

%=========================================================================
%                  Find s-shock such that p^{d}_t = 0.025
%=========================================================================

rhos     = param(22);
sigmas   = param(23);
s_star   = param(24);

shock    = 0.025;        

f        = @(x)  shock-(((exp(s_star+rhos*x)./(1+exp(s_star+rhos*x)))));

shock    = fzero(f,1);

%=========================================================================
%               Generate decomposition by simulation
%=========================================================================

% 1) Initialize state variables

e           = zeros(2000,3);
l           = -10*ones(2000,1);
initial     = zeros(6,1);
[state,obs,STATE] = simul(param,gamma',bounds,policies,TT,ss,s,V,e,l,initial);
initial     = STATE(:,end-1);
initial(4)  = 0;
initial(6)  = 0;

% 2) Draw structural shocks

randn('state',10);   % Set seed for simulations

M           = 10000;      % Number of simulations
e           = randn(3,M);
l           = logistic(1,M,0,1);

% 3) Simulate R_{K,t+1} and Lambda_{t,t+1} (See footnote to Table 4)

kret_nos    = zeros(M,1);
sdf_nos     = zeros(M,1);
kret_s      = zeros(M,1);
sdf_s       = zeros(M,1);

parfor m=1:M  
    
            e1       = [zeros(2,3);e(:,m)'];
        
            l1       = [0;0;l(m)];

% 3.a) Simulate High s_{t}

[state,obs] = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,[initial(1:5);shock]);

kret_s(m)   = obs.ret_cap(2:end)';
sdf_s(m)    = obs.sdf(2:end)';

% 3.b) Simulate Low s_{t}

[state,obs] = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,initial);

kret_nos(m) = obs.ret_cap(2:end)';
sdf_nos(m)  = obs.sdf(2:end)';

end

%=========================================================================
%                           Save results
%=========================================================================

save Matfiles/decomposition kret_s sdf_s kret_nos sdf_nos

