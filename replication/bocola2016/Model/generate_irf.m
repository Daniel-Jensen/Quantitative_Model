%==========================================================================
%                    Generate Impulse Response Functions
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
%                          Load Starting Guess 
%=========================================================================

tic

load model_solution_mean

%=========================================================================
%                        Get Policy Functions
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
%                Generate Impulse Response Functions
%=========================================================================

disp('                                                                  ');    
disp('        Simulate model and initialize state variables             ');    
disp('                                                                  ');    

randn('state',10);   % Set seed for simulations

% 1) Initialize state variables

e           = zeros(2000,3);
l           = -10*ones(2000,1);
initial     = zeros(6,1);
[state,obs,STATE] = simul(param,gamma',bounds,policies,TT,ss,s,V,e,l,initial);
initial     = STATE(:,end-1);
initial(2)  = 0;
initial(4)  = 0;
initial(6)  = 0;

% 2) Draw structural shocks

M           = 10000;
T           = 30;
e           = randn(T,3,M);

% 3) Simulate endogenous variables with and without the s-shock

bond_nos    = zeros(M,T-1);
nworth_nos  = zeros(M,T-1);
rfree_nos   = zeros(M,T-1);
cons_nos    = zeros(M,T-1);
inv_nos     = zeros(M,T-1);
gdp_nos     = zeros(M,T-1);
lab_nos     = zeros(M,T-1);
lev_nos     = zeros(M,T-1);
mult_nos    = zeros(M,T-1);
risk_nos    = zeros(M,T-1);
liqui_nos   = zeros(M,T-1);
premia_nos  = zeros(M,T-1);
state_nos   = zeros(T-1,6,M);

bond_s      = zeros(M,T-1);
nworth_s    = zeros(M,T-1);
rfree_s     = zeros(M,T-1);
cons_s      = zeros(M,T-1);
inv_s       = zeros(M,T-1);
gdp_s       = zeros(M,T-1);
lab_s       = zeros(M,T-1);
lev_s       = zeros(M,T-1);
mult_s      = zeros(M,T-1);
risk_s      = zeros(M,T-1);
liqui_s     = zeros(M,T-1);
premia_s    = zeros(M,T-1);
state_s     = zeros(T-1,6,M);

disp('                                                                  ');    
disp('            Compute IRFs for selected variables                   ');    
disp('                                                                  ');    

parfor m=1:M     
    
            e_ns        = e(:,:,m);
            e_s         = e_ns;
            e_s(2,3)    = shock/sigmas;            
           l1           = -10*ones(1,T); % this guarantees that there are no
                                         % defaults on path
    
% 3.a) Simulate without the shock
           
[state,obs,STATE]     = simul(param,gamma',bounds,policies,TT,ss,s,V,e_ns,l1,initial);

bond_nos(m,:)   = obs.q(1:end)';
nworth_nos(m,:) = obs.networth(1:end)';
rfree_nos(m,:)  = obs.R(1:end)';
cons_nos(m,:)   = obs.cons(1:end)';
inv_nos(m,:)    = obs.inv(1:end)';
gdp_nos(m,:)    = obs.gdp(1:end)';
lab_nos(m,:)    = obs.lab(1:end)';
state_nos(:,:,m)= STATE(:,1:end-1)';
mult_nos(m,:)   = obs.mult(1:end)';
risk_nos(m,:)   = obs.risk(1:end)';
lev_nos(m,:)    = obs.lev(1:end)';
liqui_nos(m,:)  = obs.liqui(1:end)';
premia_nos(m,:) = obs.premia(1:end)';

%3.b) Simulate with the shock

[state,obs,STATE] = simul(param,gamma',bounds,policies,TT,ss,s,V,e_s,l1,initial);

bond_s(m,:)     = obs.q(1:end)';
nworth_s(m,:)   = obs.networth(1:end)';
rfree_s(m,:)    = obs.R(1:end)';
cons_s(m,:)     = obs.cons(1:end)';
inv_s(m,:)      = obs.inv(1:end)';
gdp_s(m,:)      = obs.gdp(1:end)';
lab_s(m,:)      = obs.lab(1:end)';
state_s(:,:,m)  = STATE(:,1:end-1)';
mult_s(m,:)     = obs.mult(1:end)';
risk_s(m,:)     = obs.risk(1:end)';
lev_s(m,:)      = obs.lev(1:end)';
liqui_s(m,:)    = obs.liqui(1:end)';
premia_s(m,:)   = obs.premia(1:end)';

end

disp('                                                                  ');    
disp('                    Compute Risk Premia                           ');    
disp('                                                                  ');    


% 3.c) Simulate R_{K,t+1}, R_{t}, and Lambda_{t,t+1} with and without the
% s-shock

initial_s      = mean(squeeze(state_s(2:end,:,:)),3);
initial_nos    = mean(squeeze(state_nos(2:end,:,:)),3);

M              = 10000;
e              = randn(3,3,M);
l              = logistic(1,M,0,1);

kret_nos       = zeros(M,T-2);
sdf_nos        = zeros(M,T-2);
R_nos          = zeros(M,T-2);
kret_s         = zeros(M,T-2);
sdf_s          = zeros(M,T-2);
R_s            = zeros(M,T-2);


for j=1:T-2
 
    fprintf('Iteration: %4.0f of %4.0f\n', j,T-2)

    in_s     = initial_s(j,:)';
    in_nos   = initial_nos(j,:)';

 parfor m=1:M  
            e1          = e(:,:,m);
            e1(1:2,:)   = 0;
           l1           = [-10,-10,l(m)]; 
            
[state,obs]     = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,in_s);
kret_s(m,j)     = obs.ret_cap(2:end)';
sdf_s(m,j)      = obs.sdf(2:end)';
R_s(m,j)        = obs.R(1)';

[state,obs]     = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,in_nos);

kret_nos(m,j)   = obs.ret_cap(2:end)';
sdf_nos(m,j)    = obs.sdf(2:end)';
R_nos(m,j)      = obs.R(1)';

  end

end

%=========================================================================
%                                Save Results 
%=========================================================================

bond_price      = mean(log(bond_s),1)-mean(log(bond_nos),1);
networth        = mean(nworth_s,1)-mean(nworth_nos,1);
cons            = mean(cons_s,1)-mean(cons_nos,1);
inv             = mean(inv_s,1)-mean(inv_nos,1);
gdp             = mean(gdp_s,1)-mean(gdp_nos,1);
lab             = mean(lab_s,1)-mean(lab_nos,1);
liq_premia      = median(liqui_s-liqui_nos,1);
risk_premia     = zeros(1,T-2);

for j=1:T-2
    
    A_s   = cov(kret_s(:,j),-sdf_s(:,j)./mean(sdf_s(:,j)));
    A_nos = cov(kret_nos(:,j),-sdf_nos(:,j)./mean(sdf_nos(:,j)));
    
    risk_premia(j) = A_s(1,2)-A_nos(1,2);
    
end
    
prob_def = exp(s_star+initial_s(2:end,6))./(1+exp(s_star+initial_s(2:end,6)));

time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))


save Matfiles/irf_bench prob_def risk_premia liq_premia kret_s kret_nos...
     sdf_s sdf_nos R_s R_nos bond_price networth cons inv gdp lab 

