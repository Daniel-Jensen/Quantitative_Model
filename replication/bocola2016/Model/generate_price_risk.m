%==========================================================================
%                         Generate Price of Risk
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

tic

load model_solution_mean

rhoz       = param(6);   
rhog       = param(16);
rhos       = param(22);
sigmas     = param(23);
s_star     = param(24);

%=========================================================================
%           Load filtered states of Italian economy in 2011:Q4 
%=========================================================================

load state_italy

initial_ita    = median(squeeze(filtered_states(:,:,end)),2);
initial_ita(2) = initial_ita(2)/rhoz;
initial_ita(4) = initial_ita(4)/rhog;
initial_ita(6) = initial_ita(6)/rhos;

%=========================================================================
%                          Get policy functions
%=========================================================================

[policies,expectations,residuals]    = model_policies(param,gamma,bounds,...
                                       EXP,TT,coll_points,ss,s,V,rec);

%=========================================================================
%                  Find s-shock such that p^{d}_t = 0.05
%=========================================================================

shock    = 0.05;
f        = @(x)  shock-(((exp(s_star+rhos*x)./(1+exp(s_star+rhos*x)))));
shock    = fzero(f,1);

%=========================================================================
%               Generate decomposition by simulation
%=========================================================================

disp('                                                                  ');    
disp('                     Start Simulations                            ');    
disp('                                                                  ');    

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

M           = 10000;
e           = randn(3,M);
l           = logistic(1,M,0,1);

% 3) Simulate R_{K,t+1}, Lambda_{t,t+1}, mu_{t+1}, alp_{t+1}, \Delta c_{t+1} 
%    (See Table A-4)

kret        = zeros(M,3);
sdf         = zeros(M,3);
mult        = zeros(M,3);
alp         = zeros(M,3);
cons        = zeros(M,3);
counter     = 499;

for m=1:M  
    
            e1       = [zeros(2,3);e(:,m)'];        
            l1       = [0;0;l(m)];

%3.a) Simulate at the Ergodic Mean

[state,obs]   = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,initial);

kret(m,1)     = obs.ret_cap(2:end)';
sdf(m,1)      = obs.sdf(2:end)';
alp(m,1)      = obs.alp(2)';
cons(m,1)     = obs.cons_growth(2);
mult(m,1)     = obs.mult(2);

% 3.b) Simulate at p^{d}_{t}=0.05

[state,obs]   = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,[initial(1:5);shock]);

kret(m,2)     = obs.ret_cap(2:end)';
sdf(m,2)      = obs.sdf(2:end)';
mult(m,2)     = obs.mult(2)';
alp(m,2)      = obs.alp(2)';
cons(m,2)     = obs.cons_growth(2);

% 3.c) Simulate at p^{d}_{t}=0.05 and low Net-worth

[state,obs]   = simul_lown(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,[initial(1:5);shock]);

kret(m,3)     = obs.ret_cap(2:end)';
sdf(m,3)      = obs.sdf(2:end)';
mult(m,3)     = obs.mult(2)';
alp(m,3)      = obs.alp(2)';
cons(m,3)     = obs.cons_growth(2);

% 3.4) Simulate at filtered states

[state,obs]   = simul(param,gamma',bounds,policies,TT,ss,s,V,e1,l1,initial_ita);

kret(m,4)     = obs.ret_cap(2:end)';
sdf(m,4)      = obs.sdf(2:end)';
mult(m,4)     = obs.mult(2)';
alp(m,4)      = obs.alp(2)';
cons(m,4)     = obs.cons_growth(2);

if counter == 499
      fprintf('Iteration: %4.0f of %4.0f\n', m,M)
counter = 0;
end

counter = counter+1;

end

%=========================================================================
%                           Compute statistics
%=========================================================================

mu1   = zeros(M,1);
mu2   = zeros(M,1);
mu3   = zeros(M,1);
mu4   = zeros(M,1);

mu1(mult(:,1)>0) = 1;
mu2(mult(:,2)>0) = 1;
mu3(mult(:,3)>0) = 1;
mu4(mult(:,4)>0) = 1;

m1      = mu1.*alp(:,1);
m1      = sort(m1,'descend');
[a,b]   = min(m1);
E_mu(1) = mean(m1(1:b-1));

m2      = mu2.*alp(:,2);
m2      = sort(m2,'descend');
[a,b]   = min(m2);
E_mu(2) = mean(m2(1:b-1));

m3      = mu3.*alp(:,3);
m3      = sort(m3,'descend');
[a,b]   = min(m3);
E_mu(3) = mean(m3(1:b-1));

m4      = mu4.*alp(:,4);
m4      = sort(m4,'descend');
[a,b]   = min(m4);
E_mu(4) = mean(m4(1:b-1));

var_c      = var(cons).^(1/2);
var_alp    = var(alp).^(1/2);
var_sdf(1) = 400*var(sdf(:,1)./mean(sdf(:,1)));
var_sdf(2) = 400*var(sdf(:,2)./mean(sdf(:,2)));
var_sdf(3) = 400*var(sdf(:,3)./mean(sdf(:,3)));
var_sdf(4) = 400*var(sdf(:,4)./mean(sdf(:,4)));

time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))

%=========================================================================
%                           Save results
%=========================================================================

save Matfiles/price_risk var_sdf E_mu var_alp var_c
