%==========================================================================
%                     Long Term Refinancing Operation
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
%                            Load model
%=========================================================================

tic

load model_solution_mean

load policies_ltro

%=========================================================================
%                        Get Policy Functions
%=========================================================================

policies_noltro    = model_policies(param,gamma,bounds,EXP,TT,coll_points,...
                     ss,s,V,rec);
                                 
%=========================================================================
%     Select states consistent with definition of financial recession
%=========================================================================

disp('                                                                  ');    
disp('   Simulate model and select states under ''financial recession'' ');    
disp('                                                                  ');    


STATE             = zeros(6,1);
l                 = -15*ones(400,1);
Nsim              = 30000;

counter           = 0;
initial           = zeros(6,1);

for k=1:Nsim

   counter = counter+1;

  randn('state',k^(2)); 
  
e           = randn(400,4);

[state,obs,state_transformed] = simul(param,gamma',bounds,policies_noltro,...
                                TT,ss,s,V,e,l,initial);

for j=3:400-1
    if obs.gdp_growth(j)<mean(obs.gdp_growth(3:end))-2.5*var(obs.gdp_growth(3:end))^(1/2) &&...
            obs.premia(j)>mean(obs.premia(3:end))+2.5*var(obs.premia(3:end)...
            )^(1/2) && max(state_transformed(:,j)<bounds(:,1))==0 && max(state_transformed(:,j)>bounds(:,2))==0
        
        STATE = [STATE,state_transformed(:,j)];         
    end
end

if counter == 1000
      fprintf('Iteration: %4.0f of %4.0f\n', k,Nsim)
counter = 0;
end


end

STATE = STATE(:,2:end);

%=========================================================================
%    Add states from p(S^{2011:Q4}|Theta,data) and define objects
%=========================================================================

disp('                                                                  ');    
disp('       Incorporate states from p(S^{2011:Q4}|Theta,data)          ');    
disp('                                                                  ');    

load Matfiles/state_italy

filt_states   = squeeze(filtered_states(:,:,end));
M1            = size(filt_states,2);
M2            = size(STATE,2)/10;

step          = floor(M1/M2);

for j=1:M2
    
    STATE = [STATE,filt_states(:,step*j)];
    
end

M2              = M2*10;
M               = size(STATE,2); % Number of states
N               = 100;           % Number of simulations
Nsim            = 1000;          % Number of simulations for risk premia

gdp_pol         = zeros(M,N);
gdp_nopol       = zeros(M,N);
premia_pol      = zeros(M,N);
premia_nopol    = zeros(M,N);
risk_pol        = zeros(M,N);
risk_nopol      = zeros(M,N);
liqui_pol       = zeros(M,N);
liqui_nopol     = zeros(M,N);
kret            = zeros(M,Nsim);
sdf             = zeros(M,Nsim);
risk_premia     = zeros(M,1);
delta           = zeros(M,1);

m1              = zeros(2,1);
m1(2)           = 1;
m               = m1;

%=========================================================================
%        Generate paths for each state with and without LTROs
%=========================================================================

rhoz       = param(6);   
sigmaz     = param(7);   
g_star     = param(15);
rhog       = param(16);
sigmag     = param(17);
rhos       = param(22);
sigmas     = param(23);

disp('                                                                  ');    
disp('       Simulate impact effect with and without LTROs              ');    
disp('                                                                  ');    

counter    = 0;
e1         = randn(3,Nsim);
l1         = logistic(Nsim,1,0,1);

for k=1:M
  counter = counter+1;  
  initial = STATE(:,k); 
  
parfor j=1:Nsim
  
            e2       = [zeros(2,3);e1(:,j)'];        
            l2       = [0;0;l1(j)];
            
[st,obs]      = simul(param,gamma',bounds,policies_noltro,TT,ss,s,V,e2,l2,initial);
kret(k,j)     = obs.ret_cap(2:end)';
sdf(k,j)      = obs.sdf(2:end)';

end
  
parfor j=1:N

randn('state',j^(2))
e                       = randn(2,3);
l                       = -10*ones(N_pol+1,1);

[state_pol,obs_pol]     = simul_ltro(param,gamma',bounds,policies,TT,ss,s,V,e,l,initial,m);
[state_nopol,obs_nopol] = simul(param,gamma',bounds,policies_noltro,TT,ss,s,V,e,l,initial);

gdp_pol(k,j)            = obs_pol.gdp;
gdp_nopol(k,j)          = obs_nopol.gdp;
premia_pol(k,j)         = obs_pol.premia;
premia_nopol(k,j)       = obs_nopol.premia;
risk_pol(k,j)           = obs_pol.risk;
risk_nopol(k,j)         = obs_nopol.risk;
liqui_pol(k,j)          = obs_pol.liqui;
liqui_nopol(k,j)        = obs_nopol.liqui;

end
    
if counter == 100
      fprintf('Iteration: %4.0f of %4.0f\n', k,M)
counter = 0;
end
end

%=========================================================================
%                   Define effects and save results
%=========================================================================

for k=1:M   
    A = cov(kret(k,:),sdf(k,:)/mean(sdf(k,:)));
    risk_premia(k)   = -A(1,2);
    delta(k)         = risk_premia(k)/(risk_premia(k)+mean(liqui_nopol(k,:)));
end
    
impact_gdp   = mean(gdp_pol-gdp_nopol,2);
impact_ret   = mean(premia_pol-premia_nopol,2);
impact_risk  = mean(-risk_pol+risk_nopol,2);
impact_liqui = mean(liqui_pol-liqui_nopol,2);

time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))

save Matfiles/ltro delta impact_gdp impact_ret impact_risk impact_liqui STATE M2

