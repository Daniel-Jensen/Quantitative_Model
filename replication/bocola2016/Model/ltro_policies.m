%==========================================================================
%           Longer Term Refinancing Operations: Policy Functions
%                         
% 
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 12/07/2015
%==========================================================================

%=========================================================================
%                              Housekeeping
%=========================================================================

close all
clc
delete *.asv

path('Smolyak Files',path);
path('Simulation Files',path);
path('Solution Files',path);
path('Matfiles',path);

%=========================================================================
%                         Set Model Parameters
%=========================================================================

tic
disp('                                                                  ');
disp('          SOLVE FOR POLICY FUNCTIONS UNDER LTRO'                   );
disp('                                                                  ');

state       = 6;                % number of state variables
order       = [3;3;3;3;3;3];    % order of chebyshev's polynomial

%=========================================================================
%                         Load Model Solution
%=========================================================================

load model_solution_mean

%=========================================================================
%                           Define Objects     
%=========================================================================

N_g         = size(coll_points,2);
N_pol       = 12;

risk_free   = zeros(N_g,N_pol+1,2);
consumption = zeros(N_g,N_pol+1,2);
bond_prices = zeros(N_g,N_pol+1,2);
alpha       = zeros(N_g,N_pol+1,2);
PREMIUM     = zeros(N_g,N_pol+1,2);
RISK_comp   = zeros(N_g,N_pol+1,2);
MULT_comp   = zeros(N_g,N_pol+1,2);

%=========================================================================
%        Find Policy Functions at t=0 (Period before the policy)  
%=========================================================================

consumption(:,1,1) = gamma(1:N_g)';
risk_free(:,1,1)   = gamma(N_g+1:2*N_g)';
alpha(:,1,1)       = gamma(2*N_g+1:3*N_g)';
bond_prices(:,1,1) = gamma(3*N_g+1:4*N_g)';

consumption(:,1,2) = gamma(4*N_g+1:5*N_g)';
risk_free(:,1,2)   = gamma(5*N_g+1:6*N_g)';
alpha(:,1,2)       = gamma(6*N_g+1:7*N_g)';
bond_prices(:,1,2) = gamma(7*N_g+1:8*N_g)';

[pol,expectations,residuals,prem] = model_policies(param,gamma,bounds,...
                                         EXP,TT,coll_points,ss,s,V,rec);

PREMIUM(:,1,1)      = pol.premia(:,1);
RISK_comp(:,1,1)    = pol.risk_comp(:,1);
MULT_comp(:,1,1)    = pol.mult_comp(:,1);

PREMIUM(:,1,2)      = pol.premia(:,2);
RISK_comp(:,1,2)    = pol.risk_comp(:,2);
MULT_comp(:,1,2)    = pol.mult_comp(:,2);

%=========================================================================
%              Find policy functions in the last period      
%=========================================================================

disp('                                                                  ');
disp('          Find policy functions at t=T'                            );
disp('                                                                  ');

m       = ((exp(ss(5))^(0.3))*(0.318^(1-0.3)))*0.4;

obj     = @(Theta) residual_ltro_lastperiod(param,Theta,bounds,coll_points,...
                           TT,s,ss,EXP,V,m,expectations,prem,rec);
               
gamma   = parsolve(obj,gamma,1);

consumption(:,end,1) = gamma(1:N_g)';
risk_free(:,end,1)   = gamma(N_g+1:2*N_g)';
alpha(:,end,1)       = gamma(2*N_g+1:3*N_g)';
bond_prices(:,end,1) = gamma(3*N_g+1:4*N_g)';

consumption(:,end,2) = gamma(4*N_g+1:5*N_g)';
risk_free(:,end,2)   = gamma(5*N_g+1:6*N_g)';
alpha(:,end,2)       = gamma(6*N_g+1:7*N_g)';
bond_prices(:,end,2) = gamma(7*N_g+1:8*N_g)';

[residuals1,expectations,PREM,risk,MU,prem,pol] = model_ltro_lastperiod_policies(...
    param,gamma,bounds,coll_points,TT,s,ss,EXP,V,m,expectations,prem,rec);

PREMIUM(:,end,1)      = PREM(:,1);
RISK_comp(:,end,1)    = pol.risk_comp(:,1);
MULT_comp(:,end,1)    = pol.mult_comp(:,1);

PREMIUM(:,end,2)      = PREM(:,2);
RISK_comp(:,end,2)    = pol.risk_comp(:,2);
MULT_comp(:,end,2)    = pol.mult_comp(:,2);

%=========================================================================
%               Finf policy functions for t=2,..,T-1
%=========================================================================

disp('                                                                  ');
disp('          Find policy functions for t=T-1,..,2'                    );
disp('                                                                  ');

for j=1:N_pol-2
 disp('                                                                  ');   
 fprintf('Iteration: %4.0f of %4.0f\n',j,N_pol-2)
       
    obj                = @(Theta) residual_model_ltro(param,Theta,bounds,...
                        coll_points,TT,s,ss,EXP,V,expectations,prem,rec);
    
gamma                  = parsolve(obj,gamma,.8);

consumption(:,end-j,1) = gamma(1:N_g)';
risk_free(:,end-j,1)   = gamma(N_g+1:2*N_g)';
alpha(:,end-j,1)       = gamma(2*N_g+1:3*N_g)';
bond_prices(:,end-j,1) = gamma(3*N_g+1:4*N_g)';

consumption(:,end-j,2) = gamma(4*N_g+1:5*N_g)';
risk_free(:,end-j,2)   = gamma(5*N_g+1:6*N_g)';
alpha(:,end-j,2)       = gamma(6*N_g+1:7*N_g)';
bond_prices(:,end-j,2) = gamma(7*N_g+1:8*N_g)';

[residuals,expectations,PREM,risk,MU,prem,pol] = model_ltro_policies(param,...
              gamma,bounds,coll_points,TT,s,ss,EXP,V,expectations,prem,rec);

PREMIUM(:,end-j,1)      = PREM(:,1);
RISK_comp(:,end-j,1)    = pol.risk_comp(:,1);
MULT_comp(:,end-j,1)    = pol.mult_comp(:,1);

PREMIUM(:,end-j,2)      = PREM(:,2);
RISK_comp(:,end-j,2)    = pol.risk_comp(:,2);
MULT_comp(:,end-j,2)    = pol.mult_comp(:,2);

end

%=========================================================================
%                 Finf policy functions for t = 1
%=========================================================================

disp('                                                                  ');
disp('          Find policy functions for t=1'                           );
disp('                                                                  ');

obj                  = @(Theta) residual_model_ltro_firstperiod(param,Theta,...
                       bounds,coll_points,TT,s,ss,EXP,V,expectations,prem,rec,m);
               
gamma              = parsolve(obj,gamma,.5);

consumption(:,2,1) = gamma(1:N_g)';
risk_free(:,2,1)   = gamma(N_g+1:2*N_g)';
alpha(:,2,1)       = gamma(2*N_g+1:3*N_g)';
bond_prices(:,2,1) = gamma(3*N_g+1:4*N_g)';

consumption(:,2,2) = gamma(4*N_g+1:5*N_g)';
risk_free(:,2,2)   = gamma(5*N_g+1:6*N_g)';
alpha(:,2,2)       = gamma(6*N_g+1:7*N_g)';
bond_prices(:,2,2) = gamma(7*N_g+1:8*N_g)';

[residuals,expectations1,PREM,risk,MU,reject,pol] = model_ltro_firstperiod_policies(param,...
           gamma,bounds,coll_points,TT,s,ss,EXP,V,expectations,prem,rec,m);

PREMIUM(:,2,1)      = PREM(:,1);
RISK_comp(:,2,1)    = pol.risk_comp(:,1);
MULT_comp(:,2,1)    = pol.mult_comp(:,1);

PREMIUM(:,2,2)      = PREM(:,2);
RISK_comp(:,2,2)    = pol.risk_comp(:,2);
MULT_comp(:,2,2)    = pol.mult_comp(:,2);

policies.consumption  = consumption;
policies.risk_free    = risk_free;
policies.bond_prices  = bond_prices;
policies.alpha        = alpha;
policies.premium      = PREMIUM;
policies.risk_comp    = RISK_comp;
policies.mult_comp    = MULT_comp;

time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))

%=========================================================================
%                         Save results
%=========================================================================

save Matfiles/policies_ltro policies m N_pol reject



