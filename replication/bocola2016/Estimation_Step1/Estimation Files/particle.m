function [liki,N_eff,filt_state,observation,distr] = particle(param,policies,bounds,V,s,data,Nparticle)

% This function runs the particle filter for our model 
%
% Inputs: param (vector collecting the structural parameters of the model), 
% policies (structure collecting the model's policy functions), s (structure 
% used to evaluate Chebyshev's polynomials for an arbitrary point in the 
% state space), V (matrix that rotate the state variables),  data (matrix 
% collecting the time series on gdp growth and the Lagrange multiplier), 
% Nparticle (Number of particles)
%
% Output: liki (log-likelihood function), N_eff (number of effective particles)
% filt_state (filtered state variables), observation (filtered observables), 
% distr (weights for particles)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

options=optimset('Display','off','MaxFunEvals',300000,'MaxIter',10000,'LargeScale','off');

gamz     = param(5);
rhoz     = param(6);
sigmaz   = param(7);
rhog     = param(16);
sigmag   = param(17);

% Define Objects

Time          = size(data,1);
N_eff         = zeros(Time,1);
liki          = zeros(Time,1);
observation   = zeros(2,Nparticle,Time);
distr         = zeros(Nparticle,Time);
filt_state    = zeros(5,Nparticle,Time);
Sigma         = inv(diag(0.25*var(data)));
weights1      = ones(1,Nparticle);
weights2      = ones(1,Nparticle);

initial_state = zeros(5,1)*ones(1,Nparticle);
gdp_last      = zeros(1,Nparticle);
e_cent        = zeros(Time,2);

t             = 1;
 
    
while t<=Time

f             = @(ee) proposal_distr(data(t,:)',bounds,param,[0,0;ee],mean(initial_state,2),mean(gdp_last,2),policies,Sigma,s);
        
  e_cent(t,:) = fminunc(f,zeros(1,2),options);

inn           = (ones(Nparticle,1)*squeeze(e_cent(t,:)))'+randn(2,Nparticle);
  
state         = [initial_state(1,:);rhoz*initial_state(2,:)+sigmaz*inn(1,:);...
                initial_state(3,:);rhog*initial_state(4,:)+sigmag*inn(2,:);initial_state(5,:)];
                       
y             = (2*state-(bounds(:,1)+bounds(:,2))*ones(1,Nparticle))./((bounds(:,2)-...
                 bounds(:,1))*ones(1,Nparticle));

y             = min(y,1);
y             = max(y,-1);
TTT           = T(y,s);

kap           = policies.kap*TTT;
prom          = policies.prom*TTT;
bprime        = policies.debt*TTT;
gdp           = policies.gdp*TTT;
gdp_growth    = gamz+state(2,:)+gdp-gdp_last;
mult          = max(policies.mult*TTT,0);
 
filt_state(:,:,t) = state;

   obs            = [gdp_growth;mult];
   XX             = [kap',prom'];
   X_tilde        = XX*V;
   
   state(1,:)     = X_tilde(:,1)';
   state(3,:)     = X_tilde(:,2)';
   state(5,:)     = bprime;
   
v                 = obs'-ones(Nparticle,1)*data(t,:);

parfor j=1:Nparticle 
  
weights1(j)       = (((1/(2*pi)^(2))^(1/2))*det(inv(Sigma))^(-1/2)*exp(-(1/2)*...
                    v(j,:)*Sigma*v(j,:)'));
    
weights2(j)       = weights1(j)*(exp(-(1/2)*inn(:,j)'*inn(:,j))...
                    /exp(-(1/2)*(inn(:,j)-(e_cent(t,:)'))'*(inn(:,j)-(e_cent(t,:)'))));    
end

liki(t)               = log(sum(weights1)/Nparticle);    
weights               = weights2./(sum(weights2));
initial_state         = resmpl(state,weights);
N_eff(t)              = 1/sum(weights.^2);
observation(:,:,t)    = obs;
gdp_last              = resmpl(gdp,weights);
distr(:,t)            = weights';
t = t+1;

end

end

