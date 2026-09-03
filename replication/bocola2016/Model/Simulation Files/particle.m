function [filt_state,observation,distr,N_eff] = particle(data,param,gam,bounds,s,V,policies,Nparticle,quarters)

% This function runs the particle filter to the model and obtains estimates 
% for the path of model's state variables. See Section 5.2 and Appendix B
% for a description of the procedure.
%
% Inputs: data (matrix collecting the series in the measurement equations),
% param (vector collecting structural parameters), gamma (numerical solution 
% of the model), bounds (matrix collecting upper and
% lower bounds for the endogenous variables in the Smolyak grid),
% TT (matrix collecting the Chebyshev's polynomials evaluated at the 
% collocation points), s (structure used to evaluate
% Chebyshev's polynomials for an arbitrary point in the state space), ss 
% (deterministic steady state), V (matrix that rotate the state variables)
% policies (structure collecting the model's policy functuions), Nparticle
% (number of particles), quarters (vector collecting quarters in sample)
%
% Output: filt_state (filtered state variables), observation (filtered
% observables), distr (probabilities assigned to each particle), N_eff 
% (number of effective particles), liki (log-likelihood function) 
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

options=optimset('Display','off','MaxFunEvals',300000,'MaxIter',10000,'LargeScale','off');

N_g        = size(gam,1)/(2*4);
gam        = gam(1:4*N_g);

gamz       = param(5);
rhoz       = param(6);   
sigmaz     = param(7);   
rhog       = param(16);
sigmag     = param(17);
rhos       = param(22);
sigmas     = param(23);
s_star     = param(24);

% Define Objects
 
Time         = size(data,1);
N_eff        = zeros(Time,1);
distr        = zeros(Nparticle,Time);
observation  = zeros(3,Nparticle,Time);
filt_state   = zeros(6,Nparticle,Time);
Sigma        = inv(diag(0.05*var(data)));
weights1     = ones(1,Nparticle);
weights2     = ones(1,Nparticle);

% Initialize the state variables

initial_state      = zeros(6,1)*ones(1,Nparticle);
initial_state(6,:) = s_star;
gdp_last           = zeros(1,Nparticle);
e_cent             = zeros(Time,3);
t                  = 1;
counter            = 1;

while t<=Time

    f         = @(ee) proposal_distr(data(t,:)',bounds,param,[0,0,0;ee],...
                mean(initial_state,2),mean(gdp_last,2),policies,Sigma,s);

  e_cent(t,:) = fminunc(f,zeros(1,3),options);
     
randn('state',t^(2))

inn           = (ones(Nparticle,1)*squeeze(e_cent(t,:)))'+randn(3,Nparticle);

state         = [initial_state(1,:);rhoz*initial_state(2,:)+sigmaz*inn(1,:);...
                initial_state(3,:);rhog*initial_state(4,:)+sigmag*inn(2,:);...
                initial_state(5,:);rhos*initial_state(6,:)+sigmas*inn(3,:)];

filt_state(:,:,t)  = state;

y             = (2*state-(bounds(:,1)+bounds(:,2))*ones(1,Nparticle))./((bounds(:,2)-...
                 bounds(:,1))*ones(1,Nparticle));

y             = min(y,1);
y             = max(y,-1);
TTT           = T(y,s);

kap           = policies.kap(:,1)'*TTT;
prom          = policies.prom(:,1)'*TTT;
debt          = policies.debt(:,1)'*TTT;
gdp           = policies.gdp(:,1)'*TTT;
gdp_growth    = gamz+state(2,:)+gdp-gdp_last;
mult          = max(policies.mult(:,1)'*TTT,0);
def_prob      = exp(state(6,:)+s_star)./(1+exp(state(6,:)+s_star));

   obs        = [gdp_growth;mult;def_prob];
   XX         = [kap',prom'];
   X_tilde    = XX*V;
   state(1,:) = X_tilde(:,1)';
   state(3,:) = X_tilde(:,2)';
   state(5,:) = debt;

    v         = obs'-ones(Nparticle,1)*data(t,:);

parfor j=1:Nparticle 
    
weights1(j)       = (((1/(2*pi)^(2))^(1/2))*det(inv(Sigma))^(-1/2)*exp(-(1/2)*...
                    v(j,:)*Sigma*v(j,:)'));
    
weights2(j)       = weights1(j)*(exp(-(1/2)*inn(:,j)'*inn(:,j))...
                    /exp(-(1/2)*(inn(:,j)-(e_cent(t,:)'))'*(inn(:,j)-(e_cent(t,:)'))));    

end

weights               = weights2./(sum(weights2));
initial_state         = resmpl(state,weights);
filt_state(:,:,t)     = resmpl(filt_state(:,:,t),weights);
N_eff(t)              = 1/sum(weights.^2);
observation(:,:,t)    = obs;
gdp_last              = resmpl(gdp,weights);
distr(:,t)            = weights';

     
 if   counter == 5 || t==Time
    disp(['                      Period:',num2str(quarters(t))]); 
    fprintf('Effective number of Particles        %4.3f\n', N_eff(t))
    disp('                                                                  ');       
    fprintf('Observable        Data          Model\n')
    fprintf('----------        ----          -----\n')
    fprintf('GDP Growth        %4.3f         %4.3f\n', data(t,1)*100,100*sum(squeeze(observation(1,:,t))'.*distr(:,t)))
    fprintf('Multiplier        %4.3f         %4.3f\n', data(t,2),sum(squeeze(observation(2,:,t))'.*distr(:,t)))
    fprintf('Default Prob      %4.3f         %4.3f\n', data(t,3)*100,100*sum(squeeze(observation(3,:,t))'.*distr(:,t)))
    disp('                                                                  ');        
    counter = 0;
 end
 
    t = t+1;
    counter = counter+1;

end

end

