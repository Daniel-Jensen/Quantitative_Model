function [liki,gamma,N_eff,EXP,a] = dsgeliki(Theta,gamma,bounds,coll_points,TT,s,exp_g,V,data,Nparticle)

% This function evaluates the likelihood of the model without default risk 
%
% Inputs: Theta (vector collecting the structural parameters to be estimated), 
% gamma (numerical solution of the model, previous draw), bounds (matrix 
% collecting upper and lower bounds for the endogenous variables in the 
% Smolyak grid), coll_points (matrix collecting the collocation points 
% obtained with Smolyak method), TT (matrix collecting the Chebyshev's 
% polynomials evaluated at the collocation points), s (structure used to evaluate 
% Chebyshev's polynomials for an arbitrary point in the state space), exp_g 
% (matrix collecting coefficients for integral precomputation of government
% spending shock), V (matrix that rotate the state variables) data 
% (matrix collecting the time series on gdp growth and the Lagrange multiplier), 
% Nparticle (Number of particles)
%
% Output: liki (log-likelihood function), gamma (numerical solution of the
% model at Theta), N_eff (Number of effective particles), EXP (matrix 
% collecting coefficients for integral precomputation)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

[param,ss] = model_param(Theta);

rhoz     = param(6);
sigmaz   = param(7);

EXP      = precomp_integral_filter(coll_points,bounds,rhoz,sigmaz,s,exp_g);

obj      = @(Th) residual_nodefault(param,Th,bounds,coll_points,TT,ss,s,EXP,V);

[gamma,a]  = parsolve(obj,gamma,.5);
 
if a==1

   policies  = model_nodefault_policies(param,gamma,bounds,coll_points,TT,ss,s,EXP,V);
    
[liki,N_eff] = particle(param,policies,bounds,V,s,data,Nparticle);

liki         = sum(liki);

else
    
    liki = -100000000000000; 
    N_eff = 0;
       
end


end

