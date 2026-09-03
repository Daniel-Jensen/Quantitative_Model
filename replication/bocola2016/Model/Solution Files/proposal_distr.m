function [obj] = proposal_distr(data,bounds,param,e,initial,gdp_last,policies,Sigma,s)

% This routine evaluates the objective function in equation (A.20) in
% Online Appendix B.1.1. It is necessary to compute the candidate density
% for the auxiliary particle filter
%
% Inputs: data (matrix collecting the series in the measurement equations),
% param (vector collecting structural parameters), bounds (matrix collecting upper and
% lower bounds for the endogenous variables in the Smolyak grid),
% initial (value for the state variables), s (structure used to evaluate
% Chebyshev's polynomials for an arbitrary point in the state space), 
% policies (structure collecting the model's policy functuions), e (initial
% guess for the structural shocks), Sigma (variance-covariance matrix of
% measurement errors).
%
% Output: obj (value of the objective function), 
%
% Author: Luigi Bocola     lbocola@northwestern.edu
% Date  : 09/07/2015

sim_data   = simul_filter(bounds,param,e,initial,policies,s);

sim_data(1)= (sim_data(1)-gdp_last);

err        = data-sim_data;

obj        = err'*Sigma*err + e(2,:)*e(2,:)';

end

