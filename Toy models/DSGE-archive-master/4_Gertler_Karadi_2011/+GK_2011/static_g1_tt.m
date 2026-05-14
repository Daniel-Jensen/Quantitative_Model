function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 9);

T = GK_2011.static_resid_tt(T, y, x, params);

T(6) = (1-params(15))*getPowerDeriv(y(1)-y(1)*params(15),(-params(4)),1);
T(7) = getPowerDeriv(y(5)*y(12)*y(34),params(6),1);
T(8) = 1/params(10);
T(9) = (params(8)-1+params(14)*T(4)*(T(4)-1)-T(4)*(T(4)-1)*y(24)*params(14))*(params(8)-1+params(14)*T(4)*(T(4)-1)-T(4)*(T(4)-1)*y(24)*params(14));

end
