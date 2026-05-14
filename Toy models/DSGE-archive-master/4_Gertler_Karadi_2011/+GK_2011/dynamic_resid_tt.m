function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 11);

T(1) = params(23)*(y(13)*y(4)*y(42))^params(6);
T(2) = y(11)^(1-params(6));
T(3) = params(12)*y(11)^(1/params(5));
T(4) = (y(27)+params(20))/(params(20)+y(5))-1;
T(5) = params(13)*(y(27)+params(20))/(params(20)+y(5));
T(6) = ((params(20)+y(50))/(y(27)+params(20)))^2;
T(7) = params(13)*y(32)*T(6);
T(8) = (params(20)+y(50))/(y(27)+params(20))-1;
T(9) = y(49)/params(10);
T(10) = T(9)*y(32)*params(14)*(T(9)-1);
T(11) = params(8)-1+params(14)*y(25)/params(10)*(y(25)/params(10)-1)-T(10)*y(47)/y(22);

end
