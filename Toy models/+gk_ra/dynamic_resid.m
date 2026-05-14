function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = gk_ra.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(23, 1);
    residual(1) = (T(1)) - (params(1)*(1+y(19))*T(2));
    residual(2) = (params(4)*y(10)^(1/params(3))) - (T(1)*y(11));
    residual(3) = (y(12)) - (T(4)*T(5));
    residual(4) = (y(11)) - (y(12)*(1-params(5))/y(10));
    residual(5) = (y(17)) - (y(12)*params(5)/y(1));
    residual(6) = (y(16)) - (y(14)/y(1));
    residual(7) = (y(13)) - (y(1)*(1-params(6))+y(1)*T(6));
    residual(8) = (y(15)) - (1/T(7));
    residual(9) = (1+y(18)) - ((y(17)+(1-params(6))*y(15))/y(2));
    residual(10) = (y(21)) - (params(12)*(params(10)+(1-params(10))*params(11)*y(34))*(y(33)-y(19)));
    residual(11) = (y(22)) - ((1+y(19))*params(12)*(params(10)+(1-params(10))*params(11)*y(34)));
    residual(12) = (y(23)) - (y(22)/(params(11)-y(21)));
    residual(13) = (y(13)*y(15)) - (y(23)*y(24));
    residual(14) = (y(24)) - ((1-params(10))*(1+y(20))*y(5)+y(25));
    residual(15) = (y(20)) - (y(3)+y(4)*(y(18)-y(3)));
    residual(16) = (y(25)) - (y(12)*params(13));
    residual(17) = (y(26)) - ((1+y(3))*y(6)+y(28)-y(27));
    residual(18) = (y(27)) - (params(19)+params(18)*(y(6)-params(17)));
    residual(19) = (y(28)) - ((1-params(14))*params(16)+params(14)*y(7)+x(it_, 1));
    residual(20) = (log(y(29))) - (params(15)*log(y(8))+x(it_, 2));
    residual(21) = (y(12)) - (y(28)+y(9)+y(14));
    residual(22) = (y(30)) - (y(18)-y(3));
    residual(23) = (y(31)) - (y(13)*y(15)/y(24));

end
