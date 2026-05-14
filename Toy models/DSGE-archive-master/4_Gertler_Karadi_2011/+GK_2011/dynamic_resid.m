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
    T = GK_2011.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(35, 1);
    residual(1) = (y(22)) - (y(9)+y(26));
    residual(2) = (y(22)) - (T(1)*T(2));
    residual(3) = (y(10)) - (T(3)/y(31));
    residual(4) = (y(10)) - ((1-params(6))*y(24)*y(22)/y(11));
    residual(5) = (y(32)*y(48)) - (1);
    residual(6) = (y(23)) - ((1+y(8))/y(25));
    residual(7) = (y(20)) - (y(27)+y(4)*y(42));
    residual(8) = (y(27)) - (y(26)-y(20)*y(42)*y(14));
    residual(9) = (y(22)*params(6)*y(24)/y(13)) - (y(4)*y(42)*(params(21)+params(22)*(y(13)-1)));
    residual(10) = (y(14)) - (params(2)+params(21)*(y(13)-1)+params(22)/2*(y(13)-1)^2);
    residual(11) = (y(12)) - (y(42)*(y(22)*params(6)*y(24)/(y(4)*y(42))+(1-y(14))*y(28))/y(6));
    residual(12) = (y(28)) - (1+params(13)/2*T(4)^2+T(4)*T(5)-T(7)*T(8));
    residual(13) = (1/y(24)) - (params(8)/T(11));
    residual(14) = (y(43)) - ((1-params(16))*(params(11)+params(17)*(y(25)-params(10))+params(18)*(params(8)/(params(8)-1)-1/y(24)))+y(8)*params(16));
    residual(15) = (params(7)*y(15)) - (y(32)*y(46)*(y(45)-y(48)));
    residual(16) = (y(16)) - (1-params(3)+params(3)*y(17));
    residual(17) = (y(17)) - (y(48)*y(32)*y(46)/(1-y(15)));
    residual(18) = (y(18)) - (params(3)*(y(6)*(y(12)-y(23))*y(3)+y(23)*y(2))+y(6)*y(3)*params(9));
    residual(19) = (y(17)*y(18)) - (y(28)*params(7)*y(19));
    residual(20) = (y(18)+y(21)) - (y(28)*y(19));
    residual(21) = (y(20)) - (y(19));
    residual(22) = (y(30)) - (y(28)*y(19)/y(18));
    residual(23) = (y(32)) - (params(1)*y(51)/y(31));
    residual(24) = (y(33)) - (y(46)*params(1)*y(51)/y(31));
    residual(25) = (y(31)) - ((y(9)-params(15)*y(1))^(-params(4))-params(1)*params(15)*(y(44)-y(9)*params(15))^(-params(4)));
    residual(26) = (y(29)) - (y(45)-y(23));
    residual(27) = (log(y(42))) - (params(24)*log(y(7))-x(it_, 1));
    residual(28) = (y(34)) - (log(y(42)));
    residual(29) = (y(35)) - (log(y(22)));
    residual(30) = (y(36)) - (log(y(9)));
    residual(31) = (y(37)) - (log(y(26)));
    residual(32) = (y(38)) - (log(y(20)));
    residual(33) = (y(39)) - (log(y(11)));
    residual(34) = (y(40)) - (log(y(28)));
    residual(35) = (y(41)) - (log(y(18)));

end
