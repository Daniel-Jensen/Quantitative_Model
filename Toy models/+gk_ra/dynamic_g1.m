function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = gk_ra.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(23, 36);
g1(1,9)=T(8);
g1(1,32)=(-(params(1)*(1+y(19))*getPowerDeriv(y(32),(-params(2)),1)));
g1(1,19)=(-(params(1)*T(2)));
g1(2,9)=(-(y(11)*T(8)));
g1(2,10)=params(4)*getPowerDeriv(y(10),1/params(3),1);
g1(2,11)=(-T(1));
g1(3,10)=(-(T(4)*getPowerDeriv(y(10),1-params(5),1)));
g1(3,12)=1;
g1(3,1)=(-(T(5)*y(29)*getPowerDeriv(y(1),params(5),1)));
g1(3,29)=(-(T(3)*T(5)));
g1(4,10)=(-((-(y(12)*(1-params(5))))/(y(10)*y(10))));
g1(4,11)=1;
g1(4,12)=(-((1-params(5))/y(10)));
g1(5,12)=(-(params(5)/y(1)));
g1(5,1)=(-((-(y(12)*params(5)))/(y(1)*y(1))));
g1(5,17)=1;
g1(6,1)=(-((-y(14))/(y(1)*y(1))));
g1(6,14)=(-(1/y(1)));
g1(6,16)=1;
g1(7,1)=(-(1-params(6)+T(6)));
g1(7,13)=1;
g1(7,16)=(-(y(1)*params(8)*getPowerDeriv(y(16),1-params(7),1)));
g1(8,15)=1;
g1(8,16)=(-((-(params(8)*(1-params(7))*getPowerDeriv(y(16),(-params(7)),1)))/(T(7)*T(7))));
g1(9,2)=(-((-(y(17)+(1-params(6))*y(15)))/(y(2)*y(2))));
g1(9,15)=(-((1-params(6))/y(2)));
g1(9,17)=(-(1/y(2)));
g1(9,18)=1;
g1(10,33)=(-(params(12)*(params(10)+(1-params(10))*params(11)*y(34))));
g1(10,19)=params(12)*(params(10)+(1-params(10))*params(11)*y(34));
g1(10,21)=1;
g1(10,34)=(-((y(33)-y(19))*params(12)*(1-params(10))*params(11)));
g1(11,19)=(-(params(12)*(params(10)+(1-params(10))*params(11)*y(34))));
g1(11,22)=1;
g1(11,34)=(-((1+y(19))*params(12)*(1-params(10))*params(11)));
g1(12,21)=(-(y(22)/((params(11)-y(21))*(params(11)-y(21)))));
g1(12,22)=(-(1/(params(11)-y(21))));
g1(12,23)=1;
g1(13,13)=y(15);
g1(13,15)=y(13);
g1(13,23)=(-y(24));
g1(13,24)=(-y(23));
g1(14,20)=(-((1-params(10))*y(5)));
g1(14,5)=(-((1-params(10))*(1+y(20))));
g1(14,24)=1;
g1(14,25)=(-1);
g1(15,18)=(-y(4));
g1(15,3)=(-(1-y(4)));
g1(15,20)=1;
g1(15,4)=(-(y(18)-y(3)));
g1(16,12)=(-params(13));
g1(16,25)=1;
g1(17,3)=(-y(6));
g1(17,6)=(-(1+y(3)));
g1(17,26)=1;
g1(17,27)=1;
g1(17,28)=(-1);
g1(18,6)=(-params(18));
g1(18,27)=1;
g1(19,7)=(-params(14));
g1(19,28)=1;
g1(19,35)=(-1);
g1(20,8)=(-(params(15)*1/y(8)));
g1(20,29)=1/y(29);
g1(20,36)=(-1);
g1(21,9)=(-1);
g1(21,12)=1;
g1(21,14)=(-1);
g1(21,28)=(-1);
g1(22,18)=(-1);
g1(22,3)=1;
g1(22,30)=1;
g1(23,13)=(-(y(15)/y(24)));
g1(23,15)=(-(y(13)/y(24)));
g1(23,24)=(-((-(y(13)*y(15)))/(y(24)*y(24))));
g1(23,31)=1;

end
