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
    T = GK_2011.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(35, 52);
g1(1,9)=(-1);
g1(1,22)=1;
g1(1,26)=(-1);
g1(2,11)=(-(T(1)*getPowerDeriv(y(11),1-params(6),1)));
g1(2,13)=(-(T(2)*params(23)*y(4)*y(42)*T(14)));
g1(2,4)=(-(T(2)*params(23)*T(14)*y(13)*y(42)));
g1(2,22)=1;
g1(2,42)=(-(T(2)*params(23)*y(13)*y(4)*T(14)));
g1(3,10)=1;
g1(3,11)=(-(params(12)*getPowerDeriv(y(11),1/params(5),1)/y(31)));
g1(3,31)=(-((-T(3))/(y(31)*y(31))));
g1(4,10)=1;
g1(4,11)=(-((1-params(6))*y(24)*(-y(22))/(y(11)*y(11))));
g1(4,22)=(-((1-params(6))*y(24)*1/y(11)));
g1(4,24)=(-((1-params(6))*y(22)/y(11)));
g1(5,48)=y(32);
g1(5,32)=y(48);
g1(6,23)=1;
g1(6,25)=(-((-(1+y(8)))/(y(25)*y(25))));
g1(6,8)=(-(1/y(25)));
g1(7,4)=(-y(42));
g1(7,20)=1;
g1(7,27)=(-1);
g1(7,42)=(-y(4));
g1(8,14)=y(42)*y(20);
g1(8,20)=y(42)*y(14);
g1(8,26)=(-1);
g1(8,27)=1;
g1(8,42)=y(20)*y(14);
g1(9,13)=(-(y(22)*params(6)*y(24)))/(y(13)*y(13))-y(4)*y(42)*params(22);
g1(9,4)=(-(y(42)*(params(21)+params(22)*(y(13)-1))));
g1(9,22)=params(6)*y(24)/y(13);
g1(9,24)=y(22)*params(6)/y(13);
g1(9,42)=(-(y(4)*(params(21)+params(22)*(y(13)-1))));
g1(10,13)=(-(params(21)+params(22)/2*2*(y(13)-1)));
g1(10,14)=1;
g1(11,12)=1;
g1(11,14)=(-(y(42)*(-y(28))/y(6)));
g1(11,4)=(-(y(42)*(-(y(42)*y(22)*params(6)*y(24)))/(y(4)*y(42)*y(4)*y(42))/y(6)));
g1(11,22)=(-(y(42)*params(6)*y(24)/(y(4)*y(42))/y(6)));
g1(11,24)=(-(y(42)*y(22)*params(6)/(y(4)*y(42))/y(6)));
g1(11,6)=(-((-(y(42)*(y(22)*params(6)*y(24)/(y(4)*y(42))+(1-y(14))*y(28))))/(y(6)*y(6))));
g1(11,28)=(-(y(42)*(1-y(14))/y(6)));
g1(11,42)=(-((y(22)*params(6)*y(24)/(y(4)*y(42))+(1-y(14))*y(28)+y(42)*(-(y(4)*y(22)*params(6)*y(24)))/(y(4)*y(42)*y(4)*y(42)))/y(6)));
g1(12,5)=(-(params(13)/2*T(16)*2*T(4)+T(5)*T(16)+T(4)*params(13)*T(16)));
g1(12,27)=(-(params(13)/2*2*T(4)*1/(params(20)+y(5))+T(5)*1/(params(20)+y(5))+T(4)*params(13)*1/(params(20)+y(5))-(T(8)*params(13)*y(32)*(-(params(20)+y(50)))/((y(27)+params(20))*(y(27)+params(20)))*2*(params(20)+y(50))/(y(27)+params(20))+T(7)*(-(params(20)+y(50)))/((y(27)+params(20))*(y(27)+params(20))))));
g1(12,50)=T(8)*params(13)*y(32)*2*(params(20)+y(50))/(y(27)+params(20))*1/(y(27)+params(20))+T(7)*1/(y(27)+params(20));
g1(12,28)=1;
g1(12,32)=T(8)*params(13)*T(6);
g1(13,22)=(-((-(params(8)*(-(T(10)*(-y(47))/(y(22)*y(22))))))/(T(11)*T(11))));
g1(13,47)=(-((-(params(8)*(-(T(10)*1/y(22)))))/(T(11)*T(11))));
g1(13,24)=(-1)/(y(24)*y(24));
g1(13,25)=(-((-(params(8)*((y(25)/params(10)-1)*params(14)*T(15)+params(14)*y(25)/params(10)*T(15))))/(T(11)*T(11))));
g1(13,49)=(-((-(params(8)*(-(y(47)/y(22)*(y(32)*params(14)*(T(9)-1)*T(15)+T(9)*y(32)*params(14)*T(15))))))/(T(11)*T(11))));
g1(13,32)=(-((-(params(8)*(-(y(47)/y(22)*T(9)*params(14)*(T(9)-1)))))/(T(11)*T(11))));
g1(14,24)=(-((1-params(16))*params(18)*(-((-1)/(y(24)*y(24))))));
g1(14,25)=(-((1-params(16))*params(17)));
g1(14,8)=(-params(16));
g1(14,43)=1;
g1(15,45)=(-(y(32)*y(46)));
g1(15,15)=params(7);
g1(15,46)=(-(y(32)*(y(45)-y(48))));
g1(15,48)=y(32)*y(46);
g1(15,32)=(-(y(46)*(y(45)-y(48))));
g1(16,16)=1;
g1(16,17)=(-params(3));
g1(17,15)=(-(y(48)*y(32)*y(46)/((1-y(15))*(1-y(15)))));
g1(17,46)=(-(y(32)*y(48)/(1-y(15))));
g1(17,17)=1;
g1(17,48)=(-(y(32)*y(46)/(1-y(15))));
g1(17,32)=(-(y(48)*y(46)/(1-y(15))));
g1(18,12)=(-(params(3)*y(6)*y(3)));
g1(18,2)=(-(y(23)*params(3)));
g1(18,18)=1;
g1(18,3)=(-(params(3)*y(6)*(y(12)-y(23))+y(6)*params(9)));
g1(18,23)=(-(params(3)*(y(2)+y(6)*(-y(3)))));
g1(18,6)=(-(y(3)*params(9)+params(3)*(y(12)-y(23))*y(3)));
g1(19,17)=y(18);
g1(19,18)=y(17);
g1(19,19)=(-(y(28)*params(7)));
g1(19,28)=(-(params(7)*y(19)));
g1(20,18)=1;
g1(20,19)=(-y(28));
g1(20,21)=1;
g1(20,28)=(-y(19));
g1(21,19)=(-1);
g1(21,20)=1;
g1(22,18)=(-((-(y(28)*y(19)))/(y(18)*y(18))));
g1(22,19)=(-(y(28)/y(18)));
g1(22,28)=(-(y(19)/y(18)));
g1(22,30)=1;
g1(23,31)=(-((-(params(1)*y(51)))/(y(31)*y(31))));
g1(23,51)=(-(params(1)/y(31)));
g1(23,32)=1;
g1(24,46)=(-(params(1)*y(51)/y(31)));
g1(24,31)=(-(y(46)*(-(params(1)*y(51)))/(y(31)*y(31))));
g1(24,51)=(-(y(46)*params(1)/y(31)));
g1(24,33)=1;
g1(25,1)=(-((-params(15))*T(12)));
g1(25,9)=(-(T(12)-params(1)*params(15)*(-params(15))*T(13)));
g1(25,44)=params(1)*params(15)*T(13);
g1(25,31)=1;
g1(26,45)=(-1);
g1(26,23)=1;
g1(26,29)=1;
g1(27,7)=(-(params(24)*1/y(7)));
g1(27,42)=1/y(42);
g1(27,52)=1;
g1(28,34)=1;
g1(28,42)=(-(1/y(42)));
g1(29,22)=(-(1/y(22)));
g1(29,35)=1;
g1(30,9)=(-(1/y(9)));
g1(30,36)=1;
g1(31,26)=(-(1/y(26)));
g1(31,37)=1;
g1(32,20)=(-(1/y(20)));
g1(32,38)=1;
g1(33,11)=(-(1/y(11)));
g1(33,39)=1;
g1(34,28)=(-(1/y(28)));
g1(34,40)=1;
g1(35,18)=(-(1/y(18)));
g1(35,41)=1;

end
