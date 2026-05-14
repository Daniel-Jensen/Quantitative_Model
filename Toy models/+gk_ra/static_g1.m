function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = gk_ra.static_g1_tt(T, y, x, params);
end
g1 = zeros(23, 23);
g1(1,1)=T(7)-params(1)*(1+y(11))*T(7);
g1(1,11)=(-(T(1)*params(1)));
g1(2,1)=(-(y(3)*T(7)));
g1(2,2)=params(4)*getPowerDeriv(y(2),1/params(3),1);
g1(2,3)=(-T(1));
g1(3,2)=(-(T(3)*getPowerDeriv(y(2),1-params(5),1)));
g1(3,4)=1;
g1(3,5)=(-(T(4)*y(21)*getPowerDeriv(y(5),params(5),1)));
g1(3,21)=(-(T(2)*T(4)));
g1(4,2)=(-((-(y(4)*(1-params(5))))/(y(2)*y(2))));
g1(4,3)=1;
g1(4,4)=(-((1-params(5))/y(2)));
g1(5,4)=(-(params(5)/y(5)));
g1(5,5)=(-((-(y(4)*params(5)))/(y(5)*y(5))));
g1(5,9)=1;
g1(6,5)=(-((-y(6))/(y(5)*y(5))));
g1(6,6)=(-(1/y(5)));
g1(6,8)=1;
g1(7,5)=1-(1-params(6)+T(5));
g1(7,8)=(-(y(5)*params(8)*getPowerDeriv(y(8),1-params(7),1)));
g1(8,7)=1;
g1(8,8)=(-((-(params(8)*(1-params(7))*getPowerDeriv(y(8),(-params(7)),1)))/(T(6)*T(6))));
g1(9,7)=(-(((1-params(6))*y(7)-(y(9)+(1-params(6))*y(7)))/(y(7)*y(7))));
g1(9,9)=(-(1/y(7)));
g1(9,10)=1;
g1(10,10)=(-(params(12)*(params(10)+(1-params(10))*params(11)*y(15))));
g1(10,11)=params(12)*(params(10)+(1-params(10))*params(11)*y(15));
g1(10,13)=1;
g1(10,15)=(-((y(10)-y(11))*params(12)*(1-params(10))*params(11)));
g1(11,11)=(-(params(12)*(params(10)+(1-params(10))*params(11)*y(15))));
g1(11,14)=1;
g1(11,15)=(-((1+y(11))*params(12)*(1-params(10))*params(11)));
g1(12,13)=(-(y(14)/((params(11)-y(13))*(params(11)-y(13)))));
g1(12,14)=(-(1/(params(11)-y(13))));
g1(12,15)=1;
g1(13,5)=y(7);
g1(13,7)=y(5);
g1(13,15)=(-y(16));
g1(13,16)=(-y(15));
g1(14,12)=(-((1-params(10))*y(16)));
g1(14,16)=1-(1-params(10))*(1+y(12));
g1(14,17)=(-1);
g1(15,10)=(-y(15));
g1(15,11)=(-(1-y(15)));
g1(15,12)=1;
g1(15,15)=(-(y(10)-y(11)));
g1(16,4)=(-params(13));
g1(16,17)=1;
g1(17,11)=(-y(18));
g1(17,18)=1-(1+y(11));
g1(17,19)=1;
g1(17,20)=(-1);
g1(18,18)=(-params(18));
g1(18,19)=1;
g1(19,20)=1-params(14);
g1(20,21)=1/y(21)-params(15)*1/y(21);
g1(21,1)=(-1);
g1(21,4)=1;
g1(21,6)=(-1);
g1(21,20)=(-1);
g1(22,10)=(-1);
g1(22,11)=1;
g1(22,22)=1;
g1(23,5)=(-(y(7)/y(16)));
g1(23,7)=(-(y(5)/y(16)));
g1(23,16)=(-((-(y(5)*y(7)))/(y(16)*y(16))));
g1(23,23)=1;

end
