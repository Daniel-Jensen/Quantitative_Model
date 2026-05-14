function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(7, 1);
end
[T_order, T] = gk_ra.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(23, 1);
    residual(1) = (T(1)) - (params(1)*(1+y(34))*T(2));
    residual(2) = (params(4)*y(25)^(1/params(3))) - (T(1)*y(26));
    residual(3) = (y(27)) - (T(4)*T(5));
    residual(4) = (y(26)) - (y(27)*(1-params(5))/y(25));
    residual(5) = (y(32)) - (y(27)*params(5)/y(5));
    residual(6) = (y(31)) - (y(29)/y(5));
    residual(7) = (y(28)) - (y(5)*(1-params(6))+y(5)*T(6));
    residual(8) = (y(30)) - (1/T(7));
    residual(9) = (1+y(33)) - ((y(32)+(1-params(6))*y(30))/y(7));
    residual(10) = (y(36)) - (params(12)*(params(10)+(1-params(10))*params(11)*y(61))*(y(56)-y(34)));
    residual(11) = (y(37)) - ((1+y(34))*params(12)*(params(10)+(1-params(10))*params(11)*y(61)));
    residual(12) = (y(38)) - (y(37)/(params(11)-y(36)));
    residual(13) = (y(28)*y(30)) - (y(38)*y(39));
    residual(14) = (y(39)) - ((1-params(10))*(1+y(35))*y(16)+y(40));
    residual(15) = (y(35)) - (y(11)+y(15)*(y(33)-y(11)));
    residual(16) = (y(40)) - (y(27)*params(13));
    residual(17) = (y(41)) - ((1+y(11))*y(18)+y(43)-y(42));
    residual(18) = (y(42)) - (params(19)+params(18)*(y(18)-params(17)));
    residual(19) = (y(43)) - ((1-params(14))*params(16)+params(14)*y(20)+x(1));
    residual(20) = (log(y(44))) - (params(15)*log(y(21))+x(2));
    residual(21) = (y(27)) - (y(43)+y(24)+y(29));
    residual(22) = (y(45)) - (y(33)-y(11));
    residual(23) = (y(46)) - (y(28)*y(30)/y(39));
end
