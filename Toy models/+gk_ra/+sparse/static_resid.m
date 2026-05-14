function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(6, 1);
end
[T_order, T] = gk_ra.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(23, 1);
    residual(1) = (T(1)) - (T(1)*params(1)*(1+y(11)));
    residual(2) = (params(4)*y(2)^(1/params(3))) - (T(1)*y(3));
    residual(3) = (y(4)) - (T(3)*T(4));
    residual(4) = (y(3)) - (y(4)*(1-params(5))/y(2));
    residual(5) = (y(9)) - (y(4)*params(5)/y(5));
    residual(6) = (y(8)) - (y(6)/y(5));
    residual(7) = (y(5)) - (y(5)*(1-params(6))+y(5)*T(5));
    residual(8) = (y(7)) - (1/T(6));
    residual(9) = (1+y(10)) - ((y(9)+(1-params(6))*y(7))/y(7));
    residual(10) = (y(13)) - (params(12)*(params(10)+(1-params(10))*params(11)*y(15))*(y(10)-y(11)));
    residual(11) = (y(14)) - ((1+y(11))*params(12)*(params(10)+(1-params(10))*params(11)*y(15)));
    residual(12) = (y(15)) - (y(14)/(params(11)-y(13)));
    residual(13) = (y(5)*y(7)) - (y(15)*y(16));
    residual(14) = (y(16)) - (y(16)*(1-params(10))*(1+y(12))+y(17));
    residual(15) = (y(12)) - (y(11)+y(15)*(y(10)-y(11)));
    residual(16) = (y(17)) - (y(4)*params(13));
    residual(17) = (y(18)) - ((1+y(11))*y(18)+y(20)-y(19));
    residual(18) = (y(19)) - (params(19)+params(18)*(y(18)-params(17)));
    residual(19) = (y(20)) - ((1-params(14))*params(16)+y(20)*params(14)+x(1));
    residual(20) = (log(y(21))) - (log(y(21))*params(15)+x(2));
    residual(21) = (y(4)) - (y(20)+y(1)+y(6));
    residual(22) = (y(22)) - (y(10)-y(11));
    residual(23) = (y(23)) - (y(5)*y(7)/y(16));
end
