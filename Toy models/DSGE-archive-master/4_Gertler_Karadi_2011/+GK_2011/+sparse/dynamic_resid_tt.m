function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 11
    T = [T; NaN(11 - size(T, 1), 1)];
end
T(1) = params(23)*(y(40)*y(12)*y(69))^params(6);
T(2) = y(38)^(1-params(6));
T(3) = params(12)*y(38)^(1/params(5));
T(4) = (y(54)+params(20))/(params(20)+y(19))-1;
T(5) = params(13)*(y(54)+params(20))/(params(20)+y(19));
T(6) = ((params(20)+y(89))/(y(54)+params(20)))^2;
T(7) = params(13)*y(59)*T(6);
T(8) = (params(20)+y(89))/(y(54)+params(20))-1;
T(9) = y(87)/params(10);
T(10) = T(9)*y(59)*params(14)*(T(9)-1);
T(11) = params(8)-1+params(14)*y(52)/params(10)*(y(52)/params(10)-1)-T(10)*y(84)/y(49);
end
