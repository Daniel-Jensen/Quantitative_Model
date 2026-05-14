function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 6
    T = [T; NaN(6 - size(T, 1), 1)];
end
T(1) = y(1)^(-params(2));
T(2) = y(5)^params(5);
T(3) = y(21)*T(2);
T(4) = y(2)^(1-params(5));
T(5) = params(8)*y(8)^(1-params(7))+params(9);
T(6) = params(8)*(1-params(7))*y(8)^(-params(7));
end
