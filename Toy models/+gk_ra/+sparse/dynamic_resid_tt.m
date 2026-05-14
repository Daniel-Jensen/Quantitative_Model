function [T_order, T] = dynamic_resid_tt(y, x, params, steady_state, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 7
    T = [T; NaN(7 - size(T, 1), 1)];
end
T(1) = y(24)^(-params(2));
T(2) = y(47)^(-params(2));
T(3) = y(5)^params(5);
T(4) = y(44)*T(3);
T(5) = y(25)^(1-params(5));
T(6) = params(8)*y(31)^(1-params(7))+params(9);
T(7) = params(8)*(1-params(7))*y(31)^(-params(7));
end
