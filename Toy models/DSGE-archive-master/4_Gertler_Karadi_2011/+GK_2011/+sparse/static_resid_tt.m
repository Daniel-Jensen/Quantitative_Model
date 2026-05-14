function [T_order, T] = static_resid_tt(y, x, params, T_order, T)
if T_order >= 0
    return
end
T_order = 0;
if size(T, 1) < 5
    T = [T; NaN(5 - size(T, 1), 1)];
end
T(1) = params(23)*(y(5)*y(12)*y(34))^params(6);
T(2) = y(3)^(1-params(6));
T(3) = params(12)*y(3)^(1/params(5));
T(4) = y(17)/params(10);
T(5) = (y(1)-y(1)*params(15))^(-params(4));
end
