function [T_order, T] = static_g1_tt(y, x, params, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = GK_2011.sparse.static_resid_tt(y, x, params, T_order, T);
T_order = 1;
if size(T, 1) < 9
    T = [T; NaN(9 - size(T, 1), 1)];
end
T(6) = (1-params(15))*getPowerDeriv(y(1)-y(1)*params(15),(-params(4)),1);
T(7) = getPowerDeriv(y(5)*y(12)*y(34),params(6),1);
T(8) = 1/params(10);
T(9) = (params(8)-1+params(14)*T(4)*(T(4)-1)-T(4)*(T(4)-1)*y(24)*params(14))*(params(8)-1+params(14)*T(4)*(T(4)-1)-T(4)*(T(4)-1)*y(24)*params(14));
end
