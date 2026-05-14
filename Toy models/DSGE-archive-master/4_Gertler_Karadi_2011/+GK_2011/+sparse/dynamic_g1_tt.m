function [T_order, T] = dynamic_g1_tt(y, x, params, steady_state, T_order, T)
if T_order >= 1
    return
end
[T_order, T] = GK_2011.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
T_order = 1;
if size(T, 1) < 16
    T = [T; NaN(16 - size(T, 1), 1)];
end
T(12) = getPowerDeriv(y(36)-params(15)*y(1),(-params(4)),1);
T(13) = getPowerDeriv(y(71)-y(36)*params(15),(-params(4)),1);
T(14) = getPowerDeriv(y(40)*y(12)*y(69),params(6),1);
T(15) = 1/params(10);
T(16) = (-(y(54)+params(20)))/((params(20)+y(19))*(params(20)+y(19)));
end
