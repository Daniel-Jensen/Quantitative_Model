function [y, T, residual, g1] = static_5(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(3, 1);
  T(2)=params(8)-1+params(14)*y(17)/params(10)*(y(17)/params(10)-1)-y(17)/params(10)*(y(17)/params(10)-1)*y(24)*params(14);
  residual(1)=(1/y(16))-(params(8)/T(2));
  residual(2)=(y(35))-((1-params(16))*(params(11)+params(17)*(y(17)-params(10))+params(18)*(params(8)/(params(8)-1)-1/y(16)))+y(35)*params(16));
  residual(3)=(y(15))-((1+y(35))/y(17));
if nargout > 3
    g1_v = NaN(7, 1);
g1_v(1)=(-1)/(y(16)*y(16));
g1_v(2)=(-((1-params(16))*params(18)*(-((-1)/(y(16)*y(16))))));
g1_v(3)=(-((-(params(8)*((y(17)/params(10)-1)*params(14)*1/params(10)+params(14)*y(17)/params(10)*1/params(10)-((y(17)/params(10)-1)*y(24)*params(14)*1/params(10)+y(17)/params(10)*y(24)*params(14)*1/params(10)))))/(T(2)*T(2))));
g1_v(4)=(-((1-params(16))*params(17)));
g1_v(5)=(-((-(1+y(35)))/(y(17)*y(17))));
g1_v(6)=1-params(16);
g1_v(7)=(-(1/y(17)));
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 3, 3);
end
end
