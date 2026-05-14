function [y, T, residual, g1] = dynamic_4(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(2, 1);
  residual(1)=(y(42))-(params(19)+params(18)*(y(18)-params(17)));
  residual(2)=(y(41))-((1+y(11))*y(18)+y(43)-y(42));
if nargout > 3
    g1_v = NaN(3, 1);
g1_v(1)=1;
g1_v(2)=1;
g1_v(3)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 2, 2);
end
end
