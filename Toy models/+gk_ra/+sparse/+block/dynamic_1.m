function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(43)=(1-params(14))*params(16)+params(14)*y(20)+x(1);
end
