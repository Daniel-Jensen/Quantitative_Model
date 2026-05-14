function [y, T] = static_1(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(20)=1;
  y(24)=params(1);
end
