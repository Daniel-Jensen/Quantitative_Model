function [EXP] = precomp_integral_filter(coll_points,bounds,rhoz,sigmaz,s,coeff_g)

N = size(coll_points,2);

[points,weights]= GaussHermite(15);

EXP  = ones(N,N);

poli_z = zeros(N,15);

for k=1:N

x   = coll_points(:,k);

for j=1:15

    z         = rhoz*x(2)+sigmaz*points(j);
    
    y         =(2*x-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    
    y(2)      =(2*z-(bounds(2,1)+bounds(2,2)))/(bounds(2,2)-bounds(2,1));
   
    poli_z(:,j) = T(y,s);
    
end

coeff_z = ((1/sum(weights))*(poli_z*weights)+0.00000001)./(poli_z(:,8)+0.00000001);

EXP(:,k) = coeff_z.*coeff_g(:,k);
   
end

EXP = min(max(EXP,-1),1);


end

