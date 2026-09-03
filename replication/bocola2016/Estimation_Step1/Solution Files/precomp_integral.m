function [EXP,c_g] = precomp_integral(coll_points,bounds,rhog,sigmag,rhoz,sigmaz,s)

N = size(coll_points,2);

[points,weights]= GaussHermite(15);

EXP    = ones(N,N);
poli_z = zeros(N,15);
poli_g = zeros(N,15);
c_g    = zeros(N,1);

for k=1:N

x   = coll_points(:,k);

for j=1:15

    z         = rhoz*x(2)+sigmaz*points(j);
    
    y         =(2*x-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    
    y(2)      =(2*z-(bounds(2,1)+bounds(2,2)))/(bounds(2,2)-bounds(2,1));
   
    poli_z(:,j) = T(y,s);
    
end

coeff_z = ((1/sum(weights))*(poli_z*weights)+0.00000001)./(poli_z(:,8)+0.00000001);

x   = coll_points(:,k);

for j=1:15

    g         = rhog*x(4)+sigmag*points(j);
    
    y         =(2*x-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
    
    y(4)      =(2*g-(bounds(4,1)+bounds(4,2)))/(bounds(4,2)-bounds(4,1));
   
    poli_g(:,j) = T(y,s);
    
end

coeff_g  = (((1/sum(weights))*poli_g*weights)+0.00000001)./((poli_g(:,8)+0.00000001));

c_g(:,k) = coeff_g;

EXP(:,k) = coeff_z.*coeff_g;
   
end

EXP = min(max(EXP,-1),1);


end

