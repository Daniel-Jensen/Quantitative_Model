function [y] = simul_filter(bounds,param,e,initial,policies,s)

M        = size(e,1);
gamz     = param(5);
rhoz     = param(6); 
sigmaz   = param(7); 
rhog     = param(16);
sigmag   = param(17);

state    = initial*ones(1,M);

for j=2:M
    
 state(2,j)  =  rhoz*state(2,j-1)+sigmaz*e(j,1);
 state(4,j)  =  rhog*state(4,j-1)+sigmag*e(j,2);
 y           =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
 y           =  max(y,-1);
 y           =  min(y,1);
 TT          =  T(y,s); 
 gdp         =  policies.gdp*TT;
 mult        =  max(policies.mult*TT,0);

end

y            = [state(2,j)+gamz+gdp;mult];

end
