function [state,obs] = simul_ltro(param,gamma,bounds,policies,TT,ss,s,V,e,l,initial,m)

N_g     = size(gamma,1)/(2*4);

M       = size(e,1);

alpha    = param(1);
delta    = param(2);
beta     = param(3);
nu       = param(4); 
gamz     = param(5);
rhoz     = param(6); 
sigmaz   = param(7); 
chi      = param(8);
psi      = param(9);
omega    = param(10); 
lambda   = param(11);
csi      = param(12);     
a1       = param(13);
a2       = param(14);
g_star   = param(15);
rhog     = param(16);
sigmag   = param(17);
pi       = param(18);
iota     = param(19);
t_star   = param(20);
gamma_t  = param(21);
rhos     = param(22);
sigmas   = param(23);
s_star   = param(24);

state       = initial*ones(1,M);
net_worth   = zeros(M,1);
bond_price  = ones(M,1);
stock_price = zeros(M,1);
consumption = zeros(M,1);
output      = zeros(M,1);
investment  = zeros(M,1);
ytm         = zeros(M,1);
premia      = zeros(M,1);
risk        = zeros(M,1);
mult        = zeros(M,1);
liqui       = zeros(M,1);
return_cap  = zeros(M,1);
return_bonds= zeros(M,1);
sdf         = zeros(M,1);
RR          = zeros(M,1);
alp         = zeros(M,1);
labor       = zeros(M,1);
leverage    = zeros(M,1);
start       = zeros(M,1);
gdp_growth  = zeros(M,1);
cons_last   = exp(ss(1));

for j=2:M
 
 A                 =  [state(1,j),state(3,j)]*inv(V);
 P_tom             =  exp(A(2)+ss(7));
 P_tomnew          =  P_tom-m(j);
 A(2)              =  log(P_tomnew)-ss(7);
 A                 =  A*V;
 state(3,j)        =  A(2);
 state(6,j)        =  rhos*state(6,j-1) + sigmas*e(j,3);
 state(2,j)        =  rhoz*state(2,j-1) + sigmaz*e(j,1);
 state(4,j)        =  rhog*state(4,j-1) + sigmag*e(j,2);
   
 y           =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
 y           =  max(y,-1);
 y           =  min(y,1);
 TTT         =  T(y,s);
 
  if l(j)> -s_star-state(6,j-1)
 
 gamma_1           = (policies.consumption(:,j,2)'/TT)';
 gamma_2           = (policies.risk_free(:,j,2)'/TT)';
 gamma_4           = (policies.bond_prices(:,j,2)'/TT)';   
 gamma_liqui       = policies.mult_comp(:,j,2);  
 gamma_risk        = policies.risk_comp(:,j,2);      
 cons              = exp(ss(1)+gamma_1'*TTT);
 q                 = exp(ss(4)+gamma_4'*TTT); 
 R                 = exp(ss(2)+gamma_2'*TTT);
 rec               = param(end);
 
  else
 
 gamma_1           = (policies.consumption(:,j,1)'/TT)';
 gamma_2           = (policies.risk_free(:,j,1)'/TT)';
 gamma_4           = (policies.bond_prices(:,j,1)'/TT)';  
 gamma_liqui       = policies.mult_comp(:,j,1);  
 gamma_risk        = policies.risk_comp(:,j,1);  
 
 cons              = exp(ss(1)+gamma_1'*TTT);
 q                 = exp(ss(4)+gamma_4'*TTT);
 R                 = exp(ss(2)+gamma_2'*TTT);
 
 rec         = 0;
 
  end
  
 liqui(j)    = max(gamma_liqui'*TTT,0);
 risk(j)     = gamma_risk'*TTT; 
 premia(j)   = liqui(j)-risk(j);

 XX         = ([state(1,j);state(3,j)]'*inv(V));
 state(1,j) = XX(1);
 state(3,j) = XX(2);
 
 
 sdf(j) = beta*((cons/cons_last)*exp(-(state(2,j)+gamz)))*((1-psi)+psi*alp(j));
 
 lab    =  ((1-alpha)*(((exp(ss(5)+state(1,j))*(exp(-(state(2,j)+gamz))))^(alpha)))/...
           (chi*cons))^(1/((1/nu)+alpha));
       
gdp     =  (lab^(1-alpha))*(exp(ss(5)+state(1,j)).*exp(-(state(2,j)+gamz)))^(alpha);
 
inve    =  max((1-exp(state(4,j)+g_star))*gdp-cons,0.001);

k_tom   =  ((1-delta)*exp(ss(5)+state(1,j)) + (a1*(exp((state(2,j)+gamz))*...
                  (inve/exp(ss(5)+state(1,j))))^(1-csi)+a2)*exp(ss(5)+...
                  state(1,j)))*exp(-(state(2,j)+gamz));

Q       =  (1/((1-csi)*a1))*((exp((state(2,j)+gamz))*(inve/exp(ss(5)+state(1,j))))).^(csi);

k_ret   =  (1-delta)*Q+(alpha*gdp/exp(ss(5)+state(1,j)))*exp(state(2,j)+gamz);

b_ret   =  (1-rec)*(pi+(1-pi)*(iota+q));

b_tom   =  ((1-rec)*(pi+(1-pi)*(q+iota))*exp(ss(end)+state(5,j))*exp(-(state(2,j)+gamz))+...
           exp(state(4,j)+g_star).*gdp - exp(t_star+gamma_t*(log(exp(ss(end)+state(5,j))*(1-rec)))))./q;

n_tom   =  (psi*(k_ret*exp(ss(5)+state(1,j)) + b_ret.*exp(ss(end)+state(5,j)) -...
           exp(ss(7)+state(3,j))) + omega*(Q*exp(ss(5)+state(1,j))+...
           (1-rec)*q*exp(ss(end)+state(5,j))))*exp(-(state(2,j)+gamz));

p_tom   =  R.*(Q*k_tom+q*b_tom-n_tom); 

XX      = [log(k_tom)-ss(5);log(p_tom)-ss(7)]'*V;

state(1,j+1) = XX(1);  
state(3,j+1) = XX(2);
state(5,j+1) = log(b_tom)-ss(end);

net_worth(j)   = log(n_tom);
leverage(j)    = ((Q*k_tom+q*b_tom))/n_tom;
bond_price(j)  = q;
ytm(j)         = ((((iota+((1-q)/(1/pi)))./((1+q)/2))));
stock_price(j) = Q;
consumption(j) = log(cons);
RR(j)          = R;
labor(j)       = log(lab);
output(j)      = log(gdp);
gdp_growth(j)  = gamz + state(2,j) + (output(j)-output(j-1));
investment(j)  = log(inve);
return_cap(j)  = k_ret./stock_price(j-1);
return_bonds(j)= b_ret./bond_price(j-1);

cons_last     = cons;

end

state          = state(:,2:end);
obs.cons       = consumption(2:end);
obs.gdp        = output(2:end);
obs.inv        = investment(2:end);
obs.networth   = net_worth(2:end);
obs.q          = bond_price(2:end);
obs.Q          = stock_price(2:end);
obs.ytm        = ytm(2:end);
obs.premia     = premia(2:end);
obs.risk       = risk(2:end);
obs.mult       = max(mult(2:end),0);
obs.ret_cap    = return_cap(2:end);
obs.ret_bond   = return_bonds(2:end);
obs.R          = RR(2:end);
obs.sdf        = sdf(2:end);
obs.lab        = labor(2:end);
obs.alp        = alp(2:end);
obs.lev        = leverage(2:end);
obs.gdp_growth = gdp_growth(2:end);
obs.sdf        = sdf(2:end);
obs.liqui      = liqui(2:end);

end
