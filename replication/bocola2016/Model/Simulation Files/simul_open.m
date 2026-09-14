function [state,obs,dd] = simul_open(bounds,param,e,l,gamma,policies,TT,ss,initial,s)

N_g     = size(gamma,1)/(2*5);

M       = size(e,1);

% Set parameters
alpha    = param(1);
delta    = param(2);
beta     = param(3);
sigma    = param(4);   
nu       = param(5);  
rhoz     = param(6);
sigmaz   = param(7);
chi      = param(8);
Theta    = param(9);
omega    = param(10); 
lambda   = param(11);
tau      = param(12);     
a1       = param(13);
a2       = param(14);
g_star   = param(15);
rhog     = param(16);
sigmag   = param(17);
pi       = param(18);
iota     = param(19);
t_star   = param(20);
gamma_t  = param(21);
psi      = param(22);
rhos     = param(23);
sigmas   = param(24);
s_star   = param(25);
rec      = param(26);

state       = initial*ones(1,M);
net_worth   = zeros(M,1);
bond_price  = ones(M,1);
stock_price = zeros(M,1);
consumption = zeros(M,1);
output      = zeros(M,1);
investment  = zeros(M,1);
leverage    = zeros(M,1);
ytm         = zeros(M,1);
premia      = zeros(M,1);
risk        = zeros(M,1);
mult        = zeros(M,1);
return_cap  = zeros(M,1);
return_bonds= zeros(M,1);
risk_comp   = zeros(M,1);
mult_comp   = zeros(M,1);
sdf         = zeros(M,1);
RR          = zeros(M,1);
alp         = zeros(M,1);
labor       = zeros(M,1);

 gamma_1     = (gamma(1:N_g)'/TT)';
 gamma_2     = (gamma(N_g+1:2*N_g)'/TT)';
 gamma_3     = (gamma(2*N_g+1:3*N_g)'/TT)';
 gamma_4     = (gamma(3*N_g+1:4*N_g)'/TT)';
 gamma_5     = (gamma(4*N_g+1:5*N_g)'/TT)';
 gamma_1_def =  (gamma(5*N_g+1:6*N_g)'/TT)';
 gamma_2_def = (gamma(6*N_g+1:7*N_g)'/TT)';
 gamma_3_def = (gamma(7*N_g+1:8*N_g)'/TT)';
 gamma_4_def = (gamma(8*N_g+1:9*N_g)'/TT)';
 gamma_5_def = (gamma(9*N_g+1:10*N_g)'/TT)';
 
 gamma_prem  = policies.premia(:,1);
 gamma_mult  = policies.mult(:,1);
 gamma_risk  = policies.risk(:,1);
 gamma_prem_def  = policies.premia(:,2);
 gamma_mult_def  = policies.mult(:,2);
 gamma_risk_def  = policies.risk(:,2);
 
 dd          = zeros(M,1);
  
for j=2:M
  
 state(7,j)  =  rhos*state(7,j-1) + sigmas*e(j,3);

 state(2,j)  =  rhoz*state(2,j-1)+sigmaz*e(j,1);
 
 state(4,j)  =  rhog*state(4,j-1)+sigmag*e(j,2);
 
 y           =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));

 y           =  max(y,-1);
 
 y           =  min(y,1);
 
 TT          =  T(y,s);
 
 if l(j) > - s_star-state(7,j)

 cons        = exp(ss(1)+gamma_1_def'*TT);
 alp(j)      = exp(ss(2)+gamma_2_def'*TT);
 q           = exp(ss(3)+gamma_3_def'*TT);
 R_wc        = exp(ss(4)+gamma_4_def'*TT);
 B_for       = ss(5)+gamma_5_def'*TT;
 rec         = param(end);
 
mult(j)      = max(gamma_mult_def'*TT,0);
premia(j)    = gamma_prem_def'*TT;
risk(j)      = gamma_risk_def'*TT;
l            = -10000*ones(10000,1);
dd(j)        = 1;

 else
     
 cons        = exp(ss(1)+gamma_1'*TT);
 alp(j)      = exp(ss(2)+gamma_2'*TT);
 q           = exp(ss(3)+gamma_3'*TT);
 R_wc        = exp(ss(4)+gamma_4'*TT);
 B_for       = ss(5)+gamma_5'*TT;
 
 rec         = 0;
 
mult(j)       = max(gamma_mult'*TT,0);
premia(j)     = gamma_prem'*TT;
risk(j)       = gamma_risk'*TT;

 end
  
 lab         =  min((((1-alpha).*(((exp(state(2,j)).^(1-alpha)).*(exp(ss(6)+state(1,j)).^(alpha))))...
           ./(chi*((1-psi)+psi*R_wc)))).^(1/((1/nu)+alpha)),1);

gdp          =  log(((exp(state(2,j)).*lab).^(1-alpha)).*exp(ss(6)+state(1,j)).^(alpha));

R            =  (1/beta)+0.01*(B_for/exp(gdp));

inv          =  max((1-exp(state(4,j)+g_star)).*exp(gdp)-cons + (B_for/R)-state(6,j),0.001);

state(1,j+1) =  log((1-delta)*exp(ss(6)+state(1,j)) + (a1*(inv./exp(ss(6)+state(1,j))).^(1-tau)+a2).*exp(ss(6)+state(1,j)));

Q            =  (1/((1-tau)*a1))*(inv./exp(ss(6)+state(1,j))).^(tau);

k_ret        =  (1-delta)*Q+(alpha*exp(gdp)./exp(ss(6)+state(1,j)));

state(5,j+1) =  log(((1-rec)*(pi+(1-pi)*(q+iota)).*exp(ss(end-1)+state(5,j))+exp(state(4,j)+...
                g_star).*exp(gdp) - exp(t_star+gamma_t*(ss(end-1)+state(5,j))*(1-rec)))./q);

b_ret        =  (1-rec)*(pi+(1-pi)*(iota+q));

N_tom        =  max(Theta*(k_ret.*exp(ss(6)+state(1,j)) + b_ret.*exp(ss(end-1)+state(5,j)) -exp(ss(8)+state(3,j))) + ...
                omega*(Q.*exp(ss(6)+state(1,j))+q.*exp(ss(end-1)+state(5,j))*(1-rec) + psi*(1-alpha)*exp(gdp)));

leverage(j)  = ((Q.*exp(state(1,j+1))+q.*exp(state(5,j+1))+psi*(1-alpha)*gdp))./N_tom;

state(3,j+1) =  log(R.*(Q.*exp(state(1,j+1))+q.*exp(state(5,j+1))-N_tom) - psi*(1-alpha)*(R_wc-R)*exp(gdp)); 

state(1,j+1) = state(1,j+1) - ss(6);

state(3,j+1) = state(3,j+1) - ss(8);

state(5,j+1) = state(5,j+1) - ss(end-1);

state(6,j+1) = B_for - ss(end);

net_worth(j) = log(N_tom);

bond_price(j) = q;

ytm(j)        = ((((iota+((1-q)/(1/pi)))./((1+q)/2))));

stock_price(j)= Q;

consumption(j) = log(cons);

labor(j)       = log(lab);

output(j)     = gdp;

investment(j) = log(inv);

return_cap(j)   = k_ret./stock_price(j-1);

return_bonds(j) = b_ret./bond_price(j-1);

RR(j)         = R_wc;

sdf(j)        = beta*(((exp(consumption(j))-chi*(exp(labor(j))^(1+(1/nu)))/(1+...
               (1/nu)))^(-sigma))/((exp(consumption(j-1))-chi*(exp(labor(j-1))^...
               (1+(1/nu)))/(1+(1/nu)))^(-sigma)))*((1-Theta)+Theta*alp(j));

risk_comp(j)  = policies.risk_comp(:,1)'*TT;
mult_comp(j)  = policies.mult_comp(:,1)'*TT;
end

state        = state(:,2:end);
obs.cons     = consumption(2:end);
obs.gdp      = output(2:end);
obs.inv      = investment(2:end);
obs.networth = net_worth(2:end);
obs.q        = bond_price(2:end);
obs.Q        = stock_price(2:end);
obs.ytm      = ytm(2:end);
obs.premia   = premia(2:end);
obs.risk     = risk(2:end);
obs.mult     = mult(2:end);
obs.ret_cap  = return_cap(2:end);
obs.ret_bonds = return_bonds(2:end);
obs.R       = RR(2:end);
obs.sdf     = sdf(2:end);
obs.risk_comp = risk_comp(2:end);
obs.mult_comp = mult_comp(2:end);
obs.lab     = labor(2:end);
obs.alp     = alp(2:end);
obs.lev     = leverage(2:end);

end
