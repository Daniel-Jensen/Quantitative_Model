function [state,obs,STATE,dd,start] = simul_lown(param,gamma,bounds,policies,TT,ss,s,V,e,l,initial)

% This function generates a simulation from the model 
%
% Inputs: param (vector collecting structural parameters), gamma (numerical
% solution of the model), bounds (matrix collecting upper and lower bounds 
% for the endogenous variables in the Smolyak grid), policies (structure 
% collecting the model's policy functions), TT (matrix collecting 
% the Chebyshev's polynomials evaluated at the collocation points), ss 
% (deterministic steady state), s (structure used to evaluate Chebyshev's 
% polynomials for an arbitrary point in the state space), V (matrix that 
% rotate the state variables), e (matrix collecting structural shocks), l
% (vector collecting realization of logistic), initial (vector collecting 
% initial conditions for state variables).
%
% Output: state (realization of state variables), obs (structure collecting 
% the realization for several endogenous variables)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

N_g     = size(TT,1);

M        = size(e,1);
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
STATE       = initial*ones(1,M);
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
cons_growth = zeros(M,1);

gamma_1          = (gamma(1:N_g)'/TT)';
gamma_2          = (gamma(N_g+1:2*N_g)'/TT)';
gamma_3          = (gamma(2*N_g+1:3*N_g)'/TT)';
gamma_4          = (gamma(3*N_g+1:4*N_g)'/TT)';
gamma_1_def      = (gamma(4*N_g+1:5*N_g)'/TT)';
gamma_2_def      = (gamma(5*N_g+1:6*N_g)'/TT)';
gamma_3_def      = (gamma(6*N_g+1:7*N_g)'/TT)';
gamma_4_def      = (gamma(7*N_g+1:8*N_g)'/TT)';
gamma_liqui      = policies.mult_comp(:,1);
gamma_risk       = policies.risk_comp(:,1);
gamma_mult       = policies.mult(:,1);
gamma_liqui_def  = policies.mult_comp(:,2);
gamma_risk_def   = policies.risk_comp(:,2);
gamma_mult_def   = policies.mult(:,2);
cons_last        = exp(ss(1));
dd               = zeros(M,1);

for j=2:M
    
 A                 =  [state(1,j),state(3,j)]*inv(V);
 P_tom             =  exp(A(2)+ss(7));
 P_tomnew          =  P_tom+10;
 A(2)              =  log(P_tomnew)-ss(7);
 A                 =  A*V;
 state(3,j)        =  A(2);
 state(1,j)        =  A(1);    
 state(2,j)  =  rhoz*state(2,j-1)  + sigmaz*e(j,1); 
 state(4,j)  =  rhog*state(4,j-1)  + sigmag*e(j,2);
 state(6,j)  =  rhos*state(6,j-1)  + sigmas*e(j,3);

 y           =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
 y           =  max(y,-1);
 y           =  min(y,1);
 TTT         =  T(y,s);
 
 start(j)    =  y(3);
 
 STATE(:,j)  = state(:,j);
 
 if l(j) > -s_star-state(6,j-1)
 
 cons        = exp(ss(1)+gamma_1_def'*TTT); 
 R           = exp(ss(2)+gamma_2_def'*TTT);
 alp(j)      = exp(ss(3)+gamma_3_def'*TTT);
 q           = exp(ss(4)+gamma_4_def'*TTT);
 rec         = param(end);
mult(j)      = max(gamma_mult_def'*TTT,0);
liqui(j)     = max(gamma_liqui_def'*TTT,0);
risk(j)      = gamma_risk_def'*TTT;
premia(j)    = liqui(j)-risk(j);
l            = -10000*ones(10000,1);
dd(j)        = 1;

 else
     
 cons        = exp(ss(1)+gamma_1'*TTT);
 R           = exp(ss(2)+gamma_2'*TTT);
 alp(j)      = exp(ss(3)+gamma_3'*TTT);
 q           = exp(ss(4)+gamma_4'*TTT);
 
 rec         = 0;
mult(j)      = max(gamma_mult'*TTT,0);
liqui(j)     = max(gamma_liqui'*TTT,0);
risk(j)      = gamma_risk'*TTT;
premia(j)    = liqui(j)-risk(j);

 end

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
cons_growth(j) = gamz + state(2,j) + (consumption(j)-consumption(j-1));
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
obs.cons_growth= cons_growth(2:end);
obs.sdf        = sdf(2:end);
obs.liqui      = liqui(2:end);

end
