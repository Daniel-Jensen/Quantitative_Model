function [state,obs] = simul_nodefault(param,gamma,bounds,policies,TT,ss,s,V,e,initial)

% This function generates a simulation from the model without default risk
%
% Inputs: param (vector collecting structural parameters), gamma (numerical
% solution of the model), bounds (matrix collecting upper and lower bounds 
% for the endogenous variables in the Smolyak grid), policies (structure 
% collecting the model's policy functions), TT (matrix collecting 
% the Chebyshev's polynomials evaluated at the collocation points), ss 
% (deterministic steady state), s (structure used to evaluate Chebyshev's 
% polynomials for an arbitrary point in the state space), V (matrix that 
% rotate the state variables), e (matrix collecting structural shocks),
% initial (vector collecting initial conditions for state variables).
%
% Output: state (realization of state variables), obs (structure collecting 
% the realization for several endogenous variables)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

N_g     = size(gamma,1)/4;

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

state       = initial*ones(1,M);
net_worth   = zeros(M,1);
bond_price  = ones(M,1);
stock_price = zeros(M,1);
consumption = zeros(M,1);
output      = zeros(M,1);
gdp_growth  = zeros(M,1);
investment  = zeros(M,1);
leverage    = zeros(M,1);
ytm         = zeros(M,1);
premia      = zeros(M,1);
risk        = zeros(M,1);
mult        = zeros(M,1);
return_cap  = zeros(M,1);
return_bonds= zeros(M,1);
sdf         = zeros(M,1);
RR          = zeros(M,1);
alp         = zeros(M,1);
labor       = zeros(M,1);

gamma_1  = gamma(1:N_g,:)'/TT;
gamma_2  = gamma(N_g+1:2*N_g,:)'/TT;
gamma_3  = gamma(2*N_g+1:3*N_g,:)'/TT;
gamma_4  = gamma(3*N_g+1:4*N_g,:)'/TT;
STAT     = zeros(2,M);

%=========================================================================
%                SIMULATE REALIZATION FROM MODEL
%=========================================================================

for j=2:M

 state(2,j)  =  rhoz*(state(2,j-1))+sigmaz*e(j,1); 
 state(4,j)  =  rhog*state(4,j-1)+sigmag*e(j,2);
 
 y           =  (2*state(:,j)-(bounds(:,1)+bounds(:,2)))./(bounds(:,2)-bounds(:,1));
 y           =  max(y,-1);
 y           =  min(y,1);
 TTT         =  T(y,s);
 cons        =  exp(ss(1)+gamma_1*TTT);
 R           =  exp(ss(2)+gamma_2*TTT);
 alp(j)      =  exp(ss(3)+gamma_3*TTT);
 q           = exp(ss(4)+gamma_4*TTT);
 mult(j)     = max(policies.mult*TTT,0);

 
 XX         = ([state(1,j);state(3,j)]'*inv(V));
 
 state(1,j) = XX(1);
 state(3,j) = XX(2);
 
 lab    =  ((1-alpha)*(((exp(ss(5)+state(1,j))*(exp(-(state(2,j)+gamz))))^(alpha)))/...
           (chi*cons))^(1/((1/nu)+alpha));
       
gdp     =  (lab^(1-alpha))*(exp(ss(5)+state(1,j)).*exp(-(state(2,j)+gamz)))^(alpha);
 
inve    =  max((1-exp(state(4,j)+g_star))*gdp-cons,0.001);

k_tom   =  ((1-delta)*exp(ss(5)+state(1,j)) + (a1*(exp((state(2,j)+gamz))*...
                  (inve/exp(ss(5)+state(1,j))))^(1-csi)+a2)*exp(ss(5)+...
                  state(1,j)))*exp(-(state(2,j)+gamz));

Q       =  (1/((1-csi)*a1))*((exp((state(2,j)+gamz))*(inve/exp(ss(5)+state(1,j))))).^(csi);

k_ret   =  (1-delta)*Q+(alpha*gdp/exp(ss(5)+state(1,j)))*exp(state(2,j)+gamz);

b_tom   =  ((pi+(1-pi)*(q+iota))*exp(ss(end)+state(5,j))*exp(-(state(2,j)+gamz))+...
           exp(state(4,j)+g_star).*gdp - exp(t_star+gamma_t*(ss(end)+state(5,j))))./q;

b_ret   =  (pi+(1-pi)*(iota+q));

n_tom   =  (psi*(k_ret*exp(ss(5)+state(1,j)) + b_ret.*exp(ss(end)+state(5,j)) -...
           exp(ss(7)+state(3,j))) + omega*(Q*exp(ss(5)+state(1,j))+...
           q*exp(ss(end)+state(5,j))))*exp(-(state(2,j)+gamz));

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
labor(j)       = log(lab);
output(j)      = policies.gdp*TTT;
gdp_growth(j)  = gamz + state(2,j) + (output(j)-output(j-1));
investment(j)  = log(inve);
return_cap(j)  = k_ret./stock_price(j-1);
return_bonds(j)= b_ret./bond_price(j-1);
end

%=========================================================================
%                         ORGANIZE RESULTS
%=========================================================================

state          = state(:,3:end);
obs.cons       = consumption(3:end);
obs.gdp        = output(3:end);
obs.inv        = investment(3:end);
obs.networth   = net_worth(3:end);
obs.q          = bond_price(3:end);
obs.Q          = stock_price(3:end);
obs.ytm        = ytm(3:end);
obs.premia     = premia(3:end);
obs.risk       = risk(3:end);
obs.mult       = mult(3:end);
obs.ret_cap    = return_cap(3:end);
obs.ret_bond   = return_bonds(3:end);
obs.R          = RR(3:end);
obs.sdf        = sdf(3:end);
obs.lab        = labor(3:end);
obs.alp        = alp(3:end);
obs.lev        = leverage(3:end);
obs.gdp_growth = gdp_growth(3:end);

end
