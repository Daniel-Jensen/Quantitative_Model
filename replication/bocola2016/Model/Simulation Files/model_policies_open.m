function [policies,expectations,prem,residuals] = model_policies(param,gamma,bounds,EXP,TT,coll_points,ss,s)

N_g      = size(coll_points,2);

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
y_ss     = (exp(ss(6))^(alpha))*(0.33^(1-alpha));
def      = [0,1];

gam_1     = zeros(N_g,2);
gam_2     = zeros(N_g,2);
gam_3     = zeros(N_g,2);
gam_4     = zeros(N_g,2);
gam_5     = zeros(N_g,2);

exp_11     = zeros(2,N_g);
exp_22     = zeros(2,N_g);
exp_33     = zeros(2,N_g);
exp_44     = zeros(2,N_g);

mu        = zeros(N_g,2);
cc        = zeros(N_g,2);
ll        = zeros(N_g,2);
qq        = zeros(N_g,2);
aa        = zeros(N_g,2);
RR        = zeros(N_g,2);
RR_wc     = zeros(N_g,2);

QQ        = zeros(N_g,2);

x         = coll_points;

GDP       = zeros(N_g,2);
INV       = zeros(N_g,2);
DEBT      = zeros(N_g,2);
B_f       = zeros(N_g,2);
KAP       = zeros(N_g,2);
PROM      = zeros(N_g,2);
PREMIUM   = zeros(N_g,2);
YTM       = zeros(N_g,2);
prem      = zeros(2,N_g);
residuals = zeros(2,N_g,5);

x         = coll_points;

for j=1:2

gam     = gamma(1+(j-1)*N_g*5:N_g*5+(j-1)*N_g*5);    
gamma_1 = gam(1:N_g);
gamma_2 = gam(N_g+1:2*N_g);
gamma_3 = gam(2*N_g+1:3*N_g);
gamma_4 = gam(3*N_g+1:4*N_g);
gamma_5 = gam(4*N_g+1:5*N_g);
        
cons    =  exp(ss(1)+gamma_1);
alp     =  exp(ss(2)+gamma_2);
q       =  exp(ss(3)+gamma_3);
R_wc    =  exp(ss(4)+gamma_4);
B_for   =  ss(5) + gamma_5;

lab     =  min((((1-alpha).*(((exp(x(2,:)).^(1-alpha)).*(exp(ss(6)+x(1,:)).^(alpha))))...
           ./(chi*((1-psi) + psi*R_wc)))).^(1/((1/nu)+alpha)),1);
       
gdp     =  ((exp(x(2,:)).*lab).^(1-alpha)).*exp(ss(6)+x(1,:)).^(alpha);

R       =  (1/beta) + 0.01*(B_for./gdp);

inv     =  max((1-exp(x(4,:)+g_star)).*gdp-cons + B_for./R - x(6,:),0.01);

K_tom   =  (1-delta)*exp(ss(6)+x(1,:)) + (a1*(inv./exp(ss(6)+x(1,:))).^(1-tau)+a2).*exp(ss(6)+x(1,:));

Q       =  (1/((1-tau)*a1))*(inv./exp(ss(6)+x(1,:))).^(tau);

k_ret   =  (1-delta)*Q+(alpha*gdp./exp(ss(6)+x(1,:)));

B_tom   =  ((1-def(j)*rec).*(pi+(1-pi)*(q+iota)).*exp(ss(end-1)+x(5,:))+exp(x(4,:)+...
           g_star).*gdp - exp(t_star+gamma_t*(ss(end-1)+x(5,:)).*(1-def(j)*rec)))./q;

b_ret   =  (1-def(j)*rec).*(pi+(1-pi)*(iota+q));

N_tom   =  max(Theta*(k_ret.*exp(ss(6)+x(1,:)) + b_ret.*exp(ss(end-1)+x(5,:)) -exp(ss(8)+x(3,:))) + ...
            omega*(Q.*exp(ss(6)+x(1,:))+q.*exp(ss(end-1)+x(5,:)).*(1-def(j)*rec)+psi.*(1-alpha).*gdp),0.4);

P_tom   =  R.*(Q.*K_tom+q.*B_tom+psi*(1-alpha)*gdp-N_tom) - psi*R_wc.*(1-alpha).*gdp; 

x_tom   = [log(K_tom)-ss(6);rhoz*coll_points(2,:);log(P_tom)-ss(8);rhog*coll_points(4,:);...
          log(B_tom)-ss(end-1);B_for-ss(end);rhos*coll_points(7,:)];

y_tom   = (2*x_tom-(bounds(:,1)+bounds(:,2))*ones(1,N_g))./((bounds(:,2)-...
           bounds(:,1))*ones(1,N_g));

y_tom   = max(y_tom,-1);

y_tom   = min(y_tom,1);

TT1     = T(y_tom,s);

exp_11(j,:)  = log(beta*(cons-((chi*lab.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma));

exp_22(j,:)  = log(exp(exp_11(j,:)).*alp);

exp_33(j,:)  = log((exp(exp_11(j,:)).*((1-Theta)+Theta*alp).*k_ret));

exp_44(j,:)  = log((exp(exp_11(j,:)).*((1-Theta)+Theta*alp).*b_ret));

prem(j,:)    = (k_ret);

gam_1(:,j)   = ((exp_11(j,:)/TT)*(TT1.*EXP))';

gam_2(:,j)   = ((exp_22(j,:)/TT)*(TT1.*EXP))';

gam_3(:,j)   = ((exp_33(j,:)/TT)*(TT1.*EXP))';

gam_4(:,j)   = ((exp_44(j,:)/TT)*(TT1.*EXP))';

PREMIUM(:,j) = ((prem(j,:)/TT)*(TT1.*EXP))';

cc(:,j)      = cons';

ll(:,j)      = lab';

qq(:,j)      = q';

RR(:,j)      = R';

RR_wc(:,j)   = R_wc';

aa(:,j)      = alp';

mu(:,j)      = (N_tom./(lambda*(Q.*K_tom+q.*B_tom +psi*(1-alpha)*gdp)));

QQ(:,j)      = Q';

GDP(:,j)     = log(gdp)';

INV(:,j)     = log(inv)';

DEBT(:,j)    = log(B_tom)';

KAP(:,j)     = log(K_tom)';

PROM(:,j)    = log(P_tom)';

YTM(:,j)     = ((((iota+((1-q')/(1/pi)))./((1+q')/2))))-(R'-1);

B_f(:,j)     = B_for;

end

Sprob            = (exp(s_star+coll_points(end,:))./(1+exp(s_star+coll_points(end,:))))';

exp_1            = ((1-Sprob).*exp(gam_1(:,1)) + Sprob.*exp(gam_1(:,2)))*ones(1,2);

exp_2            = ((1-Sprob).*exp(gam_2(:,1)) + Sprob.*exp(gam_2(:,2)))*ones(1,2);

exp_3            = ((1-Sprob).*exp(gam_3(:,1)) + Sprob.*exp(gam_3(:,2)))*ones(1,2);

exp_4            = ((1-Sprob).*exp(gam_4(:,1)) + Sprob.*exp(gam_4(:,2)))*ones(1,2);

mu               = max(1-(((1-Theta)+Theta*RR.*(exp_2./(cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma))).*mu),0);

PREMIUM          = ((((1-Sprob).*PREMIUM(:,1)+ Sprob.*PREMIUM(:,2))*ones(1,2))./QQ);

SDF              = ((1-Theta)./RR)+Theta*(exp_2./(cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma));

RISK             = (exp_3./QQ)./(cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma) - SDF.*PREMIUM;

R_impl           = (exp_1./(cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma)).^(-1);

alp_impl         = ((1-Theta)+Theta*R_impl.*(exp_2./(cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma)))./(1-mu);

q_impl           = (exp_4./((cc-((chi*ll.^(1+(1/nu)))/(1+(1/nu)))).^(-sigma)))./(lambda*mu+(1-mu).*(alp_impl));

R_wc_impl        = R_impl + lambda*mu./SDF;

residuals(:,:,1) = log(RR./R_impl)';
residuals(:,:,2) = log(aa./alp_impl)';
residuals(:,:,3) = log(qq./q_impl)';
residuals(:,:,4) = (PREMIUM-(((lambda*mu+(1-mu).*(alp_impl))-RISK)./SDF))';
residuals(:,:,5) = log(RR_wc./R_wc_impl)';

residuals        = reshape(residuals,N_g*2*5,1);

PREMIUM          = (PREMIUM-RR);
MULT_comp        = lambda*mu./(SDF);
RISK_comp        = RISK./SDF;



policies.mult         = (mu'/TT)';
policies.premia       = (PREMIUM'/TT)';
policies.risk         = (RISK'/TT)';
policies.gdp          = ((GDP-log(y_ss))'/TT)';
policies.inv          = (INV'/TT)';
policies.ytm          = (YTM'/TT)';
policies.q            = (qq'/TT)';
policies.Q            = (QQ'/TT)';
policies.debt         = ((DEBT-ss(end-1))'/TT)';
policies.prom         = ((PROM-ss(8))'/TT)';
policies.kap          = ((KAP-ss(6))'/TT)';
policies.R            = (RR'/TT)';
policies.alp          = (aa'/TT)';
policies.cons         = (cc'/TT)';
policies.lab          = (ll'/TT)';
policies.R_wc         = (RR_wc'/TT)';
policies.mult_comp    = (MULT_comp'/TT)';
policies.risk_comp    = (RISK_comp'/TT)';
policies.Bf           = (B_f'/TT)';

expectations          = [exp_11',exp_22',exp_33',exp_44'];

end

