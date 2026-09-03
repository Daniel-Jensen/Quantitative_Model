function [residuals,expectations,PREMIUM,RISK_comp,MULT_comp,prem1,policies] = model_ltro_lastperiod_policies(param,gamma,bounds,coll_points,TT,s,ss,EXP,V,m,expectations,prem,eps)

N_g      = size(coll_points,2);

alpha        = param(1);
delta        = param(2);
beta         = param(3);
nu           = param(4); 
gamz         = param(5);
rhoz         = param(6);   
chi          = param(8);
psi          = param(9);
omega        = param(10); 
lambda       = param(11);
csi          = param(12);     
a1           = param(13);
a2           = param(14);
g_star       = param(15);
rhog         = param(16);
pi           = param(18);
iota         = param(19);
t_star       = param(20);
gamma_t      = param(21);
rhos         = param(22);
s_star       = param(24);
rec          = param(25)*ones(1,N_g);
rec(1:71)    = param(25)/eps;
def          = [0,1];

gam_1     = zeros(N_g,2);
gam_2     = zeros(N_g,2);
gam_3     = zeros(N_g,2);
gam_4     = zeros(N_g,2);
gam_5     = zeros(N_g,2);

exp_11    = zeros(N_g,2)';
exp_22    = zeros(N_g,2)';
exp_33    = zeros(N_g,2)';
exp_44    = zeros(N_g,2)';
exp_55    = zeros(N_g,2)';

mu        = zeros(N_g,2);
cc        = zeros(N_g,2);
qq        = zeros(N_g,2);
aa        = zeros(N_g,2);
RR        = zeros(N_g,2);
QQ        = zeros(N_g,2);

residuals = zeros(2,N_g,4);

xx        = (([coll_points(1,:);coll_points(3,:)])'*inv(V))';

x         = coll_points;
x(1,:)    = xx(1,:);
x(3,:)    = xx(2,:);
x(2,:)    = x(2,:)+gamz;

GDP       = zeros(N_g,2);
INV       = zeros(N_g,2);
DEBT      = zeros(N_g,2);
KAP       = zeros(N_g,2);
PROM      = zeros(N_g,2);
PREMIUM   = zeros(N_g,2);
YTM       = zeros(N_g,2);
prem1     = zeros(2,N_g);

for j=1:2

gam     = gamma(1+(j-1)*N_g*4:N_g*4+(j-1)*N_g*4);    
    
gamma_1 = gam(1:N_g);
gamma_2 = gam(N_g+1:2*N_g);
gamma_3 = gam(2*N_g+1:3*N_g);
gamma_4 = gam(3*N_g+1:4*N_g);
        
cons    =  exp(ss(1)+gamma_1);
R       =  exp(ss(2)+gamma_2);
alp     =  exp(ss(3)+gamma_3);
q       =  exp(ss(4)+gamma_4);
        
%=========================================================================
%                COMPUTE ENDOGENOUS VARIABLES GIVEN GUESS
%=========================================================================

lab     =  (((1-alpha).*(((exp(ss(5)+x(1,:))).*exp(-x(2,:))).^(alpha))...
           ./(chi*(cons)))).^(1/((1/nu)+alpha));
       
gdp     =  (lab.^(1-alpha)).*(exp(ss(5)+x(1,:)).*exp(-x(2,:))).^(alpha);

inve    =  max((1-exp(x(4,:)+g_star)).*gdp-cons,0.01);

K_tom   =  ((1-delta)*exp(ss(5)+x(1,:)) + (a1*(exp(x(2,:)).*(inve./exp(ss(5)+...
            x(1,:)))).^(1-csi)+a2).*exp(ss(5)+x(1,:))).*exp(-x(2,:));

Q       =  (1/((1-csi)*a1))*((exp(x(2,:)).*(inve./exp(ss(5)+x(1,:))))).^(csi);

k_ret   =  (1-delta)*Q+(alpha*gdp./exp(ss(5)+x(1,:))).*exp(x(2,:));

B_tom   =  ((1-def(j)*rec).*(pi+(1-pi)*(q+iota)).*exp(ss(end)+x(5,:)).*exp(-x(2,:))+exp(x(4,:)+...
           g_star).*gdp - exp(t_star+gamma_t*log(exp((ss(end)+x(5,:))).*(1-def(j)*rec))))./q;

b_ret   =  (1-def(j)*rec).*(pi+(1-pi)*(iota+q));

N_tom   =  max((psi*(k_ret.*exp(ss(5)+x(1,:)) + b_ret.*exp(ss(end)+x(5,:)) -...
           exp(ss(7)+x(3,:))) + omega*(Q.*exp(ss(5)+x(1,:))+q.*exp(ss(end)+...
           x(5,:)).*(1-def(j)*rec))).*exp(-x(2,:)),0.65);

P_tom   =  R.*(Q.*K_tom+q.*B_tom-N_tom); 

%=========================================================================
%                         COMPUTE EXPECTATIONS
%=========================================================================

X       = [log(K_tom)'-ss(5),log(P_tom)'-ss(7)];

X_tilde = (X*V)';

x_tom   = [X_tilde(1,:);rhoz*coll_points(2,:);X_tilde(2,:);rhog*coll_points(4,:);...
          log(B_tom)-ss(end);rhos*coll_points(6,:)];

y_tom   = (2*x_tom-(bounds(:,1)+bounds(:,2))*ones(1,N_g))./((bounds(:,2)-...
           bounds(:,1))*ones(1,N_g));

y_tom   = max(y_tom,-1);
y_tom   = min(y_tom,1);

TT_exp  = T(y_tom,s);

exp_11(j,:)  = log(beta*(exp(-x(2,:))).*cons.^(-1));
exp_22(j,:)  = log(exp(exp_11(j,:)).*alp);
exp_33(j,:)  = log((exp(exp_11(j,:)).*((1-psi)+psi*alp).*k_ret));
exp_44(j,:)  = log((exp(exp_11(j,:)).*((1-psi)+psi*alp).*b_ret));
exp_55(j,:)  = log(alp.*m);

exp_1        = expectations(:,j)';
exp_2        = expectations(:,j+2)';
exp_3        = expectations(:,j+4)';
exp_4        = expectations(:,j+6)';

gam_1(:,j)   = ((exp_1/TT)*(TT_exp.*EXP))';
gam_2(:,j)   = ((exp_2/TT)*(TT_exp.*EXP))';
gam_3(:,j)   = ((exp_3/TT)*(TT_exp.*EXP))';
gam_4(:,j)   = ((exp_4/TT)*(TT_exp.*EXP))';

PREMIUM(:,j) = exp(((prem(j,:)/TT)*(TT_exp.*EXP)))';

cc(:,j)      = cons';
qq(:,j)      = q';
RR(:,j)      = R';
aa(:,j)      = alp';
mu(:,j)      = ((N_tom-m)./(lambda*(Q.*K_tom+q.*B_tom)));
QQ(:,j)      = Q';
GDP(:,j)     = log(gdp)';
INV(:,j)     = log(inve)';
DEBT(:,j)    = log(B_tom)';
KAP(:,j)     = log(K_tom)';
PROM(:,j)    = log(P_tom)';
YTM(:,j)     = ((((iota+((1-q')/(1/pi)))./((1+q')/2))))-(R'-1);

end


Sprob            = (exp(s_star+coll_points(end,:))./(1+exp(s_star+coll_points(end,:))))';

exp_1            = ((1-Sprob).*exp(gam_1(:,1)) + Sprob.*exp(gam_1(:,2)))*ones(1,2);
exp_2            = ((1-Sprob).*exp(gam_2(:,1)) + Sprob.*exp(gam_2(:,2)))*ones(1,2);
exp_3            = ((1-Sprob).*exp(gam_3(:,1)) + Sprob.*exp(gam_3(:,2)))*ones(1,2);
exp_4            = ((1-Sprob).*exp(gam_4(:,1)) + Sprob.*exp(gam_4(:,2)))*ones(1,2);

mu               = max(1-(((1-psi)+psi*RR.*(exp_2.*cc)).*mu),0);
PREMIUM          = ((((1-Sprob).*PREMIUM(:,1)+ Sprob.*PREMIUM(:,2))*ones(1,2))./QQ);
SDF              = ((1-psi)./RR)+psi*(exp_2.*cc);
RISK             = (exp_3./QQ).*cc - SDF.*PREMIUM;
R_impl           = (exp_1.*cc).^(-1);
alp_impl         = ((1-psi)+psi*R_impl.*(exp_2.*cc))./(1-mu);
cons_impl        = (lambda*mu+(1-mu).*(alp_impl))./(exp_4./(qq));

residuals(:,:,1) = log(RR./R_impl)';
residuals(:,:,2) = log(aa./alp_impl)';
residuals(:,:,3) = log(cc./cons_impl)';
residuals(:,:,4) = log((PREMIUM./(((lambda*mu+(1-mu).*(alp_impl))-RISK)./SDF)))';

residuals        = reshape(residuals,N_g*2*4,1);

expectations     = [exp_11',exp_22',exp_33',exp_44',exp_55'];

PREMIUM          = (PREMIUM-RR);
MULT_comp        = R_impl.*(lambda*mu./(alp_impl.*(1-mu)));
RISK_comp        = R_impl.*(RISK./(alp_impl.*(1-mu)));
gdp_ss           = log((0.318^(1-alpha))*(exp(ss(5))*exp(-gamz))^(alpha));

policies.mult         = (mu'/TT)';
policies.premia       = (PREMIUM'/TT)';
policies.risk         = (RISK'/TT)';
policies.gdp          = ((GDP-gdp_ss)'/TT)';
policies.inv          = (INV'/TT)';
policies.ytm          = (YTM'/TT)';
policies.q            = (qq'/TT)';
policies.Q            = (QQ'/TT)';
policies.debt         = ((DEBT-ss(end))'/TT)';
policies.prom         = ((PROM-ss(7))'/TT)';
policies.kap          = ((KAP-ss(5))'/TT)';
policies.mult_comp    = (MULT_comp'/TT)';
policies.risk_comp    = (RISK_comp'/TT)';

end

