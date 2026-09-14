function [residuals] = residual_model_ltro(param,gamma,bounds,coll_points,TT,s,ss,EXP,V,expectations,prem,eps)

N_g      = size(coll_points,2);

alpha        = param(1);
delta        = param(2);
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

PREMIUM   = zeros(N_g,2);
mu        = zeros(N_g,2);
mu1       = zeros(N_g,2);
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

exp_1        = expectations(:,j)';
exp_2        = expectations(:,j+2)';
exp_3        = expectations(:,j+4)';
exp_4        = expectations(:,j+6)';
exp_5        = expectations(:,j+8)';

gam_1(:,j)   = ((exp_1/TT)*(TT_exp.*EXP))';
gam_2(:,j)   = ((exp_2/TT)*(TT_exp.*EXP))';
gam_3(:,j)   = ((exp_3/TT)*(TT_exp.*EXP))';
gam_4(:,j)   = ((exp_4/TT)*(TT_exp.*EXP))';
gam_5(:,j)   = ((log(exp(exp_1).*exp(exp_5))/TT)*(TT_exp.*EXP))';

PREMIUM(:,j) = exp(((prem(j,:)/TT)*(TT_exp.*EXP)))';

cc(:,j)      = cons';
qq(:,j)      = q';
RR(:,j)      = R';
aa(:,j)      = alp';
mu(:,j)      = ((N_tom)./(lambda*(Q.*K_tom+q.*B_tom)));
mu1(:,j)     = 1./(lambda*(Q.*K_tom+q.*B_tom));
QQ(:,j)      = Q';

end

Sprob            = (exp(s_star+coll_points(end,:))./(1+exp(s_star+coll_points(end,:))))';

exp_1            = ((1-Sprob).*exp(gam_1(:,1)) + Sprob.*exp(gam_1(:,2)))*ones(1,2);
exp_2            = ((1-Sprob).*exp(gam_2(:,1)) + Sprob.*exp(gam_2(:,2)))*ones(1,2);
exp_3            = ((1-Sprob).*exp(gam_3(:,1)) + Sprob.*exp(gam_3(:,2)))*ones(1,2);
exp_4            = ((1-Sprob).*exp(gam_4(:,1)) + Sprob.*exp(gam_4(:,2)))*ones(1,2);
exp_5            = ((1-Sprob).*exp(gam_5(:,1)) + Sprob.*exp(gam_5(:,2)))*ones(1,2);

PREMIUM          = ((((1-Sprob).*PREMIUM(:,1)+ Sprob.*PREMIUM(:,2))*ones(1,2))./QQ);
SDF              = ((1-psi)./RR)+psi*(exp_2.*cc);
RISK             = (exp_3./QQ).*cc - SDF.*PREMIUM;
R_impl           = (exp_1.*cc).^(-1);

residuals(:,:,1) = log(RR./R_impl)';

mu               = max(1-(((1-psi)+psi*R_impl.*(exp_2./cc.^(-1))).*mu) + psi*(exp_5./cc.^(-1)).*mu1,0);
alp_impl         = ((1-psi)+psi*R_impl.*(exp_2.*cc))./(1-mu);
cons_impl        = (lambda*mu+(1-mu).*(alp_impl))./(exp_4./(qq));

residuals(:,:,2) = log(aa./alp_impl)';
residuals(:,:,3) = log(cc./cons_impl)';
residuals(:,:,4) = log((PREMIUM./(((lambda*mu+(1-mu).*(alp_impl))-RISK)./SDF)))';

residuals        = reshape(residuals,N_g*2*4,1);

end

