function [residuals] = euler_errors(param,gamma,bounds,TT,EXP,coll_points,ss,expectations,prem,s,V)

N_g      = size(coll_points,2);
alpha    = param(1);
delta    = param(2);
beta     = param(3);
nu       = param(4); 
gamz     = param(5);
rhoz     = param(6);   
chi      = param(8);
psi      = param(9);
omega    = param(10); 
lambda   = param(11);
csi      = param(12);     
a1       = param(13);
a2       = param(14);
g_star   = param(15);
rhog     = param(16);
pi       = param(18);
iota     = param(19);
t_star   = param(20);
gamma_t  = param(21);
rhos     = param(22);
s_star   = param(24);
rec      = param(25);
def      = [0,1];

sigmaz   = param(7);   
sigmag   = param(17);
sigmas   = param(23);

residuals = zeros(2,N_g,4);

NN = 389;

for jj=1:size(coll_points,2)
    
    for j=1:2


x       =  min(max(coll_points(:,jj),bounds(:,1)),bounds(:,2));
y       =  (2*x-(bounds(:,1)+bounds(:,2)))./((bounds(:,2)-...
            bounds(:,1)));
y       = min(max(y,-1),1);      
xx      = (([coll_points(1,jj);coll_points(3,jj)])'*inv(V))';
x(1,:)  = xx(1,:);
x(3,:)  = xx(2,:);
x(2,:)  = x(2,:)+gamz;
gam     = gamma(1+(j-1)*NN*4:NN*4+(j-1)*NN*4);    
gamma_1 = gam(1:NN)/TT;
gamma_2 = gam(NN+1:2*NN)/TT;
gamma_3 = gam(2*NN+1:3*NN)/TT;
gamma_4 = gam(3*NN+1:4*NN)/TT;
        
TTT     = T(y,s);

cons    =  exp(ss(1)+gamma_1*TTT);
R       =  exp(ss(2)+gamma_2*TTT);
alp     =  exp(ss(3)+gamma_3*TTT);
q       =  exp(ss(4)+gamma_4*TTT);

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

x_tom   = [X_tilde(1,:);rhoz*(x(2,:)-gamz);X_tilde(2,:);rhog*x(4,:);...
          log(B_tom)-ss(end);rhos*x(6,:)];

y_tom   = (2*x_tom-(bounds(:,1)+bounds(:,2)))./((bounds(:,2)-...
           bounds(:,1)));

y_tom   = max(y_tom,-1);
y_tom   = min(y_tom,1);

TT_exp  = T(y_tom,s);

[EXP]   = precomp_integral_def(coll_points(:,jj),bounds,rhog,sigmag,rhoz,sigmaz,rhos,sigmas,s);

exp_1 = expectations(:,1:2);
exp_2 = expectations(:,3:4);
exp_3 = expectations(:,5:6);
exp_4 = expectations(:,7:8);

gam_1(:,j)   = ((exp_1(:,j)'/TT)*(TT_exp.*EXP(:,1)))';
gam_2(:,j)   = ((exp_2(:,j)'/TT)*(TT_exp.*EXP(:,1)))';
gam_3(:,j)   = ((exp_3(:,j)'/TT)*(TT_exp.*EXP(:,1)))';
gam_4(:,j)   = ((exp_4(:,j)'/TT)*(TT_exp.*EXP(:,1)))';

PREMIUM(:,j) = ((prem(j,:)/TT)*(TT_exp.*EXP))';

cc(:,j)      = cons';
qq(:,j)      = q';
RR(:,j)      = R';
aa(:,j)      = alp';
mu(:,j)      = (N_tom./(lambda*(Q.*K_tom+q.*B_tom)));
QQ(:,j)      = Q';

    end

Sprob            = (exp(s_star+coll_points(end,jj))./(1+exp(s_star+coll_points(end,jj))))';

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

residuals(:,jj,1) = log10(abs(1-RR./R_impl))';
residuals(:,jj,2) = log10(abs(1-aa./alp_impl))';
residuals(:,jj,3) = log10(abs(1-cc./cons_impl))';
residuals(:,jj,4) = log10(abs(PREMIUM-(((lambda*mu+(1-mu).*(alp_impl))-RISK)./SDF)))';

end 

end

