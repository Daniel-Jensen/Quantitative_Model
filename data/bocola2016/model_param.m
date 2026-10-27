function [param,ss] = model_param(param)

% This function performs the reparametrization of the model 
%
% Inputs: param (vector collecting the parameters estimated in Step 1)
%
% Output: param (vector collecting all model parameters), ss (vector collecting
% the steady state of the model)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

mu      = param(1);
psi     = param(2);
csi     = param(3);
gamz    = param(4);
rhoz    = param(5);  
sigmaz  = param(6); 
alpha   = 0.3;
nu      = 2;
gamma_t = 1;
g_star  = log(0.198);
rhog    = 0.89;
sigmag  = 0.012;
R       = 1.003;
lev     = 5;
pi      = 0.0546;

beta    = exp(gamz)/R;
alp     = (1-psi)/((1-mu)-psi);
lambda  = alp/lev;
exc_ret = (lambda/alp)*(1/beta)*(mu/(1-mu));
R_k     = R+exc_ret;
R_b     = R_k;
h       = 0.318;
k       = @(x) ((alpha/(R_k-(1-x))*exp((1-alpha)*gamz))^(1/(1-alpha)))*h;
ff      = @(x) (1-x)*exp(-gamz) + (exp(-alpha*gamz)*0.213*(h/k(x))^(1-alpha)-1);
delta   = fzero(ff,0.01);
k       = k(delta);
y       = (h^(1-alpha))*(k*exp(-gamz))^(alpha);
i       = log(k*(1-(1-delta)*exp(-gamz)));
c       = y*(1-exp(g_star))-exp(i);
iota    = ((R_b-pi)/(1-pi))-1;
q       = 1;
 f2     = q-(pi+(1-pi)*(q+iota))*exp(-gamz);
options=optimset('Display','off','MaxIter',10000000,'MaxFunEvals',30000000,'Algorithm',...
                 'levenberg-marquardt','UseParallel','always');   % Option to display output
 f      = @(t) [(f2*t(1)) - (exp(g_star)*exp(log(y))) + (t(1)^(gamma_t) * exp(t(2)));0.076-t(1)./(t(1)+k)];
 b      = fsolve(f,[3;log(0.076)],options);
 t_star = b(2);
 b      = b(1);
t      = t_star + gamma_t*log(b);
W      = ((1-alpha)*y/h);
n      = (k+b)/lev;
chi    = W*(h^(-1/nu))/c;
dz     = gamz;
g      = 0;
h      = log(h);
omega  = n*(1-psi*(exc_ret*lev+R)*exp(-gamz))/((k+b)*exp(-gamz));
k      = log(k);
c      = log(c);
R      = log(R);
n      = log(n);
alp    = log(alp);
b      = log(b);
Q      = 0;
q      = 0;
a1     = (exp(gamz)*(1-(1-delta)*exp(-gamz)))^(csi)/(1-csi);
a2     = exp(gamz)*(1-(1-delta)*exp(-gamz))-a1*(exp(gamz)*(1-(1-delta)*exp(-gamz)))^(1-csi);
P      = log(exp(R).*(exp(Q).*exp(k)+ exp(q).*exp(b)-exp(n)));

param   = [alpha;delta;beta;nu;gamz;rhoz;sigmaz;chi;psi;omega;lambda;csi;a1;a2;...
           g_star;rhog;sigmag;pi;iota;t_star;gamma_t];

ss      = [c;R;alp;q;k;dz;P;g;b];



end
