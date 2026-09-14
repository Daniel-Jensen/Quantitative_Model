function [residuals] = residual_nodefault(param,gamma,bounds,coll_points,TT,ss,s,EXP,V)

% This function generates the vector of residuals (left hand side of
% equations A.1-A.4 in the Online Appendix) evaluated at the collocation
% points for the model without default risk 
%
% Inputs: param (vector collecting structural parameters), gamma (guess for
% the numerical solution of the model), bounds (matrix collecting upper and
% lower bounds for the endogenous variables in the Smolyak grid),
% coll_points (matrix collecting the collocation points obtained with
% Smolyak method), TT (matrix collecting the Chebyshev's polynomials
% evaluated at the collocation points), ss (deterministic steady state), 
% s (structure used to evaluate Chebyshev's polynomials for an arbitrary 
% point in the state space), EXP (matrix collecting coefficients for integral 
% precomputation), V (matrix that rotate the state variables).
%
% Output: residuals (left hand side of equations A.1-A.4 in the Online
% Appendix evaluated at the collocation points)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

N_g     = size(coll_points,2);

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
residuals = zeros(N_g,4);

xx      = (([coll_points(1,:);coll_points(3,:)])'*inv(V))';

x       = coll_points;
x(1,:)  = xx(1,:);
x(3,:)  = xx(2,:);
x(2,:)  = x(2,:)+gamz;

gamma_1 = gamma(1:N_g);
gamma_2 = gamma(N_g+1:2*N_g);
gamma_3 = gamma(2*N_g+1:3*N_g);
gamma_4 = gamma(3*N_g+1:4*N_g);
        
cons    =  exp(ss(1)+gamma_1);   % Guess for consumption at collocation
R       =  exp(ss(2)+gamma_2);   % Guess for risk free rate 
alp     =  exp(ss(3)+gamma_3);   % Guess for bankers' marginal value of wealth
q       =  exp(ss(4)+gamma_4);   % Guess for bond prices

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

B_tom   =  ((pi+(1-pi)*(q+iota)).*exp(ss(end)+x(5,:)).*exp(-x(2,:))+exp(x(4,:)+...
           g_star).*gdp - exp(t_star+gamma_t*(ss(end)+x(5,:))))./q;

b_ret   =  (pi+(1-pi)*(iota+q));

N_tom   =  max((psi*(k_ret.*exp(ss(5)+x(1,:)) + b_ret.*exp(ss(end)+x(5,:)) -...
           exp(ss(7)+x(3,:))) + omega*(Q.*exp(ss(5)+x(1,:))+q.*exp(ss(end)+...
           x(5,:)))).*exp(-x(2,:)),0.65);

P_tom   =  R.*(Q.*K_tom+q.*B_tom-N_tom); 

%=========================================================================
%                         COMPUTE EXPECTATIONS
%=========================================================================

X       = [log(K_tom)'-ss(5),log(P_tom)'-ss(7)];
X_tilde = (X*V)';

x_tom   = [X_tilde(1,:);rhoz*coll_points(2,:);X_tilde(2,:);rhog*coll_points(4,:);...
          log(B_tom)-ss(end)];

y_tom   = (2*x_tom-(bounds(:,1)+bounds(:,2))*ones(1,N_g))./((bounds(:,2)-...
           bounds(:,1))*ones(1,N_g));

y_tom   = max(y_tom,-1);
y_tom   = min(y_tom,1);

TT_exp  = T(y_tom,s);

exp_1   = log(beta*(exp(-x(2,:))).*(cons.^(-1)));
exp_2   = log(exp(exp_1).*alp);
exp_3   = log(exp(exp_1).*((1-psi)+psi*alp).*k_ret);
exp_4   = log(exp(exp_1).*((1-psi)+psi*alp).*b_ret);
prem    = log((k_ret));

exp_1   = exp(((exp_1/TT)*(TT_exp.*EXP))');
exp_2   = exp(((exp_2/TT)*(TT_exp.*EXP))');
exp_3   = exp(((exp_3/TT)*(TT_exp.*EXP))');
exp_4   = exp(((exp_4/TT)*(TT_exp.*EXP))');
PREMIUM = exp(((prem/TT)*(TT_exp.*EXP))');

%=========================================================================
%                           COMPUTE RESIDUALS
%=========================================================================

mu           = (N_tom./(lambda*(Q.*K_tom+q.*B_tom)))';
mu           = max(1-(((1-psi)+psi*R'.*(cons'.*exp_2)).*mu),0);

PREMIUM      = PREMIUM./Q';
SDF          = ((1-psi)./R')+psi*(cons'.*exp_2);
RISK         = cons'.*(exp_3./Q') - SDF.*PREMIUM;
R_impl       = (cons'.*exp_1).^(-1);
alp_impl     = ((1-psi)+psi*R_impl.*(cons'.*exp_2))./(1-mu);
cons_impl    = (lambda*mu+(1-mu).*(alp_impl))./(exp_4./(q'));

residuals(:,1) = log(cons'./cons_impl);
residuals(:,2) = log(R'./R_impl);
residuals(:,3) = log(alp'./alp_impl);
residuals(:,4) = PREMIUM-(((lambda*mu+(1-mu).*(alp_impl))-RISK)./SDF);

residuals      = reshape(residuals,N_g*4,1);

end

