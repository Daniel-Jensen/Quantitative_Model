%==========================================================================
%                      Model Solution: Posterior Draws
%                         
% 
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 12/07/2015
%==========================================================================

%=========================================================================
%                              Housekeeping
%=========================================================================

clear all
close all
clc
delete *.asv

path('Smolyak Files',path);
path('Simulation Files',path);
path('Solution Files',path);
path('Matfiles',path);

%=========================================================================
%                       Set model parameters
%=========================================================================

state       = 6;                % number of state variables

order       = [3;3;3;3;3;3];    % order of chebyshev's polynomial

%=========================================================================
%                        Load starting guess 
%=========================================================================

load model_nodefault_posterior

load param_defprocess

param_step1 = param;

%=========================================================================
%                           Define Objects
%=========================================================================

N_gnodef           = size(coll_points,2);
Ndraw              = size(param_step1,2);
Ncoeff             = size(model_nodef,1);

coll_points_nodef  = coll_points;
param_step1        = param;

param_est         = zeros(25,Ndraw);
gamma_est         = zeros(3112,Ndraw);
EXP_est           = zeros(389,389,Ndraw);
ss_est            = zeros(9,Ndraw);
haircut           = zeros(Ndraw,1);
rec               = zeros(Ndraw,1);

disp('                                                                  ');
disp('          Solve Model for different parametrizations'              );
disp('                                                                  ');

for k=1:Ndraw

disp('                                                                  ');
fprintf('Draw:      %4.0f of %4.0f\n', k, size(param_step1,2))
disp('                                                                  ');

[param,ss]    = model_param(param_step1(:,k));
param(22:23)  = param_step2(2:3,k);
param(24)     = param_step2(1,k);

rhoz       = param(6);   
sigmaz     = param(7);   
g_star     = param(15);
rhog       = param(16);
sigmag     = param(17);
rhos       = param(22);
sigmas     = param(23);

%=========================================================================
%                        Arrange initial guess
%=========================================================================

gamma    = model_nodef(:,k);

gamma_1  = gamma(1:N_gnodef);
gamma_2  = gamma(N_gnodef+1:2*N_gnodef);
gamma_3  = gamma(2*N_gnodef+1:3*N_gnodef);
gamma_4  = gamma(3*N_gnodef+1:4*N_gnodef);

%=========================================================================
%                         Generate collocation points
%=========================================================================

bounds(6,:)         = [-4.75,4.75];
[coll_points_new s] = smolyakapprox_grid(state,order,bounds(:,1),bounds(:,2));
cons                = zeros(size(coll_points_new,2),1);
R                   = zeros(size(coll_points_new,2),1);
alp                 = zeros(size(coll_points_new,2),1);
q                   = zeros(size(coll_points_new,2),1);

%=========================================================================
%        Use solution to model without default as initial guess 
%=========================================================================


for j=1:size(coll_points_new,2)

    [a,b]   = min(sum(abs(coll_points_new(1:5,j)*ones(1,size(coll_points_nodef,2))...
                  -coll_points_nodef),1));    
    cons(j) = gamma_1(b);
    R(j)    = gamma_2(b);
    alp(j)  = gamma_3(b);
    q(j)    = gamma_4(b);

end

coll_points  = coll_points_new;

gamma   = [cons;R;alp;q];
gamma   = [gamma;gamma]';

%=========================================================================
%             Generate coefficient matrix for expectations    
%=========================================================================

[EXP]        = precomp_integral_def(coll_points,bounds,rhog,sigmag,rhoz,...
               sigmaz,rhos,sigmas,s);

%=========================================================================
%                        Generate polynomials    
%=========================================================================

N_g     = size(coll_points,2);
x       = coll_points;
y       = (2*x-(bounds(:,1)+bounds(:,2))*ones(1,N_g))./((bounds(:,2)-...
           bounds(:,1))*ones(1,N_g));
TT      = T(y,s);
      
%=========================================================================
%                            Solve model
%=========================================================================

disp('                                                                  ');
disp('                       Model Solution'                             );
disp('                                                                  ');


y       = [0,0.15,0.3,0.35,0.4,0.45,0.5,0.55];
x       = 1;
rec(k)  = 0;
tic

for j=1:8
  
    l = 0;
    
while l==0
    
fprintf('Iteration: %4.0f of %4.0f\n', j,8)

param(25)     = y(j);
obj           = @(Theta) residual_model(param,Theta,bounds,coll_points,TT,s,ss,EXP,V,2);
[gamma1 a]    = parsolve_post(obj,gamma,x);
disp('                                                                  ');

if a==1
    gamma  = gamma1;
    rec(k) = y(j);
    l      = 1;
else
    x = max(x-0.1,0.05);
end

end

end

eps    = linspace(1.8,1,5);
smooth = 2;
x      = 1;

for j=1:5
    
    l = 0;
    
while l==0
        
    fprintf('Tapering: %4.0f of %4.0f\n', j,5)    
obj           = @(Theta) residual_model(param,Theta,bounds,coll_points,TT,s,ss,EXP,V,eps(j));
[gamma1 b]    = parsolve_post(obj,gamma,x);

if b==1
    gamma   = gamma1;
    smooth  = eps(j);
    l       = 1;
else
    x = max(x-0.1,0.05);
end
end
end
    
time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))

param_est(:,k)    = param;
gamma_est(:,k)    = gamma';
EXP_est(:,:,k)    = EXP;
ss_est(:,k)       = ss;
haircut(k)        = smooth;

%=========================================================================
%                             Save Results
%=========================================================================

save Matfiles/model_solution_posterior param_est s gamma_est bounds coll_points...
     EXP_est ss_est TT V k haircut rec

end

 
