%==========================================================================
%                      Model Solution: Posterior Mean
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
%                       Set Model Parameters 
%=========================================================================

state       = 6;                % number of state variables

order       = [3;3;3;3;3;3];    % order of chebyshev's polynomial

%=========================================================================
%                        Load Starting Guess 
%=========================================================================

load model_nodefault_mean

param(22)  = 0.95;
param(23)  = 0.63;
param(24)  = -7.06;

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

N_g              = size(coll_points,2);

gamma_1          = gamma(1:N_g);
gamma_2          = gamma(N_g+1:2*N_g);
gamma_3          = gamma(2*N_g+1:3*N_g);
gamma_4          = gamma(3*N_g+1:4*N_g);

%=========================================================================
%                         Generate collocation points
%=========================================================================

bounds(6,:)         = [-4.35,4.35];
[coll_points_new s] = smolyakapprox_grid(state,order,bounds(:,1),bounds(:,2));
cons                = zeros(size(coll_points_new,2),1);
R                   = zeros(size(coll_points_new,2),1);
alp                 = zeros(size(coll_points_new,2),1);
q                   = zeros(size(coll_points_new,2),1);

%=========================================================================
%        Use solution to model without default as initial guess 
%=========================================================================

for j=1:size(coll_points_new,2)

    [a,b]   = min(sum(abs(coll_points_new(1:5,j)*ones(1,size(coll_points,2))...
                  -coll_points),1));    
    cons(j) = gamma_1(b);
    R(j)    = gamma_2(b);
    alp(j)  = gamma_3(b);
    q(j)    = gamma_4(b);

end

gamma   = [cons;R;alp;q];
gamma   = [gamma;gamma]';

coll_points  = coll_points_new;

%=========================================================================
%             Generate coefficient matrix for expectations    
%=========================================================================

[EXP]   = precomp_integral_def(coll_points,bounds,rhog,sigmag,rhoz,...
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

rec    = [0,0.3,0.35,0.4,0.45,0.5,0.55];

tic

for j=1:7
  
fprintf('Iteration: %4.0f of %4.0f\n', j,7)

param(25)     = rec(j);
obj           = @(Theta) residual_model(param,Theta,bounds,coll_points,TT,s,ss,EXP,V,2);
[gamma a]     = parsolve(obj,gamma,1);
disp('                                                                  ');

if a==0
    gamma=NaN;
break
end

end

eps = linspace(1.8,1,5);

for j=1:5
    fprintf('Tapering: %4.0f of %4.0f\n', j,5)    
obj           = @(Theta) residual_model(param,Theta,bounds,coll_points,TT,s,ss,EXP,V,eps(j));
[gamma1 a]    = parsolve(obj,gamma,1);

if a==1
    gamma = gamma1;
    rec   = eps(j);
else
    break
end
end
time=toc;
fprintf('Time for model solution: %4.3f minutes\n', (time/60))

%=========================================================================
%                             Save Results
%=========================================================================

save Matfiles/model_solution_mean param s gamma bounds coll_points EXP ss TT V rec

