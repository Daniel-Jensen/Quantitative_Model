%==========================================================================
%                         Compute Euler Errors
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
%                           Load model solution 
%=========================================================================

load model_solution_mean

%=========================================================================
%                 Get policy functions and simulate data
%=========================================================================

randn('state',10);

N                 = 5000;

[policies,expectations,residuals,prem]    = model_policies(param,gamma,bounds,...
                                            EXP,TT,coll_points,ss,s,V,rec);

initial           = zeros(6,1);
l                 = -15*ones(N,1);
e                 = randn(N,3);

[state,obs,STATE] = simul(param,gamma',bounds,policies,TT,ss,s,V,e,l,initial);

%=========================================================================
%                   Compute Euler equation errors
%=========================================================================

coll_points             = STATE;
errors                  = zeros(size(coll_points,2),4);

parfor j=1:size(coll_points,2)
x              = coll_points(:,j);    
[err]          = euler_errors(param,gamma,bounds,TT,EXP,x,ss,expectations,prem,s,V);
errors(j,:)    = err(1,1,:);
end

%=========================================================================
%                              Save results
%=========================================================================

save Matfiles/euler_errors errors

 