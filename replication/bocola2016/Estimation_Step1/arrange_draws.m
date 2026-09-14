%==========================================================================
%                       Arrange Posterior Draws
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
path('Solution Files',path);
path('Simulation Files',path);
path('Estimation Files',path);
path('Matfiles',path);

%=========================================================================
%               Load model solution and posterior draws
%=========================================================================

load model_nodefault

load param_step1

%=========================================================================
%                  Select subset of posterior draws
%=========================================================================

Nsim        = 10;
param       = zeros(6,Nsim);
model_nodef = zeros(size(gamma,2),Nsim);
steps       = floor(size(Thetasim,2)/Nsim);

for j=1:Nsim
    param(:,j)       = Thetasim(:,j*steps);
    model_nodef(:,j) = model_sol(j*steps,:)';
end

%=========================================================================
%                          Save results
%=========================================================================

save Matfiles/model_nodefault_posterior param model_nodef V bounds coll_points 

    
    
    
    