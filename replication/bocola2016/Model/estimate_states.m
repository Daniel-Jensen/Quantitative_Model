%==========================================================================
%                 Estimate state vector using particle filter
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
%                         Load model solution
%=========================================================================

load model_solution_mean

%=========================================================================
%                        Obtain policy functions
%=========================================================================

[policies,expectations,residuals]    = model_policies(param,gamma,bounds,...
                                       EXP,TT,coll_points,ss,s,V,rec);

%=========================================================================
%                             Load data
%=========================================================================

load data                                   % Load data for particle filter

quarters        = data(1:end,1);
data            = data(1:end,2:4);
data            = [data(:,2),data(:,1),data(:,3)];

%=========================================================================
%                 Estimate path for model state variables
%=========================================================================

Time           = size(data(end-7:end,1),1);
Nparticle      = 500000;

disp('                                                                  ');
disp('                   Starting Particle Filter'                       );
disp('                                                                  ');

gam = gamma';

[filt_state,observation,distr,N_eff] = particle(data,param,gam,bounds,...
                                       s,V,policies,Nparticle,quarters);

%=========================================================================
%                            Get Filtered States
%=========================================================================

M                 = 1000;
steps             = floor(Nparticle/M);
filtered_states   = zeros(6,M,Time);


for k=1:M
    
    filtered_states(:,k,:) = squeeze(filt_state(:,k*steps,end-7:end));
    
end

%=========================================================================
%                            Save Results 
%=========================================================================

save Matfiles/state_italy filtered_states
