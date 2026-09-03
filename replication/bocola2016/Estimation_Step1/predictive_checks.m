%==========================================================================
%                           Predictive Checks
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
%                       Set Model Parameters
%=========================================================================

state       = 5;                % Number of state variables

order       = [3;3;3;3;3];      % Order of Chebyshev's polynomial 
                               

%=========================================================================
%                       Load Posterior Draws  
%=========================================================================

load model_nodefault
load param_step1

Thetasim             = Thetasim(:,5000:end);
model_sol            = model_sol(5000:end,:);

%=========================================================================
%                          Simulate Model
%=========================================================================

load data

data = [data(:,3),data(:,2)];

Time       = 100;
N          = 100;
steps      = int32(size(Thetasim,2)/N)-1;
counter    = 0;

multiplier = zeros(Time-1,N);
output     = zeros(Time-1,N);

disp('                                                                  ');
disp('                                                                  ');    
disp('                   Starting Simulations                           ');
disp('                                                                  ');    
disp('                                                                  ');    

% Initialize


parfor j=1:N
    
    randn('state',j^(2));   % Set seed for simulations

    tic
    
 [param,ss] = model_param(Thetasim(:,j*steps));
    
    rhoz    = param(6);
    sigmaz  = param(7);
    rhog    = param(16);
    sigmag  = param(17);

    EXP     = precomp_integral(coll_points,bounds,rhog,sigmag,rhoz,sigmaz,s);

    gam     = model_sol(j*steps,:);
    
policies    = model_nodefault_policies(param,gam,bounds,coll_points,TT,ss,s,EXP,V);

e           = randn(Time+4,2);

initial     = zeros(5,1);

[state,obs] = simul_nodefault(param,gam',bounds,policies,TT,ss,s,V,e,initial);

output(:,j)       = obs.gdp_growth(end-Time+1:end-1);
multiplier(:,j)   = obs.mult(end-Time+1:end-1);

time = toc;

end

%=========================================================================
%                   Construct Predictive Checks
%=========================================================================

mean_mult     = mean(multiplier,1);
stdev_mult    = var(multiplier).^(1/2);
skew_mult     = skewness(multiplier);
acorr_mult    = zeros(N,1);
cycl_mult     = zeros(N,1);
kurt_mult     = kurtosis(multiplier);

mean_gdp      = mean(output,1);
stdev_gdp     = var(output).^(1/2);
skew_gdp      = skewness(output,1);
acorr_gdp     = zeros(N,1);
kurt_gdp      = kurtosis(output,1);

for j=1:N

acorr_mult(j) = corr(multiplier(2:end,j),multiplier(1:end-1,j));
acorr_gdp(j)  = corr(output(2:end,j),output(1:end-1,j));
cycl_mult(j)  = corr(multiplier(2:end,j),output(2:end,j));

end

%=========================================================================
%                             Save Results
%=========================================================================

save Matfiles/predictive_checks_step1 mean_mult stdev_mult skew_mult... 
     kurt_mult acorr_mult cycl_mult mean_gdp stdev_gdp skew_gdp kurt_gdp acorr_gdp


