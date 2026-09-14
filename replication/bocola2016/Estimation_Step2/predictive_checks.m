%==========================================================================
%            PREDICTIVE CHECKS: SOVEREIGN DEFAULT PROBABILITIES
%
% Author: Luigi Bocola     lbocola@northwestern.edu
% Date  : 12/07/2015
%==========================================================================

%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

clear all
close all
clc
delete *.asv

path('Estimation Files',path);

%=========================================================================
%                       LOAD AND PREPARE DATA
%=========================================================================

load Matfiles/def_prob

def_probabilities = physical(6:end,2);

T   = size(def_probabilities,1);

load Matfiles/param_step2

%=========================================================================
%                     POSTERIOR PREDICTIVE CHECKS
%=========================================================================

M    = 1000;

T    = 100;

mean_prob   = zeros(M,1);
stdev_prob  = zeros(M,1);
skew_prob   = zeros(M,1);
kurt_prob   = zeros(M,1);
corr_prob   = zeros(M,1);
median_prob = zeros(M,1);

for m=1:M
    
    Theta = Thetasim(:,m*Nsim/M);

    s     = Theta(1)*ones(T,1);
    
    eps   = randn(T,1);
    
    for t=2:T
        
        s(t) = (1-Theta(2))*Theta(1)+Theta(2)*s(t-1)+Theta(3)*eps(t);
                
    end
        
      prob         = exp(s)./(1+exp(s));
      
      mean_prob(m) = mean(prob);    
      stdev_prob(m)= var(prob)^(1/2); 
      skew_prob(m) = skewness(prob);
      kurt_prob(m) = kurtosis(prob);      
      corr_prob(m) = corr(prob(2:end),prob(1:end-1));      
      median_prob(m) = median(prob);
      
end
      
%=========================================================================
%                              SAVE RESULTS
%=========================================================================
    
save Matfiles/predictive_checks_def def_probabilities mean_prob stdev_prob skew_prob kurt_prob corr_prob median_prob

    






