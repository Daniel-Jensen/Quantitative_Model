%==========================================================================
%                                 TABLE 3
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

path('C:\Users\Luigi\Dropbox\projects\The Pass-Through of Sovereign Risk\Final version_JPE\Replication Files\Matfiles',path);
path('Files',path);

%=========================================================================
%                        Load posterior draws 
%=========================================================================

load predictive_checks_def

%=========================================================================
%                    Compute posterior statistics   
%=========================================================================

pred_median  = [median(def_probabilities),median(median_prob),prctile(median_prob,5),prctile(median_prob,95)]*100;
pred_mean    = [mean(def_probabilities),median(mean_prob),prctile(mean_prob,5),prctile(mean_prob,95)]*100;
pred_stdev   = [var(def_probabilities).^(1/2),median(stdev_prob),prctile(stdev_prob,5),prctile(stdev_prob,95)]*100;
pred_corr    = [corr(def_probabilities(2:end),def_probabilities(1:end-1)),median(corr_prob),prctile(corr_prob,5),prctile(corr_prob,95)];
pred_skew    = [skewness(def_probabilities),median(skew_prob),prctile(skew_prob,5),prctile(skew_prob,95)];
pred_kurt    = [kurtosis(def_probabilities),median(kurt_prob),prctile(kurt_prob,5),prctile(kurt_prob,95)];

%=========================================================================
%                          Table 3   
%=========================================================================

disp('                                                                  ');
disp('                                                                  ');    
disp(['                   TABLE 3']);
disp('                                                                  ');
fprintf('Statistic    data       median    5th pct     95th pct\n')
fprintf('---------    ----       ------    -------     -------\n')
fprintf('Median       %4.2f       %4.2f      %4.2f       %4.2f\n', pred_median)
fprintf('Mean         %4.2f       %4.2f      %4.2f       %4.2f\n', pred_mean)
fprintf('St Dev       %4.2f       %4.2f      %4.2f       %4.2f\n', pred_stdev)
fprintf('Acorr        %4.2f       %4.2f      %4.2f       %4.2f\n', pred_corr)
fprintf('Skew         %4.2f       %4.2f      %4.2f       %4.2f\n', pred_skew)
fprintf('Kurt         %4.2f     %4.2f      %4.2f       %4.2f\n', pred_kurt)


