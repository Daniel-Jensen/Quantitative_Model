%==========================================================================
%                                 TABLE 5
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
%                            Load output losses
%=========================================================================

load output_losses

%=========================================================================
%                            Construct objects
%=========================================================================

[Time,N,M]    = size(gdp_gs);

output_losses = squeeze(mean(cumsum(gdp_gs-gdp_gns,1),2))*400;

ol_forecasts  = squeeze(mean(cumsum(gdp_for_gs-gdp_for_gns,1),2))*400;

%=========================================================================
%                            Table 5 
%=========================================================================

disp('                                                                  ');
disp('                       TABLE 4                                    ');
disp('                                                                  ');
fprintf('                  Mean         20th Pctl       80th Pctl\n')
fprintf('2011:Q3          %4.2f          %4.2f           %4.2f\n',[mean(output_losses(end-1,:)),prctile(output_losses(end-1,:),20),prctile(output_losses(end-1,:),80)])
fprintf('2011:Q4          %4.2f          %4.2f           %4.2f\n',[mean(output_losses(end,:)),prctile(output_losses(end,:),20),prctile(output_losses(end,:),80)])
fprintf('2012:Q1          %4.2f          %4.2f           %4.2f\n',[mean(ol_forecasts(2,:)),prctile(ol_forecasts(2,:),20),prctile(ol_forecasts(2,:),80)])
fprintf('2012:Q2          %4.2f          %4.2f           %4.2f\n',[mean(ol_forecasts(3,:)),prctile(ol_forecasts(3,:),20),prctile(ol_forecasts(3,:),80)])
fprintf('2012:Q3          %4.2f          %4.2f           %4.2f\n',[mean(ol_forecasts(4,:)),prctile(ol_forecasts(4,:),20),prctile(ol_forecasts(4,:),80)])
fprintf('2012:Q4          %4.2f          %4.2f           %4.2f\n',[mean(ol_forecasts(5,:)),prctile(ol_forecasts(5,:),20),prctile(ol_forecasts(5,:),80)])
disp('                                                                  ');
disp('                                                                  ');
