%==========================================================================
%                              TABLE A-4
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

path('C:\Users\Luigi\Dropbox\projects\The Pass-Through of Sovereign Risk\Final version_JPE\Replication Files\Matfiles',path);
path('Files',path);

%=========================================================================
%                   Load statistics for price of risk
%=========================================================================

load price_risk

%=========================================================================
%                             TABLE A-4
%=========================================================================

disp('                                                                  ');
disp('                       TABLE A-4                                    ');
disp('                                                                  ');
disp('                    Ergodic      p=.05     p=.05,Low N   Filtered  ');
disp(['var_sdf:             ', num2str(round(var_sdf*100)./100)]);
disp(['E[alp|mu>0]:          ', num2str(round(E_mu*100)./100)]);
disp(['var_alp:             ', num2str(round(100*var_alp*100)./100)]);
disp(['var_cons:            ', num2str(round(100*var_c*100)./100)]);
disp('                                                                  ');
disp('                                                                  ');
