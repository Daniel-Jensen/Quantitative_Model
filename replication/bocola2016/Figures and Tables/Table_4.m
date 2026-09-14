%==========================================================================
%                                 TABLE 4
%                         
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 12/07/2015
%==========================================================================

%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

clear all
close all
clc
delete *.asv

path('C:\Users\Luigi\Dropbox\projects\The Pass-Through of Sovereign Risk\Final version_JPE\Replication Files\Matfiles',path);
path('Files',path);

%=========================================================================
%                    Load risk premia decomposition
%=========================================================================

load decomposition

%=========================================================================
%                         Compute statistics
%=========================================================================

A               = -cov(sdf_nos(:,1)./mean(sdf_nos(:,1)),kret_nos(:,1));
risk_premia_nos = A(1,2)*400;
A               = -cov(sdf_s(:,1)./mean(sdf_s(:,1)),kret_s(:,1));
risk_premia_s   = A(1,2)*400;

qrisk_nos =  var(kret_nos(:,1))^(1/2);
prisk_nos =  var(sdf_nos(:,1)./mean(sdf_nos(:,1)))^(1/2);
corr_nos  =  corr(kret_nos(:,1),-sdf_nos(:,1)./mean(sdf_nos(:,1)));

qrisk_s   =  var(kret_s(:,1))^(1/2);
prisk_s   =  var(sdf_s(:,1)./mean(sdf_s(:,1)))^(1/2);
corr_s    =  corr(kret_s(:,1),-sdf_s(:,1)./mean(sdf_s(:,1)));

%=========================================================================
%                            Table 4 
%=========================================================================

disp('                                                                  ');
disp('                       TABLE 4                                    ');
disp('                                                                  ');
fprintf('                  Low s_{t}         High s_{t}\n')
fprintf('Risk Premia         %4.3f           %4.3f\n',[risk_premia_nos,risk_premia_s])
fprintf('Price of risk       %4.3f           %4.3f\n',[prisk_nos,prisk_s])
fprintf('Quantity of risk    %4.3f           %4.3f\n',[qrisk_nos,qrisk_s])
fprintf('Correlation         %4.3f           %4.3f\n',[corr_nos,corr_s])
disp('                                                                  ');
disp('                                                                  ');

