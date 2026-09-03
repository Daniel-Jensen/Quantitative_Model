%==========================================================================
%                                 TABLE 2
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
%                    Load posterior draws: STEP 1  
%=========================================================================

load param_step1 

mean_mu   = mean(100*Thetasim(1,:));
mean_psi  = mean(Thetasim(2,:));
mean_csi  = mean(Thetasim(3,:));
mean_gamz = mean(400*Thetasim(4,:));
mean_rhoz = mean(Thetasim(5,:));
mean_sigz = mean(100*Thetasim(6,:));

lb_mu   = prctile(100*Thetasim(1,:),5,2);
lb_psi  = prctile(Thetasim(2,:),5,2);
lb_csi  = prctile(Thetasim(3,:),5,2);
lb_gamz = prctile(400*Thetasim(4,:),5,2);
lb_rhoz = prctile(Thetasim(5,:),5,2);
lb_sigz = prctile(100*Thetasim(6,:),5,2);

ub_mu   = prctile(100*Thetasim(1,:),95,2);
ub_psi  = prctile(Thetasim(2,:),95,2);
ub_csi  = prctile(Thetasim(3,:),95,2);
ub_gamz = prctile(400*Thetasim(4,:),95,2);
ub_rhoz = prctile(Thetasim(5,:),95,2);
ub_sigz = prctile(100*Thetasim(6,:),95,2);

%=========================================================================
%                    Load posterior draws: STEP 2  
%=========================================================================

load param_step2 

mean_sstar   = mean(Thetasim(1,:),2);
mean_rhos    = mean(Thetasim(2,:),2);
mean_sigmas  = mean(Thetasim(3,:),2);

lb_sstar     = prctile(Thetasim(1,:),5,2);
lb_rhos      = prctile(Thetasim(2,:),5,2);
lb_sigmas    = prctile(Thetasim(3,:),5,2);
       
ub_sstar     = prctile(Thetasim(1,:),95,2);
ub_rhos      = prctile(Thetasim(2,:),95,2);
ub_sigmas    = prctile(Thetasim(3,:),95,2);

%=========================================================================
%                              Table 2  
%=========================================================================
disp('                                                                  ');    
disp(['                 TABLE 2 ']);
disp('                                                                  ');
disp(['                UPPER PANEL']);
disp('                                                                   ');
fprintf('para     mean      5th pct     95th pct\n')
fprintf('----     ------    -------     -------\n')
fprintf('mu_ss    %4.2f      %4.2f       %4.2f\n', mean_mu, lb_mu, ub_mu)
fprintf('psi      %4.2f      %4.2f       %4.2f\n', mean_psi, lb_psi, ub_psi)
fprintf('csi      %4.2f      %4.2f       %4.2f\n', mean_csi, lb_csi, ub_csi)
fprintf('gammaz   %4.2f      %4.2f       %4.2f\n', mean_gamz, lb_gamz, ub_gamz)
fprintf('rhoz     %4.2f      %4.2f       %4.2f\n', mean_rhoz, lb_rhoz, ub_rhoz)
fprintf('sigz     %4.2f      %4.2f       %4.2f\n', mean_sigz, lb_sigz, ub_sigz)
disp('                                                                  ');
disp(['                LOWER PANEL']);
disp('                                                                  ');
fprintf('para     mean      5th pct     95th pct\n')
fprintf('----     ------    -------     -------\n')
fprintf('s_star   %4.2f      %4.2f       %4.2f\n', mean_sstar,lb_sstar,ub_sstar)
fprintf('rhos     %4.2f      %4.2f       %4.2f\n', mean_rhos,lb_rhos,ub_rhos)
fprintf('sigs     %4.2f      %4.2f       %4.2f\n', mean_sigmas,lb_sigmas,ub_sigmas)

