%==========================================================================
%                                FIGURE 7
%                         
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

%=========================================================================
%                             LOAD FILE
%=========================================================================

load irf_bench
load irf_open


premia = +risk_premia+liq_premia(1:28);

%=========================================================================
%                        FIGURE IRFs BENCHMARK
%=========================================================================

figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(2,3,1), plot(100*prob_def(1:27),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Probability of a default','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
grid on
box off


subplot(2,3,2), plot(400*premia(1:27),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
plot(irf_premia_open(1:27)*400,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
legend('Closed Economy','Open Economy')
set(gca,'FontSize',16,'FontWeight','bold');
title('Excess returns','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
ylim([0,1.25])
grid on
box off

subplot(2,3,3), plot(inv(1:27)*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
plot(irf_inv_open(1:27)*100,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Investment','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
ylim([-2.5,0.5])
grid on
box off

subplot(2,3,4), plot(gdp(1:27)*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
plot(irf_gdp_open(1:27)*100,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Output','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
ylim([-0.3,0.1])
grid on
box off

subplot(2,3,5), plot(lab(1:27)*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on
plot(irf_lab_open(1:27)*100,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Labor','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
ylim([-0.4,0.1])
grid on
box off

subplot(2,3,6), plot(cons(1:27)*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
plot(irf_cons_open(1:27)*100,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Consumption','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([0,30])
ylim([-0.4,0.4])
grid on
box off
