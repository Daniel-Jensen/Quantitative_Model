%==========================================================================
%                               Figure 3
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

%=========================================================================
%                             Load File
%=========================================================================

load irf_bench

premia = risk_premia+liq_premia(1:28);

%=========================================================================
%                             Figure 3 
%=========================================================================

figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(2,2,1), plot(100*prob_def,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Probability of a default','FontSize',19,'FontWeight','bold','interpreter','latex')
ylim([0,3])
xlim([0,30])
grid on
box off

subplot(2,2,2), plot(100*bond_price,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Government bond prices','FontSize',19,'FontWeight','bold','interpreter','latex')
ylim([-18,0])
xlim([0,30])
grid on
box off

subplot(2,2,3), plot(100*networth,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Bank net worth','FontSize',19,'FontWeight','bold','interpreter','latex')
ylim([-10,5])
xlim([0,30])
grid on
box off

subplot(2,2,4), plot(400*premia,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
plot(risk_premia*400,'-o','LineWidth',3,'Color',[0.66,0,0.01]), hold on
plot(liq_premia(1:28)*400,'--','LineWidth',3,'Color','k'), hold on
legend('Excess returns','Risk premia','Liquidity premia')
set(gca,'FontSize',16,'FontWeight','bold');
title('Excess returns','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0,30])
ylim([0,0.6000001])
grid on
box off

