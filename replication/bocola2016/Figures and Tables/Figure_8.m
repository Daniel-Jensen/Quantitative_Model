%==========================================================================
%                               Figure 8
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
%                             Load file
%=========================================================================

load ltro

%=========================================================================
%                             Figure 8
%=========================================================================


figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(2,2,1),
scatter(delta(1:M2),impact_gdp(1:M2)*100,10,[0,0.0500,0.7000]), hold on
scatter(delta(M2+1:end),impact_gdp(M2+1:end)*100,10,[0.66,0,0.01],'filled','s'), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Output','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
xlim([0 1])
ylim([0 0.7]) 
grid on
box off

subplot(2,2,2),
scatter(delta(1:M2),impact_ret(1:M2)*400,10,[0,0.0500,0.7000]), hold on
scatter(delta(M2+1:end),impact_ret(M2+1:end)*400,10,[0.66,0,0.01],'filled','s'), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Excess returns','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
xlim([0 1])
ylim([-1.5 0.5]) 
grid on
box off

subplot(2,2,3),
scatter(delta(1:M2),impact_liqui(1:M2)*400,10,[0,0.0500,0.7000]), hold on
scatter(delta(M2+1:end),impact_liqui(M2+1:end)*400,10,[0.66,0,0.01],'filled','s'), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Liquidity premium','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
xlim([0 1])
ylim([-1.5 0]) 
grid on
box off

subplot(2,2,4),
scatter(delta(1:M2),impact_risk(1:M2)*400,10,[0,0.0500,0.7000]), hold on
scatter(delta(M2+1:end),impact_risk(M2+1:end)*400,10,[0.66,0,0.01],'filled','s'), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Risk premium','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
xlim([0 1])
ylim([-0.2 0.6]) 
grid on
box off
