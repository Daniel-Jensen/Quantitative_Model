%==========================================================================
%                                FIGURE 5
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
%                              Load file
%=========================================================================

load irf_bench

%=========================================================================
%                              Figure 5
%=========================================================================

figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(2,2,1), plot(inv*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
set(gca,'FontSize',16,'FontWeight','bold');
title('Investment','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0,30])
ylim([-2.5,0.500001])
grid on
box off

subplot(2,2,2), plot(gdp*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
set(gca,'FontSize',16,'FontWeight','bold');
title('Output','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0,30])
ylim([-0.300001,0.01])
grid on
box off

subplot(2,2,3), plot(lab*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
set(gca,'FontSize',16,'FontWeight','bold');
title('Labor','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0,30])
ylim([-0.4000001,0.10001])
grid on
box off

subplot(2,2,4), plot(cons*100,'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on 
set(gca,'FontSize',16,'FontWeight','bold');
title('Consumption','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0,30])
ylim([-0.1000001,0.4000001])
grid on
box off
