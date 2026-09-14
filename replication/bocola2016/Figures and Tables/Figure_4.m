%==========================================================================
%                             FIGURE 4 
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

load decomposition

%=========================================================================
%                Densities: price and quantity of risk
%=========================================================================

[k1 k2]             = ksdensity(kret_nos(:,1));
[k1_def k2_def]     = ksdensity(kret_s(:,1),'width',0.005);

[sdf1 sdf2]         = ksdensity((sdf_nos(:,1)./mean(sdf_nos(:,1))));
[sdf1_def sdf2_def] = ksdensity(sdf_s(:,1)./(mean(sdf_s(:,1))));

%=========================================================================
%                             Figure 4
%=========================================================================

figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(1,2,1), plot(sdf2,sdf1/sum(sdf1),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on
plot(sdf2_def,sdf1_def/sum(sdf1_def),'LineWidth',3,'Marker','o','Color',[0.66,0,0.01]), hold on
legend('Low $s_{t}$', 'High $s_{t}$')
set(gca,'FontSize',16,'FontWeight','bold');
title('Density of $\frac{\hat{\Lambda}_{t,t+1}}{\mathbf{E}_{t}[\hat{\Lambda}_{t,t+1}]}$','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0.8 1.8])
grid on
box off

subplot(1,2,2), plot(k2,k1./sum(k1),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on
plot(k2_def,k1_def./sum(k1_def),'LineWidth',3,'Marker','o','Color',[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
title('Density of $R_{K,t+1}$','FontSize',19,'FontWeight','bold','interpreter','latex')
xlim([0.9 1.05])
grid on
box off

