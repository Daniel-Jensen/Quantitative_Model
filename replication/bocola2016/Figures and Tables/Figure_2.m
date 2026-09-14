%==========================================================================
%                               FIGURE 2
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
%                             Load file
%=========================================================================

load predictive_checks_step1

load data_fig1

gdp_data  = data(1:end-4,3);
mult_data = data(1:end-4,2);

%=========================================================================
%                             Figure 2
%=========================================================================

figure('Position',[20,20,900,600],'Name',...
    '','Color','w')

subplot(2,3,1),
notBoxPlot([mean_gdp',mean_mult'],[1 4]), hold on %',parametersmodelWAGE(i,:)',parametersmodelINFL(i,:)',parametersmodelFFR(i,:)']), hold on
scatter([1,4],[mean(gdp_data),mean(mult_data)],70,[0.66,0,0.01],'filled')
title('Mean','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold','xtick',1:4,'xticklabel',{'GDP Growth','','','Multiplier'});
grid on
box off

subplot(2,3,2),
notBoxPlot([stdev_gdp',stdev_mult'],[1 4]), hold on
scatter([1,4],[var(gdp_data).^(1/2),var(mult_data).^(1/2)],70,[0.66,0,0.01],'filled')
title('Standard deviation','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold','xtick',1:4,'xticklabel',{'GDP Growth','','','Multiplier'});
grid on
box off

subplot(2,3,3),
notBoxPlot([acorr_gdp,acorr_mult],[1 4]), hold on
scatter([1,4],[corr(gdp_data(2:end),gdp_data(1:end-1)),corr(mult_data(2:end),mult_data(1:end-1))],70,[0.66,0,0.01],'filled')
title('Autocorrelation','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold','xtick',1:4,'xticklabel',{'GDP Growth','','','Multiplier'});
grid on
box off

subplot(2,3,4),
notBoxPlot([skew_gdp',skew_mult'],[1 4]), hold on
scatter([1,4],[skewness(gdp_data),skewness(mult_data)],70,[0.66,0,0.01],'filled')
title('Skewness','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold','xtick',1:4,'xticklabel',{'GDP Growth','','','Multiplier'});
ylim([-3,4])
grid on
box off

subplot(2,3,5),
notBoxPlot([kurt_gdp',kurt_mult'],[1 4]), hold on
scatter([1,4],[kurtosis(gdp_data),kurtosis(mult_data)],70,[0.66,0,0.01],'filled')
title('Kurtosis','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold','xtick',1:4,'xticklabel',{'GDP Growth','','','Multiplier'});
ylim([2,10])
grid on
box off

subplot(2,3,6),
notBoxPlot([cycl_mult],1), hold on
scatter(1,corr(gdp_data,mult_data),70,[0.66,0,0.01],'filled')
title('Corr ($\mu_{t}$, GDP Growth$_{t}$)','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
xlim([-0.5,2.5])
set(gca,'FontSize',14,'FontWeight','bold','xtick',-0.5:2.5,'xticklabel',{'','','',''});
grid on
box off
