%==========================================================================
%                               Figure A-1
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
%                             Load Euler errors
%=========================================================================

load euler_errors

%=========================================================================
%                                 Figure A-1
%=========================================================================

a = sort(errors(:,1));

figure('Position',[20,20,900,600],'Color','w')                                                   

subplot(2,2,1), hist(a,100), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
vline(mean(a),'--r','Mean')
title('Risk Free Rate','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([-8,0])
grid on
box off

a = sort(errors(:,2));

subplot(2,2,2), hist(a,100), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
vline(mean(a),'--r')
title('Marginal Value of Wealth','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([-8,0])
grid on
box off

a = sort(errors(:,3));

subplot(2,2,3), hist(a,100), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
vline(mean(a),'--r')
title('Consumption','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([-8,0])
grid on
box off

a = sort(errors(:,4));

subplot(2,2,4), hist(a,100), hold on ;%[0.66,0,0.01]), hold on
set(gca,'FontSize',16,'FontWeight','bold');
vline(mean(a),'--r')
title('Price of Government Bonds','FontSize',19,'FontWeight','bold','Interpreter','Latex')
xlim([-8,0])
grid on
box off





