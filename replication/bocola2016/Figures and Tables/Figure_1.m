%==========================================================================
%                               FIGURE 1
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

load data_fig1

time = data(:,1);
mult = data(:,2);
gdp  = data(:,3)*100;

%=========================================================================
%                             Figure 1
%=========================================================================

figure('Position',[20,20,900,600],'Name',...
    '','Color','w')

[AX,H1,H2]=plotyy(time,gdp,time,mult);
set(AX,{'ycolor'},{'black';'black'})
set(AX(2),'FontSize',16,'FontWeight','bold','Ylim',[-.0025,0.0125],'Xlim',[2002 2013])
set(AX(1),'FontSize',16,'FontWeight','bold','Xlim',[2002 2013])
set(H1(1),'LineWidth',3,'Linestyle','-','Marker','o','Color',[0.66,0,0.01])
set(H2(1),'LineWidth',3,'Color',[0,0.0500,0.7000])
legend('GDP growth','Lagrange multiplier (right axis)')
grid on
box off 

