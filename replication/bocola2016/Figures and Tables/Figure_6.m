%==========================================================================
%                               FIGURE 6 
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

load counterfactual

%=========================================================================
%                       Construct variables to plot
%=========================================================================

[Time,M,Nsim,N] = size(kret_s);

risk_premia       = zeros(Time,M,N);
liqui_premia      = zeros(Time,M,N);
quantity_risk_s   = zeros(Time,M,N);
quantity_risk_ns  = zeros(Time,M,N);
price_risk_s      = zeros(Time,M,N);
price_risk_ns     = zeros(Time,M,N);

for j=1:Time
    for m=1:M
        for l=1:N
       AA  = cov(squeeze(kret_s(j,m,:,l)),squeeze(sdf_s(j,m,:,l))./...
            mean(squeeze(sdf_s(j,m,:,l)))); 
       AA1 = cov(squeeze(kret_ns(j,m,:,l)),squeeze(sdf_ns(j,m,:,l))./...
            mean(squeeze(sdf_ns(j,m,:,l)))); 
       risk_premia(j,m,l) = -AA(1,2)+AA1(1,2);
       
       ret = median(squeeze(kret_s(j,m,:,l)-R_s(j,m,:,l))-(squeeze(kret_ns...
             (j,m,:,l)-R_ns(j,m,:,l))));
         
       liqui_premia(j,m,l) = ret-risk_premia(j,m,l);
 
       quantity_risk_s(j,m,l)   = var(squeeze(kret_s(j,m,:,l))).^(1/2);
       quantity_risk_ns(j,m,l)  = var(squeeze(kret_ns(j,m,:,l))).^(1/2);
       price_risk_s(j,m,l)      = var(squeeze(sdf_s(j,m,:,l))./mean(squeeze...
                                  (sdf_s(j,m,:,l)))).^(1/2);
       price_risk_ns(j,m,l)     = var(squeeze(sdf_ns(j,m,:,l))./mean(squeeze...
                                  (sdf_ns(j,m,:,l)))).^(1/2);      
        end
    end
end

%=========================================================================
%                       Prepare Variables
%=========================================================================

liqui       = squeeze(median(liqui_premia,2));
liqui_50    = median(liqui,2);
liqui_5     = prctile(liqui,5,2);
liqui_20    = prctile(liqui,20,2);
liqui_80    = prctile(liqui,80,2);
liqui_95    = prctile(liqui,95,2);

risk        = squeeze(median(risk_premia,2));
risk_50     = median(risk,2);
risk_5      = prctile(risk,5,2);
risk_20     = prctile(risk,20,2);
risk_80     = prctile(risk,80,2);
risk_95     = prctile(risk,95,2);

premia       = squeeze(median(liqui_premia+risk_premia,2));
premia_50    = median(premia,2);
premia_5     = prctile(premia,5,2);
premia_20    = prctile(premia,20,2);
premia_80    = prctile(premia,80,2);
premia_95    = prctile(premia,95,2);

price_risk_s     = squeeze(median(price_risk_s,2));
price_risk_ns    = squeeze(median(price_risk_ns,2));
quantity_risk_s  = squeeze(median(quantity_risk_s,2));
quantity_risk_ns = squeeze(median(quantity_risk_ns,2));

%=========================================================================
%                             Figure 6
%=========================================================================

f1= figure('Position',[20,20,900,600],'Name',...
         '','Color','w');

time = (2010:.25:2011.75)';

subplot(2,3,1),
shadedplot(time',400*premia_95',400*premia_5',[0.9 0.9 0.9],'k'), hold on
shadedplot(time',400*premia_80',400*premia_20',[0.8 0.8 0.8],'k'), hold on
plot(time,premia_50*400,'LineWidth',3,'Color',[0.66,0,0.01]), hold on
axis([2010 2012 0 1])
set(gca,'XTick',[2010:.50:2012])
set(gca,'XTickLabel',['2010:Q1';'2010:Q3';'2011:Q1';'2011:Q3';'2012:Q1';])
title('Excess returns','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',13,'FontWeight','bold');
grid on
box off

subplot(2,3,2),
shadedplot(time',400*liqui_95',400*liqui_5',[0.9 0.9 0.9],'k'), hold on
shadedplot(time',400*liqui_80',400*liqui_20',[0.8 0.8 0.8],'k'), hold on
plot(time,liqui_50*400,'LineWidth',3,'Color',[0.66,0,0.01]), hold on
set(gca,'XTick',[2010:.50:2012])
set(gca,'XTickLabel',['2010:Q1';'2010:Q3';'2011:Q1';'2011:Q3';'2012:Q1';])
title('Liquidity premia','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',13,'FontWeight','bold');
axis([2010 2012 0 .8])
grid on
box off

subplot(2,3,3),
shadedplot(time',400*risk_95',400*risk_5',[0.9 0.9 0.9],'k'), hold on
shadedplot(time',400*risk_80',400*risk_20',[0.8 0.8 0.8],'k'), hold on
set(gca,'XTick',[2010:.50:2012])
set(gca,'XTickLabel',['2010:Q1';'2010:Q3';'2011:Q1';'2011:Q3';'2012:Q1';])
plot(time,risk_50*400,'LineWidth',3,'Color',[0.66,0,0.01]), hold on
title('Risk premia','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',14,'FontWeight','bold');
axis([2010 2012 0 .8])
grid on
box off

subplot(2,3,4),
plot(time(end-7:end),median(price_risk_s,2),'LineWidth',3,'Marker','o','Color',[0.66,0,0.01]), hold on
plot(time(end-7:end),median(price_risk_ns,2),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on
set(gca,'XTick',[2010:.50:2012])
set(gca,'XTickLabel',['2010:Q1';'2010:Q3';'2011:Q1';'2011:Q3';'2012:Q1';])
title('$\sigma_{t}\left[\frac{\hat{\Lambda}_{t,t+1}}{\mathbf{E}_{t}[\hat{\Lambda}_{t,t+1}]}\right]$','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',13,'FontWeight','bold');
axis([2010 2012 0 0.1])
grid on
box off

subplot(2,3,6),
plot(time(end-7:end),median(quantity_risk_s,2),'LineWidth',3,'Marker','o','Color',[0.66,0,0.01]), hold on
plot(time(end-7:end),median(quantity_risk_ns,2),'LineWidth',3,'Color',[0,0.0500,0.7000]), hold on
set(gca,'XTick',[2010:.50:2012])
set(gca,'XTickLabel',['2010:Q1';'2010:Q3';'2011:Q1';'2011:Q3';'2012:Q1';])
title('$\sigma_{t}\left[R_{K,t+1}\right]$','FontSize',19,'FontWeight','bold','Interpreter','Latex');    
set(gca,'FontSize',13,'FontWeight','bold');
axis([2010 2012 0.0025 0.01])
grid on
box off



