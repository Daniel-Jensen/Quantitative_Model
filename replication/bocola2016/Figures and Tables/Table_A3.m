%==========================================================================
%                            TABLE A-3
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
%             LOAD portfolios, sdf and Fama-French factors
%=========================================================================

load portfolios_beta 
load portfolio_industrysize
load sdf
load fama_french

%=========================================================================
%                                Options
%=========================================================================

reply  = input('Which portfolios? (1 for size/industry, 2 for beta, 3 for both):  ');
disp('                                                                  ');
disp('                                                                  ');
reply2  = input('Last Period? :  ');
disp('                                                                  ');
disp('                                                                  ');

   if reply==3 

       ret_capm  = [returns_quarterly,ret_capm];
       ret_nolev = [returns_quarterly,ret_nolev];
       ret_lev   = [returns_quarterly,ret_lev]; 
       N_port    = size(ret_lev,2);
       
   elseif reply==1
       
       ret_lev   = returns_quarterly; 
       ret_capm  = returns_quarterly;
       ret_nolev = returns_quarterly;
       N_port    = size(ret_lev,2);
   
   elseif reply==2
       
       N_port    = size(ret_lev,2);

   end

   
   [a,last_period] = min(abs(reply2-sdf(:,1)));
       
%=========================================================================
%        Compute betas of portfolios from time series regressions
%=========================================================================

sdf_data       = sdf;

[a first]      = min(abs(1999-sdf_data(:,1)));
[a last]       = min(abs(reply2-sdf_data(:,1)));

% Pricing Kernel

sdf            = log(sdf_data(first+1:last+1,2));
sdf_nolev      = log(sdf_data(first+1:last+1,4));

T              = size(sdf,1);

% Fama French Factors

ff_1           = factors(first:last,4);
ff_2           = factors(first:last,2);
ff_3           = factors(first:last,3);

% Compute Betas of Benchmark Specification

X              = [ones(T-1,1),sdf(1:end-1)];
Y              = sdf(2:end);
Phi            = inv(X'*X)*X'*Y;
sdf_inn        = Y-X*Phi;
X              = [ones(T-1,1),sdf_inn];
Y              = (ret_lev(first+1:last,:)-1);
Phi            = inv(X'*X)*X'*Y;
beta_lev       = Phi(2,:)';
res_blev       = Y-X*Phi;

% Compute Betas of No Lev Specification

X              = [ones(T-1,1),sdf_nolev(1:end-1)];
Y              = sdf_nolev(2:end);
Phi            = inv(X'*X)*X'*Y;
sdf_inn_nolev  = Y-X*Phi;
X              = [ones(T-1,1),sdf_inn_nolev];
Y              = (ret_nolev(first+1:last,:)-1);
Phi            = inv(X'*X)*X'*Y;
beta_nolev     = Phi(2,:)';
res_bnolev      = Y-X*Phi;

% Compute Betas of CAPM

X              = [ones(T,1),ff_1];
Y              = (ret_capm(first:last,:)-1);
Phi            = inv(X'*X)*X'*Y;
beta_capm      = Phi(2,:)';
res_bcapm      = Y-X*Phi;

% Compute Betas of 3 Factors FF

X              = [ones(T,1),ff_1,ff_2,ff_3];
Y              = (ret_capm(first:last,:)-1);
Phi            = inv(X'*X)*X'*Y;
beta_ff        = Phi(2:end,:)';
res_bff        = Y-X*Phi;

%=========================================================================
%               Cross-sectional regression: benchmark
%=========================================================================

exc_ret       = (ret_lev(1:last_period,:)-sdf_data(1:last_period,3)*ones(1,N_port))*1; 

T             = size(exc_ret,1); 

% Estimate Benchmark Specification

Y             = mean(exc_ret,1)';
X             = [ones(N_port,1),-beta_lev];
B_bench       = inv(X'*X)*X'*Y;
res_bench     = Y-X*B_bench;
R2_bench      = 1-var(res_bench)/var(Y);

% Shanken Standard Errors

c       = B_bench(2)'*(var(sdf_inn)^(-1))*B_bench(2);
Sig_f   = [0,0;0,var(sdf_inn)];
A       = (inv(X'*X)*X')*cov(res_bench)*(inv(X'*X)*X')';
std_err = ((1/T)*((1+c)*A+Sig_f)).^(1/2);

% Compute Chi^{2} Test Statistic

V       = ((eye(N_port)-X*inv(X'*X)*X')*cov(res_bench)*(eye(N_port)-X*inv(X'*X)*X')*(1+c));
chi_2   = res_bench'*ginv(V)*res_bench;
p_value = chi2cdf(chi_2,N_port-1);

% Report Results

disp('                                                                  ');
disp('                                                                  ');    
disp(['                      BENCHMARK                                  ']);
disp('                                                                  ');
fprintf('         Point Estimate    Standard Error\n')
fprintf('----     -------------     --------------\n')
fprintf('b_0         %4.2f              %4.2f\n', B_bench(1)*400, std_err(1,1)*400)
fprintf('b_1         %4.2f              %4.2f\n', B_bench(2)*400, std_err(2,2)*400)
fprintf('----     -------------     --------------\n')
fprintf('      Diagnostics\n')
fprintf('----     -------------\n')
fprintf('R2          %4.2f\n', R2_bench)
fprintf('MAPE        %4.2f\n', mean(abs(res_bench*400)))
fprintf('chi2        %4.2f\n', chi_2)
fprintf('p-value     %4.2f\n', (1-p_value))

%=========================================================================
%               Cross-sectional regression: CAPM
%=========================================================================

exc_ret        = 1*(ret_capm(1:last_period,:)-sdf_data(1:last_period,3)*ones(1,N_port));

% Estimate CAPM

Y              = mean(exc_ret,1)';
X             = [ones(N_port,1),beta_capm];
B_capm        = inv(X'*X)*X'*Y;
res_capm      = Y-X*B_capm;
R2_capm       = 1-var(res_capm)/var(Y);

% Shanken Standard Errors

c       = B_capm(2)'*(var(ff_1)^(-1))*B_capm(2);
Sig_f   = [0,0;0,var(ff_1)];
A       = (inv(X'*X)*X')*cov(res_capm)*(inv(X'*X)*X')';
std_err = ((1/T)*((1+c)*A+Sig_f)).^(1/2);

% Compute Chi^{2} Test Statistic

V       = ((eye(N_port)-X*inv(X'*X)*X')*cov(res_capm)*(eye(N_port)-X*inv(X'*X)*X')*(1+c));
chi_2   = res_capm'*ginv(V)*res_capm;
p_value = chi2cdf(chi_2,N_port-1);

% Report Results

disp('                                                                  ');
disp('                                                                  ');    
disp(['                       CAPM                                     ']);
disp('                                                                  ');
fprintf('         Point Estimate    Standard Error\n')
fprintf('----     -------------     --------------\n')
fprintf('b_0         %4.2f              %4.2f\n', B_capm(1)*400, std_err(1,1)*400)
fprintf('b_1         %4.2f              %4.2f\n', B_capm(2)*400, std_err(2,2)*400)
fprintf('----     -------------     --------------\n')
fprintf('      Diagnostics\n')
fprintf('----     -------------\n')
fprintf('R2          %4.2f\n', R2_capm)
fprintf('MAPE        %4.2f\n', mean(abs(res_capm*400)))
fprintf('chi2        %4.2f\n', chi_2)
fprintf('p-value     %4.2f\n', (1-p_value))

%=========================================================================
%               Cross-sectional regression: 3 factors FF
%=========================================================================

% Estimate FF

X             = [ones(N_port,1),beta_ff];
B_ff          = inv(X'*X)*X'*Y;
res_ff        = Y-X*B_ff;
R2_ff         = 1-var(Y-X*B_ff)/var(Y);

% Shanken Standard Errors

Sig_f   = cov([ff_1,ff_2,ff_3]);
c       = B_ff(2:end)'*(inv(Sig_f))*B_ff(2:end);
Sig_f2  = [zeros(3,1),Sig_f];
Sig_f2  = [zeros(1,4);Sig_f2];
A       = (inv(X'*X)*X')*cov(res_ff)*(inv(X'*X)*X')';
std_err = diag(((1/T)*((1+c)*A+Sig_f2)).^(1/2));

% Compute Chi^{2} Test Statistic

V       = ((eye(N_port)-X*inv(X'*X)*X')*cov(res_ff)*(eye(N_port)-X*inv(X'*X)*X')*(1+c));
chi_2   = res_ff'*ginv(V)*res_ff;
p_value = chi2cdf(chi_2,N_port-3);

% Report Results

disp('                                                                  ');    
disp(['                       Three-factor FF                           ']);
disp('                                                                  ');
fprintf('         Point Estimate    Standard Error\n')
fprintf('----     -------------     --------------\n')
fprintf('b_0         %4.2f              %4.2f\n', B_ff(1)*400, std_err(1)*400)
fprintf('b_1         %4.2f              %4.2f\n', B_ff(2)*400, std_err(2)*400)
fprintf('b_2         %4.2f              %4.2f\n', B_ff(3)*400, std_err(3)*400)
fprintf('b_3         %4.2f              %4.2f\n', B_ff(4)*400, std_err(4)*400)
fprintf('----     -------------     --------------\n')
fprintf('      Diagnostics\n')
fprintf('----     -------------\n')
fprintf('R2          %4.2f\n', R2_ff)
fprintf('MAPE        %4.2f\n', mean(abs(res_ff*400)))
fprintf('chi2        %4.2f\n', chi_2)
fprintf('p-value     %4.2f\n', (1-p_value))

%=========================================================================
%               Cross-sectional regression: No Lev
%=========================================================================

exc_ret        = 1*(ret_nolev(1:last_period,:)-sdf_data(1:last_period,3)*ones(1,N_port));

% Estimate No Lev Specification

Y             = mean(exc_ret,1)';
X             = [ones(N_port,1),-beta_nolev];
B_nolev       = inv(X'*X)*X'*Y;
res_nolev     = Y-X*B_nolev;
R2_nolev      = 1-var(res_nolev)/var(Y);

% Shanken Standard Errors

c       = B_nolev(2)'*(var(sdf_inn_nolev)^(-1))*B_nolev(2);
Sig_f   = [0,0;0,var(sdf_inn_nolev)];
A       = (inv(X'*X)*X')*cov(res_nolev)*(inv(X'*X)*X')';
std_err = ((1/T)*((1+c)*A+Sig_f)).^(1/2);

% Compute Chi^{2} Test Statistic

V       = ((eye(N_port)-X*inv(X'*X)*X')*cov(res_nolev)*(eye(N_port)-X*inv(X'*X)*X')*(1+c));
chi_2   = res_nolev'*ginv(V)*res_nolev;
p_value = chi2cdf(chi_2,N_port-1);

% Report Results

disp('                                                                  ');    
disp(['                       No Leverage                              ']);
disp('                                                                  ');
fprintf('         Point Estimate    Standard Error\n')
fprintf('----     -------------     --------------\n')
fprintf('b_0         %4.2f              %4.2f\n', B_nolev(1)*400, std_err(1,1)*400)
fprintf('b_1         %4.2f              %4.2f\n', B_nolev(2)*400, std_err(2,2)*400)
fprintf('----     -------------     --------------\n')
fprintf('      Diagnostics\n')
fprintf('----     -------------\n')
fprintf('R2          %4.2f\n', R2_nolev)
fprintf('MAPE        %4.2f\n', mean(abs(res_nolev*400)))
fprintf('chi2        %4.2f\n', chi_2)
fprintf('p-value     %4.2f\n', (1-p_value))

