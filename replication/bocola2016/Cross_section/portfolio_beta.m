%==========================================================================
%                   CONSTRUCT POTFOLIOS SORTED BY BETA
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

path('Data',path);
path('Matfiles',path);

%=========================================================================
%                              LOAD DATA
%=========================================================================

[S S1] = xlsread('compustat_99-11.xlsx');      %read the data from excel

%=========================================================================
%                             CLEAN DATA
%=========================================================================

comp_id       = S(:,1);
sic           = int32(S(:,4)/100);
time_v        = S(:,5);
adj_price     = S(:,6);
share_out     = S(:,7);
volume        = S(:,8);
stock_price   = S(:,9);
return_factor = S(:,10);

company_id    = 0;

time          = zeros(100,800);
prices        = zeros(100,100);
adj           = zeros(100,100);
return_f      = zeros(100,100);
vol           = zeros(100,100);
siz           = zeros(100,100);

l = 1;
k = 1;

for j=1:size(time_v,1)-1;
    
    if comp_id(j+1)==comp_id(j) && time_v(j)>time(max(l-1,1),k)      
        
        prices(l,k)   = stock_price(j);
        adj(l,k)      = adj_price(j);
        return_f(l,k) = return_factor(j);
        vol(l,k)      = volume(j);
        siz(l,k)      = share_out(j)*stock_price(j);
        codes(l,k)    = sic(j);
        time(l,k)     = time_v(j);       
        l = l+1;

    elseif comp_id(j+1)>comp_id(j) %last month
        prices(l,k)     = stock_price(j);
        adj(l,k)        = adj_price(j);
        return_f(l,k)   = return_factor(j);
        time(l,k)       = time_v(j);
        vol(l,k)        = volume(j);
        siz(l,k)        = share_out(j)*stock_price(j);
        codes(l,k)      = sic(j);
        company_name(k) = S1(j,3);
        company_id(k)   = comp_id(j);
        l=1;
        k = k+1;
    end
      
end

time = time(1:end,1:k);

prices(isnan(prices))     =0;
adj(isnan(adj))           =0;
return_f(isnan(return_f)) =0;
vol(isnan(vol))           =0;
siz(isnan(siz))           =0;
codes(isnan(codes))       =0;

%=========================================================================
% COMPUTE RETURNS FOR FIRMS THAT ARE IN THE PANEL FOR THE PERIOD 2003-2011
%=========================================================================

T                   = size(time,1);

prices_monthly      = zeros(T,10);
adj_monthly         = zeros(T,10);
retf_monthly        = zeros(T,10);
vol_monthly         = zeros(T,10);
size_monthly        = zeros(T,10);
sic_monthly         = zeros(T,10);

l=1;

for k=1:size(prices,2)
    
    initial_t = time(1,k);
    final_t   = max(time(:,k));
    
    if initial_t<=20030201 && final_t>=20111131
        
       [a b] = min(abs(initial_t-time(:,4)));
       [e f] = min(abs(final_t-time(:,4)));      
       [c d] = max(time(:,4));
       
      if (f-b)==d-1
       company_sample(l)     = company_name(k);
       prices_monthly(b:f,l) = prices(1:d,k);
       adj_monthly(b:f,l)    = adj(1:d,k);
       retf_monthly(b:f,l)   = return_f(1:d,k);
       vol_monthly(b:f,l)    = vol(1:d,k);
       size_monthly(b:f,l)   = siz(1:d,k);
       codes_monthly(b:f,l)  = codes(1:d,k);
       l     = l+1;
      end
       
    end
    
end

returns_monthly            = zeros(T-1,10);
volume                     = nanmedian(vol_monthly,1);

l = 0;
g = 0;

for k=1:size(prices_monthly,2)
    
    if max(prices_monthly(:,k))>0 && volume(k)>prctile(volume,10)

        
       l = l+1;
       size_monthly(:,l)     = size_monthly(:,k);
       codes_monthly(:,l)    = codes_monthly(:,k);
       
    for j=2:T
        
        if  prices_monthly(min(j-1,216),k)>0 && prices_monthly(min(j,216),k)>0             
            returns_monthly(j-1,l) = log(((prices_monthly(j,k)/adj_monthly(j,k))*...
            retf_monthly(j,k))./((prices_monthly(j-1,k)/adj_monthly(j-1,k))*retf_monthly(j-1,k)));         
        end
        
    end
    end
    
end

size_monthly    = size_monthly(:,1:l);
codes_monthly   = codes_monthly(:,1:l);
returns_monthly = 1+[returns_monthly];
returns_monthly = [ones(1,size(returns_monthly,2));returns_monthly];

%=========================================================================
%                  AVERAGE RETURNS OVER THE QUARTER
%=========================================================================

[T,M]   = size(returns_monthly);
T       = T/3;
returns_quarterly = zeros(T,M);
size_quarterly    = zeros(T,M);

for j=1:int32(T)
    returns_quarterly(j,:) = prod(returns_monthly(1+(j-1)*3:3+(j-1)*3,:),1);
    size_quarterly(j,:)    = mean(size_monthly(1+(j-1)*3:3+(j-1)*3,:),1);
end

market_return    = zeros(T,1);

% value-weighted market return

for j=1:T
    
market_return(j) = (size_quarterly(j,:)./sum(size_quarterly(j,:)))*returns_quarterly(j,:)';    
        
end

%=========================================================================
%                          COMPUTE BETAS
%=========================================================================

load sdf
load fama_french

sdf_data       = sdf;

time_quartsdf  = (1999:.25:2012.75)';
time_quartret  = (1999:.25:2011.75)';
time_quartff   = (1999:.25:2012.75)';

[a first]      = min(abs(1999-time_quartsdf));
[a last]       = min(abs(2011.5-time_quartsdf));

[a first_ret]  = min(abs(1999-time_quartret));
[a last_ret]   = min(abs(2011.5-time_quartret));

[a first_ff]  = min(abs(1999-time_quartff));
[a last_ff]   = min(abs(2011.5-time_quartff));

% Model Stochastic discount factor

sdf            = log(sdf_data(first+1:last+1,2));
sdf_nolev      = log(sdf_data(first+1:last+1,4));

% Fama French Factors

m_ret          = market_return(first_ret:last_ret);
ff_1           = factors(first_ff:last_ff,2);
ff_2           = factors(first_ff:last_ff,3);

% Compute Betas of Benchmark Model

T              = size(sdf,1);

X              = [ones(T-1,1),sdf(1:end-1)];
Y              = sdf(2:end);
Phi            = inv(X'*X)*X'*Y;
sdf_inn        = Y-X*Phi;
X              = [ones(T-1,1),sdf_inn];
Y              = returns_quarterly(first_ret+1:last_ret,:)-1;
Phi            = inv(X'*X)*X'*Y;
beta_lev       = Phi(2,:)';

% Compute Betas of CAPM

X              = [ones(T,1),m_ret];
Y              = returns_quarterly(first_ret:last_ret,:)-1;
Phi            = inv(X'*X)*X'*Y;
beta_capm      = Phi(2,:)';

% Compute Betas of No Lev Specification

X              = [ones(T-1,1),sdf_nolev(1:end-1)];
Y              = sdf_nolev(2:end);
Phi            = inv(X'*X)*X'*Y;
sdf_inn_nolev  = Y-X*Phi;
X              = [ones(T-1,1),sdf_inn_nolev];
Y              = returns_quarterly(first_ret+1:last_ret,:)-1;
Phi            = inv(X'*X)*X'*Y;
beta_nolev     = Phi(2,:)';


%=========================================================================
%                       CONSTRUCT PORTFOLIOS
%=========================================================================

% Sort Benchmark 

[betas_sort,index] = sort(beta_lev);        % Sort by beta
returns_sorted     = returns_quarterly(:,index);
size_sorted        = size_quarterly(:,index);


N_port  = 10;     % Number of Portfolios

ret_lev = zeros(size(returns_sorted,1),N_port);
M       = size(betas_sort,1);
m       = int32(M/N_port);
k       = 0;

for j=1:N_port     % Aggregate stocks and compute quarterly value-weighted returns

    k = k+1;
    
weights        = mean(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),1)./mean((...
                 sum(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),2)*ones(1,size(betas_sort(1+m*(j-1):min(m+m*(j-1),M)),1))),1);    
ret_lev(:,k)   = (weights*returns_sorted(:,1+m*(j-1):min(m+m*(j-1),M))')';

end

% Sort CAPM 

[betas_sort,index] = sort(beta_capm);        % Sort by beta
returns_sorted     = returns_quarterly(:,index);
size_sorted        = size_quarterly(:,index);



ret_capm = zeros(size(returns_sorted,1),N_port);
M        = size(betas_sort,1);
m        = int32(M/N_port);
k        = 0;

for j=1:N_port     % Aggregate stocks and compute quarterly value-weighted returns

    k = k+1;
    
weights        = mean(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),1)./mean((...
                 sum(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),2)*ones(1,size(betas_sort(1+m*(j-1):min(m+m*(j-1),M)),1))),1);    
ret_capm(:,k)   = (weights*returns_sorted(:,1+m*(j-1):min(m+m*(j-1),M))')';

end

% Sort No leverage 

[betas_sort,index] = sort(beta_nolev);        % Sort by beta
returns_sorted     = returns_quarterly(:,index);
size_sorted        = size_quarterly(:,index);


ret_nolev = zeros(size(returns_sorted,1),N_port);
M         = size(betas_sort,1);
m         = int32(M/N_port);
k         = 0;

for j=1:N_port     % Aggreagate stocks and compute quarterly value-weighted returns

    k = k+1;
    
weights        = mean(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),1)./mean((...
                 sum(size_sorted(:,1+m*(j-1):min(m+m*(j-1),M)),2)*ones(1,size(betas_sort(1+m*(j-1):min(m+m*(j-1),M)),1))),1);    
ret_nolev(:,k) = (weights*returns_sorted(:,1+m*(j-1):min(m+m*(j-1),M))')';

end

%=========================================================================
%                       SAVE PORTFOLIOS
%=========================================================================

save Matfiles/portfolios_beta ret_lev ret_capm ret_nolev



