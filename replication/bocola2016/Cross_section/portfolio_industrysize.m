%==========================================================================
%               CONSTRUCT POTFOLIOS BY INDUSTRY AND SIZE
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

volume                     = nanmedian(vol_monthly,1);
returns_monthly            = zeros(T-1,10);

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
%                    COMPUTE QUARTERLY MARKET RETURN
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
%                  PARTITION BY TWO DIGIT SIC CODES
%=========================================================================

manu=0;
retail=0;
finance=0;
serv=0;
nat_res=0;

for j=1:size(returns_monthly,2)
    
    if codes_monthly(150,j)<=30
        
        manu = manu+1;
        
        ret_manu(:,manu)  = returns_monthly(:,j);
        size_manu(:,manu) = size_monthly(:,j);        
        comp_manu(manu)   = company_sample(j);
           
    elseif codes_monthly(150,j)>=40 && codes_monthly(150,j)<60 || codes_monthly(150,j)>=70 && codes_monthly(150,j)<90
        serv = serv+1;

        ret_serv(:,serv)  = returns_monthly(:,j);
        size_serv(:,serv) = size_monthly(:,j);
        comp_serv(serv)   = company_sample(j);

    elseif codes_monthly(150,j)>=60 && codes_monthly(150,j)<70 
       
        finance = finance+1;
        ret_fin(:,finance)  = returns_monthly(:,j);
        size_fin(:,finance) = size_monthly(:,j);
        comp_fin(finance)   = company_sample(j);
      
    end
end

%=========================================================================
%                   PARTITION BY SIZE QUINTILES
%=========================================================================

% Manufacturing

ret_manu1 = zeros(size(size_manu,1),1);
ret_manu2 = zeros(size(size_manu,1),1);
ret_manu3 = zeros(size(size_manu,1),1);
ret_manu4 = zeros(size(size_manu,1),1);
ret_manu5 = zeros(size(size_manu,1),1);

size_manu1 = zeros(size(size_manu,1),1);
size_manu2 = zeros(size(size_manu,1),1);
size_manu3 = zeros(size(size_manu,1),1);
size_manu4 = zeros(size(size_manu,1),1);
size_manu5 = zeros(size(size_manu,1),1);

k = 0;
l = 0;
m = 0;
n = 0;
o = 0;

for j=1:size(size_manu,2)
    
    if median(size_manu(:,j),1) <= prctile(median(size_manu,1),20)
   
    k = k+1;    
        
    ret_manu1(:,k)  = ret_manu(:,j);
    size_manu1(:,k) = size_manu(:,j);
        
    elseif median(size_manu(:,j),1) > prctile(median(size_manu,1),20) && median(size_manu(:,j),1) <= prctile(median(size_manu,1),40)

    l = l+1;
    
    ret_manu2(:,l) = ret_manu(:,j);
    size_manu2(:,l) = size_manu(:,j);
    
    elseif median(size_manu(:,j),1) > prctile(median(size_manu,1),40) && median(size_manu(:,j),1) <= prctile(median(size_manu,1),60)

    m = m+1;
    
    ret_manu3(:,m)  = ret_manu(:,j);
    size_manu3(:,m) = size_manu(:,j);

    elseif median(size_manu(:,j),1) > prctile(median(size_manu,1),60) && median(size_manu(:,j),1) <= prctile(median(size_manu,1),80)

    n = n+1;
    
    ret_manu4(:,n)  = ret_manu(:,j);
    size_manu4(:,n) = size_manu(:,j);
    
    else
        
    o = o+1;
    
    ret_manu5(:,o) = ret_manu(:,j);
    size_manu5(:,o) = size_manu(:,j);
        
    end
    
end

% Services

ret_serv1 = zeros(size(size_serv,1),1);
ret_serv2 = zeros(size(size_serv,1),1);
ret_serv3 = zeros(size(size_serv,1),1);
ret_serv4 = zeros(size(size_serv,1),1);
ret_serv5 = zeros(size(size_serv,1),1);

size_serv1 = zeros(size(size_serv,1),1);
size_serv2 = zeros(size(size_serv,1),1);
size_serv3 = zeros(size(size_serv,1),1);
size_serv4 = zeros(size(size_serv,1),1);
size_serv5 = zeros(size(size_serv,1),1);

k = 0;
l = 0;
m = 0;
n = 0;
o = 0;

for j=1:size(size_serv,2)
    
    if median(size_serv(:,j),1) <= prctile(median(size_serv,1),20)
   
    k = k+1;    
        
    ret_serv1(:,k)  = ret_serv(:,j);
    size_serv1(:,k) = size_serv(:,j);
    
    elseif median(size_serv(:,j),1) > prctile(median(size_serv,1),20) && median(size_serv(:,j),1) <= prctile(median(size_serv,1),40)

    l = l+1;
    
    ret_serv2(:,l) = ret_serv(:,j);
    size_serv2(:,l) = size_serv(:,j);
    
    elseif median(size_serv(:,j),1) > prctile(median(size_serv,1),40) && median(size_serv(:,j),1) <= prctile(median(size_serv,1),60)

    m = m+1;
    
    ret_serv3(:,m)  = ret_serv(:,j);
    size_serv3(:,m) = size_serv(:,j);

    elseif median(size_serv(:,j),1) > prctile(median(size_serv,1),60) && median(size_serv(:,j),1) <= prctile(median(size_serv,1),80)

    n = n+1;
    
    ret_serv4(:,n)  = ret_serv(:,j);
    size_serv4(:,n) = size_serv(:,j);
    
    else
        
    o = o+1;
    
    ret_serv5(:,o) = ret_serv(:,j);
    size_serv5(:,o) = size_serv(:,j);
        
    end
    
end

% Finance

ret_fin1 = zeros(size(size_fin,1),1);
ret_fin2 = zeros(size(size_fin,1),1);
ret_fin3 = zeros(size(size_fin,1),1);
ret_fin4 = zeros(size(size_fin,1),1);
ret_fin5 = zeros(size(size_fin,1),1);

size_fin1 = zeros(size(size_fin,1),1);
size_fin2 = zeros(size(size_fin,1),1);
size_fin3 = zeros(size(size_fin,1),1);
size_fin4 = zeros(size(size_fin,1),1);
size_fin5 = zeros(size(size_fin,1),1);

k = 0;
l = 0;
m = 0;
n = 0;
o = 0;

for j=1:size(size_fin,2)
    
    if median(size_fin(:,j),1) <= prctile(median(size_fin,1),20)
   
    k = k+1;    
        
    ret_fin1(:,k)  = ret_fin(:,j);
    size_fin1(:,k) = size_fin(:,j);
        
    elseif median(size_fin(:,j),1) > prctile(median(size_fin,1),20) && median(size_fin(:,j),1) <= prctile(median(size_fin,1),40)

    l = l+1;
    
    ret_fin2(:,l) = ret_fin(:,j);
    size_fin2(:,l) = size_fin(:,j);
    
    elseif median(size_fin(:,j),1) > prctile(median(size_fin,1),40) && median(size_fin(:,j),1) <= prctile(median(size_fin,1),60)

    m = m+1;
    
    ret_fin3(:,m)  = ret_fin(:,j);
    size_fin3(:,m) = size_fin(:,j);

    elseif median(size_fin(:,j),1) > prctile(median(size_fin,1),60) && median(size_fin(:,j),1) <= prctile(median(size_fin,1),80)

    n = n+1;
    
    ret_fin4(:,n)  = ret_fin(:,j);
    size_fin4(:,n) = size_fin(:,j);
    
    else
        
    o = o+1;
    
    ret_fin5(:,o) = ret_fin(:,j);
    size_fin5(:,o) = size_fin(:,j);
        
    end
    
end

%=========================================================================
%                       VALUE WEIGHTED RETURNS
%=========================================================================

returns_fin1        = zeros(T,1);
returns_fin2        = zeros(T,1);
returns_fin3        = zeros(T,1);
returns_fin4        = zeros(T,1);
returns_fin5        = zeros(T,1);

returns_manu1        = zeros(T,1);
returns_manu2        = zeros(T,1);
returns_manu3        = zeros(T,1);
returns_manu4        = zeros(T,1);
returns_manu5        = zeros(T,1);

returns_serv1        = zeros(T,1);
returns_serv2        = zeros(T,1);
returns_serv3        = zeros(T,1);
returns_serv4        = zeros(T,1);
returns_serv5        = zeros(T,1);


k=0;
l=11;

for j=1:size(ret_fin1,1)

    l = l+1;
    
    if l==12
        k = k+1;
        l = 0;
        weights_fin1 = mean(size_fin1(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_fin1(1+12*(k-1):12+12*(k-1),:),2)));
        weights_fin2 = mean(size_fin2(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_fin2(1+12*(k-1):12+12*(k-1),:),2)));
        weights_fin3 = mean(size_fin3(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_fin3(1+12*(k-1):12+12*(k-1),:),2)));
        weights_fin4 = mean(size_fin4(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_fin4(1+12*(k-1):12+12*(k-1),:),2)));
        weights_fin5 = mean(size_fin5(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_fin5(1+12*(k-1):12+12*(k-1),:),2)));
        
        weights_manu1 = mean(size_manu1(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_manu1(1+12*(k-1):12+12*(k-1),:),2)));
        weights_manu2 = mean(size_manu2(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_manu2(1+12*(k-1):12+12*(k-1),:),2)));
        weights_manu3 = mean(size_manu3(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_manu3(1+12*(k-1):12+12*(k-1),:),2)));
        weights_manu4 = mean(size_manu4(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_manu4(1+12*(k-1):12+12*(k-1),:),2)));
        weights_manu5 = mean(size_manu5(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_manu5(1+12*(k-1):12+12*(k-1),:),2)));

        weights_serv1 = mean(size_serv1(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_serv1(1+12*(k-1):12+12*(k-1),:),2)));
        weights_serv2 = mean(size_serv2(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_serv2(1+12*(k-1):12+12*(k-1),:),2)));
        weights_serv3 = mean(size_serv3(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_serv3(1+12*(k-1):12+12*(k-1),:),2)));
        weights_serv4 = mean(size_serv4(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_serv4(1+12*(k-1):12+12*(k-1),:),2)));
        weights_serv5 = mean(size_serv5(1+12*(k-1):12+12*(k-1),:),1)./mean((nansum(size_serv5(1+12*(k-1):12+12*(k-1),:),2)));

    end
    
 % Finance   
 
 returns_fin1(j)    = weights_fin1*(ret_fin1(j,:))';
 returns_fin2(j)    = weights_fin2*(ret_fin2(j,:))';
 returns_fin3(j)    = weights_fin3*(ret_fin3(j,:))';
 returns_fin4(j)    = weights_fin4*(ret_fin4(j,:))';
 returns_fin5(j)    = weights_fin5*(ret_fin5(j,:))';

 % Manufactoring
 
 returns_manu1(j)    = weights_manu1*(ret_manu1(j,:))';
 returns_manu2(j)    = weights_manu2*(ret_manu2(j,:))';
 returns_manu3(j)    = weights_manu3*(ret_manu3(j,:))';
 returns_manu4(j)    = weights_manu4*(ret_manu4(j,:))';
 returns_manu5(j)    = weights_manu5*(ret_manu5(j,:))';

 % Services 
 
 returns_serv1(j)    = weights_serv1*(ret_serv1(j,:))';
 returns_serv2(j)    = weights_serv2*(ret_serv2(j,:))';
 returns_serv3(j)    = weights_serv3*(ret_serv3(j,:))';
 returns_serv4(j)    = weights_serv4*(ret_serv4(j,:))';
 returns_serv5(j)    = weights_serv5*(ret_serv5(j,:))';

end

%=========================================================================
%                  AVERAGE RETURNS OVER THE QUARTER
%=========================================================================

[T,M]   = size(returns_monthly);
T       = T/3;
returns_quarterly = zeros(T,4);

for j=1:int32(T)
    
    returns_quarterly(j,1) = prod(returns_fin1(1+(j-1)*3:3+(j-1)*3,:),1);
    returns_quarterly(j,2) = prod(returns_fin2(1+(j-1)*3:3+(j-1)*3,:),1);
    returns_quarterly(j,3) = prod(returns_fin3(1+(j-1)*3:3+(j-1)*3,:),1);
    returns_quarterly(j,4) = prod(returns_fin4(1+(j-1)*3:3+(j-1)*3,:),1);
    returns_quarterly(j,5) = prod(returns_fin5(1+(j-1)*3:3+(j-1)*3,:),1);

    returns_quarterly(j,6) = prod(returns_manu1(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,7) = prod(returns_manu2(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,8) = prod(returns_manu3(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,9) = prod(returns_manu4(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,10) = prod(returns_manu5(1+(j-1)*3:3+(j-1)*3,:),1);    

    returns_quarterly(j,11) = prod(returns_serv1(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,12) = prod(returns_serv2(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,13) = prod(returns_serv3(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,14) = prod(returns_serv4(1+(j-1)*3:3+(j-1)*3,:),1);    
    returns_quarterly(j,15) = prod(returns_serv5(1+(j-1)*3:3+(j-1)*3,:),1);    
    
end

%=========================================================================
%                       SAVE PORTFOLIOS
%=========================================================================

save Matfiles/portfolio_industrysize returns_quarterly

