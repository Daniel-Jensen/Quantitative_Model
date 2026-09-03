%==========================================================================
%                    TABLE A-3: CONFIDENCE INTERVAL FOR R2 
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

path('C:\Users\Luigi\Dropbox\projects\The Pass-Through of Sovereign Risk\Final version_JPE\Replication Files\Matfiles',path);
path('Files',path);

%=========================================================================
%                              LOAD DATA
%=========================================================================

load portfolios_beta
load portfolio_industrysize
load sdf

%=========================================================================
%                                OPTIONS
%=========================================================================

disp('                                                                  ');
disp('                                                                  ');
reply2  = input('Which specification? (1 for Bench, 2 for CAPM, 3 No Lev):');
disp('                                                                  ');
disp('                                                                  ');
reply4  = input('Sample R2 of Cross-sectional Regression? :  ');
disp('                                                                  ');
disp('                                                                  ');

reply  = 3;
reply3 = 2007.75;

if reply2==1
    ret_beta = ret_lev;
    
elseif reply2==2 
    ret_beta = ret_capm;
    
elseif reply2 ==3
    ret_beta = ret_nolev;
end


   if reply==3 

       port      = [ret_beta,returns_quarterly];
       
       
   elseif reply==1
       
       port      = returns_quarterly; 
   
   elseif reply==2
       
       port      = ret_beta;

   end

   [a,last_period] = min(abs(reply3-sdf(:,1)));

%=========================================================================
%          STEP 1: Compute weights reproducing a given R2
%=========================================================================

N_port    = size(port,2);

port      = port - sdf(1:end-4,3)*ones(1,N_port);

port      = port(1:last_period,:);

Time      = size(port,1);

Sigma     = cov(port);

N  = 10;    % Number of R2

K  = 20;    % Number of Repetitions for first step

r2 = linspace(0.01,0.40,N);

x  = 0;

factors = zeros(N_port,N);
weights = zeros(N_port,N);

for k=1:K
    
    for j=1:N
        
           while x==0
    
    C = chol(Sigma)'*randn(N_port,1);
    
    X = [ones(N_port,1),C];
    
    Y = mean(port,1)';
    
    B = inv(X'*X)*X'*Y;
    
    res = Y-X*B;
    
    R2  = (1-var(res)/var(Y));
        
    y   = abs(R2-r2(j));
    
    if y<0.0015
        
        x =1;
        
    end
 
    
           end       
         
   factors     = C;
w              = inv(Sigma)*factors;
w              = w./(ones(N_port,1)*sum(w,1));
weights(:,j,k) = w;
   
   x = 0;
    end
     x = 0;
     
end

% Find weights


%=========================================================================
%          STEP 2: Bootstrap sample R2 given population R2
%=========================================================================
   
Nsim        = 15000;
r2_sampling = zeros(N,Nsim,K);


for k=1:K
    for j=1:N
        parfor m=1:Nsim
                
[Y idx]= datasample(port,Time,1);

factor = (weights(:,j,k)'*Y')';

X      = [ones(Time,1),factor];

B      = inv(X'*X)*X'*Y;

betas  = B(2,:)';
    
X      = [ones(N_port,1),betas];

Y      = mean(Y,1)';

B      = inv(X'*X)*X'*Y;

res    = Y-X*B;

R2     = 1-(var(res)/var(Y));

r2_sampling(j,m,k) = R2;

        end
    end
end

%=========================================================================
%                   REPORT CONFIDENCE INTERVAL
%=========================================================================

r2_samp      = zeros(N,K*Nsim);

for j=1:N
    
r2_samp(j,:) = reshape(r2_sampling(j,:,:),K*Nsim,1)';

end

r2_80                = prctile(r2_samp,80,2);
[x y1]               = min(abs(reply4-r2_80));

disp('                                                                  ');
disp('CONFIDENCE INTERVAL');
disp('                                                                  ');
disp(['LOWER BOUND:            ', num2str(r2(y1))]);







    
    
    
    
    

