%==========================================================================
%                       Estimation of nonlinear model
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

path('Smolyak Files',path);
path('Solution Files',path);
path('Simulation Files',path);
path('Estimation Files',path);
path('Matfiles',path);

%=========================================================================
%                               Load data
%=========================================================================

load data            % Load Lagrange multiplier and GDP growth

data                 = [data(:,3),data(:,2)];
Time                 = size(data,1);

%=========================================================================
%             Load model solution and candidate density
%=========================================================================

load model_nodefault

load candidate

N_g = size(gamma,2);

rhoz     = param(6);
sigmaz   = param(7);
rhog     = param(16);
sigmag   = param(17);

[EXP,exp_g] = precomp_integral(coll_points,bounds,rhog,sigmag,rhoz,sigmaz,s);

%=========================================================================
%               Define objects for Metropolis Hastings
%=========================================================================

% Define Objects for MH Algorithm

nparticle       = 100000;               % Number of Particles 
Nsim            = 10000;                % Number of Posterior Draws
Thetasim        = zeros(6,Nsim);        % Posterior Draws
logposterior    = zeros(Nsim,1);        % Logposterior
obs             = zeros(2,Time,Nsim);   % Filtered Observations
state           = zeros(5,Time,Nsim);   % Filtered States
model_sol       = zeros(Nsim,N_g);      % Model Solution
Effective_par   = zeros(Time,Nsim);     % Effective Particles
keep            = zeros(Nsim,1);
m               = 2;                    % Current draw
counter         = 4;                    % Counter for intermediate output
accept          = 0;                    % Indicator if draw is accepted
c               = 0.005;

%Initialize the MH Algorithm at posterior mode

Thetasim(:,1)    = Theta;
model_sol(1,:)   = gamma;
model_exp(1,:,:) = EXP;
[liki]           = dsgeliki(Theta,gamma,bounds,coll_points,TT,s,exp_g,V,data,nparticle);
obj              = liki+prior(Theta);
logposterior(1)  = obj;

%=========================================================================
%                   Metropolis-Hastings algorith, 
%=========================================================================

disp('                                                                  ');
disp('                                                                  ');    
disp('         Starting Metropolis Hastings Algorithm...                ');
disp('                                                                  ');    
disp('                                                                  ');    

while m<=Nsim
   
    tic
    
    Thetacand    =  mvnrnd(Thetasim(:,m-1),c*Sigma);  % Sample a candidate draw    
    lnprior      =  prior(Thetacand);                 % Prior
   
    if lnprior<=-10000000 
    Thetasim(:,m)       = Thetasim(:,m-1);
    logposterior(m)     = logposterior(m-1);    
    model_sol(m,:)      = model_sol(m-1,:);
    obs(:,:,m)          = obs(:,:,m-1);
    state(:,:,m)        = state(:,:,m-1);
    Effective_par(:,m)  = Effective_par(:,m-1);
    else

[liki,gamma,N_eff,EXP,a]  = dsgeliki(Thetacand,model_sol(m-1,:),bounds,coll_points,TT,s,...
                            exp_g,V,data,nparticle);
                      
keep(m) = a;                      
        
objcand       = liki+lnprior;

alpha         = min(1,exp(objcand-obj)); % MH acceptance probability

u = rand(1);

if u<=alpha                        % Accept the Draw
    Thetasim(:,m)       = Thetacand';
    accept              = accept+1;
    obj                 = objcand;   
    logposterior(m)     = obj;    
    model_sol(m,:)      = gamma;
    Effective_par(:,m)  = N_eff;
    
else
    
    Thetasim(:,m)       = Thetasim(:,m-1);
    logposterior(m)     = logposterior(m-1);    
    model_sol(m,:)      = model_sol(m-1,:);
    Effective_par(:,m)  = Effective_par(:,m-1);

end

    end
    
acceptancerate     = accept/m;

counter            = counter+1;

time = toc;

if counter==5
disp('                 Information on the Chain                         ');    
disp('                                                                  ');    
fprintf('  Draw Number                     = %5.2f\n', m);
fprintf('  Acceptance Rate                 = %5.2f\n', acceptancerate);
fprintf('  Log-posterior                   = %5.2f\n', logposterior(m));
fprintf('  Number of Effective Particles   = %5.2f\n', mean(Effective_par(:,m)));
fprintf('  Time for posterior draw         = %5.2f\n', time);
disp('                                                                  ');    
fprintf('para    current        mean    std\n')
fprintf('----   ------------   ----    ----\n')
fprintf('mu_ss   %4.3f         %4.3f    %4.3f\n', Thetasim(1,m)*100, mean(Thetasim(1,1:m)*100), var(Thetasim(1,1:m)*100)^(1/2))
fprintf('psi     %4.3f         %4.3f    %4.3f\n', Thetasim(2,m), mean(Thetasim(2,1:m)), var(Thetasim(2,1:m))^(1/2))
fprintf('csi     %4.3f         %4.3f    %4.3f\n', Thetasim(3,m), mean(Thetasim(3,1:m)), var(Thetasim(3,1:m))^(1/2))
fprintf('gammaz  %4.3f         %4.3f    %4.3f\n', Thetasim(4,m)*400, mean(Thetasim(4,1:m))*400, var(Thetasim(4,1:m)*400)^(1/2))
fprintf('rhoz    %4.3f         %4.3f    %4.3f\n', Thetasim(5,m), mean(Thetasim(5,1:m)), var(Thetasim(5,1:m))^(1/2))
fprintf('sigmaz  %4.3f         %4.3f    %4.3f\n', Thetasim(6,m)*100, mean(Thetasim(6,1:m))*100, var(Thetasim(6,1:m)*100)^(1/2))
disp('                                                                  ');    
disp('                                                                  ');

save Matfiles/param_step1_prov Thetasim obj m logposterior model_sol accept Effective_par  

counter = 0;

end

m = m+1;

end    
        
%=========================================================================
%                         Save Results
%=========================================================================

save Matfiles/param_step1 Thetasim logposterior model_sol  


    
    
    
    