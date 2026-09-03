%==========================================================================
%                        ESTIMATION: STEP 2
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

path('Estimation Files',path);

%=========================================================================
%                       LOAD AND PREPARE DATA
%=========================================================================

load Matfiles/def_prob

def_probabilities = physical(6:end,2);

T   = size(def_probabilities,1);

%=========================================================================
%                    GENERATE CANDIDATE DENSITY
%=========================================================================

options = optimset('Display','iter','FunValCheck','on','MaxFunEvals',300000,'MaxIter',1000);

data    = log((def_probabilities)./(1-(def_probabilities)));

obj     = @(Theta) -(objective_kalman(Theta,data,1)+prior(Theta,1));

param   = [mean(data),2,log(0.1)];

param   = fminunc(obj,param,options);

param(2) = exp(param(2))/(1+exp(param(2)));

param(3) = exp(param(3));

mode    = param';

obj     = @(Theta) -(objective_kalman(Theta,data,0)+prior(Theta,0));

Sigma   = nhess(obj,mode);

Sigma   = inv(Sigma);

%=========================================================================
%                   2) METROPOLIS HASTINGS
%=========================================================================

obj   = objective_kalman(mode,data,0)+prior(mode,0);

Nsim  = 100000;

c     = 0.05;

Thetasim = mode*ones(1,Nsim);

liki     = obj*ones(Nsim,1);

accept   = 0;

counter  = 0;

state    = zeros(T+1,Nsim);

for i=1:Nsim
    
    cand     = mvnrnd(Thetasim(:,i),c*Sigma);

    priorcand = prior(cand,0);

    if cand(2)>0.99
        obj_cand = -1000000;
    else
    [obj_cand,statepredi] = objective_kalman(cand,data,0);
    end
    
    alpha    = min(1,exp(obj_cand+priorcand-obj));
 
           u = rand(1);

if u<=alpha

    state(:,i+1)      = statepredi;
    Thetasim(:,i+1)   = cand;
    accept            = accept+1;
    obj               = obj_cand+priorcand;
    liki(i+1)         = obj_cand+priorcand;

else
  
Thetasim(:,i+1)   = Thetasim(:,i);

    liki(i+1)     = obj;   

state(:,i+1)      = state(:,i);

end

acceptancerate     = accept/i;

counter            = counter + 1;

if counter==5000
disp('                                                                  ');
disp(['                               DRAW NUMBER:', num2str(i)]         );
disp('                                                                  ');
disp('                                                                  ');    
disp(['                           ACCEPTANCE RATE:', num2str(acceptancerate)]);
disp('                                                                  ');
disp('                                                                  ');    
disp('                            RECURSIVE AVERAGES                    ');
disp('                                                                  ');
disp('S_STAR       RHO_S          SIGMA_S');
disp(num2str(mean(Thetasim(:,1:i-1)')));
disp('                                                                  ');
disp('                                                                  ');    
disp(['                           RECURSIVE LIKI:', num2str(mean(liki(1:i)))]);
disp('                                                                  ');
disp('                                                                  ');    

counter = 0;

end

end

%=========================================================================
%                              SAVE RESULTS
%=========================================================================

save Matfiles/param_step2 Thetasim Nsim

param_step2 = zeros(3,10);
steps       = floor(size(Thetasim,2)/10);

for j=1:10
    param_step2(:,j) = Thetasim(:,j*steps);
end

save Matfiles/param_defprocess param_step2

