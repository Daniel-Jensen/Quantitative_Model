% - Estimation_Step2
% 09/06/2015

There are 2 files in the folder. 

1) estimation_step2.m             : It estimates equation (12) in the main text using a Random Walk Metropolis Hastings. It loads the data from 
                                    "def_prob.mat" in the subfolder "Matfiles". The original series can be obtained from "Data.xls". The draws for 
                                    the model parameters are saved in "param_step2.mat" in the subfolder "Matfiles". 

2) predictive_checks.m            : It computes the posterior predictive checks reported in Table 3. It loads posterior draws 
                                    "param_step2.mat". The predictive checks are saved in "predictive_checks_def.mat" in the subfolder "Matfiles".

These files use a number of routines collected in the subfolder "Estimation Files". The most relevant files are

a) objective_kalman.m             : It takes as inputs the vector of parameters ("param"), the times series of default probabilities ("data"), and an indicator
                                    equal to 1 if we are imposing constraints on the parameter vector ("constr"). It returns the log-likelihood function ("obj") and
                                    the filtered state ("statepredi").
                               
b) prior.m                        : It takes as inputs the vector of structural parameters ("param"), and an indicator equal to 1 if we are imposing constraints 
                                    on the parameter vector ("constr"). It returns the log-prior ("lnprior"). 
