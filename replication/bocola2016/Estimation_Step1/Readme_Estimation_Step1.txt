% - Estimation_Step1
% 09/30/2015

There are 3 files in the main folder. The files are designed to exploits the Parallel toolbox in Matlab. For Matlab versions earlier than 2014, the parallel 
toolbox needs to be activated by adding the line "matlabpool" at the beginning of each file. Approximate running times are reported for a 8 dual cores 
Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz and 32 GB of RAM.

1) estimation_step1.m             : It estimates the model without sovereign risk using a Random Walk Metropolis Hastings. It loads the data from 
                                    "data.mat" in the subfolder "Matfiles". The original series can be obtained from "Data.xls". The draws for 
                                    the model parameters are saved in "param_step1.mat" in the subfolder "Matfiles". Approximate running time: 1 min.
                                    per posterior draw. 

2) predictive_checks.m            : It computes the posterior predictive checks reported in Figure 2. It loads posterior draws "param_step1.mat". The 
                                    predictive checks are saved in "predictive_checks.mat", in the subfolder "Matfiles". 

3) arrange_draws.m                : It loads the draws for model parameters, select a subset, a save the results in "model_nodefault_posterior".
                                  
These files use a number of procedures collected in the sub-folders "Solution Files", "Simulation Files", and "Estimation Files". The most relevant files are

a) residual_nodefault.m         :  It is the main file used in the numerical solution of the model. It takes as inputs the model's structural parameters 
                                   ("param"), an initial guess for the model solution ("gamma"), the bounds for the model state variables ("bounds"), the collocation 
                                   points ("coll_points"), the Chebyshev's polynomials evaluated at the collocation points ("TT"), a structure that is 
                                   necessary to evaluate polynomials ("s"), a vector collecting the deterministic steady-state values for model's variables 
                                   ("ss"), the matrix for integral precomputation ("EXP"), and a matrix that rotates the state variables ("V"). It returns the 
                                   vector of residuals evaluated at the collocation points. See the Online Appendix A for additional information.   

b) simul.m                      : It generates a simulation from the model. It takes as inputs the vector collecting structural parameters ("param), the numerical
                                   solution of the model ("gamma"), the bounds for the model state variables ("bounds"), the Chebyshev's polynomials evaluated at the collocation points
                                   ("TT"), a vector collecting the deterministic steady state of the model ("ss"), a structure used to evaluate Chebyshev's polynomials ("s"), a matrix 
                                   collecting structural shocks ("e"), and a vector collecting initial conditions for the state variables ("initial"). It returns simulations for state 
                                   variables ("state"), a structure collecting the simulation for several endogenous variables ("obs"), and the untransformed state variables. 

c) model_nodefault_policies.m   : Takes as inputs the vector of structural parameters ("param"), elements of the model solution ("gamma","bounds","TT",
                                  "coll_points","ss","s"), the matrix for the pre-computation of integrals ("EXP"), and a matrix that rotates the state variables ("V"). 
                                  It returns the coefficients parametrizing the model's policy functions. These coefficients are collected in the structure "policies".  

d) dsgeliki.m                   : It is the main file used in the estimation of the model. It takes as inputs a vector collecting the structural parameters to be estimated (Theta), 
                                  a guess for the numerical solution of the model at Theta ("gamma"), elements of the model solution ("bounds","coll_points","TT","s","exp_g","V"), a
                                  matrix collecting the data used in estimation ("data"), and the number of particles ("Nparticle"). It returns the evaluation of the log-likelihood at
                                  Theta ("liki"), the numerical solution of the model at Theta ("gamma"), the number of effective particles ("N_eff"), and the matrix collecting coefficients 
                                  for integral precomputation ("EXP"). 

The folder "Smolyak Files" collects routines used for the construction of the Smolyak Grid and of the Chebyshev's polynomials. These routines are modifications of codes written by Grey Gordon, 
available at https://sites.google.com/site/greygordon/code.


