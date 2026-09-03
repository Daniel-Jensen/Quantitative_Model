% - Model 
% 09/05/2015

There are 12 main files in the folder. Most of the files are designed to exploits the Parallel toolbox in Matlab. For Matlab versions earlier than 2014, 
the parallel toolbox needs to be activated by adding the line "matlabpool" at the beginning of the file. Approximate running times are reported for a 8 
dual cores Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz and 32 GB of RAM.
 
1) model_solution_mean.m          : It generates the solution to the model when parameters are at their posterior median. The model solution is stored
                                    in "model_solution_mean.mat". It uses the file "model_nodefault_mean.mat". Approximate running time: 30 mins.
                                    
2) model_solution_posterior.m     : It generates the solution to the model for multiple draws of model's parameters. The model solution is stored
                                    in "model_solution_posterior.mat". It uses the file "model_nodefault_posterior.mat" and "param_defprocess.mat". 
                                    These two files are generated in folder "Estimation_step1" and "Estimation_step2". Approximate running time: 30 mins 
                                    per draw.

3) decomposition_riskpremia.m     : It generates the file "decomposition.mat" used in the construction of Table 4 and Figure 4. It uses the matfile
                                    "model_solution_median.mat". Approximate running time: 2 mins.

4) generate_irf.m                 : It generates the file "irf_bench.mat" used in the construction of Figure 3 and Figure 5. It uses the matfile
                                    "model_solution_median.mat". Approximate running time: 35 mins.

5) counterfactual_outputlosses.m  : It generates the file "output_losses.mat" used in the construction of Table 5. The file uses 
                                    "model_solution_posterior.mat". Approximate running time: 6 mins per draw.

6) counterfactual_returns.m       : It generates the file "counterfactual.mat" used in the construction of Figure 6. The file uses 
                                    "model_solution_posterior.mat". Approximate running time: 50 mins per draw.

7) estimate_states.m              : It generates the file "state_italy.mat" used in the construction of Figure 8 and Table A-4. It uses the file
                                    "model_solution_median.mat". Approximate running time: 3 mins.

8) generate_price_risk.m          : It generates the matfiles "price_risk.mat" used in the construction of Table A-4. It uses the matfile
                                    "model_solution_median.mat" and "state_italy.mat". Approximate running time: 25 mins.

9) ltro_policies.m                : It generates the file "policies_ltro.mat". This file is used by "ltro_experiment.mat" to generate the LTROs experiment
                                    underlying Figure 8. The files uses "model_solution_median.mat". Approximate running time: 160 minutes.

10) ltro_experiment.m             : It generates the file "ltro.mat" used in the construction of Figure 8. It uses "model_solution_median.mat" and 
                                    "policies_ltro.mat". Approximate running time: 19 hours.

11) generate_irf_open.m           : It generates the file "irf_open.mat" used in the construction of Figure 7. Approximate running time: 3 hours.

12) generate_euler_errors.m       : It generates the matfiles "euler_errors.mat" used in the construction of Figure A-1. It uses the file
                                    "model_solution_median.mat". Approximate running time: 3 mins.

These files use a number of procedures collected in the sub-folders "Solution Files" and in "Simulation Files". The most relevant files are

a) residual_model.m               : See the description of "residual_nodefault.m" in Estimation_Step1.txt in the Folder "Estimation_Step1".

b) simul.m                        : See the description of "simul.m" in Estimation_Step1.txt in the Folder "Estimation_Step1".

c) model_policies.m               : See the description of "model_policies.m" in Estimation_Step1.txt in the Folder "Estimation_Step1".

d) particle.m                     : See the description of "particle.m" in Estimation_Step1.txt in the Folder "Estimation_Step1".

The folder "Smolyak Files" collects routines used for the construction of the Smolyak Grid and of the Chebyshev's polynomials. These routines are minor modifications of 
codes written by Grey Gordon, available at https://sites.google.com/site/greygordon/code.







