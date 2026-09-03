% - Cross_section
% 09/05/2015

There are 2 files in the folder. 

1) portfolio_beta.m               : It aggregates stocks of the Italian exchange into 10 portfolios sorted by the decile of the estimated distribution of
                                    betas. The betas are computed with respect to the model implied stochastic discount factor (benchmark), the 
                                    market return (CAPM), and the model stochastic discount factor with psi=0 (No leverage). The data are obtained from 
                                    "compustat_99-11.xls" in the subfolder "Data". The file output is "portfolio_beta.mat" in the subfolder "Matfiles".

2) portfolio_industrysize.m       : It aggregates stocks of the Italian exchange into 15 portfolios sorted by industry and size. Stocks are first sorted 
                                    by industry ("Manufacturing", "Services" and "Financial"). Within each of these categories, stocks are sorted into 
                                    quintiles of the size distribution. Return for each portfolio is value-weighted. The data are obtained from 
                                    "compustat_99-11.xls" in the subfolder "Data". The file output is "portfolio_industrysize.mat" in the subfolder 
                                    "Matfiles".

