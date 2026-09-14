function [lnprior] = prior(Theta)

% This function evaluates the prior of model's parameters
%
% Inputs: Theta (vector collecting the structural parameters to be estimated)
%
% Output: lnprior (log-prior evaluated at Theta)
%
% Author: Luigi Bocola     luigi.bocola@northwestern.edu
% Date  : 09/06/2015

     if Theta(1)<0 || Theta(2)>=0.999 || Theta(4)<=0 || Theta(4)>0.002 || Theta(6) < 0 

         lnprior = -10000000;
         
     else
         
         para1 = [0.5,.35];
         para2 = [0.25,0.25];
             a = (1-para1).*para1.^2./para2.^2 - para1;
             b = a.*(1./para1 - 1);
            P1 = normpdf(Theta(1)*100,0.96,25); 
            P2 = normpdf(Theta(2),0.97,0.15);
            P3 = betapdf(Theta(3),a(1),b(1));
            P4 = normpdf(Theta(4)*400,1.00,1);
            P5 = betapdf(Theta(5),a(2),b(2));
            P6 = exp(lnpdfig(Theta(6)*100,.75,2));
            
            lnprior = log(P1)+log(P2)+log(P3)+log(P4)+log(P5)+log(P6);

     end

end

function y = lnpdfig(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);
end




