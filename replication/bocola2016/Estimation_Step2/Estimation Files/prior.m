function [lnprior] = prior(Theta,constr)

if constr==1
    
    Theta(2) = exp(Theta(2))/(1+exp(Theta(2)));
    
    Theta(3) = exp(Theta(3));
    
end

     if Theta(2)>=0.999

         lnprior = -10000000;
         
     else
         
         para1 = 0.5;
         para2 = 0.3;
             a = (1-para1).*para1.^2./para2.^2 - para1;
             b = a.*(1./para1 - 1);
            P1 = betapdf(Theta(2),a(1),b(1));
            P2 = normpdf(Theta(1),-7,5);
            P3 = exp(lnpdfig(Theta(3),.75,4));

            
            lnprior = log(P1)+log(P2)+log(P3);

     end

end

function y = lnpdfig(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);
end




