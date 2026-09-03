function [r] = logistic(n,T,mu,sigma)

p=rand(n,T);

r=log(p./(1-p)).*sigma+mu;

end

