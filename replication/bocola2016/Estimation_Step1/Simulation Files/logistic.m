function [r] = logistic(n,mu,sigma)

p=rand(n,1);

r=log(p./(1-p)).*sigma+mu;

end

