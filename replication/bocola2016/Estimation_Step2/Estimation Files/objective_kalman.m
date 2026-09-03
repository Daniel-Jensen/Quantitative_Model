function [obj,statepredi] = objective_kalman(param,data,constr)

if constr==1
    
    param(2) = exp(param(2))/(1+exp(param(2)));
    
    param(3) = exp(param(3));
    
end

A        = param(1);

x        = 0;

Phi      = param(2);

R        = param(3);

B        = 1;

Se       = 1;

H        = 0;

y        = data(:,1);

[liki,measurepredi,statepredi,varstatepredi] = kalman(A,B,H,R,Se,Phi,y-x);

obj      = sum(liki);

end

