% Evaluates the cheby poly of order "order" at the points z.
% Each row of z is a new point, each elt.  
function [T] = ChebEvalOneDim(order,z)
    if (order==0)
        T = 1; %ones(size(z,1),1); %Should just be able to return 1;
        return
    elseif (order==1)
        T = z;
        return
    else
        Tm2 = ones(size(z,1),1);
        Tm1 = z;
        for ind = 2:order-1
            T = 2*z.*Tm1 - Tm2;
            Tm2 = Tm1;
            Tm1 = T;
        end
        T = 2*z.*Tm1 - Tm2;
    end
end