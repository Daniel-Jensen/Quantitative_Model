% Given an integer n, construct the set of n extrema of the Chebyshev 
% polynomials zeta_j = -cos(pi*(j-1)/(n-1)), j = 1,n
function [set j]= smolyakG(n)
    if (n<1) 
        error('smolyakG: n must be >=1')
    elseif (n==1)
        set = 0;
    else       
        j = (1:n)';
        set = -cos(pi*(j-1)/(n-1));
        set(abs(set)<1d-12) = 0;
        set(abs(set-1.d0)<1d-12) = 1;
        set(abs(set+1.d0)<1d-12) = -1;
    end
end