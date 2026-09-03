function [val]=choose(n,k)
    %n!/(k!*(n-k)!) = n*(n-1)*...*max(k,n-k)/(min(k,n-k))!
    if (k>n | k<0) 
        error('k must be in [0,n]')
    end
       
    if (k>=n-k)
        val = prod(k+1:n)/factorial(n-k);
    else
        val = prod(n-k+1:n)/factorial(k);        
    end
%     
%     val2 = factorial(n)/factorial(k)/factorial(n-k);
%     if (val-val2~=0)
%         error(' ')
%     end
end