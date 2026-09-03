function [x,a] = parsolve_post(fnhandle,x0,cc,tol,maxcount)

% parsolve was written by Ryan Decker, University of Maryland Economics.
% Version 2.0, November 2012
% fnhandle must be a valid function handle for the function of which you desire to find roots.
% x0 is the starting value guess.
% tol is an optional argument for the solve tolerance; default is 1e-10.
% maxcount is an optional argument for timing out the solver; default is
% 1000.

d = 1;

if nargin<5
    maxcount = 150;
end

if nargin<4
    tol = 1e-6;
end

count = 0;

gap = 10;

m = length(x0);

x = x0;

eps = 1e-6;

epsvec = [eye(m)*eps,zeros(m,1)];

df = zeros(m,m+1);

J = zeros(m,m);

while gap>tol && count<=maxcount
    
    count = count+1;
    
    parfor i = 1:(m+1)
		
     epstemp    = epsvec(:,i)';
		
        df(:,i) = feval(fnhandle,x+epstemp);
  
    end
    
    for i=1:m
        
        J(:,i) = (df(:,i)-df(:,m+1))/eps;
    
    end
             
     x = x - cc*(J\df(:,m+1))';
       
    gap = sum(df(:,m+1).^2);
    fprintf('Norm of residual function: %4.3g\n', gap)

    if gap>100
        d=0;
        break
    end

end
    
   
    b = isreal(x);
    c = min(isfinite(x));
    
    if b==0 || c==0 || d==0
        a = 0;
    else
        a = 1;
    end
    
if count<maxcount && a==1 
    disp('parsolve completed because the function value is below tolerance.')
    a=1;
else
    disp('parsolve stopped because iterations exceeded the maximum specified.')
    a=0;
end
    
    
    
end