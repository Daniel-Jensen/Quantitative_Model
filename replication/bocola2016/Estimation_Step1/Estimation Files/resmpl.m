function y=resmpl(x,w)
% Y=RESMPL(X,W)
% resamples from X with weights W, with replacement
% based on Kitagawa's (1996) deterministic resampling algorithm
%
% INPUT
%	X : presample, (nx by nsmpl)
%   W : weights, (1 by nsmpl)
%
% OUTPUT
%	Y : resample, (nx by nsmpl)

% Sungbae An
% Created: 08/06/2004
% Updated: 10/06/2004

% dimenstion
[nx,nsmpl] = size(x);

% initialize
ind = zeros(1,nsmpl);

% construct CDF
c = cumsum(w);

% draw a starting point
u = (rand-1+(1:nsmpl))/nsmpl;

% start at the bottom of the CDF
j=1;

for i=1:nsmpl
	% move along the CDF
	while (u(i)>c(j))
		j=j+1;
    end
	% assign index
	ind(i) = j;
end

y=x(:,ind);
