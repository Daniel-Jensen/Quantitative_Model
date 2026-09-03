function [s] = smolyakapprox_step2(fx,s)
%smolyakapprox_step2 Smolyak Polynomial Construction 
%   [s] = smolyakapprox_step2(fx,s)
%   Given a list of values and s returned by step1, computes a
%   polynomial approximation.  If fx is a matrix, then it creates a
%   polynomial approximation treating each column as its own set of
%   function values.
% 
%   Grey Gordon 2011.  
%   "But God chose the foolish things of the world to shame the wise; 
%    God chose the weak things of the world to shame the strong." -1Cor12:47

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

    leni = s.leni;

    % Check dims
    if (size(fx,1)~=length(s.redo)) 
        error('fx dims are incorrect')
    end

    % Unpack fx 
    s.f = fx(s.undo,:);

    % Construct the polynomial coefficients
    s.theta = NaN(size(s.f));
    for iind = 1:leni
        for lind = s.Lb(iind):s.Ub(iind)
            tmp = (s.T(s.Lb(iind):s.Ub(iind),lind-s.Lb(iind)+1)./s.cjprod(s.Lb(iind):s.Ub(iind)))';
            s.theta(lind,:) = s.const(iind)/s.clprod(lind)*tmp*s.f(s.Lb(iind):s.Ub(iind),:);
        end 
    end
    
    % Precompute the constant times the coefficient
    s.constTimesTheta = s.theta;
    for iind = 1:leni
        s.constTimesTheta(s.Lb(iind):s.Ub(iind),:) = s.constVec(s.ibar(iind)-s.q+1)*s.constTimesTheta(s.Lb(iind):s.Ub(iind),:);
    end
    
end



