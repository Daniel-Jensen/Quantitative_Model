function [T] = T(y,s)
 %y=x';
%  y=[1;0;-1;1];
 d = s.d;
 y = y';

 T = ones(size(y,1),max(s.M_mup1),d);

 for dind = 1:d
        T(:,2,dind) = y(:,dind);
 end 
    
    % q<=ibar<=d+mu.  So ibar=d+mu is largest.  
    %Largest any order Cheby can be is the largest any component of m(i) can be.  
    %This occurs when i = 1 except for one component. 
    % => icomp = d+mu-(d-1) = mu+1 => i = mu+1 => m(i) = m(mu+1);
    % NOTE: I precompute this now as s.M_mup1
    
    twoz = 2*y;
    for oind = 3:max(s.M_mup1)
        for dind = 1:d
            T(:,oind,dind) = twoz(:,dind).*T(:,oind-1,dind) - T(:,oind-2,dind);
        end
    end

    %For each possible value of l, compute the product across i of
    %T(:,li) where li is a component of l
    
    a=1;
    for j=1:d
        a=a.*s.k1(:,j);
    end
    Tprod = ones(size(y,1),sum(a));
   
    for dind = 1:d
        for lind = 1:sum(a)
     
            if (s.l(lind,dind)>1) 
                Tprod(:,lind) = Tprod(:,lind).*T(:,s.l(lind,dind),dind);
            end
        end
    end


   T = Tprod';

   
end

