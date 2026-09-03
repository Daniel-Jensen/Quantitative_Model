function [x s] = smolyakapprox_step1AS(d,mu,a,b)
%d=4;
%mu=[0;2;0;1];
%a=[-1;-1;-1;-1];
%b=[1;1;1;1];
%smolyakapprox_step1 Smolyak Grid Construction 
%   [x s] = smolyakapprox_step1(d,mu,a,b) constructs the Smolyak
%   grid points for a hypercube in R^d defined by bounds a,b with level of 
%   approximation mu.  It outputs two elements, x and s.  x is
%   the collocation points that your function should be evaluated.
%   s is a structure containing information used by
%   smolyakapprox_step2 and smolyakapprox_step3.
 


    q = max(d,mu+1);        
    mu_m=max(mu);
    s.mu_m=mu_m;
    s.q = q;
    s.d = d;
    s.mu = mu;
    s.a = a(1:d);
    s.b = b(1:d);
    s.M_mup1 = smolyakM(mu+1); %This constant is used in step3
    
    if (length(a)>d), warning('length of a is longer than dim d'), end
    if (length(b)>d), warning('length of b is longer than dim d'), end
    
    %Construct necessary coefficients using fx values.
    % First, enumerate all necessary theta
    % -For all i satisfying q<=ibar<=d+mu, 
    % -need { l | l(1) is in 1...m(i1), l(2) is in 1...m(i1)}
    
    %Enumerate all such i, then all such m(i), then all such l
    %Determine all i that fit the criterion q<=|i|<=d+mu
    tmp=[];
 
        for ibar = min(d,min(q)):d+mu_m
            enum = smolyakEnumerate(d,ibar-d);        
            enum = enum+1;  
            
            for k=1:size(enum,1)
                tst=0;
                for l=1:d
                    if enum(k,l)<= mu(l)+1
                        tst=tst+1;
                    end 
                end
                if tst==d
                    tmp = [enum(k,:);tmp];
                end
            end
            
        end
        tmp = unique(tmp,'rows');
        leni = size(tmp,1);
        s.leni = leni;
        s.i = tmp;
        s.ibar = sum(s.i,2);
    
    %up there check and see if it has to start from biggest q!!!!
   

    
    
    
    s.k = smolyakM(s.i);
    s.k1 = smolyakM1(s.i);
    s.Lb = NaN(s.leni,1);
    s.Ub = NaN(s.leni,1);  
    Lc = NaN(s.leni,1);
    Uc = NaN(s.leni,1);
    tmpInt = 0;
    tmpInt1= 0;
    for iind = 1:leni        
        s.Lb(iind) = tmpInt + 1;
        Lc(iind)= tmpInt1+1;
        tmpInt = tmpInt + prod(s.k(iind,:));
        tmpInt1= tmpInt1 + prod(s.k1(iind,:));
        s.Ub(iind) = tmpInt;
        Uc(iind)= tmpInt1;
    end

   s.z = NaN(s.Ub(leni),d);
   s.j = NaN(s.Ub(leni),d);
    
    for iind = 1:leni
        ztmp = [];
        jtmp = [];
        for dind = 1:d
            ztmp = cartprod(ztmp,smolyakG(s.k(iind,dind)));
            jtmp = cartprod(jtmp,(1:s.k(iind,dind))');
        end
        s.z(s.Lb(iind):s.Ub(iind),:) = ztmp;
        s.j(s.Lb(iind):s.Ub(iind),:) = jtmp;
    end
    s.l = s.j;
    s.l=unique(s.l,'rows');
    for dind = 1:d
        s.x(:,dind) = (s.z(:,dind)+1.d0)*(b(dind)-a(dind))/2.d0 + a(dind);
    end    
     

    x=unique(s.x,'rows');
    
    x = x';
%end