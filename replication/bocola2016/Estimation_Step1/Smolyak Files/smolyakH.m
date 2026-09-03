% Given the dimensionality of the problem (d) and an order of approximation
% mu>=1, construct the Smolyak grid
function [grid] = smolyakH(d,mu,a,b)
    
    ibar = d+mu; 
    
    %Find all combinations of i in Z^d_{++} s.t. |i| = ibar where |i| = i1
    %+i2 + i3 + ... + id
            
    %The complete enumeration will be given my a matrix of dimension [ x d]
    %Complete enumeration can be done in the following way.    
    enum = smolyakEnumerate(d,mu);    
    enum = unique(enum,'rows') + 1;           

    q = max(d,mu+1);
    tmp=[];
    for ibar = q:d+mu
        enum = smolyakEnumerate(d,ibar-d);        
        enum = enum+1;        
        tmp = [enum;tmp];
    end     
    tmp = unique(tmp,'rows');
    enum = tmp;
    
%     disp('enum')
%     enum
%     disp('size enum')
%     disp(size(enum))
%     
    %Check enumeration
%     if (any(ibar~=sum(enum,2)) )
%         error('enum wrong')
%     end 
        
    %Compute all the grid points associated with the enumeration
    %Note: a row of the enumeration gives a
    %smolyakG(m(row(1)))xsmolyakG(m(row(2)))xsmolyakG(m(row(3)))x...xsmolyakG(m(row(d)))    
    smolyakind = 1;    
    for enumind = 1:size(enum,1)
        C = [];
        Cmap = [];
        enumrow = enum(enumind,:);            
        % Enumrow gives an order for sets
        for j = 1:d                       
            % For a given enumrow, we need to choose an element from 
            Gj = (1:smolyakM(enumrow(j)))';
            Cmap = cartprod(Cmap,Gj);
            
            Gj = smolyakG(smolyakM(enumrow(j)));
            C = cartprod(C,Gj);
           
        end
        if any(size(C)~=size(unique(C,'rows')))
            error('C not unique')
        end
        if any(size(Cmap)~=size(unique(Cmap,'rows')))
            error('Cmap not unique')
        end
        
        grid.val(smolyakind:smolyakind+size(C,1)-1,:) = C;        
        grid.order(smolyakind:smolyakind+size(C,1)-1,:) = repmat(smolyakM(enumrow),[size(C,1) 1]);
        grid.index(smolyakind:smolyakind+size(C,1)-1,:) = Cmap;
        
        smolyakind = smolyakind+size(C,1);
    end
       
end