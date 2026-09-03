%Returns cartesian prod of set and set B
function C = cartprod(A,B)
    if (isempty(B))
        C = A;
        return
    elseif (isempty(A))
        C = B;
        return
    end
%     C = NaN([size(A,1)*size(B,2) size(A,2)+size(B,2)]);
%     i = 0;
%     for a = 1:size(A,1)
%         for b = 1:size(B,1)
%             i = i + 1;
%             C(i,:) = [A(a,:), B(b,:)];            
%         end
%     end
    i = 1;
    bin = size(B,1);
    C = NaN([size(A,1)*size(B,2) size(A,2)+size(B,2)]);
    for a = 1:size(A,1)
        tmp = [repmat(A(a,:),[size(B,1) 1]), B];        
        C(i:i+bin-1,:) = tmp;
        i = i+bin;
    end
end 
