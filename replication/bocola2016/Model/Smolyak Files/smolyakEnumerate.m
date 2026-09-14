%Given a  d (dimensionality of state), enumerates all possible ways to 
%add one mu times to d (where adding one can occur in any of d positions).  
function [enum] = smolyakEnumerate(d,mu)
    if (mu==0)
        enum = zeros([1 d]);
        return
    end
    if (mu>=1) 
        enum_mum1 = smolyakEnumerate(d,mu-1);
        m = size(enum_mum1,1);
        %Given all previous enumerations, I can add one in d places.
        for i = 1:d
            enum(1+(i-1)*m:i*m,:) = enum_mum1;            
            enum(1+(i-1)*m:i*m,i) = enum(1+(i-1)*m:i*m,i) + 1;
        end
        return
    end 
end