% Given an integer i, deliver the function m(i) = 2^(i-1) + 1
function o = smolyakM(i)        
    o = NaN(size(i));
    if sum(i<1)>0
        error('i must be >= 1')
    else
        o(i==1) = 1;
        o(i>1) =  2.^(i(i>1)-1) + 1;
    end   
end