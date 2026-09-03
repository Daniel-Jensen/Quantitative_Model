% Given an integer i, deliver the function m(i) = 2^(i-2)
function o = smolyakM1(i) 
%i=5;
    o = NaN(size(i));
    if sum(i<1)>0
        error('i must be >= 1')
    else
        o(i==1) = 1;
        o(i==2) = 2;
        o(i>2) =  2.^(i(i>2)-2);
    end   
end