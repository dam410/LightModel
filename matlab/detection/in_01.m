function [result] = in_01(T)
        check = T>=0 & T<=1;
        result = prod(check);
end

