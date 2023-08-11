function [sourceMult] = sourceScalingFDTD(order)
    switch order
        case 2
            sourceMult = 0.00797032030276674519;
        case 6
            sourceMult = 0.00797032030276674519;
        otherwise
            error('scheme order not supported')
    end    
end