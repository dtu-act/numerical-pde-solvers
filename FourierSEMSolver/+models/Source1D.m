classdef Source1D        
    properties
        F, x0
    end
    
    methods
        function obj = Source1D(F,x0)
            obj.F = F;
            obj.x0 = x0;
        end
    end
end

