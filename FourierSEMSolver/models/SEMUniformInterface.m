classdef SEMUniformInterface        
    properties
        dx, Nk, loc, P_order, uniform_distr
    end
    
    methods
        function this = SEMUniformInterface(dx,Nk,P_order,uniform_distr,loc)
            this.uniform_distr = uniform_distr;
            this.P_order = P_order;
            this.dx = dx;
            this.Nk = Nk;
            this.loc = loc;
        end        
    end
end

