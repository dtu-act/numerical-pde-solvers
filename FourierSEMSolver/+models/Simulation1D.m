classdef Simulation1D < matlab.mixin.Copyable
    properties (SetAccess = private)        
        p_current, % pressure
        p_prev,    % pressure
        p_history, % pressures for time interpolation        
        n % time index
        
        src, solver_type, custom
    end
    
    properties (SetAccess = public)                
        domain
    end
    
    methods(Access = protected)
        %https://se.mathworks.com/help/matlab/ref/matlab.mixin.copyable-class.html
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a shallow copy of all properties
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % Make a deep copy of the Domain object
            %cpObj.domain = copy(obj.domain);
            %cpObj.custom = copy(obj.custom);
        end
    end
   
    methods
        function self = Simulation1D(solver_type, domain, src, custom, n, p_current, p_prev)
            import utilsDD.*
            import models.*

            if nargin ~= 4 && nargin < 7
                error('not enough args')
            end
            
            if nargin == 4
                N = length(domain.x1d);
                self.p_current = zeros(1,N);
                self.p_prev = zeros(1,N);
                self.n = 0;
            else
                N = length(p_current);
                self.p_current = p_current;
                self.p_prev = p_prev;
                self.n = n;
            end
            
            history_length = 5;            
            self.p_history = zeros(N, history_length);
            
            self.domain = domain;
            
            if isempty(src)
                self.src = Source1D(zeroSource(domain.x1d),0);
            else
                self.src = src;
            end
            
            self.solver_type = solver_type;
            self.custom = custom;
        end
        
        % NOTE: history is not interpolated/redistributed
        function [self,cs1] = remesh(self, xx_out, endpoint_dzero)
            x = self.domain.x1d;
            
            if endpoint_dzero
                cs1 = spline(x,[0,self.p_current,0]);
                pc = ppval(cs1,xx_out);
                self.p_current = pc;
                
                cs2 = spline(x,[0,self.p_prev,0]);
                pp = ppval(cs2,xx_out);
                self.p_prev = pp;
            else
                cs1 = spline(x,self.p_current);
                pc = ppval(cs1,xx_out);
                self.p_current = pc;
                
                cs2 = spline(x,self.p_prev);
                pp = ppval(cs2,xx_out);
                self.p_prev = pp;
            end      
                        
            self.domain.x1d = xx_out;
        end
        
        function [self] = setHistory(self, history)
            N = length(self.p_current);
            
            if (size(history,1) ~= N)
                error("Spatial dimension doesn't agree")
            end
            
            self.p_history = history;
            
            self.p_current(:) = history(:,1);
            self.p_prev(:) = history(:,2);
        end
        
        function [self] = update(self, p_current, p_prev)
            
            self.p_current = p_current;
            self.p_prev = p_prev;
            self.n = self.n+1;
            
            self.p_history(:,2:end) = self.p_history(:,1:end-1);
            self.p_history(:,1) = self.p_current'; % TODO check: in case more pressure values are present 
                                                   %             we only keep the first in history                        
        end
        
        function [self] = adjust(self, p_current, p_prev, adjust_history)                        
            self.p_current = p_current;
            self.p_prev = p_prev;
            
            if nargin == 3 || adjust_history
                self.p_history(:,1) = self.p_current';
                self.p_history(:,2) = self.p_prev';
            end
        end
        
        function [self] = adjust_current(self, p_current)                        
            self.p_current = p_current;            
            self.p_history(:,1) = self.p_current';
        end
        
        function [self] = adjust_prev(self, p_prev)
            self.p_prev = p_prev;            
            self.p_history(:,2) = self.p_prev';
        end
    end
end

