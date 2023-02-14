function [out_function] = createBoundRectDetectFunc(bounding_box, tol)
%
%    Purpose:  creates a function for detect node types on a 
%              rectangular domain (staircase type only): 
%              -1 is outer node
%               0 is boundary
%               1 is inner node
%
%    Return:   f(real, real) -> int
%

    xmin = bounding_box(1,1);
    xmax = bounding_box(2,1);
    ymin = bounding_box(1,2);
    ymax = bounding_box(2,2);
    
    function [val] = boundaryDetectionFunc(x,y)

        lrbound = @(x) abs(x - xmin) < tol || abs(x - xmax) < tol;
        ulbound = @(y) abs(y - ymax) < tol || abs(y - ymin) < tol;
        
        ulinner = @(y) y - ymin + tol > 0 && y - ymax - tol < 0;
        lrinner = @(x) x - xmin + tol > 0 && x - xmax - tol < 0;
        
        if lrbound(x) && ulinner(y)
            % left/right boundary and within upper lower boundaries
            val = 0;
            return
        end
        
        if ulbound(y) && lrinner(x)
            % upper/lower node and within left/right boundaries
            val = 0;
            return
        end
        
        % inner?
        if lrinner(x) && ulinner(y)
            val = 1;
            return
        end

        % outer node
        val = -1;
    end

    out_function = @boundaryDetectionFunc;
end