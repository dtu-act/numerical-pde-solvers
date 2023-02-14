function [p_interp] = interpElementwise(P, x1d, xx_out, V, conn_table, etov, y)
    p_interp = zeros(1,length(xx_out));
    
    for i=1:length(xx_out)
        x_i = xx_out(i);
        
        %% INTERPOLATION OPERATOR
        nodes = x1d(1:P:end);
        elem_indxs = findElementIndex1D(nodes,x_i,etov);
        if elem_indxs == 0
            fprintf('index: %i, x: %f\n', i, x_i)
            error('ERROR: only interpolation is supported (vertice not found inside any elements)')
        end

        idx = conn_table(elem_indxs,:); % node indices in the receiving element (non-periodic mesh)
        VXel = x1d([idx(1), idx(end)])'; % first and last element node. NOTE: ordering of nodes is: [first, last, polynomial_interpol]
        r = map2elem1D(x_i,VXel); % map from physical coordinates to standard element
        if (r < -1 || r > 1)
            error('r not within range')
        end
        [V1DGpoint] = spectral.Vandermonde1D(P,r); % construct interpolation from vandermonde matrix
        I1pG = V1DGpoint*inv(V); % construct inter-element interpoliation operation.
        
        if size(y,1) > 1
            disp('hello')
        end
        
        p_interp(i) = I1pG*y(idx)';
    end
end

function [r] = map2elem1D(x,VX)
    % Map to standard 1D element [-1,1]
    r = -1 + 2*(x-VX(1))/(VX(2) - VX(1));
    
%     A = [-VX(1)+VX(2];
%     b = [2*xr - VX(2)];
% 
%     rs = A\b;
% 
%     r = rs(1);
end

function idx = findElementIndex1D(nodes,x,etov)
    idx = 0;
    for kk = 1:length(etov)
        vi1 = etov(kk,1);
        vi2 = etov(kk,2);        
    
        if x >= nodes(vi1) && x <= nodes(vi2)        
            idx = kk;
            break
        end
    end

    if idx == 0
        disp('Element not found')
    end
end