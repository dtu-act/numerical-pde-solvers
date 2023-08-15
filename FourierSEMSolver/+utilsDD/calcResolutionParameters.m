function [dx1,dx2,dt1,dt2,nmodes1,nmodes2] = calcResolutionParameters(adjust, dx, c, cfl, nmodes)
    if nargin == 4
        nmodes = 0;
    end
    
    switch adjust
        case "none"
            dx1 = dx;
            dx2 = dx;
            
            nmodes1 = nmodes;
            nmodes2 = nmodes;
            
            dt1 = dx1/(c*cfl);
            dt2 = dt1;  
        case "spatialAdjustment"
            dx1 = dx*2;
            dx2 = dx;
            
            nmodes1 = nmodes/2;
            nmodes2 = nmodes;
            
            dt1 = min(dx1,dx2)/(c*cfl);
            dt2 = dt1;
        case "temporalAdjustment"
            dx1 = dx;
            dx2 = dx;
            
            nmodes1 = nmodes;
            nmodes2 = nmodes;
            
            dt1 = dx1/(c*cfl);
            dt2 = dx2/(c*cfl)/2; % twice the resolution
        case 'spatialTemporalAdjustment'
            dx1 = dx*2;
            dx2 = dx;
            
            nmodes1 = nmodes/2;
            nmodes2 = nmodes;
            
            % Temporal resolution also differs due spatial resolution
            dt1 = dx1/(c*cfl);
            dt2 = dx2/(c*cfl);
        otherwise
            error("Compare type not implemented")
    end  
end