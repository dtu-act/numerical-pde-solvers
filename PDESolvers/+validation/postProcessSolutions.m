function [P_out,Pref_out,C] = postProcessSolutions(P,Pref,XY,bbox,xy0,r0,adjust_energy)
    if nargin == 6
        adjust_energy = true;
    end
    
    % remove area around point source
    fd = @(p) min( -drectangle(p,bbox(1,1),bbox(2,1),bbox(1,2),bbox(2,2)), dcircle(p,xy0(1),xy0(2),r0));
    Proi = validation.extractROI(P,XY,fd);
    Pref_roi = validation.extractROI(Pref,XY,fd);

    % validate
    if adjust_energy
        C = validation.calcNormConstantMax(Proi,Pref_roi);
    else
        C=1;
    end

    P_out = Proi*C;
    Pref_out = Pref_roi;
end