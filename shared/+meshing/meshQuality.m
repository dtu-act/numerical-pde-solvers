function q = meshQuality(etov)
% - Equilateral triangles have q = 1.
% - Degenerate triangles have q = 0.
% - "Good" triangles we defined as having q > 0.5 (rule of thumb).
%
%   From notes "The Spectral/hp-Finite Element Method for Partial
%   Differential Equations" by by Allan P. Engsig-Karup
%
    % Compute side lengths of triangles
    lx = VX(etov(:,[1 2 3])) - VX(etov(:,[3 1 2]));
    ly = VY(etov(:,[1 2 3])) - VY(etov(:,[3 1 2]));
    l = sqrt(lx.^2+ly.^2);
    a = l(:,1); b = l(:,2); c = l(:,3);
    % Compute triangle measures
    q = (b+c-a).*(c+a-b).*(a+b-c)./(a.*b.*c);
end