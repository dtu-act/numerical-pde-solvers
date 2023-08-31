function [pqhat] = sourceAtBoundary(i, XY, k, domain, boundary_type)
    switch domain.source.type
        case models.SourceType.Line
            switch boundary_type
                case models.BoundaryCondition.Pressure
                    % pressure source: Pi => pi
                    pqhat = domain.source.F(i+1,:).*ones(size(XY,1),1);
                case models.BoundaryCondition.Velocity
                    % velocity source: NABLA Pi => [dpi_x,dpi_y]
                    pqhat = domain.source.dF(i+1,:).*ones(size(XY,1),1);
                otherwise
                    error('Not implemented: source/boundary type not supported')
            end
        case models.SourceType.PointSource
            switch boundary_type
                case models.BoundaryCondition.Pressure
                    pqhat = zeros(size(XY,1),1);
                    pqhat(:) = domain.source.F(XY(:,1),XY(:,2),k);
                case models.BoundaryCondition.Velocity
                    pqhat = zeros(size(XY,1),2);
                    % https://se.mathworks.com/matlabcentral/answers/93677-how-can-i-evaluate-the-derivatives-of-a-bessel-function-at-different-points
                    pqhat(:,:) = domain.source.dF(XY(:,1),XY(:,2),k);
                otherwise
                    error('Not implemented: source/boundary type not supported')
            end
        case models.SourceType.None
             switch boundary_type
                case models.BoundaryCondition.Pressure
                    pqhat = zeros(size(XY,1),1);
                case models.BoundaryCondition.Velocity
                    pqhat = zeros(size(XY,1),2);
                otherwise
                    error('Not implemented: source/boundary type not supported')
             end            
        otherwise
            error('Not implemented: source/boundary type not supported')
    end
end