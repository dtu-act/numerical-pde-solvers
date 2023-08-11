function [pointSource2D] = PointSource(r0,Q)
%POINTSOURCE2D Initialize a SourceModel struct
% source_type: SourceType enum
% source: Union{Matrix{T},Vector{T},PointSource2D,Function} where T <: Number

    pointSource2D = struct(...
        'type', models.SourceType.PointSource,...
        'r0', r0,...
        'Q', Q);
end