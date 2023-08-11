function [sourceModel] = SourceModel(source_type, F, r0)
%SOURCEMODEL Initialize a SourceModel struct

    if nargin == 2
        r0 = 'NA';
    end
    
    sourceModel = struct(...
        'type', source_type,...
        'F', F,...
        'r0', r0);
end