function [sourceModel] = SourceModelNone()
%SOURCEMODEL Initialize a SourceModel struct with type 'None'    
    sourceModel = struct('type', models.SourceType.None);
end