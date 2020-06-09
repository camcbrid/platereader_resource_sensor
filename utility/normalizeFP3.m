function outstruct = normalizeFP3(instruct,normfields)
%outstruct = normalizeFP3(instruct,normfields)
%normalize fluorescence data by OD in the fields contained in the cell
%array normfields.

if nargin < 2
    %fields to divide by OD
    normfields = {'BFP','GFP','RFP','BFPdiff','GFPdiff','RFPdiff'};
end

%normalize fluorescence data based on growth rate or by number of cells (OD)
cellnames = fieldnames(instruct);
%loop through cell names
for k = 1:length(cellnames)
    
    data = instruct.(cellnames{k});
    ODfield = data.ODfield{1};
    dataout = data;
    
    %props = fieldnames(instruct.(cellnames{k}));
    for ii = 1:length(normfields)
        %loop through properties of each cell
        if isprop(data,normfields{ii})
            %if fluorescence data, normalize by growth rate and OD
            %normalize FP per cell based on OD
            dataout.([normfields{ii},'OD']) = data.(normfields{ii})./data.(ODfield);
        end
    end
    %copy to output
    outstruct.(cellnames{k}) = dataout;
end
