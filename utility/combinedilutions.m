function outstruct = combinedilutions(datastruct)
%HELP need to fix for RSexpdata objects

outstruct = struct;

if iscell(datastruct)
    %loop across each element in the cell array
    for ii = 1:length(datastruct)
        %concatenate fields of each cell in cellstruct
        outstruct = catstructs(datastruct{ii},outstruct);
    end
else
    outstruct = datastruct;
end


function outstruct = catstructs(instruct,outstruct)
%outstruct = catstructs(instruct,outstruct)
%concatenate fields within the two structures instruct and outstuct

if isa(instruct,'RSexpdata')
    names = fieldnames(instruct);
    %loop across subfields
    for ii = 1:length(names)
        %recursion
        if isfield(outstruct,names{ii})
            obj = RSexpdata;
            obj = catstructs(instruct.(names{ii}),outstruct.(names{ii}));
            outstruct.(names{ii}) = obj;
        else
            
            outstruct.(names{ii}) = RSexpdata;
            outstruct.(names{ii}) = catstructs(instruct.(names{ii}),[]);
        end
    end
elseif isstruct(instruct)
    names = fieldnames(instruct);
    %loop across subfields
    for ii = 1:length(names)
        %recursion
        if isfield(outstruct,names{ii})
            outstruct.(names{ii}) = catstructs(instruct.(names{ii}),outstruct.(names{ii}));
        else
            outstruct.(names{ii}) = catstructs(instruct.(names{ii}),[]);
        end
    end
elseif ischar(instruct)
    %append character arrays into cell array
    outstruct = [outstruct;{instruct}];
elseif isnumeric(instruct)
    %concatenate numeric data along the first dimension
    outstruct = cat(1,outstruct,instruct);
else
    outstruct = instruct;
end


