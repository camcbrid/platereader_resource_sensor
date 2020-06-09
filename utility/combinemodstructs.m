function modstructout = combinemodstructs(modstruct1,modstruct2,varargin)
%modstructout = combinemodstructs(modstruct1,modstruct2,varargin)
%combine structs containing the RSmodule objects. Appends fields together.

if ~isempty(varargin)
    for jj = 1:length(varargin)
        modstructout = combinemodstructs(modstruct1,modstruct2);
        modstructout = combinemodstructs(modstructout,varargin{jj});
    end
else
    %init
    modstructout = modstruct1;
    cellnames2 = fieldnames(modstruct2);
    
    %append data from modstruct2 onto modstruct1.
    for ii = 1:length(cellnames2)
        if ~isfield(modstructout,cellnames2{ii})
            modstructout.(cellnames2{ii}) = modstruct2.(cellnames2{ii});
        else
            %append objects
            warning(['module names overlap, appending for cell ',cellnames2{ii}])
            mod1 = modstructout.(cellnames2{ii});
            mod2 = modstruct2.(cellnames2{ii});
            modstructout.(cellnames2{ii}) = appendRSmods(mod1,mod2);
        end
    end
end


function modout = appendRSmods(mod1,mod2)
%append RSmodule objects

modout = mod1;
metafields = {'u','FPout','containingmods','isResourceSensor','isalone'};
mod2names = fieldnames(mod2);

%loop through fields
for ii = 1:length(mod2names)
    if any(strcmp(metafields,mod2names{ii}))
        %just copy fields
        if iscell(modout.(mod2names{ii}))
            if ~isempty(setdiff(modout.(mod2names{ii}),mod2.(mod2names{ii})))
                warning(['metafield ',mod2names{ii},' does not match; appending fields'])
                modout.(mod2names{ii}) = [modout.(mod2names{ii}),mod2.(mod2names{ii})];
            end
        else
            if all(size(modout.(mod2names{ii})) ~= size(mod2.(mod2names{ii}))) || ...
                ~all(modout.(mod2names{ii}) == mod2.(mod2names{ii}))
                warning(['metafield ',mod2names{ii},' does not match; appending fields'])
                modout.(mod2names{ii}) = {modout.(mod2names{ii}),mod2.(mod2names{ii})};
            end
        end
    else
        %append fields
        if all(size(modout.(mod2names{ii})) ~= size(mod2.(mod2names{ii}))) || ...
            all(modout.(mod2names{ii}) ~= mod2.(mod2names{ii}))
            modout.(mod2names{ii}) = [modout.(mod2names{ii}), mod2.(mod2names{ii})];
        end
    end
end
