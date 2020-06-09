function outstruct = plotGRvsQ(modstruct,expdatastruct,ploton,figh)
%outstruct = plotGRvsQ(modstruct,expdatastruct,ploton,figh)
%get growth rate and total resource demand for all modules in modstruct

if nargin < 3
    ploton = true;
end
if nargin < 4 || isempty(figh)
    figh = figure;
end

%init
modnames = fieldnames(modstruct);
outstruct = struct;

%loop through modules
for ii = 1:length(modnames)
    modout = RSmodule;
    mod1 = modstruct.(modnames{ii});
    
    %get total resource demand as sum of individual resource demands
    if ~mod1.isalone
        modsinside = mod1.containingmods;
        Qtot = 0;
        Qtotvar = 0;
        for jj = 1:length(modsinside)
            %get Q for submodules
            Qtot = Qtot + mean(modstruct.(modsinside{jj}).Q);
            Qtotvar = Qtotvar^2 + sum(modstruct.(modsinside{jj}).Qstd).^2/...
                numel(modstruct.(modsinside{jj}).Qstd);
        end
        modout.Q = Qtot;
        modout.Qstd = sqrt(Qtotvar);
    else
        modout.Q = mean(mod1.Q);
        modout.Qstd = sqrt(sum(mod1.Qstd.^2)./numel(mod1.Qstd));
    end
    
    %copy growth rate info
    gr0 = expdatastruct.(modnames{ii}).growthrate;
    modout.growthrate = mean(gr0);
    modout.growthratestd = std(gr0);
    
    %output
    outstruct.(modnames{ii}) = modout;
end

if ploton
    
    [gr,grstd,Q,Qstd] = deal(zeros(length(modnames),1));
    for k = 1:length(modnames)
        %get data
        mod2 = outstruct.(modnames{k});
        gr(k) = mod2.growthrate;
        grstd(k) = mod2.growthratestd;
        Q(k) = mod2.Q;
        Qstd(k) = mod2.Qstd;
    end
    %plot
    figure(figh); clf;
	plot(gr,Q,'o','linewidth',2); hold on;
    errorbar(gr,Q,Qstd,'.k','vertical','linewidth',1.0)
    errorbar(gr,Q,grstd,'.k','horizontal','linewidth',1.0)
    set(gca,'fontsize',14,'linewidth',2)
    xlabel('Growthrate [1/hr]')
    ylabel('Resource Demand [ ]')
    xlim([0,inf])
    ylim([0,inf])
    box off
end
