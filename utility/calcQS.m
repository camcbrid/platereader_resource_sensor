function [modsout,RSmodsout] = calcQS(modstruct,RSmodules,ploton,fighy,fighyRS,fighQ,fighS)
%[modsout,RSmodsout] = calcQS(modstruct,RSmodules,ploton,fighy,fighyRS,fighQ,fighS)
%calculate resource demand Q and sensitivity S for arbitrary modules under
%the perturbation from a resource sensor.

if nargin < 3
    ploton = false;
end
if ploton
    if nargin < 4 || isempty(fighy)
        fighy = figure;
    elseif ishandle(fighy); clf(fighy);
    end
    if nargin < 5 || isempty(fighy)
        fighyRS = figure;
    elseif ishandle(fighyRS); clf(fighyRS);
    end
    if nargin < 6 || isempty(fighQ)
        fighQ = figure;
    elseif ishandle(fighQ); clf(fighQ);
    end
    if nargin < 7 || isempty(fighS)
        fighS = figure;
    elseif ishandle(fighS); clf(fighS);
    end
end

%find groups of non-RS modules that share properties
modsalone = findprop(modstruct,'isalone');
modstogether = findprop(modstruct,'isalone',false);
normalmods = findprop(modstruct,'isResourceSensor',false);  %not resource sensors
RSmods = findprop(modstruct,'isResourceSensor',true);

%experiments that are modules together and one module isn't a resource sensor
normalmodstogeth = intersect(normalmods,modstogether);
normalmodsalone = intersect(normalmods,modsalone);
RSmodsalone = intersect(RSmods,modsalone);
RSmodstructplt = struct;

%loop through perturbed modules to calculate Q and S for each under the
%perturbation from the Resource Sensor
for ii = 1:length(normalmodstogeth)
    %get the modules in the pair
    pairmods = modstruct.(normalmodstogeth{ii});
    containedmods = pairmods.containingmods;
    
    RSmodname = [];
    modname = [];
    %loop through contained modules to find which modules are resource
    %sensors and which are normal
    for jj = 1:length(containedmods)
        %exactly one module needs to be a resource sensor
        if modstruct.(containedmods{jj}).isResourceSensor
            RSmodname = containedmods{jj};
        else
            modname = containedmods{jj};
        end
    end
    
    %need at least one resource sensor and one normal module
    if ~isempty(modname) && ~isempty(RSmodname)
        mod1 = modstruct.(modname);
        RSmod1 = RSmodules.(RSmodname);
        
        %get fluorescence channels corresponding to outputs
        mod_FPoutname = mod1.FPout;
        RS_FPoutname = RSmod1.FPout;
        
        %get module and resource sensor outputs
        modalone = mod1.(mod_FPoutname{1});
        modalonestd = mod1.([mod_FPoutname{1},'std']);
        modtogeth = pairmods.(mod_FPoutname{1});
        modtogethstd = pairmods.([mod_FPoutname{1},'std']);
        %for resource sensor
        RSalone = RSmod1.(RS_FPoutname{1});
        RSalonestd = RSmod1.([RS_FPoutname{1},'std']);
        RStogeth = pairmods.(RS_FPoutname{1});
        RStogethstd = pairmods.([RS_FPoutname{1},'std']);
        
        %get resource sensor Q
        QRS = RSmod1.Q;
        QRSstd = RSmod1.Qstd;
        
        %find Q and S for module from data
        [Q,S,Qstd,Sstd] = calcQS0(modalone,RSalone,modtogeth,RStogeth,QRS,...
            modalonestd,RSalonestd,modtogethstd,RStogethstd,QRSstd);
        
        %output updated modules
        mod1.Q = [mod1.Q, Q];
        mod1.Qstd = [mod1.Qstd, Qstd];
        mod1.S = [mod1.S, S];
        mod1.Sstd = [mod1.Sstd, Sstd];
        mod1.y = modalone;
        mod1.ystd = modalonestd;
        mod1.perturbname = [mod1.perturbname, {RSmodname}];
        mod1.perturby = [mod1.perturby, modtogeth];
        mod1.perturbystd = [mod1.perturbystd, modtogethstd];
        modstruct.(modname) = mod1;
        
        %output for resource sensor
        %only for plotting
        RSmod2 = modstruct.(RSmodname);
        RSmod2.Q = [RSmod2.Q, QRS];
        RSmod2.Qstd = [RSmod2.Qstd, QRSstd];
        RSmod2.S = [RSmod2.S, 1];
        RSmod2.Sstd = [RSmod2.Sstd, 0];
        RSmod2.y = RSalone;
        RSmod2.ystd = RSalonestd;
        RSmod2.perturbname = [RSmod2.perturbname, {modname}];
        RSmod2.perturby = [RSmod2.perturby, RStogeth];
        RSmod2.perturbystd = [RSmod2.perturbystd, RStogethstd];
        RSmodstructplt.(RSmodname) = RSmod2;
        %total output, keeping initial perturbation data from resource
        %sensors
        RSmod1.perturbname = [RSmod1.perturbname, {modname}];
        RSmod1.perturby = [RSmod1.perturby, RStogeth];
        RSmod1.perturbystd = [RSmod1.perturbystd, RStogethstd];
        modstruct.(RSmodname) = RSmod1;
    end
end

%output data only for modules that are not resource sensors
modsout = struct;
for k = 1:length(normalmodsalone)
    modsout.(normalmodsalone{k}) = modstruct.(normalmodsalone{k});
end
for q = 1:length(RSmodsalone)
    RSmodsout.(RSmodsalone{q}) = modstruct.(RSmodsalone{q});
end

if ploton
    %plot mean and std of RS sensor outputs y and Q
    normmodnames = fieldnames(modsout);
    [yvec,ystdvec,Qvec,Qvecstd,Svec,Svecstd,bh1,bh2,ebh2] = ...
        deal(cell(length(normmodnames),1));
    %loop through resource sensors to reshape data for plotting
    for k = 1:length(normmodnames)
        mod2 = modsout.(normmodnames{k});
        yalone = mod2.y;
        yalonestd = mod2.ystd;
        yperturb = mod2.perturby;
        yperturbstd = mod2.perturbystd;
        %store normalized output for plotting
        yvec{k} = [yalone; yperturb(:)]/yalone;
        ystdvec{k} = [yalonestd; yperturbstd(:)]/yalone;
        %store resource demand for plotting
        Qvec{k} = mod2.Q;
        Qvecstd{k} = mod2.Qstd;
        Svec{k} = mod2.S;
        Svecstd{k} = mod2.Sstd;
    end
    
    %reshape data for resource sensor outputs when perturbed by modules for
    %plotting
    RSmodnames = fieldnames(RSmodstructplt);
    [RSyvec,RSystdvec,bh3,ebh3] = deal(cell(length(RSmodnames),1));
    for p = 1:length(RSmodnames)
        RSmod2 = RSmodstructplt.(RSmodnames{p});
        RSyalone = RSmod2.y;
        RSyalonestd = RSmod2.ystd;
        RSyperturb = RSmod2.perturby;
        RSyperturbstd = RSmod2.perturbystd;
        RSyvec{p} = [RSyalone; RSyperturb(:)]/RSyalone;
        RSystdvec{p} = [RSyalonestd; RSyperturbstd(:)]/RSyalone;
    end
    
    %plot module outputs
    figure(fighy);
    xprev = 0;
    for m = 1:length(normmodnames)
        %get x coordinates to plot bars
        xvec = (1:length(yvec{m})) + xprev;
        xprev = xvec(end);
        %plot bars
        bh2{m} = bar(xvec,yvec{m}); hold on
        set(bh2{m},'edgecolor','none','barwidth',0.6);
        ebh2{m} = errorbar(xvec,yvec{m},ystdvec{m},'.k','linewidth',1.5);
        %xticklabels(RSnames);
        set(gca,'fontsize',16);
    end
    ylabel('$y$','interpreter','latex')
    box off
    set(gca,'linewidth',2)
    xlim([0.3,xprev+0.8])
    
    %plot module outputs
    figure(fighyRS);
    xprev2 = 0;
    for m = 1:length(RSmodnames)
        %get x coordinates to plot bars
        xvec2 = (1:length(RSyvec{m})) + xprev2;
        xprev2 = xvec2(end);
        %plot bars
        bh3{m} = bar(xvec2,RSyvec{m}); hold on
        set(bh3{m},'edgecolor','none','barwidth',0.6);
        ebh3{m} = errorbar(xvec2,RSyvec{m},RSystdvec{m},'.k','linewidth',1.5);
        %xticklabels(RSnames);
        set(gca,'fontsize',16);
    end
    ylabel('$y$ RS','interpreter','latex')
    box off
    set(gca,'linewidth',2)
    xlim([0.3,xprev2+0.8])
    
    %plot Q for module
    figure(fighQ);
    xprev3 = 0;
    for l = 1:length(normmodnames)
        xvec3 = (1:length(Qvec{l})) + xprev3;
        xprev3 = xvec3(end);
        %make sepperate bar objects for each resoruce sensor module to
        %change colors independently
        bh1{l} = bar(xvec3,Qvec{l}); hold on
        set(bh1{l},'edgecolor','none','barwidth',0.6);
        ebh2{l} = errorbar(xvec3,Qvec{l},Qvecstd{l},'.k','linewidth',1.5);
        set(gca,'fontsize',16)
        xticklabels(normmodnames);
    end
    ylabel('$Q$','Interpreter','Latex');
    box off
    set(gca,'linewidth',2)
    xlim([0.3,xprev3+0.8])
    
    %plot Q for module
    figure(fighS);
    xprev4 = 0;
    for l = 1:length(normmodnames)
        xvec4 = (1:length(Svec{l})) + xprev4;
        xprev4 = xvec4(end);
        %make sepperate bar objects for each resoruce sensor module to
        %change colors independently
        bh1{l} = bar(xvec4,Svec{l}); hold on
        set(bh1{l},'edgecolor','none','barwidth',0.6);
        ebh2{l} = errorbar(xvec4,Svec{l},Svecstd{l},'.k','linewidth',1.5);
        set(gca,'fontsize',16)
        xticklabels(normmodnames);
    end
    ylabel('$S$','Interpreter','Latex');
    xlim([0.3,xprev4+0.8])
    box off
    set(gca,'linewidth',2)
    
    %NOTE: green color = [163,252,67]/255
    %blue color = [125,195,251]/255
    %red color = [251,122,122]/255
end


function [Q,S,Qstd,Sstd] = calcQS0(modalone,RSalone,modtogeth,RStogeth,QRS,...
    modalonestd,RSalonestd,modtogethstd,RStogethstd,QRSstd)
%calculate Q and S according to [CITE]

%module resource demand
Q = ((RSalone./RStogeth)-1)*(1+QRS);
%module sensitivity = Fhat
S = modtogeth*RSalone*(1 + QRS)./...
    (modalone*(RStogeth + QRS*(RSalone - RStogeth)));

%calculate partial derivatives for standard error calculation
dQdRSa = (1+QRS)/RStogeth;                 %dQ/dRSalone
dQdRSt = -(1+QRS)*(RSalone./RStogeth^2);   %dQ/dRStogether
dQdQRS = ((RSalone./RStogeth)-1);          %dQ/dQRS
dSdma = -modtogeth*RSalone*(1+QRS)./...
    (modalone^2*(RStogeth+QRS*(RSalone-RStogeth)));     %dS/dmodalone
dSdmt = RSalone*(1+QRS)./...
    (modalone*(RStogeth+QRS*(RSalone-RStogeth)));       %dS/dmodtogether
dSdRSa = -modtogeth*RStogeth*(1+QRS)*(QRS-1)/...
    (modalone*(RStogeth+QRS*(RSalone-RStogeth))^2);     %dS/dRSalone
dSdRSt = modtogeth*RSalone*(1+QRS)*(QRS-1)/...
    (modalone*(RStogeth+QRS*(RSalone-RStogeth))^2);     %dS/dRStogether
dSdQRS = modtogeth*RSalone*(RSalone-2*RStogeth)/...
    (modalone*(RStogeth+QRS*(RSalone-RStogeth))^2);     %dS/dQRS

%calculate standard error
Qstd = sqrt(dQdRSa^2*RSalonestd^2 + dQdRSt^2*RStogethstd^2 + dQdQRS^2*QRSstd^2);
Sstd = sqrt(dSdma^2*modalonestd^2 + dSdmt^2*modtogethstd^2 + ...
    dSdRSa^2*RSalonestd^2 + dSdRSt^2*RStogethstd^2 + dSdQRS^2*QRSstd^2);
