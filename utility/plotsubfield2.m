function [figh,axh] = plotsubfield2(cellstruct, xfield, yfield,figh, cellfields)
%[figh,axh] = plotsubfield2(cellstruct, xfield, yfield,figh, cellfields)
%plot one subfield for each field in a struct. datafield must match to a
%valid field for all substructs contined in cellstruct

%input verification

if nargin < 4 || isempty(figh)
    figh = figure;
end

fntsze = 12;

%get fields for each cell type
celltypes = fieldnames(cellstruct);
if nargin > 4 && ~isempty(cellfields)
    cellnames = intersect(cellfields,celltypes,'stable');
else; cellnames = celltypes; 
end

n = length(cellnames);                  %number of subplots (fields) to make
m = ceil(9/16*sqrt(n));                 %number of rows of subplots

%plot each struct's subfield on one plot
figure(figh);
axh = cell(n,1);
for jj = 1:n
    %make subplot
    axh{jj} = subplot(m,ceil(n/m),jj);
    %keep hold on if already set
    if ishold; hold on; end
    %plot
    plot(cellstruct.(cellnames{jj}).(xfield), cellstruct.(cellnames{jj}).(yfield),'linewidth',2);
    title(strrep(cellnames{jj},'_','-'))
    %put ylabels only on farthest left plots
    if mod(jj,ceil(n/m)) == 1; ylabel(strrep(yfield,'_',' ')); end
    %put ylabels only on bottom plots
    if n - jj < ceil(n/m); xlabel(strrep(xfield,'_',' ')); end
    set(gca,'fontsize',fntsze)
    hold off
end
