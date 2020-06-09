function plotstruct(cellstruct,xfield,yfield,figh,cellfields)
%plot two subfields from a struct for all fields

if nargin < 4 || isempty(figh)
    figh = figure;
end
if nargin < 5 || isempty(cellfields)
    cellnames = fieldnames(cellstruct);
else
    cellnames = cellfields;
end

figure(figh); clf;
for ii = 1:length(cellnames)
    %make plots
    plot(cellstruct.(cellnames{ii}).(xfield),...
        cellstruct.(cellnames{ii}).(yfield),'linewidth',2);
    hold on
end
legend(cellnames,'Location','Best')
hold off
xlabel(xfield)
ylabel(yfield)
set(gca,'fontsize',12)