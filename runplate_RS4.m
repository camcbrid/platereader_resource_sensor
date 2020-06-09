addpath utility
%initialize for data
filenamecell = {'data 2-4-20\J0_exp_2-4-2020_1'};
% %match wells on plate with experimental conditions
indstruct = cell(length(filenamecell),1);
indstruct{1}.B = [2,2; 2,3; 2,4];
indstruct{1}.G = [3,2; 3,3; 3,4];
indstruct{1}.R = [4,2; 4,3; 4,4];
indstruct{1}.BG = [5,2; 5,3; 5,4];
indstruct{1}.BR = [6,2; 6,3; 6,4];
indstruct{1}.GR = [7,2; 7,3; 7,4];
indstruct{1}.M9 = [2,5; 3,5; 4,5];
%tell which fields should be zero to subtract background
nonestruct = struct;
nonestruct.OD = {'M9'};
nonestruct.BFP = {'G','R','GR'};
nonestruct.GFP = {'B'};%,'R','BR'};
nonestruct.RFP = {'B','G','BG'};
[cellholder,basalstruct] = deal(cell(length(filenamecell),1));
%windows where each cell type are in steady state
ssstruct.B = [6.1,7.8];
ssstruct.G = [4.0,5.8];
ssstruct.R = [6.0,7.1];
ssstruct.BG = [6.5,8.0];
ssstruct.BR = [5.1,7.8];
ssstruct.GR = [4.8,6.0];
%metadata struct for matching module properties
metadata.B.FPout = {'BFPdiffOD'};
metadata.B.containingmods = {'B'};
metadata.B.isResourceSensor = true;
metadata.B.isalone = true;
metadata.R.FPout = {'RFPdiffOD'};
metadata.R.containingmods = {'R'};
metadata.R.isResourceSensor = true;
metadata.R.isalone = true;
metadata.G.FPout = {'GFPdiffOD'};
metadata.G.containingmods = {'G'};
metadata.G.isalone = true;
metadata.G.isResourceSensor = false;
%modules together
metadata.BR.FPout = {'BFPdiffOD','RFPdiffOD'};
metadata.BR.containingmods = {'B','R'};
metadata.BR.isResourceSensor = true;
metadata.BR.isalone = false;
metadata.BG.FPout = {'BFPdiffOD','GFPdiffOD'};
metadata.BG.containingmods = {'B','G'};
metadata.BG.isResourceSensor = false;
metadata.BG.isalone = false;
metadata.GR.FPout = {'GFPdiffOD','RFPdiffOD'};
metadata.GR.containingmods = {'G','R'};
metadata.GR.isResourceSensor = false;
metadata.GR.isalone = false;

%-----------------------------------
%analyze data
for ii = 1:length(filenamecell)
    %init data and load data from file
    platedata = RSexpdata;
    platedata.samplename = filenamecell{ii};
    platedata = loaddata(platedata,filenamecell{ii});
    platedata = filtplatedata3(platedata,[]);
    %repackage data into individual fields for each experimental condition
    cellstruct0 = bulkdata2cells2(platedata,indstruct{ii});
    %subtract the basal fluorescence and OD off
    [cellstruct0, basalstruct{ii}] = subtractbasal5(cellstruct0,nonestruct);
    %take time derivative of raw fluorescent data
    cellstruct0 = diffcellfields(cellstruct0,{'BFP','GFP','RFP'},'time');
    %fit growthrate
    cellstruct0 = fitgrowthrate3(cellstruct0,false,false,1);
    %combine batches
    cellholder{ii} = cellstruct0;
end

%combine batches into one
cellstruct = combinedilutions(cellholder{1});       %note: this may return a struct of structs, FIX
%divide fluorescence by OD and growthrate
cellstruct = normalizeFP3(cellstruct,{'BFP','GFP','RFP','BFPdiff','GFPdiff','RFPdiff'});
%average over time intervals to get steady state
modulestruct = findss4(cellstruct,ssstruct,{'BFPdiffOD','GFPdiffOD','RFPdiffOD'},@median,false);
%add metadata to modules
modulestruct = addmetadata(modulestruct,metadata);
%calculate resource demand for resource sensor
RSmodules = calcRS2(modulestruct,false,1,2);
%find Q and S for non-resource sensor modules
[modulesout,RSmodules] = calcQS(modulestruct,RSmodules,false,3,4,5,6);
%combine modules together
modulesout = combinemodstructs(RSmodules,modulesout,modulestruct);
%plot growth rate and resource demand
foo = plotGRvsQ(modulesout,cellstruct,true,7);

return

%make predictions of combos from measurements alone
out = predictRS3(modulestruct,modulesout,true,7);

%plot?
figure(1);
plotsubfield2(cellstruct,'OD','BFP',8,{'B','G','R','BG','BR','GR'})
