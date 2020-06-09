classdef RSexpdata
    properties
        samplename
        time                                %units of hours
        BFP
        GFP
        RFP
        OD
        timestep (1,1) double = 5/60        %measurements every 5 mins
        FPfields cell = {'BFP','GFP','RFP'}
        ODfield cell = {'OD'}
        growthrate
        BFPdiff
        GFPdiff
        RFPdiff
        BFPOD
        RFPOD
        GFPOD
        BFPdiffOD
        GFPdiffOD
        RFPdiffOD
        totalQ
    end
    methods
        function obj = loaddata(obj,filename)
            platedata = importplate([filename,'.txt']);
            %load([filename,'.mat'],'platedata')
            obj.OD = platedata.OD600_600;
            obj.RFP = platedata.RFP1_584_619; %RFP1_584_607;
            obj.GFP = platedata.GFP1_485_530; %GFP1_465_498;
            obj.BFP = platedata.GFP1_400_460; %GFP1_518_540;
            %add time
            numsamples = max(structfun(@(x) size(x,3), platedata));
            timevec = obj.timestep*(0:(numsamples-1));
            if isempty(obj.time)
                obj.time = timevec(:);
            else
                timeprev = obj.time(end);
                obj.time = [obj.time(:); timevec + timeprev + obj.timestep];
            end
        end
        function obj = catRSexpdata(obj,obj2)
            
            
        end
    end
end