classdef RSmodule
    properties
        u = []
        FPout cell
        containingmods cell
        BFPdiffOD
        GFPdiffOD
        RFPdiffOD
        BFPdiffODstd
        GFPdiffODstd
        RFPdiffODstd
        y
        ystd
        Q
        Qls
        perturbname
        perturby
        perturbystd
        S
        Sstd
        Qstd
        Qlsstd
        isResourceSensor logical = false
        isalone logical
        growthrate
        growthratestd
    end
    methods
        function plotmoduleIO(RSmodule,varargin)
            if length(RSmodule.u) > 1
                plot(RSmodule.u,RSmodule.y,varargin{:})
                if ~isempty(RSmodule.ystd)
                    hold on;
                    errorbar(RSmodule.u,RSmodule.ystd,'.k')
                end
            else
                bar(RSmodule.y,varargin{:})
                errorbar(RSmodule.ystd,'.k')
            end
        end
    end
end