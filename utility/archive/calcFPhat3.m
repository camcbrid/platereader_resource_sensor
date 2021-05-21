function [FPhat,FPhatstd] = calcFPhat3(FPalone,FPtogether,FPalonestd,FPtogetherstd)
%calculates ratios of fluorescent concentrations.

%calculate FPhat
FPhat = FPtogether./FPalone;

%estimate error bars assuming uncertianty is dominanted by variance across
%population samples and assuming var(FPtogether), var(FPalone) are uncorrelated
FPhatstd = abs(FPhat).*sqrt( (FPtogetherstd./FPtogether).^2 + ...
    (FPalonestd./FPalone).^2);