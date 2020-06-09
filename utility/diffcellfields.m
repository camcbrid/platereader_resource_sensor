function outstruct = diffcellfields(cellstruct,difffields,timefield,n)
%outstruct = diffcellstruct(cellstruct,difffields,timefield)
%apply smooth differentation to difffields in cellstruct with respect to
%timefield. difffields may be a cell array containing the fields to be
%differentaited. timefield must be a string.

cellnames = fieldnames(cellstruct);
outstruct = cellstruct;

if nargin < 4 || isempty(n)
    n = 3;          %moving average window size
end

%loop through cells
for ii = 1:length(cellnames)
    %loop through fields to be differentiated
    for jj = 1:length(difffields)
        %get fields
        FP = cellstruct.(cellnames{ii}).(difffields{jj});
        time = cellstruct.(cellnames{ii}).(timefield);
        %take filtered numerical derivative with respect to time
        dFPdt0 = mydiff(FP,time);
        %moving average filter
        dFPdt = moveavefilt(dFPdt0,n);
        %output
        outstruct.(cellnames{ii}).([difffields{jj},'diff']) = dFPdt;
    end
end


function dfdh = mydiff(f,h)
%run 9-pt Savitzky-Golay differentiation on x with respect to t
%f is vector of inputs, h is step size. h can be vector of evenly spaced
%points where f is sampled.

%find step size
if length(h) > 1
    h2 = abs(mean(h(2:end) - h(1:end-1)));
else; h2 = h;
end

% 9-pt Savitzky-Golay filtered first derivative (1st/2nd order polynomial fit)
a = 60*h2;
b = [-4,-3,-2,-1,0,1,2,3,4];
dfdh = filter(-b,a,f);


function y = moveavefilt(x,n)
%run zero phase moving average filter with window of n points
a = 1;
b = 1/n*ones(1,n);
if nnz(isnan(x)) > 0
    y = x;
else
    %zero phase moving average filter
    y = filtfilt(b,a,x);
end
