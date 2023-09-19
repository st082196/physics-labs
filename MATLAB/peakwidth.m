function width = peakwidth(x,y,cutoff)
%Evaluate peak width
%
%   x [vector]
%   x data
%
%   y [vector]
%   y data
%
%   cutoff [0.5 (default) | scalar value in the range (0,1)]
%   Fraction of the peak height at the boundaries of the width

if nargin < 3
    cutoff = 0.5;
end
if ~issorted(x)
    [x,i] = sort(x);
    y = y(i);
end

cutoff = cutoff*max(y);
included_points = y > cutoff;

lb = x(1);
i = find(included_points,1,'first');
if i ~= 1
    lb = interp1(y([i-1, i]), x([i-1, i]), cutoff);
end

ub = x(end);
i = find(included_points,1,'last');
if i ~= numel(y)
    ub = interp1(y([i, i+1]), x([i, i+1]), cutoff);
end

width = ub - lb;