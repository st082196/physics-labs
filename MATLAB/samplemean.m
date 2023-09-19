function [x, dx] = samplemean(sample,e,level)
%Calculate mean of the sample and error of the mean
%
%   sample [vector | matrix | cell array of vectors]
%   Measurement data
%
%   e [0 (default) | scalar | vector]
%   Systematic error
%
%   level [0.95 (default) | scalar value in the range (0,1)]
%   Confidence level
%
%If sample is a matrix, x and dx are column vectors containing the means
%and errors corresponding to each row.
%If sample is a cell array, x and dx are matrices containing the means and
%errors corresponding to each cell.

if nargin < 3
    level = 0.95;
    if nargin < 2
        e = 0;
    end
end

if iscell(sample)
    e = e.*ones(size(sample));
    [x, dx] = deal(zeros(size(sample)));
    for i = 1:numel(sample)
        n = numel(sample{i});
        x(i) = mean(sample{i});
        dx(i) = sqrt(tinv(0.5*(1+level),n-1).^2.*var(sample{i})/n + e(i).^2);
    end
elseif isvector(sample)
    n = numel(sample);
    x = mean(sample);
    dx = sqrt(tinv(0.5*(1+level),n-1).^2.*var(sample)/n + e.^2);
else
    n = size(sample,2);
    x = mean(sample,2);
    dx = sqrt(tinv(0.5*(1+level),n-1).^2.*var(sample,0,2)/n + e.^2);
end