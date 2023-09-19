function varargout = matsplit(A,dim)
%Split matrix elements into separate variables.
%
%   A [any type of multidimensional array]
%   Matrix to be split
%
%   dim (optional) [positive integer | positive vector of integers]
%   Dimension preserved when splitting

if nargin==1
    varargout = num2cell(A);
else
    varargout = num2cell(A,dim);
end