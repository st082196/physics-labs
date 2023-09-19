function printangle(x,dx,label)
%Print x and dx to Command Window formatted as angle values and
%error margins in degrees-minutes-seconds representation
%
%   x [scalar | column vector | n-by-3 matrix]
%   Angles in degrees or degrees-minutes-seconds representation
%
%   dx (optional) [scalar | column vector | n-by-3 matrix]
%   Error margins in degrees or degrees-minutes-seconds representation
%
%   label (optional) [string]
%   Label of x
%
%If x and dx are scalars, prints out 'label = dd°mm′ss″ ± dd°mm′ss″'
%If x and dx are vectors, prints out label followed by multiple lines
%of 'dd°mm′ss″ ± dd°mm′ss″'.
%If dx is specified as 0, the corresponding x is printed without error
%margins.
%If dx contains single angle value while x has multiple values,
%dx is assumed to be the same for every x and is printed separately.

% Check input arguments
if ~(ismatrix(x) && any(size(x,2) == [1, 3]))
    error('x must be scalar, column vector or 3-column matrix');
elseif nargin >= 2 && ~(ismatrix(dx) && any(size(dx,2) == [1, 3]))
    error('dx must be scalar, column vector or 3-column matrix');
elseif nargin >= 2 && ~any(size(dx,1) == [1, size(x,1)])
    error('dx must have 1 or the same number of rows as x');
elseif nargin >= 2 && any(dx < 0,'all')
    error('dx must contain non-negative values');
end

if nargin < 2
    dx = 0;
end
if size(x,2) == 1
    x = degrees2dms(x);
end
if size(dx,2) == 1
    dx = degrees2dms(dx);
end

% Format output to Command Window
fprintf('\n');
if size(x,1) == 1
    if nargin >= 3
        fprintf('%s = ',label);
    end
else
    fprintf('%s\n',label);
end
for i = 1:size(x,1)
    if any(sign(x(i,:)) == -1)
        fprintf('-');
    end
    fprintf('%d°%02d′%02.0f″',abs(x(i,:)));
    if size(dx,1) == size(x,1) && any(dx(i,:) ~= 0)
        fprintf(' ± %d°%02d′%02.0f″',dx(i,:));
    end
    fprintf('\n');
end
if size(dx,1) ~= size(x,1) && any(dx ~= 0)
    fprintf('\n∆');
    if nargin >= 3
        fprintf('%s',label);
    end
    fprintf(' = %d°%02d′%02.0f″\n',dx);
end