function printresult(x,dx,mode,label,units,p)
%Print the data to Command Window formatted as quantity with
%value, uncertainty, exponent, label and units
%
%   x [scalar | vector]
%   Values
%
%   dx (optional) [scalar | vector]
%   Values, the meaning of which depends on the specified mode
%
%   mode (optional)
%   ['(uncertainty)' (default) | 'pmuncertainty' | 'significant' | 'decimals']
%   Output mode
%   (uncertainty): dx is uncertainty of x.
%      Can also be specified as '()'.
%      x will be rounded to the order of dx and printed as 'x(sd)'.
%   pmuncertainty: dx is uncertainty of x.
%      Can also be specified as 'pm'.
%      x will be rounded to the order of dx and printed as 'x ± dx'.
%   significant: dx is number of significant digits of x.
%      Can also be specified as 's'.
%      x will be rounded to dx significant digits when printed.
%   decimals: dx are the rounding orders of x.
%      Can also be specified as 'd'.
%      x will be rounded to dx digits in relation to decimal point
%      when printed.
%
%   label (optional) [character vector | string]
%   Label of x. Specify as [] for no label
%
%   units (optional) [character vector | string]
%   Units of x. Specify as [] for no units
%
%   p [0 (default) | integer scalar | integer vector | 'auto' | 'auto-unique']
%   Exponent of scientific notation
%   auto: p is selected automatically and is the same for every x
%   auto-unique: p is selected automatically and is unique for every x
%
%If x are scalars, prints out 'label = x units'.
%If x are vectors, prints out header string 'label, units'
%followed by multiple lines of 'x'.
%x is rounded according to specified dx and output mode.
%If mode is specified as 'pmuncertainty', ' ± dx' is added to every x
%when printing.
%If type is specified as '(uncertainty)', '(sd)' is added to every x
%when printing, where sd is 1 or 2 significant digits of dx.
%
%If dx is specified as scalar while x is vector, dx is assumed to be
%the same for every x.
%
%If dx is uncertainty, it is rounded to:
%- 1 significant digit if first two significant digits are 15 or more;
%- 2 significant digits if first two significant digits are 14 or less;
%and x is rounded to the same order as dx.
%
%If p isn't zero, scientific notation is used.

% Initialize unspecified arguments with default values
if nargin < 6
    p = 0;
    if nargin < 5
        units = [];
        if nargin < 4
            label = [];
            if nargin < 3
                mode = '()';
                if nargin < 2
                    dx = 0;
                end
            end
        end
    end
end

% Check input arguments
if any(strcmpi(mode,["(uncertainty)","()"]))
    mode = '()';
elseif any(strcmpi(mode,["pmuncertainty","pm"]))
    mode = 'pm';
elseif any(strcmpi(mode,["significant","s"]))
    mode = 's';
elseif any(strcmpi(mode,["decimals","d"]))
    mode = 'd';
else
    error(['[type] must be either ',...
        '''(uncertainty)'', ''pmuncertainty'', ''significant'', or ''decimals''']);
end
if isempty(x)
    return;
elseif ~isvector(x)
    error('[x] must be scalar or vector');
elseif ~(isscalar(dx) || all(size(dx) == size(x)))
    error('[dx] must be scalar or vector of the same length as x');
elseif any(strcmpi(mode,["pm","()"])) && any(dx <= 0)
    error('Uncertainty must be positive');
elseif ~(isvector(p) && all(p == round(p)) || any(strcmpi(p,["auto","auto-unique"])))
    error('[p] must be an integer scalar, ''auto'' or ''auto-unique''');
end
[label,units] = convertStringsToChars(label,units);

% Transform dx into a vector in case if it's a scalar
dx = dx.*ones(size(x));
% Choose the common exponent for all x if p is 'auto' or 'auto-unique'
if strcmpi(p,'auto')
    p = round(median(floor(log10(abs(x)))));
    if abs(p) < 3
        p = 0;
    end
elseif strcmpi(p,'auto-unique')
    p = floor(log10(abs(x)));
    p(abs(p) < 3) = 0;
end
% Transform p into a vector in case if it's a scalar
if isscalar(p) && ~isscalar(x)
    is_p_unique = false;
    p = p.*ones(size(x));
else
    is_p_unique = true;
end
% Set the rounding order of each element of x
if any(strcmpi(mode,["pm","()"]))
    r = -floor(log10(dx)) + (dx.*10.^-floor(log10(dx)) < 1.5) + p;
elseif strcmpi(mode,'s')
    r = dx - 1 - floor(log10(abs(x))) + p;
elseif strcmpi(mode,'d')
    r = dx + p;
end

% Format output to Command Window
fprintf('\n');
if isscalar(x)
    if ischar(label)
        fprintf('%s = ',label);
    end
else
    if ischar(label)
        fprintf('%s',label);
        if ischar(units) || ~is_p_unique && p(1) ~= 0
            fprintf(',');
        end
    end
    if ~is_p_unique && p(1) ~= 0
        fprintf(' ×10^%d',p(1));
    end
    if ischar(units)
        fprintf(' %s',units);
    end
    if ischar(label) || ischar(units) || ~is_p_unique && p(1) ~= 0
        fprintf('\n');
    end
end
for i = 1:numel(x)
    if strcmpi(mode,'pm') && is_p_unique && p(i) ~= 0
        fprintf('(');
    end
    if r(i) > 0
        fprintf('%.*f', r(i), x(i)*10^-p(i));
        if strcmpi(mode,'pm')
            fprintf(' ± %.*f', r(i), dx(i)*10^-p(i));
        elseif strcmpi(mode,'()')
            fprintf('(%.0f)', dx(i)*10^(-p(i)+r(i)));
        end
    else
        fprintf('%.f', round(x(i)*10^-p(i),r(i)));
        if strcmpi(mode,'pm')
            fprintf(' ± %.f', round(dx(i)*10^-p(i),r(i)));
        elseif strcmpi(mode,'()')
            fprintf('(%.0f)', round(dx(i)*10^-p(i),r(i)));
        end
    end
    if is_p_unique && p(i) ~= 0
        if strcmpi(mode,'pm')
            fprintf(')');
        end
        fprintf('⋅10^%d',p(i));
    end
    if i ~= numel(x)
        fprintf('\n');
    end
end
if isscalar(x)
    if ischar(units)
        fprintf(' %s',units);
    end
end
fprintf('\n');