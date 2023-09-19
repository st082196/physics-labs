function coeffvariances = coeffvariances(fitresult,level)
%Variances of coefficient values of cfit or sfit object
%
%   fitresult [cfit | sfit]
%   Fit object
%
%   level [0.95 (default) | scalar value in the range (0,1)]
%   Confidence level

if nargin == 1
    level = 0.95;
end
ci = confint(fitresult,level);
coeffvariances = (ci(2,:) - ci(1,:))./2;