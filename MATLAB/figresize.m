function figresize(dim,length)
%Resize the current figure to the specified width or height
%while preserving aspect ratio
%
%   dim ['width' | 'height']
%   Dimension that will be adjusted to specified length
%
%   length [scalar]
%   Final length of specified dimension (including magins for text labels),
%   in figure's units

% Get current figure and screen properties
position = get(gcf,'Position');
aspect_ratio = get(gca,'PlotBoxAspectRatio');
aspect_ratio = aspect_ratio(1)/aspect_ratio(2);
set(gca,'PlotBoxAspectRatioMode','auto');
inset = get(gca,'TightInset');
width_inset = (inset(1) + inset(3))*position(3);
height_inset = (inset(2) + inset(4))*position(4);
set(groot,'units',get(gcf,'Units'));
screen_size = get(groot,'ScreenSize');

% Calculate new figure width and height
if strcmpi(dim,'width')
    width = length;
    height = (width - width_inset)/aspect_ratio + height_inset;
elseif strcmpi(dim,'height')
    height = length;
    width = (height - height_inset)*aspect_ratio + width_inset;
else
    error('First input argument must be ''width'' or ''height''');
end

% Set new figure size and position
position = [...
    (screen_size(3) - width)/2 ...
    screen_size(4) - height - 3 ...
    width ...
    height];
set(gcf,'Position',position);