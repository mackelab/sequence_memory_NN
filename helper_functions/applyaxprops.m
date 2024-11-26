function [axhandle,labhandle,leghandle]=applyaxprops(axhandle,labhandle,leghandle)

if nargin<1
    axhandle=gca;
end
% properties for axis
props1={'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'TickLength',[0.01 0.0],...
'Color',[1 1 1]};

set(axhandle,props1{:});

if nargin==2
    
    
    
% properties for labels and legend
props2={'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal'};
set(labhandle,props2{:});
end


if nargin==3
set(leghandle,props2{:});
end