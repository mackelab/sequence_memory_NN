function [fhandle]=applyfigprops(fhandle)


% properties for axis
props1={'Color',[1 1 1],'paperorientation','landscape'};

set(fhandle(:),props1{:});

