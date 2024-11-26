function [h]=lineatzero(a,xval, varargin)
%function [h]=lineatzero(a,xval)
%print straight line starting at x=0 from minimal to maximum yvalue
% unless x is specified in xval, a is axis handle
if nargin<1
    a=gca;
    xval=0;
    col='r';
   
end

if nargin<=2
    xval=0;
    col=[0 0 0];
end

if nargin>2
    col=varargin{1};
end
% multiple axes k

[y]=ylim(a);
maxi=max(y);
mini=min(y);

axes(a);
h=line([xval xval],[mini maxi],'color',col);

