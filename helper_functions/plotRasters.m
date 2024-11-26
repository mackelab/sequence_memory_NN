function [axhandle,linehandles]=plotRasters(data,varargin)
% Plot spike train rasters. input: data(trials,timebins)
%
LW=2; % linewidth for spikes as vertical lines
data=data';
% parse and remove this function's input params from varargin
%params.x = [0 1000];
params.density = 1.5;
%params.horColor = [.9 .9 .9];
params.horColor = [];
params.linewidth=5;
params.x=[];

f = fieldnames(params);

for i = 1:length(f)
    ndx = strmatch(f{i},{varargin{1:2:end}});
    if ~isempty(ndx)
        params.(varargin{ndx}) = varargin{ndx+1};
        varargin(ndx:ndx+1) = [];
    end
end

% find spikes and map indices to times in given interval
[x,y] = find(data);
x = (x-1) / size(data,1) * diff(params.x) + params.x(1);

% plot horizontal lines?
n = size(data,2);
if ~isempty(params.horColor)
    plot(repmat(params.x',1,n),repmat((1:n),2,1),'Color',params.horColor,'linewidth',1);
end

% plot spikes as short vertical lines
Y = [y-params.density/2, y + params.density/2]';
X = [x, x]';
hold on

linehandles=plot(X,Y,'k','LineWidth',LW,varargin{:});



set(gca,'XTick',params.x,'XLim',params.x, ...
        'YTick',[1 n],'YTicklabel',[1,n],'YLim',[0 n+1+params.density])
    box off
    ylim([params.density/2,n+1-params.density/2])
    axhandle=gca;
    
    
