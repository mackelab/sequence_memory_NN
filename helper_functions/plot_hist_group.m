function handles=plot_hist_group(g1,bins,edgv,col,norm,varargin)
% handles=plot_hist_group(g1,bins,edgv,col,norm)
% input g1 = matrix w/ timestamps, edgv = binvector
% plots histograms of 1 group (y axis relative count)
% per group returns handles to bars, patches, axis and figure
% default bins =11; norm=1 (normalize to maximum elements per bin

if ~exist('norm','var')
    norm=1;
end


if nargin==1
    bins=11;
    edgv=linspace(min([g1]),max([g1]),bins);
    col=[0 0 0];
end

if nargin==2
    edgv=linspace(min([g1]),max([g1]),bins);
    col=[0 0 0];
end
if nargin==3
    
    col=[0 0 0];
end
if isempty(edgv)
    edgv=linspace(min([g1]),max([g1]),bins);
end
if nargin==6
    
    normHz=varargin{1}{1};
    ntrials=varargin{1}{2};
end

if ~exist('normHz','var')
    normHz=0;
    
end

if ~exist('ntrials','var')
    ntrials=[];
end

g1=iscolumnvector(g1);
[n1,bin]=histc(g1(:),edgv);
if norm
n1=n1./numel(g1);
end


if normHz
    binsize=abs(edgv(1)-edgv(2));
    n1=(n1/ntrials).*(1000/binsize);
end

ha(1)=bar(edgv,n1,1,'histc');

h=findobj(gca,'type','patch');
set(h,'Facecolor','none','Linewidth',1);
set(h(1),'EdgeColor',col);


handles.bars=ha;
handles.patches=h;
handles.axis=gca;
handles.fig=gcf;
handles.edgv=edgv;
handles.ydata=n1;
