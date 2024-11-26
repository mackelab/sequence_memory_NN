function fig3_plot_reb(datafolder,datafiles)
% provide foldername for data files and cellarray containing filenames
% function requires circ_stats toolbox

plotwf=1; % option to plot waveform or not (default = 1)
Nunits=numel(datafiles);
posx_hist=linspace(0.03,0.9,Nunits);
freqs=[1.5000 1.7466 2.0336 2.3679 2.7571 3.2103 3.7380 4.3524 5.0678 5.9008 6.8707 8.0000];
% size bar spk histograms
sizeplot=[0.066 0.07];

fh=figure;
applyfigprops(fh);
set(gcf,'position',[379        -244        1891         872]);


unitfreqs=[];

for u=1:Nunits

    load(strcat(datafolder,datafiles{u}));
    
    
    axpos=[posx_hist(u) 0.71 sizeplot;...
        posx_hist(u) 0.63 sizeplot;...
        posx_hist(u) 0.56 sizeplot;...
        posx_hist(u) 0.48 sizeplot];
    
    [axh]=plothistograms(spikes, phases, sortinfo, axpos);
    
    % waveform plots
    if plotwf
    waveformpos=[posx_hist(u) 0.2 0.1 0.15];
    [axhwf]=plotwaveforms(wvf,idx,waveformpos);
    end
    unitfreqs(u,1)=freqs(freqindex);
    unitfreqs(u,2)=(freqindex);
    
    
       th=annotation('textbox','EdgeColor','none','position',[posx_hist(u) 0.05 0.15 0.05],...
        'String',{[num2str(u)]});
end

unitfreqs

function [axhwf]=plotwaveforms(wvf,idx,waveformpos)

axhwf=axes('position',waveformpos);
density_plot(wvf, idx);
%ylabel('microV');
%xlabel('mS');

            
function [axh]=plothistograms(spikes, phases, sortinfo, axpos)
% plot phase triggered spiking histograms

oneylim=1;
col_s=p_colors('yellow_red');
normalizedata=0;
pos=(unique(sortinfo(:,3)));
spkmean=[];
meanphasespike=[];
stdmeanphasespike=[];
phaseatspikes_pos=cell(numel(pos),1);

% go through positions
for j=1:numel(pos)
    tempmean=[];
    ix=sortinfo(:,3)==pos(j);
    
    
    temp=phases(ix,:);
    spk=spikes(ix,:);
    temp=temp(:);
    spk=spk(:);
    
    
    %% calculate histogram counts for phase vectors
    % length of window
    winsize=size(spikes,2);
    phaseatspike=temp(logical(spk));
    [binvector,NN,CN]=computephasehistogram(phaseatspike,1);
    % sort all phases into phasebins
    [N,bins]=histc(temp,binvector);
    % get rid of last (0) bin because always 0
    N=N(1:end-1);
    
    
    % number of spikes / number of phases *1000 is instantaneous firing
    % rate in Hz normalized by number of phases
    xx=(CN./N)*1000;
    
    % prop spikes per bin
    spkmean(j,:)=xx./sum(xx); 
    
    % take mean of phases@spike
    meanphasespike(j)=circ_mean(phaseatspike);
    stdmeanphasespike(j)=circ_std(phaseatspike);
    
    % collect spikes at phase
    phaseatspikes_pos{j}=phaseatspike;
    
    
end % end positions



if normalizedata
    % normalization by max across all conditions , preserving spike reate
    % differences between conditions
    spkmean=spkmean./(max(spkmean(:)));
    
    % normalization by mean per condition, not preserving spikerate differences
    [y,i]=max(spkmean,[],2)
    m=repmat(y,1,size(spkmean,2))
    spkmean=spkmean./m;
    
end


pistep=pi/(pi/abs(circ_dist(binvector(1),binvector(2))));
xval=[binvector(1:end-1) pi:pistep:3*pi];
xval=xval(1:end-1);
strpos={'pos1';'pos2';'pos3';'pos4'};
ylimmax=(max(spkmean,[],2)+0.2.*(max(spkmean,[],2)));
%%
for s=1:numel(pos)
    axh(s)=axes('position',axpos(s,:));
    hb=bar(xval,[spkmean(s,:) spkmean(s,:)],'Facecolor',col_s{s},'Edgecolor',col_s{s},'Barwidth',1.2);
    
    
    axis tight;
    set(gca,'Xtick',[]);
    yl=get(gca,'Ylim');
    xl=get(gca,'Xlim');
    hold on;
    
    
    % add horizontal line barwidth wide
    % where mean direction is
    if oneylim
        x12=[meanphasespike(s)-stdmeanphasespike(s) meanphasespike(s)+stdmeanphasespike(s)];
        y12=[1.15*yl(2) 1.15*yl(2)];
        plot(x12,y12,'-','Color', col_s{s});
        plot(meanphasespike(s),1.15*yl(2),'o','MarkerSize',5,'Markeredgecolor',col_s{s},'Markerfacecolor',col_s{s});
    else
        x12=[meanphasespike(s)-stdmeanphasespike(s) meanphasespike(s)+stdmeanphasespike(s)];
        y12=[0.95*ylimmax(s) 0.95*ylimmax(s)];
        plot(x12,y12,'-','Color', col_s{s});
        plot(meanphasespike(s),0.95*ylimmax(s),'o','MarkerSize',5,'Markeredgecolor',col_s{s},'Markerfacecolor',col_s{s});
    end
    
    axis tight;
    
end

set(gca,'Xtick',[linspace(-pi,pi,5)],'XtickLabel',{'-pi','-pi/2','0','pi/2','pi'});

if oneylim % one ylim axis for all
    set(axh, 'Ylim',[0 max(spkmean(:))+0.2*max(spkmean(:))]); % set equal axes for all
    mini=min(spkmean(:))-0.15*min(spkmean(:));
    if mini>0
        set(axh, 'Ylim',[mini max(spkmean(:))+0.15*max(spkmean(:))]);
    end
    
else
    
    for Y=1:size(ylimmax,1)
        set(axh(Y), 'Ylim',[0 ylimmax(Y)]); % set individual axes
        mini=min(spkmean(:))-0.15*min(spkmean(:));
        if mini>0
            set(axh(Y), 'Ylim',[mini ylimmax(Y)]);
        end
        
    end
end

set(axh, 'Xlim',[xl(1) 2*pi]);


function phase_color = p_colors(palette)
% function phase_color = p_colors(palette)
switch palette
    case 'yellow_red'
        
        phase_color{1}=[255,201,70]./255;
        
        phase_color{2}=[253,141,33]./255;
        %phase_color = 'r';
        
        phase_color{3}=[227,26,28]./255;
        %phase_color = 'g';
        
        phase_color{4}=[142,23,15]./255;
        %phase_color = 'k';
        
    case 'green_blue'
        phase_color={[161 218 180]./255; [65 182 196]./255; [34 94 168]./255; [10 30 69]./255};
    case 'uni_blue'
        phase_color={[239,243,255]./255; [189,215,231]./255; [107,174,214]./255; [33,113,181]./255};
        % HSL 225, 100%, 97
    case 'myblue'
        phase_color=[107,174,214]./255;
    case 'mygray'
        phase_color=[0.3 0.3 0.3];
end


function applyfigprops(fhandle)
props1={'Color',[1 1 1],'paperorientation','landscape'};
set(fhandle(:),props1{:});


function [binvector,NN,varargout]=computephasehistogram(phaseatspike,sm,varargin)
% function [binvector,NN]=computephasehistogram(phaseatspike,sm)
% computes phase histogram of spikes  at phases for inputvector of phases @
% spikes

if nargin<2
    sm=1;
end

if nargin==2
    varargout=cell(1,1);
end

% binning vector for spiking data
% each bin will be pi/binsize wide
binsize=17; % (with 17 bins difference between bins is pi/8)
binvector=linspace(-pi,pi,binsize);
[NN,BINS]=histc(phaseatspike,binvector);
NN=NN(1:end-1);

% for plotting only, does not affect stats
if sm
    % create box filter
    % window size be approximately 1/10 of original binning
    % filterlength
    filtwidth=round(binsize/4);
    u=[zeros(1, round(filtwidth/2)) ones(1, filtwidth) zeros(1, round(filtwidth/2))];
    
    % convolve 2 cycles of data
    N2=[NN;NN];
    CN=conv(N2,u,'same');
    CN=CN(1:binsize-1);
    % transform counts back to original units of NN
    CN=CN./sum(u);
    % normalize by count to get probabilities
    varargout{1}=CN;
end

function density_plot(spikes, idx)
%%function density_plot(spikes, idx)

color_ = [0 0 0];
if ~exist('idx', 'var')
    wavs = spikes;
else
    wavs = spikes(idx,:);
end
mean_spike = mean(wavs);
std_spike = std(wavs);
lbound=floor(min(min(wavs)));  %Lowest point of the figure
ubound=ceil(max(max(wavs)));  %Highest point of the figure
vps=size(spikes,2);   %Values per spike after interpolation. 64 without interpolation

%Make the 2D histogram
ybins=linspace(lbound,ubound,150);  %Vector of bins in vertical direction
ybinSize=ybins(2)-ybins(1);         %Size of a bin in vertical direction
numXbins=vps;                       %Number of bins in horizontal direction
numYbins=length(ybins);             %Number of bins in vertical direction
n=zeros(numXbins,numYbins);         %Preallocate 2D histogram

%Bin count
for k=1:numXbins
    for j=1:numYbins
        n(k,j)= sum ( wavs(:,k) <= ybins(j)+ybinSize/2 & wavs(:,k) > ybins(j)-ybinSize/2);
    end
end

%Creates colormaps for each cluster
maxN=max(max(n));   %determine amount of colorsteps in colormap (color resolution)
colorMap=zeros(maxN+1,3);   %preallocate space for the colormap. 3-column RGB vector
colorMap(1,:)=[1 1 1];   %first value of colormap set to white. 
cBuffer=zeros(1,maxN);
for e=1:3
    cBuffer(1:ceil(maxN/2))=color_(e)-linspace(0,1,ceil(maxN/2));  % in equal steps, increment from color of choice (myColors{i}) towards black
    cBuffer(cBuffer < 0)=0;     % RGB are values between 0 and 1. can't go below black (=0).
    cBuffer((ceil(maxN/2)+1):maxN)=linspace(0,1,maxN-ceil(maxN/2));  % in equal steps, increment from color of choice (myColors{i}) towards black
    cBuffer(cBuffer > 1)=1;     % RGB are values between 0 and 1. can't go above white
    colorMap(2:end,e)=cBuffer;   % assign cBuffer colormap to i in colorMaps
end

% remove extreme outliers in order to keep color resolution
cutoff = 5*std(reshape(n,numel(n),1)); %magic cutoff for too high bin values
n(n>cutoff) = cutoff; %replace n with n without too high bin values
% pcolor(n'),shading interp % plot the 2D histogram
xvals = [-19:44]*1000/(2^15);
try
    pcolor(xvals, ybins, n');
catch
    pcolor(n')
end

shading interp
% xlim([0 numXbins]);
% ylim([0 numYbins]);
xlabel('ms'); 
ylabel('\muV');