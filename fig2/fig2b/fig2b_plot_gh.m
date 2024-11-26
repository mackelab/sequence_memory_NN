function fig2b_plot_gh(datafilepath)
%function fig2b_plot(datafilepath)

load(datafilepath);
avgtimeparams=loadparams_avgtimes;
nunits=size(sortcat,1);

% new order for plotting hippocampus (2) Enthorinal Cortex (4) {Parahippocampal (3) and
% Amygdala (1)
orderreg=[2 4 3 1];
figure;
nfreqs=numel(anaparams.freqs);

% make double axis in same spot to plot on top of eachother
ah(1)=axes('Position',[0.100    0.7800    0.1800    0.2200]);
ah(2)=axes('Position',[0.100    0.5300    0.1800    0.2200]);
ah(3)=axes('Position',[0.100    0.2800    0.1800    0.2200]);
ah(4)=axes('Position',[0.100    0.030    0.1800    0.2200]);
ahs(1)=axes('Position',[0.100    0.7800    0.1800    0.2200]);
ahs(2)=axes('Position',[0.100    0.5300    0.1800    0.2200]);
ahs(3)=axes('Position',[0.100    0.2800    0.1800    0.2200]);
ahs(4)=axes('Position',[0.100    0.030    0.1800    0.2200]);
set(ah(:),'Color','none');
set(ahs(:),'Color','none');
set(ahs(:),'XColor','none');
set(ahs(:),'Yaxislocation','right');
set(gcf,'position',[109   126   991   781]);

% define ylims for single data per region
ylims_single(1,:)=[-4 4];
ylims_single(2,:)=[-2.5 2.5];
ylims_single(3,:)=[-2.5 2.5];
ylims_single(4,:)=[-2.2 2.2];

% define ylims for mean data per region
ylims_mean(1,:)=[-0.7 0.7];
ylims_mean(2,:)=[-2.5 2.5];
ylims_mean(3,:)=[-1.5 1.5];
ylims_mean(4,:)=[-.75 .75];

% line plots Z power change delay, also separate for all siteregions
% for different regions
region{1}='HPC';
region{2}='PHC';
region{3}='EC';
region{4}='AM';

for k=1:numel(orderreg)
    temp=TFdata(:,avgtimeparams.twindelay(1):avgtimeparams.twindelay(2),index(sortcat==orderreg(k)));
    
    axes(ah(k));
    plot(anaparams.freqs, squeeze((mean((temp),2))),'color',[0.7 0.7 0.7]);
    xlim([-5 66]);
    ylim(ylims_single(k,:));
    
    % plot mean on top using different axis
    hold on;
    axes(ahs(k));
    set(ahs(k),'Color','none');
    set(ahs(k),'XColor','none');
    set(ahs(k),'Yaxislocation','right');
    plot(anaparams.freqs,  mean(squeeze((mean((temp),2))),2),'color','k');
    xlim([-5 66]);
    ylim(ylims_mean(k,:));
    ylabel(region{k});
    set(ahs(k),'Xtick',[]);
    set(ahs(k),'Xticklabel',[]);
    set(ahs(k),'Yticklabel',[]);
    set(ah(k),'Xtick',anaparams.freqs(2:4:end));
    set(ah(k),'Xticklabel',[]);
    set(ahs(k),'Color','none');
    set(ah(k),'Color','none');
    set(ahs(1),'Ylim',ylims_mean(k,:));
end

applyfigprops(gcf);

