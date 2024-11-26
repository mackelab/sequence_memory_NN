function fig2c_plot_gh(datafilepath)
%function fig2c_plot(datafilepath)

load(datafilepath);
avgtimeparams=loadparams_avgtimes;

% plot of average Zscore during DELAY sorted by site region als 2d plot
TFmeandelay=(squeeze(mean(TFdata(:,avgtimeparams.twindelay(1):avgtimeparams.twindelay(2),:),2)))';
orderreg=[2 4 3 1];
figure;
set(gcf,'position',[726   580   248   606]);
ax=gca;
set(ax,'position',[0.1 0.05 0.75 0.8]);

% rearrange data so that it appears in order of line plots
indices=[];
cats=[];
for k=1:numel(orderreg)
    locindex=index(sortcat==orderreg(k));
    cats=[cats;repmat(orderreg(k),numel(locindex),1)];
    indices=[indices;locindex];
end

imagesc(anaparams.freqs,[],flipud(TFmeandelay(indices,:))),...
    ,colormap jet,axis tight;hold on;set(gca,'Ydir','normal');

set(gca,'xtick',2:10:100);
set(gca,'Color','none');
hold all;
%linehandle=hline(find(diff(sortcat)),'k');
linehandle=hline(find(diff(flipud(cats))),'k');

% make xticklabels
set(gca,'Xticklabel',{'',num2str(round(anaparams.freqs(12),1)),'','','',...
    num2str(round(anaparams.freqs(52),1)),'','','',num2str(round(anaparams.freqs(100),1))});
box off;
set(gca,'yticklabel','');
ylabel('channel count');
xlabel('frequency (Hz)');
axhcb=colorbar('position',[0.9 0.50 0.02 0.25]);
