function fig2a_plot_gh(datafilepath)

load(datafilepath,'data','anaparams','plotparams_spk');
%% create new axes
figure;
set(gcf,'position',[870   860   464   284]);
lineprops={'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',0.1};
%loc_ha=axes('position',[0.28 0.45 0.18 0.2]);
% tevents to plot, plotparams are loaded from datafile
nfreqs=numel(anaparams.freqs);

xtimes=repmat(plotparams_spk.m,2,1);
plotparams_spk.p=plotparams_spk.p-min(xtimes(:));
xtimes=xtimes-min(xtimes(:));
xtimes=repmat(plotparams_spk.s_on,2,1);
ytimes=repmat([1 nfreqs]',1, size(xtimes,2));

% sdat == pref stimulus
%imagesc(plotparams_spk.tvec+500,[],TF.sdat);% plot only 500ms baseline
imagesc(plotparams_spk.tvec+500,[],data);% plot only 500ms baseline
hold on;
line(xtimes,ytimes, lineprops{:});
line([plotparams_spk.p_on plotparams_spk.p_on],[1 nfreqs],lineprops{:});
%plotparams.p
set(gca,'Ydir','normal');
lh=lineatzero;set(lh,lineprops{:});
ch=colorbar;
xt=get(gca,'Ytick');
set(gca,'Yticklabel',num2str(round((anaparams.freqs(xt)'))));
ylabel('frequency');
xlabel('time from 1st stimulus onset in s');
%xlim([-500 4890+500]);
set(gca,'Xtick',[0 400 800 1200 plotparams_spk.p_on]);
set(gca,'Xticklabels',{'0','0.4','.8','.12','4'});
set(gca,'Box','off');
caxis([-8.2 8.2]);
colormap jet;

loc_ha=gca;
loc_ch=ch;

applyfigprops(gcf);
applyaxprops(gca);
