function fig1e_plot(datafilepath)
% function fig1e_plot(datafilepath)

load(datafilepath,'data');
convdat=data;
%convdata=timexpositionsxunits
%lineprops={'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',0.1};
figure;
stimpos=1:4;
tvec=[-1000:1:size(convdat,1)-1000-1];
col=cell2mat(p_colors('yellow_red')');
%loc_ha=axes('position',[0.27 0.53 0.16 0.12]);
set(gca,'Color','none');
s_on = [0 400 800 1200];
p_on = 3900;

for i=1:numel(stimpos)
    locdat=squeeze(convdat(:,i,:));
    data=mean(locdat,2);
    lpos(i)=plot(tvec,data,'color',col(i,:),'linewidth',1);
    upper=std(locdat,[],2)./sqrt(size(locdat,2))+data;
    lower=data-std(locdat,[],2)./sqrt(size(locdat,2));
    hold on;
    plot(tvec,upper','color',col(i,:),'linewidth',0.2,'Linestyle','-');
    hold on;
    plot(tvec,lower','color',col(i,:),'linewidth',0.2,'Linestyle','-');
    hold on;
end
hs=vline(s_on,'k-');
hs(1).Color=[0.7 0.7 0.7];
hs(2).Color=[0.7 0.7 0.7];
hs(3).Color=[0.7 0.7 0.7];
hs(4).Color=[0.7 0.7 0.7];
hp=vline(p_on,'k-');
hp.Color=[0.7 0.7 0.7];
%legend(lpos,{'P1';'P2';'P3';'P4'}, 'location','best');
%legend boxoff;
% xlabel('ms from stim onset');
% only show first 500ms (instead of 1000ms) pre stim
curxlim=get(gca,'Xlim');
curylim=get(gca,'Ylim');

xlim([-500 tvec(end)]);
ylim([-.3 7]);
set(gca,'Xlim',([curxlim(1)+500   tvec(end)]));
set(gca,'Xtick',[0:400:4000]);
set(gca,'Xticklabels',{'0','0.4','0.8','1.2','4'});
set(gca,'Fontsize',12);
xlabel('time from 1st stimulus onset (s)');
ylabel('Zscored Spikerate');
set(gca,'Color','none');
set(gca, 'Box','off');
set(gca,'YColor','k');
applyfigprops(gcf);
applyaxprops(gca);


