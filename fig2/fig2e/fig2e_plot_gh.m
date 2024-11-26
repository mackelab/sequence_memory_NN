function fig2e_plot_gh(datafilepath)
% this function requires the circ_stats toolbox and the gramm toolbox for
% matlab
% plots figures of fig2e and does stats
% requires gramm plotting toolbox

load(datafilepath);
%%
numbins=numel(params.bincenters(1:end-1));
nfreqs=numel(params.freqs);
nunits=217;
xax=linspace(-pi,pi,numbins)-1/2*pi;
binwidth=xax(1)-xax(2);

% part 1
% individual fits to get units kappa estimates, compare kappa between
% between preferred stimulus (stimulus in encoded sequence) vs.
% not
locdatps=squeeze(mean(data{1},1));
locdatnps=squeeze(mean(data{2},1));

thps=[];
kapps=[];
thnps=[];
kapnps=[];

for u=1:nunits
    [thps(u),kapps(u)]=circ_vmpar(xax,locdatps(:,u),binwidth);
    [thnps(u),kapnps(u)]=circ_vmpar(xax,locdatnps(:,u),binwidth);
end

kappa_units_ps_nps=[kapps' kapnps'];
%run stats
%stats=runstatsonkappa(kappa_units_ps_nps,Rates);


%% plotting - barplot
clear g
figure;
group=repmat([1 2],size(kappa_units_ps_nps,1),1);
group=group(:);
g(1,1)=gramm('x',group,'y', kappa_units_ps_nps(:),'color',group);
g(1).stat_boxplot('dodge',0.8,'width',1.1);
g(1,1).no_legend();
g(1,1).set_names('x',{''},'y','kappa', 'color', {''},'lightness',{''});
g.draw();
ax=gca;
set(ax,'Xticklabel',{'pS';'npS'},'Xtick',[1 2]);

g(1).results.stat_boxplot(1).box_handle.FaceAlpha=0.5;
g(1).results.stat_boxplot(1).outliers_handle.Marker='none';
g(1).results.stat_boxplot(2).outliers_handle.Marker='none';
g(1).facet_axes_handles.YLim=[-0.2 0.8];
g(1).results(1).stat_boxplot(1).box_handle.FaceColor=[0 0.5 1];
g(1).results(1).stat_boxplot(2).box_handle.FaceColor=[0.5 0.5 0.5];
g(1).set_layout_options('position',[0.05    0.1000    0.20    0.80]);
g(1).axe_property('YLim',[-0.2 0.8],'XLim',[-0.5 2.5],'YTick',[0:0.2:0.8],'XTick',[1 2]);
g(1).facet_axes_handles.XTick=[1 2];
g(1).facet_axes_handles.XTickLabel={'pS';'npS'};

% lineplot % second panel of figures: line plot with modulation of rate as function of
% phase, use bootstrapping to get CIs
% all
datup=[];
datupn=[];

% concat across all theta freqs to produce line plots
for t=1:nfreqs
    datup=cat(2,datup,squeeze(data{1}(t,:,:)));
    datupn=cat(2,datupn,squeeze(data{2}(t,:,:)));
end
YLIM=[min([datup(:);datupn(:)])-(0.001)*max([datup(:);datupn(:)]) max([datup(:);datupn(:)])+(0.001)*max([datup(:);datupn(:)])];
group=repmat([1 2],size(datup,2),1);
group=group(:);
g(1,2)=gramm('x',1:size(datup,1)*2,'y', [[datup;datup] [datupn ;datupn]]','color',group);
g(1,2).stat_summary('type','bootci');
g(1,2).no_legend();
g(1,2).set_layout_options('position',[0.3    0.100    0.3500    0.6000]);
g(1,2).set_color_options('n_color',2,'map',[0 0.5 1;0.5 0.5 0.5]);
g(1,2).set_names('x',{''},'y','prop. per bin', 'color', {'stimulus'},'lightness',{''});
g.draw();
g(1,2).facet_axes_handles.YLim=[1.0181e-04 0.00025]
g(2).facet_axes_handles.XTick=numbins;
g(2).facet_axes_handles.XTickLabel={'pref. phase'};
%g(1).set_title({'Kappa and phase modulation (centered at pref. phase) during delay'; 'for pS in encoded sequence or not (npS)'});
set(gcf,'position',[395   603   788   569]);

% collect stats on confidence intervals
stats.yci.pS=g(2).results.stat_summary(1).yci;
stats.yci.npS=g(2).results.stat_summary(2).yci;


function stats=runstatsonkappa(kappa_units_ps_nps,Rates)
% statistics, compare median kappa values , all units
[pkappa_allunitsPSNPS,h,tempstats]=signrank(kappa_units_ps_nps(:,1),kappa_units_ps_nps(:,2));
Zkappa_allunitsPSNPS=tempstats.zval;

stats.pkappa_allunitsPSNPS=pkappa_allunitsPSNPS;
stats.Zkappa_allunitsPSNPS=Zkappa_allunitsPSNPS;

% now controlling for rate differences!
rateDelay=Rates{2};
rateDelayNPS=Rates{3};
% first test difference in rates between baseline and delay
% pref stim in sequence (Rates{2}) vs not (Rates{3}
[stats.pPSNPSDelay,h,tempstats]=signrank(rateDelay,rateDelayNPS);
stats.ZPSNPSDelay=tempstats.zval;

% CONTROL 1 lower ptile group where there is no differnce in spiking
diffRateLP=rateDelay-rateDelayNPS;
[sortdiff,sortix]=sort(diffRateLP);
p50=prctile(diffRateLP,50);
index_l=diffRateLP<p50;
index_h=diffRateLP>p50;
[stats.pkappa_allunitsPSNPS_l,h,tempstats]=signrank(kappa_units_ps_nps(index_l,1),kappa_units_ps_nps(index_l,2));
stats.Zkappa_allunitsPSNPS_l=tempstats.zval;
N_l=sum(index_l);

% CONTROL 2 sequentially removing units
sortdiff=fliplr(sortdiff);
sortix=fliplr(sortix);
sortedDelayNPS=rateDelayNPS(sortix);
sortedDelay=rateDelay(sortix);
%
% remove units one by one until you have two groups with no sign
% differences in spkrate 
for kk=1:size(diffRateLP,2)
    locp=signrank(sortedDelayNPS,sortedDelay);
    % if sign difference still then remove unit from rest
    if locp<0.05
        sortedDelayNPS=sortedDelayNPS(2:end);
        sortedDelay=sortedDelay(2:end);
        sortix=sortix(2:end);
    else
        break
    end
end

[stats.pkappa_allunitsPSNPS_rm,h,tempstats]=signrank(kappa_units_ps_nps(sortix,1),kappa_units_ps_nps(sortix,2));
stats.Zkappa_allunitsPSNPS_rm=tempstats.zval;
N_rm=numel(sortix);


