function fig1b_plot(datafilepath)
% code for plotting Figure 1b
% for plotting you will need to have gramm installed for Matlab
% https://github.com/piermorel/gramm
% load performance and RT data

%         disp(['Mean performance is ' num2str(mean(pc))])
%         disp(['STD performance is ' num2str(std(pc))])
%
%         disp(['Mean RT is ' num2str(mean(mrt(:,1)))]);
%         disp(['STD RT is ' num2str(std(mrt(:,1)))]);

load(datafilepath,'pc','mrt','patients','patindex');
perf=pc';
RT=mrt(:,1);

% compare mean performance to chance
% 224 is number of trials, 0.25 is chance performance
pval=binopdf(mean(perf),224,0.25);

% correlation between perf and RT
[rho, pvalrho]=corr(perf,RT,'type','Spearman');

% RT for correct and incorrect trials
RT_p=[mrt(:,2);mrt(:,3)]; 
[nsubs]=size(mrt,1);
group1=[ones(nsubs,1);ones(nsubs,1)*2];

% define colors
green=[161 218 180]./255;
blue=[33,113,181]./255;
gray=[0.5 0.5 0.5];
gray0=[0.3 0.3 0.3];
red=[181,11,33]./255;


figure;
set(gcf,'position',[702        1150         551         188]);

% plot median RTs for correct vs incorrect trials

g(1,1)=gramm('y',RT_p./1000,'x',group1,'color',group1);
% boxplot
g(1,1).stat_boxplot('width',0.7,'notch',true,'dodge',0.7);
g(1,1).geom_point();
%g(1).axe_property('YLim',[-0.02 0.4],'XLim',[-0.02 0.4]);
%g.set_continuous_color('LCH_colormap',[20 80 ; 40 30 ; 260 260 ]);
%g(1,1).set_color_options('chroma',-20,'lightness',50);
% g(1,1).set_color_options('lightness',50,'chroma',[90] )%, 'hue_range',[25 385]); % grey 'chroma',285,
g(1,1).set_color_options('map',[blue;red],'n_color',2,'n_lightness',1)

% position axis
% g(1,1).set_layout_options('Position',[0.38 0.15 0.22 0.22]);%,...
g(1,1).no_legend();
g(1,1).axe_property('Xlim',[-.5 4.5],'Ylim',[2.000 9.000]);
g(1,1).set_names('y','RT (s)','x',{'performance'});
%g(1,1).set_color_options('hue_range',[-60 60],'chroma',40,'lightness',90);

% plot pc's across subjects

g(1,2)=gramm('y',perf,'x',[ones(size(perf,1),1)],'color',patindex);

g(1,2).geom_point();
g(1,2).set_color_options('hue_range',[-60 60],'chroma',40,'lightness',90);
% make all dots different colors for different subjects
%g(1,2).set_color_options('map','brewer1');
g(1,2).axe_property('Xlim',[0.5 3],'Ylim',[0 1.2]);
g(1,2).no_legend();
g(1,2).set_names('y', 'percentage correct','x',{'performance'});

% scatter plot median RTs vs pc correct
g(1,3)=gramm('x',RT./1000,'y',perf);
g(1,3).set_color_options('chroma',0,'lightness',30);

g(1,3).stat_glm('geom','area','disp_fit',false);
%g.draw();
%snapnow;
g(1,3).geom_point();

g(1,3).set_color_options('map',[gray],'n_color',1,'n_lightness',1)
g(1,3).set_names('y', 'percentage correct','x','RT (s)');
g(1,3).axe_property('Ylim',[0.2 1.2],'Xlim',[1.500 8.000]);
g.draw();

% edit properties
%colormap(g(2).facet_axes_handles,'jet');
colormap(g(2).facet_axes_handles,[0.5 0.5 0.5]);
colormap(g(1).facet_axes_handles,[0.3 0.3 0.3;0.5 0.5 0.5]);
set(g(2).results.geom_point_handle,'SizeData',30);
g(1).results.stat_boxplot(1).box_handle.FaceAlpha=0.3;
g(1).results.stat_boxplot(2).box_handle.FaceAlpha=0.3;



set(g(1).facet_axes_handles,'XtickLabels',{});

set(g(2).facet_axes_handles,'XtickLabels',{});

set(g(1).facet_axes_handles,'Ytick',[3.000 4.000 5.000 6 7 8]);
set(g(3).facet_axes_handles,'Xtick',[3.000 4.000 5.000 6]);
set(g(2).facet_axes_handles,'Ytick',[0.4 0.6 0.8 1]);

set(g(3).facet_axes_handles,'Ytick',[0.4 0.6 0.8 1]);
set(g(2).facet_axes_handles,'Ytick',[0.4 0.6 0.8 1]);
set(g(2).facet_axes_handles,'YLim',[[0.2 1.2]]);
set(g(3).facet_axes_handles,'XLim',[[1.5 8]]);
% add median performance across subjects
axes(g(2).facet_axes_handles);
plot(1,median(perf),'ok','Markersize',11);



% add hlines for chance level
axes(g(2).facet_axes_handles)
hold on
hline(0.25);

axes(g(3).facet_axes_handles)
hold on
hline(0.25);