function fig4_plote_pop_gh(datafilepath)
% function fig4_plote_pop_gh(datafilepath)
% this function requires the gramm plotting toolbox from matlab

load(datafilepath);
fh(1)=figure;
axh=plot_popperf(svm,params)
stats=getstats(svm);
applyaxprops(axh);
plot_popperf_region(svm,params)

function axh=plot_popperf(svm,params)
phase_colors=cell2mat(p_colors('uni_blue'));

edgv=linspace(min([1-svm.results_dpu.performance(:);1-svm.results_dpl.performance(:)]),...
max([1-svm.results_dpu.performance(:);1-svm.results_dpl.performance(:)]),15);

[nperms,ngroups]=size(svm.results_dpl.performance);


c{1}='original';
c{2}='shuffled';
g1=[repmat(c(1),params.nperms,1);repmat(c(2),params.nperms,1)];
g2=[repmat(c(1),params.nperms*2,1);repmat(c(2),params.nperms*2,1)];
group1=[repmat(c(1),nperms*ngroups,1);repmat(c(2),ngroups*nperms,1)];

edgv2=linspace(min(1-svm.results_dpu.performance(:))-max(1-svm.results_dpl.performance(:)), ...
    max(1-svm.results_dpu.performance(:))-min(1-svm.results_dpl.performance(:)),15);

fh=figure(1);
set(fh,'position',[326   704   697   607]);
set(fh,'position',[326         158        1522        1153]);

gr(1,1)=gramm('x', 1-svm.results_dpu.performance(:),'color',g1);
gr(1,1).stat_bin('geom','overlaid_bar','normalization','probability','edges',edgv); %Overlaid bar automatically changes bar coloring to transparent
gr(1,1).set_title('dp_upper');
gr(1,1).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
gr(1,1).set_names('x','performance','y','probability', 'color', {''});
gr(1,1).set_layout_options('Position',[0.0662 0.5579 0.3 0.3]);


gr(1,2)=gramm('x', 1-svm.results_dpl.performance(:),'color',g1);
gr(1,2).stat_bin('geom','overlaid_bar','normalization','probability','edges',edgv); %Overlaid bar automatically changes bar coloring to transparent
gr(1,2).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
gr(1,2).set_names('x','performance','y','probability', 'color', {''});
gr(1,2).no_legend();
gr(1,2).set_layout_options('Position',[0.5662 0.5579 0.3 0.3])
gr(1,2).set_title('dp_lower');

gr(2,1)=gramm('x', g2,'y', [1-svm.results_dpu.performance(:);1-svm.results_dpl.performance(:)],'color',[g1;g1]);
gr(2,1).stat_boxplot('dodge',0.5,'width',0.7);
gr(2,1).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
gr(2,1).set_names('x','dp_lower/dp_upper','y','probability', 'color', {''});
gr(2,1).set_title('figure4e_leftpanel');
gr(2,1).set_layout_options('Position',[0.0662 0.0414 0.3 0.3]);

gr(2,2)=gramm('x',[1-svm.results_dpu.performance(:,1);1-svm.results_dpl.performance(:,1)] , 'y',[1-svm.results_dpu.performance(:,2);1-svm.results_dpl.performance(:,2)] ,'color',g1);
gr(2,2).geom_point('alpha',0.05);
gr(2,2).stat_cornerhist('edges',edgv2,'aspect',0.5);
gr(2,2).geom_abline();
gr(2,2).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
gr(2,2).set_layout_options('Position',[[0.5662 0.0579 0.3 0.3]]);
gr(2,2).axe_property('YLim',[0.15 0.5],'XLim',[0.15 0.5]);
gr(2,2).set_names('x','performance','y','performance', 'color', {''});

gr.draw();

gr(2,1).facet_axes_handles.XLim=[-1 4];
gr(1,1).facet_axes_handles.YLim=[ 0 0.35];
gr(1,2).facet_axes_handles.YLim=[ 0 0.35];




gr(2,1).results.stat_boxplot(1).box_handle.FaceAlpha=0.5;
%gr(2,1).results.stat_boxplot(1).outliers_handle.Marker='none';
gr(2,1).results.stat_boxplot(2).box_handle.FaceAlpha=0.5;
%gr(2,1).results.stat_boxplot(2).outliers_handle.Marker='none';
gr(2,1).axe_property('YLim',[-0.02 0.4],'XLim',[0.1 0.6]);


% add vline to histograms
axes(gr(1,1).facet_axes_handles);
vline(0.25,'k:');
axes(gr(1,2).facet_axes_handles);
vline(0.25,'k:');
axh=gca;



function plot_popperf_region(svm,params)


region_label={'hi','ec','am','phc'};
nregions=size(region_label,1);
[nperms,ngroups]=size(svm.results_allfreqs_hiboth.performance);
g2=[];
g1=[];
perf=[];
Zval=[];
locp=[];
% assemble data for plotting
for r=1:numel(region_label)
g2=[g2;ones(params.nperms*ngroups,1)*r]; % shuf vs non shuff label
perf_region=eval(strcat('svm.results_allfreqs_' ,region_label{r} ,'both.performance'));
perf=[perf; 1-perf_region(:)];

% stats per region
[locp(r),h(r),locstats]=ranksum(1-perf_region(:,2),1-perf_region(:,1));
Zval(r)=locstats.zval;

g1=repmat([ones(params.nperms,1);ones(params.nperms,1)*2],1,numel(region_label));
g1=g1(:);
end  

%% plotting
figure;
gr(1,1)=gramm('x', g2,'y', perf,'color',g2, 'lightness', g1);
gr(1,1).stat_boxplot('dodge',1.5,'width',1.7);
gr(1,1).set_names('x','','y','classification performance', 'color', {''});
gr(1,1).no_legend();
gr(1,1).set_color_options('map','d3_20');

%gr(1,1).set_layout_options('Position',[0.0662 0.0414 0.3 0.3]);
gr(1,1).facet_axes_handles.XLim=[-.75 5.25];
gr(1,1).facet_axes_handles.YLim=[0.1 0.6];
gr.draw();

%gr(1,1).results.stat_boxplot(2).outliers_handle.Marker='none';
%gr(1,1).results.stat_boxplot(3).outliers_handle.Marker='none';
%gr(1,1).results.stat_boxplot(4).outliers_handle.Marker='none';
%gr(1,1).results.stat_boxplot(7).outliers_handle.Marker='none';
gr(1,1).facet_axes_handles.XTick=[1 2 3 4];
gr(1,1).facet_axes_handles.XTickLabels=region_label;
set(gr(1,1).facet_axes_handles, 'Ylim',[0.15 0.4]);
%gr(1,1).facet_axes_handles.Position=[0.1321 0.0828 0.22 0.22];
     
% add labels of regions


function stats=getstats(svm)
% assemble stats from svm output
[locp(1),temp,locstats1]=ranksum(1-svm.results_dpu.performance(:,2),1-svm.results_dpu.performance(:,1));
[locp(2),temp,locstats2]=ranksum(1-svm.results_dpl.performance(:,2),1-svm.results_dpl.performance(:,1));
stats.locp=locp;
stats.temp=temp;
stats.locstats1=locstats1;
stats.locstats2=locstats2;


% rebuttal
% report median performance
median(1-svm.results_dpu.performance(:,1))
median(1-svm.results_dpl.performance(:,1))

% test against chance level 0.25
[p_st_u,h,stats_st_u]=signtest(1-svm.results_dpu.performance(:,1),0.25)
[p_st_l,h,stats_st_l]=signtest(1-svm.results_dpl.performance(:,1),0.25)

% tests for individual regions
[p_region(1),h,statshi]=signtest(1-svm.results_allfreqs_hiboth.performance(:,1),0.25)
[p_region(2),h,statsam]=signtest(1-svm.results_allfreqs_amboth.performance(:,1),0.25)
[p_region(3),h,statsec]=signtest(1-svm.results_allfreqs_ecboth.performance(:,1),0.25)
[p_region(4),h,statsphc]=signtest(1-svm.results_allfreqs_phcboth.performance(:,1),0.25)



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
        % HSL 225?, 100%, 97
    case 'myblue'
        phase_color=[107,174,214]./255;
     case 'mygray'
        phase_color=[0.3 0.3 0.3];
end

function [axhandle,labhandle,leghandle]=applyaxprops(axhandle,labhandle,leghandle)

if nargin<1
    axhandle=gca;
end
% properties for axis
props1={'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'TickLength',[0.01 0.0],...
'Color',[1 1 1]};

set(axhandle,props1{:});

if nargin==2
    
    
    
% properties for labels and legend
props2={'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal'};
set(labhandle,props2{:});
end


if nargin==3
set(leghandle,props2{:});
end


function [fhandle]=applyfigprops(fhandle)


% properties for axis
props1={'Color',[1 1 1],'paperorientation','landscape'};

set(fhandle(:),props1{:});