function fig4_plote_su_gh(datafilepath)
% function fig4_plote_su_gh(datafilepath)
% this function requires the gramm plotting toolbox for matlab

load(datafilepath);

[nperms,ngroups,nunits_upper]=size(svm.perf_upper);
[nperms,ngroups,nunits_lower]=size(svm.perf_lower);

% stats compare means across units 
% mean across units perf

% upper % report
% mean across shuffles permtest
du1=squeeze(mean(svm.perf_upper))';
pu=1-mean((du1(:,1)<du1(:,2)));
% mean across units % reported in rebuttal (new)
du=squeeze(mean(svm.perf_upper,3));
pu_shuff=1-mean((du(:,1)<du(:,2)));
if pu_shuff==0
    pu_shuff=1/nperms;
end

% comparison shuffles original vs permuted (old, in paper)
% upper selectivity units
[psr_u_shuff, h,statssr_u_shuff]=signrank(du(:,1), du(:,2));
% test against chance level - reported in rebuttal (new)
[psr_u, h,statssr_u]=signrank(1-du(:,1),0.25)
% is same

% lower % report
% mean across shuffles permtest
dl1=squeeze(mean(svm.perf_lower))';
pl=1-mean((dl1(:,1)<dl1(:,2)));

% mean across units % reported in rebuttal (new)
dl=squeeze(mean(svm.perf_lower,3));
pl_shuff=1-mean((dl(:,1)<dl(:,2)));
if pl_shuff==0
    pl_shuff=1/nperms;
end

% test against chance level across units 
[psr_l, h,statssr_l]=signrank(1-dl(:,1),0.25)
% comparison shuffles (means across units) original vs permuted 
% lower selectivity units
[psr_l_shuff, h,statssr_l_shuff]=signrank(dl(:,1), dl(:,2));
% is same 

allp(:,1)=1-[reshape(squeeze(svm.perf_upper(:,1,:)),nperms*nunits_upper,1)];%; reshape(squeeze(perf_lower(:,1,:)),nperms*nunits_lower,1)];
allp(:,2)=1-[reshape(squeeze(svm.perf_upper(:,2,:)),nperms*nunits_upper,1)];% ; reshape(squeeze(perf_lower(:,2,:)),a*c,1)];

allp_l(:,1)=1-reshape(squeeze(svm.perf_lower(:,1,:)),nperms*nunits_lower,1);
allp_l(:,2)=1-reshape(squeeze(svm.perf_lower(:,2,:)),nperms*nunits_lower,1);

alldata=[allp(:);allp_l(:)];

c{1}='original';
c{2}='shuffled';

% groups shuffled vs non shuffled
groups_u=[repmat(c(1),nunits_upper*nperms,1);repmat(c(2),nunits_upper*nperms,1)];
groups_l=[repmat(c(1),nunits_lower*nperms,1);repmat(c(2),nunits_lower*nperms,1)];


% group label upper vs lower
group1=[[ones(nunits_upper*nperms,1);ones(nunits_lower*nperms,1)*2]];%; [ones(nunits_upper*nperms,1);ones(nunits_lower*nperms,1)*2]];
% group labels shuffles vs original data
group2=[ones(nunits_upper*nperms*ngroups,1); ones(nunits_lower*nperms*ngroups,1)*2];
edgv=linspace(min([allp(:);allp_l(:)]),max([allp(:);allp_l(:)]),15);
phase_colors=cell2mat(p_colors('uni_blue'));
edgv2=linspace(min(allp(:))-max(allp_l(:)), max(allp(:))-min(allp_l(:)),15);

% update to two colors: grey for shuffles, blue for non shuffled
fh(2)=figure;
set(2,'position',[97    27   954   678]);

g(1,1)=gramm('x', allp(:),'color',groups_u);
g(1,1).stat_bin('geom','overlaid_bar','normalization','probability','edges',edgv); %Overlaid bar automatically changes bar coloring to transparent
g(1,1).set_title('dp_u');
g(1,1).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
g(1,1).set_names('x','performance','y','probability', 'color', {''});
g(1,1).no_legend();
g(1,1).set_layout_options('Position',[0.0662 0.5579 0.35 0.35]);

g(1,2)=gramm('x', allp_l(:),'color',groups_l);
g(1,2).stat_bin('geom','overlaid_bar','normalization','probability','edges',edgv); %Overlaid bar automatically changes bar coloring to transparent
g(1,2).set_title('dp_l');
g(1,2).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
g(1,2).set_names('x','performance','y','probability', 'color', {''});
g(1,2).no_legend();
g(1,2).set_layout_options('Position',[0.5662 0.5579 0.35 0.35])


g(2,1)=gramm('x', group2,'y', [allp(:);allp_l(:)],'color',[group1;group1]);
g(2,1).stat_boxplot('dodge',0.3,'width',0.5);
g(2,1).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
g(2,1).set_names('x','','y','probability', 'color', {''});
g(2,1).no_legend();
g(2,1).set_layout_options('Position',[0.0662 0.0414 0.35 0.35]);
g(2,1).set_names('x','','y','probability', 'color', {'blue:original';'gray:shuffled'});

%
g(2,2)=gramm('x',[allp(:,1);allp_l(:,1)] , 'y',[allp(:,2);allp_l(:,2)] ,'color',group1);
%g(2,2)=gramm('x',[ppu(:,1);ppl(:,1)] , 'y',[ppu(:,2);ppl(:,2)] ,'color',groupn);




g(2,2).geom_point('alpha',0.05);
g(2,2).stat_cornerhist('edges',edgv2,'aspect',0.5);
g(2,2).geom_abline();
g(2,2).set_names('x','non-shuffled', 'y','shuffled','color',{});
g(2,2).set_color_options('map',[phase_colors(4,:);[0.5 0.5 0.5]],'n_color',2,'n_lightness',1);
g(2,2).set_names('x','performance','y','probability', 'color', {'blue:original';'gray:shuffled'});
g(2,2).set_layout_options('Position',[[0.5662 0.0579 0.35 0.35]]);

g.draw();


% set some props after
gr(2,1).facet_axes_handles.Position=[0.1321 0.0828 0.22 0.22];
g(2,1).facet_axes_handles.Position=[0.1321 0.0828 0.22 0.22];
g(2,1).facet_axes_handles.XLim=[-1 4];
g(2,2).facet_axes_handles.XLim=[0.1 0.6];
g(2,2).facet_axes_handles.YLim=[0.1 0.6];
set(g(2,1).facet_axes_handles,'XtickLabel',{},'Xtick',[1 2]);
g(1,1).facet_axes_handles.YLim=[ 0 0.35];
g(1,2).facet_axes_handles.YLim=[ 0 0.35];
g(2,1).results.stat_boxplot(1).box_handle.FaceAlpha=0.5;
g(2,1).results.stat_boxplot(1).outliers_handle.Marker='none';
g(2,1).results.stat_boxplot(2).box_handle.FaceAlpha=0.5;
g(2,1).results.stat_boxplot(2).outliers_handle.Marker='none';
g(2,1).axe_property('YLim',[-0.02 0.4],'XLim',[0.1 0.6]);

% set point size
set([g(2,2).results.geom_point_handle],'MarkerSize',2.8);

% add vline to histograms
axes(g(1,1).facet_axes_handles);
vline(0.25,'k:');
axes(g(1,2).facet_axes_handles);
vline(0.25,'k:');
g(2,1).facet_axes_handles.Position=[0.1321 0.0828 0.22 0.22];






