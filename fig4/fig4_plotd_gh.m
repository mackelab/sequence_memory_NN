function [stats]=fig4_plotd_gh(datafilepath)
% this function requires the gramm toolbox for plotting
% Cave there is a small plotting bug: figure needs to be adjusted in size
% manually for the cornerhistogram to be aligned!


load(datafilepath);
[Nunits,nreps,nfreqs]=size(data.vex_shuff);

% index to get vex for max power in delay per unit
index=sub2ind(size(data.vex),1:Nunits, params.fixpower');
newdata_power=data.vex(index);
newdatashuff_power=squeeze(median(data.vex_shuff,2));
newdatashuff_power=newdatashuff_power(index);

% effect size
for u=1:size(data.vex,1)
    dummyshuff(u,:)=data.vex_shuff(u,:,params.fixpower(u));
end
d=[];
for r=1:nreps
    d(:,r)=(newdata_power'-squeeze(dummyshuff(:,r)))./(std([newdata_power';squeeze(dummyshuff(:,r))])); % mean of d across units same as hedges g!!
end

% same as d just checking! 
% compute effect size measure per unit and shuffle
hg_ra=[];
for ra=1:nreps
    dummyshuff=squeeze(data.vex_shuff(:,ra,:));
    dummyshuff=dummyshuff(index);
    [st_r]=mes(data.vex(index),dummyshuff,{'hedgesg'},'missVal','pairwise','isDep',1);
    hg_ra(ra,:)=st_r.hedgesg;
end
stats.hg_ra=hg_ra;

phase_colors=cell2mat(p_colors('green_blue'));
figure('position',[300   741   532   577]);
% prep data for plotting
x=[newdata_power(:)];
y=[newdatashuff_power(:)];
group=repmat([1],Nunits,1);
group=group(:);
numbins=45;
edgv=linspace(-0.3,0.3,numbins);

clear g
g(1)=gramm('x',x,'y',y);
g(1).geom_point();
g(1).axe_property('YLim',[-0.02 0.4],'XLim',[-0.02 0.4]);
g(1,1).set_color_options('map',[0.6 0.6 0.6],'n_color',1,'n_lightness',1)
g(1).stat_cornerhist('edges',edgv);
g(1).geom_abline();
g(1,1).no_legend();
g.draw();

% collect axis handles for corner histograms
ch_axh1=g.results.stat_cornerhist.child_axe_handle;

% update to only plot again upper 50% percentile of all data and show that
% separately in histograms
up50=prctile([x;y],50);
ii(:,1)=x>=up50;
ii(:,2)=y>=up50;
index=sum(ii,2)>0;

g.update('x',x(index),'y',y(index));
g(1).geom_point();
g(1).axe_property('YLim',[-0.02 0.4],'XLim',[-0.02 0.4]);

g(1,1).set_color_options('map',phase_colors(3,:),'n_color',1,'n_lightness',1)
g(1).stat_cornerhist('edges',edgv);
g(1).geom_abline();
g(1,1).no_legend();
g(1,1).set_names('y','vex-shuffled trials','x' ,'vex - original');
g.set_layout_options('Position',[0.5 0.1 0.2 0.2]);
g.draw();

ch_axh2=g.results.stat_cornerhist.child_axe_handle;
g.facet_axes_handles.Position=[0.5 0.05 0.35 0.35];
set(ch_axh1,'Position',[310 150 400 80],'Visible','off');
set(ch_axh2,'Position',[310 150 400 80]);

%% stats
% reported in paper
[h,pval,ci,statstt]=ttest(newdata_power,newdatashuff_power);
[psr,hsr]=signrank(newdata_power,newdatashuff_power);

% test only upper ptile, i.e. units with high vex
[h_up,p_up,x,statsup]=ttest(newdata_power(index),newdatashuff_power(index));
[psr_up,hsr_up]=ttest(newdata_power(index),newdatashuff_power(index),'tail','right');

stats.h=h;
stats.p=pval;
stats.ci=ci;
stats.T=statstt.tstat;
stats.h_up=h_up;
stats.p_up=p_up;
stats.Tup=statsup.tstat;



% helper functions
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
