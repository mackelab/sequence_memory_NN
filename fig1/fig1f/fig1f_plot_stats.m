function fig1f_plot_stats(datafilepath)
% function fig1f_plot_stats(datafilepath)
% same as fig1e_data.mat

load(datafilepath,'data');
% baseline 1000ms, create timevector
%%
tvec=((1:size(data,1)) - 1000)-1;
% load times for averagaing windows
ptimes=loadparams_avgtimes; % load times for averaging
npos=4;
%%
for i=1:npos
    spkZ_sample(:,i)=squeeze(mean(data((ptimes.twins(i,1):ptimes.twins(i,2))+abs(tvec(1)),i,:)));
end
% for delay
for i=1:npos
    spkZ_delay(:,i)=squeeze(mean(data((ptimes.twindelay(1):ptimes.twindelay(2))+abs(tvec(1)),i,:)));
end

% for panel
for i=1:npos
    spkZ_panel(:,i)=squeeze(mean(data((ptimes.twinpanel(1):ptimes.twinpanel(2))+abs(tvec(1)),i,:)));
end

% boxplot figure
phasecolors=cell2mat(p_colors('yellow_red')');
gh=figure;
set(gh, 'position',[953 1011 625 277]);
subplot(1,3,1);
boxplot(spkZ_sample,'colors',phasecolors);
ylim([-2 15]);
ylabel('Zscored Spikerate');
xlabel('Position');
o=findobj(gca,'tag','Outliers');
set(o,'marker','none');
subplot(1,3,2);
boxplot(spkZ_delay,'colors',phasecolors);
o=findobj(gca,'tag','Outliers');
set(o,'marker','none');
ylabel('Zscored Spikerate');
xlabel('Position');
subplot(1,3,3);
boxplot(spkZ_panel,'colors',phasecolors);
o=findobj(gca,'tag','Outliers');
set(o,'marker','none');
ylim([-2 10]);
ylabel('Zscored Spikerate');
xlabel('Position');
applyfigprops(gcf);
applyaxprops(gca);

%% stats
% compare spikerates between positions per timewindow
% non parametric anova
[pkw_au(1,1),b,c]=kruskalwallis(spkZ_sample,[],'off');
Chi_au(1,1)=b{2,5};
[pkw_au(2),b,c]=kruskalwallis(spkZ_delay,[],'off');
Chi_au(2)=b{2,5};
[pkw_au(3),b,c]=kruskalwallis(spkZ_panel,[],'off');
Chi_au(3)=b{2,5};



%  anova1
[pa_au(1),b,c]=anova1(spkZ_sample,[],'off');
F_au(1)=b{2,5};
[pa_au(2),b,c]=anova1(spkZ_delay,[],'off');
F_au(2)=b{2,5};
[pa_au(3),b,c]=anova1(spkZ_panel,[],'off');
F_au(3)=b{2,5};


% print to screen
sprintf('Evoked Response F and p are %0.5g and %0.5g\n Chi2 and p are %0.5g and %0.5g.',F_au(1),pa_au(1),Chi_au(1),pkw_au(1))
sprintf('Delay F and p are %0.5g and %0.5g\n Chi2 and p are %0.5g and %0.5g.',F_au(2),pa_au(2),Chi_au(2),pkw_au(2))
sprintf('Panel F and p are %0.5g and %0.5g\n Chi2 and p are %0.5g and %0.5g.',F_au(3),pa_au(3),Chi_au(3),pkw_au(3))