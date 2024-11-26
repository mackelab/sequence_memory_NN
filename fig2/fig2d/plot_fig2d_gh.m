function [stats]=plot_fig2d_gh(datafilepath)
% same plot,  spike rates equalized 
% reported in paper: stats.ZRayleigh/stats.pRayleigh for estimation of
% non-uniformity per region,
% stats_ds.PAnglesRegions/stats_ds.pAnglesRegions comparing directions
% between regions per task phase

load(datafilepath);
data=ANG_ds;
% colors for lines
CO{1}=[0 0.3 0.8];% nPS
CO{2}=[0.8 0.8 0.8]; %  Baseline
CO{3}=[0.2 0.2 0.2]; %  PS

        YLIMS(1,:)=[0.0002 0.00035];
        YLIMS(2,:)=[0.00018 0.00033];
        YLIMS(3,:)=[0.00008 0.00018];
        YLIMS(4,:)=[0.00008 0.00018];
        
posis{1}=[0.3500    0.8200    0.1200    0.18;...
    0.500    0.8200    0.1200    0.18;...
    0.6500    0.8200    0.1200    0.18 ];
posis{2}=[0.3500    0.5600    0.1200    0.18;...
    0.500    0.5600    0.1200    0.18;...
    0.6500    0.5600    0.1200    0.18 ];
posis{3}=[0.3500    0.3400    0.1200    0.18;...
    0.500    0.3400    0.1200    0.18;...
    0.6500    0.3400    0.1200    0.18 ];
posis{4}=[0.3500    0.1200    0.1200    0.18;...
    0.500    0.1200    0.1200    0.18;...
    0.6500    0.1200    0.1200    0.18 ];

nregions=4;
nfreqs=12;
region{1}='HPC';
region{2}='PHC';
region{3}='AM';
region{4}='EC';
regionlabel=[];

figure('position',[1251 114 226 1095]);
datain_ps=[];
datain_nps=[];
datain_baseline=[];
for k=1:nregions
            subplot(4,1,k);
            title(region{k});
            locix=IX{k};
            angles{1}=data{1}(locix,:); % nps
            angles{2}=[data{2}(locix,:)]; % baseline
            angles{3}=[data{3}(locix,:)]; % ps
            [axh]=lineplotangles(angles, 1, 0.6, CO, 'global',0);
            
            [stats.pRayleigh(k,1),stats.ZRayleigh(k,1)]=circ_rtest(angles{3}(:)); % delay, ps
            [stats.pRayleigh(k,2),stats.ZRayleigh(k,2)]=circ_rtest(angles{2}(:)); % baseline
            %[stats.pRayleigh(k,3),stats.ZRayleigh(k,3)]=circ_rtest([angles{3}(:); angles{1}(:)]); % delay, ps and nps 
             [stats.pRayleigh(k,3),stats.ZRayleigh(k,3)]=circ_rtest(angles{1}(:)); % delay, nps

            % save delay mean phases per population to do comp between
            % populations
            datain_ps=[datain_ps;angles{3}(:)];
            datain_nps=[datain_nps;angles{1}(:)];
            datain_baseline=[datain_baseline;angles{2}(:)];
            
            regionlabel=[regionlabel;ones(numel(locix)*nfreqs,1)*k];
end

% compare all units preferred phase
[stats.pRayleigh_allunits(1),stats.ZRayleigh_allunits(1)]=circ_rtest(data{3}(:)); % PS delay
[stats.pRayleigh_allunits(2),stats.ZRayleigh_allunits(2)]=circ_rtest(data{2}(:)); % baseline
stats.pRayleigh=correctsimes(stats.pRayleigh,nfreqs*numel(region)*2);

% compare angles (direction) between regions
[stats.pAnglesRegions(1), medang, stats.PAnglesRegions(1)] = circ_cmtest(datain_ps,regionlabel); % delay ps
[stats.pAnglesRegions(2), medang, stats.PAnglesRegions(2)] = circ_cmtest(datain_baseline,regionlabel); % baseline
[stats.pAnglesRegions(3), medang, stats.PAnglesRegions(3)] = circ_cmtest(datain_nps,regionlabel); % nps
stats.pAnglesRegions=correctsimes(stats.pAnglesRegions',nfreqs*numel(region)*3);


labelstring{1}='blue:delay NPS'
labelstring{2}='gray:Baseline'
labelstring{3}='black:delay PS'
rh=annotation('textbox',[.1 .93 .35 .05],'String',labelstring,'Edgecolor','none');

function [p_kappa]= get_kappastats(ANG, IX)

for r=1:nregions
    
    d1=ANG{1}(IX{r},:); %nps
    b=ANG{2}(IX{r},:); % baseline
    d2=ANG{3}(IX{r},:); % ps


    % estimate kappas from population 
    [th(1),ka(1)]=circ_vmpar(d1(:));
    [th(2),ka(2)]=circ_vmpar(b(:));
    [th(3),ka(3)]=circ_vmpar(d2(:));


    % comparison groups
    ph1=[d1;b]; % nps vs baseline
    ph2=[d2;b]; % ps vs baseline
    n=size(ph1,1);
    la=repmat([1 2],n/2,1);
    la=la(:);

    %
    nperms=1999;
    n=randpermmatrix(la(:),nperms);
    for r=1:nperms

        % nps
        loc1=ph1(n(:,r)==1,:);
        % baseline
        loc2=ph1(n(:,r)==2,:);
        % ps
        loc3=ph2(n(:,r)==1,:);
        
        [a,kar(1,r)]=circ_vmpar(loc1(:));
        [a,kar(2,r)]=circ_vmpar(loc2(:));
        [a,kar(3,r)]=circ_vmpar(loc3(:));


    end

    p(1,r)=mean((kar(1,:)-kar(2,:))>ka(1)-ka(2));% nps
    p(2,r)=mean((kar(3,:)-kar(2,:))>ka(3)-ka(2));% ps
    
end


% 
% figure;
% plot_hist_group(kar(1,:)-kar(2,:));
% hold on;
% vline(ka(1)-ka(2));

p=mean((kar(1,:)-kar(2,:))>ka(1)-ka(2))

function [axh]=lineplotangles(angles,d,radiuslimit, cols,norm,plotlines)
% function [axh]=lineplotangles(angles,d,radiuslimit,cols,norm)
% angles in radians, d is whether to plot direction only with normalized
% length (d=1), or real vector length (d=0)
% r radius limits
% cols is cell with ncolors entries for color of plotting (eg cols{1}=[0 0 1]) default provided
% is argument omitted, otherwise chose ncolors for nlines to plot
% last argument is normalization to max ('local' within each condition,
% 'global', across all angles, default: 'local'

% check dimensionality of data, needs to be 1dvector
for k=1:max(size(angles))
    angles{k}=angles{k}(:);
end


if nargin<2
    d=1;
end

if nargin<3
    radiuslimit=0.75;
end

if isempty(radiuslimit)
    radiuslimit=1;
end

if nargin<4
cols{1}=[0 0 1];
cols{2}=[1 0 0];
cols{3}=[0 1 0];
cols{4}=[0 0 0];
cols{5}=[1 0 1];
cols{6}=[0 1 1];
cols{7}=[1 0 1];
cols{8}=[0.5 0.5 0.5];

end

if nargin<5
    norm='local';
end

if nargin<6
    plotlines=0;
end

if isempty(cols)
    cols{1}=[0 0 1];
cols{2}=[1 0 0];
cols{3}=[0 1 0];
cols{4}=[0 0 0];
cols{5}=[1 0 1];
cols{6}=[0 1 1];
cols{7}=[1 0 1];
cols{8}=[0.5 0.5 0.5];
end


% params
th=linspace(-pi,pi+0.00001,20);


for j=1:max(size(angles))
    if isnumeric(angles)
        data=angles(:,j);
        R=histc(angles(:),th);R(end)=R(1);
    end
    if iscell(angles)
        data=angles{j}(:);
        angtemp=cell2mat(cellfun(@(c) c', angles, 'uni',0));
        R=histc(angtemp,th);
        R(end)=R(1);
    end
    %[t,r]=rose(data);
    [r]=histc(data,th);
    r(end)=r(1);
    % normalization to maximum
    switch norm
        % local normalization per condition
        case 'local'
            r=r./max(r);
        % global normalization across condition
        case 'global'
            r=r./max(R);
    end
    

    
    t = 0 : .01 : 2 * pi;
    P = polar(t, radiuslimit * ones(size(t)));
    set(P, 'Visible', 'off');
    
    hold on;
    phan=polar(th,r','k')
    
    set(phan,'Color',cols{j});
    if ~plotlines
    if j==1
        hold all;
        h = findall(gca,'type','line');
        h(h == phan) = [];
        delete(h);
    end
    end
    
    hold on;

end

cm=circ_mean([angles{1}(:);angles{3}(:)]);
hold on;
hp(j)=polar([0 cm],[0 radiuslimit]);
set(hp(j),'Color',[0 0.3 0.7],'Linewidth',1);

axh.ax=gca;


%% replace text 

tex = findall(gca,'type','text');
tx=get(tex,'String');
% delete all text exept for labels of 0, 90, 180 , 270
% and replace with units of pi
ii=strcmp({'0'},tx);
set(tex(ii),'String','0');

ii=strcmp({'90'},tx);
set(tex(ii),'String','pi/2');

ii=strcmp({'180'},tx);
set(tex(ii),'String','pi');

ii=strcmp({'270'},tx);
set(tex(ii),'String','3/2pi');

ii=strcmp({'0'},tx)|strcmp({'90'},tx)|strcmp({'180'},tx)|strcmp({'270'},tx)|strcmp({'  0.5'},tx)|strcmp({'  1'},tx)...
    |strcmp({'  0.2'},tx)|strcmp({'  0.4'},tx)|strcmp({'  0.6'},tx)|strcmp({'  0.8'},tx);
set(tex(~ii),'String','');




