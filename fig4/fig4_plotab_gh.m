function [stats]=fig4_plotab_gh(datafilepath)

load(datafilepath)
[fh(1)]=plotpanela(data.vex,data.vex_shuff,params);
[stats,fh(2)]=plotpanelb(data,params);

%[stats2,st]=geteffectsize(data)
%[fh]=ploteffectsize(st,params)

function [fh]=plotpanela(vex,vex_shuff,params)

%%
 [Nunits,nreps,nfreqs]=size(vex_shuff);
fh=figure;
set(gcf,'position',[72   821   314   474]);
imagesc([(vex(params.sortindex,:))]);
%imagesc([(vex(in,:))]);
clim_upper=mean(prctile(vex(:), [95 99]));
set(gca,'Clim',[0 clim_upper]);
set(gca,'Clim',[0 0.22]);
colorbar;
hold on;
hline(find(diff(params.sortcat))+1,'r');
%hline(find(diff(cat))+1,'r');

set(gca, 'Ytick',[],'Xtick', [1:2:nfreqs],'Xticklabel',num2str(round(params.freqs(1:2:nfreqs),1,'decimals')'));
% add text annotations for regions
txt={'HPC';'EC';'PHC';'AM'}; 
locs=[0.02 0.82;0.02 0.52;0.02 0.4;0.02 0.18];
for r=1:numel(txt)
annotation('textbox',[locs(r,1),locs(r,2) 0.1 0.1],'string',txt{r},'EdgeColor','none');
end
applyaxprops(gca);
applyfigprops(gcf);


%% panel b
function [stats,fh]=plotpanelb(data,params)
%%
pcrit=0.05; % for plotting
[Nunits,nreps,nfreqs]=size(data.vex_shuff);
permdata=squeeze(median(data.vex_shuff,2)); % across shuffles

% ttest
for j=1:nfreqs, [h(j),ptt(j)]=ttest(data.vex(:,j), ...
        permdata(:,j),'tail','right');
end
% signrank
for j=1:nfreqs, [psr(j),hsr(j)]=signrank(data.vex(:,j), ...
        permdata(:,j),'tail','right');
end
% shuffles: p values based on mean vex across units and mean across units for each shuffle 
mvex=mean(data.vex);
mvexshuff=squeeze(mean(data.vex_shuff,1));
p_shuff=1-mean(repmat(mvex,size(data.vex_shuff,2),1)>mvexshuff);

stats.ptt=ptt;
stats.h=h;
stats.psr=psr;
stats.hsr=hsr;
stats.p_shuff=p_shuff;

% plotting
fh(1)=figure;
offset2=max(mean(data.vex(:)))+0.02;

hfake=errorbar(1:nfreqs, squeeze(mean(permdata)),std(permdata)./sqrt(Nunits));
hold on;
h_all=errorbar(1:nfreqs,mean(data.vex),std(data.vex)./sqrt(size(data.vex,1)));

% hold on;
% plot(find(offset2.*(ptt<0.05)), nonzeros((offset2.*(ptt<0.05))),'*k');
% plot(find(offset2.*(ptt<0.01))+0.1, nonzeros((offset2.*(ptt<0.01))),'*k');

ylim([0.05 0.1]);
set(gca,'Xtick',1:2:nfreqs);
%set(gca,'Xticklabel',num2str(round(params.freqs(1:2:nfreqs),1,'decimals')'));
xlabel('frequency');
ylabel('variance explained');
set(h_all, 'Color', 'k','Marker','o','Color','k','Markerfacecolor','k','MarkerSize',8);
set(hfake(:),'Color', [0.8 0.8 0.8],'Marker','o','Color',[0.8 0.8 0.8],'Markerfacecolor',...
    [0.8 0.8 0.8],'MarkerSize',8);

legend({'vex-original','vex-shuffled'},'Fontsize',6);
legend(gca,'boxoff');


applyaxprops(gca);
applyfigprops(gcf);

function [stats,st]=geteffectsize(data)
pcrit=0.05; % for plotting
[Nunits,nreps,nfreqs]=size(data.vex_shuff);
a=data.vex;
b=data.vex_shuff;

%% hedges g
hg=[];
hgra=[];
% hedges g ci
hg_ci=[];
hgra_ci=[];
d=zeros(Nunits,nreps,nfreqs);
% compute standardized difference between correct and incorrect values
for f=1:12 % freqs

    for r=1:nreps
        % compute the same using the mes toolbox
        % this does not take into account NaN / Inf values
        [st]=mes(a(:,f),squeeze(b(:,r,f)),{'hedgesg'},'missVal','listwise','isDep',1);
        % save effect size measure and ci
        hg(f,r)=st.hedgesg;
        hg_ci(:,f,r)=st.hedgesgCi;
        d(:,r,f)=(a(:,f)-squeeze(b(:,r,f)))./(std([a(:,f);squeeze(b(:,r,f))])); % mean of d across units same as hedges g!!
    end

    
end
stats.p=1-mean(hg>0,2);
stats.d=d;
stats.hg=hg;

% compute effect size measure per unit and shuffle
hg_ra=[];
for ra=1:nreps
    [st_r]=mes(data.vex,squeeze(data.vex_shuff(:,ra,:)),{'hedgesg'},'missVal','pairwise','isDep',1);
    hg_ra(ra,:)=st_r.hedgesg;
end
stats.hg_ra=hg_ra;


function [fh]=ploteffectsize(st,params)

%% one imagesc plot using displaying effect sizes
% one line plot displaying effect size
[Nunits,nreps,nfreqs]=size(st.d);
dat=squeeze(mean(st.d,2));
fh(1)=figure;%3
set(gcf,'position',[72   821   314   474]);
imagesc(dat(params.sortindex,:));
set(gca,'Clim',[-2 2]);
colorbar;
hold on;
hline(find(diff(params.sortcat))+1,'r');
set(gca, 'Ytick',[],'Xtick', [1:2:nfreqs],'Xticklabel',num2str(round(params.freqs(1:2:nfreqs),1,'decimals')'));
% add text annotations for regions
txt={'HPC';'EC';'PHC';'AM'}; 
locs=[0.02 0.82;0.02 0.52;0.02 0.4;0.02 0.18];
for r=1:numel(txt)
annotation('textbox',[locs(r,1),locs(r,2) 0.1 0.1],'string',txt{r},'EdgeColor','none');
end
applyaxprops(gca);
applyfigprops(gcf);

fh(2)=figure;%4
%subplot(1,2,1);
errorbar(1:nfreqs,mean(st.hg_ra),...
    std(st.hg_ra),'Color', 'b','Marker','o','Markerfacecolor','b','MarkerSize',8);
ylim([-0.3 0.3]);
legend({'hedges g'},'Fontsize',6);
applyaxprops(gca);
applyfigprops(gcf);



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

function cp=correctsimes(pvals,varargin)
% function cp=correctsimons(p,N)
% ncomp = size(p,1), 
% doesn't work on exactly same p vals!!!

%pvals=iscolumnvector(pvals);

[ncomp, ncases]=size(pvals);
cp=zeros(ncomp, ncases);

% simons correction factors
cr=repmat((ncomp./[1:ncomp])',1,ncases);

% sort p values
[sp,ix]=sort(pvals);

% multiply with correction factor
temp_p=sp.*cr;

% bring back in correct order per comparison
% with loop, no easy way?!
%cp(ix)=(temp_p);
for k=1:size(sp, 2)
    cp(ix(:,k),k)=temp_p(:,k);
end


%cp=iscolumnvector(cp);
