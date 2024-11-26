function [stats]=fig4_plotc(datafilepath)

load(datafilepath)
[nunits,nfreqs,nshuff]=size(data.vex_correct_shuff);

% difference to shuffled vex
% correct trials
mcor_shuff=squeeze(median(data.vex_correct_shuff,3));
sdcor_shuff=squeeze(std(data.vex_correct_shuff,[],3));
nVEXcor=(data.vex_correct-mcor_shuff)./sdcor_shuff;
% incorrect trials
mincor_shuff=squeeze(median(data.vex_incorrect_shuff,3));
sdincor_shuff=squeeze(std(data.vex_incorrect_shuff,[],3));
nVEXincor=(data.vex_incorrect-mincor_shuff)./sdincor_shuff;

%% figure 1
figure;
ya=prctile([nVEXcor(:);nVEXincor(:)],75);
errorbar(1:nfreqs, nanmean(nVEXincor),...
   nanstd(nVEXincor)./sqrt(nunits),'Color', [0.8 0.8 0.8],'Marker','o','Markerfacecolor',[0.8 0.8 0.8],'MarkerSize',8);
set(gca, 'Xtick', [1:2:nfreqs],'Xticklabel',num2str(round(params.freqs(1:2:nfreqs),1,'decimals')'));
hold on;
errorbar(1:nfreqs,nanmean(nVEXcor),...
    nanstd(nVEXcor)./sqrt(nunits),'Color', 'b','Marker','o','Markerfacecolor','b','MarkerSize',8);


legend({'incor trials','cor trials'});
legend(gca,'boxoff');
ylabel('norm. vex cor / vex incor ');
xlabel('freqs (Hz)');

ylim([-ya ya]);
ylim([-0.7 0.7]);
applyfigprops(gcf);
applyaxprops(gca);

stats=doshuffletest(nVEXcor,nVEXincor)
[fh]=plotnormdiff(data,params);

function [fh]=plotnormdiff(data,params)
[nunits,nfreqs,nshuff]=size(data.vex_correct_shuff);
for ra=1:nshuff
    [st_r]=mes(data.vex_correct,data.vex_correct_shuff(:,:,ra),{'hedgesg'},'missVal','pairwise','isDep',1);
    hg_rac(ra,:)=st_r.hedgesg;
end

for ra=1:nshuff
    [st_r]=mes(data.vex_incorrect,data.vex_incorrect_shuff(:,:,ra),{'hedgesg'},'missVal','pairwise','isDep',1);
    hg_rai(ra,:)=st_r.hedgesg;
end

pval_permdiff=1-mean(hg_rac>hg_rai);
%% figure 2 - hedges g for correct and incorrect 
fh=figure;

ya=prctile([hg_rac(:);hg_rai(:)],99);
errorbar(1:nfreqs,mean(hg_rai),...
    std(hg_rai),'Color', [0.8 0.8 0.8],'Marker','o','Markerfacecolor',[0.8 0.8 0.8],'MarkerSize',8);
set(gca, 'Xtick', [1:2:nfreqs],'Xticklabel',num2str(round(params.freqs(1:2:nfreqs),1,'decimals')'));
hold on;
errorbar(1:nfreqs,mean(hg_rac),...
    std(hg_rac),'Color', 'b','Marker','o','Markerfacecolor','b','MarkerSize',8);



legend({'incor trials','cor trials'});
legend(gca,'boxoff');
ylabel('hedges g cor/incor');
xlabel('freqs (Hz)');
ylim([-ya ya]*3);
ylim([-0.7 0.7]);
applyfigprops(gcf);
applyaxprops(gca);

function [stats]=doshuffletest(nVEXcor,nVEXincor)
%% shuffle test between correct / incorrect
[nunits,nfreqs]=size(nVEXcor);
a=nVEXcor;
b=nVEXincor;
d=[];
nperms=1999;
pm=randpermmatrix(1:size(a,1)*2,nperms);
d_r=[];

% hedges g
hg=[];
hgra=[];
% hedges g ci
hg_ci=[];
hgra_ci=[];

% compute standardized difference between correct and incorrect values
for f=1:12 % freqs
%d(:,f)=(a(:,f)-b(:,f))./(std([a(:,f);b(:,f)])); % mean of d across units same as hedges g!!
% compute the same using the mes toolbox
% this does not take into account NaN / Inf values
[st]=mes(a(:,f),b(:,f),{'hedgesg'},'missVal','listwise','isDep',1);
% save effect size measure and ci
hg(f)=st.hedgesg;

end
%
% compute nperms permutations effect sizes / standardized differences
for ra=1:nperms
    for f=1:12 % use random shuffle for each frequency 
        tempdat=[a(:,f);b(:,f)];
        tempdat=tempdat(pm(:,ra)); % randomize group label
        %d_r(:,f,ra)=(tempdat(1:size(a,1))-tempdat(size(a,1)+1:end))./(std(tempdat(:)));
        % create 2 groups with shuffled labels, calculate effect size from shuffled labels as control
        a_r=tempdat(1:size(a,1));
        b_r=tempdat(size(a,1)+1:end);
        % NaN and Inf are not taken into account using this comparison
        [st_r]=mes(a_r,b_r,{'hedgesg'},'missVal','pairwise','isDep',1);
        % save effect size measure for random shufflings
        hgra(f,ra)=st_r.hedgesg;
        hgra_ci(:,f,ra)=st_r.hedgesgCi;
    end
end
%
% calculate pvalue from permutation test per frequency
% for both own standardized diff and from mes measure - should be same!
for f=1:12
%pval(f)=1-mean(mean(d(:,f))>squeeze(mean(d_r(:,f,:))));
pval_hg(f)=1-mean((hg(f)>hgra(f,:)));
end

% output:
stats.pval_hg=pval_hg; % pval hg from mes 
stats.hg=hg; % save original effect size between conditions
stats.hgra=mean(hgra,2); % save mean effect size across shuffles



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