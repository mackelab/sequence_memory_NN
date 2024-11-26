function [stats]=fig4_statsghij_gh(datafilepath)
% function [stats]=fig4_statsghij_gh(datafilepath)

load(datafilepath);
[nunits,npos]=size(X);
[anchoredphases , ref_phase] = anchor_phases(X,'mean');
[order,correctInverse,correctReverse,sortedphases] = spk_getorderfromphase(anchoredphases,'circular','mean');


% stats on order
[h(1),p(1), chi2stat(1),df] = prop_test([sum(correctInverse) 1/6*nunits],[nunits nunits],'correct');
[h(2),p(2), chi2stat(2),df] = prop_test([sum(correctReverse) 1/6*nunits],[nunits nunits],'correct');

stats.prop_order=[mean(correctInverse) mean(correctReverse)];
stats.prop_test.h=h;
stats.prop_test.p=p;
stats.prop_test.chi2stat=chi2stat;

%% compute time differences between positions based on taken frequeny for ordered phases histograms
% 0 (adjacent stimuli)
phaseshift(:,1)=abs(rad2deg(circ_dist(sortedphases(:,1),sortedphases(:,2))));
phaseshift(:,2)=abs(rad2deg(circ_dist(sortedphases(:,2),sortedphases(:,3))));
phaseshift(:,3)=abs(rad2deg(circ_dist(sortedphases(:,3),sortedphases(:,4))));
for k=1:size(phaseshift,2)
    [timedelay1(:,k)]=phaseshift2timedelay(phaseshift(:,k), t.fpower);
end

% 1 (one in between)
phaseshift2(:,1)=abs(rad2deg(circ_dist(sortedphases(:,1),sortedphases(:,3))));
phaseshift2(:,2)=abs(rad2deg(circ_dist(sortedphases(:,2),sortedphases(:,4))));

for k=1:size(phaseshift2,2)
    [timedelay2(:,k)]=phaseshift2timedelay(phaseshift2(:,k), t.fpower);
end

% 2 (2 in between, only one condition)
phaseshift3(:,1)=abs(rad2deg(circ_dist(sortedphases(:,1),sortedphases(:,4))));
[timedelay3]=phaseshift2timedelay(phaseshift3, t.fpower);

stats.med_timedelay=median(timedelay1);
stats.med_phaseshift=median(phaseshift);

[P_timedelay1,ANOVATAB_timedelay1,STATS_timedelay1] = kruskalwallis(timedelay1,[],'off');
stats.P_timedelay1=P_timedelay1;
stats.ANOVATAB_timedelay1=ANOVATAB_timedelay1;
stats.STATS_timedelay1=STATS_timedelay1;

% test differences between phase differences between pairs
[pval1, med1, P1] =circ_cmtest(phaseshift(:,1),phaseshift(:,2));
[pval2, med2, P2] =circ_cmtest(phaseshift(:,1),phaseshift(:,3));
[pval3, med3, P3] =circ_cmtest(phaseshift(:,2),phaseshift(:,3));
stats.P_phasediff=[pval1 pval2 pval2];
stats.m_phasediff=[med1 med2 med3];
stats.TestP_phasediff=[P1 P2 P3];

function [diffphases , ref_phase] = anchor_phases(phases,mode)
% function [diffphases , ref_phase] = anchor_phases(phases,mode)
%given an array of phases (of size n trials by nphases), 'anchor' the
%phases, i.e. calculate the phase difference with respect to some reference
%phase
%mode tells us how we 'anchor' the phases
%mode = 'first'% means first position is set to 0 (and all phases are
%nonnegative)
%mode = 'mean' anchors phases relative to the circular mean , second dim
%mode = 'circgap' sets as the first phase the one which has the biggest gap
%from the previous phase
%default is 'mean'
if nargin == 1 mode ='mean'; end

[N,M]=size(phases);

switch mode
    case 'first'
        ref = phases(:,1);
    case 'mean'
        ref = circ_mean(phases,[],2);
    case {'gap','circgap'}
        for n=1:N
            locphases= mod(phases(n,:),2*pi);
            %find biggest gap-- for that, sort the data, and append first
            %elemtent +2*pi to also take 'wraparound' gap into account:
            [sortphases,sortindex]= sort(locphases,'ascend');
            diffphases = [sortphases]-[sortphases(end)-2*pi,sortphases(1:end-1)];
            [biggap,ind_biggap_sorted]=max(diffphases);
            %ind_biggap_sorted gives us the index of the angle (after
            %sorting) which is preceeded by the biggest gap
            ref(n,1)=sortphases(ind_biggap_sorted);
        end
end

diffphases=circ_dist(phases,repmat(ref,1,M));
ref_phase = ref;

function [order,correctI,correctR,varargout] = spk_getorderfromphase(diffphases,mode1,mode2)
% function [order,correctI,correctR,varargout] = spk_getorderfromphase(diffphases,mode1,mode2)
% input = diffphases i.e. recentered phases from function anchor_phases.m
% mode1 = 'circular', i.e. relative order counts (relative ordering of
% phases within the circle
% mode1 = 'absolute', absolute 0 point
% compute prop. of forward (1  2  3 4] and reverse orderings [4 3 2 1]
% mode2 = 'circgap', 'first', 'mean' for anchoring phases
% varargout{1} = sorted phases

if nargin<2
    mode1='circular';
    mode2='mean';
end

[N,n]= size(diffphases);

switch mode2
    case {'mean','circgap','gap'}
        [sortphases, order] = sort(diffphases,2,'ascend');
    case {'first'}
        [sortphases, order] = sort(diffphases(:,2:4),2,'ascend');
        sortphases=[diffphases(:,1) sortphases(:,1:3)];
        order=[ones(size(sortphases,1),1) order+1];
end


switch mode1
    case 'circular'
        correctI = zeros(N,1);
        for i=1:n
            locsort = circshift([1:n]',i)';
            correctI = max(correctI,min(order == repmat(locsort,N,1),[],2));
        end
        correctR = zeros(N,1);
        vec=[4 3 2 1];
        for i=1:n
            locsort = circshift(vec',i)';
            correctR = max(correctR,min(order == repmat(locsort,N,1),[],2));
        end

    case 'absolute'
        correctI = min(order == repmat([1:n],N,1),[],2);
        correctR = min(order == repmat([4 3 2 1],N,1),[],2);
end

varargout{1}=sortphases;


function [timedelay]=phaseshift2timedelay(phaseshift, frequency)
% computes the timedelay given a phaseshift (deg) , i.e. timedifference (in
% ms)
% based on phasedifference

phaseshift=iscolumnvector(phaseshift);
frequency=iscolumnvector(frequency);
timedelay=(phaseshift./(360.*frequency)).*1000;


function [h,p, chi2stat,df] = prop_test(X , N, correct)

% [h,p, chi2stat,df] = prop_test(X , N, correct)
%
% A simple Chi-square test to compare two proportions
% It is a 2 sided test with alpha=0.05
%
% Input:
% X = vector with number of success for each sample (e.g. [20 22])
% N = vector of total counts for each sample (e.g. [48 29])
% correct = true/false : Yates continuity correction for small samples?
%
% Output:
% h = hypothesis (H1/H0)
% p = p value
% chi2stat= Chi-square value
% df = degrees of freedom (always equal to 1: 2 samples)
%
% Needs chi2cdf from the Statistics toolbox
% Inspired by prop.test() in "R" but much more basic
%
% Example: [h,p,chi]=prop_test([20 22],[48 29], true)
% The above example tests if 20/48 differs from 22/29 using Yate's correction

if (length(X)~= 2)||(length(X)~=length(N))
    disp('Error: bad vector length')
elseif (X(1)>N(1))|| (X(2)>N(2))
    disp('Error: bad counts (X>N)')
else
    df=1; % 2 samples

    % Observed data
    n1 = X(1);
    n2 = X(2);
    N1 = N(1);
    N2 = N(2);

    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);

    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];

    if correct == false
        % Standard Chi-square test
        chi2stat = sum((observed-expected).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    else
        % Yates continuity correction
        chi2stat = sum((abs(observed - expected) - 0.5).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    end

    h=0;
    if p<0.05
        h=1;
    end
end



