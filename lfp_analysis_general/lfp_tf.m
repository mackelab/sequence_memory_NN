function [Z,pw,ang]=lfp_tf(lfpdata,anaparams)
% [Z,pw,cs]=lfp_tf(lfpdata,anaparams)
% this function takes time!!!!!
% function calculated time-frequency spectograms based on
% anaparams.analysis and specified paramteres in anaparams
% lfpdata = trialsxtimepoints(ms,1000Hz)
% default anaparams from anaparams=lfp_defaults_analysis_seq(analysis)

cs = lfp_wavelet_decomposition_seq(lfpdata, anaparams.Fs, anaparams.freqs,[],numel(anaparams.freqs));
pw=abs(cs);
ang=angle(cs);
clear cs;

if anaparams.Zscoringbaseline
    % Zscoring
    % Z-Score using baseline values (first 500ms)
    mix=find(anaparams.t==0);
    mb=mean(mean(pw(:,:,1:mix),3),2);
    % mean baseline
    mb=(repmat(mb,1,size(pw,2),size(pw,3)));
    % std baseline
    sdb=std(std(pw(:,:,1:mix),[],3),[],2);
    sdb=repmat(sdb,1,size(pw,2),size(pw,3));
    Z=(pw-mb)./sdb;
else
    Z=[];
end

