function [cdata]=spk_cutdata(data,twindow)
% cuts out timestamps in between twindow(1) and twindow(2)
% input data (cell array which one trial per cell
% twindow ms , default is post fixation onset; need to change of data is
% aligned to different event
% expinfo

cdata=cellfun(@(data) data(data>=twindow(1)&data<=twindow(2)),data,...
    'Uniformoutput',false);

