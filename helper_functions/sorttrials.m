function [sortdata,sortinfo]=sorttrials(data,expinfo, stimids, stimpos)

% sort data according to stimulus conditions
% input is
% 1. cellarray of spktimes where each cell corresponds to one trial with spktimes
% 2. sort variables 'stimids', sortvalues 1-8 (e.g. [1], [1 8])
% stimpos, [1 4]
% if you omit 'stimpos' by default trials are sorted  by all positions 1-4
% if data is not a cell array, but matrix, function sorts data in
% trialdimension
% if data is empty, function just returns the sorted trial indices and
% empty sortdata var

%%

ntrials=size(expinfo.stimuli.stimpos,1);
sortinfo=[];

if isempty(stimpos)
    for i=1:numel(stimids)
        [rows,cols]=find(expinfo.stimuli.stimpos==stimids(i));
        sortinfo=[sortinfo;[rows repmat(stimids(i),size(rows,1),1) cols]];
    end
else
    
    for i=1:numel(stimids)
        [rows,cols]=find(expinfo.stimuli.stimpos==stimids(i));
        for j=1:numel(stimpos)
            % [trialidx stimid stimpos]
            sortinfo=[sortinfo;[rows(cols==stimpos(j)) repmat(stimids(i),...
                size(rows(cols==stimpos(j)))) cols(cols==stimpos(j))]];
            
        end
    end
end



if ~isempty(data)
    %% sort data according to new trial indices
    if iscell(data)
        sortdata=data(sortinfo(:,1));
    else
        dims=size(data);
        sortdim=find(dims==ntrials);
        if numel(dims)==1
            sortdata=data(sortinfo(:,1));
        elseif numel(dims)==2
            if sortdim==1
                sortdata=data(sortinfo(:,1),:);
            elseif sortdim==2
                sortdata=data(:,sortinfo(:,1));
            end
        elseif numel(dims)==3
            if sortdim==1
                sortdata=data(sortinfo(:,1),:,:);
            elseif sortdim==2
                sortdata=data(:,sortinfo(:,1),:);
            elseif sortdim==3
                sortdata=data(:,:,sortinfo(:,1));
            end
        end
    end
else
    sortdata=data;
end


