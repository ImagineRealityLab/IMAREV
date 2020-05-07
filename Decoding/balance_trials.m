function idx = balance_trials(y,samplemethod,ntrials)
% function idx = balance_trials(y,samplemethod,ntrials)
% balances the number of trials per class for classification
%
%  INPUT:   y             = trials x 1 labels
%           balancemethod = 'downsample' or 'upsample'
%           ntrials (optional) = specify the number of trials per class
% OUTPUT:   idx           = cell structure with balanced indices per class
%
% Developed by Nadine Dijkstra 2016
% 


% for repeatability
rng(1,'twister');

if strcmp(samplemethod,'upsample') % upsample other classes
    [~,b] = max(histc(y,unique(y))); % class with most trials
    nclasses = numel(unique(y));
    if nargin < 3
    ntrials = size(find(y == b),1);
    end
    idx = cell(nclasses,1);
    
    for c = 1:nclasses
        indices = find(y == c);
        while numel(indices)<ntrials
            ind = indices;
            indices = [indices;ind(randi(numel(ind)))];
        end
        idx{c} = indices;
        clear indices
    end
    
elseif strcmp(samplemethod,'downsample') % downsample other classes
    [~,b] = min(histc(y,unique(y))); % class with least trials
    nclasses = numel(unique(y));
    if nargin < 3
        ntrials = size(find(y == b),1);
    end
    idx = cell(nclasses,1);
    
    for c = 1:nclasses
        indices = find(y == c);
        while numel(indices)>ntrials
            indices(randi(numel(indices))) = []; % randomly delete one
        end
        idx{c} = indices;
        clear indices
    end
end
