function [X,Y,trialnumbers] = prepDataDecoding(X,classIdx,balancemethod,Trl)
% function [X,Y,trialnumbers] = prepDataDecoding(X,classIdx,balancemethod,Trl)
% 
% Prepares data for decoding. Selects the trials for each class and
% balances the number of trials per class based on balancemethod. 
%
%         X             = data matrix with trials x features
%         classIdx      = cell array of logical trial indices per class
%         balancemethod = 'upsample' or 'downsample' to balance classes,
%         Trl           = optional field, cell array of trial numbers per class 
%         when using for cross-validation, use 'downsample' to prevent
%         same trials going in train and test folds
%
%         See also BALANCE_TRIALS
%
%         Created by Nadine Dijkstra, March. 2017

if ischar(balancemethod) % balance the trials or not?
    balance = 1;
else 
    balance = 0;
end

if nargin < 4
    Trl = 0;
end

% get data per class
nCond        = numel(classIdx);
classData    = cell(nCond,1);
labels       = cell(nCond,1);

% make three dimensions
nDim = length(size(X));
if nDim < 3
    X(:,:,2) = 0;
end

for c = 1:nCond
classData{c} = (X(classIdx{c},:,:));
labels{c}    = zeros(size(find(classIdx{c}==1),1),1)+c;
end

% concatenate over classes
X            = cat(1,classData{:});
Y            = cat(1,labels{:});
if iscell(Trl)
    trialnumbers = cat(1,Trl{:});
else 
    trialnumbers = [];
end

% balance the classes 
if balance
    idx          = balance_trials(Y,balancemethod);
    X            = X(cell2mat(idx),:,:);
    Y            = Y(cell2mat(idx),:,:);
    if iscell(Trl)
    trialnumbers = trialnumbers(cell2mat(idx),:,:);
    end
end
Y            = Y == 1;