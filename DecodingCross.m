function DecodingCross(cfg)
% function DecodingCross(subjectdata)
% trains on one data set and decodes the other one. If data files are the
% same, perform n-fold cross-validation.
%
% INPUT: cfg = configuration structure with the following fields:
%        cfg.dataFile        = 2 x 1 cell containing location of train{1} and test{2} data.
%        cfg.outputDir       = where to store the results
%        cfg.conIdx          = 2 x 2 cell with trigger values per dataset per class
%        cfg.channels        = on which channels to do the decoding, 'MEG',
%        selects MEG channels
%        cfg.Permutation     = 'true' or 'false' whether to permute labels, default is 'false'
%        cfg.nPermutations   = if Permutation is 'true', number of permutations
%        cfg.Visualising     = 'true' or 'false', whether to plot results, default is 'false'
%        cfg.Vividnessrelate = 'true' or 'false', whether to correlate activations with vividness ratings per trial, default is 'false'
%
%        cfg.decfg = configuration structure with the following fields:
%        cfg.decfg.gamma        = regularization parameter LDA classification
%        cfg.decfg.nMeanS       = n samples sliding window classification
%        cfg.decfg.samplemethod = method to balance classes,
%                    'downsample' or 'upsample'
%
%        subjectdata = structure containing subject specific information,
%                      such as location of data, output folder etc.
%
%        See also PREPDATADECODING, BALANCE_TRIALS, DECODINGCROSSTIME
%
%        Created by Nadine Dijkstra, April 2017


% outputDir
outputDir = cfg.outputDir;
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Fill in defaults cfg struct.
if ~isfield(cfg, 'Permutation'); cfg.Permutation = false; end
if ~isfield(cfg, 'Visualising'); cfg.Visualising = false; end
if ~isfield(cfg, 'Vividnessrelate'); cfg.Vividnessrelate = false; end
if ~isfield(cfg, 'appName'); cfg.appName = []; end
if ~isfield(cfg, 'channels'); cfg.channels = 'MEG'; end
if ~isfield(cfg, 'zscore'); cfg.zscore = false; end
if ~isfield(cfg, 'demean'); cfg.demean = false; end

%% Get the data
% load the training data
load(cfg.dataFiles{1})
cfgS             = [];
cfgS.channel     = cfg.channels;
dataStruct{1}    = ft_selectdata(cfgS,data);

% check whether test and train are the same data set
if strcmp(cfg.dataFiles{1},cfg.dataFiles{2})
    sameData = true; else sameData = false;
end

switch sameData
    case false % load the test data
        load(cfg.dataFiles{2})
        dataStruct{2}    = ft_selectdata(cfgS,data);
end
clear data

nDataSets = numel(dataStruct);

%% Prepare data for train and test set
decfg        = cfg.decfg;
X            = cell(nDataSets,1);
Y            = cell(nDataSets,1);
trialnumbers = cell(nDataSets,1);
dataNames    = cell(nDataSets,1);
time         = cell(nDataSets,1);

for d = 1:nDataSets
    
    % name of the data set for later saving
    [~,dataNames{d},~] = fileparts(cfg.dataFiles{d});
    time{d}            = dataStruct{d}.time{1,1};
    
    % determine size dimensions
    nChannels        = size(dataStruct{d}.trial{1,1},1);
    nTrials          = numel(dataStruct{d}.trial);
    nSamples         = size(dataStruct{d}.trial{1,1},2);
    
    % get class indices and trialnumbers
    nCon        = size(cfg.conIdx,2);
    classIdx    = cell(nCon,1);
    Trl         = cell(nCon,1);
    for cl = 1:nCon
        classIdx{cl} = ismember(dataStruct{d}.trialinfo,cfg.conIdx{d,cl});
        Trl{cl} = dataStruct{d}.trialnumbers(classIdx{cl}); % get trial numbers
    end
    
    % extract data in convenient format
    X{d} = reshape(cat(1,dataStruct{d}.trial{:}),[nChannels,nTrials,nSamples]);
    X{d} = permute(X{d},[2,1,3]);
    
    % select and balance trials
    [X{d},Y{d},trialnumbers{d}] = prepDataDecoding(X{d},classIdx,decfg.samplemethod,Trl);
    
    % zscore
    if cfg.zscore && d == 2 % remove time-locked info from test set
        X{d} = zscore(X{d});       
    end
    
    % demean per time point and per feature
    if cfg.demean
        for t = 1:size(X{d},3)
            for f = 1:size(X{d},2)
                X{d}(:,f,t) = X{d}(:,f,t) - squeeze(mean(X{d}(:,f,t),1));
            end
        end
    end
end

% clean up
clear dataStruct nDataSets nSamples nTrials nChannels

% name for saving
switch sameData
    case true
        if ~isempty(cfg.appName)
            name = sprintf('%s to %s_%s',dataNames{1},dataNames{1},cfg.appName);
        else
            name = sprintf('%s to %s',dataNames{1},dataNames{1});
        end
    case false
        if ~isempty(cfg.appName)
            name = sprintf('%s to %s_%s',dataNames{1},dataNames{2},cfg.appName);
        else
        name = sprintf('%s to %s',dataNames{1},dataNames{2});
        end
end

% permute or not
switch cfg.Permutation
    
    case false
        
    %% Do the decoding    
        switch sameData
            case true % apply cross-validation
                trueClass = Y{1};
                folds = createFolds(cfg,Y{1});
                Xhat  = zeros(size(X{1},3),size(X{1},3),size(X{1},1));
                for f = 1:cfg.nFold
                    fprintf('Running fold %d of %d \n',f,cfg.nFold)

                    testidx = folds{f}; trainidx = setdiff(1:numel(Y{1}),folds{f});
                    
                    % select train and test trials
                    labels{1} = Y{1}(trainidx,1);
                    labels{2} = Y{1}(testidx,1);
                    data{1}   = X{1}(trainidx,:,:);
                    data{2}   = X{1}(testidx,:,:);
                    
                    % decode per fold
                    Xhat(:,:,testidx) = decodingCrossTime(decfg,data,labels);
                end        
            case false
                trueClass = Y{2};
                folds = createFolds(cfg,Y{2}); % divide test data into folds
                Xhat  = zeros(size(X{1},3),size(X{2},3),size(X{2},1));                
                for f = 1:cfg.nFold
                    fprintf('Running fold %d of %d \n',f,cfg.nFold)
                    
                    testidx = folds{f}; testidxTnum = trialnumbers{2}(testidx);
                    trainidx = ~ismember(trialnumbers{1},testidxTnum); % deselect perception train data from same trial numbers as imagery test data
                    
                    % select train and test trials
                    labels{1} = Y{1}(trainidx,1);
                    labels{2} = Y{2}(testidx,1);
                    data{1}   = X{1}(trainidx,:,:);
                    data{2}   = X{2}(testidx,:,:);   
                    
                    % decode per fold
                    Xhat(:,:,testidx) = decodingCrossTime(decfg,data,labels);
                    
                end
        end
        %% Do the classification and save data
        % Classify
        class       = (Xhat > 0);
        Accuracy    = mean(class == permute(repmat(trueClass,[1,size(Xhat,1),size(Xhat,2)]),[2,3,1]),3);
        
        % Mean discriminative channel
        m0 = squeeze(mean(Xhat(:, :, trueClass==0), 3));
        m1 = squeeze(mean(Xhat(:, :, trueClass==1), 3));
        
        % p-values
        [~, p] = ttest2(Xhat(:, :, trueClass==0), Xhat(:, :, trueClass==1), 'dim', 3);
        p = squeeze(p);
        
        % save everything
        save(fullfile(outputDir,name),'Accuracy','m0','m1','p','Xhat','cfg','trialnumbers','time','trueClass','-v7.3')
        
    %% Permute the labels    
    case true
        Accuracy = [];
        Ytrue = Y;
        for per = 1:cfg.nPermutations
            fprintf('Permutation %d out of %d \n',per,cfg.nPermutations)
            
            % permute the labels of the training set
            Y     = Ytrue;
            Y{1}  = Ytrue{1}(randperm(size(Ytrue{1},1)),1);
            
            switch sameData
                case true % apply cross-validation
                    trueClass = Y{1};
                    folds = createFolds(cfg,Y{1});
                    Xhat  = zeros(size(X{1},3),size(X{1},3),size(X{1},1));
                    for f = 1:cfg.nFold
                        fprintf('Running fold %d of %d \n',f,cfg.nFold)
                        testidx = folds{f}; trainidx = setdiff(1:numel(Y{1}),folds{f});
                        
                        % select train and test trials
                        labels{1} = Y{1}(trainidx,1);
                        labels{2} = Y{1}(testidx,1);
                        data{1}   = X{1}(trainidx,:,:);
                        data{2}   = X{1}(testidx,:,:);
                        
                        % decode per fold
                        Xhat(:,:,testidx) = decodingCrossTime(decfg,data,labels);
                    end
                case false
                    trueClass = Y{2};
                    folds = createFolds(cfg,Y{2}); % divide test data into folds
                    Xhat  = zeros(size(X{1},3),size(X{2},3),size(X{2},1));
                    for f = 1:cfg.nFold
                        fprintf('Running fold %d of %d \n',f,cfg.nFold)
                        
                        testidx = folds{f}; testidxTnum = trialnumbers{2}(testidx);
                        trainidx = ~ismember(trialnumbers{1},testidxTnum); % deselect perception train data from same trial numbers as imagery test data
                        
                        % select train and test trials
                        labels{1} = Y{1}(trainidx,1);
                        labels{2} = Y{2}(testidx,1);
                        data{1}   = X{1}(trainidx,:,:);
                        data{2}   = X{2}(testidx,:,:);
                        
                        % decode per fold
                        Xhat(:,:,testidx) = decodingCrossTime(decfg,data,labels);
                        
                    end
            end
            class            = (Xhat > 0);
            Accuracy(per,:,:)= mean(class == permute(repmat(trueClass,[1,size(Xhat,1),size(Xhat,2)]),[2,3,1]),3);
            
        end
        save(fullfile(outputDir,[sprintf('Permuted(%d)',cfg.nPermutations) name]),'Accuracy','cfg','time','-v7.3','Xhat','trueClass')      
end


%% Visualisation

% visualising
switch cfg.Visualising
    case true
        figure;
        corrAcc = Accuracy;
        
        subplot(3,1,1); imagesc(time{2},time{1},corrAcc); colorbar
        title('Accuracy')
        xlabel('Imagery'); ylabel('Perception'); axis xy; axis image
        
        subplot(3,1,2); imagesc(testTime,trainTime,log10(p)); colorbar
        title(sprintf('%d versus %d log10 p',class1Num,class2Num))
        xlabel('Imagery'); ylabel('Perception'); axis xy; axis image
        
        subplot(3,1,3); imagesc(testTime,trainTime,m1-m0); colorbar
        title(sprintf('%d versus %d m1-m0',class1Num,class2Num))
        xlabel('Imagery'); ylabel('Perception'); axis xy; axis image
end

%% Relate to vividness
switch cfg.Vividnessrelate
    
    case true
        % get the vividness ratings
        load(fullfile(root,subjectdata.outputDir,'Behaviour','vividness'))
        
        % change to correct-class activation
        Activation  = Xhat;
        Activation(:, :, trueClass==0) = Activation(:, :, trueClass==0)*-1;
        
        % correlate to vividness rating
        viv  = vividness(trialnumbers{1});
        r    = zeros(size(Xhat,1),size(Xhat,1));
        pval = zeros(size(Xhat,1),size(Xhat,1));
        for s = 1:size(Xhat,1)
        [r(:,s),pval(:,s)] = corr(squeeze(Activation(s,:,:))',viv);
        end
        
        % show the things
        alpha = 0.05;
        figure;
        pcorr = r;
        pcorr(pval>alpha) = NaN;
        imagesc(time{1},time{1},pcorr); colorbar
        axis xy; axis image
        title(sprintf('Correlation vividness p < %.3f',alpha))
        xlabel('Imagery'); ylabel('Perception')
        
end

