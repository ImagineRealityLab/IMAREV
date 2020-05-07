function Xhat = decodingCrossTime(cfg,X,Y)
% function decodingCrossTime(cfg,X,Y)
%
% Trains a classifier on the data in X{1} and classifies the data in X{2}.
% 
%     INPUT: X{n}  = a trials x features x sample points matrix. 
%            Y{n}  = a trials x 1 vector containing the class labels
%            cfg   = a configuration structure with the fields:
%                .nMeanS = amount of sample points to average over, default
%                is 0.
%                .gamma  = the shrinkage regularisation parameter for the
%                LDA classification
%    OUTPUT: Xhat  = sample point x sample point x test trials matrix of
%                 classifier activations
%
%    See also TRAIN_LDA, DECODE_LDA
%
%    Created by Nadine Dijkstra April 2017

nSamplesTrain = size(X{1},3);
nSamplesTest  = size(X{2},3);
nTrialsTest   = size(X{2},1);

Xhat          = zeros(nSamplesTrain,nSamplesTest,nTrialsTest);

nMeanS        = cfg.nMeanS;


for s1 = 1:nSamplesTrain
    
    % define the training set
    if s1 <= nMeanS/2 || s1 >= nSamplesTrain - (nMeanS/2)
        Xhat(s1,:,:) = NaN;
    else
        train = squeeze(mean(X{1}(:,:,round(s1 - nMeanS/2):round(s1 + nMeanS/2)),3));
        if sum(isnan(train(1,:))) > 1
                Xhat(s1,:,:) = NaN;
        else
        
        if mod(s1,100) == 0
        fprintf('\t Training on sample %d out of %d \r',s1,nSamplesTrain);
        end
        
        % train the decoder
        decoder = train_LDA(cfg, Y{1}, train');
        
        for s2 = 1:nSamplesTest
            
            % define the testing set
            if s2 <= nMeanS/2 || s2 >= nSamplesTest - (nMeanS/2)
                Xhat(s1,s2,:) = NaN;
            else
                
                test = squeeze(mean(X{2}(:,:,round(s2 - nMeanS/2):round(s2 + nMeanS/2)),3));
                
                % decoding
                Xhat(s1,s2,:) = decode_LDA(cfg, decoder, test');

            end
        end
        end
    end
end

