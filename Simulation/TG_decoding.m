function acc = TG_decoding(cfg,labels,data)

% data is subjects x trials x timepoints x features

nSubs = size(data,1);
nTrls = size(data,2);
nTPs  = size(data,3);

folds     = createFolds(cfg,labels);
acc       = zeros(nSubs,nTPs,nTPs);
%Yhat      = zeros(nSubs,nTrls,nTPs,nTPs);
for sub = 1:nSubs
   
    fprintf('Decoding subject %d out of %d \n',sub,nSubs);
    X = squeeze(data(sub,:,:,:));
    Yhat = zeros(nTrls,nTPs,nTPs);
    
    % run n fold cross-validation
    for f = 1:cfg.nFold
       
        fprintf('\t Running fold %d out of %d \n',f,cfg.nFold)
        testidx = folds{f}; trainidx = setdiff(1:length(labels),folds{f});
        trainY = labels(trainidx);        
        
        % run over time points
        for s1 = 1:nTPs % train 
            
            if mod(s1,100)==0;
                fprintf('\t \t Training on timepoint %d out of %d \n',s1,nTPs)
            end
            
            trainX = squeeze(X(trainidx,s1,:));             
            decoder = train_LDA(cfg, trainY, trainX');
            
            for s2 = 1:nTPs % test
                
                testX  = squeeze(X(testidx,s2,:));
                Yhat(testidx,s1,s2) = decode_LDA(cfg,decoder,testX');
                %Yhat(sub,testidx,s1,s2) = decode_LDA(cfg,decoder,testX');
            end
        end     
        
    end  
    
    acc(sub,:,:) = mean(Yhat>0 == repmat(labels,1,nTPs,nTPs),1);
    %acc(sub,:,:) = mean(squeeze(Yhat(sub,:,:,:))>0 == repmat(labels,1,nTPs,nTPs),1);
    
end