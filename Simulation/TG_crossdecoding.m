function [acc,react] = TG_crossdecoding(labels,data1,data2)

cfg = []; 
cfg.gamma = 0.05;

nSubs = size(data1,1);
nTrls = size(data1,2);

trainTP = size(data1,3);
testTP  = size(data2,3);

acc = zeros(nSubs,trainTP,testTP);
react = zeros(nSubs,nTrls,trainTP);
react_filter = zeros(nSubs,nTrls,trainTP);

for sub = 1:nSubs  
   fprintf('Crossdecoding sub %d out of %d \n', sub,nSubs)
   
   % decoding - generate distance to hyperplane for each trial and time
   % point
   Yhat = zeros(trainTP,testTP,nTrls);
     
   for t1 = 1:trainTP
       if mod(t1,3) == 0; fprintf('\t Time point %d out of %d \n',t1,trainTP); end
       
       for t2 = 1:testTP
           
           trainX = squeeze(data1(sub,:,t1,:));
           testX  = squeeze(data2(sub,:,t2,:));
           
           decoder = train_LDA(cfg,labels,trainX'); 
           Yhat(t1,t2,:)    = decode_LDA(cfg,decoder,testX'); 
       end
       
   end
   
    % convert to evidence correct class
   Yhat(:,:,labels==0) = Yhat(:,:,labels==0)*-1;
  
   acc(sub,:,:) = squeeze(mean(Yhat>0,3)); % accuracy
   
   % calculate reactivation
   for t = 1:trainTP
       
       yhat  = squeeze(Yhat(t,:,:));         
       
       % find peak
       for tr = 1:nTrls
           [~,t1] = max(yhat(:,tr));
           react(sub,tr,t) = t1;
       end   
       
%        %low pass filter
%        yhat(isnan(yhat)) = 0;
%        yhat = lowpass(yhat,30,300);       
%        
%        % find peak
%        for tr = 1:nTrls
%            [~,t1] = max(yhat(:,tr));
%            react_filter(sub,tr,t) = t1;
%        end   
   end
end

