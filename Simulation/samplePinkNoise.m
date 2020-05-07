function data = samplePinkNoise (X,T,nSubs,nTPP,nTrls,nChanInArea,nBrainAreas,noiseW,noiseP)

data = noiseW * randn(nSubs,nTPP,nTrls*2,nChanInArea*nBrainAreas);
for subj = 1:nSubs % we use one channel of real data per subject, for simplicity
    % use some data
    x = reshape(zscore(reshape(X(:,subj),T(1),length(T))),T(1)*length(T),1);
    options = struct();
    options.order = 1;
    for j = 1:nChanInArea*nBrainAreas
        Tnew = nTrls*2*nTPP+1000;
        y = simmar(x,T,Tnew,options);
        y = reshape(y(501:end-500),nTPP,(Tnew-1000)/nTPP);
        data(subj,:,:,j) = squeeze(data(subj,:,:,j)) + noiseP * zscore(y);
    end
end
data = permute(data,[1 3 2 4]);


% % Create data with pink noise
% noiseW = 0.1; noiseP = 0.5; 
% dataP = noiseW * randn(nTrls*2*nTPP+1000,nSubs,nChanInArea*nBrainAreas);
% for sub = 1:nSubs
%     for j = 1:nChanInArea*nBrainAreas
%         dataP(:,sub,j) = smooth(dataP(:,sub,j),smooth_par);
%     end
% end
% dataP = noiseP * dataP(501:end-500,:,:);

end