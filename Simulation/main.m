clear all;
rng(1); 
restoredefaultpath;
addpath(genpath('/vol/ccnlab1/naddij/Analyses/HMM-MAR-master'))
simulationDir = fullfile('/vol/ccnlab1/naddij/Analyses/MEGREV/Simulation');
cd(simulationDir);

% contains data X, time by sensors; and T with the lenghts of the trials
real_data = load('example_data','X','T');


%% Simulations for eLife revision (2020) MEGREV 
% Simulations should show:
% - Whether trial by trial distance measures can dissociate between same
%   and reversed order of reactivations
% - Do that better than decoding accuracy (King & Dehaene, 2014) 
% - The effect of temporal jitter between trials 


%% Settings
Fs     = 300; 
nTrls  = 100; % per class - both for I and P 
nSubs  = 10; % 
nClass = 2;
nTPP   = 19; % perception time points
nTPI   = 600; % imagery time points
nBrainAreas = 5; 
nChanInArea = 4;  

labels    = [ones(nTrls,1); zeros(nTrls,1)];

% time vectors 
timeP     = 0:(1/Fs):(nTPP-1)*(1/Fs);
timeI     = 0:(1/Fs):(nTPI-1)*(1/Fs);

%% Model the feedforward flow

models = zeros(nChanInArea*nBrainAreas,nTPP,nClass); 
tPA = 7; % # time points each brain area is active
OL  = 2; % how much each brain area overlaps
for s = 1:nBrainAreas    
    idx = s+((s-1)*OL):s+((s-1)*OL)+tPA-1; 
    idx(idx>nTPP) = nTPP;
    ids = (1:nChanInArea) + nChanInArea*(s-1);
    models(ids,idx,1) = repmat(rand(nChanInArea,1),1,length(idx));
    models(ids,idx,2) = repmat(rand(nChanInArea,1),1,length(idx));
end

figure(1);
subplot(2,2,1); imagesc(squeeze(models(:,:,1))); title('Faces'); xlabel('Model');
ylabel('sensor'); caxis([-2 2]); colorbar
subplot(2,2,2); imagesc(squeeze(models(:,:,2))); title('Houses'); xlabel('Model');
ylabel('sensor'); caxis([-2 2]); colorbar
subplot(2,2,3); imagesc(corr(squeeze(models(:,:,1)))); axis xy; 
xlabel('Model'); ylabel('Model'); colorbar; title('Correlation faces')
subplot(2,2,4); imagesc(corr(squeeze(models(:,:,2)))); axis xy; 
xlabel('Model'); ylabel('Model'); colorbar; title('Correlation houses')


time_factor = 10; % slow down for imagery
models_slow = zeros(nChanInArea*nBrainAreas,time_factor,nTPP,nClass); 
for t = 1:time_factor, models_slow(:,t,:,:) = models; end
models_slow = reshape(models_slow,[nChanInArea*nBrainAreas,time_factor*nTPP,nClass]); 
figure(2);
subplot(2,1,1); imagesc(squeeze(models_slow(:,:,1))); title('Faces'); xlabel('Model');
ylabel('sensor'); caxis([-2 2]);
subplot(2,1,2); imagesc(squeeze(models_slow(:,:,2))); title('Houses'); xlabel('Model');
ylabel('sensor'); caxis([-2 2]);


%% Simulate the data, pink noise based on real data
noiseW = 0.1; noiseP = 1; 

% perception data
dataP = samplePinkNoise (real_data.X,real_data.T,nSubs,nTPP,nTrls,nChanInArea,nBrainAreas,noiseW,noiseP);
for sub = 1:nSubs
    for t = 1:nTrls*2
        if t < 101; c = 1; else; c = 2; end
        
        dataP(sub,t,:,:) = squeeze(dataP(sub,t,:,:)) + squeeze(models(:,:,c))';
    end
end
  
% imagery data 
jitM       = 0.8; % normally distributed jitter, mean onset is after 0.8 seconds
jitSD_vals = 0:0.1:0.8; 
njit       = length(jitSD_vals);
accRO      = cell(njit,1);
accSO      = cell(njit,1);
reactSO    = cell(njit,1);
reactRO    = cell(njit,1);
for j = 1:njit % check different SD's
    
    jitSD = jitSD_vals(j);
    
    fprintf('RUNNING JITTER %d OUT OF %d \n',j,njit)
    
    % same order
    dataI_SO = samplePinkNoise (real_data.X,real_data.T,nSubs,nTPI,nTrls,nChanInArea,nBrainAreas,noiseW,noiseP);
    Onsets_SO = zeros(nSubs,nTrls*2);
    
    for sub = 1:nSubs
        for t = 1:nTrls*2
            if t < nTrls+1; c = 1; else; c = 2; end
            
            % add random jitter per trial
            onset = round((jitM + jitSD.*randn(1))/(1/Fs));
            if onset < 1; onset = 1; elseif onset > nTPI-size(models_slow,2); onset = nTPI-size(models_slow,2); end
            Onsets_SO(sub,t) = onset;
            
            dataI_SO(sub,t,onset:onset+(size(models_slow,2))-1,:) = squeeze(dataI_SO(sub,t,onset:onset+(size(models_slow,2))-1,:)) ...
                + squeeze(models_slow(:,:,c))';
        end
    end
    
    % reversed order
    dataI_RO = samplePinkNoise (real_data.X,real_data.T,nSubs,nTPI,nTrls,nChanInArea,nBrainAreas,noiseW,noiseP);
    Onsets_RO = zeros(nSubs,nTrls*2);
    for sub = 1:nSubs
        for t = 1:nTrls*2
            if t < nTrls+1; c = 1; else; c = 2; end
            
            % add random jitter per trial
            onset = round((jitM + jitSD.*randn(1))/(1/Fs));
            if onset < 1; onset = 1; elseif onset > nTPI-size(models_slow,2); onset = nTPI-size(models_slow,2); end
            Onsets_RO(sub,t) = onset;
            
            dataI_RO(sub,t,onset:onset+(size(models_slow,2))-1,:) = squeeze(dataI_RO(sub,t,onset:onset+(size(models_slow,2))-1,:)) ...
                + flip(squeeze(models_slow(:,:,c)),2)';
        end
    end
    
    % decoding
    [accSO{j},reactSO{j}] = TG_crossdecoding(labels,dataP,dataI_SO);
    [accRO{j},reactRO{j}] = TG_crossdecoding(labels,dataP,dataI_RO);
end

save('results','jitSD_vals','reactRO','reactSO','accRO','accSO','-v7.3')


%% Plot for different jitter values

% calculate reactivations based on accuracy (instead of distance per trial)
njit = length(jitSD_vals);
react_accSO = zeros(njit,nSubs,nTPP);
react_accRO = zeros(njit,nSubs,nTPP);
for j = 1:njit    
    accRO{j} = accRO{j};
    accSO{j} = accSO{j};
    reactRO{j} = reactRO{j};
    reactSO{j} = reactSO{j};
    for sub = 1:nSubs        
        SO = squeeze(accSO{j}(sub,:,:));
        [~,react_accSO(j,sub,:)] = max(SO,[],2);
        
        RO = squeeze(accRO{j}(sub,:,:));
        [~,react_accRO(j,sub,:)] = max(RO,[],2);
    end
end


% plot the things
slopes_acc_RO = zeros(njit,3); slopes_acc_SO = zeros(njit,3);
slopes_LDA_RO = zeros(njit,3); slopes_LDA_SO = zeros(njit,3);
c = 1; plotting = 1; 
for j = 2:2:6
   
    figure(1);
    subplot(2,3,c); imagesc(timeI,timeP,squeeze(mean(accSO{j},1)));
    axis xy; ylabel('Perception training time'); xlabel('Imagery testing time');
    caxis([0.45 0.6]); title(sprintf('SO SD: %.2f',jitSD_vals(j)));
    
    subplot(2,3,c+3); imagesc(timeI,timeP,squeeze(mean(accRO{j},1)));
    axis xy; ylabel('Perception training time'); xlabel('Imagery testing time');
    caxis([0.45 0.6]); title(sprintf('RO SD: %.2f',jitSD_vals(j)));
    colormap('hot')
    
    % calculate slopes
    tmp1 = reshape(reactSO{j},size(reactSO{j},1)*size(reactSO{j},2),size(reactSO{j},3));
    tmp2 = reshape(reactRO{j},size(reactRO{j},1)*size(reactRO{j},2),size(reactRO{j},3));
    
    % from LDA
    x = [ones(nSubs*nTrls*2*nTPP,1) repmat(timeP',nSubs*nTrls*2,1)];
    y1 = reshape(squeeze(tmp1.*(1/Fs))',nSubs*nTrls*2*nTPP,1);
    y2 = reshape(squeeze(tmp2.*(1/Fs))',nSubs*nTrls*2*nTPP,1);
    [B1_LDA,Bint1_LDA] = regress(y1,x); [B2_LDA,Bint2_LDA] = regress(y2,x);
    line_SO_LDA = x*B1_LDA; line_RO_LDA = x*B2_LDA;
    slopes_LDA_SO(j,:) = [B1_LDA(2) Bint1_LDA(2,:)];
    slopes_LDA_RO(j,:) = [B2_LDA(2) Bint2_LDA(2,:)];
    
    % accuracy
    x = [ones(nSubs*nTPP,1) repmat(timeP',nSubs,1)];
    y1 = reshape(squeeze(react_accSO(j,:,:).*(1/Fs))',nSubs*nTPP,1);
    y2 = reshape(squeeze(react_accRO(j,:,:).*(1/Fs))',nSubs*nTPP,1);
    [B1_acc,Bint1_acc] = regress(y1,x); [B2_acc,Bint2_acc] = regress(y2,x);
    line_SO_acc = x*B1_acc; line_RO_acc = x*B2_acc;
    slopes_acc_SO(j,:) = [B1_acc(2) Bint1_acc(2,:)];
    slopes_acc_RO(j,:) = [B2_acc(2) Bint2_acc(2,:)];
    
    if plotting
    figure(2);
    subplot(2,3,c); 
    plotCI(tmp1'*(1/Fs),timeP,'CI','b','b','over'); hold on
    plot(timeP,squeeze(mean(mean(reactSO{j},2),1))*(1/Fs),'b','LineWidth',2);
    hold on; plot(timeP,line_SO_LDA(1:nTPP),'k','LineWidth',2);
    xlabel('Perception time (s)'); ylabel('Imagery reactivation time (s)');
    title('Same order'); ylim([0.9 1.15])

    subplot(2,3,c+3); tmp2 = reshape(reactRO{j},size(reactRO{j},1)*size(reactRO{j},2),size(reactRO{j},3));
    plotCI(tmp2'*(1/Fs),timeP,'CI','b','b','over'); hold on
    plot(timeP,squeeze(mean(mean(reactRO{j},2),1))*(1/Fs),'b','LineWidth',2); 
    hold on; plot(timeP,line_RO_LDA(1:nTPP),'k','LineWidth',2);
    xlabel('Perception time (s)'); ylabel('Imagery reactivation time (s)');
    title('Reversed order'); ylim([0.9 1.15])
    
    figure(3);
    subplot(2,3,c)
    plotCI(squeeze(react_accSO(j,:,:))'*(1/Fs),timeP,'CI','b','b','over'); hold on
    plot(timeP,squeeze(mean(react_accSO(j,:,:),2)*(1/Fs)),'b','LineWidth',2); ylim([0.5 2])
    hold on; plot(timeP,line_SO_acc(1:nTPP),'k','LineWidth',2);
    
    subplot(2,3,c+3)
    plotCI(squeeze(react_accRO(j,:,:))'*(1/Fs),timeP,'CI','b','b','over'); hold on
    plot(timeP,squeeze(mean(react_accRO(j,:,:),2)*(1/Fs)),'b','LineWidth',2); ylim([0.5 2])
    hold on; plot(timeP,line_RO_acc(1:nTPP),'k','LineWidth',2);
    
    % plot the slopes
    figure(4); 
    subplot(2,3,c);     
    plot(1,[Bint1_acc(2,1) B1_acc(2) Bint1_acc(2,2)],'bo'); hold on;
    plot(2,[Bint2_acc(2,1) B2_acc(2) Bint2_acc(2,2)],'go'); hold on;
    plot([0 3],[0 0],'k--'); xlim([0 3])
    
    subplot(2,3,c+3); 
    plot(1,[Bint1_LDA(2,1) B1_LDA(2) Bint1_LDA(2,2)],'bo'); hold on;
    plot(2,[Bint2_LDA(2,1) B2_LDA(2) Bint2_LDA(2,2)],'go'); hold on;
    plot([0 3],[0 0],'k--'); xlim([0 3])
    end
    c = c+1;   
    
    
end

legend(cellstr(string(jitSD_vals(1:6))))

