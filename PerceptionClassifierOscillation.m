function PerceptionClassifierOscillation(cfg)
% function PerceptionClassifierOscillation(cfg)

outputDir = cfg.outputDir;
if ~exist(outputDir,'dir'); mkdir(outputDir); end

%% Define different phases based on extrema
idx = cfg.extrema;
phases = zeros(length(idx)-1,2); % start end for each phase
for e = 1:length(idx)-1
    phases(e,1) = idx(e); phases(e,2) = idx(e+1);
end

%% Get reactivation for all phases
nPhases = length(phases);
S = [];
L = cell(nPhases,nPhases);

%% Run through subjects
nsubjects = length(cfg.subjects);

if ~exist(fullfile(cfg.outputDir,[cfg.outputName '.mat']),'file')
    
    for sub = 1:nsubjects
        
        fprintf('Processing subject %s \n',cfg.subjects{sub});
        
        % get the cross decoding
        load(fullfile(cfg.root,cfg.subjects{sub},cfg.data),'Xhat','trueClass')
        
        % convert distance to correct class
        Xhat(:,:,trueClass==0) = Xhat(:,:,trueClass==0)*-1;
        nTrials  = size(Xhat,3);
        
        for te = 1:nPhases % test on each phases
            xhat = Xhat(:,phases(te,1):phases(te,2),:);
            
            for tr = 1:nPhases % train on each phase
                %if tr == te
                %    L{te,tr} = NaN;
                %else
                xhat2 = xhat(phases(tr,1):phases(tr,2),:,:);
                npoints = size(xhat2,1);
                
                % settings
                react_time = zeros(nTrials,npoints);
                
                % calculate the reactivation time of each time point
                for l = 1:npoints
                    x  = squeeze(xhat2(l,:,:));
                    
                    % find peak
                    for t = 1:nTrials
                        [~,t1] = max(x(:,t));
                        
                        react_time(t,l) = t1; % save it
                    end
                end
                
                
                L{te,tr} = [L{te,tr}; react_time]; % concatenate reactivations
                clear react_time xhat2
                %end
            end
        end
        S = [S; ones(nTrials,1)*sub]; % save subject id
        clear Xhat trueClass
    end
    
    % save results
    save(fullfile(cfg.outputDir,cfg.outputName),'L','S','cfg')
    
else
    load(fullfile(cfg.outputDir,cfg.outputName),'L','S')
end

%% Determine time-line and slope per phase
nTrials = size(L{1},1);
fs   = (cfg.time(end)-cfg.time(1))/length(cfg.time);
slopes = zeros(nPhases,nPhases,2);
X = zeros(nPhases,length(phases(1,1):phases(end,end)),1);
Y = zeros(nPhases,length(phases(1,1):phases(end,end)),nTrials);
for te = 1:nPhases
    
    % time vector for testing phase
    t = cfg.time(phases(te,1):phases(te,2));
    
    for tr = 1:nPhases
        % time vector for training phase
        x = cfg.time(phases(tr,1):phases(tr,2));
        X(te,phases(tr,1):phases(tr,2)) = x;
        
        % reactivation score for phase
        y = L{te,tr};%mean(L{te,tr},1);
        y = (y.*fs) + t(1); % convert to s
        
        Y(te,phases(tr,1):phases(tr,2),:) = y';
        
        % calculate slope
        slopes(te,tr,:) = regress(mean(y,1)',[ones(length(x),1), x']);
        %end
    end
end


% figure; % feedback and feedforward
FW = 1:2:nPhases;
FB = 2:2:nPhases;

% plot reactivation curves per test phase
figure;
for p = 1:nPhases
    idx = phases(p,1):phases(p,2);
    subplot(nPhases,1,p)
    if ismember(p,FW); c = 'b'; else; c = 'r'; end
    plot(X(p,cfg.tStart:end),mean(Y(p,cfg.tStart:end,:),3),c) ;
    hold on; plot(X(p,idx),mean(Y(p,idx,:),3),c,'LineWidth',2)
    xlabel('Training time (s)'); ylabel('Reactivation time');
    title(sprintf('Testing period %.3f to %.3f',cfg.time(phases(p,1)),cfg.time(phases(p,2))))
end

figure;
subplot(2,1,1);
imagesc(squeeze(slopes(:,:,2))); colorbar
title('Slopes'); caxis([-0.2 0.2]); axis xy;
subplot(2,1,2)
normSlopes = zscore(squeeze(slopes(:,:,2)));
imagesc(normSlopes); colorbar
title('Normalized slopes'); caxis([-2 2]); axis xy;

figure;
bar([mean(diag(normSlopes,-1)) mean(diag(normSlopes,0)) mean(diag(normSlopes,1))])
ylabel('Normalizes slopes');
set(gca,'XTickLabels',{'Previous phase','Current phase','Next phase'})
