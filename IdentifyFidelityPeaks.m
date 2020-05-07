function IdentifyFidelityPeaks(cfg)

if ~isfield(cfg,'permute'); cfg.permute = false; end

%% Calculate fidelity peaks
% some initialisation
nsubjects = length(cfg.subjects);
npoints = length(cfg.p_timepoints);
if ~exist(fullfile(cfg.outputDir,[cfg.outputName '.mat']),'file')
    S     = [];
    L     = [];
    X     = [];
    trl   = [];
    tidx = cfg.time >= cfg.timeIdx(1) & cfg.time <= cfg.timeIdx(2);
    
    % load data per subject and calculate reactivation per trial
    for sub = 1:nsubjects
        
        fprintf('Processing subject %s \n',cfg.subjects{sub});
        
        % get the cross decoding
        load(fullfile(cfg.dataDir,cfg.subjects{sub},cfg.dataSet),'Xhat','trueClass','trialnumbers')
        
        % only select after cue and before vividness instr text
        Xhat = Xhat(:,tidx,:);
        nTrials  = size(Xhat,3);
                
        % convert distance to correct class
        Xhat(:,:,trueClass==0) = Xhat(:,:,trueClass==0)*-1; 
                       
        % settings        
        react_time = zeros(nTrials,npoints);
        Xhat_lagged = zeros(nTrials,npoints);
        
        % calculate the reactivation time of each time point
        for l = 1:npoints
            point = cfg.p_timepoints(l);
            xhat  = squeeze(Xhat(point,:,:));   
            
            % low pass filter
            if isfield(cfg,'lpfilter')
                fprintf('Filtering perception point %d out of %d \n',l,npoints);
                xhat(isnan(xhat)) = 0;
                xhat = lowpass(xhat,cfg.lpfilter,300);
            end
            
            % find peak
            for t = 1:nTrials
                [~,t1] = max(xhat(:,t)); 
                
                react_time(t,l) = t1; % save it
                Xhat_lagged(t,l) = Xhat(point,t1,t); % get the distance
            end
        end
        
        S = [S; ones(nTrials,1)*sub]; % save subject id
        L = [L; react_time]; % concatenate reactivations
        X = [X; Xhat_lagged]; % and distances
        %trl = [trl; trialnumbers];       
        clear Xhat trueClass react_time Xhat_lagged trialnumbers
    end
    
    % save results
    if ~exist(cfg.outputDir,'dir'); mkdir(cfg.outputDir); end
    save(fullfile(cfg.outputDir,cfg.outputName),'L','S','X','trl','cfg');
    
else
    load(fullfile(cfg.outputDir,cfg.outputName),'L','S','X')
end

pTime = cfg.pTime(cfg.p_timepoints);

%% Put in long format and add
if cfg.outputMatrix    
    S_long_trial = [];
    for s = 1:nsubjects
        S_long_trial = [S_long_trial;ones(sum(S(:,1)==s)*npoints,1)*s];
    end
    L_long_trial = reshape(L'*(1/300),length(S_long_trial),1);
    Trial = [];
    for sub = 1:nsubjects
        trl = [];
        for t = 1:sum(S==sub)
            tmp = ones(npoints,1)*t;
            trl = [trl; tmp]; clear tmp
        end
        Trial = [Trial;trl]; clear trl
    end
    PTidx = repmat(pTime',length(L),1);
    m     = [S_long_trial, Trial, PTidx, L_long_trial];
    save(fullfile(cfg.outputDir,cfg.outputName),'m','-append')
end

%% Get vividness and append
if cfg.vividness
    V = []; R = [];
    for sub = 1:nsubjects
        load(fullfile(cfg.dataDir,cfg.subjects{sub},cfg.dataSet),'trialnumbers')
        load(fullfile('/vol/ccnlab-scratch1/naddij/ImaMEG/',cfg.subjects{sub},'Behaviour','vividness'),'vividness')
        load(fullfile('/vol/ccnlab-scratch1/naddij/ImaMEG/',cfg.subjects{sub},'Behaviour','responses'),'B')
        vividness = vividness(trialnumbers{2});
        B = B(trialnumbers{2},3);
        R = [R; B]; 
        V = [V;vividness]; clear vividness trialnumbers B
    end
    save(fullfile(cfg.outputDir,cfg.outputName),'V','R','-append')
end

%% Bootstrap
if cfg.bootstrap
   nBootstrap = 10000;
   nTrials    = size(L,1);
   
   btstrp = zeros(nBootstrap,size(L,2));
   for b = 1:nBootstrap
       if mod(b,100) == 0; fprintf('Btstrp %d out of %d \n',b,nBootstrap); end
       tmp = zeros(nTrials,size(L,2));
       for t = 1:nTrials
          idx = randi(nTrials); % pick a random trl
          tmp(t,:) = L(idx,:);
       end
       btstrp(b,:) = mean(tmp,1); % average      
   end
   
   % calculate 95% CI
   tmp = sort(btstrp,'ascend');
   idx(1) = nBootstrap*0.025; idx(2) = nBootstrap*0.975;
   z1 = tmp(idx(2),:);
   z2 = tmp(idx(1),:);
   save(fullfile(cfg.outputDir,cfg.outputName),'z1','z2','btstrp','-append') 
    
end

%% Plot results
if cfg.plotting
    load(fullfile(cfg.outputDir,cfg.outputName),'z1','z2')
    if ~exist('z1','var'); error('First bootstrap to get CI! \n'); end
    
    figure(1); mReact = mean(L*(1/300),1);
    h = fill([pTime';flipud(pTime')],[z1'*(1/300);flipud(z2'*(1/300))],[0.5 0.5 0.5],'LineStyle','none');
    set(h,'facealpha',.5);
    hold on; plot(pTime,mReact,'-o','Color','k','LineWidth',2);   
    xlabel('Perceptual time point');
    ylabel('Reactivation during imagery')
    
    if cfg.vividness
       figure(2);
       [r,p] = corr(X,V);
       plot(pTime,r); hold on; plot(pTime(p<0.01),r(p<0.01),'*')
    end
end
