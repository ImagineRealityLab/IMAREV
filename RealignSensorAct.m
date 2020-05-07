function RealignSensorAct(cfg)
% function RealignSensorAct(cfg)(cfg)
%
% cfg.dataName = name of the group effect
% cfg.dataDir  = dir of the the data file
% cfg.time     = time point to plot activity of
% cfg.clim     = limits of the color axis

load(fullfile(cfg.dataDir,cfg.subjects{1},cfg.dataFolder,cfg.dataName))
cfgS = []; cfgS.channel = 'MEG'; data = ft_selectdata(cfgS,data);

nChannels = length(data.label); clear data

%% Realign the data
% load the reactivation
load([cfg.reactivations '.mat'],'L','trl','S','X');
trl = trl(2:2:end); % select only imagery trial numbers
nPoints = size(L,2);
time    = cfg.pTime;

for sub = 1:length(cfg.subjects)
    
    outputDir = fullfile(cfg.output,cfg.subjects{sub},cfg.outputDir);
    
    if ~exist(fullfile(outputDir,'dataIM.mat'),'file')
        
        fprintf('Realigning data for subject %s \n',cfg.subjects{sub})
        
        % load the sensor data
        load(fullfile(cfg.dataDir,cfg.subjects{sub},cfg.dataFolder,cfg.dataName))
        cfgS = []; cfgS.channel = 'MEG'; data = ft_selectdata(cfgS,data);
        
        % select the reactivations
        react = L(S==sub,:); %amplitude = X(S==sub,:);
        tmp = find(cfg.tidx); react = react+tmp(1)-1; % make timing compatible
        clear tmp
        
        % select the same trials
        trl_sub = trl{sub};
        nTrials = length(trl_sub);
        sensorData = zeros(nTrials,size(data.trial{1},1),size(data.trial{1},2));
        for t = 1:nTrials
            idx = ismember(data.trialnumbers,trl_sub(t));
            sensorData(t,:,:) = data.trial{idx};
        end
        
        % realign the trials
        realigned_sensorData = zeros(nTrials,size(react,2),size(sensorData,2));
        for t=1:nTrials
            for p=1:size(react,2)
                realigned_sensorData(t,p,:) = sensorData(t,:,react(t,p));%*abs(amplitude(t,p));
            end
        end
        clear sensorData
        
        % put in data struct and safe
        if ~exist(outputDir,'dir'); mkdir(outputDir); end
        data.trialinfo = ones(nTrials,1); data.trialinfo(nTrials/2+1:nTrials) = 2;
        data = rmfield(data,{'sampleinfo','trial','time'}); data.trialnumbers = trl_sub;
        for t = 1:nTrials
            data.trial{t} = squeeze(realigned_sensorData(t,:,:))';
            data.time{t}  = time;
        end
        save(fullfile(outputDir,'dataIM'),'data');
        clear realigned_sensorData data
        
    elseif ~exist(fullfile(outputDir,'dataIM_planar.mat'),'file')
        load(fullfile(outputDir,'dataIM'),'data')
        
        fprintf('Calculating planar for %s \n',cfg.subjects{sub})
        
        % calculate condition diff
        cfgS = []; cfgS.trials = data.trialinfo == 1;
        data_face = ft_timelockanalysis(cfgS,data);
        cfgS.trials = data.trialinfo == 2;
        data_house = ft_timelockanalysis(cfgS,data);
        
        cfgS = []; cfgS.parameter = 'avg'; cfgS.operation = 'x1-x2';
        diff       = ft_math(cfgS,data_face,data_house);
        
        % average data
        cfgS = [];
        data_avg   = ft_timelockanalysis(cfgS,data);
        
        % calculate the planar
        cfgS = []; cfgS.feedback = 'no'; cfgS.method = 'template';
        cfgS.neighbours = ft_prepare_neighbours(cfgS,data);
        cfgS.planarmethod = 'sincos'; cfgS.parameter = 'avg';
        planar = ft_megplanar(cfgS,data_avg);
        planarDiff = ft_megplanar(cfgS,diff);
        
        % combine
        cfgS = []; data = ft_combineplanar(cfgS,planar);
        save(fullfile(outputDir,'dataIM_planar'),'data');
        
        cfgS = []; data = ft_combineplanar(cfgS,planarDiff);
        save(fullfile(outputDir,'dataIM_planar_diff'),'data')
        
        clear data diff data_face data_house planar planarDiff
    end
    
    
end

%% Plot group means
MeanAct  = zeros(length(cfg.subjects),nPoints,nChannels);
kind     = 'planar'; 
for sub = 1:length(cfg.subjects)
    outputDir = fullfile(cfg.output,cfg.subjects{sub},cfg.outputDir);
    if strcmp(kind,'planar')
        load(fullfile(outputDir,'dataIM_planar_diff.mat'),'data');
        MeanAct(sub,:,:) = data.avg';
    elseif strcmp(kind,'axial')
        load(fullfile(outputDir,'dataIM'),'data');
        cfgS = []; cfgS.keeptrials = 'yes';
        data = ft_timelockanalysis(cfgS,data);
        data.avg = squeeze(mean(data.trial(data.trialinfo==1,:,:)))-...
            squeeze(mean(data.trial(data.trialinfo==2,:,:)));
        MeanAct(sub,:,:) = data.avg';
    end
    clear data
end

% plot
diff = squeeze(mean(MeanAct,1))';
N    = size(diff,2);
bin  = [1 7 12 18 22 26];%[1:3:N-3;3:3:N];% % bin some time points together
time = cfg.pTime;
times = cell(length(bin),1);
for t = 1:length(bin); times{t} = int2str(time(bin(t))*1000); end

% put into data struct
load(fullfile(cfg.dataDir,cfg.subjects{1},cfg.dataFolder,cfg.dataName))
cfgS = []; cfgS.channel = 'MEG'; data = ft_selectdata(cfgS,data);
data = rmfield(data,{'trialnumbers','trial','sampleinfo','trialinfo','time'});
data.dimord = 'chan';

figure;
map = makeColorMaps('redblue');
map = flipud(map);
map2 = map(256/2:end,:);
for t = 1:length(bin)
    
    %T = mean(diff(:,bin(1,t):bin(2,t)),2);
    T = diff(:,bin(t));
    data.avg = T;
    
    subplot(2, ceil(length(bin)/2),t)
    cfgP = [];
    cfgP.layout = 'CTF275.lay';
    cfgP.marker      = 'on';
    cfgP.highlight   = 'off';
    cfgP.style       = 'straight';
    if strcmp(kind,'planar')
        cfgP.zlim = [1e-13 2.5e-13];%[-quantile(diff(:),.99) quantile(diff(:),.99)];
    elseif strcmp(kind,'axial')
        cfgP.zlim = [-5e-14 5e-14];
    end
    cfgP.comment = 'no';
    ft_topoplotER(cfgP,data)
    
    if strcmp(kind,'planar')
        colormap(map2);
    elseif strcmp(kind,'axial')
        colormap(map);
    end
    
    title(times{t})
    %title(sprintf('%.0f to %.0f',time(bin(1,t)),time(bin(2,t))))
    %colorbar
    
end