function PlotPerceptionModels(cfg)
% function PlotPerceptionModels(cfg)

%% Load data
nsubjects = length(cfg.subjects);
load(fullfile(cfg.dataDir,cfg.subjects{1},'CleanData',cfg.sensorData),'data');
nTime = length(data.time{1}); clear data

sensorData = zeros(nsubjects,270,nTime);
sourceDataL = zeros(nsubjects,76,nTime);
sourceDataR = zeros(nsubjects,76,nTime);

for sub = 1:nsubjects
    
    fprintf('PROCESSING SUBJECT %s \n',cfg.subjects{sub})
    
    outputDir = fullfile(cfg.outputDir,cfg.subjects{sub},'PerceptionModels');
    if ~exist(outputDir,'dir'); mkdir(outputDir); end
    
    load(fullfile(cfg.dataDir,cfg.subjects{sub},'CleanData',cfg.sensorData),'data');
    cfgS = []; cfgS.channel = 'MEG'; cfgS.keeptrials = 'yes';
    data = ft_timelockanalysis(cfgS,data);
    
    sensorData(sub,:,:) = data.avg;
    
    % source reconstruction
   if ~exist(fullfile(outputDir,'dataFP_atlasRegions.mat'),'file') && cfg.sourceRec
        load(fullfile(cfg.dataDir,cfg.subjects{sub},...
            cfg.sourceFolder,'headmodel'),'vol')
        load(fullfile(cfg.dataDir,cfg.subjects{sub},...
            cfg.sourceFolder,'sourcespace'),'sourcespace')
        vol    = ft_convert_units(vol,'mm');
        
        % forward solution
        cfgS           = [];
        cfgS.grad      = data.grad;
        cfgS.grid.pos  = sourcespace.pos;
        cfgS.grid.inside = 1:size(sourcespace.pos,1);
        if strcmp(cfg.subjects{sub},'S15'); cfgS.channel = {'MEG','-MLT31'};
        else; cfgS.channel   = {'MEG'}; end
        cfgS.headmodel = vol;
        cfgS.backproject = 'no';
        leadfield     = ft_prepare_leadfield(cfgS);
        
        % inverse solution
        cfgS = [];
        cfgS.lambda = 0.01;                 % Regularization; 0.01 corresponds to '1%' in FieldTrip
        cfgS.numPerm = 1000;                 % 100 is fine for testing, but you want more for the real analysis
        cfgS.feedback = 100;                  % provide feedback (soonest) every X seconds
        cfgS.allowNan = 'no';              % Use nanmean, nancov, etc.?
        cfgS.lf = leadfield;
        cfgS.cov.method = 'fixed';
        cfgS.cov.window = 'all';
        cfgS.returnSqrt = 'no';
        cfgS.locationRescale = 'no';
        
        source = sourceAnalysisOneCondition(cfgS, data);
        source.tri = sourcespace.tri;
        
        % correct for positivity and depth bias
        source.avg.pow2 = (source.avg.pow - source.mPerm)./source.mPerm;
        source.avg.pow2(source.avg.pow2 < 0) = 0;
        source.avg.pow2 = sqrt(source.avg.pow2);
        
        % turn into atlas regions
        load(fullfile(outputDir,'dataFP_atlasRegions'),'source')
        load(fullfile(cfg.dataDir,cfg.subjects{sub},cfg.sourceFolder,...
            'dataFP_atlasRegions'),'labIdxL','labIdxR');        
        
        % left-sided areas
        actPerAreaL = zeros(max(labIdxL),size(source.avg.pow,2));
        for i = 1:max(labIdxL) % for every label, get the averaged activity. Exclude 0. This is something I gave to no label assigned, so this is not informative anyway
            % get grid points with a given label
            selGrids = (labIdxL==i);
            selGrids = logical([selGrids; zeros(length(labIdxR), 1)]); % pad, because the first half is the left side.
            actPerAreaL(i,:) = nanmean(source.avg.pow2(selGrids,:),1);
        end
        clear selGrids
        
        % right-sided areas
        actPerAreaR = zeros(max(labIdxL),size(source.avg.pow,2));
        for i = 1:max(labIdxR) % for every label, get the averaged activity.  Exclude 0. This is something I gave to no label assigned, so this is not informative anyway
            % get grid points with a given label
            selGrids = (labIdxR==i);
            selGrids = logical([zeros(length(labIdxL), 1); selGrids]); % pad, because the first half is the left side.
            actPerAreaR(i,:) =  nanmean(source.avg.pow2(selGrids,:),1);
        end
        clear selGrids
        
        % save 
        save(fullfile(outputDir,'dataFP_atlasRegions'),'source','actPerAreaL','actPerAreaR');
        clear source labIdx labIdxL labIdxR nAreas idx 
       
    else
        load(fullfile(outputDir,'dataFP_atlasRegions'),'actPerAreaL','actPerAreaR')
    end
    sourceDataL(sub,:,:) = actPerAreaL;
    sourceDataR(sub,:,:) = actPerAreaR;
    clear data actPerAreaL actPerAreaR
end


%% Plotting

% get some data structures and info
load(fullfile(cfg.dataDir,cfg.subjects{1},'CleanData',cfg.sensorData),'data');
cfgS = []; cfgS.channel = 'MEG'; data = ft_timelockanalysis(cfgS,data);
data.avg = squeeze(mean(sensorData,1));

load(fullfile(cfg.outputDir,cfg.subjects{1},'PerceptionModels','dataFP_atlasRegions'),'source')
load(fullfile(cfg.dataDir,cfg.subjects{1},cfg.sourceFolder,...
            'dataFP_atlasRegions'),'labIdxL','labIdxR')

% for the left side
for a = 1:length(labIdxL)
    if labIdxL(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
        source.avg.pow2(a,:) = NaN;
    else
        source.avg.pow2(a,:) = squeeze(mean(sourceDataL(:,labIdxL(a), :),1));
    end
end

% for the right side
for a = 1:length(labIdxR)
    if labIdxR(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
        source.avg.pow2(a+length(labIdxL),:) = NaN;
    else
        source.avg.pow2(a+length(labIdxL),:) = squeeze(mean(sourceDataR(:,labIdxR(a), :),1));
    end
end
        
bnd.pnt = source.pos;
bnd.tri = source.tri;

% plot the things
nTime = length(cfg.timePoints);
idx   = zeros(nTime,1);
for t = 1:nTime
   [~,idx(t)] = min(abs(cfg.timePoints(t)-data.time)); % find closest match    
end

map = makeColorMaps('redblue');
map = flipud(map);
map2 = map(256/2:end,:);
for t = 1:nTime
    
    % plot the sensor data
    figure(1)
    subplot(2,ceil(nTime/2),t)
    cfgP = [];
    cfgP.layout = 'CTF275.lay';
    cfgP.xlim   = [cfg.timePoints(t) cfg.timePoints(t)];
    cfgP.zlim = [-1.2e-13 1.2e-13];
    cfgP.marker = 'off';
    cfgP.comment = ' ';
    ft_topoplotER(cfgP,data); 
    
    colormap(map); colorbar
    title(string(cfg.timePoints(t)))
    caxis(cfgP.zlim)
    
    % plot the source data
    figure(2)
    T = source.avg.pow2(:,idx(t));
    subplot(2,ceil(nTime/2),t)
    ft_plot_mesh(bnd, 'vertexcolor',T,'edgecolor','none');  
    view([-25 -10]) % right hemisphere
    camlight;
    lighting gouraud;
    material dull;  
    colormap(map2); 
    caxis([2 10]); colorbar
    title(string(cfg.timePoints(t)))
    
end