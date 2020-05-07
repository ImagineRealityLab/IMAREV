function RealignSourceReconstruction(cfg0)
% function RealignSourceReconstruction(cfg0)

%% Sourcereconstruction for all subjects
% save name
if strcmp(cfg0.kind,'activation')
    saveName = 'dataIM_sourceAct.mat';
elseif strcmp(cfg0.kind,'condifference')
    saveName = 'dataIM_sourceConDiff.mat';
end
load(fullfile(cfg0.dataDir,cfg0.subjects{1},cfg0.dataName),'data');
nPoints = size(data.trial{1},2); clear data

nsubjects = length(cfg0.subjects);
sourceDataL = zeros(nsubjects,76,nPoints);
sourceDataR = zeros(nsubjects,76,nPoints);

for sub = 1:nsubjects
    
    
    % save dir
    saveDir = fullfile(cfg0.dataDir,cfg0.subjects{sub},...
        fileparts(cfg0.dataName));
    
    if ~exist(fullfile(saveDir,saveName),'file')
        
        fprintf('SOURCERECONSTRUCTING SUBJECT %s \n',cfg0.subjects{sub})
        
        % get the sensor data
        load(fullfile(cfg0.dataDir,cfg0.subjects{sub},cfg0.dataName),'data')
        
        % load the source information
        load(fullfile(cfg0.sourceRecDir,cfg0.subjects{sub},...
            cfg0.sourceRecFolder,'headmodel'),'vol')
        load(fullfile(cfg0.sourceRecDir,cfg0.subjects{sub},...
            cfg0.sourceRecFolder,'sourcespace'),'sourcespace')
        vol    = ft_convert_units(vol,'mm');
        
        % timelock on sensor data to get covariance
        %             for t = 1:length(data.trialinfo)
        %                 data.time{t} = linspace(0,0.0903,28); % correct time axis
        %             end
        cfg  = [];
        cfg.channel    = {'MEG'};
        cfg.keeptrials = 'yes';
        tlckAll        = ft_timelockanalysis(cfg,data);
        
        % forward solution
        cfg           = [];
        cfg.grad      = tlckAll.grad;
        cfg.grid.pos  = sourcespace.pos;
        cfg.grid.inside = 1:size(sourcespace.pos,1);
        if strcmp(cfg0.subjects{sub},'S15'); cfg.channel = {'MEG','-MLT31'};
        else; cfg.channel   = {'MEG'}; end
        cfg.headmodel = vol;
        cfg.backproject = 'no';
        leadfield     = ft_prepare_leadfield(cfg);
        
        % inverse solution
        cfg = [];
        cfg.lambda = 0.01;                 % Regularization; 0.01 corresponds to '1%' in FieldTrip
        cfg.numPerm = 1000;                 % 100 is fine for testing, but you want more for the real analysis
        cfg.feedback = 5;                  % provide feedback (soonest) every X seconds
        cfg.allowNan = 'no';              % Use nanmean, nancov, etc.?
        cfg.lf = leadfield;
        cfg.cov.method = 'fixed';
        cfg.cov.window = 'all';
        cfg.returnSqrt = 'no';
        cfg.locationRescale = 'no';
        
        % reconstruct mean activation or difference between conditions
        if strcmp(cfg0.kind,'condifference')
            design = data.trialinfo == 1;
            source = sourceAnalysisTwoConditions(cfg, tlckAll, design);
        elseif strcmp(cfg0.kind,'activation')
            source = sourceAnalysisOneCondition(cfg, tlckAll);
        end
        
        source.tri = sourcespace.tri;
        
        % correct for positivity and depth bias
        source.avg.pow2 = (source.avg.pow - source.mPerm)./source.mPerm;
        source.avg.pow2(source.avg.pow2 < 0) = 0;
        source.avg.pow2 = sqrt(source.avg.pow2);
        
        % save source reconstruction
        save(fullfile(saveDir,saveName),'source')
        
        % plot
        if isfield(cfg0,'plot')
            bnd.pnt = source.pos;
            bnd.tri = source.tri;
            m = source.avg.pow2(:,end);
            
            figure; ft_plot_mesh(bnd,'vertexcolor',m); colorbar;
            title(cfg0.dataName); view([-45 0])
            pause(3); close
        end
        
        clear data vol sourcespace tlckAll leadfield bnd
    else
        load(str2fullfile(saveDir,saveName),'source')
    end
    
    if ~exist(fullfile(saveDir,'dataIM_atlas_sourceConDiff.mat'),'file')
        fprintf('PUTTING SOURCE DATA IN ATLAS FOR %s \n',cfg0.subjects{sub})
        
        % turn into atlas regions
        load(fullfile(cfg0.sourceRecDir,cfg0.subjects{sub},cfg0.sourceRecFolder,...
            'dataIM_atlasRegions'),'labIdxL','labIdxR');
        
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
        
        % save atlas
        if strcmp(cfg0.kind,'activation')
            save(fullfile(saveDir,'dataIM_atlas_sourceAct'),'actPerAreaL','actPerAreaR')
        elseif strcmp(cfg0.kind,'condifference')
            save(fullfile(saveDir,'dataIM_atlas_sourceConDiff'),'actPerAreaL','actPerAreaR')
        end
    else
        load(fullfile(saveDir,'dataIM_atlas_sourceConDiff'),'actPerAreaL','actPerAreaR')
        
        sourceDataL(sub,:,:) = actPerAreaL;
        sourceDataR(sub,:,:) = actPerAreaR;
        clear sourceAtlas labIdx labIdxL labIdxR nAreas idx act
    end
    
end


%% Average and plot over participants
sourceDataL(cfg0.nData==1,:,:) = [];
sourceDataR(cfg0.nData==1,:,:) = [];

% get some plotting info
load(fullfile(cfg0.dataDir,cfg0.subjects{1},fileparts(cfg0.dataName),...
    saveName),'source');
load(fullfile(cfg0.sourceRecDir,cfg0.subjects{1},cfg0.sourceRecFolder,...
    'dataIM_atlasRegions'),'labIdxL','labIdxR');

% put them in a source structure
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

% plot the stuff
bnd.pnt = source.pos; %bnd.pnt = bnd.pnt(1:length(bnd.pnt)/2,:);
bnd.tri = source.tri; %bnd.tri = bnd.tri(1:length(bnd.tri)/2,:);

% define plot points
time = source.time;
nPlots  = length(cfg0.plot);
[~,bin] = min(abs(time'-cfg0.plot));

times = cell(nPlots,1);
for t = 1:length(bin); times{t} = int2str(time(bin(t))*1000); end

% plot them
figure;
for t = 1:nPlots
    
    T = source.avg.pow2(:,bin(t)); T = T(1:length(T)/2);
    
    subplot(2, ceil(length(bin)/2),t)
    ft_plot_mesh(bnd, 'vertexcolor',T,'edgecolor','none');
    %view([-92 -40]) % ventral view
    %view([-90 90]); % dorsal view
    view([-25 -10]) % right hemisphere
    %view([-150 -10]) % left hemishere
    %view([-180 15]) % left frontal
    camlight;
    lighting gouraud;
    material dull;
    title(times{t})
    map = makeColorMaps('redblue');
    map = flipud(map);
    map = map(256/2:end,:);
    colormap(map);
    caxis([1 6]);
    %colorbar
    
end


%% Plot specific atlas areas

% get atlas parcel names
addpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));
[~,~,colortable] = read_annotation('/vol/ccnlab-scratch1/naddij/ImaMEG/FreesurferOutput/S01/label/lh.aparc.a2009s.annot');

areas{1} = 'Pole_occipital';
areas{2} = {'G_oc-temp_lat-fusifor','S_oc-temp_lat',...
    'S_oc-temp_med_and_Lingual'};
nAreas = length(areas); idx = cell(nAreas,1);
figure; nSources = zeros(nAreas,1);
map = makeColorMaps('redblue');
for a = 1:nAreas
    idx{a} = find(ismember(colortable.struct_names,areas{a}));
    T = [labIdxL;labIdxR]; T(~ismember(T,idx{a})) = 0;
    if a == 1
        T(ismember(T,idx{a})) = 1;
    elseif a == 2
        T(ismember(T,idx{a})) = -1;
    end
    nSources(a) = sum(T);
    
    %T = T(1:length(T)/2);
    
    subplot(2,2,a)
    ft_plot_mesh(bnd, 'vertexcolor',T,'edgecolor','none');
    view([-25 -10]) % right hemisp
    camlight;
    lighting gouraud;
    material dull;
    colormap(map); caxis([-1 1])
    title(areas{a})
end

subplot(2,2,3:4);
sourceData = (sourceDataL+sourceDataR)./2;
tidx = 1:22;
load(cfg0.reactivations,'L','S'); L = mean(L(:,tidx),1)*(1/300);
imaTime = time(tidx)*-1.171+2.045;

activation = zeros(nsubjects,nAreas,length(tidx));%time));
colors = ['b','r'];
for a = 1:nAreas
    
    for sub = 1:nsubjects
        tmp = squeeze(sourceData(sub,idx{a},tidx));
        if size(tmp,2) > 1
            [~,act] = pca(tmp');
        else
            [~,act] = pca(tmp);
        end
        activation(sub,a,:) = act(:,1);
        %activation(sub,a,:) = mean(tmp(sub,:,:),2);
    end
    
    %[~,act] = pca(source.avg.pow2(ind,:)');
    
    plotCI(squeeze(activation(:,a,tidx))',imaTime,'SEM',colors(a),colors(a),'over');  hold on
    plot(imaTime,squeeze(mean(activation(:,a,tidx),1)),colors(a),'LineWidth',2); hold on
    [~,b] = max(activation(a,tidx));
end
legend('area1','area2')


save(fullfile(cfg0.dataDir,'GroupResults','SourceActivation','realignedIM'),'sourceData','time')

%% Compare realigned to aligned activation

real_act = mean(sourceData(:,idx{1},tidx),3); % average over time
real_act(:,2) = mean(mean(sourceData(:,idx{2},tidx),3),2);

% get un-realigned
load('/vol/ccnlab-scratch1/naddij/ImaMEG/GroupResults/SourceReconstruction/dataIM_atlasRegions',...
    'actLeft','actRight')
actunreal = (actLeft+actRight)./2;
sourceData = permute(actunreal,[3,1,2]); time = linspace(-0.2,4.5,1410);
save(fullfile(cfg0.dataDir,'GroupResults','SourceActivation','IM'),'sourceData','time')


[~,b] = max(squeeze(nanmean(nanmean(actunreal,1),3)));
tidx2 = b-length(tidx)/2:b+length(tidx)/2-1;

unreal_act = squeeze(mean(actunreal(idx{1},tidx2,:),2));
unreal_act(:,2) = mean(mean(actunreal(idx{2},tidx2,:),1),2);

act = [unreal_act, real_act];
err = std(act)./sqrt(nsubjects);

figure;
subplot(2,2,1:2)
bar(mean(act)); set(gca,'XTickLabels',{'EVC-ur','IT-ur','EVC-r','IT-r'})
for b = 1:length(err) % add SEM
    hold on;
    plot([b b],[mean(act(:,b)) mean(act(:,b))+err(b)],'k','LineWidth',2)
end


% put in a source structure
source_real = source; source_unreal = source;
source_real.avg.pow2 = []; source_unreal.avg.pow2 = [];
for a = 1:length(labIdxL)
    if labIdxL(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
        source_real.avg.pow2(a,:) = NaN;
        source_unreal.avg.pow2(a,:) = NaN;
    else
        source_real.avg.pow2(a,:) = squeeze(mean(mean(sourceDataL(:,labIdxL(a), tidx),3),1));
        source_unreal.avg.pow2(a,:) = squeeze(mean(mean(actLeft(labIdxL(a), tidx2,:),2),3));
    end
end

% for the right side
for a = 1:length(labIdxL)
    if labIdxL(a) == 0 % this seems to happen in desikan. In that case, give no value to this grid point, as it has no defined label.
        source_real.avg.pow2(a+length(labIdxL),:) = NaN;
        source_unreal.avg.pow2(a+length(labIdxL),:) = NaN;
    else
        source_real.avg.pow2(a+length(labIdxL),:) = squeeze(mean(mean(sourceDataR(:,labIdxR(a), tidx),3),1));
        source_unreal.avg.pow2(a+length(labIdxL),:) = squeeze(mean(mean(actRight(labIdxR(a), tidx2,:),2),3));
    end
end

bnd.pnt = source.pos; %bnd.pnt = bnd.pnt(1:length(bnd.pnt)/2,:);
bnd.tri = source.tri; %bnd.tri = bnd.tri(1:length(bnd.tri)/2,:);

map = makeColorMaps('redblue');
map = flipud(map);
map = map(256/2:end,:);

subplot(2,2,3)
T = source_unreal.avg.pow2;
ft_plot_mesh(bnd, 'vertexcolor',T,'edgecolor','none');
view([-25 -10]) % right hemisp
camlight;
lighting gouraud;
material dull;
colormap(map); caxis([0 1]);colorbar
title('unrealigned')

subplot(2,2,4)
T = source_real.avg.pow2;
ft_plot_mesh(bnd, 'vertexcolor',T,'edgecolor','none');
view([-25 -10]) % right hemisp
camlight;
lighting gouraud;
material dull;
colormap(map); caxis([2 5]); colorbar
title('realigned')

% also for perception
load('/vol/ccnlab-scratch1/naddij/ImaMEG/GroupResults/SourceReconstruction/dataFP_atlasRegions',...
    'actLeft','actRight')
actPerc = (actLeft+actRight)./2; tidx3 = 61:181;%82:103;
actPerception = zeros(nsubjects,nAreas,length(tidx3));
t = linspace(-0.2,1,360);
pTime = t(tidx3);
for a = 1:nAreas
    for sub = 1:nsubjects
        tmp = squeeze(actPerc(idx{a},tidx3,sub));
        if size(tmp,2) > 1
            [~,act] = pca(tmp');
        else
            [~,act] = pca(tmp);
        end
        actPerception(sub,a,:) = act(:,1);    %
    end
    
    plotCI(squeeze(actPerception(:,a,:))',pTime,'SEM',colors(a),colors(a),'over');  hold on
    plot(pTime,squeeze(mean(actPerception(:,a,:),1)),colors(a),'LineWidth',2); hold on
    
end

sourceData = permute(actPerc,[3,1,2]); time = linspace(-0.2,1,360);
save(fullfile(cfg0.dataDir,'GroupResults','SourceActivation','FP'),'sourceData','time')
