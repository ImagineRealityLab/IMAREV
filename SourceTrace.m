function SourceTrace(cfg)
% function SourceTrace(cfg)

%% Preprocess data
% load the data
load(fullfile(cfg.dataDir,cfg.data),'sourceData','time');
nSubjects = size(sourceData,1);

% select time window

tidx = time >= cfg.tidx(1) & time <= cfg.tidx(2);
time = time(tidx);
sourceData = sourceData(:,:,tidx);

% convert to imagery time for realigned IM data
if strcmp(cfg.data,'realignedIM')
    time = time*-1.171+2.045;
    time = flip(time); % make right order    
end

% get the area names and indices
addpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));
[~,~,colortable] = read_annotation('/vol/ccnlab-scratch1/naddij/ImaMEG/FreesurferOutput/S01/label/lh.aparc.a2009s.annot');
rmpath(genpath('/vol/ccnlab1/naddij/fieldtrip-20170801/external/freesurfer'));

% get activation for specified areas
nAreas     = length(cfg.area);
idx        = cell(nAreas,1);
activation = zeros(nSubjects,nAreas,length(time));
for a = 1:nAreas
    idx{a} = find(ismember(colortable.struct_names,cfg.area{a}));
    
    for sub = 1:nSubjects
        tmp = squeeze(sourceData(sub,idx{a},:));
        if size(tmp,2) > 1
            [~,act] = pca(tmp');
        else
            [~,act] = pca(tmp); % if only 1 dimension, this just zscores
        end
        
        if strcmp(cfg.data,'realignedIM')
            act = flip(act); % make right order
        end
        
        activation(sub,a,:) = act(:,1);
    end 
    
end

%% Calculate cross-correlation
act = squeeze(mean(activation,1));
[cr,lags] =  xcorr(act(1,:),act(2,:));    

fs = abs(time(end)-time(1))/length(time);
lags = lags*fs;

% bootstrap
nBootstrap = 10000;
cr_btstrp  = zeros(nBootstrap,1);
for b = 1:nBootstrap
    tmp = zeros(size(activation));
    for sub = 1:nSubjects
        idx = randi(nSubjects);
        tmp(sub,:,:) = activation(idx,:,:);
    end
    act_b = squeeze(mean(tmp,1));
    [~,cr_btstrp(b)] = max(abs(xcorr(act_b(1,:), act_b(2,:))));
end    

[~,cr_I] = max(abs(cr));

% calculate 95% CI
tmp = sort(abs(cr_btstrp),'ascend');
idx(1) = nBootstrap*0.025; idx(2) = nBootstrap*0.975;
CI = [tmp(idx(1)) tmp(idx(2))];

% plot results
figure;
subplot(2,1,1);
plotCI(squeeze(activation(:,1,:))',time,'SEM','b','b','over');
hold on
plotCI(squeeze(activation(:,2,:))',time,'SEM','r','r','over');
legend('area 1','area 2'); title('Source traces');
xlabel('Time (s)'); ylabel('Stimulus information');

subplot(2,1,2);
plot(lags(cr_I),0,'go')
hold on;
plot(lags(CI),[0 0],'g','LineWidth',2);
xlim([lags(1) lags(end)])
title('Cross correlation');
xlabel('Shift area 2 wrt area 1'); 
fprintf('Lag is %.4f CI from %.4f to %.4f \n',lags(cr_I),lags(CI(1)),lags(CI(2)))


