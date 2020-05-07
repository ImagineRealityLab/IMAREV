function LinkOscillatorReactivationRawData(cfg)


%% get the reactivation data
load(fullfile(cfg.output,cfg.dataSet),'L','S');
real_data = L; clear L
load(fullfile(cfg.output,cfg.controlSet),'L');
control_data = L; clear L

%% Create hanning taper
time = linspace(-0.2,1,360);
tidx  = 1:241; % -0.2 to 0.6 - same as freq analysis
n = length(tidx);
w = hann(n); % make hanning taper
time = time(tidx);
sf = round(1/((time(end)-time(1))/n));
f = sf*(0:n-1)/n;

%% Loop over subjects to get FFTs
nsubjects = length(cfg.subjects);
fftReac = zeros(nsubjects,n); fftRaw = zeros(nsubjects,n,270);
fftCntr = zeros(nsubjects,n);
for sub = 1:nsubjects
    
    
    fprintf('Processing subject %s \n',cfg.subjects{sub})
    
    %% Preprocess raw data
    load(fullfile(cfg.output,cfg.subjects{sub},cfg.decodingData),'trialnumbers');
    load(fullfile(cfg.root,cfg.subjects{sub},cfg.rawData),'data')
    
    % time lock
    cfgT = []; cfgT.channel = 'MEG';
    cfgT.trials = find(ismember(data.trialnumbers,trialnumbers{1}));
    data = ft_timelockanalysis(cfgT,data);
    
    % convert to planar
    cfgS = []; cfgS.feedback = 'no'; cfgS.method = 'template';
    cfgS.neighbours = ft_prepare_neighbours(cfgS,data);
    cfgS.planarmethod = 'sincos';
    planar = ft_megplanar(cfgS,data);
    cfgS = []; data = ft_combineplanar(cfgS,planar);
    
    % demean and hanning
    rawData = data.avg(:,tidx)-mean(data.avg(:,tidx),2);
    rawData = rawData'.*w;
    
    %% Preprocess reactivations
    % real data
    rlDat = mean(real_data(S==sub,tidx),1);
    rlDat = rlDat-mean(rlDat);
    rlDat = rlDat.*w';
    
    % control data
    cnDat = mean(control_data(S==sub,tidx),1);
    cnDat = cnDat-mean(cnDat);
    cnDat = cnDat.*w';
    
    %% Get FFT
    for c = 1:size(rawData,2) % loop over channels
        fftRaw(sub,:,c) = fft(rawData(:,c)); % raw data
    end
    
    fftReac(sub,:) = fft(rlDat);
    fftCntr(sub,:) = fft(cnDat);
    
end

%% Calculate PSD and CSD
nChannels = size(rawData,2);
nsubjects = size(fftRaw,1);

% for the real data
coh_real = zeros(n,nChannels);
angle_diff_real = zeros(n,nChannels);
csd_real = zeros(n,nsubjects,nChannels);
for c = 1:nChannels
    signals = cat(3,fftReac',squeeze(fftRaw(:,:,c))');
    psd = 2.*abs(signals).^2./(nsubjects.^2);
    csd_real(:,:,c) = 2.*(signals(:,:,1).*conj(signals(:,:,2)))./(nsubjects.^2);
    sumpsd = squeeze(sum(psd,2));
    sumcsd = squeeze(sum(csd_real(:,:,c),2));
    coh_real(:,c) = abs(sumcsd ./ sqrt(sumpsd(:,1) .* sumpsd(:,2)));
    angle_diff_real(:,c) = angle(sumcsd(:,1))*(180/pi);
end

% for the control data
coh_cntr = zeros(n,nChannels);
angle_diff_cntr = zeros(n,nChannels);
csd_cntr = zeros(n,nsubjects,nChannels);
for c = 1:nChannels
    signals = cat(3,fftCntr',squeeze(fftRaw(:,:,c))');
    psd = 2.*abs(signals).^2./(nsubjects.^2);
    csd_cntr(:,:,c) = 2.*(signals(:,:,1).*conj(signals(:,:,2)))./(nsubjects.^2);
    sumpsd = squeeze(sum(psd,2));
    sumcsd = squeeze(sum(csd_cntr(:,:,c),2));
    coh_cntr(:,c) = abs(sumcsd ./ sqrt(sumpsd(:,1) .* sumpsd(:,2)));
    angle_diff_cntr(:,c) = angle(sumcsd(:,1))*(180/pi);
end

%% Plot the results
data = rmfield(data,{'time'});
data.dimord = 'chan';

cfgP = [];
cfgP.layout = 'CTF275.lay';
cfgP.parameter = 'avg';
cfgP.zlim  = [0 0.5];%'maxmin';
cfgP.marker = 'on';
cfgP.highlight = 'off';
cfgP.comment = 'no';
cfgP.colorbar = 'no';

foi = 11.2;
[~,fidx] = min(abs(f-foi));

figure(1);
% coherence
map = makeColorMaps('dusk');
map = flipud(map);
subplot(2,2,1);
data.avg = coh_real(fidx,:)';
ft_topoplotER(cfgP,data);
title('coherence - real')

subplot(2,2,2);
data.avg = coh_cntr(fidx,:)';
ft_topoplotER(cfgP,data);
title('coherence - control')

subplot(2,2,3);
data.avg = coh_real(fidx,:)'-coh_cntr(fidx,:)';
cfgP.zlim  = [-0.3 0.3];
%cfgP.colorbar = 'west';
ft_topoplotER(cfgP,data);
title('real - control');


% angle
figure;
subplot(1,2,1);
data.avg = angle_diff_real(fidx,:)';
data.dimord = 'chan';
ft_topoplotER(cfgP,data);
title('angle difference - real')
colormap(hsv); %caxis([-180 180])

subplot(1,2,2);
data.avg = angle_diff_cntr(fidx,:)';
data.dimord = 'chan';
ft_topoplotER(cfgP,data);
colormap(hsv); %caxis([-180 180])
title('angle difference - control')

%% Do stats
% compare real and shuffled csd
nPerm = 10000;
permutations = zeros(nChannels,nPerm);
true_distr   = zeros(nChannels,1);
for p = 1:nPerm
    perm = randperm(nsubjects*2); % random permutation
    if mod(p,100)==0
        fprintf('Permutation %d out of %d \n',p,nPerm)
    end
    for c = 1:nChannels
        real = abs(squeeze(csd_real(fidx,:,c)));
        cntr = abs(squeeze(csd_cntr(fidx,:,c)));
        
        true_distr(c,1) = mean(real-cntr);
        
        dat = [real,cntr]; dat = dat(perm);
        permutations(c,p) = mean(dat(1:nsubjects)-dat(nsubjects+1:nsubjects*2));
    end
end

pVals        = sum(permutations >= true_distr,2)./nPerm;

figure(1);
subplot(2,2,4);
data.avg = log10(pVals);
cfgP.zlim = [-2 0];%'maxmin';
cfgP.colorbar = 'east';
ft_topoplotER(cfgP,data);
colormap(map);
