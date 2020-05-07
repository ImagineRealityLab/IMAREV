function DetermineShape(cfg)

fitshapes = false;

outputDir = cfg.outputDir;
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% load the data 
load(fullfile(cfg.root,cfg.dataSet),'L');
real_data = L; clear L
load(fullfile(cfg.root,cfg.controlSet),'L');
control_data = L; clear L

% settings
sf = 300;
time = linspace(-0.2,1,360);
win  = time >= cfg.window(1) & time <= cfg.window(2);
f    = 11.2;

%% Emperical mode decomposition
data = real_data;

numExtr = 8;
X = squeeze(mean(data(:,win),1))*(1/sf);
[imf,residual,info] = emd(X,'Interpolation','pchip');
[~,pickedMode] = min(abs(info.NumExtrema-numExtr));

% plot results
nModes = size(imf,2);
figure; 
subplot(nModes+2,1,1);
plot(time(win),X,'LineWidth',2); title('Signal');
for m = 1:nModes
    subplot(nModes+2,1,m+1)
    plot(time(win),imf(:,m)); 
    title(sprintf('IMF %d',m))
end
subplot(nModes+2,1,m+2);
plot(time(win),residual)
title('Residual')

figure;
hht(imf,sf); 

figure; tWin = time(win); sig = (X-mean(X))./2;
subplot(3,1,1)
plot(tWin,sig,'b-o','LineWidth',2); hold on
plot(tWin,imf(:,pickedMode),'r','LineWidth',2)
xlm = xlim;

% find maximum and minimum
[iHi,iLo] = findextrema(imf(:,pickedMode));
extrema = sort([1;iHi;iLo]);
hold on; plot(tWin(extrema),imf(extrema,pickedMode),'y*')

% cut into peaks and throughs and check slopes
subplot(3,1,2)
nExtr = length(extrema);
slopes = zeros(nExtr-1,2);
for e = 1:nExtr-1
    idx = extrema(e):extrema(e+1);
    segment = X(idx);
    slopes(e,:) = [ones(1,length(idx));tWin(idx)]'\segment';
    bar(tWin(round(median(idx))),slopes(e,2),0.02); hold on
end
xlim(xlm);

if imf(extrema(1),pickedMode) > imf(extrema(2),pickedMode)
    down = slopes(1:2:nExtr-1,2); up = slopes(2:2:nExtr-1,2);
else
    up = slopes(1:2:nExtr-1,2); down = slopes(2:2:nExtr-1,2);
end
subplot(3,1,3); bar([mean(down) mean(up)]); ylim([-1.2 1.2]);
set(gca,'XtickLabel',{'downwards','upwards'})
downwards_mean = mean(down); upwards_mean = mean(up);

save(fullfile(outputDir,'extrema'),'extrema');

% plot slopes
figure;
plot(tWin,X,'b-o','LineWidth',2); hold on
plot(tWin(extrema),X(extrema),'y*','LineWidth',4); hold on;
for s = 1:nExtr-1
ywin = tWin(extrema(s):extrema(s+1))*slopes(s,2)+slopes(s,1);
plot(tWin(extrema(s):extrema(s+1)),ywin,'r','LineWidth',2); hold on
end
%% Bootstrap
nBootstrap = 10000;
nTrials    = size(data,1);
downwards  = zeros(nBootstrap,1);
upwards    = zeros(nBootstrap,1);
for b = 1:nBootstrap
    
    if mod(b,100) == 0
        fprintf('Bootstrap %d out of %d \n',b,nBootstrap)
    end
    
    tmp = zeros(nTrials,sum(win));
    for t = 1:nTrials
        tmp(t,:) = data(randi(nTrials),win);
    end
    
    X = squeeze(mean(tmp,1))*(1/sf);
    [imf,~,info] = emd(X,'Interpolation','pchip');
    [~,pickedMode] = min(abs(info.NumExtrema-numExtr));
    
    % find maximum and minimum
    [iHi,iLo] = findextrema(imf(:,pickedMode));
    extrema = sort([1;iHi;iLo]);
    
    % cut into peaks and throughs and estimate slopes
    nExtr = length(extrema);
    slopes = zeros(nExtr-1,1);
    for e = 1:nExtr-1
        idx = extrema(e):extrema(e+1);
        segment = X(idx);
        B = [ones(1,length(idx));tWin(idx)]'\segment';
        slopes(e) = B(2);
    end
    
    % divide into downwards and upwards
    if imf(extrema(1),pickedMode) > imf(extrema(2),pickedMode)
        down = slopes(1:2:nExtr-1); up = slopes(2:2:nExtr-1);
    else
        up = slopes(1:2:nExtr-1); down = slopes(2:2:nExtr-1);
    end
    
    downwards(b,1) = mean(down);
    upwards(b,1)   = mean(up); 
    clear up down extrema
    
end

downwards = sort(downwards,'ascend');
CI_down   = [downwards(0.025*nBootstrap) downwards(0.975*nBootstrap)];
upwards = sort(upwards,'ascend');
CI_up   = [upwards(0.025*nBootstrap) upwards(0.975*nBootstrap)];

boxplot([downwards,upwards],'symbol',''); ylim([-3 3])

% absolute difference
absDiff = sort(abs(downwards)-abs(upwards));
CI      = [absDiff(0.025*nBootstrap) absDiff(0.975*nBootstrap)]; % crosses zero so


