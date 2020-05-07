function OscillatorReactivation(cfg)

outputDir = fullfile(cfg.root,cfg.output);
if ~exist(outputDir,'dir'); mkdir(outputDir); end

% get the data
load(fullfile(cfg.root,cfg.dataSet),'L');
real_data = L; clear L
load(fullfile(cfg.root,cfg.controlSet),'L');
control_data = L; clear L

%% Wavelet decomposition
sf = 300;
t = -0.5:(1/sf):0.5;
w = 2*pi*10;
s = (1*2*pi)/w;
time = linspace(-0.2,1,360);

% create wavelet
mw = exp((-t.^2)/(2*s^2)+1i*w*t);
mw = mw(50:250);
% plot(real(mw))

%% Evoked oscillation - real data
% throw out last data point
L = real_data(:,1:end-1);
n = size(L,2);
time = time(1:n);

if ~exist(fullfile(outputDir,'EvokedOscillation.mat'),'file')
    
    [mX_Evo,btstrp_Evo,z1_Evo,z2_Evo] = EvokedOscillations(L,mw,1);
    
    w1 = ((n-size(btstrp_Evo,2))/2);
    tp = time(w1+1:end-w1);
    
    % save
    save(fullfile(outputDir,'EvokedOscillation'),'mX_Evo','btstrp_Evo','z1_Evo','z2_Evo','tp')
else
    % load
    load(fullfile(outputDir,'EvokedOscillation'),'mX_Evo','z1_Evo','z2_Evo','tp')
end

% plot the things
figure(1); subplot(2,1,1)
h = fill([tp';flipud(tp')],[z1_Evo';flipud(z2_Evo')],[0 0 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(tp,mX_Evo,'Color','b','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');

%% Evoked oscillation - control data
% throw out last data point
L = control_data(:,1:end-1);
n = size(L,2);
time = time(1:n);

if ~exist(fullfile(outputDir,'EvokedOscillationPerm.mat'),'file')
    
    [mX_Evo,btstrp_Evo,z1_Evo,z2_Evo] = EvokedOscillations(L,mw,1);
    
    w1 = ((n-size(btstrp_Evo,2))/2);
    tp = time(w1+1:end-w1);
    
    % save
    save(fullfile(outputDir,'EvokedOscillationPerm'),'mX_Evo','btstrp_Evo','z1_Evo','z2_Evo','tp')
else
    % load
    load(fullfile(outputDir,'EvokedOscillationPerm'),'mX_Evo','btstrp_Evo','z1_Evo','z2_Evo','tp')
end

% plot the things
figure(1); subplot(2,1,1); hold on
h = fill([tp';flipud(tp')],[z1_Evo';flipud(z2_Evo')],[0.5 0.5 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(tp,mX_Evo,'Color','k','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');

%% Evoked oscillation - compare
L1 = real_data(:,1:end-1);
L2 = control_data(:,1:end-1);
n = size(L1,2);
time = time(1:n);

if ~exist(fullfile(outputDir,'EvokedOscillationCompare.mat'),'file')
    
    [mX_Evo,btstrp_Evo,z1_Evo,z2_Evo] = EvokedOscillationsCompare(L1,L2,mw);
    
    w1 = ((n-size(btstrp_Evo,2))/2);
    tp = time(w1+1:end-w1);
    
    % save
    save(fullfile(outputDir,'EvokedOscillationCompare'),'mX_Evo','btstrp_Evo','z1_Evo','z2_Evo','tp')
else
    % load
    load(fullfile(outputDir,'EvokedOscillationCompare'),'mX_Evo','btstrp_Evo','z1_Evo','z2_Evo','tp')
end

figure(2); subplot(2,1,1); hold on
h = fill([tp';flipud(tp')],[z1_Evo';flipud(z2_Evo')],[0 0.5 0],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(tp,mX_Evo,'Color','g','LineWidth',2);
hold on; plot(xlim,[0 0],'Color','k','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');
xlim([-0.2 1])

% get p values 
pVal = 1-(sum(btstrp_Evo>0)/10000);
[~,qVal] = mafdr(pVal); tmp = ones(length(tp),1)*-25;
hold on; plot(tp(qVal<0.05),tmp(qVal<0.05),'*r','LineWidth',2)


%% Induced oscillation - real data
% throw out last data point
L = real_data(:,1:end-1);
n = size(L,2);
time = time(1:n);

if ~exist(fullfile(outputDir,'InducedOscillation.mat'),'file')
    
    [mX_Ind,btstrp_Ind,z1_Ind,z2_Ind] = InducedOscillations(L,mw);
    
    w1 = ((n-size(btstrp_Ind,2))/2);
    tp = time(w1+1:end-w1);
    
    % save
    save(fullfile(outputDir,'InducedOscillation'),'mX_Ind','btstrp_Ind','z1_Ind','z2_Ind','tp')
    
else
    % load
    load(fullfile(outputDir,'InducedOscillation'),'mX_Ind','btstrp_Ind','z1_Ind','z2_Ind','tp')
end

% plot the things
figure(1); subplot(2,1,2)
h = fill([tp';flipud(tp')],[z1_Ind';flipud(z2_Ind')],[0 0 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(tp,mX_Ind,'Color','b','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Induced oscillation');

%% Induced oscillation - control data
% throw out last data point
L = control_data(:,1:end-1);
n = size(L,2);
time = time(1:n);

if ~exist(fullfile(outputDir,'InducedOscillationPerm.mat'),'file')
    
    [mX_Ind,btstrp_Ind,z1_Ind,z2_Ind] = InducedOscillations(L,mw);
    
    w1 = ((n-size(btstrp_Ind,2))/2);
    tp = time(w1+1:end-w1);
    
    % save
    save(fullfile(outputDir,'InducedOscillationPerm'),'mX_Ind','btstrp_Ind','z1_Ind','z2_Ind','tp')
    
else
    % load
    load(fullfile(outputDir,'InducedOscillationPerm'),'mX_Ind','btstrp_Ind','z1_Ind','z2_Ind','tp')
end

% plot the things
figure(1); subplot(2,1,2); hold on
h = fill([tp';flipud(tp')],[z1_Ind';flipud(z2_Ind')],[0.5 0.5 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(tp,mX_Ind,'Color','k','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Induced oscillation');


%% Frequency power - real data 

% do hanning taper before fft focus on first 400 ms 
idx = 1:241; % -0.2 to 0.6 s 
L = real_data(:,idx); 
L = L-mean(L,2); % remove mean
w = hann(size(L,2)); % make hanning taper
L = L.*w'; % multiply signal by hanning taper

n = size(L,2);
time = linspace(-0.2,1,360);
time = time(idx);

% Evoked power
if ~exist(fullfile(outputDir,'EvokedPower.mat'),'file')
    [mX_power_evo,f,btstrp_power_evo,z1_power_evo,z2_power_evo] = EvokedPower(L,sf,1);
    
    save(fullfile(outputDir,'EvokedPower'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
else
    load(fullfile(outputDir,'EvokedPower'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
end

% plot it
figure(3); subplot(2,1,1); hold on
h = fill([f';flipud(f')],[z1_power_evo';flipud(z2_power_evo')],[0 0 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(f,mX_power_evo,'Color','b','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation'); xlim([2 50])


% induced power
if ~exist(fullfile(outputDir,'InducedPower.mat'),'file')    
    [mX_power_ind,btstrp_power_ind,z1_power_ind,z2_power_ind,f] = InducedPower(L,sf);
    
    save(fullfile(outputDir,'InducedPower'),'mX_power_ind','btstrp_power_ind','z1_power_ind','z2_power_ind','f')
else
    load(fullfile(outputDir,'InducedPower'),'mX_power_ind','btstrp_power_ind','z1_power_ind','z2_power_ind','f')
end

% plot it
figure(3); subplot(2,1,2); 
h = fill([f';flipud(f')],[z1_power_ind';flipud(z2_power_ind')],[0 0 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(f,mX_power_ind,'Color','b','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation'); xlim([2 50])


%% Frequency power - control data 
% do hanning taper before fft focus on first 400 ms 
idx = 1:241; % -0.2 to 0.6 s 
L = control_data(:,idx); 
L = L-mean(L,2); % remove mean
w = hann(size(L,2)); % make hanning taper
L = L.*w'; % multiply signal by hanning taper


n = size(L,2);
time = linspace(-0.2,1,360);
time = time(idx);

% Evoked power
if ~exist(fullfile(outputDir,'EvokedPowerPerm.mat'),'file')
    [mX_power_evo,f,btstrp_power_evo,z1_power_evo,z2_power_evo] = EvokedPower(L,sf,1);
    
    save(fullfile(outputDir,'EvokedPowerPerm'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
else
    load(fullfile(outputDir,'EvokedPowerPerm'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
end

% plot it
figure(3); subplot(2,1,1); hold on
h = fill([f';flipud(f')],[z1_power_evo';flipud(z2_power_evo')],[0.5 0.5 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(f,mX_power_evo,'Color','k','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');
xlim([2 50])


% induced power
if ~exist(fullfile(outputDir,'InducedPowerPerm.mat'),'file')    
    [mX_power_ind,btstrp_power_ind,z1_power_ind,z2_power_ind,f] = InducedPower(L,sf);
    
    save(fullfile(outputDir,'InducedPowerPerm'),'mX_power_ind','btstrp_power_ind','z1_power_ind','z2_power_ind','f')
else
    load(fullfile(outputDir,'InducedPowerPerm'),'mX_power_ind','btstrp_power_ind','z1_power_ind','z2_power_ind','f')
end

% plot it
figure(3); subplot(2,1,2); hold on
h = fill([f';flipud(f')],[z1_power_ind';flipud(z2_power_ind')],[0.5 0.5 0.5],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(f,mX_power_ind,'Color','k','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');
xlim([2 50])

%% Frequency power - compare

% do hanning taper before fft focus on first 400 ms 
idx = 1:241; % -0.2 to 0.6 s 
n   = length(idx);
w   = hann(n); % make hanning taper

L1 = real_data(:,idx); 
L1 = L1-mean(L1,2); % remove mean
L1 = L1.*w'; % multiply signal by hanning taper

L2 = control_data(:,idx); 
L2 = L2-mean(L2,2); % remove mean
L2 = L2.*w'; % multiply signal by hanning taper

time = linspace(-0.2,1,360);
time = time(idx);

% Evoked power
if ~exist(fullfile(outputDir,'EvokedPowerCompare.mat'),'file')
    [mX_power_evo,btstrp_power_evo,z1_power_evo,z2_power_evo,f] = EvokedPowerCompare(L1,L2,sf);
    
    save(fullfile(outputDir,'EvokedPowerCompare'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
else
    load(fullfile(outputDir,'EvokedPowerCompare'),'mX_power_evo','btstrp_power_evo','z1_power_evo','z2_power_evo','f');
end


% plot it
figure(2); subplot(2,1,2); 
h = fill([f';flipud(f')],[z1_power_evo';flipud(z2_power_evo')],[0 0.5 0],'LineStyle','none');
set(h,'facealpha',.5);
hold on; plot(f,mX_power_evo,'Color','g','LineWidth',2);
xlabel('Time (s)'); ylabel('Power');
title('Evoked oscillation');
hold on; plot(xlim,[0 0],'k','LineWidth',2);
xlim([2 50])

% get p values 
pVal = 1-(sum(btstrp_power_evo>0)/10000);
[~,qVal] = mafdr(pVal); tmp = ones(length(f),1)*-0.5;
hold on; plot(f(qVal<0.025),tmp(qVal<0.025),'*r','LineWidth',2)