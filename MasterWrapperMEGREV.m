% Check Linde-Domingo et al., Nature Communications (2019) paper 
restoredefaultpath;

% add paths
addpath(genpath('/vol/ccnlab1/naddij/Analyses/MEGREV/FinalCode'))
%addpath(genpath('/vol/ccnlab1/naddij/Analyses/ImaMEG'))
%addpath(genpath('/vol/ccnlab1/naddij/Analyses/Utilities'))
%addpath(genpath('/vol/ccnlab1/naddij/Analyses/Decoding'))

addpath('/vol/ccnlab1/naddij/fieldtrip-20170801');
addpath('/vol/ccnlab1/naddij/fieldtrip-20170801/qsub');
%addpath('/vol/optdcc/fieldtrip-latest/fieldtrip/');
%addpath('/vol/optdcc/fieldtrip-latest/fieldtrip/qsub');

% set defaultsc
ft_defaults

% directories
root = '/vol/ccnlab-scratch1/naddij/ImaMEG'; % datadir - download data from http://hdl.handle.net/11633/di.dcc.DSC_2017.00072_245
output = '/vol/ccnlab-scratch1/naddij/MEGREV'; % where results will be stored 
if ~exist(output,'dir'); mkdir(output); end
cd(output)

% subjects
subjects = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S15','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'};
nsubjects = length(subjects);

%% Plot perception models
cfg = [];
cfg.timePoints = [0.07 0.09 0.11 0.13];
cfg.sensorData = 'dataFP';
cfg.sourceRec  = 1;
cfg.dataDir    = root;
cfg.sourceFolder = 'SourceReconstruction';
cfg.subjects   = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'}; % skip S15 because they had a bad sensor 
cfg.outputDir  = output;

PlotPerceptionModels(cfg)

%% Cross-decoding perception -> imagery
for sub = 1:nsubjects
cfg = [];
cfg.dataFiles{1}       = fullfile(root,subjects{sub},'CleanData','dataFP');
cfg.dataFiles{2}       = fullfile(root,subjects{sub},'CleanData','dataIM');
cfg.outputDir          = fullfile(output,subjects{sub},'cross_decoding');
cfg.nFold              = 5;


cfg.conIdx{1,1}        = [101:108,201:208];
cfg.conIdx{1,2}        = [109:116,209:216]; 
cfg.conIdx{2,1}        = [121:128,221:228];
cfg.conIdx{2,2}        = [129:136,229:236];

cfg.channels           = 'MEG';
cfg.demean             = true;
cfg.zscore             = true;

cfg.decfg.gamma        = 0.05;
cfg.decfg.samplemethod = 'downsample';
cfg.decfg.nMeanS       = 0; % average over x samplescfg.Permutation = 'true' or 'false' whether to permute labels, default is 'false'

DecodingCross(cfg)

end

%% Get the reactivations
cfg = [];
cfg.subjects     = {'S01'};%subjects;
cfg.p_timepoints = 82:100; % 70 to 130 ms
cfg.pTime        = linspace(-0.2,1,360); % perception time vector
cfg.dataSet      = fullfile('cross_decoding','dataFP to dataIM');
cfg.dataDir      = output;
cfg.outputDir    = fullfile(output,'GroupResults','reactivation_times');
cfg.outputName   = 'Reactivation_timing';
cfg.time         = linspace(-0.2,4.5,1410); % imagery time vector
cfg.timeIdx      = [0 4]; % look only at imagery time period cue - vividness rating
cfg.outputMatrix = true;
cfg.plotting     = true;
cfg.bootstrap    = true;

IdentifyFidelityPeaks(cfg);


%% Get reactivations for entire perception period
cfg = [];
cfg.subjects     = subjects;
cfg.p_timepoints = 1:360;
cfg.pTime        = linspace(-0.2,1,360);
cfg.dataSet      = fullfile('cross_decoding','dataFP to dataIM');
cfg.dataDir      = output;
cfg.outputDir    = fullfile(output,'GroupResults','reactivation_times');
cfg.outputName   = 'Reactivation_timing_1to360';
cfg.time         = linspace(-0.2,4.5,1410);
cfg.timeIdx      = [0 4];
cfg.outputMatrix = true;
cfg.bootstrap    = true;
cfg.plotting     = true;

IdentifyFidelityPeaks(cfg);

%% Train on random classifier
for sub = 1:nsubjects
cfg = [];
cfg.dataFiles{1}       = fullfile(root,subjects{sub},'CleanData','dataFP');
cfg.dataFiles{2}       = fullfile(root,subjects{sub},'CleanData','dataIM');
cfg.outputDir          = fullfile(output,subjects{sub},'ShuffledPerception');
cfg.nFold              = 5;

cfg.conIdx{1,1}        = [101:108,201:208];
cfg.conIdx{1,2}        = [109:116,209:216]; 
cfg.conIdx{2,1}        = [121:128,221:228];
cfg.conIdx{2,2}        = [129:136,229:236];

cfg.channels           = 'MEG';
cfg.appName            = 'perm';

cfg.decfg.gamma        = 0.05;
cfg.decfg.samplemethod = 'downsample';
cfg.decfg.nMeanS       = 0; % average over x samplescfg.Permutation     = 'true' or 'false' whether to permute labels, default is 'false'

cfg.Permutation        = true; % shuffle training labels
cfg.nPermutations      = 1;

DecodingCross(cfg)
end 


%% Get reactivations of random classifier
cfg = [];
cfg.subjects     = subjects;
cfg.p_timepoints = 82:100;
cfg.pTime        = linspace(-0.2,1,360);
cfg.dataSet      = fullfile('ShuffledPerception','Permuted(1)dataFP to dataIM_perm');
cfg.dataDir      = output;
cfg.outputDir    = fullfile(output,'GroupResults','reactivation_times');
cfg.outputName   = 'Reactivation_timing_perm';
cfg.time         = linspace(-0.2,4.5,1410);
cfg.timeIdx      = [0 4];
cfg.outputMatrix = true;
cfg.plotting     = true;
cfg.bootstrap    = true;

IdentifyFidelityPeaks(cfg);

% entire period
cfg.p_timepoints = 1:360;
cfg.outputName   = 'Reactivation_timing_perm_1to360';
IdentifyFidelityPeaks(cfg);


%% Realign data on reactivation times

% realign sensor data
cfg            = [];
cfg.subjects   = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'}; % skip S15 because they had a bad sensor 
cfg.dataName   = 'dataIM';
cfg.dataDir    = root;
cfg.dataFolder = '/CleanData';
cfg.reactivations = fullfile(output,'GroupResults','reactivation_times','Reactivation_timing_1to360');
cfg.time       = linspace(-0.2,4.5,1410);
cfg.tidx       = cfg.time >= 0 & cfg.time <= 4;
cfg.output     = output;
cfg.outputDir  = 'RealignedData_entireperiod';
cfg.pTime      = linspace(0.07,0.8,219);

RealignSensorAct(cfg)

% source reconstruction of realigned data
cfg            = [];
cfg.reactivations = fullfile(output,'GroupResults','ReactivationTimes','Reactivation_timing_1to360');
cfg.subjects   = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'}; % skip S15 because they had a bad sensor 
cfg.dataName   = 'RealignedData_entireperiod/dataIM';
cfg.dataDir    = output;
cfg.sourceRecDir = root;
cfg.sourceRecFolder = 'SourceReconstruction';
cfg.kind       = 'condifference'; % activation or condifference
cfg.plot       = [0.15 0.14 0.13 0.12 0.11 0.1 0.09 0.08];

RealignSourceReconstruction(cfg)

%% Look at alpha in reactivations
% calculate evoked and induced power from real and control data set
cfg = [];
cfg.root = output;
cfg.dataSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_1to360');
cfg.controlSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_perm_1to360');
cfg.output = fullfile('GroupResults','OscillationReactivation');

OscillatorReactivation(cfg);

% find the extrema of the oscillation
cfg = [];
cfg.root = output;
cfg.dataSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_1to360');
cfg.controlSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_perm_1to360');
cfg.plot = true;
cfg.window = [0.07 0.47]; % only when there is alpha
cfg.outputDir = fullfile(output,'GroupResults','EMD_nofilter');

DetermineShape(cfg)


% Link reactivation alpha to alpha in raw data
cfg = [];
cfg.subjects   = {'S01','S03','S04','S06','S07','S08','S09','S11','S14','S16','S17','S18','S19','S20','S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'}; % skip S15 because they had a bad sensor
cfg.root         = root;
cfg.output       = output;
cfg.decodingData = fullfile('cross_decoding','dataFP to dataFP.mat');
cfg.dataSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_1to360');
cfg.controlSet = fullfile('GroupResults','reactivation_times',...
    'Reactivation_timing_cc_perm_1to360');
cfg.rawData      = fullfile('CleanData','dataFP.mat');
LinkOscillatorReactivationRawData(cfg)

% Do a similar analysis for perception
cfg = [];
cfg.root    = output; 
cfg.data    = 'CrossDecoding/dataFP to dataFP';
cfg.subjects = subjects;
cfg.outputDir = fullfile(output,'GroupResults','PerceptionOscillation');
cfg.tStart  = 1;
foi = 10; lphase = (1/foi)/2; tphase = lphase/(1/300); % for 10 hz
%load(fullfile(output,'GroupResults','EMD','extrema')); % use extrema defined in earlier analysis
cfg.extrema = cfg.tStart:tphase:360; %extrema+cfg.tStart-1;%cfg.tStart:tphase:200;%cfg.tStart:tphase:241;%
cfg.time    = linspace(-0.2,1,360);
cfg.outputName = 'PercReactPhases10Hz1to360';

PerceptionClassifierOscillation(cfg)


%% Source activation traces

% for realigned imagery
cfg = [];
cfg.dataDir = fullfile(output,'GroupResults','SourceActivation');
cfg.data    = 'realignedIM';
cfg.tidx    = [0.07 0.143];
cfg.area{1} = 'Pole_occipital';
cfg.area{2} = {'G_oc-temp_lat-fusifor','S_oc-temp_lat',... 
	 'S_oc-temp_med_and_Lingual'};
 
SourceTrace(cfg); 

% for perception
cfg = [];
cfg.dataDir = fullfile(output,'GroupResults','SourceActivation');
cfg.data    = 'FP';
cfg.tidx    = [0 0.4];
cfg.area{1} = 'Pole_occipital';
cfg.area{2} = {'G_oc-temp_lat-fusifor','S_oc-temp_lat',... 
	 'S_oc-temp_med_and_Lingual'};
cfg.extrema = fullfile(output,'GroupResults','EMD','extrema');
 
SourceTrace(cfg); 

% for unrealigned imagery
cfg = [];
cfg.dataDir = fullfile(output,'GroupResults','SourceActivation');
cfg.data    = 'IM';
cfg.tidx    = [0.5 1.5];
cfg.area{1} = 'Pole_occipital';
cfg.area{2} = {'G_oc-temp_lat-fusifor','S_oc-temp_lat',... 
	 'S_oc-temp_med_and_Lingual'};
 
SourceTrace(cfg);
