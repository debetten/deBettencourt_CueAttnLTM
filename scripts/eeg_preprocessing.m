function eeg_preprocessing(subj_name,port_codes)
%function eeg_preprocessing(subj_name,port_codes)
%
% Accomplishes the preprocessing steps for EEG analysis
%
% 1. Reads in BrainVision EEG files
% 2. Saves it as a matlab file using EEGLAB
% 3. Calculates ERPs preprocessing
%
% Originally written by JJF & DS, edited by MdB & SDW 06/2017
%
% Inputs:
% - subj_name: string file corresponding to subject name (e.g., 0517171_rtPreStim02)
% - port_codes: matrix containing port codes of interest (e.g., [30 31])
%
% all details to change are specified in the function EEG_Settings

fprintf('Preprocessing EEG data for subject:\t%s\n',subj_name)

%% filepaths

%environments
if strcmp(computer,'MACI64') %running on Megan's desktop
    proj_dir = '/Users/megan/Documents/projects/rtPreStim02';
else %running on the RA PC
    proj_dir = 'C:/Desktop/Stephanie/rtPreStim02';
end

% add folder with eeg functions
eeg_functions_dir = fullfile(proj_dir,'/scripts/functions_eeg/');
addpath(eeg_functions_dir)

% add folder with artifact detection functions
arf_functions_dir = fullfile(proj_dir,'/scripts/functions_art/');
addpath(arf_functions_dir);


%% preprocessing settings

fprintf('loading settings... \n')

settings.droppedElectrodes = {}; % electrodes to remove

% %% directory info
settings.dir.eeg_data_path = fullfile(proj_dir,'/subjects/',subj_name,'/data/eeg/');
settings.dir.eye_data_path = fullfile(proj_dir,'/subjects/',subj_name,'/data/eye/');

%segmentation settings
settings.seg.codes = port_codes; % vector of all event codes of interest
%settings.seg.codes = [10,11,21,30,31,41,50,51,60,61,71,101,102,103,104,111,112,113,114,121,122,123,131,132,133,141,142,143,151,200,201]; % vector of all event codes of interest

%timing for artifact rejection (times should be ABSOLUTE VALUES) 
settings.seg.arfPreTime = 300; % msec prior to timelock
settings.seg.arfPostTime = 1500; % msec post timelock (600 ms after array onset)

%timing stuff for building waveforms (times should be absolute values) 
settings.seg.preTime = settings.seg.arfPreTime+800; % msecs prior to timelock built in extra 800ms for time freq analyses
settings.seg.postTime = settings.seg.arfPostTime+800; % msecs post timelock

%window for baselining (time should be NEGATIVE IF PRE-=TIMELOCKING
settings.seg.baseStart = -300;  %%% if using the whole time period, just use -settings.preTime
settings.seg.baseEnd = -100; 


%% artifact rejection settings

%Noise threshold for artifact rejection
settings.arf.noiseThr = 200; % microvolts
settings.arf.noiseWin = 15; % ms (short so the peak-to-peak algorithm is selective for high-freq noise)
settings.arf.noiseStep = 50; % ms no need to check more than every 50 ms for noise

%Threshold for drift
settings.arf.driftThr = 100; %microvolts

%Step function settings for channel drop out (main cap channels sometimes have step functions when they suddenly drop out! 
% do a wide window length to avoid catching alpha!! 
settings.arf.dropoutWin = 250; %ms
settings.arf.dropoutStep = 20; % ms
settings.arf.dropoutThr = 60; % microvolts

%Step function settings for blink rejection
settings.arf.blinkWin = 150; % ms
settings.arf.blinkStep = 10; % ms
settings.arf.blinkThr = 50; % microvolts

%Step function settings for horizontal eye movements 
settings.arf.eMoveWin = 150; % ms
settings.arf.eMoveStep = 10; % ms
settings.arf.eMoveThr = 20; %microvolts

%Settings for block rejection
settings.arf.blockWin = 100; % ms
settings.arf.blockStep = 50; % ms
settings.arf.blockX = 60; % ms
settings.arf.blockY = 1; % microvolts


%% preprocess EEG data and save file (no segmentation)

tic;

%check whether EEG file exists
eegFile = [settings.dir.eeg_data_path, subj_name '_EEG.mat'];
% if exist(eegFile,'file')
%     reply = input('EEG file already exists, do you want to overwrite it? [y/n]','s');
% else
    reply = 'n';
% end

if reply == 'y' || ~exist(eegFile,'file')
    eeg = doEEG(subj_name,settings);
    fprintf('saving EEG file.. \n')
    eegName = [settings.dir.eeg_data_path,subj_name '_EEG','.mat'];
    eeg.droppedElectrodes = settings.droppedElectrodes; % only save dropped electrodes, other settings apply to segmentation
    save(eegName,'eeg');
else
    fprintf('loading EEG file.. \n')
    load(eegFile)
end

toc; % report time to process EEG data


%% do ERPs

tic;
erp = eeg;
clear eeg; % pass eeg data to erp structure

erp = doERPs(erp,settings); % doERP pipeline

fprintf('saving ERP file.. \n')
if any(ismember(port_codes,10))
    erpName = [settings.dir.eeg_data_path,num2str(subj_name),'_EEG_SEG.mat'];
else
    erpName = [settings.dir.eeg_data_path,num2str(subj_name),'_EEG_SEG_ENC.mat'];
end
erp.settings = settings; % add settings to erp struct before saving

save(erpName,'erp');
toc; % report time to doERPs


end
