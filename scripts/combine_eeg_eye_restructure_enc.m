function combine_eeg_eye_restructure_enc(subj_name,eye_track)
% Script to restructure EEG data from(Trials,chans,timepoints) to (Chans,Timepoints,Trials) for plotting with the EEG Lab plot function
% This script also concatenates the EEG data and Eye tracking data (if we have it for the subject) 
% along with the logical indices from our automated artifact rejection procedures.
%
% Inputs:
% subjectName = subject name
% eyeTrack: 1 = usable eye tracking, 0 = no eye tracking, or eye tracking unusable
%
%The script saves 2 files:
%1. A restructured .mat file (subj#_restruct_for_arf).
%2. An EEG lab file for manual marking (subj#MarkingComplete.set and .fdt


%filepaths
proj_dir = '/Users/megan/Documents/projects/rtPreStim02/';
eyeData_dir = fullfile(proj_dir,'subjects',subj_name,'data/eye');
eegData_dir = fullfile(proj_dir,'subjects',subj_name,'data/eeg');
fname_eye = [subj_name,'_EYE_SEG_ENC.mat'];
fname_eeg = [subj_name '_EEG_SEG_ENC.mat'];

% name of restructured file to be saved
fname_restruct = [eegData_dir,'/',(subj_name),'_Restruct_ENC.mat'];

%% abort if a restructured file already exists
% if exist(restructuredFilename,'file')
%     reply = input('Restructured EEG file already exists, do you want to overwrite it? [y/n]','s');
% else
%     reply = 'n';
% end
% 
% if ~strcmp(reply,'y')
%     return
% end

%% load eeg data file

fprintf('Load eeg data... ')

% load eeg data
tic; load([eegData_dir '/' fname_eeg]); toc;

% index time points of interest
startTime = -erp.settings.seg.arfPreTime;
endTime = erp.settings.seg.arfPostTime;

%get artifact window (which might be shorter than the whole ERP time
arfWindow = ismember(erp.trial.times,startTime:endTime);
nSampsEEG = sum(arfWindow);% get dimensions of data
tmpeeg = erp.trial.baselined(:,:,arfWindow);% grab eeg data of interest


%% load eye tracking data (if available)


if eye_track
    
    load([eyeData_dir '/' fname_eye]);
    
    % grab eye position in degrees of vis ang from fixation, multiply by 16
    % to covert to microvolts
    microV_H = eyeData.trial.xDeg.*16; % horizontal gaze position
    microV_V = eyeData.trial.yDeg.*16; % vertical gaze posiion
    
    %EDITED
    if eyeData.trial.nTrials ~= erp.trial.nTrials
        disp('number of trials do not match across eye tracking and eeg files')
        return   
    end

    if eyeData.sRate ~= erp.srate
        disp('EEG and eye tracking sampling rates do not matched\n')
        return
    end
    
else % otherwise, set microV variables to nans
    
    microV_H = nan(erp.trial.nTrials,nSampsEEG);
    microV_V = nan(erp.trial.nTrials,nSampsEEG);
    
end

%% ensure the sampling rates match for eyetracking and eeg data

disp('calculate nsamps\n')
if eye_track
    
    newsamps_eeg = size(tmpeeg,3);
    newsamps_eye = size(microV_H,3); %IF BINOCULAR -- USE 3, IF MONOCULAR USE 2
    
    if newsamps_eeg ~= newsamps_eye
        error('did not successfully match sampling rates of eeg and eye tracking\n')
    end
    
    nsamps = newsamps_eeg;

else
   
    nsamps = size(tmpeeg,3);
    
end


%% restructure the data so it's compatible with EEGLab

disp('restructuring\n')
if eye_track
    n_extra_elecs = 8; % 4 for eyetracker, 2 for eyetracker fixation lines, and 2 EOG fixation lines
    if size(microV_H,2)==1
        n_extra_elecs = 6;
    end
else
    n_extra_elecs =  2; % 2 EOG fixation lines
end

nelecs = erp.nChans + n_extra_elecs;

% preallocate the restructured matrix
restructured = nan(nelecs,nsamps,erp.trial.nTrials);

% put EEG data into restructured matrix
for ch = 1:erp.nChans
    for tr = 1:erp.trial.nTrials
        restructured(ch,:,tr) = squeeze(tmpeeg(tr,ch,:));
    end
end

% amount to shift by
HEOGShift = 100;
HEyeShift = 200;
VEOGShift = -200;
VEyeShift = -100;
StimTrackShift = 400;

% add eye tracking to restructured matrix  
for tr = 1:erp.trial.nTrials
    if eye_track
        % add eye tracking
        if size(microV_H,2)==2 
            restructured(nelecs-7,:,tr) = microV_H(tr,1,:);
            restructured(nelecs-6,:,tr) = microV_V(tr,1,:);
            restructured(nelecs-5,:,tr) = microV_H(tr,2,:);
            restructured(nelecs-4,:,tr) = microV_V(tr,2,:);
        else
            restructured(nelecs-5,:,tr) = microV_H(tr,1,:);
            restructured(nelecs-4,:,tr) = microV_V(tr,1,:);
        end
    end
    % add baseline lines
    restructured(nelecs-3,:,tr) = HEOGShift;
    restructured(nelecs-2,:,tr) = HEyeShift;   
    restructured(nelecs-1,:,tr) = VEOGShift;
    restructured(nelecs,:,tr) = VEyeShift;
end

% shift the position of the eye track, eog, and stim track so it's visible on the butterfly plot
% restructured(nChans-8,:,:) = restructured(nChans-8,:,:)+HEOGShift; %shift heog
% restructured(nChans-7,:,:) = restructured(nChans-7,:,:)+VEOGShift; %shift veog
% restructured(nChans-6,:,:) = restructured(nChans-6,:,:)+StimTrackShift; %shift stim track to top
% if eyeTrack
%     restructured(nChans-5,:,:) = restructured(nChans-5,:,:)+HEyeShift; %shift horizontal eye track
%     restructured(nChans-4,:,:) = restructured(nChans-4,:,:)+VEyeShift; %shift vertical eye track
% end

restructured(31,:,:) = restructured(31,:,:)+HEOGShift; %shift heog
restructured(32,:,:) = restructured(32,:,:)+VEOGShift; %shift veog
restructured(33,:,:) = restructured(33,:,:)+StimTrackShift; %shift stim track to top
if eye_track
     restructured(34,:,:) = restructured(34,:,:)+HEyeShift; %shift horizontal eye track
     restructured(35,:,:) = restructured(35,:,:)+VEyeShift; %shift vertical eye track
     
     if size(microV_H,2)==2 
         restructured(36,:,:) = restructured(36,:,:)+HEyeShift; %shift horizontal eye track
         restructured(37,:,:) = restructured(37,:,:)+VEyeShift; %shift vertical eye track
     end
end

%% save the data to the restructured file

disp('saving\n')
save([fname_restruct],'restructured','-v7.3')


%% import data into EEGLab format

% srate = sampling rate, pnts = number of samples, xmin = start of epoch (e.g., -0.3 for -300 ms)
EEG = pop_importdata('dataformat','matlab','nbchan',nelecs,'data',fname_restruct,'srate',erp.srate,'pnts',nsamps,'xmin',startTime/1000);
EEG = eeg_checkset( EEG );

%% import auto rejection information to the EEGLab file

if eye_track == 1 % eyetracking data!
    
    % BLUE: blocking, noise, or dropout in EEG
    EEG.reject(:).rejkurt = (erp.arf.blocking | erp.arf.noise | erp.arf.dropout);
    % GREEN: drift in EEG
    EEG.reject(:).rejconst = erp.arf.drift;
    % LIME GREEN: blinks, lost eyes, missing eye tracker data
    blinks = eyeData.arf.missingPupil + eyeData.arf.parserBlinks + erp.arf.blink; blinks(blinks>1) = 1; % make all 0s and 1s again
    EEG.reject(:).rejthresh = blinks;
    % RED: EyeTracker saccades or deviations from fixation
    EEG.reject(:).rejjp = eyeData.arf.saccadeX | eyeData.arf.saccadeY ;
    % PINK/PURPLE: EOG saccade reject
    % EEG.reject(:).rejfreq = erp.arf.eMove;
    
else % no eyetracking data
    
    % BLUE: blocking, noise, or dropout
    EEG.reject(:).rejkurt = (erp.arf.blocking | erp.arf.noise | erp.arf.dropout);
    % GREEN: drift
    EEG.reject(:).rejconst = erp.arf.drift;
    % LIME GREEN: blinks
    EEG.reject(:).rejthresh = erp.arf.blink;
    % RED: not assigned
    %EEG.reject(:).rejjp =
    % PINK/PURPLE: EOG saccade reject
    EEG.reject(:).rejfreq = erp.arf.eMove;
    
end

%% Save .set file with the artifact markers
EEG = eeg_checkset(EEG);
ename = [num2str(subj_name),'_Art_ENC'];
EEG = pop_saveset( EEG,'filename',ename,'filepath',eegData_dir);

