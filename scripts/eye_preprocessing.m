function eyeData = eye_preprocessing(subj_name)
%function eyeData = eye_preprocessing(subj_name)
%
% CUE
%
% OUTPUTS
% - eyeData: structure containing eye data from each trial

fprintf('Preprocessing Eye-tracking data for subject:\t%s\n',subj_name)


%% filepaths

%environments
if strcmp(computer,'MACI64') %running on Megan's desktop
    proj_dir = '/Users/megan/Documents/projects/rtPreStim02';
    packages_dir = '/Users/megan/Documents/packages';
else %running on the RA PC
    proj_dir = 'C:/Desktop/Stephanie/rtPreStim02';
end

% add folder with eeg functions
eye_functions_dir = fullfile(proj_dir,'/scripts/functions_eye/');
addpath(eye_functions_dir)

% add folder with artifact detection functions
arf_functions_dir = fullfile(proj_dir,'/scripts/functions_art/');
addpath(arf_functions_dir);

% add folder with edf mex  functions
edfmex_functions_dir = fullfile(packages_dir,'edfmex');
addpath(edfmex_functions_dir);

%% load in the preprocessing settings.

fprintf('loading eye tracking settings... \n')

%directory
settings.dir.eye_data_path = fullfile(proj_dir,'subjects',subj_name,'data/eye/');

% eye tracker settings
settings.sRate = 1000; % sampling rate (Hz)
settings.rateAcq = 1000/settings.sRate; % 500 Hz = 2 ms rateAcq
%if subjName == 
settings.recordingMode = 'ChinRest_Binocular'; % 'RemoteMode_Monocular', 'RemoteMode_Binocular', 'ChinRest_Monocular', 'ChinRest_Binocular'

% key distances
% settings.cam2screenDist % distance from monitor to eye tracker, if recorded
settings.viewDist = 75; % viewing distance (cm)

% monitor details
settings.monitor.xPixels = 1920;
settings.monitor.yPixels = 1080;
settings.monitor.pxSize = 0.0276; % 53cm width / 1920 pixels wide

%segmentation settings
settings.seg.timeLockMessage = 'CUE'; % message for time locking
settings.seg.preTime = 300;  % pre-stimulus end of segment, absolute value (ms)
settings.seg.postTime = 1500; % post-stimulus end of segment, absolute value (ms)

%artifact rejection settings - for doing our artifact rejection
settings.arf.winSize = 80; % ms --- size of sliding window that is looking for saccades
settings.arf.stepSize = 10;  % ms ---- how much the window moves by 
settings.arf.maxDeg = .5; % degrees of visual angle! if it's bigger than this reject it 


%% preprocess eye track data and save file (no segmentation)

%if ~exist([settings.dir.eye_data_path,subjName,'_EYE.mat'],'file')
    fname = ['PS2_' subj_name(3:4) '_' subj_name(7)];
    f_dir_name = [settings.dir.eye_data_path,'/',fname,'.edf'];
    if ~exist(f_dir_name,'file')
        error(['Cannot find the file at ',f_dir_name])
    end
    eye = edfmex([settings.dir.eye_data_path,'/',fname,'.edf']); % read in edf file
    
    
    % check the specified sampling rate is correct
    if settings.sRate ~= eye.RECORDINGS(1).('sample_rate')
        error('specified sampling rate does not match data file.')
    end
    
    % check the specified recording mode is correct
    recordingMode = getRecordingMode(eye);
    if ~strcmp(subj_name, '0905171_rtPreStim02')
        if ~strcmp(recordingMode,settings.recordingMode) % if these strings don't match...
            fprintf('specified recording mode does not match data file: %s vs. %s',recordingMode,settings.recordingMode)
            settings.recordingMode = recordingMode;
        end
    end
    % save the data file
    save([settings.dir.eye_data_path,subj_name,'_EYE.mat'],'eye')
    
% else
%     load([settings.dir.eye_data_path,subjName,'_EYE.mat'])
% end

%% preprocess eyetracking data and save file

% sampling rate
eyeData.sRate = getSamplingRate(eye);
eyeData.rateAcq = 1000./eyeData.sRate; % calculate rate of data acquisition (ms)

% message and codestrings
eyeData.messages = {eye.FEVENT(:).message}; % grab messages sent from experiment display (Psychtoolbox/Psychopy)
eyeData.codestrings = {eye.FEVENT(:).codestring}; % grab codestrings (includes STARTSACC,ENDSACC,STARTFIX,ENDFIX,STARTBLINK,ENDBLINK among other things)
eyeData.eventTimes = [eye.FEVENT(:).sttime]; % when events occured

% which eye was recorded on each trial
RecordedEyeVec = [eye.RECORDINGS(:).('eye')]; % the eye that was tracked (left or right)
RecordedEyeIdx = 1:2:length(RecordedEyeVec); % RECORDINGS taken at start and end of each trial, only grab the starts (i.e. odd entries)
eyeData.RecordedEye=RecordedEyeVec(RecordedEyeIdx);

% eye tracking data
eyeData.sampleTimes = [eye.FSAMPLE(:).time]; % the times at which data was sampled
eyeData.gx = [eye.FSAMPLE(:).gx]; % gaze referenced x coords
eyeData.gy = [eye.FSAMPLE(:).gy]; % gaze referenced y coords
eyeData.hx = [eye.FSAMPLE(:).hx]; % head referenced x coords
eyeData.hy = [eye.FSAMPLE(:).hy]; % head referenced y coords
eyeData.pa = [eye.FSAMPLE(:).pa]; % head referenced pupil size / area

% get distance to eye tracker if using remote mode, otherwise store vector of nans
if strcmp(settings.recordingMode,'RemoteMode_Monocular') || strcmp(settings.recordingMode,'RemoteMode_Binocular')
    eyeData.dist = (double((eye.FSAMPLE(:).hdata(3,:))))./100; % 3rd row is distance, divide by 100 to scale to cm
else
    eyeData.dist = nan(size(eyeData.sampleTimes)); % same size as eyeData.sampleTimes
end

%% Segment data

timeLockInd = strcmp(settings.seg.timeLockMessage,eyeData.messages); % index the time-locking message (e.g., 'StimOnset')
eyeData.trial.timeLockTimes = eyeData.eventTimes(timeLockInd); % times for time-locking messsage
eyeData.trial.nTrials = sum(timeLockInd); % adds up logical index to get number of trials

% throw an error if no trials were found
if eyeData.trial.nTrials == 0
    error('Did not find any trials. Did you specify the right event marker in the settings file?')
end

% save times vector for each trial
eyeData.trial.times = -settings.seg.preTime:eyeData.rateAcq:settings.seg.postTime; % time points in segment
eyeData.trial.nSamps = length(eyeData.trial.times); % expected number of samples per segment

% specify start and end times of each segment
eyeData.trial.startTimes = double(eyeData.trial.timeLockTimes) - settings.seg.preTime; % start time, ms
eyeData.trial.endTimes = double(eyeData.trial.timeLockTimes) + settings.seg.postTime;  % end time, ms

% preallocate matrices for segmented data (all the same size)
eyeData.trial.gx = nan(eyeData.trial.nTrials,2,eyeData.trial.nSamps);
eyeData.trial.gy = eyeData.trial.gx;
eyeData.trial.hx = eyeData.trial.gx;
eyeData.trial.hy = eyeData.trial.gx;
eyeData.trial.pa = eyeData.trial.gx;
eyeData.trial.dist = eyeData.trial.gx;
eyeData.trial.exist = nan(eyeData.trial.nTrials,eyeData.trial.nSamps);

% loop through trials and segment data
for t = 1:eyeData.trial.nTrials
    
    % grab the start and end of trial t
    tStart = eyeData.trial.startTimes(t); tEnd = eyeData.trial.endTimes(t);
    
    % specify window of interest
    tWindow = tStart:double(eyeData.rateAcq):tEnd; 
    
    % index times of interest with logical
    tWindowInd = ismember(eyeData.sampleTimes,tWindow);
    
    % throw an error if sampling rate is less than 500 Hz
    if eyeData.rateAcq > 1
        error('Sampling rate lower than 1000 Hz. Have not prepared the fix above for sampling freqs lower than 1000 Hz')
    end
    
    % create index of the time points that actually exist in the data (i.e., that were recorded).
    existInd = ismember(tWindow,(eyeData.sampleTimes));
    
    % determine which eye was recorded for trial t
    if eyeData.RecordedEye(t)==3
        recordedEye = 1:2;
    else
        recordedEye = eyeData.RecordedEye(t); %MDB CHECK THIS
    end
    
    % grab the relevant segment of data (from the recorded eye)
    if numel(recordedEye)<2
        eyeData.trial.gx(t,existInd) = eyeData.gx(recordedEye,tWindowInd);
        eyeData.trial.gy(t,existInd) = eyeData.gy(recordedEye,tWindowInd);
        eyeData.trial.hx(t,existInd) = eyeData.hx(recordedEye,tWindowInd);
        eyeData.trial.hy(t,existInd) = eyeData.hy(recordedEye,tWindowInd);
        eyeData.trial.pa(t,existInd) = eyeData.pa(recordedEye,tWindowInd);
        eyeData.trial.dist(t,existInd) = eyeData.dist(tWindowInd);
    else %CAN I USE BOTH EYES INDEPENDENTLY??? FOR NOW, MEAN
        eyeData.trial.gx(t,:,existInd) = (eyeData.gx(:,tWindowInd));
        eyeData.trial.gy(t,:,existInd) = (eyeData.gy(:,tWindowInd));
        eyeData.trial.hx(t,:,existInd) = (eyeData.hx(:,tWindowInd));
        eyeData.trial.hy(t,:,existInd) = (eyeData.hy(:,tWindowInd));
        eyeData.trial.pa(t,:,existInd) = (eyeData.pa(:,tWindowInd));
        eyeData.trial.dist(t,existInd) = eyeData.dist(tWindowInd);        
    end
    % save exist to the trial structure to make it easy to check where data is missing
    eyeData.trial.exist(t,:) = existInd;
    
end

% plot the missing data to alert experimenter to problems
figure; imagesc(eyeData.trial.exist);
title('Missing samples (1 = data present, 0 = data missing)')
xlabel('Samples')
ylabel('Trials')
colorbar

%% calculate eye position in degrees of vis angle from fixation

% if data collected in remote mode
if strcmp(settings.recordingMode,'RemoteMode_Monocular') || strcmp(settings.recordingMode,'RemoteMode_Binocular')
    
    % check that settings.cam2screenDist exists, then use it!
    if isfield(settings,'cam2screenDist')
        cam2screenDist = settings.cam2screenDist;
    else % otherwise estimate it from viewing dist and head2cam dist
        cam2screenDist = approx_cam2screenDist(settings.viewDist,eyeData.dist);
    end
    
    % calculate degrees of visual angle
    [eyeData.trial.xDeg,eyeData.trial.yDeg] = pix2deg_remoteMode(eyeData.trial.gx,eyeData.trial.gy,eyeData.trial.dist,...
                                                cam2screenDist,settings.monitor.xPixels,settings.monitor.yPixels,settings.monitor.pxSize);
end

% if data collected with the chin rest...
if strcmp(settings.recordingMode,'ChinRest_Monocular') || strcmp(settings.recordingMode,'ChinRest_Binocular') || strcmp(settings.recordingMode,'Tower_Binocular') 
    
    % calculate degrees of visual angle
    %[eyeData.trial.xDeg,eyeData.trial.yDeg] = pix2deg_chinRest(eyeData.trial.gx,eyeData.trial.gy,...
    %                                            settings.monitor.xPixels,settings.monitor.yPixels,settings.monitor.pxSize,settings.viewDist);
    %[degHfromFix degVfromFix] = pix2deg_chinRest(pixH,pixV,pixelsH,pixelsV,pxSize,viewDist)
    % calculate pixels from the middle of the screen
    pixHfromFix = eyeData.trial.gx-(settings.monitor.xPixels/2);
    pixVfromFix = eyeData.trial.gy-(settings.monitor.yPixels/2);
    
    % convert these values to cm to calculate degrees of visual angle
    cmHfromFix = pixHfromFix.*settings.monitor.pxSize;
    cmVfromFix = pixVfromFix.*settings.monitor.pxSize;
    
    % calculate degrees of visual angle from fixation
    eyeData.trial.xDeg = atand(cmHfromFix./settings.viewDist);
    eyeData.trial.yDeg = atand(cmVfromFix./settings.viewDist);
end


%% Artifact rejection: mark bad data

% mark bad data based on eyelink parser
% grab relevant codestrings
blinkStartInd = strcmp('STARTBLINK ',eyeData.codestrings); %space here is necessary
blinkEndInd = strcmp('ENDBLINK',eyeData.codestrings);
saccStartInd = strcmp('STARTSACC',eyeData.codestrings);
saccEndInd = strcmp('ENDSACC',eyeData.codestrings);

% grab the times for events of interest, save to arf structure
parserTimes.blinkStart = eyeData.eventTimes(blinkStartInd); % times for blink start, and so on...
parserTimes.blinkEnd = eyeData.eventTimes(blinkEndInd);
parserTimes.saccStart = eyeData.eventTimes(saccStartInd);
parserTimes.saccEnd = eyeData.eventTimes(saccEndInd);

nTrials = eyeData.trial.nTrials;

% preallocate matrices for detected saccades and blinks
rejBlMat = zeros(1,nTrials);
rejSaccMat = zeros(1,nTrials);

% loop through trials and check whether they contained artifacts
for t = 1:nTrials
    
    % get trial start and trial end
    tStart = eyeData.trial.startTimes(t);
    tEnd = eyeData.trial.endTimes(t);
    tWindow = tStart:tEnd;
    
    % was there a blink during the trial?
    if sum(ismember(tWindow,parserTimes.blinkStart))>0 || sum(ismember(tWindow,parserTimes.blinkEnd))>0
        rejBlMat(t) = 1;
    end
    
    % was there a saccade during the trial?
    if sum(ismember(tWindow,parserTimes.saccStart))>0 || sum(ismember(tWindow,parserTimes.saccEnd))>0
        rejSaccMat(t) = 1;
    end
    
end

eyeData.arf.parserBlinks = logical(rejBlMat); % logical of if there were blinks
eyeData.arf.parserSaccs = logical(rejSaccMat); % logical of if there were saccades


% run our own check for artifacts

% preallocate vectors - FOR BINOCULAR DATA ONLY!!!! NEED TO MAKE THIS FLEXIBLE
missingPupil = nan(2,eyeData.trial.nTrials);
saccadeX = missingPupil;
saccadeY = missingPupil;

% loop through trials
for t = 1:eyeData.trial.nTrials
    
    for iEye = 1:2
        %grab gaze data for current trial
        xGaze = squeeze(eyeData.trial.gx(t,iEye,:));
        yGaze = squeeze(eyeData.trial.gy(t,iEye,:));
        xDeg = squeeze(eyeData.trial.xDeg(t,iEye,:));
        yDeg = squeeze(eyeData.trial.yDeg(t,iEye,:));
        
        %mark trials where the eye tracker lost the pupil (e.g. blinks)
        mpx = zeros(size(xGaze));
        mpy = zeros(size(yGaze));
        
        mpx(xGaze > 10000) = 1;
        mpy(yGaze > 10000) = 1;
        
        if (sum(mpx)>0) || (sum(mpy)>0); % mark if missing pupil in x or y data
            missingPupil(iEye,t) = 1;
        else
            missingPupil(iEye,t) = 0;
        end
        
        %run step function as a second check for saccades
        % check xgaze
        stepX = art_step(xDeg,eyeData.rateAcq,settings.arf.stepSize,settings.arf.winSize,settings.arf.maxDeg);
        if sum(stepX) > 0
            saccadeX(iEye,t) = 1;
        else
            saccadeX(iEye,t) = 0;
        end
        
        % check ygaze
        stepY = art_step(yDeg,eyeData.rateAcq,settings.arf.stepSize,settings.arf.winSize,settings.arf.maxDeg);
        if sum(stepY) > 0
            saccadeX(iEye,t) = 1;
        else
            saccadeX(iEye,t) = 0;
        end
    end
    
end

% covert to logicals
eyeData.arf.saccadeX = logical(sum(saccadeX)==2);
eyeData.arf.saccadeY = logical(sum(saccadeY)==2);
eyeData.arf.missingPupil = logical(sum(missingPupil)==2);

% print artifact summary
fprintf(sprintf('Parser Blink Rate: %.2f \n',sum(eyeData.arf.parserBlinks)./length(eyeData.arf.parserBlinks)))
fprintf(sprintf('Parser Saccade Rate: %.2f \n',sum(eyeData.arf.parserSaccs)./length(eyeData.arf.parserSaccs)))
fprintf(sprintf('Calculated nonsensical position vals: %.2f \n',sum(eyeData.arf.missingPupil)./length(eyeData.arf.missingPupil)))
fprintf(sprintf('Calculated horizontal saccades: %.2f \n',sum(eyeData.arf.saccadeX)./length(eyeData.arf.saccadeX)))
fprintf(sprintf('Calculated vertical saccades: %.2f \n',sum(eyeData.arf.saccadeY)./length(eyeData.arf.saccadeY)))


%% Remove continuos data from eyeData structure

eyeData = rmfield(eyeData,'gx');
eyeData = rmfield(eyeData,'gy');
eyeData = rmfield(eyeData,'hx');
eyeData = rmfield(eyeData,'hy');
eyeData = rmfield(eyeData,'pa');
eyeData = rmfield(eyeData,'dist');


%% save

eyeData.settings = settings; % save settings to eyeData structure
save([settings.dir.eye_data_path,'/',num2str(subj_name),'_EYE_SEG.mat'],'eyeData')


end
