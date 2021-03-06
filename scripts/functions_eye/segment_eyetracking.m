function trial = segment_eyetracking(eyeData,settings)
% segment eye tracker data
%
% Inputs:
% eyeData: struct with relevant eye tracking variables
% settings: settings struct, which specifies segmentation parameters
%
% Outputs:
% trial: structure with segmented data and segmentation parameters
%
% Note: Past versions of this function did not properly handle missing data
% points at the start of the trial. If data points were missing, the first
% available sample was assumed to be the first wanted sample point. So
% instead of there being NaNs at the start of the trial, they ended up at
% the end of the trial. The problem has been corrected in this function.
%
% IMPORTANT FIX: sometimes the recorded time for an event marker
% doesn't line up with the times at which data was sampled. For
% example, when sampling at 500 Hz, we might have data sampled at 1000
% ms, 1002 ms, 1004 ms. However, it's possible that a time-locking
%  message was sent at 1003 ms. In this case, we would try to
% get the time points at 1001 ms, 1003 ms, 1005 ms etc. but these don't
% exist! To deal with this problem, we shift the tWindow index by 1 ms
% if no data was found. This fix isn't necessary for a sampling rate of 1000 Hz. 
% See "fix for when marker time is out of sync with sample points" below

timeLockInd = strcmp(settings.seg.timeLockMessage,eyeData.messages); % index the time-locking message (e.g., 'StimOnset')
trial.timeLockTimes = eyeData.eventTimes(timeLockInd); % times for time-locking messsage
trial.nTrials = sum(timeLockInd); % adds up logical index to get number of trials

% throw an error if no trials were found
if trial.nTrials == 0
    error('Did not find any trials. Did you specify the right event marker in the settings file?')
end

% save times vector for each trial
trial.times = -settings.seg.preTime:eyeData.rateAcq:settings.seg.postTime; % time points in segment
trial.nSamps = length(trial.times); % expected number of samples per segment

% specify start and end times of each segment
trial.startTimes = double(trial.timeLockTimes) - settings.seg.preTime; % start time, ms
trial.endTimes = double(trial.timeLockTimes) + settings.seg.postTime;  % end time, ms

% preallocate matrices for segmented data (all the same size)
trial.gx = nan(trial.nTrials,trial.nSamps);
trial.gy = trial.gx;
trial.hx = trial.gx;
trial.hy = trial.gx;
trial.pa = trial.gx;
trial.dist = trial.gx;
trial.exist = trial.gx;

% loop through trials and segment data
for t = 1:trial.nTrials
    
    % grab the start and end of trial t
    tStart = trial.startTimes(t); tEnd = trial.endTimes(t);
    
    % specify window of interest
    tWindow = tStart:double(eyeData.rateAcq):tEnd; 
    
    % index times of interest with logical
    tWindowInd = ismember(single(eyeData.sampleTimes),tWindow);
    
    % fix for when marker time is out of sync with sample points
    if eyeData.rateAcq == 2
        if sum(tWindowInd) == 0
            tWindow = tWindow-1;
            tWindowInd = ismember(single(eyeData.sampleTimes),tWindow);
        end
    end
    
    % throw an error if sampling rate is less than 500 Hz
    if eyeData.rateAcq > 2
        error('Sampling rate lower than 500 Hz. Have not prepared the fix above for sampling freqs lower than 500 Hz')
    end
    
    % create index of the time points that actually exist in the data (i.e., that were recorded).
    existInd = ismember(tWindow,single(eyeData.sampleTimes));
    
    % determine which eye was recorded for trial t
    if eyeData.RecordedEye(t)==3
        recordedEye = 1:2;
    else
        recordedEye = eyeData.RecordedEye(t); %MDB CHECK THIS
    end
    
    % grab the relevant segment of data (from the recorded eye)
    if numel(recordedEye)>2
        trial.gx(t,existInd) = eyeData.gx(recordedEye,tWindowInd);
        trial.gy(t,existInd) = eyeData.gy(recordedEye,tWindowInd);
        trial.hx(t,existInd) = eyeData.hx(recordedEye,tWindowInd);
        trial.hy(t,existInd) = eyeData.hy(recordedEye,tWindowInd);
        trial.pa(t,existInd) = eyeData.pa(recordedEye,tWindowInd);
        trial.dist(t,existInd) = eyeData.dist(tWindowInd);
    else %CAN I USE BOTH EYES INDEPENDENTLY??? FOR NOW, MEAN
        trial.gx(t,existInd) = mean(eyeData.gx(recordedEye,tWindowInd));
        trial.gy(t,existInd) = meaN(eyeData.gy(recordedEye,tWindowInd));
        trial.hx(t,existInd) = mean(eyeData.hx(recordedEye,tWindowInd));
        trial.hy(t,existInd) = mean(eyeData.hy(recordedEye,tWindowInd));
        trial.pa(t,existInd) = mean(eyeData.pa(recordedEye,tWindowInd));
        trial.dist(t,existInd) = eyeData.dist(tWindowInd);        
    end
    % save exist to the trial structure to make it easy to check where data is missing
    trial.exist(t,:) = existInd;
    
end

% plot the missing data to alert experimenter to problems
figure; imagesc(trial.exist);
title('Missing samples (1 = data present, 0 = data missing)')
xlabel('Samples')
ylabel('Trials')
colorbar