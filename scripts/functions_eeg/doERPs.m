function erp = doERPs(erp,settings)
%function erp = doERPs(erp,settings)
%
% create ERPs (segmentation, artifact rejection etc.)
%
% Written by JJF & DS, edited by MdB & SDW 06/2017
%
% INPUTS
% - erp:
% - settings:

fprintf('Processing ERPs \n')

%% Segment data (including buffer time for time-freq analyses)

fprintf('segmenting data... \n')

% Determine onset of each trial
erp.trial.times = -settings.seg.preTime:erp.rateAcq:settings.seg.postTime;% the sample times in segments
preTimeSamp = settings.seg.preTime./erp.rateAcq;         % # of samples we need to go back from event code
postTimeSamp = settings.seg.postTime./erp.rateAcq;       % # of samples we need to go forward from event
tLength = preTimeSamp + postTimeSamp + 1;                % trial length in samples
codesInd = ismember(erp.eventCodes,settings.seg.codes);  % index where these codes of interest occured           
erp.trial.nTrials = sum(codesInd);                       % calculate number of trials

% preallocate matrices
erp.trial.data = nan(erp.trial.nTrials,erp.nChans,tLength);
erp.trial.codes = nan(erp.trial.nTrials,1);

tCnt = 1;
% loop through all event codes
for ii = 1:length(erp.eventCodes)
    if codesInd(ii) % if a code of interest, grab the segment!
        
    % Determine start and stop of trial
    tStart = erp.eventTimes(ii)-preTimeSamp; % jjf: rename this??
    tEnd = erp.eventTimes(ii)+postTimeSamp;   %jjf: rename this??
    tWindow = tStart:tEnd;
    
    % get time-series data for segment
    rawTS = erp.data(:,tWindow);
    
    % save time-series data to data matrix
    erp.trial.data(tCnt,1:erp.nChans,:) = rawTS;
    
    % grab trial lable....
    erp.trial.codes(tCnt) = erp.eventCodes(ii);
    
    tCnt = tCnt + 1; % advance trial indexing counter

    end
end


%% artifact Rejection: mark bad data

fprintf('checking scalp channels for artifacts... \n');
tic

% Preallocate matrices 
erp.arf.blockingFull = nan(erp.nChans,erp.trial.nTrials);
erp.arf.noiseFull = nan(erp.nChans,erp.trial.nTrials);
erp.arf.driftFull = nan(erp.nChans,erp.trial.nTrials);
erp.arf.dropoutFull = nan(erp.nChans,erp.trial.nTrials);

erp.arf.blink = nan(1,erp.trial.nTrials);
erp.arf.eMove = nan(1,erp.trial.nTrials);
erp.arf.artifactInd = nan(1,erp.trial.nTrials);
erp.arf.grand = zeros(erp.nChans,erp.trial.nTrials);  % matrix for labeling all artifacts


%loop through channels and check for artifacts
for i = 1:erp.nChans
    
    % grab data the portion of the segment to apply artifact rejection
    arf_tois = ismember(erp.trial.times,-settings.seg.arfPreTime:settings.seg.arfPostTime);
    chanDat = squeeze(erp.trial.data(:,i,arf_tois));
    
    for t = 1:erp.trial.nTrials

        % get raw time-series for the trial and electrode
        rawTS = chanDat(t,:);
        
        % check for blocking in all channels except StimTrak
        checkChannel = ~ismember(erp.chanLabels,'StimTrak');  % specify names of channels to skip
        if checkChannel(i)        
            block = art_block(rawTS,erp.rateAcq,settings.arf.blockStep,settings.arf.blockWin,settings.arf.blockX,settings.arf.blockY);
            if sum(block) > 0
                erp.arf.blockingFull(i,t) = 1;
            else
                erp.arf.blockingFull(i,t) = 0;
            end
        end
        
        % check for noise in all scalp channels
        checkChannel = ~ismember(erp.chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip        
        if checkChannel(i)
            noise = art_ppa(rawTS,erp.rateAcq,settings.arf.noiseStep,settings.arf.noiseWin,settings.arf.noiseThr);
            if sum(noise) > 0
                erp.arf.noiseFull(i,t) = 1;
            else
                erp.arf.noiseFull(i,t) = 0;
            end
        end
        
        % check for extreme drift in all scalp channles 
        checkChannel = ~ismember(erp.chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip        
        if checkChannel(i)
            drift = art_drift(rawTS,erp.rateAcq,settings.arf.driftThr);
            if sum(drift) > 0
                erp.arf.driftFull(i,t) = 1;
            else
                erp.arf.driftFull(i,t) = 0;
            end
        end

        % check for extreme channel drop out (step function) in all scalp channles
        checkChannel = ~ismember(erp.chanLabels,{'HEOG','VEOG','StimTrak'}); % specify names of channels to skip
        if checkChannel(i)
            dropout = art_step(rawTS,erp.rateAcq,settings.arf.dropoutStep,settings.arf.dropoutWin,settings.arf.dropoutThr);
            if sum(dropout) > 0
                erp.arf.dropoutFull(i,t) = 1;
            else
                erp.arf.dropoutFull(i,t) = 0;
            end
        end
        
    end
end
toc;

%Check for blinks
fprintf('checking for blinks... \n');
tic;
erp.arf.blink = nan(1,erp.trial.nTrials);
veogDat = squeeze(erp.trial.data(:,ismember(erp.chanLabels,'VEOG'),arf_tois));
for t = 1:erp.trial.nTrials
    rawTS = veogDat(t,:);
    % Check for blinks using step function
    blink = art_step(rawTS,erp.rateAcq,settings.arf.blinkStep,settings.arf.blinkWin,settings.arf.blinkThr);
    if sum(blink) > 0
        erp.arf.blink(t) = 1;
    else
        erp.arf.blink(t) = 0;
    end
end
toc;

%check for eye movements
fprintf('checking for eye movements... \n');
tic;
erp.arf.eMove = nan(1,erp.trial.nTrials);
heogDat = squeeze(erp.trial.data(:,ismember(erp.chanLabels,'HEOG'),arf_tois));
%%%%% check for horizontal eye movements 
for t = 1:erp.trial.nTrials
        rawTS = heogDat(t,:);
        % check for eye movements using step function
        eMoveH = art_step(rawTS,erp.rateAcq,settings.arf.eMoveStep,settings.arf.eMoveWin,settings.arf.eMoveThr);     
        if sum(eMoveH) > 0
            erp.arf.eMove(t) = 1;
        else
            erp.arf.eMove(t) = 0;
        end
end
toc;

%create a vector each trials as having an artifact or not
erp.arf.noise = summarizeArtifacts(erp.arf.noiseFull);
erp.arf.blocking = summarizeArtifacts(erp.arf.blockingFull);
erp.arf.drift = summarizeArtifacts(erp.arf.driftFull);
erp.arf.dropout = summarizeArtifacts(erp.arf.dropoutFull);

%artSum = squeeze(nansum(artifactMatrix,1));
%nn = 0;
%artifacts = ~ismember(e,nn);

%loop through trials and create an index of all artifacts
for t = 1:erp.trial.nTrials
    erp.arf.artifactInd(t) = erp.arf.blocking(t) | erp.arf.noise(t) | erp.arf.blink(t) | erp.arf.eMove(t) | erp.arf.drift(t) | erp.arf.dropout(t);
end


%% save rejection statistics
    
erp.arf.totalArtProp = (sum(sum(erp.arf.artifactInd))/(erp.trial.nTrials)).*100;
erp.arf.blockingProp = (sum(sum(erp.arf.blocking))/(erp.trial.nTrials)).*100;
erp.arf.noiseProp = (sum(sum(erp.arf.noise))/(erp.trial.nTrials)).*100;
erp.arf.blinkProp =  (sum(sum(erp.arf.blink))/(erp.trial.nTrials)).*100;
erp.arf.eMoveProp = (sum(sum(erp.arf.eMove))/(erp.trial.nTrials)).*100;
erp.arf.driftProp =  (sum(sum(erp.arf.drift))/(erp.trial.nTrials)).*100;
erp.arf.dropoutProp =  (sum(sum(erp.arf.dropout))/(erp.trial.nTrials)).*100;

% print proportion of trials lost due to each kind of artifact.
fprintf('%d \tPercent trials rejected (total) \n', (round(erp.arf.totalArtProp)));
fprintf('%d \tPercent blocking \n', round(erp.arf.blockingProp));
fprintf('%d \tPercent noise \n', round(erp.arf.noiseProp));
fprintf('%d \tPercent eye movements \n', round(erp.arf.eMoveProp));
fprintf('%d \tPercent blinks \n', round(erp.arf.blinkProp));
fprintf('%d \tPercent drift \n', round(erp.arf.driftProp));
fprintf('%d \tPercent chan dropout \n', round(erp.arf.dropoutProp));


%% Baseline correction

fprintf('baselining data... \n')

% get dimensions from data
nTrials = size(erp.trial.data,1); 
nChans = size(erp.trial.data,2); 
%nSamps = size(erp.trial.data,3);

% create baseline index
bInd = ismember(erp.trial.times,settings.seg.baseStart:settings.seg.baseEnd); 

% preallocate array for baselined data
erp.trial.baselined = nan(size(erp.trial.data));
erp.trial.baselineCorrection = nan(nTrials,nChans);

for t = 1:nTrials
    for chan = 1:nChans
        dat = squeeze(erp.trial.data(t,chan,:));    % grab the raw time-series
        base = mean(dat(bInd));           % calcuate mean during baseline period
        erp.trial.baselined(t,chan,:) = dat-base;   % do baseline subtraction
        erp.trial.baselineCorrection(t,chan) = base;    % save how much the data was shifted by (for easy undoing of baseline correction)
    end
end


%% remove unwanted data for ERP file

erp.trial = rmfield(erp.trial,'data'); % just keep the baseline corrected data
erp = rmfield(erp,'data'); % ditch unsegmented data - it's saved in the other file
erp = rmfield(erp,'eventTimes'); % only relevant to unsegmented data
erp = rmfield(erp,'event'); % only relevant to unsegmented data
erp = rmfield(erp,'eventCodes'); % only relevant to unsegmented data


