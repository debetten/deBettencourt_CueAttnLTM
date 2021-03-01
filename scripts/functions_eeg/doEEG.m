function eeg = doEEG(subjName,settings)
%function eeg = doEEG(sub,data_path,settings)
%
% Accomplishes much of the necessary steps of obtaining EEG data:
% 1. loads EEG data
% 2. rereferences
% 3. organizes condition codes.
%
% Inputs:
% - subjName: subject name (e.g., '0517171_rtPreStim02')
% - settings: struct with preprocessing info, created by EEG_Settings.m
%             should contain the path to the .vhdr file in settings.dir.eeg_data_path
%
% OUTPUTS:
% - eeg: structure with the eeg data and all associated details
%
% Originally written by JJF & DS, edited by MdB & SDW 06/2017


%% read in data

fprintf('Importing Data\n');


%load in the EEG file
if strcmp(subjName,'0803171_rtPreStim02')
    data_file = '0803171.1_rtPreStim02.vhdr';
elseif strcmp(subjName,'0829171_rtPreStim02')
    data_file = '08291712_rtPreStim02.vhdr';
elseif strcmp(subjName,'1013171_rtPreStim02')
    data_file = '1013172_rtPreStim02.vhdr';
else
    data_file = [subjName '.vhdr'];
end

% Read the header file
CONF = readbvconf(settings.dir.eeg_data_path, data_file); %channel info etc.

if strcmp(subjName,'0816171_rtPreStim02')
    [EEG1, ~]=pop_loadbv(settings.dir.eeg_data_path, data_file);  % requires eeglab
    [EEG2, ~]=pop_loadbv(settings.dir.eeg_data_path, '08161712_rtPreStim02.vhdr');  % requires eeglab
    EEG = pop_mergeset(EEG1,EEG2);
elseif strcmp(subjName,'0919171_rtPreStim02')
    [EEG1, ~]=pop_loadbv(settings.dir.eeg_data_path, data_file);  % requires eeglab
    [EEG2, ~]=pop_loadbv(settings.dir.eeg_data_path, '0919172_rtPreStim02.vhdr');  % requires eeglab
    EEG = pop_mergeset(EEG1,EEG2);
elseif strcmp(subjName,'1011171_rtPreStim02')
    [EEG1, ~]=pop_loadbv(settings.dir.eeg_data_path, data_file);  % requires eeglab
    [EEG2, ~]=pop_loadbv(settings.dir.eeg_data_path, '1011172_rtPreStim02.vhdr');  % requires eeglab
    EEG = pop_mergeset(EEG1,EEG2);
else
    [EEG, ~]=pop_loadbv(settings.dir.eeg_data_path, data_file);  % requires eeglab
end

% Rename things so we have them to use later!!
eeg.data = EEG.data; % changed from eeg to erp to match josh's naming schema
eeg.srate = EEG.srate;
eeg.nChans = size(EEG.data,1); % Number of Channels
eeg.nArfChans = size(EEG.data,1) - 1; % don't count the stim trak as a channel to artifact reject ! 
eeg.pnts = size(EEG.data,2); % number of overall sample data points
eeg.rateAcq = 1000/EEG.srate; % 2 ms= Rate of Data Acquisition
eeg.event = struct( 'type', { EEG.event.type }, 'latency', {EEG.event.latency});
eeg.eventTimes = round(cell2mat({eeg.event.latency})); % Event Times
eeg.headerInfo = CONF.comment;
eeg.data_file = data_file;


% loop through channels and get labels and coordintates
for c = 1:eeg.nChans  
    
    % grab channel labels
    label = strsplit(CONF.channelinfos{c},',');
    eeg.chanLabels{c} = label{1};
       
    % grab channel coordinates
    coordinates = strsplit(CONF.coordinates{c},',');
    eeg.chanCoordinates{c} = [str2double(coordinates{1}) str2double(coordinates{2}) str2double(coordinates{3})];

end 


%% Re-reference and/or drop unwanted electrodes

% TP9 is the left mastoid, TP10 is the right mastoid (and the online referece)
leftMastoid = ismember(eeg.chanLabels,'TP9')';

l = length(settings.droppedElectrodes);

% add the offline reference (TP9) to the notIncluded vector
settings.droppedElectrodes{l+1} = 'TP9';

notIncludedInd = ismember(eeg.chanLabels,settings.droppedElectrodes)';

chanLabels_reref = eeg.chanLabels(~notIncludedInd);

% eye channels and stim track are last and have their own references so don't rereference them
skipReRef = {'HEOG','VEOG','StimTrak'};

% drop TP9 from our reReferenceSet and don't include eye channels and stim
Chan2ReRef = sum(~ismember(chanLabels_reref,skipReRef));

tmpreRef = eeg.data(~notIncludedInd,:);

% r/2
mastoidValue = eeg.data(leftMastoid,:) ./ 2;

%NOTE: This for loop is set up for HEOG,VEOG, and Stimtrack (which have their own face reference) as the last 3 channels.
%If you change the order make sure you restructure the script to skip rereferencing those channels!
for chan = 1:Chan2ReRef 
    tmpreRef(chan,:) = tmpreRef(chan,:) - mastoidValue; % a - (r/2)
end

%save chan labels for moving forward
eeg.data = tmpreRef;
eeg.chanLabels = chanLabels_reref;
eeg.nChans = sum(~notIncludedInd);
eeg.nArfChans = sum(~notIncludedInd)-1; % don't count stim trak as a chan to be looked at for arfs
eeg.chanCoordinates = eeg.chanCoordinates(~notIncludedInd);   % JJF: update this....


%% Remove spaces and letters from the BrainProducts codes

eeg.eventCodes = nan(1,length(eeg.eventTimes)); 

% preallocate event codes matrix
eeg.eventCodes(1,1)= NaN; %boundary event (i.e., start of recording)

for ii = 2 : length(eeg.eventTimes)
    % get characters 2-4 (i.e. the digits). 1 is actually stored as blank-blank-1.
    s = eeg.event(1,ii).type(2:4);
    
    % covert string to double
    s = str2double(s);
    
    % save the event code structure
    eeg.eventCodes(1,ii) = s;      
end


end


