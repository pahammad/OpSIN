% code for converting old benchmark results ('data' matrix) into new format
% ('result' structure)

clear all; close all; clc;

generalDataRootDirectory= '/Users/ahammadp/Documents/parvez-laptop/work/data/SpikeSorting/OpSIN_BenchmarkDatasets/David_Noisy/';


warp_type = 'PCR';
% min. peak height for peak-finding step
threshold_low = 0.01; 
%@20KHz ; default: 2 times the length of expected spike
spikeTemplateWidth = 180; 
% Decides how many light-stimulus windows to skip for sampling the spike templates - useful for speeding up unsupervised clustering part
LightStimSkip = 1; 
% 1 to do bandpass filtering, 0 otherwise
bandpass_prefilter_flag = 1; 

% sampling frequency
F_s = 20*10^3; 
% lower cut-off frequency for bandpass pre-filter
low_cutoff_freq = 5; 
% high cut-off frequency for bandpass pre-filter
high_cutoff_freq = 5*10^3;
lightChannelNum = 3;
% value that defines the light windows where templates should be collected
lightChannelONthresh = 0;
dataChannelNum = 1;
% max number of clusters after clustering the light-window spike candidates
max_display_clusters = 8; 
% min number of clusters after clustering the light-window spike candidates
min_display_clusters = 3; 
procrustes_scaling_flag = 0;

locs_max = 1000; % maximum number of spikes to be manually annotated %ML

%************************************************************

%-------------------------------------
%% set up the processing module flags
%-------------------------------------

fprintf('\n\n Pick the original data file:\n');
[fName directory fOpenFlag] = uigetfile(generalDataRootDirectory);

if fOpenFlag==0
    fprintf('\n File selection failed !!\n');
    return;
end
fprintf('--------------------------------------------\n');
fprintf('\nEXPECTED INPUT DATA FORMAT:\n');
fprintf('Channel-1: electrode data\n');
fprintf('Channel-2: odor stimulus presentation information\n');
fprintf('Channel-3: light stimulus presentation information\n');
fprintf('\nIf your data does not fit this format, please reformat your data!!\n');
fprintf('--------------------------------------------\n');

fprintf('Current data directory is: %s\n', directory);
fprintf('Current working file name is: %s\n', fName);
fprintf('Current light stimulus channel: %d\n', lightChannelNum);
fprintf('Current data channel: %d\n', dataChannelNum);
fprintf('--------------------------------------------\n');


% 1 to apply bandpass pre-filtering, 0 for no pre-filtering
%bandpass_prefilter_flag = input('\n Enter 1 to prefilter the data and 0 for no pre-filtering:');

fprintf('Sampling frequency is: %d Hz\n', F_s);
if bandpass_prefilter_flag ==1
    fprintf('Lower cut-off frequency is: %d Hz\n', low_cutoff_freq);
    fprintf('Upper cut-off frequency is: %d Hz\n', high_cutoff_freq);
end

Process_Flags_OK = input('\nEnter 1 to continue:');
if Process_Flags_OK ~= 1
    return;
end

%-------------------------------------
%% set up the data input
%-------------------------------------

data = importdata([directory fName]);

if size(data,1)<size(data,2)
    data = data';
end
% formatting the input data to follow the structure specified above:
if size(data,2)>3
    foo = data(:,1:3);
    clear data;
    data = foo;
    clear foo;
end

origData = data;
%-------------------------------------
%% Pre-bandpass-filter the data to remove electrical artifacts
%-------------------------------------

if bandpass_prefilter_flag == 1
    foo = bandpassmu(data(:,1), F_s, low_cutoff_freq, high_cutoff_freq);
    data(:,1) = foo;
    clear foo;
end

%-------------------------------------
%% Parse the data into trials
%-------------------------------------

% each trial is defined as the period between two light pulses
foo = diff(data(:,lightChannelNum)); % differentiate the light channel to find start/stop
[LightStartTimes junk] = find(foo>lightChannelONthresh); % time-points before light turns ON
[LightStopTimes junk] = find(foo<lightChannelONthresh); % time-points when light turns OFF
clear foo junk;
fprintf('\n Number of trials in this dataset: %d \n', size(LightStartTimes,1));

%-------------------------------------
%% Collect all the light stimulated spike candidates across trials
%-------------------------------------

lightStimCandidateCounter=1;
for iter = 1: LightStimSkip: size(LightStartTimes,1)
    fprintf('\n Now selecting candidates from Trial Number: %d', iter);
    currentTrialStartTime = LightStartTimes(iter)+1;
    currentLightStimWindow = data(currentTrialStartTime:LightStopTimes(iter),:);
    
    % find spike candidates in light stim window
    % peak-finding step (couple of ways to do this - default is Matlab's)
    [pksLight,locsLight] = findpeaks(currentLightStimWindow(:,1),'minpeakheight',threshold_low, 'minpeakdistance',spikeTemplateWidth);
    % *** setting MINPEAKDISTANCE properly is very helpful !! ***
    %[maxtab, mintab] = peakdet(v, delta, x) % an alternative to findpeaks
    %figure, plot(currentLightStimWindow(:,1:3)), hold on, plot(locsLight, pksLight,'r*'), hold off, title('Light Stim Window: Spike Candidates')
    
    
    % pool the detected spike candidates from light stim window
    for i=1:length(locsLight)
        if min(locsLight(i)+spikeTemplateWidth/2,length(currentLightStimWindow(:,1))) - max(locsLight(i)-spikeTemplateWidth/2,1)< spikeTemplateWidth
            continue
        else
            detectedSpikeCandidates(:,lightStimCandidateCounter) = currentLightStimWindow(max(locsLight(i)-spikeTemplateWidth/2,1): min(locsLight(i)+spikeTemplateWidth/2,length(currentLightStimWindow(:,1))),1);
        end
        lightStimCandidateCounter = lightStimCandidateCounter+1;
    end
end

figure, plot(detectedSpikeCandidates), title('Detected spike candidates (light stim. windows)')
pause(0.1)
fprintf('\n Number of spike candidates collected:%d \n', size(detectedSpikeCandidates,2))


%-------------------------------------
%% Calculate pair-wise warp distances between all candidates
%-------------------------------------

distMat = zeros(size(detectedSpikeCandidates,2), size(detectedSpikeCandidates,2));
fprintf('\n Now computing pairwise distances.. \n')
for i=1:size(detectedSpikeCandidates,2)
    for j=1: size(detectedSpikeCandidates,2)
        if strcmp(warp_type,'DTW') && size(detectedSpikeCandidates,2)<10
            [distMat(i,j), junk1 ,junk2] = dtw_WarpingDistance(detectedSpikeCandidates(:,i), detectedSpikeCandidates(:,j));
            clear junk1 junk2;
        else
            [distMat(i,j)] = procrustes(detectedSpikeCandidates(:,i), detectedSpikeCandidates(:,j), 'scaling', procrustes_scaling_flag);
        end
    end
    fprintf('%d, ',i)
end
fprintf('\n done! \n')

%-------------------------------------
%% Unsupervised grouping via affinity propagation
%-------------------------------------
apply_AP_flag =1;
fprintf('\n Applying affinity propagation.. \n')
k=1;
while apply_AP_flag==1
    sim = -distMat;
    % choose the median biasing statement as default
    [idx,netsim,dpsim,expref] = apcluster(sim, k*median(median(sim)));
    clustCenters = unique(idx);
    numCurrentClusters = length(unique(idx));
    %fprintf('\n k = %f,  numClusters = %d',k, numClusters);
    if (numCurrentClusters <= max_display_clusters) && (numCurrentClusters>=min_display_clusters)
        apply_AP_flag=0;
    elseif numCurrentClusters > max_display_clusters
        k = k*1.1;
    elseif numCurrentClusters < min_display_clusters
        k = k*0.9;
    end
    %fprintf('\n Number of current clusters = %d; Updating the clusters.. \n',numCurrentClusters)
end
fprintf('done! \n')

maxPlotLimit = max(detectedSpikeCandidates(:));
minPlotLimit = min(detectedSpikeCandidates(:));

figure;
for i=1:numCurrentClusters
    counter=1;
    for j=1:size(detectedSpikeCandidates,2)
        if idx(j)== clustCenters(i)
            foo(:,counter)= detectedSpikeCandidates(:,j);
            counter=counter+1;
        end
    end
    subplot(numCurrentClusters,1, i), plot(foo), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), title(sprintf('Cluster Template ID = %d, cluster member size = %d',clustCenters(i), size(foo,2))),...
        hold on, plot(detectedSpikeCandidates(:,clustCenters(i)),'k*'), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold off;
    clear foo;
end

%-------------------------------------
%% Get user input to choose the right spike template
%-------------------------------------

numChosenClusters = input('\n Enter the number of CLEAN cluster template IDs you want to pool together:');
for i=1:numChosenClusters
    ChosenClustCenters(i) = input('\n Enter the CLEAN template waveform ID:');
    spikeTemplates(:,i) = detectedSpikeCandidates(:,ChosenClustCenters(i));
end

%-------------------------------------
%% Find spike candidate locations across the trial period
%-------------------------------------

% Find potential spikes using a data-dependent threshold
detection_threshold = min(mean(spikeTemplates))+min(std(spikeTemplates));
fprintf('\n Spike candidate detection threshold estimated from light windows: %f', detection_threshold);

[pks,locs] = findpeaks(data(:,1),'minpeakheight',detection_threshold, 'minpeakdistance', spikeTemplateWidth/2);
clear pks;
if length(locs) < locs_max
    locs_max = length(locs);
end
fprintf('\n Total number of spikes to be annotated is %d\n', locs_max);

locs = locs(1:locs_max,:); % ML, limits the total number of spikes to be annotated
fprintf('\n Number of potential spike candidate locations: %d\n\n', length(locs));


%-------------------------------------
%% save results and processing parameters
%-------------------------------------
fprintf('\n\n\n Open the annotated file\n\n');

[fName directory fOpenFlag] = uigetfile(generalDataRootDirectory);
data = importdata([directory fName]);

if size(data,1)<size(data,2)
    data = data';
end
if size(data,2)<4
    fprintf('Annotation channel (ch-4) is missing !!');
    return
end

if size(data,2)>4
    foo = data(:,1:4); 
    clear data;
    data = foo; 
    clear foo;
end

result.data = data;
result.locs_max = locs_max;
result.detection_threshold = detection_threshold;
clear data

fprintf('\n Output is stored into a datastructure called result in the file : %s \n', [directory fName]);
eval(sprintf('save %s result', [directory fName]));

