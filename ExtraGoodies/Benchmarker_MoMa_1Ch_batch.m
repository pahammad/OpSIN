%************************************************************
% Mostly Manual (MoMA) Benchmarking Tool for OpSIN - Single Channel
%************************************************************
% Parvez Ahammad, Janelia Farm, HHMI;
% parvez@ieee.org
%************************************************************
%
% ALGORITHM STEPS:
% 1. Apply bandpass filter to remove noise artifacts
% 2. Find spike candidates among the light stimulation windows
% 3. Group the candidate light spikes into clusters and ask user for input about the preferred candidate cluster(s).
% 4. Store the representatives of the preferred cluster(s) as spike template(s).
% 5. Find spike candidates in the entire session via thresholding (findpeaks function).
% 6. Display finalized results sequentially, for proof-reading
% 7. Store the finalized spike results
%************************************************************
%
% EXPECTED INPUT DATA FORMAT: T x 3 matrix
% T = number of time points recorded
% first channel: electrode data
% second channel: odor
% third channel: light
%************************************************************
%
% Comments/Feedback:
% Matthieu Louis (CRG, Barcelona)
%************************************************************
%
% UPDATE NOTES:
% 06/01/2011: first version of benchmarking tool completed
% 06/03/2011: 1-step undo added, larger context window added
% 06/05/2011: modification by ML to solve 2 bugs and limit the number of
% spikes to be annotated
%************************************************************

clc; clear all; close all; pause (0.1);
tic
%************************************************************
%% set up the processing parameters
%************************************************************
generalDataRootDirectory= '/Users/ahammadp/Documents/parvez-laptop/work/data/';


% this parameter determines the acceptable first bump height in double-bump filtering scenario (default: 0.2 or 0.1)
DB_FirstBumpHeightRatio = 0.1;
% min. peak height for peak-finding step
threshold_low = 0.01;
%@20KHz ; default: 2 times the length of expected spike
spikeTemplateWidth = 180;
% Decides how many light-stimulus windows to skip for sampling the spike templates - useful for speeding up unsupervised clustering part
LightStimSkip = 1;
%@20KHz ; if OpSIN detection and User annotation are within these many time-points, we declare a win!
spikeDetectionAcceptanceRange = 40;
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
lightChannelONlevel = 1;
lightPulseWidth = 5000;
% warp type used: default is DTW for the template-matching phase
warp_type = 'DTW';

dataChannelNum = 1;
% max number of clusters after clustering the light-window spike candidates
max_display_clusters = 8;
% min number of clusters after clustering the light-window spike candidates
min_display_clusters = 3;
numHistBins = 20;
procrustes_scaling_flag = 0;

locs_max = 1000; % maximum number of spikes to be manually annotated %ML

%************************************************************

%-------------------------------------
%% set up the processing module flags
%-------------------------------------

[fName directory fOpenFlag] = uigetfile(generalDataRootDirectory);

if fOpenFlag==0
    fprintf('\n File selection failed !!\n');
    return;
end
fprintf('\n\n Running OpSIN spike-detection tool:\n');
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

foo = importdata([directory fName]);
if isstruct(foo)
    data = foo.data;
else
    data= foo;
end

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

Orig_data = data;
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
foo2 = data(:,lightChannelNum);
foo2(find(foo2>lightChannelONlevel))=lightChannelONthresh;
foo2(find(foo2<lightChannelONlevel))=lightChannelONthresh;
foo = diff(foo2); % differentiate the light channel to find start/stop
[LightStartTimes junk] = find(foo>lightChannelONthresh); % time-points before light turns ON
[LightStopTimes junk] = find(foo<lightChannelONthresh); % time-points when light turns OFF
firstPulseWidth = abs(LightStopTimes(1)-LightStartTimes(1));
clear foo foo2 junk;

% check for condition where the light is being pulsed
if firstPulseWidth/F_s < 200/F_s
    fprintf('\n Cleared the LightStopTimes .. \n')
    clear LightStopTimes
    LightStopTimes = LightStartTimes+lightPulseWidth; % time-points when light turns OFF
end

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
    if firstPulseWidth/F_s < 400/F_s
        [pksLight,locsLight] = findpeaks(currentLightStimWindow(:,1),'minpeakheight',threshold_low, 'minpeakdistance',spikeTemplateWidth/20);
    else
        [pksLight,locsLight] = findpeaks(currentLightStimWindow(:,1),'minpeakheight',threshold_low, 'minpeakdistance',spikeTemplateWidth/2);      
    end
    
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
        if size(detectedSpikeCandidates,2)<30
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
        hold on, plot(detectedSpikeCandidates(:,clustCenters(i)),'r*'), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold off;
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
%% Present each candidate to the user and store the annotated result
%-------------------------------------
close all;

figure(1);
for i=1:numChosenClusters
    counter=1;
    for j=1:size(detectedSpikeCandidates,2)
        if idx(j)== ChosenClustCenters(i)
            foo(:,counter)= detectedSpikeCandidates(:,j);
            counter=counter+1;
        end
    end
    ChosenClusters{i}=foo;
    subplot(numChosenClusters,1, i), plot(foo), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), title(sprintf('Selected Cluster Template ID = %d, cluster member size = %d',ChosenClustCenters(i), size(foo,2))),...
        hold on, plot(detectedSpikeCandidates(:,ChosenClustCenters(i)),'r*-'), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold off;
    clear foo;
end

benchmarkedSpikesCounter=0;
fprintf('\n...................................\n')
fprintf('\n Starting MoMa (Mostly Manual) benchmarking process..')
fprintf('\n Enter 1 to mark the current waveform as valid spike.')
fprintf('\n Enter -1 for Undo previous entry.')
fprintf('\n Note: Just hitting Enter means current waveform is NOT a valid spike.')
fprintf('\n...................................\n')

i=0;
%%
while i<length(locs) %error corrected by ML
    i=i+1;
    validSpikeFlag =0;
    figure(2);
    if length(data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1))< spikeTemplateWidth
        continue
    else
        curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);
        curSpikeTarget_long = data(max(locs(i)-10*spikeTemplateWidth/2,1): min(locs(i)+10*spikeTemplateWidth/2,length(data(:,1))),1);
        subplot(numChosenClusters+1,1, 1), plot(curSpikeTarget_long), axis([1 length(curSpikeTarget_long) minPlotLimit maxPlotLimit]), title(sprintf('Zoomed-out view of current candidate waveform')),...
            hold on, plot(10*spikeTemplateWidth/2+1,curSpikeTarget_long(10*spikeTemplateWidth/2+1),'r*'), axis([1 length(curSpikeTarget_long) minPlotLimit maxPlotLimit]), hold off;
        for j=1:numChosenClusters
            subplot(numChosenClusters+1,1, j+1), plot(ChosenClusters{j}), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), title(sprintf('Current Candidate Waveform #%d (Red *) + Template Cluster %d',i,j)),...
                hold on, plot(curSpikeTarget,'r*'), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold off;
        end
        validSpikeFlag = input(sprintf('#%d of %d -- Enter your input (1=valid spike, -1=undo):',i,length(locs)));
        if ~isempty(validSpikeFlag)
            if ~isnumeric(validSpikeFlag)
                i=i-1;
                continue;
            end
        end
        if validSpikeFlag == 1
            benchmarkedSpikesCounter = benchmarkedSpikesCounter+1;
            BenchmarkedSpikeLocs(benchmarkedSpikesCounter)=locs(i);
        elseif validSpikeFlag == -1
            i = i-2;
            benchmarkedSpikesCounter = benchmarkedSpikesCounter-1;
        end       
    end
    fprintf('Number of spikes selected so far: %d\n\n',benchmarkedSpikesCounter)
end

%-------------------------------------
%% save results and processing parameters
%-------------------------------------

fprintf('\n\n Finalized output and parameter settings stored in a data structure named: result');

% save an non-bandpass-filtered copy of original ephys data
data(:,1)= Orig_data(:,1);
clear Orig_data;

% adding the benchmark annotations into channel-4
fprintf('\n\n Benchmarked spike locations are stored in channel-4 of data');
foo = zeros(size(data,1),1); % error correcred by ML
foo(BenchmarkedSpikeLocs)=1;
data(:,4) = foo;
clear foo;

result.data = data;
result.locs_max = locs_max;
result.detection_threshold = detection_threshold;
clear data

dt = date;
fprintf('\n Output is stored into a datastructure called result in the file : %s \n', strcat(directory,strcat(fName,'_MoMa_BenchmarkResult.mat')));
eval(sprintf('save %s result', strcat(directory,strcat(fName,'_MoMa_BenchmarkResult.mat'))));
fName
toc
