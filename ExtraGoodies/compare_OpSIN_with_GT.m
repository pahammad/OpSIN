%************************************************************
% Code to compare OpSIN performance with a ground-truth
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
% 5A. If Double-bump spikes exist, apply two-step filtering (based on peak-height to prune down the spike candidates and center the
% double-bump shapes appropriately.
% 6. Estimate spike probability at candidate locations via multi-template matching.
% 7. Display histogram of spike probabilities, and ask user to set appropriate threshold(s).
% 8. Display the aligned spike populations divided separated by the threshold(s)
% 9. If manual annotation is available, compute statistics of correct detection
%************************************************************
%
% EXPECTED INPUT DATA FORMAT: T x 4 matrix
% T = number of time points recorded
% first channel: electrode data
% second channel: odor
% third channel: light
%************************************************************
%
% Comments/Feedback:
% David Jarriault
% Matthieu Louis (CRG, Spain)
% Vivek Jayaraman (Janelia Farm, HHMI)
%************************************************************
%
% OPEN QUESTIONS:
% 1. What is the "correct" form of the probability function to convert warp
% distances into probability estimates?
% 2. How to handle the double bump spikes in a systematic way without doing
% a hacky feature check?
%************************************************************
%
% UPDATE NOTES:
% 03/24/2011: implemented batch processing mode
% 04/05/2011: Found the bug in parametric warping function and fixed it
% 04/06/2011: Added the Procrustes distance into the framework
% 04/06/2011: tried to fix the issue with mis-localization of spikes on
% double bump shapes, but not sure if it's working well yet.
%
% 04/07/2011: batch mode ver1 released to David/Matthieu.
% 04/27/2011: batch mode ver2 (lot of small changes to incorporate feedback, 
% including allowing two different types of histograms to choose from)
%
% 04/28/2011: convention is to declare spike location on first bump in
% double bump spikes (based on David's e-mail)
%
% 05/09/2011: Fixed the plot axes, updated the spike prob. estimates to
% compute the scaling parameter adaptively
%
% 05/17/11: double-bump spikes (OSN) are handled via a two-step filtering
% approach where all the local peaks below a certain height are discarded,
% and the peak locations of rest of the spike candidates are adjusted
% appropriately.
% 
% 05/17/11: Added new plots to deal with situations where user-annotation
% does not exist
%
% 06/20/11: Major re-write to change the result data structures
%************************************************************


%clc; 
clear all; close all; pause (0.1);
tic
%************************************************************
%% set up the processing parameters
%************************************************************
%generalDataRootDirectory = '/Users/ahammadp/Documents/parvez-laptop/work/data/SpikeSorting/OpSIN_BenchmarkDatasets/David_Noisy/';
generalDataRootDirectory = '/Users/ahammadp/Documents/parvez-laptop/work/data/';


warp_type = 'PCR';
% this parameter determines the acceptable first bump height in double-bump filtering scenario (default: 0.2 or 0.1)
DB_FirstBumpHeightRatio = 0.05; 
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
dataChannelNum = 1;
% max number of clusters after clustering the light-window spike candidates
max_display_clusters = 8; 
% min number of clusters after clustering the light-window spike candidates
min_display_clusters = 3; 
numHistBins = 20;
procrustes_scaling_flag = 0;
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

%Process_Flags_OK = input('\nEnter 1 to continue:');
%if Process_Flags_OK ~= 1
%    return;
%end

%-------------------------------------
%% set up the data input
%-------------------------------------

foo = importdata([directory fName]);

data = foo.data;
detection_threshold = foo.detection_threshold;
locs_max = foo.locs_max;
clear foo;
if size(data,1)<size(data,2)
    data = data';
end
% formatting the input data to follow the structure specified above:
if size(data,2)>4
    foo = data(:,1:4);
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
%detection_threshold = min(mean(spikeTemplates))+min(std(spikeTemplates));
fprintf('\n Spike candidate detection threshold estimated from light windows: %f', detection_threshold);

[pks,locs] = findpeaks(data(:,1),'minpeakheight',detection_threshold, 'minpeakdistance', spikeTemplateWidth/2);
if length(locs) < locs_max
    locs_max = length(locs);
end
fprintf('\n Total number of spikes to be annotated is %d\n', locs_max);

locs = locs(1:locs_max,:); % ML, limits the total number of spikes to be annotated
pks = pks(1:locs_max,:); 
fprintf('\n Number of potential spike candidate locations: %d\n\n', length(locs));


%-------------------------------------
%% Adjust the peak location for double bump spikes (if they exist)
%-------------------------------------

doubleBumpAdjust_flag = input('\n\n Enter 1 if the data contains double bump spikes, 0 otherwise:');
if (doubleBumpAdjust_flag ~= 0) && (doubleBumpAdjust_flag~=1)
    fprintf('\n Wrong Entry! \n')
    doubleBumpAdjust_flag = input('\n Enter 1 if the data contains double bump spikes, 0 otherwise:');
end

if doubleBumpAdjust_flag ==1
    doubleBumpSeparation = round(spikeTemplateWidth/6);
    DBspikeTemplateWidth = max(3*doubleBumpSeparation, 2*round(spikeTemplateWidth/2));
    %doubleBumpSeparation = input('\n Enter approx. separation between double bumps (Suggested default @ 20Khz = 25):');
end

counter=1;
if doubleBumpAdjust_flag == 1
    for i=1:length(locs)
        if length(data(max(locs(i)-DBspikeTemplateWidth/2,1): min(locs(i)+DBspikeTemplateWidth/2,length(data(:,1))),1))< DBspikeTemplateWidth
            continue
        else
            curSpikeTarget = data(max(locs(i)-DBspikeTemplateWidth/2,1): min(locs(i)+DBspikeTemplateWidth/2,length(data(:,1))),1);
            [DBpks,DBlocs] = findpeaks(curSpikeTarget,'minpeakheight', detection_threshold, 'minpeakdistance', doubleBumpSeparation/2);
            if length(DBlocs) > 1
                [junk DBvalleyLocs] = findpeaks(-curSpikeTarget);
                for p=1:length(DBlocs)
                    % trying to find the valleys sitting around the DBpeak
                    [valleyID junk] = find((DBlocs(p)-DBvalleyLocs)<0); % first entry contains the valley on RHS of peak
                    if numel(valleyID)==0
                        DBpkHeights(p) = curSpikeTarget(DBlocs(p))- curSpikeTarget(DBvalleyLocs(end));
                    elseif valleyID(1)==1
                        DBpkHeights(p)= curSpikeTarget(DBlocs(p))- curSpikeTarget(DBvalleyLocs(1));
                    else
                        DBpkHeights(p) = min( curSpikeTarget(DBlocs(p))- curSpikeTarget(DBvalleyLocs(valleyID(1))), curSpikeTarget(DBlocs(p))- curSpikeTarget(DBvalleyLocs(valleyID(1)-1)) );
                    end
                end
                
                [junk MainBumpID] = max(DBpkHeights);
                if max(DBpkHeights)> DB_FirstBumpHeightRatio*(maxPlotLimit-minPlotLimit)
                    MainBumpOffset = (DBspikeTemplateWidth/2+1) - DBlocs(MainBumpID);
                    adjLocs(counter) = locs(i) - MainBumpOffset;
                    counter=counter+1;
                end

                clear DBpks DBlocs DBpkHeights DBpkSlopes MainBumpID junk
            end
        end
    end
    % adjust the peak values appropriately
    
    length(locs)   
    length(adjLocs)
    
    clear locs
    locs = adjLocs;
    foo = data(:,1);
    pks = foo(locs);
    clear foo adjLocs;   
end

%-------------------------------------
%% Compute the spike probability and target spike distance
%-------------------------------------

% closing all figures upto this point.
close all; pause(0.1)

% compute the target spike distance
targetSpikeDist = zeros(length(locs),numChosenClusters);
minDist = zeros(length(locs),1);
for i=1:length(locs)
    if length(data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1))< spikeTemplateWidth
        continue
    else
        curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);
        for j=1:numChosenClusters
            if strcmp(warp_type,'DTW')
                [targetSpikeDist(i,j), junk1,junk2] = dtw_WarpingDistance(curSpikeTarget, spikeTemplates(:,j));
                clear junk1 junk2;
            elseif strcmp(warp_type,'PCR')
                [targetSpikeDist(i,j)] = procrustes(spikeTemplates(:,j),curSpikeTarget, 'scaling', procrustes_scaling_flag);
            end
            minDist(i,1) = min(targetSpikeDist(i,:));
        end
    end
end

% adaptively estimate the spikeMatchScalingParam
spikeMatchScalingParam = mean(targetSpikeDist(:));

% compute the spike probability
spikeProbs = zeros(length(locs),1);
for i=1:length(locs)
    spikeProbs(i) = exp(- (minDist(i))/spikeMatchScalingParam);
end

figure;
subplot(2,1,1), hist(spikeProbs,numHistBins); title ('Spike Probability Histogram');
subplot(2,1,2), plot(data(:,1:3)), hold on, plot(locs, spikeProbs,'go'), title ('Estimated spike probabilities (Higher Probability = Better)'), hold off

% Normalized distance from target spike template
normTargetSpikeDist = min(targetSpikeDist,[],2) / max(min(targetSpikeDist,[],2));
figure;
subplot(2,1,1), hist(normTargetSpikeDist,numHistBins); title ('Spike Distance Histogram');
subplot(2,1,2), plot(data(:,1:3)), hold on, plot(locs, normTargetSpikeDist ,'go'), title ('Estimated spike distances (Lower Distance = Better)'), hold off

clear targetSpikeDist;

%-------------------------------------
%% choose the appropriate histogram to threshold
%-------------------------------------

fprintf('\n \n Which of the two histogram do you prefer to process?')
fprintf('\n Note-1: Pick the one with better separability between peaks')
fprintf('\n Note-2: Pick the one where most spikes in light windows (Red Rectangles) qualify')

hist_choice = input('\n Enter 1 for Spike Prob. Histogram and 2 for Spike Distance Histogram:');
if (hist_choice ~= 1) && (hist_choice~=2)
    fprintf('\n Wrong Entry! \n')
    hist_choice = input('\n Enter 1 for Spike Prob. Histogram and 2 for Spike Distance Histogram:');
end
    
% Re-plot the chosen histogram and close the other figures
close all; pause(0.1);
figure;
if hist_choice == 1
    ChosenMeasure = spikeProbs;
    subplot(2,1,1), hist(spikeProbs,numHistBins); title ('Spike Probability Histogram');
    subplot(2,1,2), plot(data(:,1)), hold on, plot(data(:,2),'g-'), hold on, plot(data(:,3),'r-'), hold on, ...
        plot(locs, spikeProbs,'go'), title ('Estimated spike probabilities (Higher Probability = Better)'), hold off
elseif hist_choice == 2
    ChosenMeasure = normTargetSpikeDist;
    subplot(2,1,1), hist(normTargetSpikeDist,numHistBins); title ('Spike Distance Histogram');
    subplot(2,1,2), plot(data(:,1)), hold on, plot(data(:,2),'g-'), hold on, plot(data(:,3),'r-'), hold on, ...
        plot(locs, normTargetSpikeDist ,'go'), title ('Estimated spike distances (Lower Distance = Better)'), hold off
end

%-------------------------------------
%% threshold the chosen histogram
%-------------------------------------

accept_HistThreshold_flag = 0;
while accept_HistThreshold_flag==0
    clear spikesAboveThresh spikesBelowThresh
    HistThreshold = input('\n \n Enter a threshold value to split the histogram:');
    i =1;
    close all; pause (0.1)
    aboveThreshCounter =1;
    belowThreshCounter=1;
    while i<= length(locs)
        %************************
        %divide the spike candidates above/below the threshold
        %************************
        curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);

        if length(curSpikeTarget) >= spikeTemplateWidth
            if ChosenMeasure(i)> HistThreshold
                spikesAboveThresh(:, aboveThreshCounter)= curSpikeTarget;
                aboveThreshCounter = aboveThreshCounter+1;
            else
                spikesBelowThresh(:, belowThreshCounter)= curSpikeTarget;
                belowThreshCounter = belowThreshCounter+1;
            end
        end
        i=i+1;
    end
    % adjust the counters
    belowThreshCounter = belowThreshCounter-1;
    aboveThreshCounter = aboveThreshCounter-1;
    
    %************************
    % Display the thresholded & aligned waveforms
    %************************
    if ~exist('spikesAboveThresh', 'var')
        spikesAboveThresh = zeros(size(curSpikeTarget),1);
    end
    if ~exist('spikesBelowThresh', 'var')
        spikesBelowThresh = zeros(size(curSpikeTarget),1);
    end
    maxPlotLimit = max( max(spikesBelowThresh(:)) , max(spikesAboveThresh(:)) );
    minPlotLimit = min( min(spikesBelowThresh(:)) , min(spikesAboveThresh(:)) );

    figure;
    if hist_choice == 1
        ChosenMeasure = spikeProbs;
        subplot(2,1,1), hist(spikeProbs,numHistBins); title ('Spike Probability Histogram');
        subplot(2,1,2), plot(data(:,1:3)), hold on, plot(locs, spikeProbs,'go'), title ('Estimated spike probabilities (Higher Probability = Better)'), hold off
    elseif hist_choice == 2
        ChosenMeasure = normTargetSpikeDist;
        subplot(2,1,1), hist(normTargetSpikeDist,numHistBins); title ('Spike Distance Histogram');
        subplot(2,1,2), plot(data(:,1:3)), hold on, plot(locs, normTargetSpikeDist ,'go'), title ('Estimated spike distances (Lower Distance = Better)'), hold off
    end
    
    figure;
    subplot(2,1,1), plot(spikesAboveThresh), axis([1 size(spikesAboveThresh,1) minPlotLimit maxPlotLimit]), title(sprintf('Candidate waveforms above %f: number = %d', HistThreshold, size(spikesAboveThresh,2)))
    subplot(2,1,2), plot(spikesBelowThresh), axis([1 size(spikesBelowThresh,1) minPlotLimit maxPlotLimit]), title(sprintf('Candidate waveforms below %f: number = %d', HistThreshold, size(spikesBelowThresh,2)))
       
    accept_HistThreshold_flag = input('\n Are you happy with the spike output? \n Enter 1 for Yes and 0 for No:');
end


%-------------------------------------
%% Compute the finalized spikes and correct detection statistics
%-------------------------------------

% User-annotated spikes
annotatedSpikeLocs = find(data(:,4)>=1);
foo = data(:,1);
annotatedSpikePks = foo(annotatedSpikeLocs);
displayOffset = 0.05*max(data(:,1));
clear foo;

% OpSIN's spike results
if hist_choice ==1 % Prob. Histogram
    ind = find(ChosenMeasure>HistThreshold);
    figure; errorbar(1:spikeTemplateWidth+1,mean(spikesAboveThresh,2)', std(spikesAboveThresh,[],2)'), title('Point-wise Mean and Std. Dev. of Finalized Spikes')
elseif hist_choice==2 % Distance Histogram
    ind = find(ChosenMeasure<=HistThreshold);
    figure; errorbar(1:spikeTemplateWidth+1,mean(spikesBelowThresh,2)', std(spikesBelowThresh,[],2)'), title('Point-wise Mean and Std. Dev. of Finalized Spikes')
end
finalSpikeLocs = locs(ind);
finalSpikePks = pks(ind);

% Adjust the user annotated spike locations if they are within user-defined
% range from OpSIN detected spikes
adjustedAnnotatedSpikeLocs = zeros(size(annotatedSpikeLocs));
for i=1: length(annotatedSpikeLocs)
    [minTimeOffset closestOpSINSPikeLoc] = min(abs(finalSpikeLocs - annotatedSpikeLocs(i)*ones(size(finalSpikeLocs))));   
    if minTimeOffset <= spikeDetectionAcceptanceRange
        adjustedAnnotatedSpikeLocs(i) = finalSpikeLocs(closestOpSINSPikeLoc);
    else
        adjustedAnnotatedSpikeLocs(i) = annotatedSpikeLocs(i);
    end
end
foo = data(:,1);
adjustedAnnotatedSpikePks = foo(adjustedAnnotatedSpikeLocs);
clear foo;

%-------------------------------------
%% Calculate precision and recall (assuming that user annotation is correct)
%-------------------------------------

fprintf('\n\n Assuming that user annotation is correct, \n')

SpikesInAgreementCounter=0;
for i=1:length(adjustedAnnotatedSpikeLocs)
    SpikesInAgreementCounter = SpikesInAgreementCounter+numel(find(adjustedAnnotatedSpikeLocs(i)==finalSpikeLocs));
end

if length(adjustedAnnotatedSpikeLocs)>0
    Precision = SpikesInAgreementCounter/length(adjustedAnnotatedSpikeLocs);
    Recall = SpikesInAgreementCounter/length(finalSpikeLocs);
    Fscore = 2*Precision*Recall/(Precision+Recall);
    
    fprintf('Precision (defined as: #SpikesInAgreement/#UserAnnotatedSpikes) = %f\n',Precision);
    fprintf('Recall (defined as: #SpikesInAgreement/#OpSINSpikes) = %f\n', Recall);
    fprintf('Fscore (defined as: 2*Precision*Recall/(Precision+Recall)) = %f\n', Fscore);
else
    fprintf('User did not annotate spike locations in this dataset.\n');
end


%-------------------------------------
%% Finalize and save the output data structure (result)
%-------------------------------------

fprintf('\n\n Finalized output and parameter settings stored in a data structure named: result');

% reload the original data so that the band-pass filtering effects can be removed.
clear foo;
foo = importdata([directory fName]);
data = foo.data;
clear foo;
if size(data,1)<size(data,2)
    data = data';
end
if size(data,2)>4
    foo = data(:,1:4); 
    clear data;
    data = foo; 
    clear foo;
end

% Create the result data structure and store output+parameters
result.data = data;
result.data(:,1)= Orig_data(:,1);
result.OpSINSpikeLocs = finalSpikeLocs';
result.detection_threshold = detection_threshold;

result.stats.Precision = Precision;
result.stats.Recall = Recall;
result.stats.Fscore = Fscore;

result.params.HistThreshold=HistThreshold;
result.params.warp_type = warp_type;
result.params.doubleBumpAdjust_flag = doubleBumpAdjust_flag;
result.params.DB_FirstBumpHeightRatio = DB_FirstBumpHeightRatio; 
result.params.threshold_low = threshold_low; 
result.params.spikeTemplateWidth = spikeTemplateWidth; 
result.params.LightStimSkip = LightStimSkip; 
result.params.spikeDetectionAcceptanceRange = spikeDetectionAcceptanceRange; 
result.params.bandpass_prefilter_flag = bandpass_prefilter_flag; 
result.params.F_s = F_s; 
result.params.low_cutoff_freq = low_cutoff_freq; 
result.params.high_cutoff_freq = high_cutoff_freq; 
result.params.lightChannelNum = lightChannelNum;
result.params.lightChannelONthresh = lightChannelONthresh;
result.params.dataChannelNum = dataChannelNum;
result.params.max_display_clusters = max_display_clusters; 
result.params.min_display_clusters = min_display_clusters; 
result.params.numHistBins = numHistBins;
result.params.procrustes_scaling_flag = procrustes_scaling_flag;

dt = date;
eval(sprintf('save %s result', [directory fName '_OpSINResult_' dt '.mat']));
fprintf('\n Output is written to file: %s \n', [directory fName '_OpSINResult_' dt '.mat']);
fName
%clear all; 
toc


