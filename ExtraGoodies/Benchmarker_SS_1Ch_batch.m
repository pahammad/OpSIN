%************************************************************
% Semi-supervised (SS) Benchmarking Tool for OpSIN
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
% 5. Find spike candidates in the entire session via amplitude thresholding (findpeaks function).
% 5A. If Double-bump spikes exist, apply two-step filtering (based on peak-height to prune down the spike candidates and center the
% double-bump shapes appropriately.
% 6. Estimate spike probability at candidate locations via multi-template matching.
% 7. Display histogram of spike probabilities, and ask user to set two conservative threshold(s).
% 8. Display the spike populations in a 3-panel display (good/unsure/bad spikes)
% 9. Guide the user manually through the unsure spikes to obtain user input															 
% 10. Store the finalized spike results
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
%
% 07/12/11: added conditions to handle both long & short light pulses
%
% 07/15/2011: Implemented two-threshold based approach, with 3-panel
% display - divides the spike candidates into Good/Unsure/Bad sets
%
%************************************************************


clc; clear all; close all; pause (0.1);
tic
%************************************************************
%% set up the processing parameters
%************************************************************
generalDataRootDirectory = '/Users/aljoschaschulze/Desktop/Dynamic odour_Or42a_ChR2/';


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
clear foo junk;

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
            [DBpks,DBlocs] = findpeaks(curSpikeTarget,'minpeakheight', detection_threshold, 'minpeakdistance', round(doubleBumpSeparation/2));
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
                %if max(DBpkHeights)> DB_FirstBumpHeightRatio*(maxPlotLimit-minPlotLimit)
                    MainBumpOffset = (DBspikeTemplateWidth/2+1) - DBlocs(MainBumpID);
                    adjLocs(counter) = locs(i) - MainBumpOffset;
                    counter=counter+1;
                %end
                
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
targetSpikeDist = zeros(length(pks),numChosenClusters);
minDist = zeros(length(pks),1);
for i=1:length(locs)
    if length(data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1))< spikeTemplateWidth
        continue
    else
        curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);
        for j=1:numChosenClusters
            if strcmp(warp_type,'DTW')
                [targetSpikeDist(i,j), junk1,junk2] = dtw_WarpingDistance(curSpikeTarget, spikeTemplates(:,j));
                clear junk1 junk2;
            elseif strcmp(warp_type,'DTW_Wav')
                [C1,L1] = wavedec(curSpikeTarget,3,'db1');
                [C2,L2] = wavedec(spikeTemplates(:,j),3,'db1');
                [targetSpikeDist(i,j), junk1,junk2] = dtw_WarpingDistance(C1(1:100), C2(1:100));
%             elseif strcmp(warp_type,'PCR')
%                 [targetSpikeDist(i,j)] = procrustes(spikeTemplates(:,j),curSpikeTarget, 'scaling', procrustes_scaling_flag);
%             elseif strcmp(warp_type,'PCR_Wav')
%                 [C1,L1] = wavedec(spikeTemplates(:,j),3,'db1');
%                 [C2,L2] = wavedec(curSpikeTarget,3,'db1');
%                 [targetSpikeDist(i,j)] = procrustes(C1,C2, 'scaling', procrustes_scaling_flag);
            end
            minDist(i,1) = min(targetSpikeDist(i,:));
        end
    end
end

% adaptively estimate the spikeMatchScalingParam
spikeMatchScalingParam = mean(targetSpikeDist(:));

% compute the spike probability
spikeProbs = zeros(length(pks),1);
for i=1:length(locs)
    spikeProbs(i) = exp(- (minDist(i))/spikeMatchScalingParam);
end

figure;
subplot(2,1,1), hist(spikeProbs,numHistBins); title ('Spike Probability Histogram');
subplot(2,1,2), plot(data(:,1:3)), hold on, plot(locs, spikeProbs,'go'), title ('Estimated spike probabilities (Higher Probability = Better)'), hold off

%-------------------------------------
%% threshold the chosen histogram
%-------------------------------------

update_HistThreshold_flag = 1;
HistThreshold_Low = 0.1;
HistThreshold_High = 0.8;

while update_HistThreshold_flag==1
    clear goodSpikes unsureSpikes badSpikes
    fprintf('\n\n Let us start with HistThreshold_Low=%f & HistThreshold_High=%f',HistThreshold_Low,HistThreshold_High)
    fprintf('\n\n Splitting the probability histogram into 3 regions(Bad, Unsure and Good)..')   

    i =1;
    goodSpikesCounter =0;
    badSpikesCounter=0;
    unsureSpikesCounter=0;
    while i<= length(locs)
        %************************
        %divide the spike candidates into good/unsure/bad
        %************************
        curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);
        
        if length(curSpikeTarget) >= spikeTemplateWidth
            if spikeProbs(i)> HistThreshold_High
                goodSpikesCounter = goodSpikesCounter+1;
                goodSpikes(:, goodSpikesCounter)= curSpikeTarget;
            elseif spikeProbs(i)< HistThreshold_Low
                badSpikesCounter = badSpikesCounter+1;
                badSpikes(:, badSpikesCounter)= curSpikeTarget;
            else
                unsureSpikesCounter = unsureSpikesCounter+1;
                unsureSpikes(:, unsureSpikesCounter)= curSpikeTarget;                
            end
        end
        i=i+1;
    end
    
    %************************
    % Display the three sets of waveforms
    %************************
    if ~exist('goodSpikes', 'var')
        goodSpikes = zeros(size(curSpikeTarget),1);
    end
    if ~exist('unsureSpikes', 'var')
        unsureSpikes = zeros(size(curSpikeTarget),1);
    end
    if ~exist('badSpikes', 'var')
        badSpikes = zeros(size(curSpikeTarget),1);
    end
    
    maxPlotLimit = max( [max(badSpikes(:)) max(unsureSpikes(:)) max(goodSpikes(:))] );
    minPlotLimit = min( [min(badSpikes(:)) min(unsureSpikes(:)) min(goodSpikes(:))] );
        
    % Display good/bad/unsure waveform sets
    figure;
    subplot(3,1,1), plot(goodSpikes), axis([1 size(goodSpikes,1) minPlotLimit maxPlotLimit]), title(sprintf('Good Spikes (waveforms above %f): Count = %d', HistThreshold_High, goodSpikesCounter))
    subplot(3,1,2), plot(unsureSpikes), axis([1 size(unsureSpikes,1) minPlotLimit maxPlotLimit]), title(sprintf('Unsure Spikes : Count = %d', unsureSpikesCounter))
    subplot(3,1,3), plot(badSpikes), axis([1 size(badSpikes,1) minPlotLimit maxPlotLimit]), title(sprintf('Bad Spikes (waveforms below %f): Count = %d', HistThreshold_Low, badSpikesCounter))
    
    update_HistThreshold_flag = input('\n Would you like to update your thresholds? \n Enter 1 for Yes and 0 for No:');
    if update_HistThreshold_flag == 1
        HistThreshold_Low = input('\n \n Pick a conservatve lower threshold value:');
        HistThreshold_High = input('\n \n Pick a conservatve upper threshold value:');
    end

end

%-------------------------------------
%% Proof-read the unsure spikes and finalize the spike events
%-------------------------------------

fprintf('\n...................................\n')
fprintf('\n Total spike candidates = %d [Good= %d, Unsure = %d, Bad = %d]', length(locs), goodSpikesCounter, unsureSpikesCounter, badSpikesCounter);
fprintf('\n...................................\n')
fprintf('\n Starting the proof-reading phase to finalize the "unsure" spike candidates..')
fprintf('\n Enter 11 to mark the presented waveform as high confidence good spike.')
fprintf('\n Enter 1 to mark the presented waveform as low confidence good spike.')
fprintf('\n Enter 22 to mark the presented waveform as high confidence bad spike.')
fprintf('\n Enter 2 to mark the presented waveform as low confidence bad spike.')
fprintf('\n Enter 0 to mark the presented waveform as unsure (default setting).')
fprintf('\n Enter -1 to Undo previous entry.')
fprintf('\n...................................\n')

spikeMarkers = zeros(length(locs));
for i=1:length(locs)
    if spikeProbs(i)> HistThreshold_High
        spikeMarkers(i)=1;
    elseif spikeProbs(i)< HistThreshold_Low
        spikeMarkers(i)=-1;
    else
        spikeMarkers(i)=0;
    end    
end

% only display the unsure spike events to the user and ask for feedback
validatedSpikeCounter=0;
i=0;
while i<length(locs)
        i=i+1;
    if spikeMarkers(i)==0
        validSpikeFlag =0;
        figure(2);
        if length(data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1))< spikeTemplateWidth
            continue
        else
            curSpikeTarget = data(max(locs(i)-spikeTemplateWidth/2,1): min(locs(i)+spikeTemplateWidth/2,length(data(:,1))),1);
            curSpikeTarget_long = data(max(locs(i)-10*spikeTemplateWidth/2,1): min(locs(i)+10*spikeTemplateWidth/2,length(data(:,1))),1);
            subplot(2,1, 1), plot(curSpikeTarget_long), axis([1 length(curSpikeTarget_long) minPlotLimit maxPlotLimit]), title(sprintf('Zoomed-out view of current "Unsure" waveform')),...
                hold on, plot(10*spikeTemplateWidth/2+1,curSpikeTarget_long(10*spikeTemplateWidth/2+1),'r*'), axis([1 length(curSpikeTarget_long) minPlotLimit maxPlotLimit]), hold off;
            subplot(2,1, 2), plot(goodSpikes), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold on, plot(curSpikeTarget,'r*'), axis([1 size(detectedSpikeCandidates,1) minPlotLimit maxPlotLimit]), hold off;
            validSpikeFlag = input('\n\nEnter your input (11= highConf good spike, 1= lowConf good spike, 22= highConf bad spike, 2= lowConf bad spike, 0=unsure (default), -1=undo):');
            if validSpikeFlag == 11
                spikeMarkers(i)=1;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            elseif validSpikeFlag == 1
                spikeMarkers(i)=0.5;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            elseif validSpikeFlag == 22
                spikeMarkers(i)=-1;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            elseif validSpikeFlag == 2
                spikeMarkers(i)=-0.5;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            elseif validSpikeFlag == 0
                spikeMarkers(i)=0;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            elseif validSpikeFlag == -1
                i = i-2;
                if i<1, i=1; end                 
            else 
                spikeMarkers(i)=0;
                validatedSpikeCounter = validatedSpikeCounter+1;
                fprintf('\nvalidatedSpikeCounter=%d', validatedSpikeCounter);
            end
        end                
    end
end

% finalize the benchmarked spike locations
count=0;
for i=1:length(locs)
    if spikeMarkers(i)> 0
        count=count+1;
        BenchmarkedSpikeLocs(count)=locs(i);
    end
end


%-------------------------------------
%% Compute the statistics on unsure spike percentage
%-------------------------------------

Separation_Quality_Score1 = (length(locs)-goodSpikesCounter-badSpikesCounter) / (length(locs)); 
Separation_Quality_Score2 = (length(locs)-length(find(spikeMarkers>0))-length(find(spikeMarkers<0)))/ (length(locs)); 

fprintf('\n\n Separation_Quality_Score1=%f\n', Separation_Quality_Score1);
fprintf('\n\n Separation_Quality_Score2=%f\n', Separation_Quality_Score2);


%-------------------------------------
%% Finalize and save the output data structure (result)
%-------------------------------------

fprintf('\n\n Finalized output and parameter settings stored in a data structure named: result');

% save an non-bandpass-filtered copy of original ephys data
data(:,1)= Orig_data(:,1);

% Create the result data structure and store output+parameters

% adding the benchmark annotations into channel-4
fprintf('\n\n Benchmarked spike locations are stored in channel-4 of data');
foo = zeros(size(data,1),1);
foo(BenchmarkedSpikeLocs)=1;
data(:,4) = foo;
clear foo clear Orig_data;

result.data = data;
result.detection_threshold = detection_threshold;

result.stats.Separation_Quality_Score1 = Separation_Quality_Score1;
result.stats.Separation_Quality_Score2 = Separation_Quality_Score2;
result.stats.numGoodSpikes = goodSpikesCounter;
result.stats.numBadSpikes = badSpikesCounter;
result.stats.numUnsureSpikes = (length(locs)-goodSpikesCounter-badSpikesCounter);
result.stats.numHighConfidenceGoodSpikes = length(find(spikeMarkers>0.5));
result.stats.numLowConfidenceGoodSpikes = length(find(spikeMarkers>0))-length(find(spikeMarkers>0.5));
result.stats.numHighConfidenceBadSpikes = length(find(spikeMarkers<0.5));
result.stats.numLowConfidenceBadSpikes = length(find(spikeMarkers<0))-length(find(spikeMarkers<0.5));


result.params.HistThreshold_High=HistThreshold_High;
result.params.HistThreshold_Low=HistThreshold_Low;
result.params.warp_type = warp_type;
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
fprintf('\n Output is stored into a datastructure called result in the file : %s \n', strcat(directory,strcat(fName,'_SS_BenchmarkResult.mat')));
eval(sprintf('save %s result', strcat(directory,strcat(fName,'_SS_BenchmarkResult.mat'))));
fName
toc

clear all; 

wwinput=input('\n Please press "y" in case you are OK to close all graphs (Default: N).','s');

if isempty(wwinput)|wwinput == 'n'
else
end

if wwinput == 'y'

   close all;

elseif wwinput == 'n'
end



cd('/Users/aljoschaschulze/Desktop/Dynamic odour_Or42a_ChR2/')

file = uigetfile('./*.mat');
load(file);

n = length(result.data)./20000;
x = [0.00005:0.00005:n];
nl = length(result.data);


figure;

hold on;
plot(x,result.data(:,1));
plot(x,0.005.*result.data(:,2),'r','LineWidth',2);
if size(result.data,2)==3
plot(x,result.data(:,3)/60,'k','LineWidth',2);
elseif size(result.data,2)==4
plot(x,result.data(:,3)/10000,'k','LineWidth',2);
plot(x,(0.02.*result.data(:,4))+0.05,'m','LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_PSTH_analysis_%%%%%%%%%%
clear spikes;

plpl=length(result.data);

numTrials = 1;

fSample = 20000; % Sampling frequency (in Hz)
tStep = 1/fSample;

% How many SDs away from normal peak-to-peak variation do we expect spikes to be?
% How large is the descent from peak that we expect for a real spike?
% Both used in detectSpikesWithTemplate()
sdThresh = 2;
posDeflectThresh = 1.5; % This param is ignored at the moment

acqGain = 1000;

% Raster display param (see createRaster)
spikeDisplayHt = 0.8; % How tall should the spike look in the raster

% Params for PSTH
binSize = 500/1000; % 500 msec resolution

filtH = .005*fSample; % 5 msec smoothing
HSIZE = [1 filtH];
SIGMA = 1;
H = fspecial('gaussian', HSIZE, SIGMA);

totTime = length(result.data)/fSample;
timeVec = 0:tStep:totTime-tStep;

bins = [0:binSize*fSample:length(timeVec)];

spikes=[];
spikes(1:plpl,1)=zeros;
spikes(:,1)=result.data(:,4);
spikes(:,2)=result.data(:,1);

% In case that I want to make the script run several trials  
%spikes = zeros(numTrials, length(ephys_trace));

filtPSTHPerTrial = zeros(numTrials, length(bins)-1);
for i = 1:numTrials
    if ~isempty(spikes)
        % Only if there are spikes, add them up in bins of specified size
        % after convolving with the smoothing filter
        psth = histc(find(spikes(:,1)), bins)*(1./binSize); % Avg. Inst. Firing Rate
        psthFilt = conv(psth(1:length(psth)-1), H);
        psthFiltOut = psthFilt(floor(filtH/2):(length(psthFilt)-floor(filtH/2)));
    else
        psthFiltOut = zeros(size(bins));
    end
    filtPSTHPerTrial(i, :) = psthFiltOut;
end

%figure;

% Plot stimulus presentation signal. Below, we assume there may be slight
% changes in stimulus presentation times (let me know if you want code to
% compensate for such jitter). If you're sure stimuli are delivered at the
% same time across trials, just plot the var "stim" from above (first
% trial only)
% subplot(3, 1, 1);
% for i=1:numTrials
%     plot(timeVec, stimulus(i,:)', 'k');
%     hold on
%     ylabel('Stimulus signal');
% end

% Reorder before displaying so first trial is on top and so on down.
% Should put that in createRaster
%subplot(2,1,1);
% figure;
% createRaster(spikes, spikeDisplayHt, fSample, 'k');
% xlabel('Time (s)');
% ylabel('Trial/trace number');
% xlim([0 length(timeVec)/fSample]);
% ylim([0 numTrials+(1-spikeDisplayHt)]);
% set(gca, 'YTick', [0.5:numTrials+1]);
% set(gca, 'YTickLabel', numTrials:-1:1);
% %H = title([filename]);


% Create PSTH by averaging across trials
filtPSTH = mean(filtPSTHPerTrial, 1);

figure;
hold on;
%subplot(2,1,2);
% To account for cases where total time isn't a perfect multiple of binSize
binsForPSTH = [binSize/2:binSize:binSize/2+length(timeVec)/fSample];
binsForPSTH = binsForPSTH(1:length(filtPSTH));

plot(result.data(:,1));
plot(0.005.*result.data(:,2),'r','LineWidth',2);
if size(result.data,2)==3
plot(result.data(:,3)/60,'y','LineWidth',2);
elseif size(result.data,2)==4
plot(result.data(:,3)/10000,'k','LineWidth',2);
plot((0.02.*result.data(:,4))+0.05,'m','LineWidth',2);
end

figure;
hold on;
plot(result.data(:,1));
plot(0.005.*result.data(:,2),'c','LineWidth',2);
if size(data,2)==3
plot(data(:,3)/20,'k','LineWidth',2);
elseif size(result.data,2)==4
plot(data(:,3)/20,'k','LineWidth',2);
plot(0.0005.*result.data(:,4),'r');
elseif size(result.data,2)==5
plot(data(:,3)/20,'k','LineWidth',2);
plot(0.0005.*result.data(:,4),'r');
plot(0.0005.*result.data(:,5),'y');
elseif size(result.data,2)==6
plot(data(:,3)/20,'k','LineWidth',2);
plot((0.02.*result.data(:,4)+0.05),'m');
plot(0.0005.*result.data(:,5),'r');
plot(0.0005.*result.data(:,6),'y');
end


figure;
hold on;

% Not plotting this at the moment. But if want SDs, uncomment below:
sdPSTH = std(filtPSTHPerTrial, 0, 1); 
errorbar(binsForPSTH, filtPSTH, sdPSTH, 'LineWidth', 1, 'Color', 'k');
plot(binsForPSTH, filtPSTH, 'LineWidth', 2, 'Color', 'k');

xlabel('Time (s)');
ylabel('Inst. Firing Rate (Hz)');
xlim([0 length(timeVec)/fSample]);

title('whatever');
