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

%Process_Flags_OK = input('\nEnter 1 to continue:');
%if Process_Flags_OK ~= 1
%    return;
%end

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
%% save results and processing parameters
%-------------------------------------
fprintf('\n\n\n Open the annotated file\n\n');

[fName directory fOpenFlag] = uigetfile(generalDataRootDirectory);
foo = importdata([directory fName]);
data = foo.data;
locs_max = foo.locs_max;
detection_threshold = foo.detection_threshold;

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

data(:,1) = origData(:,1);
result.data = data;
result.locs_max = locs_max;
result.detection_threshold = detection_threshold;

fprintf('\n Output is stored into a datastructure called result in the file : %s \n', [directory fName]);
eval(sprintf('save %s result', [directory fName]));

