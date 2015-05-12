%************************************************************
% Code to make a consensus ground-truth by taking majority vote among GTs
%************************************************************
% Parvez Ahammad, Janelia Farm, HHMI;
% parvez@ieee.org
%************************************************************
%
% ALGORITHM STEPS:
% 1. Load the main data/recording file
% 2. Ask user-input about how many GT files exist
% 3. Interactively load each of the GT files
% 4. Compute the consensus spike annotations via majority vote (if above
% 50% benchmarkers agree, it's a consensus spike)
% 5. Store the results
%************************************************************

clc; clear all; close all;

%-------------------------------------
%% Load all the GT files
%-------------------------------------

generalDataRootDirectory = '/Users/ahammadp/Documents/parvez-laptop/work/data/SpikeSorting/';

% Get user input about how many GT files exist
numGT = input('\n How many ground-truths exist for each recording?');

% Get the GT filenames
fprintf('\n Select the original data file:');
[dataFName dir junk] = uigetfile(generalDataRootDirectory);

dataFName

for i=1:numGT
    fprintf('\n Select the ground-truth(GT) file number %d',i);
    [fName{i} directory{i} fOpenFlag] = uigetfile(dir);
    fName{i}
    res{i} = importdata([directory{i} fName{i}]);
    foo(:,i) = res{i}.data(:,4);
end
bar = mean(foo,2);
bar(find(bar>=0.5))=1;
bar(find(bar<0.5))=0;


%-------------------------------------
%% Compute the consensus GT
%-------------------------------------
result.locs_max = res{1}.locs_max;
result.detection_threshold = res{1}.detection_threshold;
result.data(:,1:3) = res{1}.data(:,1:3);
result.data(:,4)= bar;

%-------------------------------------
%% Save the consensus GT
%-------------------------------------

dataFName
fprintf('\n GT comparison statistics are stored in the file : %s \n', [dir dataFName '_BenchmarkResult_consensus.mat']);
eval(sprintf('save %s result', [dir dataFName '_BenchmarkResult_consensus.mat']));
