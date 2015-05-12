%************************************************************
% Code to compare different ground-truths for spike events
%************************************************************
% Parvez Ahammad, Janelia Farm, HHMI;
% parvez@ieee.org
%************************************************************
%
% ALGORITHM STEPS:
% 1. Load the main data/recording file
% 2. Ask user-input about how many GT files exist
% 3. Interactively load each of the GT files
% 4. Compute precision/recall/F-Scores across each pair of GTs
% 5. Compute overall statistics
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
    result{i} = importdata([directory{i} fName{i}]);
end

%-------------------------------------
%% Compute comparative statistics between various GTs
%-------------------------------------
Precision = zeros(numGT,numGT);
Recall = zeros(numGT,numGT);
FScore = zeros(numGT,numGT);

counter=0;
for GTi=1:numGT
    for GTj=1:numGT
        if GTj==GTi
            continue
        else
            counter=counter+1;
            % code to compare ground-truths
            fprintf('\n Taking the user annotation %d as baseline.. \n',GTi)
            
            [junk Locs1] = findpeaks(result{GTi}.data(:,4));
            [junk Locs2] = findpeaks(result{GTj}.data(:,4));
            clear junk;
            
            % find the agreement
            SpikesInAgreementCounter=0;
            for i=1:length(Locs1)
                SpikesInAgreementCounter = SpikesInAgreementCounter+numel(find(Locs1(i)==Locs2));
            end
            % compute precision, recall and F-score
            Precision(GTi,GTj) = SpikesInAgreementCounter/length(Locs1);
            PrVector(counter) = SpikesInAgreementCounter/length(Locs1);
            
            Recall(GTi,GTj) = SpikesInAgreementCounter/length(Locs2);
            RecVector(counter) = SpikesInAgreementCounter/length(Locs2);
            
            FScore(GTi,GTj) = 2* Precision(GTi,GTj) * Recall(GTi,GTj) / (Precision(GTi,GTj)+Recall(GTi,GTj));
            FScVector(counter) = 2* Precision(GTi,GTj) * Recall(GTi,GTj) / (Precision(GTi,GTj)+Recall(GTi,GTj));
        end
    end
end

%-------------------------------------
%% Compute overall statistics
%-------------------------------------

GTstats.Precision = Precision;
GTstats.Mean_Precision = mean(PrVector);
GTstats.Std_Precision = std(PrVector);

GTstats.Recall = Recall;
GTstats.Mean_Recall = mean(RecVector);
GTstats.Std_Recall = std(RecVector);

GTstats.FScore = FScore;
GTstats.Mean_FScore = mean(FScVector);
GTstats.Std_FScore = std(FScVector);

dataFName
GTstats
fprintf('\n GT comparison statistics are stored in the file : %s \n', [dir dataFName 'GTstats.mat']);
eval(sprintf('save %s GTstats', [dir dataFName 'GTstats.mat']));
