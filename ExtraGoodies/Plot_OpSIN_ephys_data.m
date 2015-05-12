clear all; 

cd('/Volumes/Macintosh HD/PhD/Data/Ephys/Analyzed_by_OpSIN/')

file = uigetfile('./*.mat');
load(file);

n = length(result.data)./20000;
x = [0.00005:0.00005:n];
nl = length(result.data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_PSTH_analysis_%%%%%%%%%%

plpl=length(result.data);

numTrials = 1;

fSample = 20000; % Sampling frequency (in Hz)
tStep = 1/fSample;

sdThresh = 2;
posDeflectThresh = 1.5; % This param is ignored at the moment

acqGain = 1000;

% Raster display param (see createRaster)
spikeDisplayHt = 0.8; % How tall should the spike look in the raster

% Params for PSTH
binSize = 50/1000; % 50 msec resolution

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

% In case that you want to make the script run several trials  
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


% Create PSTH by averaging across trials
filtPSTH = mean(filtPSTHPerTrial, 1);

% To account for cases where total time isn't a perfect multiple of binSize
binsForPSTH = [binSize/2:binSize:binSize/2+length(timeVec)/fSample];
binsForPSTH = binsForPSTH(1:length(filtPSTH));


figure;
hold on;
plot(result.data(:,1));
plot(0.005.*result.data(:,2),'m','LineWidth',2);
if size(result.data,2)==3
%plot(result.data(:,3)/20,'k','LineWidth',2);
elseif size(result.data,2)==4
plot(result.data(:,3)/20,'k','LineWidth',2);
plot(0.01.*result.data(:,4),'c');
end


figure;
hold on;

% Not plotting this at the moment. But if want SDs, uncomment below:
sdPSTH = std(filtPSTHPerTrial, 0, 1); 
%errorbar(binsForPSTH, filtPSTH, sdPSTH, 'LineWidth', 1, 'Color', 'k');
plot(binsForPSTH, filtPSTH, 'LineWidth', 2, 'Color', 'k');

xlabel('Time (s)');
ylabel('Inst. Firing Rate (Hz)');
xlim([0 length(timeVec)/fSample]);


