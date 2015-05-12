function varargout = ProbabilityThresholdSelector(varargin)
% PROBABILITYTHRESHOLDSELECTOR MATLAB code for ProbabilityThresholdSelector.fig
%      PROBABILITYTHRESHOLDSELECTOR, by itself, creates a new PROBABILITYTHRESHOLDSELECTOR or raises the existing
%      singleton*.
%
%      H = PROBABILITYTHRESHOLDSELECTOR returns the handle to a new PROBABILITYTHRESHOLDSELECTOR or the handle to
%      the existing singleton*.
%
%      PROBABILITYTHRESHOLDSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROBABILITYTHRESHOLDSELECTOR.M with the given input arguments.
%
%      PROBABILITYTHRESHOLDSELECTOR('Property','Value',...) creates a new PROBABILITYTHRESHOLDSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProbabilityThresholdSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProbabilityThresholdSelector_OpeningFcn via varargin.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Modified from: J.P. Rueff, Aout 2004, modified Juin 2005
% Last Modified by Vivek, 4th Sept, 2011

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ProbabilityThresholdSelector_OpeningFcn, ...
    'gui_OutputFcn',  @ProbabilityThresholdSelector_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT
end

% --- Executes just before ProbabilityThresholdSelector is made visible.
function ProbabilityThresholdSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProbabilityThresholdSelector (see VARARGIN)

spikeData = varargin{1}; % Structure with lots of data

% Choose default output for ProbabilityThresholdSelector
HistThreshold_High = spikeData.HistThreshold_High;
HistThreshold_Low = spikeData.HistThreshold_Low;

varargout
handles.spikeData = spikeData;
handles.spikeData.data = []; % Don't pass around the huge and entire dataset

% Set up histogram plot
handles.HistThreshold_High = HistThreshold_High;
handles.HistThreshold_Low = HistThreshold_Low;

% Update handles structure for access in other functions
guidata(gca, handles);

[goodSpikes nGood unsureSpikes nUnsure badSpikes nBad] = thresholdChange(handles, 0);

handles.nGood = nGood;
handles.nBad = nBad;
handles.nUnsure = nUnsure;

% Set up initial plots when GUI isn't yet visible
if strcmp(get(hObject,'Visible'),'off')
    
    % Plot histogram
    hist(handles.histogramAxes, spikeData.spikeProbs, spikeData.numHistBins);
    title(handles.histogramAxes, 'Spike Probability Histogram');
    hold on;
    
    % Plot raw data with spike positions on top, and probability bars
    % IMPROVE THIS TO SHOW ONLY CHOSEN SUBSET OF DATA
    plot(handles.rawSpikesAxes, spikeData.data(:,:));
    hold on;
    plot(handles.rawSpikesAxes, spikeData.locs, spikeData.spikeProbs,'gx');
    title(handles.rawSpikesAxes, 'Estimated spike probabilities (Higher Probability = Better)');
    
    % Set up initial good and bad spike plots
    maxPlotLimit = max( [max(badSpikes(:)) max(unsureSpikes(:)) max(goodSpikes(:))] );
    minPlotLimit = min( [min(badSpikes(:)) min(unsureSpikes(:)) min(goodSpikes(:))] );
    
    % Display good/bad/unsure waveform sets
    plot(handles.goodSpikesAxes, goodSpikes);
    axis(handles.goodSpikesAxes, [1 size(goodSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.goodSpikesAxes, sprintf('Good Spikes (waveforms above %f): Count = %d', HistThreshold_High, nGood));
    
    plot(handles.unsureSpikesAxes, unsureSpikes);
    axis(handles.unsureSpikesAxes, [1 size(unsureSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.unsureSpikesAxes, sprintf('Unsure Spikes : Count = %d', nUnsure));
    
    plot(handles.badSpikesAxes, badSpikes);
    axis(handles.badSpikesAxes, [1 size(badSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.badSpikesAxes, sprintf('Bad Spikes (waveforms below %f): Count = %d', HistThreshold_Low, nBad));
end

axes(handles.histogramAxes);
handles.histogramLowLine = line([HistThreshold_Low HistThreshold_Low], [0 max(ylim)], 'Color', 'r');
set(handles.histogramLowLine, 'LineWidth', 2);
handles.histogramHighLine = line([HistThreshold_High HistThreshold_High], [0 max(ylim)], 'Color', 'g');
set(handles.histogramHighLine, 'LineWidth', 2);

axes(handles.rawSpikesAxes);
handles.rawSpikesLowLine = line([0 size(spikeData.data, 1)], [HistThreshold_Low HistThreshold_Low], 'Color', 'r');
handles.rawSpikesHighLine = line([0 size(spikeData.data, 1)], [HistThreshold_High HistThreshold_High], 'Color', 'g');

% Update handles structure for access in other functions
guidata(hObject, handles);

% Now make the GUI visible if the axes commands above haven't done it
set(hObject, 'Visible', 'on');

% Don't return until ready
uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = ProbabilityThresholdSelector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Returns current thresholds
varargout = {[handles.HistThreshold_High handles.HistThreshold_Low ...
    handles.nGood handles.nUnsure handles.nBad]};
delete(handles.figure1);
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

% The GUI is still in UIWAIT, use UIRESUME
uiresume(hObject);
end

% --- Executes on slider movement.
function GoodSlider_Callback(hObject, eventdata, handles)
% hObject    handle to GoodSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

hiThresh =  get(hObject, 'Value');
handles.HistThreshold_High = hiThresh;
if handles.HistThreshold_High < handles.HistThreshold_Low
    handles.HistThreshold_Low = max(handles.HistThreshold_High-0.01, 0);
    set(handles.BadSlider, 'Value', handles.HistThreshold_Low);
end

lowThresh = handles.HistThreshold_Low;
% Update lines on histogram plot
axes(handles.histogramAxes);
set(handles.histogramHighLine, 'XData', [hiThresh hiThresh]);
set(handles.histogramLowLine, 'XData', [lowThresh lowThresh]);

axes(handles.rawSpikesAxes);
set(handles.rawSpikesHighLine, 'YData', [hiThresh hiThresh]);
set(handles.rawSpikesLowLine, 'YData', [lowThresh lowThresh]);

[jnk1 nGood jnk1 nUnsure jnk3 nBad] = thresholdChange(handles, 1);

handles.nGood = nGood;
handles.nBad = nBad;
handles.nUnsure = nUnsure;

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function GoodSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GoodSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

initThresholdHigh = 0.8;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', initThresholdHigh);
end

% --- Executes on slider movement.
function BadSlider_Callback(hObject, eventdata, handles)
% hObject    handle to BadSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

lowThresh = get(hObject, 'Value');
handles.HistThreshold_Low = lowThresh;
if handles.HistThreshold_Low > handles.HistThreshold_High
    handles.HistThreshold_High = min(lowThresh+0.1,1);
    set(handles.GoodSlider, 'Value', handles.HistThreshold_High);
end

hiThresh = handles.HistThreshold_High;

% Update lines on histogram plot
axes(handles.histogramAxes);
set(handles.histogramLowLine, 'XData', [lowThresh lowThresh]);
set(handles.histogramHighLine, 'XData', [hiThresh hiThresh]);

axes(handles.rawSpikesAxes);
set(handles.rawSpikesHighLine, 'YData', [hiThresh hiThresh]);
set(handles.rawSpikesLowLine, 'YData', [lowThresh lowThresh]);

[jnk1 nGood jnk1 nUnsure jnk3 nBad] = thresholdChange(handles, 1);

handles.nGood = nGood;
handles.nBad = nBad;
handles.nUnsure = nUnsure;

guidata(hObject, handles);

end


% --- Executes during object creation, after setting all properties.
function BadSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BadSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

initThresholdLow = 0.1;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', initThresholdLow);
end

% --- Executes on button press in FinalizeButton.
function FinalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to FinalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the figure's handle, then destroy it

uiresume(handles.figure1);
end

% Function below does all the work, updating plots etc.
function [goodSpikes nGood unsureSpikes nUnsure badSpikes nBad] = ...
    thresholdChange(handles, updateDisplay)
% updateDisplay = 0: don't refresh data - no plots exist yet (dlg just opening)
% updateDisplay = 1: update data on existing plots

HistThreshold_High = handles.HistThreshold_High;
HistThreshold_Low = handles.HistThreshold_Low;

spikeProbs = handles.spikeData.spikeProbs;
curSpikeTargets = handles.spikeData.spikes;

goodSpikeIdxs = spikeProbs > HistThreshold_High;
badSpikeIdxs = spikeProbs < HistThreshold_Low;
unsureSpikeIdxs = spikeProbs <= HistThreshold_High & spikeProbs >= HistThreshold_Low;

%************************
% Display the three sets of waveforms
%************************
if ~goodSpikeIdxs
    goodSpikes = zeros(size(curSpikeTargets,1),1);
    nGood = 0;
else
    goodSpikes = curSpikeTargets(:, goodSpikeIdxs);
    nGood = size(goodSpikes, 2);
end

if ~unsureSpikeIdxs
    unsureSpikes = zeros(size(curSpikeTargets,1),1);
    nUnsure = 0;
else
    unsureSpikes = curSpikeTargets(:, unsureSpikeIdxs);
    nUnsure = size(unsureSpikes, 2);
end

if ~badSpikeIdxs
    badSpikes = zeros(size(curSpikeTargets,1),1);
    nBad = 0;
else
    badSpikes = curSpikeTargets(:, badSpikeIdxs);
    nBad = size(badSpikes,2);
end

if updateDisplay
    maxPlotLimit = max( [max(badSpikes(:)) max(unsureSpikes(:)) max(goodSpikes(:))] );
    minPlotLimit = min( [min(badSpikes(:)) min(unsureSpikes(:)) min(goodSpikes(:))] );

    plot(handles.goodSpikesAxes, goodSpikes);
    axis(handles.goodSpikesAxes, [1 size(goodSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.goodSpikesAxes, sprintf('Good Spikes (waveforms above %f): Count = %d', HistThreshold_High, nGood));
    
    plot(handles.unsureSpikesAxes, unsureSpikes);
    axis(handles.unsureSpikesAxes, [1 size(unsureSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.unsureSpikesAxes, sprintf('Unsure Spikes : Count = %d', nUnsure));
    
    plot(handles.badSpikesAxes, badSpikes);
    axis(handles.badSpikesAxes, [1 size(badSpikes,1) minPlotLimit maxPlotLimit]);
    title(handles.badSpikesAxes, sprintf('Bad Spikes (waveforms below %f): Count = %d', HistThreshold_Low, nBad));
end

end
