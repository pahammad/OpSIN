function varargout = DetectionThresholdSelector(varargin)
% DETECTIONTHRESHOLDSELECTOR MATLAB code for DetectionThresholdSelector.fig
%      DETECTIONTHRESHOLDSELECTOR, by itself, creates a new DETECTIONTHRESHOLDSELECTOR or raises the existing
%      singleton*.
%
%      H = DETECTIONTHRESHOLDSELECTOR returns the handle to a new DETECTIONTHRESHOLDSELECTOR or the handle to
%      the existing singleton*.
%
%      DETECTIONTHRESHOLDSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTIONTHRESHOLDSELECTOR.M with the given input arguments.
%
%      DETECTIONTHRESHOLDSELECTOR('Property','Value',...) creates a new DETECTIONTHRESHOLDSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DetectionThresholdSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DetectionThresholdSelector_OpeningFcn via varargin.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Modified from: J.P. Rueff, Aout 2004, modified Juin 2005
% Last Modified by Vivek, 4th Sept, 2011

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DetectionThresholdSelector_OpeningFcn, ...
    'gui_OutputFcn',  @DetectionThresholdSelector_OutputFcn, ...
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

% --- Executes just before DetectionThresholdSelector is made visible.
function DetectionThresholdSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DetectionThresholdSelector (see VARARGIN)

chosenSpikes = varargin{1};
clustCenter = varargin{2};
spikeTemplate = varargin{3};
minTemplates = varargin{4};
minSDTemplates = varargin{5};
init_threshold = varargin{6};

% Choose default output for DetectionThresholdSelector
varargout = init_threshold;

% Set up histogram plot
handles.minTemplates = minTemplates;
handles.minSDTemplates = minSDTemplates;
handles.detectionThreshold = init_threshold;
handles.thresh = 0.9; % Initial value

% Update handles structure for access in other functions
guidata(gca, handles);

% Set up initial plots when GUI isn't yet visible
if strcmp(get(hObject,'Visible'),'off')
    maxPlotLimit = max(chosenSpikes(:));
    minPlotLimit = min(chosenSpikes(:));

    plot(chosenSpikes);
    hold on;
    plot(spikeTemplate,'r*');
    handles.threshLine = line([0 length(spikeTemplate)], [init_threshold init_threshold]);
    set(handles.threshLine, 'Color' ,'k');
    axis([1 length(spikeTemplate) minPlotLimit maxPlotLimit]);
    title(sprintf('Cluster Template ID = %d, cluster member size = %d', ...
        clustCenter, size(chosenSpikes,2)));
end

% Update handles structure for access in other functions
guidata(hObject, handles);

% Now make the GUI visible if the axes commands above haven't done it
set(hObject, 'Visible', 'on');

% Don't return until ready
uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = DetectionThresholdSelector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Returns current thresholds
varargout = {[handles.detectionThreshold]};
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

uiresume(hObject);
end

% --- Executes on slider movement.
function ThreshSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

thresh =  get(hObject, 'Value');
handles.thresh = thresh;
handles.detectionThreshold = handles.minTemplates+thresh*handles.minSDTemplates;

% Update lines on histogram plot
axes(handles.bestClusterAxes);
set(handles.threshLine, 'YData', [handles.detectionThreshold handles.detectionThreshold]);

guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function ThreshSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThreshSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

init_thresh = 0.9;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', init_thresh);
end

% --- Executes on button press in FinalizeButton.
function FinalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to FinalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the figure's handle, then destroy it

uiresume(handles.figure1);
end
