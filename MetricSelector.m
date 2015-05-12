function varargout = MetricSelector(varargin)
% METRICSELECTOR MATLAB code for MetricSelector.fig
%      METRICSELECTOR, by itself, creates a new METRICSELECTOR or raises the existing
%      singleton*.
%
%      H = METRICSELECTOR returns the handle to a new METRICSELECTOR or the handle to
%      the existing singleton*.
%
%      METRICSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in METRICSELECTOR.M with the given input arguments.
%
%      METRICSELECTOR('Property','Value',...) creates a new METRICSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MetricSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MetricSelector_OpeningFcn via varargin.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Modified from: J.P. Rueff, Aout 2004, modified Juin 2005
% Last Modified by Vivek, 4th Sept, 2011

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MetricSelector_OpeningFcn, ...
    'gui_OutputFcn',  @MetricSelector_OutputFcn, ...
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

% --- Executes just before MetricSelector is made visible.
function MetricSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MetricSelector (see VARARGIN)

minDist = varargin{1};
corrToTemplate = varargin{2};
init_alpha = varargin{3};

% Choose default output for MetricSelector
varargout = {init_alpha};

% Set up plot
handles.minDist = minDist;
handles.corrToTemplate = corrToTemplate;
handles.alpha = init_alpha; % Initial value
handles.numHistBins = 40;
handles.spikeProbs = init_alpha*minDist + (1-init_alpha)*corrToTemplate;

% Update handles structure for access in other functions
guidata(gca, handles);

% Set up initial plots when GUI isn't yet visible
if strcmp(get(hObject,'Visible'),'off')
    % Plot histogram
    hist(handles.bestSeparatorAxes, handles.spikeProbs, handles.numHistBins);
    title(handles.bestSeparatorAxes, ...
        ['Spike Probability Histogram, Alpha = ' num2str(handles.alpha)]);
end

% Update handles structure for access in other functions
guidata(hObject, handles);

% Now make the GUI visible if the axes commands above haven't done it
set(hObject, 'Visible', 'on');

% Don't return until ready
uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = MetricSelector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Returns current thresholds
varargout = {[handles.alpha]};
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
function AlphaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to AlphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

alpha =  get(hObject, 'Value');
handles.alpha = alpha;

handles.spikeProbs = alpha*handles.minDist + (1-alpha)*handles.corrToTemplate;

% Update lines on histogram plot
axes(handles.bestSeparatorAxes);
hist(handles.bestSeparatorAxes, handles.spikeProbs, handles.numHistBins);
title(handles.bestSeparatorAxes, ...
    ['Spike Probability Histogram, Alpha = ' num2str(handles.alpha)]);

guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function AlphaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

init_alpha = 0.5;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', init_alpha);
end

% --- Executes on button press in FinalizeButton.
function FinalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to FinalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the figure's handle, then destroy it

uiresume(handles.figure1);
end
