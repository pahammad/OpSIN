function varargout = TwinPeaksSeparationSelector(varargin)
% TWINPEAKSSEPARATIONSELECTOR MATLAB code for TwinPeaksSeparationSelector.fig
%      TWINPEAKSSEPARATIONSELECTOR, by itself, creates a new TWINPEAKSSEPARATIONSELECTOR or raises the existing
%      singleton*.
%
%      H = TWINPEAKSSEPARATIONSELECTOR returns the handle to a new TWINPEAKSSEPARATIONSELECTOR or the handle to
%      the existing singleton*.
%
%      TWINPEAKSSEPARATIONSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWINPEAKSSEPARATIONSELECTOR.M with the given input arguments.
%
%      TWINPEAKSSEPARATIONSELECTOR('Property','Value',...) creates a new TWINPEAKSSEPARATIONSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TwinPeaksSeparationSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TwinPeaksSeparationSelector_OpeningFcn via varargin.
%
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Modified from: J.P. Rueff, Aout 2004, modified Juin 2005
% Last Modified by Vivek, 7th Oct, 2011

% Begin initialization code
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TwinPeaksSeparationSelector_OpeningFcn, ...
    'gui_OutputFcn',  @TwinPeaksSeparationSelector_OutputFcn, ...
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

% --- Executes just before TwinPeaksSeparationSelector is made visible.
function TwinPeaksSeparationSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TwinPeaksSeparationSelector (see VARARGIN)

chosenSpikes = varargin{1};
clustCenter = varargin{2};
spikeTemplate = varargin{3};
minTemplates = varargin{4};
minSDTemplates = varargin{5};
dBumpSep = varargin{6};

% Choose default output for TwinPeaksSeparationSelector
varargout = dBumpSep;

% Set up histogram plot
handles.minTemplates = minTemplates;
handles.minSDTemplates = minSDTemplates;
handles.dBumpSep = 20;
handles.firstPeakLoc = 20;
handles.secondPeakLoc = 40;
handles.lengthSpike = size(chosenSpikes,1);

% Update handles structure for access in other functions
guidata(gca, handles);

% Set up initial plots when GUI isn't yet visible
if strcmp(get(hObject,'Visible'),'off')
    maxPlotLimit = max(chosenSpikes(:));
    minPlotLimit = min(chosenSpikes(:));

    plot(chosenSpikes);
    hold on;
    plot(spikeTemplate,'r*');
    handles.firstPeakLine = line([handles.firstPeakLoc handles.firstPeakLoc], ...
        [minPlotLimit maxPlotLimit]);
    set(handles.firstPeakLine, 'Color' ,'g');
    handles.secondPeakLine = line([handles.secondPeakLoc handles.secondPeakLoc], ...
        [minPlotLimit maxPlotLimit]);
    set(handles.secondPeakLine, 'Color' ,'r');
    
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
function varargout = TwinPeaksSeparationSelector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Returns current thresholds
varargout = {[handles.dBumpSep]};
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
function firstPeakSlider_Callback(hObject, eventdata, handles)
% hObject    handle to firstPeakSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

firstPeakLoc =  get(hObject, 'Value');
handles.firstPeakLoc = firstPeakLoc;
if handles.firstPeakLoc >= handles.secondPeakLoc
    handles.secondPeakLoc = handles.firstPeakLoc+1;
    set(handles.secondPeakSlider, 'Value', handles.secondPeakLoc);
end
handles.dBumpSep = round(handles.secondPeakLoc-handles.firstPeakLoc);

% Update lines on spike shape plot
axes(handles.bestClusterAxes);
set(handles.firstPeakLine, 'XData', [handles.firstPeakLoc handles.firstPeakLoc]);
set(handles.secondPeakLine, 'XData', [handles.secondPeakLoc handles.secondPeakLoc]);

guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function firstPeakSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstPeakSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

init_firstPeakLoc = 20;
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', init_firstPeakLoc);
end

% --- Executes on button press in FinalizeButton.
function FinalizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to FinalizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get the figure's handle, then destroy it

uiresume(handles.figure1);
end


% --- Executes on slider movement.
function secondPeakSlider_Callback(hObject, eventdata, handles)
% hObject    handle to secondPeakSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

secondPeakLoc =  get(hObject, 'Value');
handles.secondPeakLoc = secondPeakLoc;
if handles.secondPeakLoc > handles.lengthSpike
    handles.secondPeakLoc = handles.lengthSpike;
    set(handles.secondPeakSlider, 'Value', handles.lengthSpike);
end

if handles.firstPeakLoc >= handles.secondPeakLoc
    handles.firstPeakLoc = handles.secondPeakLoc-1;
    set(handles.firstPeakSlider, 'Value', handles.firstPeakLoc);
end
handles.dBumpSep = round(handles.secondPeakLoc-handles.firstPeakLoc);

% Update lines on spike shape plot
axes(handles.bestClusterAxes);
set(handles.firstPeakLine, 'XData', [handles.firstPeakLoc handles.firstPeakLoc]);
set(handles.secondPeakLine, 'XData', [handles.secondPeakLoc handles.secondPeakLoc]);

guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function secondPeakSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to secondPeakSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

init_secondPeakLoc = 40;

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Value', init_secondPeakLoc);
end
