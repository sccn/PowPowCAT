% History
% 05/13/2021 Makoto. Waitbar bug fixed.
% 12/25/2020 Makoto. Batch version created for Pal Gunnar Larsson.

% Copyright (C) 2020, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function varargout = PowPowCAT_batch(varargin)
% POWPOWCAT_BATCH MATLAB code for PowPowCAT_batch.fig
%      POWPOWCAT_BATCH, by itself, creates a new POWPOWCAT_BATCH or raises the existing
%      singleton*.
%
%      H = POWPOWCAT_BATCH returns the handle to a new POWPOWCAT_BATCH or the handle to
%      the existing singleton*.
%
%      POWPOWCAT_BATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWPOWCAT_BATCH.M with the given input arguments.
%
%      POWPOWCAT_BATCH('Property','Value',...) creates a new POWPOWCAT_BATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PowPowCAT_batch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PowPowCAT_batch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PowPowCAT_batch

% Last Modified by GUIDE v2.5 25-Dec-2020 12:40:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PowPowCAT_batch_OpeningFcn, ...
                   'gui_OutputFcn',  @PowPowCAT_batch_OutputFcn, ...
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


% --- Executes just before PowPowCAT_batch is made visible.
function PowPowCAT_batch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PowPowCAT_batch (see VARARGIN)

% Choose default command line output for PowPowCAT_batch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PowPowCAT_batch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PowPowCAT_batch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in inputDataPopupmenu.
function inputDataPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to inputDataPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns inputDataPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from inputDataPopupmenu


% --- Executes during object creation, after setting all properties.
function inputDataPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to inputDataPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in startPushbutton.
function startPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain .set file.
[fileName, pathName] = uigetfile('*.set', 'MultiSelect', 'on');
if isempty(fileName)
    disp('Cancelled.')
    return
end

upperFreqLimit = str2num(get(handles.upperFreqLimitEdit, 'String'));
inputDataType  = get(handles.inputDataPopupmenu, 'Value');
methodType     = get(handles.methodPopupmenu, 'Value');
numIterations  = str2num(get(handles.upperFreqLimitEdit, 'String'));

waitbarHandle = waitbar(0,'Please wait...');
% Loop for the subjects.
for subjIdx = 1:length(fileName)
    
    % Plot the wait bar.
    waitbar(subjIdx/length(fileName), waitbarHandle, sprintf('%d/%d subjects done.', subjIdx, length(fileName)));
    
    % Load data.
    currentSubjName = fileName{subjIdx};
    EEG = pop_loadset('filename', currentSubjName, 'filepath', pathName);
    
    % Run PowPowCAT.
    EEG = calc_PowPowCAT(EEG, upperFreqLimit, inputDataType, methodType, numIterations);

    % Save data.
    pop_saveset(EEG, 'filename', currentSubjName, 'filepath', pathName)
end
close


% --- Executes on selection change in methodPopupmenu.
function methodPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to methodPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methodPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methodPopupmenu


% --- Executes during object creation, after setting all properties.
function methodPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upperFreqLimitEdit_Callback(hObject, eventdata, handles)
% hObject    handle to upperFreqLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperFreqLimitEdit as text
%        str2double(get(hObject,'String')) returns contents of upperFreqLimitEdit as a double


% --- Executes during object creation, after setting all properties.
function upperFreqLimitEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperFreqLimitEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function surroStatsIterEdit_Callback(hObject, eventdata, handles)
% hObject    handle to surroStatsIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of surroStatsIterEdit as text
%        str2double(get(hObject,'String')) returns contents of surroStatsIterEdit as a double


% --- Executes during object creation, after setting all properties.
function surroStatsIterEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to surroStatsIterEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
