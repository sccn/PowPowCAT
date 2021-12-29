% History
% 12/29/2021 Makoto. Debugged and updated.
% 12/19/2021 Makoto. Written for Nicholas Dogris.

% Copyright (C) 2021, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
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

function varargout = PowPowCAT_viewGroupSummary(varargin)
% POWPOWCAT_VIEWGROUPSUMMARY MATLAB code for PowPowCAT_viewGroupSummary.fig
%      POWPOWCAT_VIEWGROUPSUMMARY, by itself, creates a new POWPOWCAT_VIEWGROUPSUMMARY or raises the existing
%      singleton*.
%
%      H = POWPOWCAT_VIEWGROUPSUMMARY returns the handle to a new POWPOWCAT_VIEWGROUPSUMMARY or the handle to
%      the existing singleton*.
%
%      POWPOWCAT_VIEWGROUPSUMMARY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POWPOWCAT_VIEWGROUPSUMMARY.M with the given input arguments.
%
%      POWPOWCAT_VIEWGROUPSUMMARY('Property','Value',...) creates a new POWPOWCAT_VIEWGROUPSUMMARY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PowPowCAT_viewGroupSummary_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PowPowCAT_viewGroupSummary_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PowPowCAT_viewGroupSummary

% Last Modified by GUIDE v2.5 29-Dec-2021 10:41:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PowPowCAT_viewGroupSummary_OpeningFcn, ...
                   'gui_OutputFcn',  @PowPowCAT_viewGroupSummary_OutputFcn, ...
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


% --- Executes just before PowPowCAT_viewGroupSummary is made visible.
function PowPowCAT_viewGroupSummary_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PowPowCAT_viewGroupSummary (see VARARGIN)

% Choose default command line output for PowPowCAT_viewGroupSummary
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PowPowCAT_viewGroupSummary wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PowPowCAT_viewGroupSummary_OutputFcn(hObject, eventdata, handles) 
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



% --- Executes on button press in selectPrecomputedSetFilesPushbutton.
function selectPrecomputedSetFilesPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectPrecomputedSetFilesPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain .set file.
[fileName, pathName] = uigetfile('*.set', 'MultiSelect', 'on');
if isempty(fileName)
    disp('Cancelled.')
    return
end

subjIcMatrix       = [];
comodulogram_group = [];
psd_group          = [];
dipxyz_group       = [];
totalIcIdx         = 0;
for setIdx = 1:length(fileName)

    % Load data.
    currentSubjName = fileName{setIdx};
    EEG = pop_loadset('filename', currentSubjName, 'filepath', pathName, 'loadmode', 'info');

    comodulogram_group = cat(3, comodulogram_group, EEG.etc.PowPowCAT.covMatrix);
    psd_group          = cat(2, psd_group,          EEG.etc.PowPowCAT.meanPSD);

    % Detect the dipole with larger moment in case symmetric dual diploles are fitted.
        % Set up dummy dual dipoles.
        EEG.dipfit.model(1).posxyz(2,:) =  EEG.dipfit.model(1).posxyz(1,:);
        EEG.dipfit.model(1).posxyz(2,1) =  EEG.dipfit.model(1).posxyz(2,1)*-1;
        EEG.dipfit.model(1).momxyz(2,:) =  EEG.dipfit.model(1).momxyz(1,:);
        EEG.dipfit.model(1).momxyz(2,:) =  EEG.dipfit.model(1).momxyz(2,:)*1.1;
    
    for icIdx = 1:length(EEG.dipfit.model)
        totalIcIdx = totalIcIdx+1;
        [~, largerMomIdx] = max(sum(EEG.dipfit.model(icIdx).momxyz.^2,2));
        dipxyz_group(totalIcIdx,:) = EEG.dipfit.model(icIdx).posxyz(largerMomIdx,:);
        subjIcMatrix(totalIcIdx,:) = [setIdx icIdx];
    end
end

% Flatten the comodulogram.
flattenedComodulogram = reshape(comodulogram_group, [size(comodulogram_group,1)*size(comodulogram_group,2) size(comodulogram_group,3)])';

% Precomputation for optimizing the number of clusters.
kmeansClusterIdxMatrix = zeros(size(flattenedComodulogram,1),11);
meanWithinClusterDistance = nan(11+4,11);
for clustIdx = 1:11
    disp(sprintf('Precomputing for clustering: %d/%d done', clustIdx, 11))
    [IDX, ~, SUMD] = kmeans(flattenedComodulogram, clustIdx+4, 'emptyaction', 'singleton', 'maxiter', 5000, 'replicate', 20);
    kmeansClusterIdxMatrix(:,clustIdx) = IDX;
    numIcEntries = hist(IDX, 1:clustIdx+4);
    meanWithinClusterDistance(1:clustIdx+4, clustIdx) = SUMD./numIcEntries';
end
eva1 = evalclusters(flattenedComodulogram, kmeansClusterIdxMatrix, 'CalinskiHarabasz');
eva2 = evalclusters(flattenedComodulogram, kmeansClusterIdxMatrix, 'Silhouette');
eva3 = evalclusters(flattenedComodulogram, kmeansClusterIdxMatrix, 'DaviesBouldin');

figure
set(gcf, 'NumberTitle', 'off', 'name', 'PowPowCAT_viewGroupSummary: Optimizing the nubmer of clusters', 'color', [0.93 0.96 1])
subplot(2,2,1)
boxplot(meanWithinClusterDistance)
set(gca, 'xticklabel', 5:15)
title(sprintf('Distance to cluster centroid (smaller, better)'))
xlabel('Number of clusters')
ylabel('Squared Euclidean distance (a.u.)')
subplot(2,2,2)
plot(eva1); title('Calinski-Harabasz (larger, better)');
subplot(2,2,3)
plot(eva2); title('Silhouette (larger, better)');
subplot(2,2,4)
plot(eva3); title('Davies-Bouldin (smaller, better)');

% Attach the precomputed data structures to the GUI object to bring over.
dataPackage.kmeansClusterIdxMatrix = kmeansClusterIdxMatrix;
dataPackage.comodulogram_group     = comodulogram_group;
dataPackage.dipxyz_group           = dipxyz_group;
dataPackage.subjIcMatrix           = subjIcMatrix;
dataPackage.currentFreqs           = EEG.etc.PowPowCAT.freqs;
dataPackage.meanWithinClusterDistance = meanWithinClusterDistance;
handles.selectPrecomputedSetFilesPushbutton.UserData = dataPackage;



function numClustersEdit_Callback(hObject, eventdata, handles)
% hObject    handle to numClustersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numClustersEdit as text
%        str2double(get(hObject,'String')) returns contents of numClustersEdit as a double


% --- Executes during object creation, after setting all properties.
function numClustersEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numClustersEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in visualizeTheResultsPushbutton.
function visualizeTheResultsPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to visualizeTheResultsPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Obtain the precomptued data structures.
carriedOverData = handles.selectPrecomputedSetFilesPushbutton.UserData;
kmeansClusterIdxMatrix = carriedOverData.kmeansClusterIdxMatrix;
comodulogram_group     = carriedOverData.comodulogram_group;
dipxyz_group           = carriedOverData.dipxyz_group;
subjIcMatrix           = carriedOverData.subjIcMatrix;
currentFreqs           = carriedOverData.currentFreqs;
meanWithinClusterDistance = carriedOverData.meanWithinClusterDistance;


% Obtain the user-defined optimum number.
optimumIdx = str2num(get(handles.numClustersEdit, 'String'));
if optimumIdx < 5 | optimumIdx > 15
    error('Cluster number must be between 5 and 15.')
end

% Determine the optimum clustering parameters.
cluteringIdx = kmeansClusterIdxMatrix(:,optimumIdx-4);

% Cluster the data.
clusterMeanComodulogram = zeros(size(comodulogram_group,1), size(comodulogram_group,2), optimumIdx);
clusterDipxyz           = cell(1,optimumIdx);
for cluteringIdxIdx = 1:optimumIdx
    currentClusterIdx = find(cluteringIdx==cluteringIdxIdx);
    clusterMeanComodulogram(:,:,cluteringIdxIdx) = mean(comodulogram_group(:,:,currentClusterIdx),3);
    clusterDipxyz{cluteringIdxIdx} = dipxyz_group(currentClusterIdx,:);
end

% Reorder the clusters.
numIcVec = cellfun(@(x) size(x,1), clusterDipxyz);
[sortedNumIcVec,sortingIdx] = sort(numIcVec, 'descend');
sortedClusterMeanComodulogram = clusterMeanComodulogram(:,:,sortingIdx);
sortedClusterDipxyz           = clusterDipxyz(sortingIdx);

% Obtain the within-cluster distance.
currentClsDistance = meanWithinClusterDistance(:,optimumIdx);
currentClsDistance = currentClsDistance(~isnan(currentClsDistance));
sortedMeanClusterDistance = currentClsDistance(sortingIdx);

% Show dipole plot.
addpath('C:\GYREE\Nick\powpowcatUpdate\code')

% Obtain frequency ticks.
[~,tick1Hz] = min(abs(currentFreqs-1));
[~,tick2Hz] = min(abs(currentFreqs-2));
[~,tick4Hz] = min(abs(currentFreqs-4));
[~,tick8Hz] = min(abs(currentFreqs-8));
[~,tick13Hz] = min(abs(currentFreqs-13));
[~,tick20Hz] = min(abs(currentFreqs-20));
[~,tick30Hz] = min(abs(currentFreqs-30));
[~,tick40Hz] = min(abs(currentFreqs-40));
[~,tick80Hz] = min(abs(currentFreqs-80));
freqTicIdx = [tick1Hz tick2Hz tick4Hz tick8Hz tick13Hz tick20Hz tick30Hz tick40Hz tick80Hz];
freqLabels = [1 2 4 8 13 20 30 40 80];

outOfRangeIdx = find(freqLabels<currentFreqs(1) | freqLabels>currentFreqs(end));
freqTicIdx(outOfRangeIdx) = [];
freqLabels(outOfRangeIdx) = [];

figure
set(gcf, 'color', [0.13 0.13 0.13], 'NumberTitle', 'off', 'name', 'PowPowCAT_viewGroupSummary')
tiledLayoutHandles = tiledlayout(7,optimumIdx, 'TileSpacing', 'tight', 'Padding', 'tight'); % 7 is comodulogram and -40:20:60

for clusterIdx = 1:optimumIdx

    % Plot comodulogram.
    nexttile(clusterIdx)
    imagesc(sortedClusterMeanComodulogram(:,:,clusterIdx), [-0.8 0.8])
    set(gca, 'XColor', [1 1 1], 'YColor', [1 1 1],...
        'xtick', freqTicIdx, 'xticklabel', freqLabels,...
        'ytick', freqTicIdx, 'yticklabel', freqLabels)
    colormap('jet')
    axis xy
    axis square
    title(sprintf('Cluster %d (%d ICs)', clusterIdx, sortedNumIcVec(clusterIdx)), 'color', 'white')
    xlabel('Freq (Hz)')
    ylabel('Freq (Hz)')

    % Plot dipole densities.
    xyzList = sortedClusterDipxyz{clusterIdx};
    plotDipoleDensityCustom(xyzList, 20, 1);
    dipdensityHandle = gcf;
    dipDensityChildren = get(dipdensityHandle,'children');

    dipDensityChildren(8).Parent = tiledLayoutHandles;
    dipDensityChildren(8).Layout.Tile = optimumIdx*1+clusterIdx;

    dipDensityChildren(7).Parent = tiledLayoutHandles;
    dipDensityChildren(7).Layout.Tile = optimumIdx*2+clusterIdx;

    dipDensityChildren(6).Parent = tiledLayoutHandles;
    dipDensityChildren(6).Layout.Tile = optimumIdx*3+clusterIdx;

    dipDensityChildren(5).Parent = tiledLayoutHandles;
    dipDensityChildren(5).Layout.Tile = optimumIdx*4+clusterIdx;

    dipDensityChildren(4).Parent = tiledLayoutHandles;
    dipDensityChildren(4).Layout.Tile = optimumIdx*5+clusterIdx;

    dipDensityChildren(3).Parent = tiledLayoutHandles;
    dipDensityChildren(3).Layout.Tile = optimumIdx*6+clusterIdx;

    close(dipdensityHandle)
end

% Store the output.
PowPowCAT.subjIcMatrix       = carriedOverData.subjIcMatrix;
PowPowCAT.dipxyz_group       = carriedOverData.dipxyz_group;
PowPowCAT.cluteringIdx       = cluteringIdx;
PowPowCAT.comodulogram_group = comodulogram_group;
PowPowCAT.freqs              = carriedOverData.currentFreqs;
currentClusDist = carriedOverData.meanWithinClusterDistance(:,optimumIdx-4);
PowPowCAT.clusterDistances = currentClusDist(~isnan(currentClusDist));
assignin('base', 'PowPowCAT', PowPowCAT)
disp('PowPowCAT output is stored in the workspace.')

% If requested, export a part of the results.
if ~isempty(get(handles.folderPathEdit, 'String'))
    outputPath = get(handles.folderPathEdit, 'String');
    outputTable = table(PowPowCAT.subjIcMatrix(:,1), PowPowCAT.subjIcMatrix(:,2), PowPowCAT.cluteringIdx, 'VariableNames', {'subjID' 'subjIcIdx' 'clusterIdx'});
    writetable(outputTable, [outputPath filesep 'PowPowCAT_cluteringResult.xlsx'])
    disp(['PowPowCAT output is exported to ' outputPath filesep 'PowPowCAT_cluteringResult.xlsx'])
end



function folderPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to folderPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of folderPathEdit as text
%        str2double(get(hObject,'String')) returns contents of folderPathEdit as a double


% --- Executes during object creation, after setting all properties.
function folderPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folderPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in folderPushbutton.
function folderPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to folderPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exportPath = uigetdir;
set(handles.folderPathEdit, 'String', exportPath);