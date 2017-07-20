function varargout = ThreeA(varargin)
% THREEA MATLAB code for ThreeA.fig
%      THREEA, by itself, creates a new THREEA or raises the existing
%      singleton*.
%
%      H = THREEA returns the handle to a new THREEA or the handle to
%      the existing singleton*.
%
%      THREEA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THREEA.M with the given input arguments.
%
%      THREEA('Property','Value',...) creates a new THREEA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThreeA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThreeA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThreeA

% Last Modified by GUIDE v2.5 02-Apr-2017 22:30:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ThreeA_OpeningFcn, ...
                   'gui_OutputFcn',  @ThreeA_OutputFcn, ...
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


% --- Executes just before ThreeA is made visible.
function ThreeA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThreeA (see VARARGIN)

% Choose default command line output for ThreeA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThreeA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThreeA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
global fitnessfcn;
fitnessfcn = get(hObject,'String');

if isempty(fitnessfcn)
    errordlg('Objective function is empty !','Input Error');
    return;
else
    fitnessfcn = str2func(get(hObject,'String'));
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function starb_Callback(hObject, eventdata, handles)
% hObject    handle to starb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of starb as text
%        str2double(get(hObject,'String')) returns contents of starb as a double
global starbound;
[filename,pathname] = uigetfile({'*.xls;*.xlsx','Excel Files (*.xls, *.xlsx)';'*.mat','Mat Files'},'Select Datasets');
if strcmp(filename(length(filename) - 3:length(filename)), '.mat')
    
    temp = load(fullfile(pathname,filename));
    starbound = temp.ans;

    
    if isempty(starbound)
        errordlg('Vars bound is empty !','Input Error');
        return;        
    end
    
else
    starbound = xlsread(fullfile(pathname,filename));

    if isempty(starbound)
        errordlg('Vars bound is empty !','Input Error');
        return;        
    end
    
    
end


% --- Executes during object creation, after setting all properties.
function starb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to starb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
global nvars;
nvars = get(hObject,'String');

if isempty(nvars)
    errordlg('The number of vars is empty !','Input Error');
    return;
else
    nvars = str2double(get(hObject,'String'));
end



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
global nobjs;
nobjs = get(hObject,'String');

if isempty(nobjs)
    errordlg('The number of objs is empty !','Input Error');
    return;
else
    nobjs = str2double(get(hObject,'String'));
end

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
global bh_option;
bh_option.maxgen = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
global narc;
narc = get(hObject,'String');

if isempty(narc)
    errordlg('The archive size is empty !','Input Error');
    return;
else
    narc = str2double(get(hObject,'String'));
end

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pushbutton1 as text
%        str2double(get(hObject,'String')) returns contents of pushbutton1 as a double
global fitnessfcn;
global nvars;
global nobjs;
global starbound;
global narc;
global bh_option;

axes(handles.axes1)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')
axes(handles.axes4)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')
axes(handles.axes6)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')

if isempty(fitnessfcn)
    errordlg('Input arguments has not initialized!','Input Error');
    return;    
elseif isempty(nvars)
    errordlg('Input arguments has not initialized!','Input Error');
    return;        
elseif isempty(nobjs)
    errordlg('Input arguments has not initialized!','Input Error');
    return;        
elseif isempty(starbound)
    errordlg('Input arguments has not initialized!','Input Error');
    return;        
elseif isempty(narc)
    errordlg('Input arguments has not initialized!','Input Error');
    return;
end

if isnan(bh_option.maxgen)
    bh_option.maxgen = 600;
elseif isnan(bh_option.sizestar)
    bh_option.sizestar = 300;
elseif isnan(bh_option.els_max)
    bh_option.els_max = 0.6;
elseif isnan(bh_option.els_min)  
    bh_option.els_min = 0.1;    
end


axes(handles.axes1)

j = 0:1:narc;
i = 0:nobjs;
[x,y] = meshgrid(i,j);
plot(x,y,'k')
hold on
plot(x',y','k');

title('Pareto Front')
set(gca,'xticklabel',{''});
set(gca,'yticklabel',{''});
xlabel('Objective')
ylabel('Label number')
axis([0 nobjs 0 narc])


global h;
h = handles.axes1;
bh_option.gui = 1;
AMOBH(fitnessfcn, nvars, nobjs, starbound, narc, bh_option);








function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pushbutton2 as text
%        str2double(get(hObject,'String')) returns contents of pushbutton2 as a double
close

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global se;
global dse;
global APS;
global APF;
global nobjs;
axes(handles.axes4)
plot(1:length(se),se)
hold on 
plot(1:length(dse),dse)
xlabel('t')
ylabel('H')
title('Entropy')
legend('Entropy','Delta entropy')
axes(handles.axes6)
switch(nobjs)
    case 2
        y1 = APF(:,1);
        y2 = APF(:,2);
        plot(y1,y2,'mp', 'MarkerSize',5)
        hold on
        title('Pareto front')
        xlabel('function 1')
        ylabel('function 2')
    case 3
        y1 = APF(:,1);
        y2 = APF(:,2);
        y3 = APF(:,3);
        plot3(y1,y2,y3,'mp', 'MarkerSize',5)
        hold on
        title('Pareto front')
        xlabel('function 1')
        ylabel('function 2')
        zlabel('function 3')
    otherwise
        
end






% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global APS;
global APF;

[filename,pathname] = uiputfile({'*.xls;*.xlsx','Excel Files (*.xls, *.xlsx)';'*.mat','Mat Files'}, 'Save as','ApproximateParetoSolutionSet');
if filename==0
    return
end
if strcmp(filename(length(filename)-3:length(filename)),'.mat')
    save(fullfile(pathname,filename),'APS');
else
    xlswrite(fullfile(pathname,filename),APS); 
end
[filename,pathname] = uiputfile({'*.xls;*.xlsx','Excel Files (*.xls, *.xlsx)';'*.mat','Mat Files'}, 'Save as','ApproximateParetoFront');
if filename==0
    return
end
if strcmp(filename(length(filename)-3:length(filename)),'.mat')
    save(fullfile(pathname,filename),'APF');
else
    xlswrite(fullfile(pathname,filename),APF); 
end



% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
helpdlg('AMOBH Tool Box version 1.0, developed by Alicewithrabbit, Copyright(c) 2017, 3A Team','About')


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1




% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')
axes(handles.axes4)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')
axes(handles.axes6)
cla reset
set(gca,'xticklabel',{''});
set(gca,'xcolor','w')
set(gca,'yticklabel',{''});
set(gca,'ycolor','w')
set(gca,'zticklabel',{''});
set(gca,'zcolor','w')



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
global bh_option
bh_option.sizestar = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
global bh_option;
bh_option.els_max = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
global bh_option;
bh_option.els_min = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global bh_option;
bh_option = struct('maxgen',600,'sizestar',300,'els_max',0.6,'els_min',0.1);
