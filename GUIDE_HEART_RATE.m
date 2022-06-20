function varargout = GUIDE_HEART_RATE(varargin)
% GUIDE_HEART_RATE MATLAB code for GUIDE_HEART_RATE.fig
% This code is develop by Dr. Julio Sotelo, that together with 
% Dr. Sergio Uribe and Dr. Daniel Hurtado, have been working in Cardiac MR 
% and particularly in the quantification of 4D flow data for about 7 years 
% in the Biomedical Imaging center at Universidad Catolica of Chile
% (www.mri.cl).
% We have developed a methodology for the non-invasive quantification of 
% hemodynamics from 4D flow data sets based on Finite Elements methods. 
% This technique is unique and is possible obtain several hemodynamic 
% parameters in 3D as WSS, OSI, vorticity, helicity density, relative helicity
% density, viscouss dissipation, energy loss and kinetic energy. 
%   Author:      Dr. Julio Sotelo
%   Time-stamp:  2018-04-08 v1.0
%   E-mail:      jasotelo@uc.cl
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_HEART_RATE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_HEART_RATE_OutputFcn, ...
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
function GUIDE_HEART_RATE_OpeningFcn(hObject, eventdata, handles, varargin)
handles.heart_rate = 60;
handles.output = hObject;
guidata(hObject, handles);
uiwait(handles.figure1);
function varargout = GUIDE_HEART_RATE_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.heart_rate;
delete(handles.figure1);
function edit1_Callback(hObject, eventdata, handles)
handles.heart_rate = str2double(get(hObject,'String'));
handles.output = hObject;
guidata(hObject, handles);
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbutton1_Callback(hObject, eventdata, handles)
close(handles.figure1);
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');

    set(handles.edit1,'FontUnits','Normalized','FontSize',0.39)
    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.26)
    set(handles.text2,'FontUnits','Normalized','FontSize',0.5)

handles.output = hObject;
guidata(hObject, handles); 
