function varargout = GUIDE_CONTRAST(varargin)
% GUIDE_CONTRAST MATLAB code for GUIDE_CONTRAST.fig
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
                   'gui_OpeningFcn', @GUIDE_CONTRAST_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_CONTRAST_OutputFcn, ...
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

function GUIDE_CONTRAST_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

    handles.IPCMRA = varargin{1};
    handles.IPCMRA_new = zeros(size(handles.IPCMRA));
    handles.IPCMRA_reset = varargin{2};
    handles.min_value = min(handles.IPCMRA(:));
    handles.max_value = max(handles.IPCMRA(:));
    slider_step(1) = 1/255;
    slider_step(2) = 10/255;
    set(handles.slider1,'Value', handles.min_value/handles.max_value,'sliderstep',slider_step,'max',1,'min',0)
    slider_step(1) = 1/255;
    slider_step(2) = 10/255;
    set(handles.slider2,'Value', handles.max_value/handles.max_value,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.text5,'String',num2str(handles.min_value))
    set(handles.text6,'String',num2str(handles.max_value))
    [counts,binLocations] = imhist(uint8(handles.IPCMRA(:)),100);
    axes(handles.axes1);
    histogram(uint8(handles.IPCMRA(:)),255)
    axis([0 255 0 mean(counts)])
    ax = gca;
    ax.FontWeight = 'bold';
    ax.YColor = 'white';
    ax.XColor = 'white';
    ax.Title.String = 'Histogram';
    ax.Title.Color = 'white';
    ax.Title.FontWeight = 'bold';
    hold on
    plot([handles.min_value handles.min_value],[0 max(counts)],'-r','Linewidth',2)
    plot([handles.max_value handles.max_value],[0 max(counts)],'-m','Linewidth',2)
    hold off
    handles.id = 3;
    set(handles.uipanel1,'Visible','on')
    
guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = GUIDE_CONTRAST_OutputFcn(hObject, eventdata, handles)
if handles.id == 1

    varargout{1} = handles.output;
    finish_id = 0;
    setappdata(0,'closed_loop',finish_id);
    IPCMRA_new = handles.IPCMRA_reset;
    setappdata(0,'OUT',IPCMRA_new);
    
elseif handles.id == 2

    varargout{1} = handles.output;
    finish_id = 0;
    setappdata(0,'closed_loop',finish_id);
    IPCMRA_new = handles.IPCMRA_new;
    setappdata(0,'OUT',IPCMRA_new);
    
elseif handles.id == 3

    varargout{1} = handles.output;
    IPCMRA_new = handles.IPCMRA;
    setappdata(0,'OUT',IPCMRA_new);
    finish_id = 1;
    setappdata(0,'closed_loop',finish_id);
    delete(handles.figure1);
end

function slider1_Callback(hObject, eventdata, handles)
    pp=1/255;
    slider_step(1) = pp;
    slider_step(2) = 10/255;
    set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = 255;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.min_value = handles.slider_value;
    set(handles.text5,'String',num2str(handles.slider_value));
   
    [counts,binLocations] = imhist(uint8(handles.IPCMRA(:)),100);
    axes(handles.axes1);
    histogram(uint8(handles.IPCMRA(:)),255)
    axis([0 255 0 mean(counts)])
    ax = gca;
    ax.FontWeight = 'bold';
    ax.YColor = 'white';
    ax.XColor = 'white';
    ax.Title.String = 'Histogram';
    ax.Title.Color = 'white';
    ax.Title.FontWeight = 'bold';
    hold on
    plot([handles.min_value handles.min_value],[0 max(counts)],'-r','Linewidth',2)
    plot([handles.max_value handles.max_value],[0 max(counts)],'-m','Linewidth',2)
    hold off
handles.output = hObject;
guidata(hObject, handles);

function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider2_Callback(hObject, eventdata, handles)
    pp=1/255;
    slider_step(1) = pp;
    slider_step(2) = 10/255;
    set(handles.slider2,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = 255;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.max_value = handles.slider_value;
    set(handles.text6,'String',num2str(handles.slider_value));

    [counts,binLocations] = imhist(uint8(handles.IPCMRA(:)),100);
    axes(handles.axes1);
    histogram(uint8(handles.IPCMRA(:)),255)
    axis([0 255 0 mean(counts)])
    ax = gca;
    ax.FontWeight = 'bold';
    ax.YColor = 'white';
    ax.XColor = 'white';
    ax.Title.String = 'Histogram';
    ax.Title.Color = 'white';
    ax.Title.FontWeight = 'bold';
    hold on
    plot([handles.min_value handles.min_value],[0 max(counts)],'-r','Linewidth',2)
    plot([handles.max_value handles.max_value],[0 max(counts)],'-m','Linewidth',2)
    hold off
handles.output = hObject;
guidata(hObject, handles);

function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function pushbutton1_Callback(hObject, eventdata, handles)
handles.id = 1;
handles.output = hObject;
guidata(hObject, handles);
close(handles.figure1);

function pushbutton2_Callback(hObject, eventdata, handles)
handles.id = 2;

for n=1:size(handles.IPCMRA,3)
    handles.IPCMRA_new(:,:,n) = imadjust(uint8(handles.IPCMRA(:,:,n)),[handles.min_value/255,handles.max_value/255],[0 1]);
end

handles.output = hObject;
guidata(hObject, handles);
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

set(handles.text2,'FontUnits','Normalized','FontSize',0.55)
set(handles.text3,'FontUnits','Normalized','FontSize',0.55)
set(handles.text4,'FontUnits','Normalized','FontSize',0.55)
set(handles.text5,'FontUnits','Normalized','FontSize',0.66)
set(handles.text6,'FontUnits','Normalized','FontSize',0.66)
set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.37)
set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.37)

handles.output = hObject;
guidata(hObject, handles);
