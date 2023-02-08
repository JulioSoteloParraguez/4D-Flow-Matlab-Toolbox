function varargout = GUIDE_SECTIONS_LAPLACE(varargin)
% GUIDE_SECTIONS_LAPLACE MATLAB code for GUIDE_SECTIONS_LAPLACE.fig
%      GUIDE_SECTIONS_LAPLACE, by itself, creates a new GUIDE_SECTIONS_LAPLACE or raises the existing
%      singleton*.
%
%      H = GUIDE_SECTIONS_LAPLACE returns the handle to a new GUIDE_SECTIONS_LAPLACE or the handle to
%      the existing singleton*.
%
%      GUIDE_SECTIONS_LAPLACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_SECTIONS_LAPLACE.M with the given input arguments.
%
%      GUIDE_SECTIONS_LAPLACE('Property','Value',...) creates a new GUIDE_SECTIONS_LAPLACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_SECTIONS_LAPLACE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_SECTIONS_LAPLACE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_SECTIONS_LAPLACE

% Last Modified by GUIDE v2.5 25-Jan-2023 11:22:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_SECTIONS_LAPLACE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_SECTIONS_LAPLACE_OutputFcn, ...
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


% --- Executes just before GUIDE_SECTIONS_LAPLACE is made visible.
function GUIDE_SECTIONS_LAPLACE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE_SECTIONS_LAPLACE (see VARARGIN)

% Choose default command line output for GUIDE_SECTIONS_LAPLACE
handles.output = hObject;

path(path,'iso2mesh/')

% load data
handles.elem            = varargin{1}.elem;
handles.faces           = varargin{1}.faces;
handles.nodes           = varargin{1}.nodes;
handles.Laplace         = varargin{1}.Laplace;
handles.length_vessel   = varargin{1}.length_vessel;
handles.azimuth         = varargin{1}.azimuth;
handles.elevation       = varargin{1}.elevation;
handles.NODES_SECTION   = varargin{1}.NODES_SECTION;
handles.id_while        = varargin{1}.id_while;

% load('FE Mesh/elem.mat')
% load('FE Mesh/faces.mat')
% load('FE Mesh/nodes.mat')
% load('FE Laplace/Laplace.mat')
% load('FE Length Vessel/length_vessel.mat')

% handles.elem = elem;
% handles.faces = faces;
% handles.nodes = nodes;
% handles.Laplace = Laplace;
% handles.length_vessel = length_vessel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obtenemos los datos de los editores de texto que poseen por default
set(handles.edit1,'String',handles.azimuth)
set(handles.edit2,'String',handles.elevation)
set(handles.edit3,'String',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We show the laplace 
axes(handles.axes1);
plot(0,0)
patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','red')
hold on 
plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
hold off
axis off
daspect([1,1,1])
rotate3d
view(handles.azimuth,handles.elevation)

% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);
% UIWAIT makes GUIDE_SECTIONS_LAPLACE wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_SECTIONS_LAPLACE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

if handles.id_while==0
    id_while = handles.id_while;    
    azimuth = handles.elevation;
    elevation = handles.elevation;
    elem = handles.elem;
    faces = handles.faces;
    nodes = handles.nodes;
    Laplace = handles.Laplace;
    length_vessel = handles.length_vessel;
    NODES_SECTION = handles.NODES_SECTION;
    
    setappdata(0,'id_while',id_while);
    setappdata(0,'azimuth',azimuth);
    setappdata(0,'elevation',elevation);
    setappdata(0,'elem',elem);
    setappdata(0,'faces',faces);
    setappdata(0,'nodes',nodes);
    setappdata(0,'Laplace',Laplace);
    setappdata(0,'length_vessel',length_vessel);
    setappdata(0,'NODES_SECTION',NODES_SECTION);
else
    id_while = handles.id_while; 
    azimuth = handles.elevation;
    elevation = handles.elevation;
    elem = handles.elem;
    faces = handles.faces;
    nodes = handles.nodes;
    Laplace = handles.Laplace;
    length_vessel = handles.length_vessel;
    NODES_SECTION = handles.NODES_SECTION;

    setappdata(0,'id_while',id_while);
    setappdata(0,'azimuth',azimuth);
    setappdata(0,'elevation',elevation);
    setappdata(0,'elem',elem);
    setappdata(0,'faces',faces);
    setappdata(0,'nodes',nodes);
    setappdata(0,'Laplace',Laplace);
    setappdata(0,'length_vessel',length_vessel);
    setappdata(0,'NODES_SECTION',NODES_SECTION);
    
    delete(handles.figure1);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% obtenemos los datos de los editores de texto que poseen por default
handles.azimuth = str2num(get(handles.edit1,'String'));
handles.elevation = str2num(get(handles.edit2,'String'));

set(handles.edit1,'String',handles.azimuth);
set(handles.edit2,'String',handles.elevation);

% mostramos laplace en pantalla
axes(handles.axes1);
plot(0,0)
patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','red')
hold on 
plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
hold off
daspect([1,1,1])
axis off
view(handles.azimuth,handles.elevation)

handles.output = hObject;
guidata(hObject, handles);

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sec = str2num(get(handles.edit3,'String'));


if handles.sec > 1

    % values of length that identify the section
    values_sec = linspace(0,max(handles.length_vessel),handles.sec+1);
    ID_nodes = [1:size(handles.nodes,1)]';
    
    % laplace node selected
    ID_selected = zeros(handles.sec + 1,1);
    for n = 2:2+(handles.sec-2)
        [~,I] = min(sqrt((handles.length_vessel - values_sec(n)).^2));
        ID_selected(n) = ID_nodes(I);
    end
    % [~,I] = min(handles.length_vessel(handles.length_vessel>0));
    min_value = min(handles.Laplace(handles.Laplace>0));
    I = ID_nodes(handles.Laplace==min_value);
    ID_selected(1) = ID_nodes(I(1));
    [~,I] = max(handles.length_vessel);
    ID_selected(end) = ID_nodes(I);
    
    % section of nodes.
    handles.NODES_SECTION = cell(length(ID_selected)-1,1);
    for n = 1:length(ID_selected)-1
        handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
    end

    % verifico que numero incluyo con un color y cual con otro color
    out_id = rem([1:length(ID_selected)-1], 2) == 0;
    
    % mostramos laplace en pantalla
    axes(handles.axes1);
    plot(0,0)
    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
    hold on 
    for n = 1:length(out_id)
        if out_id(n)==0
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
        else
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
        end
    end
    plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
    plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
    hold off
    daspect([1,1,1])
    axis off
    view(handles.azimuth,handles.elevation)

else

    
    ID_nodes = [1:size(handles.nodes,1)]';
    % laplace node selected
    ID_selected = zeros(handles.sec + 1,1);
%     [~,I] = min(handles.length_vessel(handles.length_vessel>0));
%     ID_selected(1) = ID_nodes(I);
    min_value = min(handles.Laplace(handles.Laplace>0));
    I = ID_nodes(handles.Laplace==min_value);
    ID_selected(1) = ID_nodes(I(1));
    
    [~,I] = max(handles.length_vessel);
    ID_selected(end) = ID_nodes(I);
    
    % section of nodes.
    handles.NODES_SECTION = cell(length(ID_selected)-1,1);
    for n = 1:length(ID_selected)-1
        handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
    end

    % verifico que numero incluyo con un color y cual con otro color
    out_id = rem([1:length(ID_selected)-1], 2) == 0;
    
    % mostramos laplace en pantalla
    axes(handles.axes1);
    plot(0,0)
    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
    hold on 
    for n = 1:length(out_id)
        if out_id(n)==0
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
        else
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
        end
    end
    plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
    plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
    hold off
    daspect([1,1,1])
    axis off
    view(handles.azimuth,handles.elevation)
end


handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sec = str2num(get(handles.edit3,'String'));

handles.points_intersec = NaN(handles.sec-1,3);
if handles.sec > 1

    for n = 1:handles.sec-1
    
        
        % mostramos laplace en pantalla
        axes(handles.axes1);
        plot(0,0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','red')
        hold on 
        plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
        plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
        plot3(handles.points_intersec(:,1),handles.points_intersec(:,2),handles.points_intersec(:,3),'.g','markersize',20);
        hold off
        daspect([1,1,1])
        axis off
        view(handles.azimuth,handles.elevation)

        uiwait(msgbox({['Please select the interesect point #',num2str(n),', and press enter'],'Press OK or Close this message first...'}));
        dcmObject = datacursormode;
        pause
        datacursormode off
        cursor = getCursorInfo(dcmObject);
        handles.points_intersec(n,:) = [cursor.Position(1),cursor.Position(2),cursor.Position(3)];
    
        % mostramos laplace en pantalla
        axes(handles.axes1);
        plot(0,0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','red')
        hold on 
        plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
        plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
        plot3(handles.points_intersec(:,1),handles.points_intersec(:,2),handles.points_intersec(:,3),'.g','markersize',20);
        hold off
        daspect([1,1,1])
        axis off
        view(handles.azimuth,handles.elevation)

    end

    % values of length that identify the section
    ID_nodes = [1:size(handles.nodes,1)]';
    
    % laplace node selected
    ID_selected = zeros(handles.sec + 1,1);
    for n = 2:2+(handles.sec-2)
        [~,I] = min(sqrt(sum((handles.nodes - repmat(handles.points_intersec(n-1,:),size(handles.nodes,1),1)).^2,2)));
        ID_selected(n) = ID_nodes(I);
    end
%     [~,I] = min(handles.length_vessel(handles.length_vessel>0));
%     ID_selected(1) = ID_nodes(I);
    min_value = min(handles.Laplace(handles.Laplace>0));
    I = ID_nodes(handles.Laplace==min_value);
    ID_selected(1) = ID_nodes(I(1));
    
    [~,I] = max(handles.length_vessel);
    ID_selected(end) = ID_nodes(I);
    
    % section of nodes.
    handles.NODES_SECTION = cell(length(ID_selected)-1,1);
    for n = 1:length(ID_selected)-1
        handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
    end

    % verifico que numero incluyo con un color y cual con otro color
    out_id = rem([1:length(ID_selected)-1], 2) == 0;
    
    % mostramos laplace en pantalla
    axes(handles.axes1);
    plot(0,0)
    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
    hold on 
    for n = 1:length(out_id)
        if out_id(n)==0
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
        else
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
        end
    end
    plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
    plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
    hold off
    daspect([1,1,1])
    axis off
    view(handles.azimuth,handles.elevation)

else

    ID_nodes = [1:size(handles.nodes,1)]';
    % laplace node selected
    ID_selected = zeros(handles.sec + 1,1);
%     [~,I] = min(handles.length_vessel(handles.length_vessel>0));
%     ID_selected(1) = ID_nodes(I);
    min_value = min(handles.Laplace(handles.Laplace>0));
    I = ID_nodes(handles.Laplace==min_value);
    ID_selected(1) = ID_nodes(I(1));
    [~,I] = max(handles.length_vessel);
    ID_selected(end) = ID_nodes(I);
    
    % section of nodes.
    handles.NODES_SECTION = cell(length(ID_selected)-1,1);
    for n = 1:length(ID_selected)-1
        handles.NODES_SECTION{n} = ID_nodes( (handles.Laplace >= handles.Laplace(ID_selected(n))) & (handles.Laplace <= handles.Laplace(ID_selected(n+1))) );
    end

    % verifico que numero incluyo con un color y cual con otro color
    out_id = rem([1:length(ID_selected)-1], 2) == 0;
    
    % mostramos laplace en pantalla
    axes(handles.axes1);
    plot(0,0)
    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','none')
    hold on 
    for n = 1:length(out_id)
        if out_id(n)==0
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.r','markersize',12);
        else
            plot3(handles.nodes(handles.NODES_SECTION{n},1),handles.nodes(handles.NODES_SECTION{n},2),handles.nodes(handles.NODES_SECTION{n},3),'.g','markersize',12);
        end
    end
    plot3(handles.nodes((handles.Laplace==0),1),handles.nodes((handles.Laplace==0),2),handles.nodes((handles.Laplace==0),3),'.y','markersize',12);
    plot3(handles.nodes((handles.Laplace==1),1),handles.nodes((handles.Laplace==1),2),handles.nodes((handles.Laplace==1),3),'.y','markersize',12);
    hold off
    daspect([1,1,1])
    axis off
    view(handles.azimuth,handles.elevation)

end



handles.output = hObject;
guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(hObject,'waitstatus'),'waiting')
    handles.id_while = 1;
    handles.output = hObject;
    guidata(hObject, handles);
    uiresume(hObject);
else
    delete(hObject);
end
