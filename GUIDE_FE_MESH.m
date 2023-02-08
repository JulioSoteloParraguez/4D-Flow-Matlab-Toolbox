function varargout = GUIDE_FE_MESH(varargin)
% GUIDE_FE_MESH MATLAB code for GUIDE_FE_MESH.fig
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
                   'gui_OpeningFcn', @GUIDE_FE_MESH_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_FE_MESH_OutputFcn, ...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GUIDE_FE_MESH_OpeningFcn(hObject, eventdata, handles, varargin)

    handles.output = hObject;
    
    handles.cont_show = 1;
    
    path(path,'iso2mesh/')

    %%%%%%%%%%%%%%%%%%%%%%% read input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.SEG           = varargin{1}.SEG;
    handles.voxel_MR      = varargin{1}.voxel_MR;
    handles.L             = varargin{1}.L;
    handles.Lrgb          = varargin{1}.Lrgb;
    handles.NUM           = varargin{1}.NUM;
    handles.MR_PCA_FH     = varargin{1}.MR_PCA_FH; % Load velocity
    handles.MR_PCA_AP     = varargin{1}.MR_PCA_AP; % Load velocity
    handles.MR_PCA_RL     = varargin{1}.MR_PCA_RL; % Load velocity
    handles.xd            = varargin{1}.xd;
    handles.yd            = varargin{1}.yd;
    handles.zd            = varargin{1}.zd;
    handles.a             = varargin{1}.a;
    handles.b             = varargin{1}.b;
    handles.c             = varargin{1}.c;
    handles.d             = varargin{1}.d;
    handles.IPCMRA        = varargin{1}.IPCMRA;
    
    %%%%%%%%%%%%%%%%%%%%%%% read input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % increase the size of the image 
    [handles.a_n,handles.b_n,handles.c_n,handles.d_n] = size(handles.MR_PCA_FH);
    
    SEG_n                                   = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4);
    L_n                                     = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4);
    Lrgb_n                                  = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4,3);
    IPCMRA_n                                = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4);
    MR_PCA_FH_n                             = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4,handles.d_n);
    MR_PCA_AP_n                             = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4,handles.d_n);
    MR_PCA_RL_n                             = zeros(handles.a_n+4,handles.b_n+4,handles.c_n+4,handles.d_n);
    SEG_n(3:end-2,3:end-2,3:end-2)          = handles.SEG;
    L_n(3:end-2,3:end-2,3:end-2)            = handles.L;
    Lrgb_n(3:end-2,3:end-2,3:end-2,:)       = handles.Lrgb;
    IPCMRA_n(3:end-2,3:end-2,3:end-2)       = handles.IPCMRA;
    MR_PCA_FH_n(3:end-2,3:end-2,3:end-2,:)  = handles.MR_PCA_FH; % velocity information FH
    MR_PCA_AP_n(3:end-2,3:end-2,3:end-2,:)  = handles.MR_PCA_AP; % velocity information AP
    MR_PCA_RL_n(3:end-2,3:end-2,3:end-2,:)  = handles.MR_PCA_RL; % velocity information RL
    handles.SEG                             = SEG_n;
    handles.L                               = L_n;
    handles.Lrgb                            = Lrgb_n;
    handles.IPCMRA                          = IPCMRA_n;
    handles.MR_PCA_FH                       = MR_PCA_FH_n; % velocity information FH
    handles.MR_PCA_AP                       = MR_PCA_AP_n; % velocity information AP
    handles.MR_PCA_RL                       = MR_PCA_RL_n; % velocity information RL
    handles.MR_PCA_FH_ori                   = MR_PCA_FH_n; % velocity information FH original information of velocity
    handles.MR_PCA_AP_ori                   = MR_PCA_AP_n; % velocity information AP original information of velocity
    handles.MR_PCA_RL_ori                   = MR_PCA_RL_n; % velocity information RL original information of velocity
    
    [X,Y,Z] = meshgrid(0:size(handles.SEG,1)-1,0:size(handles.SEG,2)-1,0:size(handles.SEG,3)-1);
    
    handles.xd                              = X*handles.voxel_MR(1);
    handles.yd                              = Y*handles.voxel_MR(2);
    handles.zd                              = Z*handles.voxel_MR(3);
    handles.xd                              = permute(handles.xd,[2 1 3]);
    handles.yd                              = permute(handles.yd,[2 1 3]);
    handles.zd                              = permute(handles.zd,[2 1 3]);
    handles.S                               = handles.SEG(:);
    handles.R                               = squeeze(handles.Lrgb(:,:,:,1));
    handles.G                               = squeeze(handles.Lrgb(:,:,:,2));
    handles.B                               = squeeze(handles.Lrgb(:,:,:,3));
    handles.R                               = handles.R(:);
    handles.G                               = handles.G(:);
    handles.B                               = handles.B(:);
    
    handles.id_kernel                       = 0;
    handles.int_method                      = 'cubic';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cla(handles.axes1)
    axes(handles.axes1);
    data = smooth3(handles.SEG,'box',3);
    fv = isosurface(handles.xd,handles.yd,handles.zd,data,.5);
    handles.fv = fv;
    patch(fv,'FaceColor',[mean(handles.R(handles.S==1)), mean(handles.G(handles.S==1)), mean(handles.B(handles.S==1))],'EdgeColor','k');
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])
    colorbar('off')
    IMG = handles.SEG;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [handles.node_bin_aorta,handles.elem_bin_aorta] = binsurface(IMG,4);
    handles.node_bin_aorta(:,1) = handles.node_bin_aorta(:,1)*handles.voxel_MR(1)-0.5*handles.voxel_MR(1);
    handles.node_bin_aorta(:,2) = handles.node_bin_aorta(:,2)*handles.voxel_MR(2)-0.5*handles.voxel_MR(2);
    handles.node_bin_aorta(:,3) = handles.node_bin_aorta(:,3)*handles.voxel_MR(3)-0.5*handles.voxel_MR(3);
    
    
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton2,'Value',1);
    set(handles.radiobutton3,'Value',0);
    
    set(handles.text01,'visible','on','String',['Voxel Size: ',sprintf('%.2f',handles.voxel_MR(1)),' x ',sprintf('%.2f',handles.voxel_MR(2)),' x ',sprintf('%.2f',handles.voxel_MR(3))]);
    
    handles.gridsize    = mean(handles.voxel_MR)*0.5;
    handles.closesize   = 0;
    handles.elemsize    = mean(handles.voxel_MR)*0.5;
    
    set(handles.edit1,'string',num2str(handles.gridsize,3))
    set(handles.edit2,'string',num2str(handles.closesize,3))
    set(handles.edit3,'string',num2str(handles.elemsize,3))
    
    handles.iter        = 100;
    handles.useralpha   = 0.1;
    handles.userbeta    = 0.1;
    
    set(handles.edit4,'string',num2str(handles.iter,3))
    set(handles.edit5,'string',num2str(handles.useralpha,3))
    set(handles.edit6,'string',num2str(handles.userbeta,3))
    

%     if ismac
%         
%         % Code to run on Mac platform
%         
%         handles.vol_factor      = 0.86;
%         
%         set(handles.edit7,'string',num2str(handles.vol_factor,3))
%         set(handles.text11,'Visible','off')
%         set(handles.edit8,'Visible','off')
%         
%         handles.suggestedvolume = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3))*0.2; % 20
%         set(handles.text12,'string',['Suggested Volume: ',num2str(handles.suggestedvolume,3),' [cubic mm]'])
%        
%     elseif isunix
%         % Code to run on Linux platform
%     
%         handles.vol_factor      = 0.86;
%         
%         set(handles.edit7,'string',num2str(handles.vol_factor,3))
%         set(handles.text11,'Visible','off')
%         set(handles.edit8,'Visible','off')
%         
%         handles.suggestedvolume = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3))*0.2; % 20
%         set(handles.text12,'string',['Suggested Volume: ',num2str(handles.suggestedvolume,3),' [cubic mm]'])
%         
%     elseif ispc
%         % Code to run on widnows platform
%         

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial parameters for the mesh
    handles.vol_factor      = 2;
    handles.keepratio       = 0.3;
    set(handles.edit7,'string',num2str(handles.vol_factor,3),'visible','on')
    set(handles.edit8,'string',num2str(handles.keepratio,3),'visible','on')

    handles.suggestedvolume = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3))*0.2; % 20
    set(handles.text12,'style', 'text','string',['Suggested Volume, less than ',num2str(handles.suggestedvolume,3),' [mm',char(179),']'])

    set(handles.text10,'Visible','on')
    set(handles.text11,'Visible','on')

%     end
        
    
    handles.idmesh = 0;
    handles.faces = [];
    handles.nodes = [];
    handles.elem = [];
    handles.veset = [];
    handles.peak_flow = [];
    handles.cmap = [];
    handles.mags = [];
    handles.mags_vol = [];
    handles.Lrgb_vel = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot surface
    fv = handles.fv;
    
    handles.nodes1 = fv.vertices;
    handles.faces1 = fv.faces;
    
    handles.nodes = fv.vertices;
    handles.faces = fv.faces;
    
    % surface area
    d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
    d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
    d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
    s = (d1+d2+d3)/2;
    
    handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    
    
    
    set(handles.text01, 'FontUnits', 'points');
    FontS = get(handles.text01, 'FontSize');
    set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);
    
    axes(handles.axes2);
    h = histfit(handles.surface_area,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Area [mm^{2}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I))];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str(mean([d1;d2;d3]),3),' [mm]']);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.radiobutton5,'Visible','on','Value',1)
    set(handles.radiobutton6,'Visible','off','Value',0)
    set(handles.radiobutton7,'Visible','off','Value',0)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.uipanel1,'Visible','on')
    set(handles.uipanel2,'Visible','on')
    set(handles.uipanel3,'Visible','on')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.vector_size = 5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.v_size_image = handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3);
%     if handles.v_size_image<0.15
%         
%         % Code to run on Windows platform
%         handles.vol_factor      = 2;
%         handles.keepratio       = 0.3;
%         set(handles.edit7,'string',num2str(handles.vol_factor,3),'visible','on')
%         set(handles.edit8,'string',num2str(handles.keepratio,3),'visible','on')
%         set(handles.text10,'Visible','on')
%         set(handles.text11,'Visible','on')
%     end
    
    handles.id_mod = 0;
    handles.kernel_size = 3;

guidata(hObject, handles);
uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = GUIDE_FE_MESH_OutputFcn(hObject, eventdata, handles)

    varargout{1} = handles.output;
    faces = handles.faces;
    setappdata(0,'faces',faces);
    nodes = [handles.nodes(:,1)-(2*handles.voxel_MR(1)),handles.nodes(:,2)-(2*handles.voxel_MR(2)),handles.nodes(:,3)-(2*handles.voxel_MR(3))];
    setappdata(0,'nodes',nodes);
    elem = handles.elem;
    setappdata(0,'elem',elem);
    veset = handles.veset;
    setappdata(0,'veset',veset);
    peak_flow = handles.peak_flow;
    setappdata(0,'peak_flow',peak_flow);
    cmap = handles.cmap;
    setappdata(0,'cmap',cmap);
    mags = handles.mags;
    setappdata(0,'mags',mags);
    mags_vol = handles.mags_vol;
    setappdata(0,'mags_vol',mags_vol);
    Lrgb_vel = handles.Lrgb_vel;
    setappdata(0,'Lrgb_vel',Lrgb_vel);
    MR_PCA_AP = handles.MR_PCA_AP;
    setappdata(0,'MR_PCA_AP',MR_PCA_AP);
    MR_PCA_FH = handles.MR_PCA_FH;
    setappdata(0,'MR_PCA_FH',MR_PCA_FH);
    MR_PCA_RL = handles.MR_PCA_RL;
    setappdata(0,'MR_PCA_RL',MR_PCA_RL);
    MR_PCA_AP_ori = handles.MR_PCA_AP_ori; % change the direction of the velocity original
    setappdata(0,'MR_PCA_AP_ori',MR_PCA_AP_ori);
    MR_PCA_FH_ori = handles.MR_PCA_FH_ori; % change the direction of the  velocity original
    setappdata(0,'MR_PCA_FH_ori',MR_PCA_FH_ori);
    MR_PCA_RL_ori = handles.MR_PCA_RL_ori; % change the direction of the velocity original
    setappdata(0,'MR_PCA_RL_ori',MR_PCA_RL_ori);
    
    delete(handles.figure1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit1_Callback(hObject, eventdata, handles)
    
    handles.gridsize = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit1_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit2_Callback(hObject, eventdata, handles)

    handles.closesize = str2double(get(hObject,'String'));
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit2_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit3_Callback(hObject, eventdata, handles)

    handles.elemsize = str2double(get(hObject,'String'));
        
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit3_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit4_Callback(hObject, eventdata, handles)

    handles.iter = str2double(get(hObject,'String'));
        
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit4_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit5_Callback(hObject, eventdata, handles)
    
    handles.useralpha = str2double(get(hObject,'String'));
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit5_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit6_Callback(hObject, eventdata, handles)

    handles.userbeta = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit6_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit7_Callback(hObject, eventdata, handles)

    handles.vol_factor = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit7_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit8_Callback(hObject, eventdata, handles)

    handles.keepratio = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit8_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushbutton1_Callback(hObject, eventdata, handles)
   
    opt.gridsize  = handles.gridsize;
    opt.closesize = handles.closesize;
    opt.elemsize  = handles.elemsize;
    
    f = waitbar(0,'Please wait: Remeshing surface ...');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [newno,newfc] = remeshsurf(handles.nodes1,handles.faces1,opt);

    handles.faces2 = newfc(:,1:3);
    handles.nodes2 = newno(:,1:3);
    
    handles.faces = newfc(:,1:3);
    handles.nodes = newno(:,1:3);
    
    close(f)
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])
    colorbar('off')
    handles.idmesh = 1;
    set(handles.radiobutton1,'Value',1);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface area
    d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
    d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
    d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
    s = (d1+d2+d3)/2;
    
    handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    
    
    set(handles.text01, 'FontUnits', 'points');
    FontS = get(handles.text01, 'FontSize');
    set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);
    
    axes(handles.axes2);
    h = histfit(handles.surface_area,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Area [mm^{2}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I))];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str(mean([d1;d2;d3]),3),' [mm]']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.radiobutton1,'Visible','on','Value',1);
    set(handles.radiobutton5,'Visible','on','Value',1)
    set(handles.radiobutton6,'Visible','off','Value',0)
    set(handles.radiobutton7,'Visible','off','Value',0)
    set(handles.pushbutton4,'Visible','off');
    set(handles.radiobutton4,'Visible','off');
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushbutton2_Callback(hObject, eventdata, handles)

    f = waitbar(0,'Please wait: Smoothing surface ...');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% smooth-surface
    conn = meshconn(handles.faces(:,1:3),size(handles.nodes,1));
    handles.nodes = smoothsurf(handles.nodes2,[],conn,handles.iter,handles.useralpha,'laplacianhc', handles.userbeta);
    
    handles.faces3 = handles.faces(:,1:3);
    handles.nodes3 = handles.nodes(:,1:3);
    
    close(f)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])
    colorbar('off')
    handles.idmesh = 1;
    set(handles.radiobutton1,'Value',1);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface area
    d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
    d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
    d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
    s = (d1+d2+d3)/2;
    
    handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    
    
    set(handles.text01, 'FontUnits', 'points');
    FontS = get(handles.text01, 'FontSize');
    set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);
    
    axes(handles.axes2);
    h = histfit(handles.surface_area,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Area [mm^{2}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I))];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str(mean([d1;d2;d3]),3),' [mm]']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.radiobutton1,'Visible','on','Value',1);
    set(handles.radiobutton5,'Visible','on','Value',1)
    set(handles.radiobutton6,'Visible','off','Value',0)
    set(handles.radiobutton7,'Visible','off','Value',0)
    set(handles.pushbutton4,'Visible','off');
    set(handles.radiobutton4,'Visible','off');
    
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushbutton3_Callback(hObject, eventdata, handles)

    %%% surface-tetrahedral mesh
    
    f = waitbar(0,'Please wait: Generating tetrahedral mesh ...');
        
%     if handles.v_size_image<0.15
%         
%         % Code to run on Windows platform
%             
%         [node,elem,face]=s2m(handles.nodes3(:,1:3),handles.faces3(:,1:3),handles.keepratio,handles.suggestedvolume*handles.vol_factor);
% 
%         
%     else
%         if ismac
%         % Code to run on Mac platform
% 
%             arista = ((handles.suggestedvolume*12)/sqrt(2))^(1/3);
%             SuggestedArea = ((sqrt(3)/4)*mean(arista)^2);
% 
%             opt2.radbound = SuggestedArea*3;
%             opt2.angbound = 0;
%             opt2.distbound = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3));
%             opt2.reratio = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3)); 
% 
%             [node,elem,~]=s2m(handles.nodes3(:,1:3),handles.faces3(:,1:3),opt2,handles.suggestedvolume*handles.vol_factor,'cgalpoly');
%             elem = meshreorient(node(:,1:3),elem(:,1:4));
%             face = volface(elem(:,1:4));
% 
%         elseif isunix
%         % Code to run on Linux platform
% 
%             arista = ((handles.suggestedvolume*12)/sqrt(2))^(1/3);
%             SuggestedArea = ((sqrt(3)/4)*mean(arista)^2);
% 
%             opt2.radbound = SuggestedArea*3;
%             opt2.angbound = 0;
%             opt2.distbound = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3));
%             opt2.reratio = (handles.voxel_MR(1)*handles.voxel_MR(2)*handles.voxel_MR(3)); 
% 
%             [node,elem,~]=s2m(handles.nodes3(:,1:3),handles.faces3(:,1:3),opt2,handles.suggestedvolume*handles.vol_factor,'cgalpoly');
%             elem = meshreorient(node(:,1:3),elem(:,1:4));
%             face = volface(elem(:,1:4));
% 
%         elseif ispc
%         Code to run on Windows platform

         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface to tetrahedral mesh
    [node,elem,face]=s2m(handles.nodes3(:,1:3),handles.faces3(:,1:3),handles.keepratio,handles.suggestedvolume*handles.vol_factor);

%         end
%     end
    
    
    handles.faces = face(:,1:3);
    handles.nodes = node(:,1:3);
    handles.elem = elem(:,1:4);

    close(f)
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])
    colorbar('off')
    handles.idmesh = 1;
    set(handles.radiobutton1,'Value',1);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % surface area
    d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
    d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
    d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
    s = (d1+d2+d3)/2;
    
    handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    
    
    set(handles.text01, 'FontUnits', 'points');
    FontS = get(handles.text01, 'FontSize');
    set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);
    
    set(handles.axes2,'visible','on')
    axes(handles.axes2);
    h = histfit(handles.surface_area,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Area [mm^{2}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I))];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str(sqrt(xdata(I)*4/sqrt(3)),3),' [mm]'],'visible','on');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.pushbutton4,'Visible','on');
    set(handles.radiobutton1,'Visible','on');
    set(handles.radiobutton5,'Visible','on','Value',1);
    set(handles.radiobutton6,'Visible','on','Value',0);
	set(handles.radiobutton7,'Visible','on','Value',0);
    
    
handles.output = hObject;
guidata(hObject, handles);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushbutton4_Callback(hObject, eventdata, handles)
    
    % answer finite element mesh
    handles.answer1 = questdlg(  '¿Do you want to smooth the velocities field before interpolating it to the mesh?','Question','Yes','No','No');
    switch handles.answer1
        case 'Yes'
            
            prompt = {'Select the kernel size to smooth the velocities (one odd number):'};
            dlgtitle = 'Kernel size';
            definput = {'3'};
            dims = [1 40];
            handles.kernel_size = str2double(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
            handles.id_kernel = 1;
        
        case 'No'   
            
            handles.id_kernel = 0;
            
    end
    
    % answer finite element mesh
    handles.answer2 = questdlg(  '¿Which interpolation method you want to use to interpolate the data to the finite element mesh?','Question','linear','nearest','cubic','cubic');
    switch handles.answer2
        case 'linear'
            
            handles.int_method = 'linear';
        
        case 'nearest'  
            
            handles.int_method = 'nearest';
            
        case 'cubic' 
            
            handles.int_method = 'cubic';
            
    end

    f = waitbar(0,'Please wait: Interpolating the velocity field ...');
    
    [handles.a,handles.b,handles.c,handles.d] = size(handles.MR_PCA_FH);
    
    if handles.id_kernel == 1 
        for n=1:handles.d
            handles.MR_PCA_FH(:,:,:,n)=smooth3(squeeze(handles.MR_PCA_FH(:,:,:,n)),'box',handles.kernel_size);
            handles.MR_PCA_AP(:,:,:,n)=smooth3(squeeze(handles.MR_PCA_AP(:,:,:,n)),'box',handles.kernel_size);
            handles.MR_PCA_RL(:,:,:,n)=smooth3(squeeze(handles.MR_PCA_RL(:,:,:,n)),'box',handles.kernel_size);
        end
    end
    
    nodes_n = [(handles.nodes(:,2)./handles.voxel_MR(2))+handles.voxel_MR(2)/2,(handles.nodes(:,1)./handles.voxel_MR(1))+handles.voxel_MR(1)/2,(handles.nodes(:,3)./handles.voxel_MR(3))+handles.voxel_MR(3)/2];
    
    [X,Y,Z] = meshgrid(1:handles.b,1:handles.a,1:handles.c);
    X=double(X);
    Y=double(Y);
    Z=double(Z);
    handles.veset = zeros(size(nodes_n(:,1),1),3,handles.d);
        
    for m=1:handles.d
        handles.veset(:,:,m)=[  interp3(X,Y,Z,squeeze(handles.MR_PCA_AP(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method),...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_FH(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1,...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_RL(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1];
    end
    
    handles.veset = handles.veset/100;    
    handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0; %% Julio Sotelo
    
    
    cc = sum(sqrt(sum(handles.veset.^2,2)),1);
    [C,I] = max(cc);
    handles.peak_flow = I;
    handles.min_vel = min(min(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.max_vel = max(max(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.mags_vol = squeeze(sqrt(sum(handles.veset.^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
    rgb_vel = cmap_vol/255;

    handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
    
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
    end
    
    close(f)
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
    set(handles.radiobutton4,'Visible','on','Value',1);
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
%     handles.nodes = [handles.nodes(:,1)-(2*handles.voxel_MR(1)),handles.nodes(:,2)-(2*handles.voxel_MR(2)),handles.nodes(:,3)-(2*handles.voxel_MR(3))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pp=1/100;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', 5/100,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    
    set(handles.text30,'visible','on')
    set(handles.text31,'visible','on')
    set(handles.text32,'visible','on')
    set(handles.pushbutton8,'visible','on')
    set(handles.pushbutton9,'visible','on')
    set(handles.pushbutton10,'visible','on')
    set(handles.pushbutton12,'visible','on')


handles.output = hObject;
guidata(hObject, handles);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton1_Callback(hObject, eventdata, handles)
    
    rb1 = get(handles.radiobutton1,'Value');
    if rb1 == 1
        set(handles.radiobutton2,'Value',0)
        set(handles.radiobutton3,'Value',0)
        set(handles.radiobutton4,'Value',0)
        cla(handles.axes1)
        axes(handles.axes1);
        patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        colorbar('off')
        
        set(handles.slider1,'Visible','off')
        set(handles.text30,'Visible','off')
        set(handles.text31,'Visible','off')
        set(handles.text32,'Visible','off')
        set(handles.pushbutton8,'Visible','off')
        set(handles.pushbutton9,'Visible','off')
        set(handles.pushbutton10,'Visible','off')
        set(handles.pushbutton12,'Visible','off')
               
        
    elseif rb1 == 0
        
        cla(handles.axes1)
        axes(handles.axes1);
        axis off
        colorbar('off')
        
    end

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton2_Callback(hObject, eventdata, handles)
    
    rb2 = get(handles.radiobutton2,'Value');
    if rb2 == 1
        set(handles.radiobutton1,'Value',0)
        set(handles.radiobutton3,'Value',0)
        set(handles.radiobutton4,'Value',0)
        cla(handles.axes1)
        axes(handles.axes1);
        data = smooth3(handles.SEG,'box',3);
        fv = isosurface(handles.xd,handles.yd,handles.zd,data,.5);
        p1 = patch(fv,'FaceColor',[mean(handles.R(handles.S==1)), mean(handles.G(handles.S==1)), mean(handles.B(handles.S==1))],'EdgeColor','k');
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        colorbar('off')
        
        set(handles.slider1,'Visible','off')
        set(handles.text30,'Visible','off')
        set(handles.text31,'Visible','off')
        set(handles.text32,'Visible','off')
        set(handles.pushbutton8,'Visible','off')
        set(handles.pushbutton9,'Visible','off')
        set(handles.pushbutton10,'Visible','off')
        set(handles.pushbutton12,'Visible','off')
        
    elseif rb2 == 0
        cla(handles.axes1)
        axes(handles.axes1);
        axis off
        colorbar('off')
    end
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton3_Callback(hObject, eventdata, handles)
    
    rb3 = get(handles.radiobutton3,'Value');
    if rb3 == 1
        set(handles.radiobutton1,'Value',0)
        set(handles.radiobutton2,'Value',0)
        set(handles.radiobutton4,'Value',0)
        cla(handles.axes1)
        axes(handles.axes1);
        patch('Vertices', handles.node_bin_aorta, 'Faces', handles.elem_bin_aorta,'FaceColor',[0.85 0.85 0.85],'EdgeColor','k');
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        colorbar('off')
        
        set(handles.slider1,'Visible','off')
        set(handles.text30,'Visible','off')
        set(handles.text31,'Visible','off')
        set(handles.text32,'Visible','off')
        set(handles.pushbutton8,'Visible','off')
        set(handles.pushbutton9,'Visible','off')
        set(handles.pushbutton10,'Visible','off')
        set(handles.pushbutton12,'Visible','off')
        
    elseif rb3 == 0
        cla(handles.axes1)
        axes(handles.axes1);
        axis off
        colorbar('off')
    end
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figure1_CloseRequestFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'waitstatus'),'waiting')
        handles.output = hObject;
        guidata(hObject, handles);
        uiresume(hObject);
    else
        delete(hObject);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function radiobutton4_Callback(hObject, eventdata, handles)
    
    rb4 = get(handles.radiobutton4,'Value');
    if rb4 == 1
        set(handles.radiobutton1,'Value',0)
        set(handles.radiobutton2,'Value',0)
        set(handles.radiobutton3,'Value',0)
        
        set(handles.slider1,'Visible','on')
        set(handles.text30,'Visible','on')
        set(handles.text31,'Visible','on')
        set(handles.text32,'Visible','on')
        set(handles.pushbutton8,'Visible','on')
        set(handles.pushbutton9,'Visible','on')
        set(handles.pushbutton10,'Visible','on')
        set(handles.pushbutton12,'Visible','on')
        
        cla(handles.axes1)
        axes(handles.axes1);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
        handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
        currentColormap = colormap(handles.axes1,jet(128));
        [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.mags(:));
        handles.max_vel = max(handles.mags(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Velocity [m/s]';
        c.FontSize = 12;
        c.Label.FontSize = 14;
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
    elseif rb4 == 0
        cla(handles.axes1)
        axes(handles.axes1);
        axis off
    end
    
handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function figure1_SizeChangedFcn(hObject, eventdata, handles)


    set(handles.text01, 'FontUnits', 'points');
    FontS = get(handles.text01, 'FontSize');
    set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);
    
    rb5 = get(handles.radiobutton5,'Value');
    rb6 = get(handles.radiobutton6,'Value');
    rb7 = get(handles.radiobutton7,'Value');

    if rb5==1

        d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
        d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
        d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
        s = (d1+d2+d3)/2;

        handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    

        set(handles.axes2,'visible','on')
        axes(handles.axes2);
        h = histfit(handles.surface_area,500);
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.FontWeight = 'bold';
        ax.FontSize = FontS-2;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Frequency';
        ax.XLabel.String = 'Area [mm^{2}]';

        g = get(h(2));
        xdata = g.XData';
        ydata = g.YData';
        [C,I] = max(ydata);
        txt1 = ['\leftarrow ', num2str(xdata(I)),' [mm^{2}]'];
        t = text(xdata(I),max(ydata),txt1);
        t.Color = 'red';
        t.FontWeight = 'bold';
        t.FontSize = FontS;

    elseif rb6==1

        % element volume    
        handles.elem_volume = elemvolume(handles.nodes(:,1:3),handles.elem);

        set(handles.axes2,'visible','on')
        axes(handles.axes2);
        h = histfit(handles.elem_volume,500);
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.FontWeight = 'bold';
        ax.FontSize = FontS-2;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Frequency';
        ax.XLabel.String = 'Volume [mm^{3}]';

        g = get(h(2));
        xdata = g.XData';
        ydata = g.YData';
        [C,I] = max(ydata);
        txt1 = ['\leftarrow ', num2str(xdata(I)),' [mm^{3}]'];
        t = text(xdata(I),max(ydata),txt1);
        t.Color = 'red';
        t.FontWeight = 'bold';
        t.FontSize = FontS;

    elseif rb7 ==1


        % mesh quality    
        handles.quality=meshquality(handles.nodes(:,1:3),handles.elem);

        set(handles.axes2,'visible','on')
        axes(handles.axes2);
        h = histfit(handles.quality,500);
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.FontWeight = 'bold';
        ax.FontSize = FontS-2;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Frequency';
        ax.XLabel.String = 'Quality [-]';

        g = get(h(2));
        xdata = g.XData';
        ydata = g.YData';
        [C,I] = max(ydata);
        txt1 = ['\leftarrow ', num2str(xdata(I))];
        t = text(xdata(I),max(ydata),txt1);
        t.Color = 'red';
        t.FontWeight = 'bold';
        t.FontSize = FontS;

    end


    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');

    set(handles.text01,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text30,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text31,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text32,'FontUnits','Normalized','FontSize',0.6)

    set(handles.text1,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text2,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text3,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text4,'FontUnits','Normalized','FontSize',0.55)

    set(handles.text9,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text5,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text6,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text7,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text8,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text11,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text10,'FontUnits','Normalized','FontSize',0.55)
    set(handles.text12,'FontUnits','Normalized','FontSize',0.6)
    set(handles.text13,'FontUnits','Normalized','FontSize',0.6)

    set(handles.edit1,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit2,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit3,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit4,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit5,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit6,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit7,'FontUnits','Normalized','FontSize',0.5)
    set(handles.edit8,'FontUnits','Normalized','FontSize',0.5)

    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton4,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton8,'FontUnits','Normalized','FontSize',0.57)
    set(handles.pushbutton9,'FontUnits','Normalized','FontSize',0.57)
    set(handles.pushbutton10,'FontUnits','Normalized','FontSize',0.57)
    set(handles.pushbutton12,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton3,'FontUnits','Normalized','FontSize',0.24)

    set(handles.radiobutton1,'FontUnits','Normalized','FontSize',0.4)
    set(handles.radiobutton2,'FontUnits','Normalized','FontSize',0.4)
    set(handles.radiobutton3,'FontUnits','Normalized','FontSize',0.4)
    set(handles.radiobutton4,'FontUnits','Normalized','FontSize',0.4)
    set(handles.radiobutton5,'FontUnits','Normalized','FontSize',0.48)
    set(handles.radiobutton6,'FontUnits','Normalized','FontSize',0.48)
    set(handles.radiobutton7,'FontUnits','Normalized','FontSize',0.48)

handles.output = hObject;
guidata(hObject, handles);    


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rb5 = get(handles.radiobutton5,'Value');

set(handles.text01, 'FontUnits', 'points');
FontS = get(handles.text01, 'FontSize');
set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);

if rb5 == 1
    
    rb6 = get(handles.radiobutton6,'Value');
	rb7 = get(handles.radiobutton7,'Value');
    
    if rb6 == 1
        set(handles.radiobutton6,'Value',0);
    end
    if rb7 == 1
        set(handles.radiobutton7,'Value',0);
    end
    
    % surface area
    d1 = sqrt(sum((handles.nodes(handles.faces(:,1),:)-handles.nodes(handles.faces(:,2),:)).^2,2));
    d2 = sqrt(sum((handles.nodes(handles.faces(:,2),:)-handles.nodes(handles.faces(:,3),:)).^2,2));
    d3 = sqrt(sum((handles.nodes(handles.faces(:,3),:)-handles.nodes(handles.faces(:,1),:)).^2,2));
    s = (d1+d2+d3)/2;
    
    handles.surface_area = sqrt(s.*(s-d1).*(s-d2).*(s-d3));    
    
    set(handles.axes2,'visible','on')
    axes(handles.axes2);
    h = histfit(handles.surface_area,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Area [mm^{2}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I)),' [mm^{2}]'];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str(sqrt(xdata(I)*4/sqrt(3)),3),' [mm]'],'visible','on');
    
elseif rb5 == 0
    
    rb6 = get(handles.radiobutton6,'Value');
	rb7 = get(handles.radiobutton7,'Value');
    
    if (rb6 + rb7) == 0
        cla(handles.axes2)
        set(handles.axes2,'visible','off')
        set(handles.text13,'visible','off');
    end
    
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb6 = get(handles.radiobutton6,'Value');

set(handles.text01, 'FontUnits', 'points');
FontS = get(handles.text01, 'FontSize');
set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);


if rb6 == 1
    
    rb5 = get(handles.radiobutton5,'Value');
	rb7 = get(handles.radiobutton7,'Value');
    
    if rb5 == 1
        set(handles.radiobutton5,'Value',0);
    end
    if rb7 == 1
        set(handles.radiobutton7,'Value',0);
    end
    
    % element volume    
    handles.elem_volume = elemvolume(handles.nodes(:,1:3),handles.elem);
    
    set(handles.axes2,'visible','on')
    axes(handles.axes2);
    h = histfit(handles.elem_volume,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Volume [mm^{3}]';
    
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I)),' [mm^{3}]'];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'string',['Characteristic Length: ',num2str((xdata(I)*12/sqrt(2))^(1/3),3),' [mm]'],'visible','on');
    
else
    
    rb5 = get(handles.radiobutton5,'Value');
	rb7 = get(handles.radiobutton7,'Value');
    
    if (rb5 + rb7) == 0
        cla(handles.axes2)
        set(handles.axes2,'visible','off')
        set(handles.text13,'visible','off');
    end
    
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb7 = get(handles.radiobutton7,'Value');

set(handles.text01, 'FontUnits', 'points');
FontS = get(handles.text01, 'FontSize');
set(handles.text01, 'FontUnits', 'normalized','FontSize',0.6);

if rb7 == 1
    
    rb5 = get(handles.radiobutton5,'Value');
	rb6 = get(handles.radiobutton6,'Value');
    
    if rb5 == 1
        set(handles.radiobutton5,'Value',0);
    end
    if rb6 == 1
        set(handles.radiobutton6,'Value',0);
    end
    
    % mesh quality    
    handles.quality=meshquality(handles.nodes(:,1:3),handles.elem);
    
    set(handles.axes2,'visible','on')
    axes(handles.axes2);
    h = histfit(handles.quality,500);
    axis square
    ax = gca;
    ax.XColor = 'w';
    ax.YColor = 'w';
    ax.FontWeight = 'bold';
    ax.FontSize = FontS-2;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridColor = 'k';
    ax.YLabel.String = 'Frequency';
    ax.XLabel.String = 'Quality [-]';
       
    g = get(h(2));
    xdata = g.XData';
    ydata = g.YData';
    [C,I] = max(ydata);
    txt1 = ['\leftarrow ', num2str(xdata(I))];
    t = text(xdata(I),max(ydata),txt1);
    t.Color = 'red';
    t.FontWeight = 'bold';
    t.FontSize = FontS;
    
    set(handles.text13,'visible','off');
    
else
    
    rb5 = get(handles.radiobutton5,'Value');
	rb6 = get(handles.radiobutton6,'Value');
    
    if (rb5 + rb6) == 0
        cla(handles.axes2)
        set(handles.axes2,'visible','off')
        set(handles.text13,'visible','off');
    end
    
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes during object creation, after setting all properties.
function radiobutton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radiobutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function radiobutton4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --------------------------------------------------------------------
function uitoggletool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if handles.cont_show == 1 % Julio Sotelo
        rotate3d on
        handles.cont_show = handles.cont_show + 1;
    else 
        rotate3d off
        handles.cont_show = 1;
    end
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

    pp=2/100;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = 100;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.vector_size = handles.slider_value;
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
    set(handles.radiobutton4,'Visible','on','Value',1);
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
    
    
    
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    f = waitbar(0,'Please wait: Interpolating the velocity field ...');
    
    handles.MR_PCA_AP = handles.MR_PCA_AP;
    handles.MR_PCA_FH = handles.MR_PCA_FH*-1;
    handles.MR_PCA_RL = handles.MR_PCA_RL;
    
    handles.MR_PCA_AP_ori = handles.MR_PCA_AP_ori; % change the direction of the velocity original
    handles.MR_PCA_FH_ori = handles.MR_PCA_FH_ori*-1; % change the direction of the  velocity original
    handles.MR_PCA_RL_ori = handles.MR_PCA_RL_ori; % change the direction of the velocity original

    [handles.a,handles.b,handles.c,handles.d] = size(handles.MR_PCA_FH);
    nodes_n = [(handles.nodes(:,2)./handles.voxel_MR(2))+handles.voxel_MR(2)/2,(handles.nodes(:,1)./handles.voxel_MR(1))+handles.voxel_MR(1)/2,(handles.nodes(:,3)./handles.voxel_MR(3))+handles.voxel_MR(3)/2];
    
    [X,Y,Z] = meshgrid(1:handles.b,1:handles.a,1:handles.c);
    X=double(X);
    Y=double(Y);
    Z=double(Z);
    handles.veset = zeros(size(nodes_n(:,1),1),3,handles.d);
    
    for m=1:handles.d
        handles.veset(:,:,m)=[  interp3(X,Y,Z,squeeze(handles.MR_PCA_AP(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method),...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_FH(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1,...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_RL(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1];
    end
    
    handles.veset = handles.veset/100;    
    handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0; %% Julio Sotelo
    
    
    cc = sum(sqrt(sum(handles.veset.^2,2)),1);
    [C,I] = max(cc);
    handles.peak_flow = I;
    handles.min_vel = min(min(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.max_vel = max(max(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.mags_vol = squeeze(sqrt(sum(handles.veset.^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
    rgb_vel = cmap_vol/255;
    
    handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
    end

    close(f)
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
   
    
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	f = waitbar(0,'Please wait: Interpolating the velocity field ...');
    
    handles.MR_PCA_AP = handles.MR_PCA_AP;
    handles.MR_PCA_FH = handles.MR_PCA_FH;
    handles.MR_PCA_RL = handles.MR_PCA_RL*-1;
    
    handles.MR_PCA_AP_ori = handles.MR_PCA_AP_ori; % change the direction of the velocity original
    handles.MR_PCA_FH_ori = handles.MR_PCA_FH_ori; % change the direction of the  velocity original
    handles.MR_PCA_RL_ori = handles.MR_PCA_RL_ori*-1; % change the direction of the velocity original

    [handles.a,handles.b,handles.c,handles.d] = size(handles.MR_PCA_FH);
    nodes_n = [(handles.nodes(:,2)./handles.voxel_MR(2))+handles.voxel_MR(2)/2,(handles.nodes(:,1)./handles.voxel_MR(1))+handles.voxel_MR(1)/2,(handles.nodes(:,3)./handles.voxel_MR(3))+handles.voxel_MR(3)/2];
    
    [X,Y,Z] = meshgrid(1:handles.b,1:handles.a,1:handles.c);
    X=double(X);
    Y=double(Y);
    Z=double(Z);
    handles.veset = zeros(size(nodes_n(:,1),1),3,handles.d);
    
    for m=1:handles.d
        handles.veset(:,:,m)=[  interp3(X,Y,Z,squeeze(handles.MR_PCA_AP(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method),...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_FH(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1,...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_RL(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1];
    end
    
    handles.veset = handles.veset/100;    
    handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0; %% Julio Sotelo
    
    
    cc = sum(sqrt(sum(handles.veset.^2,2)),1);
    [C,I] = max(cc);
    handles.peak_flow = I;
    handles.min_vel = min(min(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.max_vel = max(max(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.mags_vol = squeeze(sqrt(sum(handles.veset.^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
    rgb_vel = cmap_vol/255;
    handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
    end
    
    close(f)
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	f = waitbar(0,'Please wait: Interpolating the velocity field ...');
    
    handles.MR_PCA_AP = handles.MR_PCA_AP*-1;
    handles.MR_PCA_FH = handles.MR_PCA_FH;
    handles.MR_PCA_RL = handles.MR_PCA_RL;
    
    handles.MR_PCA_AP_ori = handles.MR_PCA_AP_ori*-1; % change the direction of the velocity original
    handles.MR_PCA_FH_ori = handles.MR_PCA_FH_ori; % change the direction of the  velocity original
    handles.MR_PCA_RL_ori = handles.MR_PCA_RL_ori; % change the direction of the velocity original
    
    [handles.a,handles.b,handles.c,handles.d] = size(handles.MR_PCA_FH);
    nodes_n = [(handles.nodes(:,2)./handles.voxel_MR(2))+handles.voxel_MR(2)/2,(handles.nodes(:,1)./handles.voxel_MR(1))+handles.voxel_MR(1)/2,(handles.nodes(:,3)./handles.voxel_MR(3))+handles.voxel_MR(3)/2];
    
    [X,Y,Z] = meshgrid(1:handles.b,1:handles.a,1:handles.c);
    X=double(X);
    Y=double(Y);
    Z=double(Z);
    handles.veset = zeros(size(nodes_n(:,1),1),3,handles.d);
        
    for m=1:handles.d
        handles.veset(:,:,m)=[  interp3(X,Y,Z,squeeze(handles.MR_PCA_AP(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method),...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_FH(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1,...
                                interp3(X,Y,Z,squeeze(handles.MR_PCA_RL(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1];
    end
    
    handles.veset = handles.veset/100;    
    handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0; %% Julio Sotelo
    
    
    cc = sum(sqrt(sum(handles.veset.^2,2)),1);
    [C,I] = max(cc);
    handles.peak_flow = I;
    handles.min_vel = min(min(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.max_vel = max(max(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
    handles.mags_vol = squeeze(sqrt(sum(handles.veset.^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
    rgb_vel = cmap_vol/255;
    handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
    end
    
    close(f)
    
    
    cla(handles.axes1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
    
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % answer finite element mesh        
    prompt = sprintf('¿How many voxels (thickness) do you want to smooth at the wall?\n (This is only applied for vWERP):');
    dlgtitle = 'Voxels (thickness)';
    definput = {'1'};
    dims = [1 80];
    out_dialog = inputdlg(prompt,dlgtitle,dims,definput);
    kernel_radius = str2double(cell2mat(out_dialog));
    nhood = strel("sphere",kernel_radius);
    handles.id_smooth_wall = 1;
    
    J1 = imdilate(handles.SEG,nhood.Neighborhood);
    J2 = imerode(handles.SEG,nhood.Neighborhood);
    handles.SEG_WALL = double((J1-J2)>0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smooth step decide how many 
    kernel_size = handles.kernel_size;
    
    MR_PCA_FH_temp = handles.MR_PCA_FH;
    MR_PCA_AP_temp = handles.MR_PCA_AP;
    MR_PCA_RL_temp = handles.MR_PCA_RL;
    
    for n=1:size(handles.MR_PCA_FH,4)
        MR_PCA_FH_temp(:,:,:,n) = smooth3(squeeze(handles.MR_PCA_FH(:,:,:,n)),'box',kernel_size);
        MR_PCA_AP_temp(:,:,:,n) = smooth3(squeeze(handles.MR_PCA_AP(:,:,:,n)),'box',kernel_size);
        MR_PCA_RL_temp(:,:,:,n) = smooth3(squeeze(handles.MR_PCA_RL(:,:,:,n)),'box',kernel_size);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new velocities
    
    % remove values from velocities
    MR_PCA_FH_old = handles.MR_PCA_FH.*repmat((abs(handles.SEG_WALL-1)),1,1,1,size(handles.MR_PCA_FH,4));
    MR_PCA_AP_old = handles.MR_PCA_AP.*repmat((abs(handles.SEG_WALL-1)),1,1,1,size(handles.MR_PCA_FH,4));
    MR_PCA_RL_old = handles.MR_PCA_RL.*repmat((abs(handles.SEG_WALL-1)),1,1,1,size(handles.MR_PCA_FH,4));
    
    % add values from wall smoothed
    MR_PCA_FH_old = MR_PCA_FH_old + MR_PCA_FH_temp.*repmat(handles.SEG_WALL,1,1,1,size(handles.MR_PCA_FH,4));
    MR_PCA_AP_old = MR_PCA_AP_old + MR_PCA_AP_temp.*repmat(handles.SEG_WALL,1,1,1,size(handles.MR_PCA_FH,4));
    MR_PCA_RL_old = MR_PCA_RL_old + MR_PCA_RL_temp.*repmat(handles.SEG_WALL,1,1,1,size(handles.MR_PCA_FH,4));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f = waitbar(0,'Please wait: Interpolating the velocity field ...');
    
    [handles.a,handles.b,handles.c,handles.d] = size(handles.MR_PCA_FH);
    nodes_n = [(handles.nodes(:,2)./handles.voxel_MR(2))+handles.voxel_MR(2)/2,(handles.nodes(:,1)./handles.voxel_MR(1))+handles.voxel_MR(1)/2,(handles.nodes(:,3)./handles.voxel_MR(3))+handles.voxel_MR(3)/2];
    
    [X,Y,Z] = meshgrid(1:handles.b,1:handles.a,1:handles.c);
    X=double(X);
    Y=double(Y);
    Z=double(Z);
    handles.veset_new = zeros(size(nodes_n(:,1),1),3,handles.d);
        
    for m=1:handles.d
        handles.veset_new(:,:,m)=[  interp3(X,Y,Z,squeeze(MR_PCA_AP_old(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method),...
                                interp3(X,Y,Z,squeeze(MR_PCA_FH_old(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1,...
                                interp3(X,Y,Z,squeeze(MR_PCA_RL_old(:,:,:,m)),nodes_n(:,1),nodes_n(:,2),nodes_n(:,3),handles.int_method)*-1];
    end
    
    handles.veset_new = handles.veset_new/100;    
    handles.veset_new(unique(handles.faces(:)),:,:) = handles.veset_new(unique(handles.faces(:)),:,:)*0; %% Julio Sotelo
    
    
    cc = sum(sqrt(sum(handles.veset_new.^2,2)),1);
    [C,I] = max(cc);
    handles.peak_flow = I;
    handles.min_vel = min(min(sqrt(sum(handles.veset_new(1:end,:,:).^2,2))));
    handles.max_vel = max(max(sqrt(sum(handles.veset_new(1:end,:,:).^2,2))));
    handles.mags_vol = squeeze(sqrt(sum(handles.veset_new.^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset_new,1) size(handles.veset_new,3)]);
    cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
    rgb_vel = cmap_vol/255;

    handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
    
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
    end
    
    
    cla(handles.axes1)
    delete('handles.cmap')
    pause(0.1)
    axes(handles.axes1);
    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
    hold on
    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset_new(1:end,1,handles.peak_flow),handles.veset_new(1:end,2,handles.peak_flow),handles.veset_new(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
    handles.mags = squeeze(sqrt(sum(handles.veset_new(1:end,:,:).^2,2)));
    currentColormap = colormap(handles.axes1,jet(128));
    [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
    ind = reshape(ind,[size(handles.veset_new(1:end,1,1),1) size(handles.veset_new,3)]);
    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(handles.mags(:));
    handles.max_vel = max(handles.mags(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [1 1 1];
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 12;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])
    set(handles.radiobutton4,'Visible','on','Value',1);
    set(handles.radiobutton1,'Value',0);
    set(handles.radiobutton2,'Value',0);
    set(handles.radiobutton3,'Value',0);
%     handles.nodes = [handles.nodes(:,1)-(2*handles.voxel_MR(1)),handles.nodes(:,2)-(2*handles.voxel_MR(2)),handles.nodes(:,3)-(2*handles.voxel_MR(3))];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pp=1/100;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', 5/100,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    
    set(handles.text30,'visible','on')
    set(handles.text31,'visible','on')
    set(handles.text32,'visible','on')
    set(handles.pushbutton8,'visible','on')
    set(handles.pushbutton9,'visible','on')
    set(handles.pushbutton10,'visible','on')
    set(handles.pushbutton12,'visible','on')

    close(f)

    % answer finite element mesh
    handles.answer1 = questdlg(  '¿Do you agree with the smoothing result?','Question','Yes','No','No');
    switch handles.answer1
        case 'Yes'
            
            handles.MR_PCA_FH = MR_PCA_FH_old;
            handles.MR_PCA_RL = MR_PCA_RL_old;
            handles.MR_PCA_AP = MR_PCA_AP_old;
            handles.veset = handles.veset_new;
            handles.id_mod = 1;
        
        case 'No'   

            f = waitbar(0,'Please wait: Velocities has not been changed ...');

            handles.id_mod = 0;
            cc = sum(sqrt(sum(handles.veset.^2,2)),1);
            [C,I] = max(cc);
            handles.peak_flow = I;
            handles.min_vel = min(min(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
            handles.max_vel = max(max(sqrt(sum(handles.veset(1:end,:,:).^2,2))));
            handles.mags_vol = squeeze(sqrt(sum(handles.veset.^2,2)));
            currentColormap = colormap(handles.axes1,jet(128));
            [~, ~, ind] = histcounts(handles.mags_vol(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            cmap_vol = uint8(ind2rgb(ind, currentColormap) * 255);
            rgb_vel = cmap_vol/255;
        
            handles.Lrgb_vel = ones(size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),size(handles.MR_PCA_FH,4),3);
            
            MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[1,2,3]);
            MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
            xd_seg = MASK.*handles.xd;
            yd_seg = MASK.*handles.yd;
            zd_seg = MASK.*handles.zd;
            xd_seg(MASK==0) = [];
            yd_seg(MASK==0) = [];
            zd_seg(MASK==0) = [];
            pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
            [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,3));
            X_seg = MASK.*X;
            Y_seg = MASK.*Y;
            Z_seg = MASK.*Z;
            X_seg(MASK==0) = [];
            Y_seg(MASK==0) = [];
            Z_seg(MASK==0) = [];
            pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
            for n=1:length(pos_voxel(:,1))
                d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
                handles.Lrgb_vel(pos_v(n,2),pos_v(n,1),pos_v(n,3),:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
            end

            cla(handles.axes1)
            delete('handles.cmap')
            pause(0.1)
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),handles.vector_size,'Linewidth',1);
            handles.mags = squeeze(sqrt(sum(handles.veset(1:end,:,:).^2,2)));
            currentColormap = colormap(handles.axes1,jet(128));
            [~, ~, ind] = histcounts(handles.mags(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset(1:end,1,1),1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vel = min(handles.mags(:));
            handles.max_vel = max(handles.mags(:));
            handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vel handles.max_vel];
            c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
            c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.2 0.02 0.3];
            c.FontWeight = 'bold';
            c.Label.String = 'Velocity [m/s]';
            c.FontSize = 12;
            c.Label.FontSize = 14;
            c.Label.FontWeight = 'bold';
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_vel handles.max_vel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])

            close(f)
    end


handles.output = hObject;
guidata(hObject, handles);
