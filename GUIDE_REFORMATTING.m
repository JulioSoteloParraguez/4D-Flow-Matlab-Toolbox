function varargout = GUIDE_REFORMATTING(varargin)
% GUIDE_REFORMATTING MATLAB code for GUIDE_REFORMATTING.fig
%      GUIDE_REFORMATTING, by itself, creates a new GUIDE_REFORMATTING or raises the existing
%      singleton*.
%
%      H = GUIDE_REFORMATTING returns the handle to a new GUIDE_REFORMATTING or the handle to
%      the existing singleton*.
%
%      GUIDE_REFORMATTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_REFORMATTING.M with the given input arguments.
%
%      GUIDE_REFORMATTING('Property','Value',...) creates a new GUIDE_REFORMATTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_REFORMATTING_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_REFORMATTING_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_REFORMATTING

% Last Modified by GUIDE v2.5 21-Sep-2022 16:58:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_REFORMATTING_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_REFORMATTING_OutputFcn, ...
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


% --- Executes just before GUIDE_REFORMATTING is made visible.
function GUIDE_REFORMATTING_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE_REFORMATTING (see VARARGIN)

% Choose default command line output for GUIDE_REFORMATTING
handles.output = hObject;

    path(path,'IO_CODES/')

% folder_name = uigetdir([],'Load Folder...');
%     disp(folder_name)
% 
%     % READING MATLAB FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     h = msgbox({'Please wait ...','Loading structure data ...'});
%     
%     load([folder_name,'\data.mat']);



    
%     handles.VENC        = varargin{1}.VENC;
%     handles.voxel_MR    = varargin{1}.voxel_MR;
%     handles.heart_rate  = varargin{1}.heart_rate;
%     handles.MR_FFE_FH   = varargin{1}.MR_FFE_FH;
%     handles.MR_FFE_AP   = varargin{1}.MR_FFE_AP;
%     handles.MR_FFE_RL   = varargin{1}.MR_FFE_RL;
%     handles.MR_PCA_FH   = varargin{1}.MR_PCA_FH;
%     handles.MR_PCA_AP   = varargin{1}.MR_PCA_AP;
%     handles.MR_PCA_RL   = varargin{1}.MR_PCA_RL;

    handles.VENC                = varargin{1}.VENC_l;
    handles.voxel_MR            = varargin{1}.voxel_MR_l;
    handles.heart_rate          = varargin{1}.heart_rate_l;
%     handles.MR_FFE_FH           = varargin{1}.MR_FFE_FH_l;
%     handles.MR_FFE_AP           = varargin{1}.MR_FFE_AP_l;
%     handles.MR_FFE_RL           = varargin{1}.MR_FFE_RL_l;
%     handles.MR_PCA_FH           = varargin{1}.MR_PCA_FH_l;
%     handles.MR_PCA_AP           = varargin{1}.MR_PCA_AP_l;
%     handles.MR_PCA_RL           = varargin{1}.MR_PCA_RL_l;
    handles.IPCMRA_ORI_COOR     = varargin{1}.IPCMRA_ORI_COOR;
    handles.FFE_ORI_COOR        = varargin{1}.FFE_ORI_COOR;
    handles.MR_FFE_FH_ORI_COOR  = varargin{1}.MR_FFE_FH_ORI_COOR;
    handles.MR_PCA_FH_ORI_COOR  = varargin{1}.MR_PCA_FH_ORI_COOR;
    handles.MR_PCA_AP_ORI_COOR  = varargin{1}.MR_PCA_AP_ORI_COOR;
    handles.MR_PCA_RL_ORI_COOR  = varargin{1}.MR_PCA_RL_ORI_COOR;
    handles.min_voxel           = varargin{1}.min_voxel;
    handles.MR_FFE_FH           = varargin{1}.MR_FFE_FH_resize;
    handles.MR_FFE_AP           = varargin{1}.MR_FFE_AP_resize;
    handles.MR_FFE_RL           = varargin{1}.MR_FFE_RL_resize;
    handles.MR_PCA_FH           = varargin{1}.MR_PCA_FH_resize;
    handles.MR_PCA_AP           = varargin{1}.MR_PCA_AP_resize;
    handles.MR_PCA_RL           = varargin{1}.MR_PCA_RL_resize;
    handles.IPCMRA              = varargin{1}.IPCMRA_resize;
    handles.FFE                 = varargin{1}.FFE;
    handles.pos_a               = varargin{1}.pos_a;
    handles.pos_b               = varargin{1}.pos_b;
    handles.pos_c               = varargin{1}.pos_c;
    handles.corte_sagital_x     = varargin{1}.corte_sagital_x;
    handles.corte_sagital_y     = varargin{1}.corte_sagital_y;
    handles.corte_sagital_z     = varargin{1}.corte_sagital_z;
    handles.corte_coronal_x     = varargin{1}.corte_coronal_x;
    handles.corte_coronal_y     = varargin{1}.corte_coronal_y;
    handles.corte_coronal_z     = varargin{1}.corte_coronal_z;
    handles.corte_axial_x       = varargin{1}.corte_axial_x;
    handles.corte_axial_y       = varargin{1}.corte_axial_y;
    handles.corte_axial_z       = varargin{1}.corte_axial_z;
    handles.IPCMRA_1            = varargin{1}.IPCMRA_1;
    handles.IPCMRA_2            = varargin{1}.IPCMRA_2;
    handles.IPCMRA_3            = varargin{1}.IPCMRA_3;
    handles.FFE_1               = varargin{1}.FFE_1;
    handles.FFE_2               = varargin{1}.FFE_2;
    handles.FFE_3               = varargin{1}.FFE_3;
    handles.coordx              = varargin{1}.coordx;
    handles.coordy              = varargin{1}.coordy;
    handles.coordz              = varargin{1}.coordz;

    popup5_id                   = varargin{1}.popup5_id;
    popup6_id                   = varargin{1}.popup6_id;
    popup7_id                   = varargin{1}.popup7_id;
    popup8_id                   = varargin{1}.popup8_id;
    
    handles.id_while = 0;
    
    set(handles.axes1,'visible','on')
    set(handles.axes2,'visible','on')
    set(handles.axes3,'visible','on')
    set(handles.axes4,'visible','on')

    list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
    set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);

    list_string2 = {'IPCMRA','MAGNITUDE'};
    set(handles.popupmenu5,'visible','on','String',list_string2,'value',popup5_id);
    set(handles.popupmenu6,'visible','on','String',list_string2,'value',popup6_id);
    set(handles.popupmenu7,'visible','on','String',list_string2,'value',popup7_id);
    set(handles.popupmenu8,'visible','on','String',list_string2,'value',popup8_id);

    list_string3 = {'PLANES','PLANES + Sagital View','PLANES + Coronal View','PLANES + Axial View'};
    set(handles.popupmenu4,'visible','on','String',list_string3,'value',1);

    if get(handles.popupmenu5,'Value')==1
        
        axes(handles.axes1);
        axis off
        imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
        hold on
        images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes1,'Visible','off');
        
    elseif get(handles.popupmenu5,'Value')==2
        
        axes(handles.axes1);
        axis off
        imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
        hold on
        images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes1,'Visible','off');
        
    end
    
    if get(handles.popupmenu6,'Value')==1
        
        axes(handles.axes2);
        axis off
        imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
        hold on
        images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes2,'Visible','off');
        
    elseif get(handles.popupmenu6,'Value')==2
    
        axes(handles.axes2);
        axis off
        imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
        hold on
        images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes2,'Visible','off');
        
    end

    if get(handles.popupmenu7,'Value')==1

        axes(handles.axes3);
        axis off
        imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
        hold on
        images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes3,'Visible','off');
        
    elseif get(handles.popupmenu7,'Value')==2
        
        axes(handles.axes3);
        axis off
        imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
        hold on
        images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes3,'Visible','off');

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mostramos los cortes a la salida

    if get(handles.popupmenu8,'Value')==1

        axes(handles.axes4)
        plot(0,0)
        surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        hold on
        surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        hold off
        axis image
        colormap gray
        grid on
        axis on
        xlabel('F-H')
        ylabel('A-P')
        zlabel('R-L')
        view(35,45)
        daspect([1 1 1])
        axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

    elseif get(handles.popupmenu8,'Value')==2

        handles.FFE_1 = (handles.FFE_1/max(handles.FFE_1(:)))*255;
        handles.FFE_2 = (handles.FFE_2/max(handles.FFE_2(:)))*255;
        handles.FFE_3 = (handles.FFE_3/max(handles.FFE_3(:)))*255;

        axes(handles.axes4)
        plot(0,0)
        surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.FFE_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        hold on
        surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.FFE_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.FFE_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        hold off
        axis image
        colormap gray
        grid on
        axis on
        xlabel('F-H')
        ylabel('A-P')
        zlabel('R-L')
        view(35,45)
        daspect([1 1 1])
        axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

    end

    set(handles.uipanel1,'visible','on')
    set(handles.uipanel2,'visible','on')
    set(handles.uipanel3,'visible','on')
    set(handles.uipanel4,'visible','on')

    set(handles.pushbutton1,'visible','on')
    set(handles.pushbutton2,'visible','on')
    set(handles.pushbutton3,'visible','on')

    


% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);

% UIWAIT makes GUIDE_REFORMATTING wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_REFORMATTING_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if handles.id_while == 0
    varargout{1} = handles.output;
    id_while = handles.id_while;
    MR_FFE_FH_resize = handles.MR_FFE_FH;
    MR_FFE_AP_resize = handles.MR_FFE_AP;
    MR_FFE_RL_resize = handles.MR_FFE_RL;
    MR_PCA_FH_resize = handles.MR_PCA_FH;
    MR_PCA_AP_resize = handles.MR_PCA_AP;
    MR_PCA_RL_resize = handles.MR_PCA_RL;
    IPCMRA_resize = handles.IPCMRA;
    FFE = handles.FFE;
    pos_a = handles.pos_a;
    pos_b = handles.pos_b;
    pos_c = handles.pos_c;
    coordx = handles.coordx;
    coordy = handles.coordy;
    coordz = handles.coordz;

    popup5_id = get(handles.popupmenu5,'Value');
    popup6_id = get(handles.popupmenu6,'Value');
    popup7_id = get(handles.popupmenu7,'Value');
    popup8_id = get(handles.popupmenu8,'Value');

    setappdata(0,'id_while',id_while);
    setappdata(0,'MR_FFE_FH_resize',MR_FFE_FH_resize);
    setappdata(0,'MR_FFE_AP_resize',MR_FFE_AP_resize);
    setappdata(0,'MR_FFE_RL_resize',MR_FFE_RL_resize);
    setappdata(0,'MR_PCA_FH_resize',MR_PCA_FH_resize);
    setappdata(0,'MR_PCA_AP_resize',MR_PCA_AP_resize);
    setappdata(0,'MR_PCA_RL_resize',MR_PCA_RL_resize);
    setappdata(0,'IPCMRA_resize',IPCMRA_resize);
    setappdata(0,'FFE',FFE);
    setappdata(0,'pos_a',pos_a);
    setappdata(0,'pos_b',pos_b);
    setappdata(0,'pos_c',pos_c);
    setappdata(0,'coordx',coordx);
    setappdata(0,'coordy',coordy);
    setappdata(0,'coordz',coordz);
    setappdata(0,'popup5_id',popup5_id);
    setappdata(0,'popup6_id',popup6_id);
    setappdata(0,'popup7_id',popup7_id);
    setappdata(0,'popup8_id',popup8_id);

elseif handles.id_while == 1
    varargout{1} = handles.output;
    id_while = handles.id_while;
    MR_FFE_FH_resize = handles.MR_FFE_FH;
    MR_FFE_AP_resize = handles.MR_FFE_AP;
    MR_FFE_RL_resize = handles.MR_FFE_RL;
    MR_PCA_FH_resize = handles.MR_PCA_FH;
    MR_PCA_AP_resize = handles.MR_PCA_AP;
    MR_PCA_RL_resize = handles.MR_PCA_RL;
    IPCMRA_resize = handles.IPCMRA;
    FFE = handles.FFE;
    pos_a = handles.pos_a;
    pos_b = handles.pos_b;
    pos_c = handles.pos_c;
    coordx = handles.coordx;
    coordy = handles.coordy;
    coordz = handles.coordz;

    popup5_id = get(handles.popupmenu5,'Value');
    popup6_id = get(handles.popupmenu6,'Value');
    popup7_id = get(handles.popupmenu7,'Value');
    popup8_id = get(handles.popupmenu8,'Value');

    setappdata(0,'id_while',id_while);
    setappdata(0,'MR_FFE_FH_resize',MR_FFE_FH_resize);
    setappdata(0,'MR_FFE_AP_resize',MR_FFE_AP_resize);
    setappdata(0,'MR_FFE_RL_resize',MR_FFE_RL_resize);
    setappdata(0,'MR_PCA_FH_resize',MR_PCA_FH_resize);
    setappdata(0,'MR_PCA_AP_resize',MR_PCA_AP_resize);
    setappdata(0,'MR_PCA_RL_resize',MR_PCA_RL_resize);
    setappdata(0,'IPCMRA_resize',IPCMRA_resize);
    setappdata(0,'FFE',FFE);
    setappdata(0,'pos_a',pos_a);
    setappdata(0,'pos_b',pos_b);
    setappdata(0,'pos_c',pos_c);
    setappdata(0,'coordx',coordx);
    setappdata(0,'coordy',coordy);
    setappdata(0,'coordz',coordz);
    setappdata(0,'popup5_id',popup5_id);
    setappdata(0,'popup6_id',popup6_id);
    setappdata(0,'popup7_id',popup7_id);
    setappdata(0,'popup8_id',popup8_id);
    delete(handles.figure1);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

switch get(handles.popupmenu1,'Value')

    case 2

        if get(handles.popupmenu5,'Value')==1
            
            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
            hold on
            ROI = images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','c', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');


            center = [handles.pos_b,handles.pos_a,0];
            position = customWait(ROI);
            position = [round(position(1)),round(position(2)),0];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');


            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');
            
        elseif get(handles.popupmenu5,'Value')==2
            
            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
            hold on
            ROI = images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','c', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');


            center = [handles.pos_b,handles.pos_a,0];
            position = customWait(ROI);
            position = [round(position(1)),round(position(2)),0];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');


            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');
           
        end

    case 3
        if get(handles.popupmenu5,'Value')==1
            
            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes1,'Position',[handles.pos_b, handles.pos_a; handles.pos_b, round((handles.pos_a + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');

            vect_ini = [handles.pos_b, handles.pos_a] - [handles.pos_b, round((handles.pos_a + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_b)*-1;


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 0 1],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 0 1],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 0 1],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 0 1],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 0 1],"linear","crop");

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

            else

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        elseif get(handles.popupmenu5,'Value')==2
            
            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes1,'Position',[handles.pos_b, handles.pos_a; handles.pos_b, round((handles.pos_a + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');

            vect_ini = [handles.pos_b, handles.pos_a] - [handles.pos_b, round((handles.pos_a + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_b)*-1;


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 0 1],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 0 1],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 0 1],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 0 1],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 0 1],"linear","crop");

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

            else

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        end

    case 4
        
        if get(handles.popupmenu5,'Value')==1

            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes1,'Position',[handles.pos_b, handles.pos_a; round((handles.pos_b + size(handles.IPCMRA,1))/2),handles.pos_a]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');

            vect_ini = [handles.pos_b, handles.pos_a] - [round((handles.pos_b + size(handles.IPCMRA,1))/2),handles.pos_a];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_a);

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 0 1],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 0 1],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 0 1],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 0 1],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 0 1],"linear","crop");

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

            else

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        elseif get(handles.popupmenu5,'Value')==2
            
            axes(handles.axes1);
            axis off
            imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
            hold on
            images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes1,'Position',[handles.pos_b, handles.pos_a; round((handles.pos_b + size(handles.IPCMRA,1))/2),handles.pos_a]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes1,'Visible','off');

            vect_ini = [handles.pos_b, handles.pos_a] - [round((handles.pos_b + size(handles.IPCMRA,1))/2),handles.pos_a];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_a);

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 0 1],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 0 1],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 0 1],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 0 1],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 0 1],"linear","crop");

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

            else

                axes(handles.axes1);
                axis off
                imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                hold on
                images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes1,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
        end
end

if get(handles.popupmenu6,'Value')==1
    
    axes(handles.axes2);
    axis off
    imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
    hold on
    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes2,'Visible','off');
    
elseif get(handles.popupmenu6,'Value')==2
    
    axes(handles.axes2);
    axis off
    imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
    hold on
    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes2,'Visible','off');
    
end
    
if get(handles.popupmenu7,'Value')==1
    
    axes(handles.axes3);
    axis off
    imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
    hold on
    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes3,'Visible','off');
    
elseif get(handles.popupmenu7,'Value')==2
    
    axes(handles.axes3);
    axis off
    imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
    hold on
    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes3,'Visible','off');
    
end

    list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
    set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        point = [];
        point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
        point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
        point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);
        
        point = point*handles.min_voxel;
        point(1) = point(1)/handles.voxel_MR(1);
        point(2) = point(2)/handles.voxel_MR(2);
        point(3) = point(3)/handles.voxel_MR(3);
        
        p_norm = [];
        p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c+4);
        p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c+4);
        p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c+4);

        p_norm = p_norm*handles.min_voxel;
        p_norm(1) = p_norm(1)/handles.voxel_MR(1);
        p_norm(2) = p_norm(2)/handles.voxel_MR(2);
        p_norm(3) = p_norm(3)/handles.voxel_MR(3);
        
        normal = p_norm - point;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mostramos los cortes a la salida
    
        [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);
        [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);
        MR_FFE_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_AP_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_RL_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        
        for dt = 1:size(handles.MR_FFE_FH_ORI_COOR,4)
            MR_FFE_FH_2D(:,:,dt) = obliqueslice(handles.MR_FFE_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_FH_2D(:,:,dt) = obliqueslice(handles.MR_PCA_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_AP_2D(:,:,dt) = obliqueslice(handles.MR_PCA_AP_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_RL_2D(:,:,dt) = obliqueslice(handles.MR_PCA_RL_ORI_COOR(:,:,:,dt),point,normal);
        end

        if get(handles.popupmenu8,'Value')==1
            
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        elseif get(handles.popupmenu8,'Value')==2
           
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyectamos la velocidad a travez del plano
%         disp(min(MR_FFE_FH_2D(:)))
%         disp(min(MR_PCA_FH_2D(:)))
%         disp(min(MR_PCA_AP_2D(:)))
%         disp(min(MR_PCA_RL_2D(:)))
%         
%         disp(max(MR_FFE_FH_2D(:)))
%         disp(max(MR_PCA_FH_2D(:)))
%         disp(max(MR_PCA_AP_2D(:)))
%         disp(max(MR_PCA_RL_2D(:)))
        
        % normal vector
        nomal_vector = normal/norm(normal);
        
        % projection of the vector
        magnitud = MR_FFE_FH_2D;
        velocity = MR_PCA_FH_2D*nomal_vector(1) + MR_PCA_AP_2D*nomal_vector(2) + MR_PCA_RL_2D*nomal_vector(3);
        
%         figure, 
%         subplot 121, imshow(magnitud(:,:,5),[])
%         subplot 122, imshow(velocity(:,:,5),[])
        
        %%% coordinates        
        pp = zeros(size(cx,1),size(cx,2),3);
        pp(:,:,1) = cx*handles.voxel_MR(1);
        pp(:,:,2) = cy*handles.voxel_MR(2);
        pp(:,:,3) = cz*handles.voxel_MR(3);
        
        x_coor = zeros(size(cx));
        y_coor = zeros(size(cy));
        
        x_coor(:,2:end) = sqrt(sum((pp(:,2:end,:) - pp(:,1:end-1,:)).^2,3));
        y_coor(2:end,:) = sqrt(sum((pp(2:end,:,:) - pp(1:end-1,:,:)).^2,3));
        
        pixel_size = [x_coor(2,2),y_coor(2,2)];
        disp(pixe_size)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % llamamos a la funcion para escribir el archivo
        structure = create_structure(magnitud, velocity, handles.heart_rate, handles.VENC, pixel_size);
        
        setstruct = structure.setstruct;
        preview = structure.preview;
        info = structure.info;
        im = structure.im;
        
        directory = uigetdir(pwd, 'Select Directory');
        save([directory,'\structure.mat'],'setstruct','preview','info','im')
        
        f = msgbox("The data has been successfully saved","Success");
        pause(2)
        close(f)

handles.output = hObject;
guidata(hObject, handles);



% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.popupmenu5,'Value')

    case 1
        
        axes(handles.axes1);
        axis off
        imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
        hold on
        images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes1,'Visible','off');
        
    case 2
        
        axes(handles.axes1);
        axis off
        imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
        hold on
        images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes1,'Visible','off');
        
end

        
handles.output = hObject;
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

switch get(handles.popupmenu2,'Value')
    case 2

        if get(handles.popupmenu6,'Value')==1
            
            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
            hold on
            ROI = images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','cyan', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');

            center = [handles.pos_b,0,handles.pos_c];
            position = customWait(ROI);
            position = [round(position(2)),0,round(position(1))];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');

            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');
            
        elseif get(handles.popupmenu6,'Value')==2
            
            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
            hold on
            ROI = images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','cyan', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');

            center = [handles.pos_b,0,handles.pos_c];
            position = customWait(ROI);
            position = [round(position(2)),0,round(position(1))];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');

            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');
            
        end

    case 3

        if get(handles.popupmenu6,'Value')==1

            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes2,'Position',[handles.pos_c, handles.pos_b; handles.pos_c, round((handles.pos_b + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');

            vect_ini = [handles.pos_c, handles.pos_b] - [handles.pos_c, round((handles.pos_b + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_c)*-1;

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 1 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 1 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 1 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 1 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 1 0],"linear","crop");

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

            else

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        elseif get(handles.popupmenu6,'Value')==2
            
            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes2,'Position',[handles.pos_c, handles.pos_b; handles.pos_c, round((handles.pos_b + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');

            vect_ini = [handles.pos_c, handles.pos_b] - [handles.pos_c, round((handles.pos_b + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_c)*-1;

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 1 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 1 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 1 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 1 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 1 0],"linear","crop");

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

            else

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        end

    case 4

        if get(handles.popupmenu6,'Value')==1
        
            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes2,'Position',[handles.pos_c, handles.pos_b; round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_b]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');


            vect_ini = [handles.pos_c, handles.pos_b] - [round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_b];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_b);


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 1 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 1 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 1 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 1 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 1 0],"linear","crop");

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

            else

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        elseif get(handles.popupmenu6,'Value')==2
            
            axes(handles.axes2);
            axis off
            imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
            hold on
            images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            ROI = images.roi.Line(handles.axes2,'Position',[handles.pos_c, handles.pos_b; round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_b]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes2,'Visible','off');


            vect_ini = [handles.pos_c, handles.pos_b] - [round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_b];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_b);


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[0 1 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[0 1 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[0 1 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[0 1 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[0 1 0],"linear","crop");

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

            else

                axes(handles.axes2);
                axis off
                imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                hold on
                images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes2,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
        end

end

if get(handles.popupmenu5,'Value')==1
    
    axes(handles.axes1);
    axis off
    imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
    hold on
    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes1,'Visible','off');

elseif get(handles.popupmenu5,'Value')==2
   
    axes(handles.axes1);
    axis off
    imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
    hold on
    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes1,'Visible','off');
    
end
    
if get(handles.popupmenu7,'Value')==1
    
    axes(handles.axes3);
    axis off
    imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
    hold on
    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes3,'Visible','off');
    
elseif get(handles.popupmenu7,'Value')==2
    
    axes(handles.axes3);
    axis off
    imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
    hold on
    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes3,'Visible','off');

end

    list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
    set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);


handles.output = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        point = [];
        point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
        point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
        point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

        point = point*handles.min_voxel;
        point(1) = point(1)/handles.voxel_MR(1);
        point(2) = point(2)/handles.voxel_MR(2);
        point(3) = point(3)/handles.voxel_MR(3);
        
        p_norm = [];
        p_norm(1) = handles.coordx(handles.pos_a+4,handles.pos_b,handles.pos_c);
        p_norm(2) = handles.coordy(handles.pos_a+4,handles.pos_b,handles.pos_c);
        p_norm(3) = handles.coordz(handles.pos_a+4,handles.pos_b,handles.pos_c);

        p_norm = p_norm*handles.min_voxel;
        p_norm(1) = p_norm(1)/handles.voxel_MR(1);
        p_norm(2) = p_norm(2)/handles.voxel_MR(2);
        p_norm(3) = p_norm(3)/handles.voxel_MR(3);

%         disp(point)
%         disp(p_norm)
        normal = p_norm - point;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mostramos los cortes a la salida
    
        [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);
        [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);
        MR_FFE_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_AP_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_RL_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        
        for dt = 1:size(handles.MR_FFE_FH_ORI_COOR,4)
            MR_FFE_FH_2D(:,:,dt) = obliqueslice(handles.MR_FFE_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_FH_2D(:,:,dt) = obliqueslice(handles.MR_PCA_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_AP_2D(:,:,dt) = obliqueslice(handles.MR_PCA_AP_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_RL_2D(:,:,dt) = obliqueslice(handles.MR_PCA_RL_ORI_COOR(:,:,:,dt),point,normal);
        end

        if get(handles.popupmenu8,'Value')==1
            
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        elseif get(handles.popupmenu8,'Value')==2
           
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyectamos la velocidad a travez del plano
%         disp(min(MR_FFE_FH_2D(:)))
%         disp(min(MR_PCA_FH_2D(:)))
%         disp(min(MR_PCA_AP_2D(:)))
%         disp(min(MR_PCA_RL_2D(:)))
%         
%         disp(max(MR_FFE_FH_2D(:)))
%         disp(max(MR_PCA_FH_2D(:)))
%         disp(max(MR_PCA_AP_2D(:)))
%         disp(max(MR_PCA_RL_2D(:)))
        
        % normal vector
        nomal_vector = normal/norm(normal);
        
        % projection of the vector
        magnitud = MR_FFE_FH_2D;
        velocity = MR_PCA_FH_2D*nomal_vector(1) + MR_PCA_AP_2D*nomal_vector(2) + MR_PCA_RL_2D*nomal_vector(3);
        
%         figure, 
%         subplot 121, imshow(magnitud(:,:,5),[])
%         subplot 122, imshow(velocity(:,:,5),[])
        
        %%% coordinates        
        pp = zeros(size(cx,1),size(cx,2),3);
        pp(:,:,1) = cx*handles.voxel_MR(1);
        pp(:,:,2) = cy*handles.voxel_MR(2);
        pp(:,:,3) = cz*handles.voxel_MR(3);
        
        x_coor = zeros(size(cx));
        y_coor = zeros(size(cy));
        
        x_coor(:,2:end) = sqrt(sum((pp(:,2:end,:) - pp(:,1:end-1,:)).^2,3));
        y_coor(2:end,:) = sqrt(sum((pp(2:end,:,:) - pp(1:end-1,:,:)).^2,3));
        
        pixel_size = [x_coor(2,2),y_coor(2,2)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % llamamos a la funcion para escribir el archivo
        structure = create_structure(magnitud, velocity, handles.heart_rate, handles.VENC, pixel_size);
        
        setstruct = structure.setstruct;
        preview = structure.preview;
        info = structure.info;
        im = structure.im;
        
        directory = uigetdir(pwd, 'Select Directory');
        save([directory,'\structure.mat'],'setstruct','preview','info','im')
        

        f = msgbox("The data has been successfully saved","Success");
        pause(2)
        close(f)
        
handles.output = hObject;
guidata(hObject, handles);
        
        
        
        
% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6

switch get(handles.popupmenu6,'Value')
    case 1
        
        axes(handles.axes2);
        axis off
        imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
        hold on
        images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes2,'Visible','off');
        
     case 2
        
        axes(handles.axes2);
        axis off
        imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
        hold on
        images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes2,'Visible','off');
end
            
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

switch get(handles.popupmenu3,'Value')
    case 2

        if get(handles.popupmenu7,'Value')==1
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
            hold on
            ROI = images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','cyan', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');


            center = [0,handles.pos_a,handles.pos_c];
            position = customWait(ROI);
            position = [0,round(position(2)),round(position(1))];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');

            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');
        
        elseif get(handles.popupmenu7,'Value')==2
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
            hold on
            ROI = images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','cyan', 'InteractionsAllowed', 'all');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');


            center = [0,handles.pos_a,handles.pos_c];
            position = customWait(ROI);
            position = [0,round(position(2)),round(position(1))];
            movement = center - position;
            handles.IPCMRA = imtranslate(handles.IPCMRA,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.FFE = imtranslate(handles.FFE,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordx = imtranslate(handles.coordx,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordy = imtranslate(handles.coordy,movement,'FillValues',0,'OutputView','same', 'method','nearest');
            handles.coordz = imtranslate(handles.coordz,movement,'FillValues',0,'OutputView','same', 'method','nearest');

            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');
            
        end
        
    case 3

        if get(handles.popupmenu7,'Value')==1
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
            ROI = images.roi.Line(handles.axes3,'Position',[handles.pos_c, handles.pos_a; handles.pos_c, round((handles.pos_a + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');

            vect_ini = [handles.pos_c, handles.pos_a] - [handles.pos_c, round((handles.pos_a + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_c);

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[1 0 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[1 0 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[1 0 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[1 0 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[1 0 0],"linear","crop");

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

            else

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
        
        elseif get(handles.popupmenu7,'Value')==2
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
            ROI = images.roi.Line(handles.axes3,'Position',[handles.pos_c, handles.pos_a; handles.pos_c, round((handles.pos_a + size(handles.IPCMRA,2))/2)]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');

            vect_ini = [handles.pos_c, handles.pos_a] - [handles.pos_c, round((handles.pos_a + size(handles.IPCMRA,2))/2)];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,1) - handles.pos_c);

            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[1 0 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[1 0 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[1 0 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[1 0 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[1 0 0],"linear","crop");

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

            else

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
            
        end

    case 4

        if get(handles.popupmenu7,'Value')==1
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
            ROI = images.roi.Line(handles.axes3,'Position',[handles.pos_c, handles.pos_a; round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_a]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');


            vect_ini = [handles.pos_c, handles.pos_a] - [round((handles.pos_c + size(handles.IPCMRA,1))/2), handles.pos_a];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_c)*-1;


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[1 0 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[1 0 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[1 0 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[1 0 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[1 0 0],"linear","crop");

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

            else

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
        
        elseif get(handles.popupmenu7,'Value')==2
            
            axes(handles.axes3);
            axis off
            imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
            hold on
            images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
            ROI = images.roi.Line(handles.axes3,'Position',[handles.pos_c, handles.pos_a; round((handles.pos_c + size(handles.IPCMRA,3))/2), handles.pos_a]);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axtoolbar(handles.axes3,'Visible','off');


            vect_ini = [handles.pos_c, handles.pos_a] - [round((handles.pos_c + size(handles.IPCMRA,1))/2), handles.pos_a];
            position = customWait(ROI);
            vect_final = [position(1,:)] - [position(2,:)];
            angulo = acosd((dot(vect_ini,vect_final)/(norm(vect_ini)*norm(vect_final))))*sign(position(2,2) - handles.pos_c)*-1;


            if abs(angulo)<=90

                handles.IPCMRA = imrotate3(handles.IPCMRA,angulo,[1 0 0],"linear","crop");
                handles.FFE = imrotate3(handles.FFE,angulo,[1 0 0],"linear","crop");
                handles.coordx = imrotate3(handles.coordx,angulo,[1 0 0],"linear","crop");
                handles.coordy = imrotate3(handles.coordy,angulo,[1 0 0],"linear","crop");
                handles.coordz = imrotate3(handles.coordz,angulo,[1 0 0],"linear","crop");

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

            else

                axes(handles.axes3);
                axis off
                imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                hold on
                images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
                hold off
                axis image
                colormap gray
                axis off
                daspect([1 1 1])
                axtoolbar(handles.axes3,'Visible','off');

                waitfor(msgbox("The angle need to be less than 90 degree","Error","error"));

            end
        end

end

if get(handles.popupmenu5,'Value')==1
    
    axes(handles.axes1);
    axis off
    imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
    hold on
    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes1,'Visible','off');
    
elseif get(handles.popupmenu5,'Value')==2
    
    axes(handles.axes1);
    axis off
    imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
    hold on
    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes1,'Visible','off');
    
end

if get(handles.popupmenu6,'Value')==1
    
    axes(handles.axes2);
    axis off
    imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
    hold on
    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes2,'Visible','off');
    
elseif get(handles.popupmenu6,'Value')==2

    axes(handles.axes2);
    axis off
    imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
    hold on
    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axtoolbar(handles.axes2,'Visible','off');
    
end

    list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
    set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
    set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);


    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        point = [];
        point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
        point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
        point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

        point = point*handles.min_voxel;
        point(1) = point(1)/handles.voxel_MR(1);
        point(2) = point(2)/handles.voxel_MR(2);
        point(3) = point(3)/handles.voxel_MR(3);
        
        p_norm = [];
        p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b+4,handles.pos_c);
        p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b+4,handles.pos_c);
        p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b+4,handles.pos_c);

        p_norm = p_norm*handles.min_voxel;
        p_norm(1) = p_norm(1)/handles.voxel_MR(1);
        p_norm(2) = p_norm(2)/handles.voxel_MR(2);
        p_norm(3) = p_norm(3)/handles.voxel_MR(3);

%         disp(point)
%         disp(p_norm)
        normal = p_norm - point;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % mostramos los cortes a la salida
    
        [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);
        [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);
        MR_FFE_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_FH_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_AP_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        MR_PCA_RL_2D = zeros(size(IPCMRA_2D,1),size(IPCMRA_2D,2),size(handles.MR_FFE_FH_ORI_COOR,4));
        
        for dt = 1:size(handles.MR_FFE_FH_ORI_COOR,4)
            MR_FFE_FH_2D(:,:,dt) = obliqueslice(handles.MR_FFE_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_FH_2D(:,:,dt) = obliqueslice(handles.MR_PCA_FH_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_AP_2D(:,:,dt) = obliqueslice(handles.MR_PCA_AP_ORI_COOR(:,:,:,dt),point,normal);
            MR_PCA_RL_2D(:,:,dt) = obliqueslice(handles.MR_PCA_RL_ORI_COOR(:,:,:,dt),point,normal);
        end

        if get(handles.popupmenu8,'Value')==1
            
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        elseif get(handles.popupmenu8,'Value')==2
           
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            plot3(point(1),point(2),point(3),'*r')
            plot3(p_norm(1),p_norm(2),p_norm(3),'*g')
            quiver3(point(1),point(2),point(3),normal(1),normal(2),normal(3),10,'c')
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');
            
        end

%         figure, 
%         subplot 221, imshow(MR_FFE_FH_2D(:,:,3),[])
%         subplot 222, imshow(MR_PCA_FH_2D(:,:,3),[])
%         subplot 223, imshow(MR_PCA_AP_2D(:,:,3),[])
%         subplot 224, imshow(MR_PCA_RL_2D(:,:,3),[])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyectamos la velocidad a travez del plano
%         disp(min(MR_FFE_FH_2D(:)))
%         disp(min(MR_PCA_FH_2D(:)))
%         disp(min(MR_PCA_AP_2D(:)))
%         disp(min(MR_PCA_RL_2D(:)))
%         
%         disp(max(MR_FFE_FH_2D(:)))
%         disp(max(MR_PCA_FH_2D(:)))
%         disp(max(MR_PCA_AP_2D(:)))
%         disp(max(MR_PCA_RL_2D(:)))
        
        % normal vector
        nomal_vector = normal/norm(normal);
        
        % projection of the vector
        magnitud = MR_FFE_FH_2D;
        velocity = MR_PCA_FH_2D*nomal_vector(1) + MR_PCA_AP_2D*nomal_vector(2) + MR_PCA_RL_2D*nomal_vector(3);
        
%         figure, 
%         subplot 121, imshow(magnitud(:,:,3),[])
%         subplot 122, imshow(velocity(:,:,3),[])
        
        %%% coordinates        
        pp = zeros(size(cx,1),size(cx,2),3);
        pp(:,:,1) = cx*handles.voxel_MR(1);
        pp(:,:,2) = cy*handles.voxel_MR(2);
        pp(:,:,3) = cz*handles.voxel_MR(3);
        
        x_coor = zeros(size(cx));
        y_coor = zeros(size(cy));
        
        x_coor(:,2:end) = sqrt(sum((pp(:,2:end,:) - pp(:,1:end-1,:)).^2,3));
        y_coor(2:end,:) = sqrt(sum((pp(2:end,:,:) - pp(1:end-1,:,:)).^2,3));
        
        pixel_size = [x_coor(2,2),y_coor(2,2)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % llamamos a la funcion para escribir el archivo
        structure = create_structure(magnitud, velocity, handles.heart_rate, handles.VENC, pixel_size);
        
        setstruct = structure.setstruct;
        preview = structure.preview;
        info = structure.info;
        im = structure.im;
        
        directory = uigetdir(pwd, 'Select Directory');
        save([directory,'\structure.mat'],'setstruct','preview','info','im')
      
        f = msgbox("The data has been successfully saved","Success");
        pause(2)
        close(f)

handles.output = hObject;
guidata(hObject, handles);      


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7

switch get(handles.popupmenu7,'Value')

    case 1

        axes(handles.axes3);
        axis off
        imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
        hold on
        images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes3,'Visible','off');
        
    case 2

        axes(handles.axes3);
        axis off
        imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
        hold on
        images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'all');
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axtoolbar(handles.axes3,'Visible','off');
        
end

handles.output = hObject;
guidata(hObject, handles);  

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

if get(handles.popupmenu8,'Value')==1
    
    switch get(handles.popupmenu4,'Value')

        case 1

            axes(handles.axes4)
            plot(0,0)
            surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');



        case 2

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);


            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c+4);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


        case 3

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a+4,handles.pos_b,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        case 4

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b+4,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


    end

elseif get(handles.popupmenu8,'Value')==2
    
    switch get(handles.popupmenu4,'Value')

        case 1

            axes(handles.axes4)
            plot(0,0)
            surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');



        case 2

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);


            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c+4);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);
    
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


        case 3

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a+4,handles.pos_b,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        case 4

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b+4,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


    end
    
end
handles.output = hObject;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on selection change in popupmenu8.
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8
if get(handles.popupmenu8,'Value')==1
    
    switch get(handles.popupmenu4,'Value')

        case 1

            axes(handles.axes4)
            plot(0,0)
            surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');



        case 2

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);


            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c+4);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


        case 3

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a+4,handles.pos_b,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        case 4

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b+4,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [IPCMRA_2D,cx,cy,cz] = obliqueslice(handles.IPCMRA_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,handles.IPCMRA_1,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,handles.IPCMRA_2,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,handles.IPCMRA_3,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,IPCMRA_2D,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


    end

elseif get(handles.popupmenu8,'Value')==2
    
    switch get(handles.popupmenu4,'Value')

        case 1

            axes(handles.axes4)
            plot(0,0)
            surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');



        case 2

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);


            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c+4);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c+4);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);
    
            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


        case 3

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a+4,handles.pos_b,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a+4,handles.pos_b,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');

        case 4

            point = [];
            point(1) = handles.coordx(handles.pos_a,handles.pos_b,handles.pos_c);
            point(2) = handles.coordy(handles.pos_a,handles.pos_b,handles.pos_c);
            point(3) = handles.coordz(handles.pos_a,handles.pos_b,handles.pos_c);

            point = point*handles.min_voxel;
            point(1) = point(1)/handles.voxel_MR(1);
            point(2) = point(2)/handles.voxel_MR(2);
            point(3) = point(3)/handles.voxel_MR(3);

            p_norm = [];
            p_norm(1) = handles.coordx(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(2) = handles.coordy(handles.pos_a,handles.pos_b+4,handles.pos_c);
            p_norm(3) = handles.coordz(handles.pos_a,handles.pos_b+4,handles.pos_c);

            p_norm = p_norm*handles.min_voxel;
            p_norm(1) = p_norm(1)/handles.voxel_MR(1);
            p_norm(2) = p_norm(2)/handles.voxel_MR(2);
            p_norm(3) = p_norm(3)/handles.voxel_MR(3);

            normal = p_norm - point;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % mostramos los cortes a la salida

            [FFE_2D,cx,cy,cz] = obliqueslice(handles.FFE_ORI_COOR,point,normal);

            axes(handles.axes4)
            plot(0,0)
            s1 = surface(handles.corte_sagital_x,handles.corte_sagital_y,handles.corte_sagital_z,(handles.FFE_1/max(handles.FFE_1(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold on
            s2 = surface(handles.corte_coronal_x,handles.corte_coronal_y,handles.corte_coronal_z,(handles.FFE_2/max(handles.FFE_2(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            s3 = surface(handles.corte_axial_x,handles.corte_axial_y,handles.corte_axial_z,(handles.FFE_3/max(handles.FFE_3(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            alpha(s1,.5)
            alpha(s2,.5)
            alpha(s3,.5)
            surf(cx,cy,cz,(FFE_2D/max(FFE_2D(:)))*255,'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            hold off
            axis image
            colormap gray
            grid on
            axis on
            xlabel('F-H')
            ylabel('A-P')
            zlabel('R-L')
            view(35,45)
            daspect([1 1 1])
            axtoolbar(handles.axes4,{'rotate', 'restoreview'},'Visible','on');


    end
    
end
handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    handles.id_while = 1;
    handles.output = hObject;
    guidata(hObject, handles);
    uiresume(hObject);
else
    delete(hObject);
end


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    answer = questdlg('which image do you want to modify?', ...
	'Image selection', ...
	'Magnitude','MR-Angiography','MR-Angiography');
    % Handle response
    switch answer
        case 'Magnitude'
            handles.id_image = 1;

            
            handles.input1 = (handles.FFE/max(handles.FFE(:)))*255;
            handles.input2 = (handles.FFE/max(handles.FFE(:)))*255;
            id_while = 0;
            while(id_while == 0)
                GUIDE_CONTRAST(handles.input1, handles.input2);
                handles.input1 = getappdata(0,'OUT');
                id_while = getappdata(0,'closed_loop');
                handles.FFE = handles.input1;
                handles.FFE = (handles.FFE/max(handles.FFE(:)))*255;

                if get(handles.popupmenu5,'Value')==1
        
                    axes(handles.axes1);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                    hold on
                    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes1,'Visible','off');
                    
                elseif get(handles.popupmenu5,'Value')==2
                    
                    axes(handles.axes1);
                    axis off
                    imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                    hold on
                    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes1,'Visible','off');
                    
                end
                
                if get(handles.popupmenu6,'Value')==1
                    
                    axes(handles.axes2);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                    hold on
                    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes2,'Visible','off');
                    
                elseif get(handles.popupmenu6,'Value')==2
                
                    axes(handles.axes2);
                    axis off
                    imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                    hold on
                    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes2,'Visible','off');
                    
                end
        
                if get(handles.popupmenu7,'Value')==1
            
                    axes(handles.axes3);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                    hold on
                    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes3,'Visible','off');
                    
                elseif get(handles.popupmenu7,'Value')==2
                    
                    axes(handles.axes3);
                    axis off
                    imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                    hold on
                    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes3,'Visible','off');
        
                end
                
                list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
                set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
                set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
                set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);

            end

        case 'MR-Angiography'
            handles.id_image = 2;


            handles.input1 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
            handles.input2 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
            id_while = 0;
            while(id_while == 0)
                GUIDE_CONTRAST(handles.input1, handles.input2);
                handles.input1 = getappdata(0,'OUT');
                id_while = getappdata(0,'closed_loop');
                handles.IPCMRA = handles.input1;
                handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;

                if get(handles.popupmenu5,'Value')==1
        
                    axes(handles.axes1);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(:,:,handles.pos_c)))
                    hold on
                    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes1,'Visible','off');
                    
                elseif get(handles.popupmenu5,'Value')==2
                    
                    axes(handles.axes1);
                    axis off
                    imagesc(squeeze(handles.FFE(:,:,handles.pos_c)))
                    hold on
                    images.roi.Crosshair(handles.axes1,'Position',[handles.pos_b,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes1,'Visible','off');
                    
                end
                
                if get(handles.popupmenu6,'Value')==1
                    
                    axes(handles.axes2);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(handles.pos_a,:,:)))
                    hold on
                    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes2,'Visible','off');
                    
                elseif get(handles.popupmenu6,'Value')==2
                
                    axes(handles.axes2);
                    axis off
                    imagesc(squeeze(handles.FFE(handles.pos_a,:,:)))
                    hold on
                    images.roi.Crosshair(handles.axes2,'Position',[handles.pos_c,handles.pos_b],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes2,'Visible','off');
                    
                end
        
                if get(handles.popupmenu7,'Value')==1
            
                    axes(handles.axes3);
                    axis off
                    imagesc(squeeze(handles.IPCMRA(:,handles.pos_b,:)))    
                    hold on
                    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes3,'Visible','off');
                    
                elseif get(handles.popupmenu7,'Value')==2
                    
                    axes(handles.axes3);
                    axis off
                    imagesc(squeeze(handles.FFE(:,handles.pos_b,:)))    
                    hold on
                    images.roi.Crosshair(handles.axes3,'Position',[handles.pos_c,handles.pos_a],'LineWidth',1,'Color','y', 'InteractionsAllowed', 'none');
                    hold off
                    axis image
                    colormap gray
                    axis off
                    daspect([1 1 1])
                    axtoolbar(handles.axes3,'Visible','off');
        
                end
                
                list_string1 = {'...','move crosshair','angle vertical','angle horizontal'};
                set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
                set(handles.popupmenu2,'visible','on','String',list_string1,'value',1);
                set(handles.popupmenu3,'visible','on','String',list_string1,'value',1);
            end
    end

%     handles.input1 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
%     handles.input2 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
%     id_while = 0;
%     while(id_while == 0)
%         GUIDE_CONTRAST(handles.input1, handles.input2);
%         handles.input1 = getappdata(0,'OUT');
%         id_while = getappdata(0,'closed_loop');
%         handles.IPCMRA = handles.input1;
%         handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;

%         % Julio Sotelo 23-11-2018
%         if handles.id_ang == 1% Julio Sotelo 23-11-2018
%             handles.ANG = handles.IPCMRA;% Julio Sotelo 23-11-2018
%         elseif handles.id_mag == 1% Julio Sotelo 23-11-2018
%             handles.MAG = handles.IPCMRA;% Julio Sotelo 23-11-2018
%         end% Julio Sotelo 23-11-2018
% 
%         if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
%         if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
%         if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
%         axes(handles.axes1);
%         imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
%         hold on
%         plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
%         plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
%         if handles.id_seg == 1 && sum(handles.L(:))~=0
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
%         end
%         if handles.id_resizing == 1
%             rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
%         end
%         hold off
%         axis image
%         colormap gray
%         axis off
%         daspect([1 1 1])
%         axes(handles.axes2);
%         imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
%         hold on
%         plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
%         plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
%         if handles.id_seg == 1 && sum(handles.L(:))~=0
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
%         end
%         if handles.id_resizing == 1
%             rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
%         end
%         hold off
%         axis image
%         colormap gray
%         axis off
%         daspect([1 1 1])
%         axes(handles.axes3);
%         imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
%         hold on
%         plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
%         plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
%         if handles.id_seg == 1 && sum(handles.L(:))~=0
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
%         end
%         if handles.id_resizing == 1
%             rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
%         end
%         hold off
%         axis image
%         colormap gray
%         axis off
%         daspect([1 1 1])

        



%     end
%     
handles.output = hObject;
guidata(hObject, handles);
