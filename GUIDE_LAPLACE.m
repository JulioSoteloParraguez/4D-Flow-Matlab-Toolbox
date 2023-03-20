function varargout = GUIDE_LAPLACE(varargin)
% GUIDE_LAPLACE MATLAB code for GUIDE_LAPLACE.fig
%      GUIDE_LAPLACE, by itself, creates a new GUIDE_LAPLACE or raises the existing
%      singleton*.
%
%      H = returns the handle to a new GUIDE_LAPLACE or the handle to
%      the existing singleton*.
%
%      GUIDE_LAPLACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_LAPLACE.M with the given input arguments.
%
%      GUIDE_LAPLACE('Property','Value',...) creates a new GUIDE_LAPLACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_LAPLACE_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_LAPLACE_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_LAPLACE

% Last Modified by GUIDE v2.5 17-Jan-2023 17:29:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_LAPLACE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_LAPLACE_OutputFcn, ...
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


% --- Executes just before GUIDE_LAPLACE is made visible.
function GUIDE_LAPLACE_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE_LAPLACE (see VARARGIN)

% Choose default command line output for GUIDE_LAPLACE
handles.output = hObject;

    path(path,'iso2mesh/')
    path(path,'IO_CODES_MATLAB_CODES/') % cambiar
    
    handles.rot = 1;
    
    handles.id_vwerp        = varargin{1}.id_vwerp;
    
    if handles.id_vwerp == 1

        addpath(genpath('vWERP/'));
        handles.SEG             = varargin{1}.SEG;
        handles.SEG_for_vwerp   = varargin{1}.SEG_for_vwerp;
        handles.IPCMRA          = varargin{1}.IPCMRA;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.L               = varargin{1}.L;
        handles.Lrgb            = varargin{1}.Lrgb;
        handles.Lrgb_vel        = varargin{1}.Lrgb_vel;
        handles.NUM             = varargin{1}.NUM;
        handles.xd              = varargin{1}.xd;
        handles.yd              = varargin{1}.yd;
        handles.zd              = varargin{1}.zd;
        handles.a               = varargin{1}.a;
        handles.b               = varargin{1}.b;
        handles.c               = varargin{1}.c;
        handles.d               = varargin{1}.d;
        handles.slider_axes1    = varargin{1}.slider_id_axes1;
        handles.MR_FFE_FH       = varargin{1}.MR_FFE_FH;
        handles.MAG             = mean(handles.MR_FFE_FH,4); 

        handles.MR_PCA_AP       = varargin{1}.MR_PCA_AP;
        handles.MR_PCA_FH       = varargin{1}.MR_PCA_FH;
        handles.MR_PCA_RL       = varargin{1}.MR_PCA_RL;
        
        set(handles.pushbutton11,'Visible','on');
        set(handles.pushbutton10,'Visible','on');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % debo corregir el valor de slider axes, para que abra siempre el mismo
        % valor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % read mesh
        handles.faces       = varargin{1}.faces;
        handles.nodes       = varargin{1}.nodes;
        handles.elem        = varargin{1}.elem;
        handles.veset       = varargin{1}.veset;
        handles.mags_vel    = varargin{1}.mags_vel;

        handles.heart_rate  = varargin{1}.heart_rate;
        handles.id_view     = varargin{1}.id_view;
        handles.id_ipcmra   = varargin{1}.id_ipcmra;
        handles.id_mag      = varargin{1}.id_mag;

        handles.peak_flow   = varargin{1}.peak_flow;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_ipcmra == 1
            handles.id_image = 1;
        elseif handles.id_mag == 1
            handles.id_image = 2;
        else
            handles.id_image = 0;
        end
        
        handles.nodes_id            = (1:size(handles.nodes,1))'; 
        handles.elem_id             = (1:size(handles.elem,1))';
        handles.faces_id            = (1:size(handles.faces,1))';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0;


        % read id

        handles.id_mesh_inlet = varargin{1}.id_mesh_inlet;
        handles.id_mesh_outlet = varargin{1}.id_mesh_outlet;

        handles.id_inlet = varargin{1}.id_inlet;
        handles.id_outlet = varargin{1}.id_outlet;

        handles.inlet_exec = varargin{1}.inlet_exec;
        handles.outlet_exec = varargin{1}.outlet_exec;

        % id loop
        handles.id_while = 0;

        % handles.list_n = varargin{1}.list_n;
        
        
        
        if handles.inlet_exec == 0 && handles.outlet_exec == 0

            set(handles.text3,'Visible','on');
            set(handles.text4,'Visible','on');
            set(handles.text5,'Visible','on');
            set(handles.text6,'Visible','on');
            set(handles.radiobutton1,'Visible','on','Value',1);
            set(handles.radiobutton2,'Visible','on','Value',0);
            set(handles.radiobutton3,'Visible','on','Value',1);
            set(handles.radiobutton4,'Visible','on','Value',0);
            set(handles.edit2,'Visible','on');
            set(handles.edit3,'Visible','on');

            handles.werp = 2;
            handles.volInt = true;
            handles.resamp_factor = 1;
            handles.percentage_noise = 0.0;
            set(handles.pushbutton3,'String','vWERP','FontUnits','Normalized','FontSize',0.24)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % adjust the image pointer

        if handles.id_mesh_inlet == 1

            
            set(handles.popupmenu1,'Value',1)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.outlet_exec == 1
                handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
                set(handles.popupmenu2,'String',handles.list_string,'Value',6);

            else

                handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet'};
                set(handles.popupmenu2,'String',handles.list_string,'Value',6);

            end
            
            if get(handles.popupmenu2,'Value')==6
                set(handles.pushbutton5,'Visible','off');
                set(handles.pushbutton6,'Visible','on');
            end

            axes(handles.axes1);
            plot(0.0)
    %         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
            center_in = mean(handles.nodes(handles.id_inlet,1:3));
            quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')


            

            handles.id_ipcmra = 0;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 1;
            handles.id_mesh_outlet = 0;
%             handles.id_mesh_laplace = 0;
%             handles.id_centerline = 0;
%             handles.id_diameter = 0;
%             handles.id_radius = 0;
%             handles.id_axial_unit_vectors = 0;
%             handles.id_circumferential_unit_vectors = 0;
%             handles.id_wss_a = 0;
%             handles.id_wss_c = 0;
%             handles.id_axial_angle = 0;
%             handles.id_forward_vel = 0;
%             handles.id_backward_vel = 0;
%             handles.id_regurgitant_flow = 0;
%             handles.id_centerline_flow = 0;
%             handles.id_eccentricity = 0;

%             handles.id_curvature = 0; % new data Julio Sotelo
%             handles.id_ellipticity = 0; % new data Julio Sotelo
%             handles.id_flattening = 0; % new data Julio Sotelo
%     %         handles.id_circulation = 0; % new data Julio Sotelo
%             handles.id_length_vessel = 0; % new data Julio Sotelo
%             handles.id_forward_vortex = 0; % new data Julio Sotelo
%             handles.id_area = 0; % new data Julio Sotelo
%             handles.id_axial_circulation = 0; % new data Julio Sotelo
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.inlet_exec == 1 && handles.outlet_exec == 1

                set(handles.pushbutton3,'Visible','on');

            end

        elseif handles.id_mesh_outlet == 1

            set(handles.popupmenu1,'Value',1)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',7);

            if get(handles.popupmenu2,'Value')==7
                set(handles.pushbutton5,'Visible','on');
                set(handles.pushbutton6,'Visible','off');
            end
                     
            axes(handles.axes1);
            plot(0.0)
    %         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            center_out = mean(handles.nodes(handles.id_outlet,1:3));
            quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
            plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')

            
            
            handles.id_ipcmra = 0;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 1;
%             handles.id_mesh_laplace = 0;
%             handles.id_centerline = 0;
%             handles.id_diameter = 0;
%             handles.id_radius = 0;
%             handles.id_axial_unit_vectors = 0;
%             handles.id_circumferential_unit_vectors = 0;
%             handles.id_wss_a = 0;
%             handles.id_wss_c = 0;
%             handles.id_axial_angle = 0;
%             handles.id_forward_vel = 0;
%             handles.id_backward_vel = 0;
%             handles.id_regurgitant_flow = 0;
%             handles.id_centerline_flow = 0;
%             handles.id_eccentricity = 0;
% 
%             handles.id_curvature = 0; % new data Julio Sotelo
%             handles.id_ellipticity = 0; % new data Julio Sotelo
%             handles.id_flattening = 0; % new data Julio Sotelo
%     %         handles.id_circulation = 0; % new data Julio Sotelo
%             handles.id_length_vessel = 0; % new data Julio Sotelo
%             handles.id_forward_vortex = 0; % new data Julio Sotelo
%             handles.id_area = 0; % new data Julio Sotelo
%             handles.id_axial_circulation = 0; % new data Julio Sotelo

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.inlet_exec == 1 && handles.outlet_exec == 1

                set(handles.pushbutton3,'Visible','on');
                set(handles.pushbutton1,'Visible','off');
                set(handles.pushbutton2,'Visible','off');
                set(handles.popupmenu1,'Value',1);
                set(handles.popupmenu3,'Value',1);

            end

        else

            if handles.id_ipcmra == 1

                handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
                set(handles.popupmenu1,'String',handles.list_string,'Value',2);

                if handles.id_view == 1

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                    himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %                 himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.c;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 2

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.b;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 3

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.a;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                end

            elseif handles.id_mag == 1

                handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
                set(handles.popupmenu1,'String',handles.list_string,'Value',3);

                if handles.id_view == 1

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                    himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %                 himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.c;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 2

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.b;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 3

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.a;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet'};
            set(handles.popupmenu2,'String',handles.list_string);

            set(handles.pushbutton1,'Visible','on');
            if handles.inlet_exec == 1

                set(handles.pushbutton2,'Visible','on');

            end

            set(handles.pushbutton3,'Visible','off');

            handles.id_ipcmra = 1;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
%             handles.id_mesh_laplace = 0;
%             handles.id_centerline = 0;
%             handles.id_diameter = 0;
%             handles.id_radius = 0;
%             handles.id_axial_unit_vectors = 0;
%             handles.id_circumferential_unit_vectors = 0;
%             handles.id_wss_a = 0;
%             handles.id_wss_c = 0;
%             handles.id_axial_angle = 0;
%             handles.id_forward_vel = 0;
%             handles.id_backward_vel = 0;
%             handles.id_regurgitant_flow = 0;
%             handles.id_centerline_flow = 0;
%             handles.id_eccentricity = 0;
% 
%             handles.id_curvature = 0; % new data Julio Sotelo
%             handles.id_ellipticity = 0; % new data Julio Sotelo
%             handles.id_flattening = 0; % new data Julio Sotelo
%     %         handles.id_circulation = 0; % new data Julio Sotelo
%             handles.id_length_vessel = 0; % new data Julio Sotelo
%             handles.id_forward_vortex = 0; % new data Julio Sotelo
%             handles.id_area = 0; % new data Julio Sotelo
%             handles.id_axial_circulation = 0; % new data Julio Sotelo

            handles.inlet_exec = 0;
            handles.outlet_exec = 0;

            handles.id_section_loaded = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             answer = questdlg({'Do you want load the mesh inlet-outlet ids?...'},'Warning','Yes','No','No');
%             switch answer
% 
%                 case 'Yes'
% 
%                     folder_name = uigetdir([],'Load Folder...');
% 
%                     mat_inlet   = readtable([folder_name,'/Inlet.csv']);
%                     mat_outlet  = readtable([folder_name,'/Outlet.csv']);
% 
%                     handles.id_inlet    = mat_inlet.IDNodes___;
%                     handles.id_outlet   = mat_outlet.IDNodes___;
% 
%                     handles.inlet_exec = 1;
%                     handles.outlet_exec = 1;
% 
%                     handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
%                     set(handles.popupmenu2,'String',handles.list_string);
% 
%                     set(handles.pushbutton2,'Visible','on');
%                     set(handles.pushbutton3,'Visible','on');
% 
%                     handles.id_section_loaded = 1;
% 
% 
%                 case 'No'
% 
%             end

%             handles.Laplace = [];
%             handles.centerline = [];
%             handles.centerline_lapid = [];
%             handles.radius = [];
%             handles.diameter = [];
%             handles.axial_unit_vectors = [];
%             handles.circumferential_unit_vectors = [];
%             handles.WSS_A = [];
%             handles.WSS_C = [];
%             handles.mag_WSS_A = [];
%             handles.mag_WSS_C = [];
%             handles.angle_axial_direction = [];
%             handles.forward_velocity = [];
%             handles.backward_velocity = [];
%             handles.mag_forward_velocity = [];
%             handles.mag_backward_velocity = [];
%             handles.regurgitant_flow = [];
%             handles.centerline_flow = [];
%             handles.eccentricity = [];
% 
%             handles.curvature = []; % new data Julio Sotelo
%             handles.ellipticity = []; % new data Julio Sotelo
%             handles.flattening = []; % new data Julio Sotelo
%             handles.circulation = []; % new data Julio Sotelo
%             handles.length_vessel = []; % new data Julio Sotelo
%             handles.forward_vortex = []; % new data Julio Sotelo
%             handles.area = []; % new data Julio Sotelo
%             handles.axial_circulation = []; % new data Julio Sotelo


        end
        
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    else
        
        handles.SEG             = varargin{1}.SEG;
        handles.IPCMRA          = varargin{1}.IPCMRA;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.L               = varargin{1}.L;
        handles.Lrgb            = varargin{1}.Lrgb;
        handles.Lrgb_vel        = varargin{1}.Lrgb_vel;
        handles.NUM             = varargin{1}.NUM;
        handles.xd              = varargin{1}.xd;
        handles.yd              = varargin{1}.yd;
        handles.zd              = varargin{1}.zd;
        handles.a               = varargin{1}.a;
        handles.b               = varargin{1}.b;
        handles.c               = varargin{1}.c;
        handles.d               = varargin{1}.d;
        handles.slider_axes1    = varargin{1}.slider_id_axes1;
        handles.MR_FFE_FH       = varargin{1}.MR_FFE_FH;
        handles.MAG             = mean(handles.MR_FFE_FH,4); 

        handles.SEG             = varargin{1}.SEG;
        handles.MR_PCA_AP       = varargin{1}.MR_PCA_AP;
        handles.MR_PCA_FH       = varargin{1}.MR_PCA_FH;
        handles.MR_PCA_RL       = varargin{1}.MR_PCA_RL;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.vorticity       = varargin{1}.vorticity;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % debo corregir el valor de slider axes, para que abra siempre el mismo
        % valor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % read mesh
        handles.faces       = varargin{1}.faces;
        handles.nodes       = varargin{1}.nodes;
        handles.elem        = varargin{1}.elem;
        handles.veset       = varargin{1}.veset;
        handles.mags_vel    = varargin{1}.mags_vel;
        handles.WSS         = varargin{1}.WSS;
        handles.mags_wss    = varargin{1}.mags_wss;
        handles.peak_flow   = varargin{1}.peak_flow;
        handles.peak_flow_ori = varargin{1}.peak_flow_ori;
    %     handles.mesh_id     = [1:size(handles.nodes,1)]';
        handles.heart_rate  = varargin{1}.heart_rate;
        handles.id_view     = varargin{1}.id_view;
        handles.id_ipcmra   = varargin{1}.id_ipcmra;
        handles.id_mag      = varargin{1}.id_mag;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_ipcmra == 1
            handles.id_image = 1;
        elseif handles.id_mag == 1
            handles.id_image = 2;
        else
            handles.id_image = 0;
        end
        handles.nodes_id            = (1:size(handles.nodes,1))'; 
        handles.elem_id             = (1:size(handles.elem,1))';
        handles.faces_id            = (1:size(handles.faces,1))';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0;


        % read id

        handles.id_mesh_inlet = varargin{1}.id_mesh_inlet;
        handles.id_mesh_outlet = varargin{1}.id_mesh_outlet;

        handles.id_inlet = varargin{1}.id_inlet;
        handles.id_outlet = varargin{1}.id_outlet;

        handles.inlet_exec = varargin{1}.inlet_exec;
        handles.outlet_exec = varargin{1}.outlet_exec;

        % id loop
        handles.id_while = 0;

        handles.list_n = varargin{1}.list_n;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % adjust the image pointer

        if handles.id_mesh_inlet == 1

            set(handles.popupmenu1,'Value',1)

            axes(handles.axes1);
            plot(0.0)
    %         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.outlet_exec == 1

                handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
                set(handles.popupmenu2,'String',handles.list_string,'Value',4);

            else

                handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet'};
                set(handles.popupmenu2,'String',handles.list_string,'Value',4);

            end

            handles.id_ipcmra = 0;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 1;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;

            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_area = 0; % new data Julio Sotelo
            handles.id_axial_circulation = 0; % new data Julio Sotelo

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.inlet_exec == 1 && handles.outlet_exec == 1

                set(handles.pushbutton3,'Visible','on');

            end

        elseif handles.id_mesh_outlet == 1

            set(handles.popupmenu1,'Value',1)

            axes(handles.axes1);
            plot(0.0)
    %         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity Vectors','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',5);

            handles.id_ipcmra = 0;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 1;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;

            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_area = 0; % new data Julio Sotelo
            handles.id_axial_circulation = 0; % new data Julio Sotelo

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if handles.inlet_exec == 1 && handles.outlet_exec == 1

                set(handles.pushbutton3,'Visible','on');
                set(handles.pushbutton1,'Visible','off');
                set(handles.pushbutton2,'Visible','off');
                set(handles.popupmenu1,'Value',1);
                set(handles.popupmenu3,'Value',1);

            end

        else

            if handles.id_ipcmra == 1

                handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
                set(handles.popupmenu1,'String',handles.list_string,'Value',2);

                if handles.id_view == 1

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                    himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %                 himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.c;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 2

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.b;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 3

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.a;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                end

            elseif handles.id_mag == 1

                handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
                set(handles.popupmenu1,'String',handles.list_string,'Value',3);

                if handles.id_view == 1

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                    himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %                 himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.c;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 2

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.b;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                elseif handles.id_view == 3

                    handles.list_string = {'Sagital View','Axial View','Coronal View'};
                    set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    axes(handles.axes1);
                    plot(0.0)
                    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
                    hold on
                    Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
                    himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                    cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                    set(himage, 'AlphaData', cdata);
    %                 Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %                 himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %                 cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %                 cdata = double(cdata)*0.5;
    %                 set(himage, 'AlphaData', cdata);
                    hold off
                    axis image
                    set(handles.axes1,'xticklabel',[],'yticklabel',[])
                    axis off
                    colormap gray

                    % slider adjustment
                    slider_step(1) = 1/handles.a;
                    slider_step(2) = 0.1;
                    set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
                    set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity'};
            set(handles.popupmenu2,'String',handles.list_string);

            set(handles.pushbutton1,'Visible','on');
            if handles.inlet_exec == 1

                set(handles.pushbutton2,'Visible','on');

            end

            set(handles.pushbutton3,'Visible','off');

            handles.id_ipcmra = 1;
            handles.id_mag = 0;
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;

            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_area = 0; % new data Julio Sotelo
            handles.id_axial_circulation = 0; % new data Julio Sotelo

            handles.inlet_exec = 0;
            handles.outlet_exec = 0;

            handles.id_section_loaded = 0;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             answer = questdlg({'Do you want load the mesh inlet-outlet ids?...'},'Warning','Yes','No','No');
%             switch answer
% 
%                 case 'Yes'
% 
%                     folder_name = uigetdir([],'Load Folder...');
% 
%                     mat_inlet   = readtable([folder_name,'/Inlet.csv']);
%                     mat_outlet  = readtable([folder_name,'/Outlet.csv']);
% 
%                     handles.id_inlet    = mat_inlet.IDNodes___;
%                     handles.id_outlet   = mat_outlet.IDNodes___;
% 
%                     handles.inlet_exec = 1;
%                     handles.outlet_exec = 1;
% 
%                     handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
%                     set(handles.popupmenu2,'String',handles.list_string);
% 
%                     set(handles.pushbutton2,'Visible','on');
%                     set(handles.pushbutton3,'Visible','on');
% 
%                     handles.id_section_loaded = 1;
% 
% 
%                 case 'No'
% 
%             end

            handles.Laplace = [];
            handles.centerline = [];
            handles.centerline_lapid = [];
            handles.radius = [];
            handles.diameter = [];
            handles.axial_unit_vectors = [];
            handles.circumferential_unit_vectors = [];
            handles.WSS_A = [];
            handles.WSS_C = [];
            handles.mag_WSS_A = [];
            handles.mag_WSS_C = [];
            handles.angle_axial_direction = [];
            handles.forward_velocity = [];
            handles.backward_velocity = [];
            handles.mag_forward_velocity = [];
            handles.mag_backward_velocity = [];
            handles.regurgitant_flow = [];
            handles.centerline_flow = [];
            handles.eccentricity = [];

            handles.curvature = []; % new data Julio Sotelo
            handles.ellipticity = []; % new data Julio Sotelo
            handles.flattening = []; % new data Julio Sotelo
            handles.circulation = []; % new data Julio Sotelo
            handles.length_vessel = []; % new data Julio Sotelo
            handles.forward_vortex = []; % new data Julio Sotelo
            handles.area = []; % new data Julio Sotelo
            handles.axial_circulation = []; % new data Julio Sotelo


        end
    end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUIDE_LAPLACE wait for user response (see UIRESUME)
uiwait(handles.GUIDE_LAPLACE);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_LAPLACE_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.id_vwerp == 1
    if handles.id_while == 0

        varargout{1} = handles.output;

    %     Laplace = handles.Laplace;
    %     setappdata(0,'Laplace',Laplace);
    %     centerline = handles.centerline;
    %     setappdata(0,'centerline',centerline);
    %     centerline_lapid = handles.centerline_lapid;
    %     setappdata(0,'centerline_lapid',centerline_lapid);
    %     radius = handles.radius;
    %     setappdata(0,'radius',radius);
    %     diameter = handles.diameter;
    %     setappdata(0,'diameter',diameter);
    %     axial_unit_vectors = handles.axial_unit_vectors;
    %     setappdata(0,'axial_unit_vectors',axial_unit_vectors);
    %     circumferential_unit_vectors = handles.circumferential_unit_vectors;
    %     setappdata(0,'circumferential_unit_vectors',circumferential_unit_vectors);
    %     WSS_A = handles.WSS_A;
    %     setappdata(0,'WSS_A',WSS_A);
    %     WSS_C = handles.WSS_C;
    %     setappdata(0,'WSS_C',WSS_C);
    %     mag_WSS_A = handles.mag_WSS_A;
    %     setappdata(0,'mag_WSS_A',mag_WSS_A);
    %     mag_WSS_C = handles.mag_WSS_C;
    %     setappdata(0,'mag_WSS_C',mag_WSS_C);
    %     angle_axial_direction = handles.angle_axial_direction;
    %     setappdata(0,'angle_axial_direction',angle_axial_direction);
    %     forward_velocity = handles.forward_velocity;
    %     setappdata(0,'forward_velocity',forward_velocity);
    %     backward_velocity = handles.backward_velocity;
    %     setappdata(0,'backward_velocity',backward_velocity);
    %     mag_forward_velocity = handles.mag_forward_velocity;
    %     setappdata(0,'mag_forward_velocity',mag_forward_velocity);
    %     mag_backward_velocity = handles.mag_backward_velocity;
    %     setappdata(0,'mag_backward_velocity',mag_backward_velocity);
    %     regurgitant_flow = handles.regurgitant_flow;
    %     setappdata(0,'regurgitant_flow',regurgitant_flow);
    %     centerline_flow = handles.centerline_flow;
    %     setappdata(0,'centerline_flow',centerline_flow);
    %     eccentricity = handles.eccentricity;
    %     setappdata(0,'eccentricity',eccentricity);
    %     
    %     curvature = handles.curvature;% new data Julio Sotelo
    %     setappdata(0,'curvature',curvature);% new data Julio Sotelo
    %     ellipticity = handles.ellipticity;% new data Julio Sotelo
    %     setappdata(0,'ellipticity',ellipticity);% new data Julio Sotelo
    %     length_vessel = handles.length_vessel;% new data Julio Sotelo
    %     setappdata(0,'length_vessel',length_vessel);% new data Julio Sotelo
    %     flattening = handles.flattening;% new data Julio Sotelo
    %     setappdata(0,'flattening',flattening);% new data Julio Sotelo
    % %     circulation = handles.circulation;% new data Julio Sotelo
    % %     setappdata(0,'circulation',circulation);% new data Julio Sotelo
    %     forward_vortex = handles.forward_vortex;% new data Julio Sotelo
    %     setappdata(0,'forward_vortex',forward_vortex);% new data Julio Sotelo
    %     area = handles.area;% new data Julio Sotelo
    %     setappdata(0,'area',area);% new data Julio Sotelo
    %     axial_circulation = handles.axial_circulation;% new data Julio Sotelo
    %     setappdata(0,'axial_circulation',axial_circulation);% new data Julio Sotelo

        id_while = handles.id_while;
        setappdata(0,'id_while',id_while);
        id_mesh_inlet = handles.id_mesh_inlet;
        setappdata(0,'id_mesh_inlet',id_mesh_inlet);
        id_mesh_outlet = handles.id_mesh_outlet;
        setappdata(0,'id_mesh_outlet',id_mesh_outlet);
        id_view = handles.id_view;
        setappdata(0,'id_view',id_view);
        id_inlet = handles.id_inlet;
        setappdata(0,'id_inlet',id_inlet);
        id_outlet = handles.id_outlet;
        setappdata(0,'id_outlet',id_outlet);
        inlet_exec = handles.inlet_exec;
        setappdata(0,'inlet_exec',inlet_exec);
        outlet_exec = handles.outlet_exec;
        setappdata(0,'outlet_exec',outlet_exec);
        id_ipcmra = handles.id_ipcmra;
        setappdata(0,'id_ipcmra',id_ipcmra);
        id_mag = handles.id_mag;
        setappdata(0,'id_mag',id_mag);
        slider_id_axes1 = handles.slider_axes1;
        setappdata(0,'slider_id_axes1',slider_id_axes1);
        SEG_for_vwerp = handles.SEG_for_vwerp;
        setappdata(0,'SEG_for_vwerp',SEG_for_vwerp);

    elseif handles.id_while == 1

        varargout{1} = handles.output;

    %     Laplace = handles.Laplace;
    %     setappdata(0,'Laplace',Laplace);
    %     centerline = handles.centerline;
    %     setappdata(0,'centerline',centerline);
    %     centerline_lapid = handles.centerline_lapid;
    %     setappdata(0,'centerline_lapid',centerline_lapid);
    %     radius = handles.radius;
    %     setappdata(0,'radius',radius);
    %     diameter = handles.diameter;
    %     setappdata(0,'diameter',diameter);
    %     axial_unit_vectors = handles.axial_unit_vectors;
    %     setappdata(0,'axial_unit_vectors',axial_unit_vectors);
    %     circumferential_unit_vectors = handles.circumferential_unit_vectors;
    %     setappdata(0,'circumferential_unit_vectors',circumferential_unit_vectors);
    %     WSS_A = handles.WSS_A;
    %     setappdata(0,'WSS_A',WSS_A);
    %     WSS_C = handles.WSS_C;
    %     setappdata(0,'WSS_C',WSS_C);
    %     mag_WSS_A = handles.mag_WSS_A;
    %     setappdata(0,'mag_WSS_A',mag_WSS_A);
    %     mag_WSS_C = handles.mag_WSS_C;
    %     setappdata(0,'mag_WSS_C',mag_WSS_C);
    %     angle_axial_direction = handles.angle_axial_direction;
    %     setappdata(0,'angle_axial_direction',angle_axial_direction);
    %     forward_velocity = handles.forward_velocity;
    %     setappdata(0,'forward_velocity',forward_velocity);
    %     backward_velocity = handles.backward_velocity;
    %     setappdata(0,'backward_velocity',backward_velocity);
    %     mag_forward_velocity = handles.mag_forward_velocity;
    %     setappdata(0,'mag_forward_velocity',mag_forward_velocity);
    %     mag_backward_velocity = handles.mag_backward_velocity;
    %     setappdata(0,'mag_backward_velocity',mag_backward_velocity);
    %     regurgitant_flow = handles.regurgitant_flow;
    %     setappdata(0,'regurgitant_flow',regurgitant_flow);
    %     centerline_flow = handles.centerline_flow;
    %     setappdata(0,'centerline_flow',centerline_flow);
    %     eccentricity = handles.eccentricity;
    %     setappdata(0,'eccentricity',eccentricity);
    %     
    %     curvature = handles.curvature;% new data Julio Sotelo
    %     setappdata(0,'curvature',curvature);% new data Julio Sotelo
    %     ellipticity = handles.ellipticity;% new data Julio Sotelo
    %     setappdata(0,'ellipticity',ellipticity);% new data Julio Sotelo
    %     length_vessel = handles.length_vessel;% new data Julio Sotelo
    %     setappdata(0,'length_vessel',length_vessel);% new data Julio Sotelo
    %     flattening = handles.flattening;% new data Julio Sotelo
    %     setappdata(0,'flattening',flattening);% new data Julio Sotelo
    % %     circulation = handles.circulation;% new data Julio Sotelo
    % %     setappdata(0,'circulation',circulation);% new data Julio Sotelo
    %     forward_vortex = handles.forward_vortex;% new data Julio Sotelo
    %     setappdata(0,'forward_vortex',forward_vortex);% new data Julio Sotelo
    %     area = handles.area;% new data Julio Sotelo
    %     setappdata(0,'area',area);% new data Julio Sotelo
    %     axial_circulation = handles.axial_circulation;% new data Julio Sotelo
    %     setappdata(0,'axial_circulation',axial_circulation);% new data Julio Sotelo

        id_while = handles.id_while;
        setappdata(0,'id_while',id_while);
        id_mesh_inlet = handles.id_mesh_inlet;
        setappdata(0,'id_mesh_inlet',id_mesh_inlet);
        id_mesh_outlet = handles.id_mesh_outlet;
        setappdata(0,'id_mesh_outlet',id_mesh_outlet);
        id_view = handles.id_view;
        setappdata(0,'id_view',id_view);
        id_inlet = handles.id_inlet;
        setappdata(0,'id_inlet',id_inlet);
        id_outlet = handles.id_outlet;
        setappdata(0,'id_outlet',id_outlet);
        inlet_exec = handles.inlet_exec;
        setappdata(0,'inlet_exec',inlet_exec);
        outlet_exec = handles.outlet_exec;
        setappdata(0,'outlet_exec',outlet_exec);
        id_ipcmra = handles.id_ipcmra;
        setappdata(0,'id_ipcmra',id_ipcmra);
        id_mag = handles.id_mag;
        setappdata(0,'id_mag',id_mag);
        slider_id_axes1 = handles.slider_axes1;
        setappdata(0,'slider_id_axes1',slider_id_axes1);
        SEG_for_vwerp = handles.SEG_for_vwerp;
        setappdata(0,'SEG_for_vwerp',SEG_for_vwerp);

        delete(handles.GUIDE_LAPLACE);
    end
else
    
    if handles.id_while == 0

        varargout{1} = handles.output;

        Laplace = handles.Laplace;
        setappdata(0,'Laplace',Laplace);
        centerline = handles.centerline;
        setappdata(0,'centerline',centerline);
        centerline_lapid = handles.centerline_lapid;
        setappdata(0,'centerline_lapid',centerline_lapid);
        radius = handles.radius;
        setappdata(0,'radius',radius);
        diameter = handles.diameter;
        setappdata(0,'diameter',diameter);
        axial_unit_vectors = handles.axial_unit_vectors;
        setappdata(0,'axial_unit_vectors',axial_unit_vectors);
        circumferential_unit_vectors = handles.circumferential_unit_vectors;
        setappdata(0,'circumferential_unit_vectors',circumferential_unit_vectors);
        WSS_A = handles.WSS_A;
        setappdata(0,'WSS_A',WSS_A);
        WSS_C = handles.WSS_C;
        setappdata(0,'WSS_C',WSS_C);
        mag_WSS_A = handles.mag_WSS_A;
        setappdata(0,'mag_WSS_A',mag_WSS_A);
        mag_WSS_C = handles.mag_WSS_C;
        setappdata(0,'mag_WSS_C',mag_WSS_C);
        angle_axial_direction = handles.angle_axial_direction;
        setappdata(0,'angle_axial_direction',angle_axial_direction);
        forward_velocity = handles.forward_velocity;
        setappdata(0,'forward_velocity',forward_velocity);
        backward_velocity = handles.backward_velocity;
        setappdata(0,'backward_velocity',backward_velocity);
        mag_forward_velocity = handles.mag_forward_velocity;
        setappdata(0,'mag_forward_velocity',mag_forward_velocity);
        mag_backward_velocity = handles.mag_backward_velocity;
        setappdata(0,'mag_backward_velocity',mag_backward_velocity);
        regurgitant_flow = handles.regurgitant_flow;
        setappdata(0,'regurgitant_flow',regurgitant_flow);
        centerline_flow = handles.centerline_flow;
        setappdata(0,'centerline_flow',centerline_flow);
        eccentricity = handles.eccentricity;
        setappdata(0,'eccentricity',eccentricity);
        
        curvature = handles.curvature;% new data Julio Sotelo
        setappdata(0,'curvature',curvature);% new data Julio Sotelo
        ellipticity = handles.ellipticity;% new data Julio Sotelo
        setappdata(0,'ellipticity',ellipticity);% new data Julio Sotelo
        length_vessel = handles.length_vessel;% new data Julio Sotelo
        setappdata(0,'length_vessel',length_vessel);% new data Julio Sotelo
        flattening = handles.flattening;% new data Julio Sotelo
        setappdata(0,'flattening',flattening);% new data Julio Sotelo
    %     circulation = handles.circulation;% new data Julio Sotelo
    %     setappdata(0,'circulation',circulation);% new data Julio Sotelo
        forward_vortex = handles.forward_vortex;% new data Julio Sotelo
        setappdata(0,'forward_vortex',forward_vortex);% new data Julio Sotelo
        area = handles.area;% new data Julio Sotelo
        setappdata(0,'area',area);% new data Julio Sotelo
        axial_circulation = handles.axial_circulation;% new data Julio Sotelo
        setappdata(0,'axial_circulation',axial_circulation);% new data Julio Sotelo

        id_while = handles.id_while;
        setappdata(0,'id_while',id_while);
        id_mesh_inlet = handles.id_mesh_inlet;
        setappdata(0,'id_mesh_inlet',id_mesh_inlet);
        id_mesh_outlet = handles.id_mesh_outlet;
        setappdata(0,'id_mesh_outlet',id_mesh_outlet);
        id_view = handles.id_view;
        setappdata(0,'id_view',id_view);
        id_inlet = handles.id_inlet;
        setappdata(0,'id_inlet',id_inlet);
        id_outlet = handles.id_outlet;
        setappdata(0,'id_outlet',id_outlet);
        inlet_exec = handles.inlet_exec;
        setappdata(0,'inlet_exec',inlet_exec);
        outlet_exec = handles.outlet_exec;
        setappdata(0,'outlet_exec',outlet_exec);
        id_ipcmra = handles.id_ipcmra;
        setappdata(0,'id_ipcmra',id_ipcmra);
        id_mag = handles.id_mag;
        setappdata(0,'id_mag',id_mag);
        slider_id_axes1 = handles.slider_axes1;
        setappdata(0,'slider_id_axes1',slider_id_axes1);

    elseif handles.id_while == 1

        varargout{1} = handles.output;

        Laplace = handles.Laplace;
        setappdata(0,'Laplace',Laplace);
        centerline = handles.centerline;
        setappdata(0,'centerline',centerline);
        centerline_lapid = handles.centerline_lapid;
        setappdata(0,'centerline_lapid',centerline_lapid);
        radius = handles.radius;
        setappdata(0,'radius',radius);
        diameter = handles.diameter;
        setappdata(0,'diameter',diameter);
        axial_unit_vectors = handles.axial_unit_vectors;
        setappdata(0,'axial_unit_vectors',axial_unit_vectors);
        circumferential_unit_vectors = handles.circumferential_unit_vectors;
        setappdata(0,'circumferential_unit_vectors',circumferential_unit_vectors);
        WSS_A = handles.WSS_A;
        setappdata(0,'WSS_A',WSS_A);
        WSS_C = handles.WSS_C;
        setappdata(0,'WSS_C',WSS_C);
        mag_WSS_A = handles.mag_WSS_A;
        setappdata(0,'mag_WSS_A',mag_WSS_A);
        mag_WSS_C = handles.mag_WSS_C;
        setappdata(0,'mag_WSS_C',mag_WSS_C);
        angle_axial_direction = handles.angle_axial_direction;
        setappdata(0,'angle_axial_direction',angle_axial_direction);
        forward_velocity = handles.forward_velocity;
        setappdata(0,'forward_velocity',forward_velocity);
        backward_velocity = handles.backward_velocity;
        setappdata(0,'backward_velocity',backward_velocity);
        mag_forward_velocity = handles.mag_forward_velocity;
        setappdata(0,'mag_forward_velocity',mag_forward_velocity);
        mag_backward_velocity = handles.mag_backward_velocity;
        setappdata(0,'mag_backward_velocity',mag_backward_velocity);
        regurgitant_flow = handles.regurgitant_flow;
        setappdata(0,'regurgitant_flow',regurgitant_flow);
        centerline_flow = handles.centerline_flow;
        setappdata(0,'centerline_flow',centerline_flow);
        eccentricity = handles.eccentricity;
        setappdata(0,'eccentricity',eccentricity);
        
        curvature = handles.curvature;% new data Julio Sotelo
        setappdata(0,'curvature',curvature);% new data Julio Sotelo
        ellipticity = handles.ellipticity;% new data Julio Sotelo
        setappdata(0,'ellipticity',ellipticity);% new data Julio Sotelo
        length_vessel = handles.length_vessel;% new data Julio Sotelo
        setappdata(0,'length_vessel',length_vessel);% new data Julio Sotelo
        flattening = handles.flattening;% new data Julio Sotelo
        setappdata(0,'flattening',flattening);% new data Julio Sotelo
    %     circulation = handles.circulation;% new data Julio Sotelo
    %     setappdata(0,'circulation',circulation);% new data Julio Sotelo
        forward_vortex = handles.forward_vortex;% new data Julio Sotelo
        setappdata(0,'forward_vortex',forward_vortex);% new data Julio Sotelo
        area = handles.area;% new data Julio Sotelo
        setappdata(0,'area',area);% new data Julio Sotelo
        axial_circulation = handles.axial_circulation;% new data Julio Sotelo
        setappdata(0,'axial_circulation',axial_circulation);% new data Julio Sotelo

        id_while = handles.id_while;
        setappdata(0,'id_while',id_while);
        id_mesh_inlet = handles.id_mesh_inlet;
        setappdata(0,'id_mesh_inlet',id_mesh_inlet);
        id_mesh_outlet = handles.id_mesh_outlet;
        setappdata(0,'id_mesh_outlet',id_mesh_outlet);
        id_view = handles.id_view;
        setappdata(0,'id_view',id_view);
        id_inlet = handles.id_inlet;
        setappdata(0,'id_inlet',id_inlet);
        id_outlet = handles.id_outlet;
        setappdata(0,'id_outlet',id_outlet);
        inlet_exec = handles.inlet_exec;
        setappdata(0,'inlet_exec',inlet_exec);
        outlet_exec = handles.outlet_exec;
        setappdata(0,'outlet_exec',outlet_exec);
        id_ipcmra = handles.id_ipcmra;
        setappdata(0,'id_ipcmra',id_ipcmra);
        id_mag = handles.id_mag;
        setappdata(0,'id_mag',id_mag);
        slider_id_axes1 = handles.slider_axes1;
        setappdata(0,'slider_id_axes1',slider_id_axes1);

        delete(handles.GUIDE_LAPLACE);
    end
end
% Get default command line output from handles structure


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel2,'BackgroundColor',[0,0,0])
switch get(handles.popupmenu1,'Value')
    
    case 1
        cla(handles.axes1,'reset');
        set(handles.axes1,'visible','off')
        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        handles.id_ipcmra = 0;
        handles.id_mag = 0;
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        
    case 2
        
        set(handles.popupmenu2,'Value',1)
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        if handles.id_view == 1
            
            
            handles.slider_axes1 = round(size(handles.IPCMRA,3)/2); %Julio
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
                
            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_view == 2

            handles.slider_axes1 = round(size(handles.IPCMRA,2)/2); %Julio
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_view == 3
            
            handles.slider_axes1 = round(size(handles.IPCMRA,1)/2); %Julio
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        end
        
        
        % show the axes1
        set(handles.slider1,'visible','on')
        set(handles.text1,'visible','on')
        
        % id IPCMRA =1
        handles.id_image = 1;
        handles.id_ipcmra = 1;
        handles.id_mag = 0;
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;

        
     case 3
         
        set(handles.popupmenu2,'Value',1)
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        
        if handles.id_view == 1
            
            handles.slider_axes1 = round(size(handles.IPCMRA,3)/2); %Julio
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');

            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_view == 2
            
            handles.slider_axes1 = round(size(handles.IPCMRA,2)/2); %Julio

            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_view == 3
            
            handles.slider_axes1 = round(size(handles.IPCMRA,1)/2); %Julio
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');

            handles.list_string = {'Sagital View','Axial View','Coronal View'};
            set(handles.popupmenu3,'String',handles.list_string,'Value',handles.id_view);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(handles.axes1,'Visible','on');
            axis off
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        end
        
        % show the axes1
        set(handles.slider1,'visible','on')
        set(handles.text1,'visible','on')
        
        % set the id
        handles.id_image = 2;
        handles.id_ipcmra = 0;
        handles.id_mag = 1;
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;

end

handles.output = hObject;  
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.id_ipcmra == 1
        if handles.id_view == 1
            
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.c;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.c;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;
            
            
        elseif handles.id_view == 2
            
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.b;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.b;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;
        
        elseif handles.id_view == 3
            
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.a;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.a;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;
            
            
        end
        
        
    elseif handles.id_mag == 1
         
        if handles.id_view == 1
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.c;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.c;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;
            
        elseif handles.id_view == 2
            
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.b;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.b;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;
            
        elseif handles.id_view == 3
            
            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            pp=1/handles.a;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.a;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes1 = handles.slider_value;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;
            
        end
        
    elseif handles.id_mesh_vel == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),5, 'Linewidth',1);
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.mags_vel(:));
        handles.max_vel = max(handles.mags_vel(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Velocity [m/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        
    elseif handles.id_wss_a == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_A(1:end,1,handles.peak_flow),handles.WSS_A(1:end,2,handles.peak_flow),handles.WSS_A(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.WSS_A,1) size(handles.WSS_A,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_wssa = min(handles.mag_WSS_A(:));
        handles.max_wssa = max(handles.mag_WSS_A(:));
        handles.mean_wssa = (handles.min_wssa + handles.max_wssa)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_wssa handles.max_wssa];
        c.Ticks = [handles.min_wssa, (handles.min_wssa + handles.mean_wssa)/2, handles.mean_wssa, (handles.max_wssa + handles.mean_wssa)/2, handles.max_wssa];
        c.TickLabels = {num2str(handles.min_wssa,'%0.2f'), num2str((handles.min_wssa + handles.mean_wssa)/2,'%0.2f'), num2str(handles.mean_wssa,'%0.2f'), num2str((handles.max_wssa + handles.mean_wssa)/2,'%0.2f'), num2str(handles.max_wssa,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'WSS-A [N/m^{2}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
        
        
    elseif handles.id_wss_c == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.WSS_C,1) size(handles.WSS_C,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_wssc = min(handles.mag_WSS_C(:));
        handles.max_wssc = max(handles.mag_WSS_C(:));
        handles.mean_wssc = (handles.min_wssc + handles.max_wssc)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_wssc handles.max_wssc];
        c.Ticks = [handles.min_wssc, (handles.min_wssc + handles.mean_wssc)/2, handles.mean_wssc, (handles.max_wssc + handles.mean_wssc)/2, handles.max_wssc];
        c.TickLabels = {num2str(handles.min_wssc,'%0.2f'), num2str((handles.min_wssc + handles.mean_wssc)/2,'%0.2f'), num2str(handles.mean_wssc,'%0.2f'), num2str((handles.max_wssc + handles.mean_wssc)/2,'%0.2f'), num2str(handles.max_wssc,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'WSS-C [N/m^{2}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
        
    elseif handles.id_axial_angle == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        c = colorbar(handles.axes1);
        handles.min_axa = min(handles.angle_axial_direction(:));
        handles.max_axa = max(handles.angle_axial_direction(:));
        handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_axa handles.max_axa];
        c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
        c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Axial Angle [^{O}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_axa handles.max_axa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')        
        
    elseif handles.id_forward_vel == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,handles.peak_flow),handles.forward_velocity(1:end,2,handles.peak_flow),handles.forward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_fvel = min(handles.mag_forward_velocity(:));
        handles.max_fvel = max(handles.mag_forward_velocity(:));
        handles.mean_fvel = (handles.min_fvel + handles.max_fvel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_fvel handles.max_fvel];
        c.Ticks = [handles.min_fvel, (handles.min_fvel + handles.mean_fvel)/2, handles.mean_fvel, (handles.max_fvel + handles.mean_fvel)/2, handles.max_fvel];
        c.TickLabels = {num2str(handles.min_fvel,'%0.2f'), num2str((handles.min_fvel + handles.mean_fvel)/2,'%0.2f'), num2str(handles.mean_fvel,'%0.2f'), num2str((handles.max_fvel + handles.mean_fvel)/2,'%0.2f'), num2str(handles.max_fvel,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Forward Velocity [m/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

    elseif handles.id_backward_vel == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,handles.peak_flow),handles.backward_velocity(1:end,2,handles.peak_flow),handles.backward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_bvel = min(handles.mag_backward_velocity(:));
        handles.max_bvel = max(handles.mag_backward_velocity(:));
        handles.mean_bvel = (handles.min_bvel + handles.max_bvel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_bvel handles.max_bvel];
        c.Ticks = [handles.min_bvel, (handles.min_bvel + handles.mean_bvel)/2, handles.mean_bvel, (handles.max_bvel + handles.mean_bvel)/2, handles.max_bvel];
        c.TickLabels = {num2str(handles.min_bvel,'%0.2f'), num2str((handles.min_bvel + handles.mean_bvel)/2,'%0.2f'), num2str(handles.mean_bvel,'%0.2f'), num2str((handles.max_bvel + handles.mean_bvel)/2,'%0.2f'), num2str(handles.max_bvel,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Backward Velocity [m/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
        
        
%     elseif handles.id_circulation == 1
%         
%         % show the axes1
%         set(handles.axes1,'Visible','on');
%         axis off
%         
%         pp=1/size(handles.veset,3);
%         slider_step(1) = pp;
%         slider_step(2) = 0.1;
%         set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
%         maxslice = size(handles.veset,3);
%         handles.slider_value = get(hObject,'Value');
%         if handles.slider_value<=pp
%                 handles.slider_value = 1;
%         else
%                 handles.slider_value = round(handles.slider_value*maxslice);
%         end
%         handles.peak_flow = handles.slider_value;
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
%         colormap(handles.axes1,'winter');
%         c = colorbar(handles.axes1);
%         handles.min_cir = min(handles.circulation(:));
%         handles.max_cir = max(handles.circulation(:));
%         handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%         c.LimitsMode = 'manual';
%         c.Limits = [handles.min_cir handles.max_cir];
%         c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%         c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%         c.Color = [1 1 1]; % color
%         c.Location = 'manual';
%         c.Position = [0.2 0.2 0.02 0.3];
%         c.FontWeight = 'bold';
%         c.Label.String = 'Circulation [mm^{2}/s]';
%         windows_screen_size = get(0,'ScreenSize');
%         if windows_screen_size(4)<=900
%             c.FontSize = 9;
%         elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%             c.FontSize = 10;
%         else
%             c.FontSize = 11;
%         end
%         if windows_screen_size(4)<=900
%             c.Label.FontSize = 10;
%         elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%             c.Label.FontSize = 12;
%         else
%             c.Label.FontSize = 14;
%         end
%         c.Label.FontWeight = 'bold';
%         c.Label.Color = [1 1 1];% color
%         caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%         hold off
%         axis vis3d
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
%         
        
    elseif handles.id_forward_vortex == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.forward_vortex(:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        c = colorbar(handles.axes1);
        handles.min_fov = min(handles.forward_vortex(:));
        handles.max_fov = max(handles.forward_vortex(:));
        handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_fov handles.max_fov];
        c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
        c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Forward Vortex [1/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_fov handles.max_fov]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
    
    elseif handles.id_axial_circulation == 1
        
        % show the axes1
        set(handles.axes1,'Visible','on');
        axis off
        
        pp=1/size(handles.veset,3);
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = size(handles.veset,3);
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.peak_flow = handles.slider_value;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.axial_circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'autumn');
        c = colorbar(handles.axes1);
        handles.min_acir = min(handles.axial_circulation(:));
        handles.max_acir = max(handles.axial_circulation(:));
        handles.mean_acir = (handles.min_acir + handles.max_acir)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_acir handles.max_acir];
        c.Ticks = [handles.min_acir, (handles.min_acir + handles.mean_acir)/2, handles.mean_acir, (handles.max_acir + handles.mean_acir)/2, handles.max_acir];
        c.TickLabels = {num2str(handles.min_acir,'%0.2f'), num2str((handles.min_acir + handles.mean_acir)/2,'%0.2f'), num2str(handles.mean_acir,'%0.2f'), num2str((handles.max_acir + handles.mean_acir)/2,'%0.2f'), num2str(handles.max_acir,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Axial Circulation [cm^{2}/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_acir handles.max_acir]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
        
    end

handles.output = hObject;  
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.id_vwerp  == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Sagital view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if handles.id_view == 1

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [mean(PUNTOS_C),handles.slider_axes1];
        plot(center(1),center(2),'*g')
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(2),center(1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

    %     center_ori = [center(2),center(1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [point1(1:2),center(3)]; 
        p2       = [center(1:2),center(3) + sqrt(sum((center-point1).^2))]; 
        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles

            Rz = [cosd(n) , -sind(n) , 0 ; sind(n) , cosd(n) , 0 ; 0 , 0 , 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1*Rz;
            VB = (VB/norm(VB))*20;
            VG = v2;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1) , 0 ; sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; 0 , 0 , 1];

        % genero vectores cambiando r y el azul
        VB = v1*Rz;
        VB = VB/norm(VB);
        VR = normal*Rz;
        VR = VR/norm(VR);
        VG = v2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB*20;
        v2_n = VG;
        normal_n = VR*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB_n = v1_n;               
            VB_n = VB_n/norm(VB_n)*20;

            v1_n = (v1_n/norm(v1_n));
            ux = v1_n(1); 
            uy = v1_n(2); 
            uz = v1_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VG_n = (v2_n/norm(v2_n))*Raxis; 
            VG_n = (VG_n/norm(VG_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification

            Area2(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        ux = v1_n(1); 
        uy = v1_n(2); 
        uz = v1_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (v2_n/norm(v2_n))*Raxis; 
        VG_n = (VG_n/norm(VG_n))*20;

        VR_n = normal_n*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        VB_n = cross(VG_n/norm(VG_n),VR_n/norm(VR_n));
        VB_n = VB_n/norm(VB_n)*20;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);     
        handles.faceid = faceid;
        handles.cutpos = cutpos;
        Selected_element = handles.elem(elemid,:);



        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            [r,~,~] = find(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end

        [r,~,~] = find(nodesid==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        % Julio Sotelo
        handles.id_inlet = unique(handles.elem(elemid,:));
        
%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_inlet1 = r;
%         idnodes_inlet2 = handles.nodes_id(id_temp~=0);
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
%         plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_inlet = idnodes_inlet1;
% 
%             case 'Green'
% 
%               handles.id_inlet = idnodes_inlet2;
% 
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
%         distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(handles.id_inlet,1)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
%         quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
    %     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(2);
            Y = Y*handles.voxel_MR(1);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    % 
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
    %         quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_in(2) + b1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_INLET)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        close(h)
        
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        end

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  

        % we generate the next variables for vwerp
        inlet = handles.PLANE_INLET;
        n_inlet = handles.normal_inlet/norm(handles.normal_inlet);
        
        handles.inlet_for_plot = inlet;
        handles.n_inlet_for_plot = n_inlet;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Axial   view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 2

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(3),(handles.slider_axes1-1)*handles.voxel_MR(2),center(2)];
        point1 = [PUNTOS_C(1,2),(handles.slider_axes1-1)*handles.voxel_MR(2),PUNTOS_C(1,1)];
        point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     center_ori = [center(3),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),center(2)];
    %     point1 = [PUNTOS_C(1,2),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
    %     plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
    %     plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d
    %     xlabel('ejex')
    %     ylabel('ejey')
    %     zlabel('ejez')

    %     figure,
    %     imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
    %     hold on
    %     Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %     himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %     cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %     cdata = double(cdata)*0.5;
    %     set(himage, 'AlphaData', cdata);
    %     plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
    %     center = [handles.slider_axes1,mean(PUNTOS_C)];
    %     plot(center(2),center(3),'*g');
    %     hold off
    %     axis image
    %     set(handles.axes1,'xticklabel',[],'yticklabel',[])
    %     axis off
    %     colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'*c')
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1      = [point1(1), center(2), point1(3)]; 
        p2      = [center(1), center(2) + sqrt(sum((center-p1).^2)), center(3)];

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(p1(1),p1(2),p1(3),'*b')
    %     plot3(p2(1),p2(2),p2(3),'*g')
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'.c')
    %     quiver3(center(1),center(2),center(3),v1(1),v1(2),v1(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),v2(1),v2(2),v2(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),normal(1),normal(2),normal(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles


    %         Rx = [1 ,       0 ,        0 ; ... 
    %               0 , cosd(n) , -sind(n) ; ...
    %               0 , sind(n) , cosd(n)] ;

    %         Ry = [cosd(n),  0 ,  sind(n) ; ...
    %               0      ,  1 ,        0 ; ...
    %               -sind(n), 0 , cosd(n)] ;

            Rz = [cosd(n) , -sind(n), 0 ; ...
                  sind(n) , cosd(n) , 0 ; ...
                  0       ,        0, 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1;
            VB = (VB/norm(VB))*20;
            VG = v2*Rz;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1), 0 ; ...
              sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; ...
                            0       ,                      0, 1];

        % genero vectores cambiando r y el azul
        VB = v1/norm(v1);
        VB = (VB/norm(VB))*20;
        VG = (v2/norm(v2))*Rz;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rz;
        VR = VR/norm(VR)*20;

        PB = center + VB;  
        PG = center + VG;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(PB(1),PB(2),PB(3),'*b')
    %     plot3(PG(1),PG(2),PG(3),'*g')
    %     quiver3(center(1),center(2),center(3),VB(1),VB(2),VB(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VG(1),VG(2),VG(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VR(1),VR(2),VR(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

    %     handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        % Julio Sotelo
        handles.id_inlet = unique(handles.elem(elemid,:));
        
%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_inlet1 = r;
%         idnodes_inlet2 = handles.nodes_id(id_temp~=0);
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
%         plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_inlet = idnodes_inlet1;
% 
%             case 'Green'
% 
%               handles.id_inlet = idnodes_inlet2;
% 
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
%         distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(handles.id_inlet,1)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

        
%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
%         quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
    %     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(2);
            Y = Y*handles.voxel_MR(1);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    % 
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
    %         quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_in(2) + b1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_INLET)
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)
  
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        end


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  

        % we generate the next variables for vwerp
        inlet = handles.PLANE_INLET;
        n_inlet = handles.normal_inlet/norm(handles.normal_inlet);
        
        handles.inlet_for_plot = inlet;
        handles.n_inlet_for_plot = n_inlet;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Coronal view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 3

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [(handles.slider_axes1-1)*handles.voxel_MR(1),center(3),center(2)];
        point1 = [(handles.slider_axes1-1)*handles.voxel_MR(1),PUNTOS_C(1,2),PUNTOS_C(1,1)];
        point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

    %     center_ori = [mean(mean(handles.yd(handles.slider_axes1,:,:))),center(3),center(2)];
    %     point1 = [mean(mean(handles.yd(handles.slider_axes1,:,:))),PUNTOS_C(1,2),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [center(1), point1(2:3)]; 
        p2       = [center(1) + sqrt(sum((center-point1).^2)), center(2:3)]; 

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles


            Rx = [1 , 0 , 0 ; 0 , cosd(n) , -sind(n) ; 0 , sind(n), cosd(n)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v2;
            VB = (VB/norm(VB))*20;
            VG = v1*Rx;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rx = [1 , 0 , 0 ; 0 , cosd(Angle_selected1) , -sind(Angle_selected1) ; 0 , sind(Angle_selected1), cosd(Angle_selected1)];

        % genero vectores cambiando r y el azul
        VB = v2/norm(v2);
        VB = (VB/norm(VB))*20;
        VG = (v1/norm(v1))*Rx;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rx;
        VR = VR/norm(VR)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

    %     handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        % Julio Sotelo
        handles.id_inlet = unique(handles.elem(elemid,:));
        
%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_inlet1 = r;
%         idnodes_inlet2 = handles.nodes_id(id_temp~=0);
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
%         plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_inlet = idnodes_inlet1;
% 
%             case 'Green'
% 
%               handles.id_inlet = idnodes_inlet2;
% 
%         end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
%         distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(handles.id_inlet,1)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

        
%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
%         plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
%         quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
    %     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
            
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(2);
            Y = Y*handles.voxel_MR(1);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    % 
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
    %         quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (b1*center_in(2) + a1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);

    %         imlook3d(handles.SEG)
    %         imlook3d(double(BW2))
    %         imlook3d(double(OUT<=0))
    %         
    %         disp(size(handles.SEG))
    %         disp(size(double(BW2)))
    %         disp(size(double(OUT<=0)))
    %         
    %         BW2  = permute(BW2,[2,1,3]);
    %         OUT  = permute(OUT,[2,1,3]);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_INLET)
        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',6);
        end


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  
        
        % we generate the next variables for vwerp
        inlet = handles.PLANE_INLET;
        n_inlet = handles.normal_inlet/norm(handles.normal_inlet);
        
        handles.inlet_for_plot = inlet;
        handles.n_inlet_for_plot = n_inlet;

    end

else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Sagital view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if handles.id_view == 1

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [mean(PUNTOS_C),handles.slider_axes1];
        plot(center(1),center(2),'*g')
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(2),center(1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

    %     center_ori = [center(2),center(1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [point1(1:2),center(3)]; 
        p2       = [center(1:2),center(3) + sqrt(sum((center-point1).^2))]; 
        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles

            Rz = [cosd(n) , -sind(n) , 0 ; sind(n) , cosd(n) , 0 ; 0 , 0 , 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1*Rz;
            VB = (VB/norm(VB))*20;
            VG = v2;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1) , 0 ; sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; 0 , 0 , 1];

        % genero vectores cambiando r y el azul
        VB = v1*Rz;
        VB = VB/norm(VB);
        VR = normal*Rz;
        VR = VR/norm(VR);
        VG = v2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB*20;
        v2_n = VG;
        normal_n = VR*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB_n = v1_n;               
            VB_n = VB_n/norm(VB_n)*20;

            v1_n = (v1_n/norm(v1_n));
            ux = v1_n(1); 
            uy = v1_n(2); 
            uz = v1_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VG_n = (v2_n/norm(v2_n))*Raxis; 
            VG_n = (VG_n/norm(VG_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification

            Area2(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        ux = v1_n(1); 
        uy = v1_n(2); 
        uz = v1_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (v2_n/norm(v2_n))*Raxis; 
        VG_n = (VG_n/norm(VG_n))*20;

        VR_n = normal_n*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        VB_n = cross(VG_n/norm(VG_n),VR_n/norm(VR_n));
        VB_n = VB_n/norm(VB_n)*20;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);     
        handles.faceid = faceid;
        handles.cutpos = cutpos;
        Selected_element = handles.elem(elemid,:);



        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            [r,~,~] = find(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end

        [r,~,~] = find(nodesid==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_inlet1 = r;
        idnodes_inlet2 = handles.nodes_id(id_temp~=0);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
        plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_inlet = idnodes_inlet1;

            case 'Green'

              handles.id_inlet = idnodes_inlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
        distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(id_nodes_inlet,:)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
    %     plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
    %     quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])
    %     

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];

        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
        %     quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_in(2) + b1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        end

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Axial   view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 2

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(3),(handles.slider_axes1-1)*handles.voxel_MR(2),center(2)];
        point1 = [PUNTOS_C(1,2),(handles.slider_axes1-1)*handles.voxel_MR(2),PUNTOS_C(1,1)];
        point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     center_ori = [center(3),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),center(2)];
    %     point1 = [PUNTOS_C(1,2),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
    %     plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
    %     plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d
    %     xlabel('ejex')
    %     ylabel('ejey')
    %     zlabel('ejez')

    %     figure,
    %     imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
    %     hold on
    %     Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %     himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %     cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %     cdata = double(cdata)*0.5;
    %     set(himage, 'AlphaData', cdata);
    %     plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
    %     center = [handles.slider_axes1,mean(PUNTOS_C)];
    %     plot(center(2),center(3),'*g');
    %     hold off
    %     axis image
    %     set(handles.axes1,'xticklabel',[],'yticklabel',[])
    %     axis off
    %     colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'*c')
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1      = [point1(1), center(2), point1(3)]; 
        p2      = [center(1), center(2) + sqrt(sum((center-p1).^2)), center(3)];

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(p1(1),p1(2),p1(3),'*b')
    %     plot3(p2(1),p2(2),p2(3),'*g')
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'.c')
    %     quiver3(center(1),center(2),center(3),v1(1),v1(2),v1(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),v2(1),v2(2),v2(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),normal(1),normal(2),normal(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles


    %         Rx = [1 ,       0 ,        0 ; ... 
    %               0 , cosd(n) , -sind(n) ; ...
    %               0 , sind(n) , cosd(n)] ;

    %         Ry = [cosd(n),  0 ,  sind(n) ; ...
    %               0      ,  1 ,        0 ; ...
    %               -sind(n), 0 , cosd(n)] ;

            Rz = [cosd(n) , -sind(n), 0 ; ...
                  sind(n) , cosd(n) , 0 ; ...
                  0       ,        0, 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1;
            VB = (VB/norm(VB))*20;
            VG = v2*Rz;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1), 0 ; ...
              sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; ...
                            0       ,                      0, 1];

        % genero vectores cambiando r y el azul
        VB = v1/norm(v1);
        VB = (VB/norm(VB))*20;
        VG = (v2/norm(v2))*Rz;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rz;
        VR = VR/norm(VR)*20;

        PB = center + VB;  
        PG = center + VG;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(PB(1),PB(2),PB(3),'*b')
    %     plot3(PG(1),PG(2),PG(3),'*g')
    %     quiver3(center(1),center(2),center(3),VB(1),VB(2),VB(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VG(1),VG(2),VG(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VR(1),VR(2),VR(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

    %     handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_inlet1 = r;
        idnodes_inlet2 = handles.nodes_id(id_temp~=0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
        plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_inlet = idnodes_inlet1;

            case 'Green'

              handles.id_inlet = idnodes_inlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
        distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(id_nodes_inlet,:)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
    %     plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
    %     quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
    
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];

        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
        %     quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_in(2) + b1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));

        %     imlook3d(handles.SEG)
        %     imlook3d(PLANE_INLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        end

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Coronal view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 3

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [(handles.slider_axes1-1)*handles.voxel_MR(1),center(3),center(2)];
        point1 = [(handles.slider_axes1-1)*handles.voxel_MR(1),PUNTOS_C(1,2),PUNTOS_C(1,1)];
        point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

    %     center_ori = [mean(mean(handles.yd(handles.slider_axes1,:,:))),center(3),center(2)];
    %     point1 = [mean(mean(handles.yd(handles.slider_axes1,:,:))),PUNTOS_C(1,2),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [center(1), point1(2:3)]; 
        p2       = [center(1) + sqrt(sum((center-point1).^2)), center(2:3)]; 

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles


            Rx = [1 , 0 , 0 ; 0 , cosd(n) , -sind(n) ; 0 , sind(n), cosd(n)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v2;
            VB = (VB/norm(VB))*20;
            VG = v1*Rx;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rx = [1 , 0 , 0 ; 0 , cosd(Angle_selected1) , -sind(Angle_selected1) ; 0 , sind(Angle_selected1), cosd(Angle_selected1)];

        % genero vectores cambiando r y el azul
        VB = v2/norm(v2);
        VB = (VB/norm(VB))*20;
        VG = (v1/norm(v1))*Rx;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rx;
        VR = VR/norm(VR)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

    %     handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_inlet1 = r;
        idnodes_inlet2 = handles.nodes_id(id_temp~=0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_inlet1,1),handles.nodes(idnodes_inlet1,2),handles.nodes(idnodes_inlet1,3),'*r')
        plot3(handles.nodes(idnodes_inlet2,1),handles.nodes(idnodes_inlet2,2),handles.nodes(idnodes_inlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_inlet = idnodes_inlet1;

            case 'Green'

              handles.id_inlet = idnodes_inlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal inlet
        normal_in = VR_n; % normal of the slice 
        center_in = center; % center of the section inlet
        distance_in = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_inlet = handles.id_inlet(distance_in(handles.id_inlet)<10);
        center_in_n = mean(handles.nodes(id_nodes_inlet,:)); % center of the nurb of point inlet
        normal_in_n = (center_in_n - center_in)/norm(center_in_n - center_in)*20;

        angle_b_normal = acos(sum(normal_in_n.*normal_in)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_inlet = (normal_in/norm(normal_in))*-1;
        else
            handles.normal_inlet = (normal_in/norm(normal_in));
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*g')
    %     plot3(center_in(1),center_in(2),center_in(3),'*b','Linewidth',5)
    %     quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
    
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
        %  
        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_in(2),center_in(1),center_in(3),'*c','Linewidth',3)
        %     quiver3(center_in(2),center_in(1),center_in(3),handles.normal_inlet(2),handles.normal_inlet(1),handles.normal_inlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_inlet(2);
            b1 = handles.normal_inlet(1);
            c1 = handles.normal_inlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_in(2) + b1*center_in(1) + c1*center_in(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_in(2)).^2 + (Y_C-center_in(1)).^2 + (Z_C-center_in(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_INLET = double(LL == max(POINT(:).*LL(:)));

        %     imlook3d(handles.SEG)
        %     imlook3d(PLANE_INLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet','Mesh + Outlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        else
            handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Mesh + Inlet'};
            set(handles.popupmenu2,'String',handles.list_string,'Value',4);
        end

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 1;
        handles.id_mesh_outlet = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.inlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.outlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else 

            set(handles.pushbutton3,'Visible','off');

        end  

    end


end


handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.id_vwerp == 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Sagital view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if handles.id_view == 1

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [mean(PUNTOS_C),handles.slider_axes1];
        plot(center(1),center(2),'*g')
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(2),center(1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

    %     center_ori = [center(2),center(1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [point1(1:2),center(3)]; 
        p2       = [center(1:2),center(3) + sqrt(sum((center-point1).^2))]; 
        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles

            Rz = [cosd(n) , -sind(n) , 0 ; sind(n) , cosd(n) , 0 ; 0 , 0 , 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1*Rz;
            VB = (VB/norm(VB))*20;
            VG = v2;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1) , 0 ; sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; 0 , 0 , 1];

        % genero vectores cambiando r y el azul
        VB = v1*Rz;
        VB = VB/norm(VB);
        VR = normal*Rz;
        VR = VR/norm(VR);
        VG = v2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB*20;
        v2_n = VG;
        normal_n = VR*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB_n = v1_n;               
            VB_n = VB_n/norm(VB_n)*20;

            v1_n = (v1_n/norm(v1_n));
            ux = v1_n(1); 
            uy = v1_n(2); 
            uz = v1_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VG_n = (v2_n/norm(v2_n))*Raxis; 
            VG_n = (VG_n/norm(VG_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification

            Area2(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        ux = v1_n(1); 
        uy = v1_n(2); 
        uz = v1_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (v2_n/norm(v2_n))*Raxis; 
        VG_n = (VG_n/norm(VG_n))*20;

        VR_n = normal_n*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        VB_n = cross(VG_n/norm(VG_n),VR_n/norm(VR_n));
        VB_n = VB_n/norm(VB_n)*20;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);     
        handles.faceid = faceid;
        handles.cutpos = cutpos;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        
        % Julio Sotelo
        handles.id_outlet = unique(handles.elem(elemid,:));
 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_outlet1 = r;
%         idnodes_outlet2 = handles.nodes_id(id_temp~=0);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
%         plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_outlet = idnodes_outlet1;
% 
%             case 'Green'
% 
%               handles.id_outlet = idnodes_outlet2;
% 
%         end
% 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
%         distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(handles.id_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*g')
%         plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
%         quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
    
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    %      
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
    %         quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)
        
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',7);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end

        % we generate the next variables for vwerp
        outlet = handles.PLANE_OUTLET;
        n_outlet = handles.normal_outlet/norm(handles.normal_outlet);
        
        handles.outlet_for_plot = outlet;
        handles.n_outlet_for_plot = n_outlet;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Axial   view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 2

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(3),(handles.slider_axes1-1)*handles.voxel_MR(2),center(2)];
        point1 = [PUNTOS_C(1,2),(handles.slider_axes1-1)*handles.voxel_MR(2),PUNTOS_C(1,1)];
        point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     center_ori = [center(3),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),center(2)];
    %     point1 = [PUNTOS_C(1,2),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
    %     plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
    %     plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d
    %     xlabel('ejex')
    %     ylabel('ejey')
    %     zlabel('ejez')

    %     figure,
    %     imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
    %     hold on
    %     Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %     himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %     cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %     cdata = double(cdata)*0.5;
    %     set(himage, 'AlphaData', cdata);
    %     plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
    %     center = [handles.slider_axes1,mean(PUNTOS_C)];
    %     plot(center(2),center(3),'*g');
    %     hold off
    %     axis image
    %     set(handles.axes1,'xticklabel',[],'yticklabel',[])
    %     axis off
    %     colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'*c')
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1      = [point1(1), center(2), point1(3)]; 
        p2      = [center(1), center(2) + sqrt(sum((center-p1).^2)), center(3)];

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(p1(1),p1(2),p1(3),'*b')
    %     plot3(p2(1),p2(2),p2(3),'*g')
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'.c')
    %     quiver3(center(1),center(2),center(3),v1(1),v1(2),v1(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),v2(1),v2(2),v2(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),normal(1),normal(2),normal(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles


    %         Rx = [1 ,       0 ,        0 ; ... 
    %               0 , cosd(n) , -sind(n) ; ...
    %               0 , sind(n) , cosd(n)] ;

    %         Ry = [cosd(n),  0 ,  sind(n) ; ...
    %               0      ,  1 ,        0 ; ...
    %               -sind(n), 0 , cosd(n)] ;

            Rz = [cosd(n) , -sind(n), 0 ; ...
                  sind(n) , cosd(n) , 0 ; ...
                  0       ,        0, 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1;
            VB = (VB/norm(VB))*20;
            VG = v2*Rz;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1), 0 ; ...
              sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; ...
                            0       ,                      0, 1];

        % genero vectores cambiando r y el azul
        VB = v1/norm(v1);
        VB = (VB/norm(VB))*20;
        VG = (v2/norm(v2))*Rz;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rz;
        VR = VR/norm(VR)*20;

        PB = center + VB;  
        PG = center + VG;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(PB(1),PB(2),PB(3),'*b')
    %     plot3(PG(1),PG(2),PG(3),'*g')
    %     quiver3(center(1),center(2),center(3),VB(1),VB(2),VB(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VG(1),VG(2),VG(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VR(1),VR(2),VR(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

         % Julio Sotelo
        handles.id_outlet = unique(handles.elem(elemid,:));
 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_outlet1 = r;
%         idnodes_outlet2 = handles.nodes_id(id_temp~=0);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
%         plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_outlet = idnodes_outlet1;
% 
%             case 'Green'
% 
%               handles.id_outlet = idnodes_outlet2;
% 
%         end
% 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
%         distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(handles.id_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*g')
%         plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
%         quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
    
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    %      
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
    %         quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)
        
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',7);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end

        % we generate the next variables for vwerp
        outlet = handles.PLANE_OUTLET;
        n_outlet = handles.normal_outlet/norm(handles.normal_outlet);
        
        handles.outlet_for_plot = outlet;
        handles.n_outlet_for_plot = n_outlet;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Coronal view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 3

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [(handles.slider_axes1-1)*handles.voxel_MR(1),center(3),center(2)];
        point1 = [(handles.slider_axes1-1)*handles.voxel_MR(1),PUNTOS_C(1,2),PUNTOS_C(1,1)];
        point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

    %     center_ori = [mean(mean(handles.yd(handles.slider_axes1,:,:))),center(3),center(2)];
    %     point1 = [mean(mean(handles.yd(handles.slider_axes1,:,:))),PUNTOS_C(1,2),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [center(1), point1(2:3)]; 
        p2       = [center(1) + sqrt(sum((center-point1).^2)), center(2:3)]; 

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles


            Rx = [1 , 0 , 0 ; 0 , cosd(n) , -sind(n) ; 0 , sind(n), cosd(n)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v2;
            VB = (VB/norm(VB))*20;
            VG = v1*Rx;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rx = [1 , 0 , 0 ; 0 , cosd(Angle_selected1) , -sind(Angle_selected1) ; 0 , sind(Angle_selected1), cosd(Angle_selected1)];

        % genero vectores cambiando r y el azul
        VB = v2/norm(v2);
        VB = (VB/norm(VB))*20;
        VG = (v1/norm(v1))*Rx;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rx;
        VR = VR/norm(VR)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

                % Julio Sotelo
        handles.id_outlet = unique(handles.elem(elemid,:));
 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         elem_new = handles.elem;
%         elem_new(elemid,:) = []; 
% 
%         [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
%         nodesid = zeros(size(handles.nodes,1),1);
%         nodesid(conn{I}) = 1;
% 
%         while(1)
%             nodes_id_t = nodesid;  
%             r = handles.nodes_id(nodesid==1); 
%             for n=1:length(r)  
%                 nodesid(conn{r(n)}) = 1; 
%             end
%             if sum(nodes_id_t - nodesid)==0
%                 break
%             end
%         end
%         [r,~,~] = find(nodesid==1);
% 
%         id_temp = handles.nodes_id;
%         for n=1:length(r)
%             id_temp(r(n)) = 0;
%         end
% 
%         idnodes_outlet1 = r;
%         idnodes_outlet2 = handles.nodes_id(id_temp~=0);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% 
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
%         plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
% 
%         answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
%         switch answer1
% 
%             case 'Red'
% 
%               handles.id_outlet = idnodes_outlet1;
% 
%             case 'Green'
% 
%               handles.id_outlet = idnodes_outlet2;
% 
%         end
% 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
%         distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
%         id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(handles.id_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

%         figure,
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
%         hold on
%         plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*g')
%         plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
%         quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
%         hold off
%         axis vis3d
%         lighting gouraud
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1

            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

    %         [a,b,c] = size(handles.IPCMRA);
    %     
    %         IM = (handles.SEG/max(handles.SEG(:)))*255;
    %     
    %         node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
    %         cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
    %      
    %         figure,
    %         surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         hold on
    %         surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
    %         patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
    %         patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
    %         plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
    %         quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
    %         hold off
    %         colormap('gray')
    %         daspect([1 1 1])
    %         view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG_for_vwerp.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
    %         imlook3d(handles.SEG)
    %         imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)
        
        
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',7);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end
        
        % we generate the next variables for vwerp
        outlet = handles.PLANE_OUTLET;
        n_outlet = handles.normal_outlet/norm(handles.normal_outlet);
        
        handles.outlet_for_plot = outlet;
        handles.n_outlet_for_plot = n_outlet;

    end

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Sagital view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    if handles.id_view == 1

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            daspect([1 1 1])

            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
    %         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [mean(PUNTOS_C),handles.slider_axes1];
        plot(center(1),center(2),'*g')
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(2),center(1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

    %     center_ori = [center(2),center(1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
    %     point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [point1(1:2),center(3)]; 
        p2       = [center(1:2),center(3) + sqrt(sum((center-point1).^2))]; 
        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles

            Rz = [cosd(n) , -sind(n) , 0 ; sind(n) , cosd(n) , 0 ; 0 , 0 , 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1*Rz;
            VB = (VB/norm(VB))*20;
            VG = v2;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1) , 0 ; sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; 0 , 0 , 1];

        % genero vectores cambiando r y el azul
        VB = v1*Rz;
        VB = VB/norm(VB);
        VR = normal*Rz;
        VR = VR/norm(VR);
        VG = v2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB*20;
        v2_n = VG;
        normal_n = VR*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB_n = v1_n;               
            VB_n = VB_n/norm(VB_n)*20;

            v1_n = (v1_n/norm(v1_n));
            ux = v1_n(1); 
            uy = v1_n(2); 
            uz = v1_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VG_n = (v2_n/norm(v2_n))*Raxis; 
            VG_n = (VG_n/norm(VG_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification

            Area2(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        ux = v1_n(1); 
        uy = v1_n(2); 
        uz = v1_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (v2_n/norm(v2_n))*Raxis; 
        VG_n = (VG_n/norm(VG_n))*20;

        VR_n = normal_n*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        VB_n = cross(VG_n/norm(VG_n),VR_n/norm(VR_n));
        VB_n = VB_n/norm(VB_n)*20;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);     
        handles.faceid = faceid;
        handles.cutpos = cutpos;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_outlet1 = r;
        idnodes_outlet2 = handles.nodes_id(id_temp~=0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
        plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_outlet = idnodes_outlet1;

            case 'Green'

              handles.id_outlet = idnodes_outlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
        distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(id_nodes_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(id_nodes_outlet,1),handles.nodes(id_nodes_outlet,2),handles.nodes(id_nodes_outlet,3),'*g')
    %     plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
    %     quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1

            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
        %  
        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
        %     quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
        %     imlook3d(handles.SEG)
        %     imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity Vectors','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',5);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Axial   view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 2

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [center(3),(handles.slider_axes1-1)*handles.voxel_MR(2),center(2)];
        point1 = [PUNTOS_C(1,2),(handles.slider_axes1-1)*handles.voxel_MR(2),PUNTOS_C(1,1)];
        point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     center_ori = [center(3),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),center(2)];
    %     point1 = [PUNTOS_C(1,2),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
    %     plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
    %     plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d
    %     xlabel('ejex')
    %     ylabel('ejey')
    %     zlabel('ejez')

    %     figure,
    %     imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
    %     hold on
    %     Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
    %     himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
    %     cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %     cdata = double(cdata)*0.5;
    %     set(himage, 'AlphaData', cdata);
    %     plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
    %     center = [handles.slider_axes1,mean(PUNTOS_C)];
    %     plot(center(2),center(3),'*g');
    %     hold off
    %     axis image
    %     set(handles.axes1,'xticklabel',[],'yticklabel',[])
    %     axis off
    %     colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'*c')
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1      = [point1(1), center(2), point1(3)]; 
        p2      = [center(1), center(2) + sqrt(sum((center-p1).^2)), center(3)];

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(p1(1),p1(2),p1(3),'*b')
    %     plot3(p2(1),p2(2),p2(3),'*g')
    %     plot3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),'.c')
    %     quiver3(center(1),center(2),center(3),v1(1),v1(2),v1(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),v2(1),v2(2),v2(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),normal(1),normal(2),normal(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles


    %         Rx = [1 ,       0 ,        0 ; ... 
    %               0 , cosd(n) , -sind(n) ; ...
    %               0 , sind(n) , cosd(n)] ;

    %         Ry = [cosd(n),  0 ,  sind(n) ; ...
    %               0      ,  1 ,        0 ; ...
    %               -sind(n), 0 , cosd(n)] ;

            Rz = [cosd(n) , -sind(n), 0 ; ...
                  sind(n) , cosd(n) , 0 ; ...
                  0       ,        0, 1];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v1;
            VB = (VB/norm(VB))*20;
            VG = v2*Rz;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rz = [cosd(Angle_selected1) , -sind(Angle_selected1), 0 ; ...
              sind(Angle_selected1) , cosd(Angle_selected1) , 0 ; ...
                            0       ,                      0, 1];

        % genero vectores cambiando r y el azul
        VB = v1/norm(v1);
        VB = (VB/norm(VB))*20;
        VG = (v2/norm(v2))*Rz;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rz;
        VR = VR/norm(VR)*20;

        PB = center + VB;  
        PG = center + VG;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
    %     hold on
    %     plot3(center(1),center(2),center(3),'*k')
    %     plot3(PB(1),PB(2),PB(3),'*b')
    %     plot3(PG(1),PG(2),PG(3),'*g')
    %     quiver3(center(1),center(2),center(3),VB(1),VB(2),VB(3),'b','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VG(1),VG(2),VG(3),'g','Linewidth',3)
    %     quiver3(center(1),center(2),center(3),VR(1),VR(2),VR(3),'r','Linewidth',3)
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

    %     figure, % julio
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k') % julio
    %     hold on % julio

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

    %         plot3(cutat(2,1),cutat(2,2),cutat(2,3),'*k')
    %         plot3(cutat(1,1),cutat(1,2),cutat(1,3),'*b')
    %         plot3(cutat(3,1),cutat(3,2),cutat(3,3),'*g')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end
    %     hold off
    %     daspect([1 1 1])
    %     axis vis3d


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_outlet1 = r;
        idnodes_outlet2 = handles.nodes_id(id_temp~=0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
        plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_outlet = idnodes_outlet1;

            case 'Green'

              handles.id_outlet = idnodes_outlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
        distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(id_nodes_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(id_nodes_outlet,1),handles.nodes(id_nodes_outlet,2),handles.nodes(id_nodes_outlet,3),'*g')
    %     plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
    %     quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
    
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
        %  
        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
        %     quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
        %     imlook3d(handles.SEG)
        %     imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity Vectors','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',5);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %________________ Coronal view _______________________________________%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif handles.id_view == 3

        % show the image
        if handles.id_image == 1

            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id IPCMRA =1
            handles.id_ipcmra = 1;
            handles.id_mag = 0;

        elseif handles.id_image == 2

            % show the axes1
            set(handles.axes1,'Visible','on');
            axis off

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            % id MAG =1
            handles.id_ipcmra = 0;
            handles.id_mag = 1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = imline(handles.axes1);
        PUNTOS_C = wait(h);

        if handles.id_image ==1

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        elseif handles.id_image == 2

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
    %         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
    %         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
    %         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
    %         cdata = double(cdata)*0.5;
    %         set(himage, 'AlphaData', cdata);

        end

        plot(PUNTOS_C(:,1),PUNTOS_C(:,2),'r','Linewidth',3)
        center = [handles.slider_axes1,mean(PUNTOS_C)];
        plot(center(2),center(3),'*g');
        hold off
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        pause(0.5)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % wait bar for plane 

        h = waitbar(0,'Please wait...');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initial variables of the plane

        center_ori = [(handles.slider_axes1-1)*handles.voxel_MR(1),center(3),center(2)];
        point1 = [(handles.slider_axes1-1)*handles.voxel_MR(1),PUNTOS_C(1,2),PUNTOS_C(1,1)];
        point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

    %     center_ori = [mean(mean(handles.yd(handles.slider_axes1,:,:))),center(3),center(2)];
    %     point1 = [mean(mean(handles.yd(handles.slider_axes1,:,:))),PUNTOS_C(1,2),PUNTOS_C(1,1)];
    %     point2 = [center_ori(1) + sqrt(sum((center_ori-point1).^2)),center_ori(2), center_ori(3) ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cutat = [point1;center_ori;point2]; % tres puntos para generar el plano
        [~,~,~,elemid]  = qmeshcut(handles.elem,handles.nodes,ones(size(handles.nodes(:,1))),cutat);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(sqrt(sum((handles.nodes - repmat(center_ori,size(handles.nodes,1),1)).^2,2)));% julio 
        [conn,~,~] = meshconn(handles.elem(elemid,:),size(handles.nodes,1));    
        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);
        handles.id_section = r;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % center of the vessel
        center = mean(handles.nodes(handles.id_section,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p1       = [center(1), point1(2:3)]; 
        p2       = [center(1) + sqrt(sum((center-point1).^2)), center(2:3)]; 

        v1       = p1 - center;
        v2       = p2 - center;
        normal   = cross(v1,v2);

        % normalization
        v1      = v1/norm(v1)*20;
        v2      = v2/norm(v2)*20;
        normal  = normal/norm(normal)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area        = zeros(2,size(Angles,2));
        Area(1,:)   = Angles;
        cont        = 1;

        for n = Angles


            Rx = [1 , 0 , 0 ; 0 , cosd(n) , -sind(n) ; 0 , sind(n), cosd(n)];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VB = v2;
            VB = (VB/norm(VB))*20;
            VG = v1*Rx;
            VG = (VG/norm(VG))*20;
            PB = center + VB;  
            PG = center + VG;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB;center;PG]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area(2,cont) = pi*rad^2;
            cont = cont+1;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(40/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C,I] = min(Area(2,:));
        % angle selected
        Angle_selected1 = Area(1,I);
        Rx = [1 , 0 , 0 ; 0 , cosd(Angle_selected1) , -sind(Angle_selected1) ; 0 , sind(Angle_selected1), cosd(Angle_selected1)];

        % genero vectores cambiando r y el azul
        VB = v2/norm(v2);
        VB = (VB/norm(VB))*20;
        VG = (v1/norm(v1))*Rx;
        VG = (VG/norm(VG))*20;
        VR = (normal/norm(normal))*Rx;
        VR = VR/norm(VR)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %%% segunda rotacion
        v1_n = VB;
        v2_n = VG;
        normal_n = VR;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Angles      = -60:1:60;
        Area2 = zeros(2,size(Angles,2));
        Area2(1,:) = Angles;
        cont = 1;

        for n = Angles

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % genero vectores cambiando r y el azul
            VG_n = v2_n;               
            VG_n = VG_n/norm(VG_n)*20;

            v2_n = (v2_n/norm(v2_n));
            ux = v2_n(1); 
            uy = v2_n(2); 
            uz = v2_n(3);

            Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                      uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                      uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            VB_n = (v1_n/norm(v1_n))*Raxis; 
            VB_n = (VB_n/norm(VB_n))*20;

            PB_n = center + VB_n; 
            PG_n = center + VG_n;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cutat = [PB_n;center;PG_n]; % tres puntos para generar el plano

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [~,~,~,facesid]  = qmeshcut(handles.faces,handles.nodes,handles.nodes_id,cutat);% julio
            nodid = handles.faces(facesid,:); 
            nodid = unique(nodid(:));
            [~,I] = min(sqrt(sum((handles.nodes(nodid,:) - repmat(center_ori,size(handles.nodes(nodid,:),1),1)).^2,2)));% julio 

            [conn,~,~] = meshconn(handles.faces(facesid,:),size(handles.nodes,1));    % julio        
            nodes_id = zeros(size(handles.nodes,1),1);
            nodes_id(conn{nodid(I)}) = 1;

            while(1)
                nodes_id_t = nodes_id;  
                [r,~,~] = find(nodes_id==1); 
                for m=1:length(r)  
                    nodes_id(conn{r(m)}) = 1; 
                end
                if sum(nodes_id_t - nodes_id)==0
                    break
                end
            end

            [r,~,~] = find(nodes_id==1);
            handles.id_section = r;

            rad = mean(sqrt(sum((handles.nodes(handles.id_section,:) - repmat(center,length(handles.id_section),1)).^2,2)));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % area quantification
            Area2(2,cont) = pi*rad^2;
            cont = cont+1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(80/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~,I] = min(Area2(2,:));
        Angle_selected1 = Area2(1,I);
        n = Angle_selected1;

        v2_n = (v2_n/norm(v2_n));
        ux = v2_n(1); 
        uy = v2_n(2); 
        uz = v2_n(3);

        Raxis = [ cosd(n)+(ux*ux)*(1-cosd(n))  ux*uy*(1-cosd(n))-uz*sind(n) ux*uz*(1-cosd(n))+uy*sind(n)
                  uy*ux*(1-cosd(n))+uz*sind(n) cosd(n)+(uy*uy)*(1-cosd(n))  uy*uz*(1-cosd(n))-ux*sind(n)
                  uz*ux*(1-cosd(n))-uy*sind(n) uz*uy*(1-cosd(n))+ux*sind(n) cosd(n)+(uz*uz)*(1-cosd(n))];

        VG_n = (VG_n/norm(VG_n))*20;

        VB_n = (v1_n/norm(v1_n))*Raxis; 
        VB_n = (VB_n/norm(VB_n))*20;

        VR_n = (normal_n/norm(normal_n))*Raxis; 
        VR_n = VR_n/norm(VR_n)*20;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        PB = center + VB_n;  
        PG = center + VG_n;
        cutat = [PB;center;PG]; % tres puntos para generar el plano
        [cutpos,~,faceid,elemid]  = qmeshcut(handles.elem,handles.nodes,handles.nodes_id,cutat);  
        handles.cutpos = cutpos;
        handles.faceid = faceid;
        Selected_element = handles.elem(elemid,:);
        area = Area2(2,I);
        [conn,~,~]=meshconn(Selected_element,size(handles.nodes,1));

        negative_vector = [0,1,0];
        negative_vector = negative_vector/norm(negative_vector)*20;

        [~,I] = min(sqrt(sum( (repmat(center,size(handles.nodes,1),1) - handles.nodes).^2,2 )));

        nodes_id = zeros(size(handles.nodes,1),1);
        nodes_id(conn{I}) = 1;

        while(1)
            nodes_id_t = nodes_id;  
            [r,~,~] = find(nodes_id==1); 
            for n=1:length(r)  
                nodes_id(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodes_id)==0
                break
            end
        end

        [r,~,~] = find(nodes_id==1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SE = Selected_element(:);
        for n=1:length(r)
            EE = abs(Selected_element(:) - r(n));
            SE(EE==0) = 0;
        end
        Selected_element_new = reshape(SE,[size(Selected_element,1), size(Selected_element,2)]);
        elemid(sum(Selected_element_new,2)~=0) = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elem_new = handles.elem;
        elem_new(elemid,:) = []; 

        [conn,~,~] = meshconn(elem_new,size(handles.nodes,1));
        nodesid = zeros(size(handles.nodes,1),1);
        nodesid(conn{I}) = 1;

        while(1)
            nodes_id_t = nodesid;  
            r = handles.nodes_id(nodesid==1); 
            for n=1:length(r)  
                nodesid(conn{r(n)}) = 1; 
            end
            if sum(nodes_id_t - nodesid)==0
                break
            end
        end
        [r,~,~] = find(nodesid==1);

        id_temp = handles.nodes_id;
        for n=1:length(r)
            id_temp(r(n)) = 0;
        end

        idnodes_outlet1 = r;
        idnodes_outlet2 = handles.nodes_id(id_temp~=0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close(h)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(idnodes_outlet1,1),handles.nodes(idnodes_outlet1,2),handles.nodes(idnodes_outlet1,3),'*r')
        plot3(handles.nodes(idnodes_outlet2,1),handles.nodes(idnodes_outlet2,2),handles.nodes(idnodes_outlet2,3),'*g')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        answer1 = questdlg('¿which color do you want to select as inlet?','Question','Red','Green','Green');
        switch answer1

            case 'Red'

              handles.id_outlet = idnodes_outlet1;

            case 'Green'

              handles.id_outlet = idnodes_outlet2;

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % normal outlet
        normal_out = VR_n; % normal of the slice 
        center_out = center; % center of the section inlet
        distance_out = sqrt(sum((handles.nodes - repmat(center,size(handles.nodes,1),1)).^2,2));
        id_nodes_outlet = handles.id_outlet(distance_out(handles.id_outlet)<10);
        center_out_n = mean(handles.nodes(id_nodes_outlet,:)); % center of the nurb of point inlet
        normal_out_n = (center_out_n - center_out)/norm(center_out_n - center_out)*20;

        angle_b_normal = acos(sum(normal_out_n.*normal_out)/(20*20))*180/pi;

        if angle_b_normal < 90
            handles.normal_outlet = (normal_out/norm(normal_out));
        else
            handles.normal_outlet = (normal_out/norm(normal_out))*-1;
        end

    %     figure(1)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    %     hold on
    %     plot3(handles.nodes(id_nodes_outlet,1),handles.nodes(id_nodes_outlet,2),handles.nodes(id_nodes_outlet,3),'*g')
    %     plot3(center_out(1),center_out(2),center_out(3),'*b','Linewidth',5)
    %     quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
    %     hold off
    %     axis vis3d
    %     lighting gouraud
    %     daspect([1,1,1])
    %     axis off
    %     view([-34,-51])


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________ JULIO SOTELO ___________________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if handles.id_vwerp == 1
            
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,2)-1, 0:size(handles.IPCMRA,1)-1, 0:size(handles.IPCMRA,3)-1);
            X = X*handles.voxel_MR(1);
            Y = Y*handles.voxel_MR(2);
            Z = Z*handles.voxel_MR(3);

        %     [a,b,c] = size(handles.IPCMRA);
        % 
        %     IM = (handles.SEG/max(handles.SEG(:)))*255;
        % 
        %     node_n = [handles.nodes(:,2),handles.nodes(:,1),handles.nodes(:,3)];
        %     cut_p_node_n = [handles.cutpos(:,2),handles.cutpos(:,1),handles.cutpos(:,3)];
        %  
        %     figure,
        %     surface(X(:,:,round(c/2)),Y(:,:,round(c/2)),Z(:,:,round(c/2)),IM(:,:,round(c/2)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     hold on
        %     surface(squeeze(X(:,round(b/2),:)),squeeze(Y(:,round(b/2),:)),squeeze(Z(:,round(b/2),:)),squeeze(IM(:,round(b/2),:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     surface(squeeze(X(round(a/2),:,:)),squeeze(Y(round(a/2),:,:)),squeeze(Z(round(a/2),:,:)),squeeze(IM(round(a/2),:,:)),'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct')
        %     patch('faces',handles.faces,'vertices',node_n,'EdgeColor','b','FaceColor','none')
        %     patch('faces',handles.faceid,'vertices',cut_p_node_n,'EdgeColor','k','FaceColor','r')
        %     plot3(center_out(2),center_out(1),center_out(3),'*c','Linewidth',3)
        %     quiver3(center_out(2),center_out(1),center_out(3),handles.normal_outlet(2),handles.normal_outlet(1),handles.normal_outlet(3),30,'g','Linewidth',3)
        %     hold off
        %     colormap('gray')
        %     daspect([1 1 1])
        %     view(3)

            a1 = handles.normal_outlet(2);
            b1 = handles.normal_outlet(1);
            c1 = handles.normal_outlet(3);
            OUT = a1*X + b1*Y + c1*Z - (a1*center_out(2) + b1*center_out(1) + c1*center_out(3));

            SE = strel('sphere', 1);
            BW2 = imdilate(double(OUT<=0),SE);
            PLANE = handles.SEG.*(double(BW2) - double(OUT<=0));

            X_C = X.*PLANE;
            Y_C = Y.*PLANE;
            Z_C = Z.*PLANE;
            DD = sqrt((X_C-center_out(2)).^2 + (Y_C-center_out(1)).^2 + (Z_C-center_out(3)).^2);

            SEGMENT = DD.*PLANE;
            SEGMENT = SEGMENT(:);
            POINT = reshape(double(SEGMENT == min(SEGMENT(SEGMENT~=0))),[size(DD,2),size(DD,1),size(DD,3)]);
            [LL,~] = bwlabeln(PLANE,18);

            handles.PLANE_OUTLET = double(LL == max(POINT(:).*LL(:)));
        %     imlook3d(handles.SEG)
        %     imlook3d(handles.PLANE_OUTLET)

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ___________________ JULIO SOTELO ___________________________________%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
    %     patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
        hold on
        plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity Vectors','Mesh + Inlet','Mesh + Outlet'};
        set(handles.popupmenu2,'String',handles.list_string,'Value',5);

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 1;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        handles.outlet_exec = 1;
        set(handles.pushbutton1,'Visible','off');
        set(handles.pushbutton2,'Visible','off');
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu3,'Value',1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.inlet_exec == 1

            set(handles.pushbutton3,'Visible','on');

        else

            set(handles.pushbutton3,'Visible','off');

        end

    end  
end


handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.popupmenu2,'Value')

    case 1
        
        set(handles.uipanel2,'BackgroundColor',[0,0,0])
        cla(handles.axes1,'reset');
        set(handles.axes1,'visible','off')
        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_ipcmra = 0;
        handles.id_mag = 0;


    case 2
        
        set(handles.uipanel2,'BackgroundColor',[0,0,0])
        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','r','edgecolor','k')
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        
        handles.id_mesh = 1;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

     case 3
        
            
            set(handles.uipanel2,'BackgroundColor',[0,0,0])
            set(handles.popupmenu1,'Value',1)

            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            set(handles.pushbutton5,'Visible','off')
            set(handles.pushbutton6,'Visible','off')

            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'jet');
            [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vel = min(handles.mags_vel(:));
            handles.max_vel = max(handles.mags_vel(:));
            handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vel handles.max_vel];
            c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
            c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
            c.Color = [1 1 1]; % color
            c.Location = 'manual';
            c.Position = [0.2 0.2 0.02 0.3];
            c.FontWeight = 'bold';
            c.Label.String = 'Velocity [m/s]';
            windows_screen_size = get(0,'ScreenSize');
            if windows_screen_size(4)<=900
                c.FontSize = 9;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.FontSize = 10;
            else
                c.FontSize = 11;
            end
            if windows_screen_size(4)<=900
                c.Label.FontSize = 10;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.Label.FontSize = 12;
            else
                c.Label.FontSize = 14;
            end
            c.Label.FontWeight = 'bold';
            c.Label.Color = [1 1 1];% color
            caxis(handles.axes1, [handles.min_vel handles.max_vel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])


            % id mesh
            handles.id_mesh = 0;
            handles.id_mesh_vel = 1;

            % slider adjustment
            slider_step(1) = 1/size(handles.veset,3);
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
            set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

            % id mesh
            handles.id_mesh = 0;
            handles.id_mesh_vel = 1;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;

        
    case 4
        
        if handles.id_vwerp ==1

            % adjusting values for plot
            mask = handles.SEG_for_vwerp;

            % adjusting the velocities
                v = [];
                if max(abs(handles.MR_PCA_FH(:)))>10
                    v{3}.im = (handles.MR_PCA_RL.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_RL,4))/100)*-1;
                    v{2}.im = (handles.MR_PCA_FH.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_FH,4))/100)*-1;
                    v{1}.im = (handles.MR_PCA_AP.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_AP,4))/100);
                else
                    v{3}.im = (handles.MR_PCA_RL.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_RL,4)))*-1;
                    v{2}.im = (handles.MR_PCA_FH.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_FH,4)))*-1;
                    v{1}.im = (handles.MR_PCA_AP.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_AP,4)));
                end

            handles.v_for_plot = v;
            handles.mask_for_plot = mask;

            % we extract the voxelization using binsurface iso2mesh
            [node_bin,elem_bin] = binsurface(handles.mask_for_plot);

            
            axes(handles.axes1);
            plot(0.0)
            set(handles.uipanel2,'BackgroundColor',[1,1,1])
            
            
            %visu_mask_inlet_outlet(handles.v_for_plot, handles.mask_for_plot, handles.inlet_for_plot, handles.outlet_for_plot, handles.n_inlet_for_plot, handles.n_outlet_for_plot)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % visualization of the velocities vector field %%%%%%%%%%%%%%%%
            % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            % code
            time = 5;
            n = sum(sum(sum(handles.mask_for_plot)));
            pt = zeros(n,3);
            V = zeros(n,3);

            n = 0;
            for i = 1:size(handles.mask_for_plot,1)
               for j = 1:size(handles.mask_for_plot,2)
                   for k = 1:size(handles.mask_for_plot,3)
                       if(handles.mask_for_plot(i,j,k) ~= 1)
                           continue
                       end

                       n = n + 1;
                       pt(n,:) = [i, j, k];
                       V(n,:) = [handles.v_for_plot{1}.im(i,j,k,time), handles.v_for_plot{2}.im(i,j,k,time), handles.v_for_plot{3}.im(i,j,k,time)];
                   end
               end
            end

            vs = 10;
            jump = 1;
            q = quiver3(pt(1:jump:end,1),pt(1:jump:end,2),pt(1:jump:end,3),vs*V(1:jump:end,1),vs*V(1:jump:end,2),vs*V(1:jump:end,3), 'Linewidth',1);
            mag_vel = sqrt(sum(V(1:jump:end,:).^2,2));
            currentColormap = colormap(handles.axes1,'jet');
            [~, ~, ind] = histcounts(mag_vel(:), size(currentColormap, 1));
            handles.cmap = uint8(ind2rgb(ind, currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vel = min(mag_vel(:));
            handles.max_vel = max(mag_vel(:));
            handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vel handles.max_vel];
            c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
            c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
            c.Color = [0 0 0]; % color
            c.FontWeight = 'bold';
            c.Label.String = 'Velocity [m/s]';
            c.FontSize = 11;
            c.Label.FontSize = 14;
            c.Label.FontWeight = 'bold';
            c.Label.Color = [0 0 0];% color
            caxis(handles.axes1, [handles.min_vel handles.max_vel]);
            axis equal
            q.AutoScale='off';

            hold on
            
%             if handles.inlet_exec == 1 
%             
%             handles.n_inlet_for_plot = handles.n_inlet_for_plot*-1;
%             [r,c,v] = ind2sub(size((handles.inlet_for_plot)),find((handles.inlet_for_plot) == 1));
%             Xo(1) = mean(r); Xo(2) = mean(c); Xo(3) = mean(v); 
%             scatter3(Xo(1),Xo(2),Xo(3),'ro')
%             quiver3(Xo(1),Xo(2),Xo(3),vs*handles.n_inlet_for_plot(1),vs*handles.n_inlet_for_plot(2),vs*handles.n_inlet_for_plot(3),2,'r', 'linewidth',6)
% 
%             end
%             
%             if handles.outlet_exec == 1
%             [r,c,v] = ind2sub(size((handles.outlet_for_plot)),find((handles.outlet_for_plot) == 1));
%             Xo(1) = mean(r); Xo(2) = mean(c); Xo(3) = mean(v); 
%             scatter3(Xo(1),Xo(2),Xo(3),'go')
%             quiver3(Xo(1),Xo(2),Xo(3),vs*handles.n_outlet_for_plot(1),vs*handles.n_outlet_for_plot(2),vs*handles.n_outlet_for_plot(3),2,'g', 'linewidth',6)
%             
%             end
            
            % binsurface plot
            trisurf(elem_bin,node_bin(:,1)+0.5,node_bin(:,2)+0.5,node_bin(:,3)+0.5, 'facecolor','none');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % visualization of the velocities vector field %%%%%%%%%%%%%%%%
            % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            xlabel('A-P', 'Fontsize', 14,'FontWeight','bold')
            ylabel('F-H', 'Fontsize', 14,'FontWeight','bold')
            zlabel('R-L', 'Fontsize', 14,'FontWeight','bold')
            title('Vector Field','Fontsize', 16,'FontWeight','bold')
            hold off
            grid on
            daspect([1,1,1])
            view(3)

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            
            if (handles.inlet_exec == 1) && (handles.outlet_exec == 1)

                set(handles.pushbutton3,'Visible','on');

            else 

                set(handles.pushbutton3,'Visible','off');

            end  
            
%             set(handles.pushbutton3,'Visible','off')
            set(handles.pushbutton5,'Visible','off')
            set(handles.pushbutton6,'Visible','off')

        else
        
            set(handles.popupmenu1,'Value',1)
            
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])


            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 1;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;

        end
    case 5
        
        if handles.id_vwerp ==1
            
            % adjusting values for plot
            mask = handles.SEG_for_vwerp;
            handles.mask_for_plot = mask;
            vs = 10;
            % we extract the voxelization using binsurface iso2mesh
            [node_bin,elem_bin] = binsurface(handles.mask_for_plot);

            
            axes(handles.axes1);
            plot(0.0)
            set(handles.uipanel2,'BackgroundColor',[1,1,1])
           
            % binsurface plot
            trisurf(elem_bin,node_bin(:,1)+0.5,node_bin(:,2)+0.5,node_bin(:,3)+0.5, 'facecolor','none', 'edgealpha',0.5);
            hold on
            
            if handles.inlet_exec == 1 
                
                [node_inletp,elem_inletp] = binsurface(handles.PLANE_INLET);
                trisurf(elem_inletp,node_inletp(:,1)+0.5,node_inletp(:,2)+0.5,node_inletp(:,3)+0.5, 'facecolor','r');
            
                handles.n_inlet_for_plot = handles.normal_inlet;
                
                [r,c,v] = ind2sub(size((handles.inlet_for_plot)),find((handles.inlet_for_plot) == 1));
                Xo(1) = mean(r); Xo(2) = mean(c); Xo(3) = mean(v); 
                scatter3(Xo(1),Xo(2),Xo(3),'ro')
                quiver3(Xo(1),Xo(2),Xo(3),vs*handles.n_inlet_for_plot(1),vs*handles.n_inlet_for_plot(2),vs*handles.n_inlet_for_plot(3),2,'r', 'linewidth',6)

            end
            
            if handles.outlet_exec == 1
                
                [node_outletp,elem_outletp] = binsurface(handles.PLANE_OUTLET);
                trisurf(elem_outletp,node_outletp(:,1)+0.5,node_outletp(:,2)+0.5,node_outletp(:,3)+0.5, 'facecolor','g');

                handles.n_outlet_for_plot = handles.normal_outlet;
                
                [r,c,v] = ind2sub(size((handles.outlet_for_plot)),find((handles.outlet_for_plot) == 1));
                Xo(1) = mean(r); Xo(2) = mean(c); Xo(3) = mean(v); 
                scatter3(Xo(1),Xo(2),Xo(3),'go')
                quiver3(Xo(1),Xo(2),Xo(3),vs*handles.n_outlet_for_plot(1),vs*handles.n_outlet_for_plot(2),vs*handles.n_outlet_for_plot(3),2,'g', 'linewidth',6)

            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % visualization of the velocities vector field %%%%%%%%%%%%%%%%
            % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            xlabel('A-P', 'Fontsize', 14,'FontWeight','bold')
            ylabel('F-H', 'Fontsize', 14,'FontWeight','bold')
            zlabel('R-L', 'Fontsize', 14,'FontWeight','bold')
            title('Vector Field','Fontsize', 16,'FontWeight','bold')
            hold off
            grid on
            daspect([1,1,1])
            view(3)

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            if (handles.inlet_exec == 1) && (handles.outlet_exec == 1)

                set(handles.pushbutton3,'Visible','on');

            else 

                set(handles.pushbutton3,'Visible','off');

            end  
            
%             set(handles.pushbutton3,'Visible','off')
            set(handles.pushbutton5,'Visible','off')
            set(handles.pushbutton6,'Visible','off')
            
        else
            set(handles.popupmenu1,'Value',1)

            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 1;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;
        end
    case 6

        if handles.id_vwerp ==1

            set(handles.uipanel2,'BackgroundColor',[0,0,0])
            set(handles.popupmenu1,'Value',1)
            set(handles.pushbutton5,'Visible','off')
            set(handles.pushbutton6,'Visible','on')
            
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
            center_in = mean(handles.nodes(handles.id_inlet,1:3));
            quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])


            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            
            if (handles.inlet_exec == 1) && (handles.outlet_exec == 1)

                set(handles.pushbutton3,'Visible','on');

            else 

                set(handles.pushbutton3,'Visible','off');

            end  
            
%             set(handles.pushbutton3,'Visible','off')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 1;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;


        else

            set(handles.popupmenu1,'Value',1)
    
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace,'CDataMapping','Scaled')
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_lap = min(handles.Laplace(:));
            handles.max_lap = max(handles.Laplace(:));
            handles.mean_lap = (handles.min_lap + handles.max_lap)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_lap handles.max_lap];
            c.Ticks = [handles.min_lap, (handles.min_lap + handles.mean_lap)/2, handles.mean_lap, (handles.max_lap + handles.mean_lap)/2, handles.max_lap];
            c.TickLabels = {num2str(handles.min_lap,'%0.2f'), num2str((handles.min_lap + handles.mean_lap)/2,'%0.2f'), num2str(handles.mean_lap,'%0.2f'), num2str((handles.max_lap + handles.mean_lap)/2,'%0.2f'), num2str(handles.max_lap,'%0.2f')};
            c.Color = [1 1 1]; % color
            c.Location = 'manual';
            c.Position = [0.2 0.2 0.02 0.3];
            c.FontWeight = 'bold';
            c.Label.String = 'Laplace Value [-]';
            windows_screen_size = get(0,'ScreenSize');
            if windows_screen_size(4)<=900
                c.FontSize = 9;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.FontSize = 10;
            else
                c.FontSize = 11;
            end
            if windows_screen_size(4)<=900
                c.Label.FontSize = 10;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.Label.FontSize = 12;
            else
                c.Label.FontSize = 14;
            end
            c.Label.FontWeight = 'bold';
            c.Label.Color = [1 1 1];% color
            caxis(handles.axes1, [handles.min_lap handles.max_lap]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
    
            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 1;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;

        end

    case 7 % centerline

        if handles.id_vwerp ==1

            set(handles.uipanel2,'BackgroundColor',[0,0,0])
            set(handles.popupmenu1,'Value',1)
            set(handles.pushbutton5,'Visible','on')
            set(handles.pushbutton6,'Visible','off')
            
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
            hold on
            plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
            center_out = mean(handles.nodes(handles.id_outlet,1:3));
            quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            
            if (handles.inlet_exec == 1) && (handles.outlet_exec == 1)

                set(handles.pushbutton3,'Visible','on');

            else 

                set(handles.pushbutton3,'Visible','off');

            end  
            
%             set(handles.pushbutton3,'Visible','off')

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 1;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;

        else

            set(handles.popupmenu1,'Value',1)
    
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            plot3(handles.centerline(:,1),handles.centerline(:,2),handles.centerline(:,3),'-c','LineWidth',3)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])
    
            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 1;
            handles.id_diameter = 0;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;

        end

    case 8 % diameter
        
        if handles.id_vwerp ==1
            
            set(handles.uipanel2,'BackgroundColor',[1,1,1])
        
            axes(handles.axes1);
            plot(handles.time, [handles.DP_VWERP(end), handles.DP_VWERP(1:end-1)]*handles.pa2mmhg,'b--*','LineWidth',2, 'MarkerSize', 10)
            hold on
            plot(handles.time, zeros(size(handles.time)),'k-')
            hold off
            xlabel('time [s]', 'Fontsize', 14,'FontWeight','bold')
            ylabel('Relative pressure drop [mmHg]', 'Fontsize', 14,'FontWeight','bold')
            grid on
            title(['Estimated relative pressure drops, WERP'], 'Fontsize', 16,'FontWeight','bold');
        
            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')
            if (handles.inlet_exec == 1) && (handles.outlet_exec == 1)

                set(handles.pushbutton3,'Visible','on');

            else 

                set(handles.pushbutton3,'Visible','off');

            end  
%             set(handles.pushbutton3,'Visible','off')
            set(handles.pushbutton5,'Visible','off')
            set(handles.pushbutton6,'Visible','off')
            
        else

            set(handles.popupmenu1,'Value',1)

            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.diameter,'CDataMapping','Scaled')
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            handles.min_dia = min(handles.diameter(:));
            handles.max_dia = max(handles.diameter(:));
            handles.mean_dia = (handles.min_dia + handles.max_dia)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_dia handles.max_dia];
            c.Ticks = [handles.min_dia, (handles.min_dia + handles.mean_dia)/2, handles.mean_dia, (handles.max_dia + handles.mean_dia)/2, handles.max_dia];
            c.TickLabels = {num2str(handles.min_dia,'%0.2f'), num2str((handles.min_dia + handles.mean_dia)/2,'%0.2f'), num2str(handles.mean_dia,'%0.2f'), num2str((handles.max_dia + handles.mean_dia)/2,'%0.2f'), num2str(handles.max_dia,'%0.2f')};
            c.Color = [1 1 1]; % color
            c.Location = 'manual';
            c.Position = [0.2 0.2 0.02 0.3];
            c.FontWeight = 'bold';
            c.Label.String = 'Diameter [cm]';
            windows_screen_size = get(0,'ScreenSize');
            if windows_screen_size(4)<=900
                c.FontSize = 9;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.FontSize = 10;
            else
                c.FontSize = 11;
            end
            if windows_screen_size(4)<=900
                c.Label.FontSize = 10;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.Label.FontSize = 12;
            else
                c.Label.FontSize = 14;
            end
            c.Label.FontWeight = 'bold';
            c.Label.Color = [1 1 1];% color
            caxis(handles.axes1, [handles.min_dia handles.max_dia]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.slider1,'visible','off')
            set(handles.text1,'visible','off')
            set(handles.popupmenu1,'value',1)
            set(handles.popupmenu3,'value',1)
            set(handles.pushbutton1,'visible','off')
            set(handles.pushbutton2,'visible','off')


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handles.id_mesh = 0;
            handles.id_mesh_vel = 0;
            handles.id_mesh_inlet = 0;
            handles.id_mesh_outlet = 0;
            handles.id_mesh_laplace = 0;
            handles.id_centerline = 0;
            handles.id_diameter = 1;
            handles.id_radius = 0;
            handles.id_axial_unit_vectors = 0;
            handles.id_circumferential_unit_vectors = 0;
            handles.id_wss_a = 0;
            handles.id_wss_c = 0;
            handles.id_axial_angle = 0;
            handles.id_forward_vel = 0;
            handles.id_backward_vel = 0;
            handles.id_regurgitant_flow = 0;
            handles.id_centerline_flow = 0;
            handles.id_eccentricity = 0;
            handles.id_curvature = 0; % new data Julio Sotelo
            handles.id_flattening = 0; % new data Julio Sotelo
            handles.id_ellipticity = 0; % new data Julio Sotelo
            handles.id_length_vessel = 0; % new data Julio Sotelo
    %         handles.id_circulation = 0; % new data Julio Sotelo
            handles.id_forward_vortex = 0; % new data Julio Sotelo
            handles.id_ipcmra = 0;
            handles.id_mag = 0;
        
        end

    case 9 % radius

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.radius,'CDataMapping','Scaled')
        colormap(handles.axes1,'parula');
        c = colorbar(handles.axes1);
        handles.min_rad = min(handles.radius(:));
        handles.max_rad = max(handles.radius(:));
        handles.mean_rad = (handles.min_rad + handles.max_rad)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_rad handles.max_rad];
        c.Ticks = [handles.min_rad, (handles.min_rad + handles.mean_rad)/2, handles.mean_rad, (handles.max_rad + handles.mean_rad)/2, handles.max_rad];
        c.TickLabels = {num2str(handles.min_rad,'%0.2f'), num2str((handles.min_rad + handles.mean_rad)/2,'%0.2f'), num2str(handles.mean_rad,'%0.2f'), num2str((handles.max_rad + handles.mean_rad)/2,'%0.2f'), num2str(handles.max_rad,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Radius [cm]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_rad handles.max_rad]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 1;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 10 % axial unit vector


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.axial_unit_vectors(:,1),handles.axial_unit_vectors(:,2),handles.axial_unit_vectors(:,3),1,'g','LineWidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 1;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;


    case 11 % circumferential unit vector

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.circumferential_unit_vectors(:,1),handles.circumferential_unit_vectors(:,2),handles.circumferential_unit_vectors(:,3),1,'r','LineWidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 1;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 12 % wall shear stress axial


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_A(1:end,1,handles.peak_flow),handles.WSS_A(1:end,2,handles.peak_flow),handles.WSS_A(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.WSS_A,1) size(handles.WSS_A,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_wssa = min(handles.mag_WSS_A(:));
        handles.max_wssa = max(handles.mag_WSS_A(:));
        handles.mean_wssa = (handles.min_wssa + handles.max_wssa)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_wssa handles.max_wssa];
        c.Ticks = [handles.min_wssa, (handles.min_wssa + handles.mean_wssa)/2, handles.mean_wssa, (handles.max_wssa + handles.mean_wssa)/2, handles.max_wssa];
        c.TickLabels = {num2str(handles.min_wssa,'%0.2f'), num2str((handles.min_wssa + handles.mean_wssa)/2,'%0.2f'), num2str(handles.mean_wssa,'%0.2f'), num2str((handles.max_wssa + handles.mean_wssa)/2,'%0.2f'), num2str(handles.max_wssa,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'WSS-A [N/m^{2}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        % slider adjustment
        slider_step(1) = 1/size(handles.WSS_A,3);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.WSS_A,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 1;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 13 % wall shear stress circumferential


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.WSS_C,1) size(handles.WSS_C,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_wssc = min(handles.mag_WSS_C(:));
        handles.max_wssc = max(handles.mag_WSS_C(:));
        handles.mean_wssc = (handles.min_wssc + handles.max_wssc)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_wssc handles.max_wssc];
        c.Ticks = [handles.min_wssc, (handles.min_wssc + handles.mean_wssc)/2, handles.mean_wssc, (handles.max_wssc + handles.mean_wssc)/2, handles.max_wssc];
        c.TickLabels = {num2str(handles.min_wssc,'%0.2f'), num2str((handles.min_wssc + handles.mean_wssc)/2,'%0.2f'), num2str(handles.mean_wssc,'%0.2f'), num2str((handles.max_wssc + handles.mean_wssc)/2,'%0.2f'), num2str(handles.max_wssc,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'WSS-C [N/m^{2}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        % slider adjustment
        slider_step(1) = 1/size(handles.WSS_C,3);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.WSS_C,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 1;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 14 % axial angle

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        c = colorbar(handles.axes1);
        handles.min_axa = min(handles.angle_axial_direction(:));
        handles.max_axa = max(handles.angle_axial_direction(:));
        handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_axa handles.max_axa];
        c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
        c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Axial Angle [^{O}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_axa handles.max_axa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        % slider adjustment
        slider_step(1) = 1/size(handles.angle_axial_direction,3);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.angle_axial_direction,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 1;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 15 % forward velocity

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,handles.peak_flow),handles.forward_velocity(1:end,2,handles.peak_flow),handles.forward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_fvel = min(handles.mag_forward_velocity(:));
        handles.max_fvel = max(handles.mag_forward_velocity(:));
        handles.mean_fvel = (handles.min_fvel + handles.max_fvel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_fvel handles.max_fvel];
        c.Ticks = [handles.min_fvel, (handles.min_fvel + handles.mean_fvel)/2, handles.mean_fvel, (handles.max_fvel + handles.mean_fvel)/2, handles.max_fvel];
        c.TickLabels = {num2str(handles.min_fvel,'%0.2f'), num2str((handles.min_fvel + handles.mean_fvel)/2,'%0.2f'), num2str(handles.mean_fvel,'%0.2f'), num2str((handles.max_fvel + handles.mean_fvel)/2,'%0.2f'), num2str(handles.max_fvel,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Forward Velocity [m/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        % slider adjustment
        slider_step(1) = 1/size(handles.forward_velocity,3);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.forward_velocity,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 1;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 16 % backward velocity


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,handles.peak_flow),handles.backward_velocity(1:end,2,handles.peak_flow),handles.backward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_bvel = min(handles.mag_backward_velocity(:));
        handles.max_bvel = max(handles.mag_backward_velocity(:));
        handles.mean_bvel = (handles.min_bvel + handles.max_bvel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_bvel handles.max_bvel];
        c.Ticks = [handles.min_bvel, (handles.min_bvel + handles.mean_bvel)/2, handles.mean_bvel, (handles.max_bvel + handles.mean_bvel)/2, handles.max_bvel];
        c.TickLabels = {num2str(handles.min_bvel,'%0.2f'), num2str((handles.min_bvel + handles.mean_bvel)/2,'%0.2f'), num2str(handles.mean_bvel,'%0.2f'), num2str((handles.max_bvel + handles.mean_bvel)/2,'%0.2f'), num2str(handles.max_bvel,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Backward Velocity [m/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        % slider adjustment
        slider_step(1) = 1/size(handles.backward_velocity,3);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.backward_velocity,3),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 1;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 17 % regurgitant flow

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.regurgitant_flow,'CDataMapping','Scaled')
        colormap(handles.axes1,'summer');
        c = colorbar(handles.axes1);
        handles.min_rfl = min(handles.regurgitant_flow(:));
        handles.max_rfl = max(handles.regurgitant_flow(:));
        handles.mean_rfl = (handles.min_rfl + handles.max_rfl)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_rfl handles.max_rfl];
        c.Ticks = [handles.min_rfl, (handles.min_rfl + handles.mean_rfl)/2, handles.mean_rfl, (handles.max_rfl + handles.mean_rfl)/2, handles.max_rfl];
        c.TickLabels = {num2str(handles.min_rfl,'%0.2f'), num2str((handles.min_rfl + handles.mean_rfl)/2,'%0.2f'), num2str(handles.mean_rfl,'%0.2f'), num2str((handles.max_rfl + handles.mean_rfl)/2,'%0.2f'), num2str(handles.max_rfl,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Regurgitant Flow [%]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_rfl handles.max_rfl]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 1;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 18 % centerline flow


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        plot3(handles.centerline_flow(:,1),handles.centerline_flow(:,2),handles.centerline_flow(:,3),'-r','LineWidth',3)
        plot3(handles.centerline(:,1),handles.centerline(:,2),handles.centerline(:,3),'-c','LineWidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 1;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;


    case 19 % eccentricity

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.eccentricity,'CDataMapping','Scaled')
        colormap(handles.axes1,'winter');
        c = colorbar(handles.axes1);
        handles.min_ecc = min(handles.eccentricity(:));
        handles.max_ecc = max(handles.eccentricity(:));
        handles.mean_ecc = (handles.min_ecc + handles.max_ecc)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_ecc handles.max_ecc];
        c.Ticks = [handles.min_ecc, (handles.min_ecc + handles.mean_ecc)/2, handles.mean_ecc, (handles.max_ecc + handles.mean_ecc)/2, handles.max_ecc];
        c.TickLabels = {num2str(handles.min_ecc,'%0.2f'), num2str((handles.min_ecc + handles.mean_ecc)/2,'%0.2f'), num2str(handles.mean_ecc,'%0.2f'), num2str((handles.max_ecc + handles.mean_ecc)/2,'%0.2f'), num2str(handles.max_ecc,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Eccentricity [%]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_ecc handles.max_ecc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 1;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 20 % curvature

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)

        if min(handles.curvature(:)) == 0 && max(handles.curvature(:)) ==0 % julio
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none')
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        else
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.curvature,'CDataMapping','Scaled')
            colormap(handles.axes1,'hot');
            c = colorbar(handles.axes1);
            handles.min_cur = min(handles.curvature(:));
            handles.max_cur = max(handles.curvature(:));
            handles.mean_cur = (handles.min_cur + handles.max_cur)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_cur handles.max_cur];
            c.Ticks = [handles.min_cur, (handles.min_cur + handles.mean_cur)/2, handles.mean_cur, (handles.max_cur + handles.mean_cur)/2, handles.max_cur];
            c.TickLabels = {num2str(handles.min_cur,'%0.2f'), num2str((handles.min_cur + handles.mean_cur)/2,'%0.2f'), num2str(handles.mean_cur,'%0.2f'), num2str((handles.max_cur + handles.mean_cur)/2,'%0.2f'), num2str(handles.max_cur,'%0.2f')};
            c.Color = [1 1 1]; % color
            c.Location = 'manual';
            c.Position = [0.2 0.2 0.02 0.3];
            c.FontWeight = 'bold';
            c.Label.String = 'Curvature [1/m]';
            windows_screen_size = get(0,'ScreenSize');
            if windows_screen_size(4)<=900
                c.FontSize = 9;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.FontSize = 10;
            else
                c.FontSize = 11;
            end
            if windows_screen_size(4)<=900
                c.Label.FontSize = 10;
            elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
                c.Label.FontSize = 12;
            else
                c.Label.FontSize = 14;
            end
            c.Label.FontWeight = 'bold';
            c.Label.Color = [1 1 1];% color
            caxis(handles.axes1, [handles.min_cur handles.max_cur]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        end

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 1; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 21 % ellipticity

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.flattening,'CDataMapping','Scaled')
        colormap(handles.axes1,'parula');
        c = colorbar(handles.axes1);
        handles.min_fla = min(handles.flattening(:));
        handles.max_fla = max(handles.flattening(:));
        handles.mean_fla = (handles.min_fla + handles.max_fla)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_fla handles.max_fla];
        c.Ticks = [handles.min_fla, (handles.min_fla + handles.mean_fla)/2, handles.mean_fla, (handles.max_fla + handles.mean_fla)/2, handles.max_fla];
        c.TickLabels = {num2str(handles.min_fla,'%0.2f'), num2str((handles.min_fla + handles.mean_fla)/2,'%0.2f'), num2str(handles.mean_fla,'%0.2f'), num2str((handles.max_fla + handles.mean_fla)/2,'%0.2f'), num2str(handles.max_fla,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Flattening [-]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_fla handles.max_fla]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 1; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

    case 22 % ellipticity

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.ellipticity,'CDataMapping','Scaled')
        colormap(handles.axes1,'parula');
        c = colorbar(handles.axes1);
        handles.min_ell = min(handles.ellipticity(:));
        handles.max_ell = max(handles.ellipticity(:));
        handles.mean_ell = (handles.min_ell + handles.max_ell)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_ell handles.max_ell];
        c.Ticks = [handles.min_ell, (handles.min_ell + handles.mean_ell)/2, handles.mean_ell, (handles.max_ell + handles.mean_ell)/2, handles.max_ell];
        c.TickLabels = {num2str(handles.min_ell,'%0.2f'), num2str((handles.min_ell + handles.mean_ell)/2,'%0.2f'), num2str(handles.mean_ell,'%0.2f'), num2str((handles.max_ell + handles.mean_ell)/2,'%0.2f'), num2str(handles.max_ell,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Ellipticity [-]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_ell handles.max_ell]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 1; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;



    case 23 % length_vessel

        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.length_vessel,'CDataMapping','Scaled')
        colormap(handles.axes1,'parula');
        c = colorbar(handles.axes1);
        handles.min_len = min(handles.length_vessel(:));
        handles.max_len = max(handles.length_vessel(:));
        handles.mean_len = (handles.min_len + handles.max_len)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_len handles.max_len];
        c.Ticks = [handles.min_len, (handles.min_len + handles.mean_len)/2, handles.mean_len, (handles.max_len + handles.mean_len)/2, handles.max_len];
        c.TickLabels = {num2str(handles.min_len,'%0.2f'), num2str((handles.min_len + handles.mean_len)/2,'%0.2f'), num2str(handles.mean_len,'%0.2f'), num2str((handles.max_len + handles.mean_len)/2,'%0.2f'), num2str(handles.max_len,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Length of the Vessel [cm]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_len handles.max_len]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 1; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;

%     case 24 % Circulation
%         
%         
%         set(handles.popupmenu1,'Value',1)
%         
%         axes(handles.axes1);
%         plot(0.0)
%         patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
%         colormap(handles.axes1,'winter');
%         c = colorbar(handles.axes1);
%         handles.min_cir = min(handles.circulation(:));
%         handles.max_cir = max(handles.circulation(:));
%         handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%         c.LimitsMode = 'manual';
%         c.Limits = [handles.min_cir handles.max_cir];
%         c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%         c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%         c.Color = [1 1 1]; % color
%         c.Location = 'manual';
%         c.Position = [0.2 0.2 0.02 0.3];
%         c.FontWeight = 'bold';
%         c.Label.String = 'Circulation [mm^{2}/s]';
%         windows_screen_size = get(0,'ScreenSize');
%         if windows_screen_size(4)<=900
%             c.FontSize = 9;
%         elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%             c.FontSize = 10;
%         else
%             c.FontSize = 11;
%         end
%         if windows_screen_size(4)<=900
%             c.Label.FontSize = 10;
%         elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%             c.Label.FontSize = 12;
%         else
%             c.Label.FontSize = 14;
%         end
%         c.Label.FontWeight = 'bold';
%         c.Label.Color = [1 1 1];% color
%         caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%         hold off
%         axis vis3d
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
%         
%         
%         % slider adjustment
%         slider_step(1) = 1/size(handles.circulation,2);
%         slider_step(2) = 0.1;
%         set(handles.slider1,'Value', handles.peak_flow/size(handles.circulation,2),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
%         set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         handles.id_mesh = 0;
%         handles.id_mesh_vel = 0;
%         handles.id_mesh_inlet = 0;
%         handles.id_mesh_outlet = 0;
%         handles.id_mesh_laplace = 0;
%         handles.id_centerline = 0;
%         handles.id_diameter = 0;
%         handles.id_radius = 0;
%         handles.id_axial_unit_vectors = 0;
%         handles.id_circumferential_unit_vectors = 0;
%         handles.id_wss_a = 0;
%         handles.id_wss_c = 0;
%         handles.id_axial_angle = 0;
%         handles.id_forward_vel = 0;
%         handles.id_backward_vel = 0;
%         handles.id_regurgitant_flow = 0;
%         handles.id_centerline_flow = 0;
%         handles.id_eccentricity = 0;
%         handles.id_curvature = 0; % new data Julio Sotelo
%         handles.id_flattening = 0; % new data Julio Sotelo
%         handles.id_ellipticity = 0; % new data Julio Sotelo
%         handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 1; % new data Julio Sotelo
%         handles.id_forward_vortex = 0; % new data Julio Sotelo
%         handles.id_ipcmra = 0;
%         handles.id_mag = 0;    
% 
%         set(handles.popupmenu1,'value',1)
%         set(handles.popupmenu3,'value',1)
%         set(handles.pushbutton1,'visible','off')
%         set(handles.pushbutton2,'visible','off')

    case 24 % Forward Vortex


        set(handles.popupmenu1,'Value',1)

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.forward_vortex(:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        c = colorbar(handles.axes1);
        handles.min_fov = min(handles.forward_vortex(:));
        handles.max_fov = max(handles.forward_vortex(:));
        handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_fov handles.max_fov];
        c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
        c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Forward Vortex [1/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_fov handles.max_fov]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        % slider adjustment
        slider_step(1) = 1/size(handles.forward_vortex,2);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.forward_vortex,2),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 1; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;   

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')

    case 25 % area

        set(handles.popupmenu1,'Value',1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.area,'CDataMapping','Scaled')
        colormap(handles.axes1,'parula');
        c = colorbar(handles.axes1);
        handles.min_are = min(handles.area(:));
        handles.max_are = max(handles.area(:));
        handles.mean_are = (handles.min_are + handles.max_are)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_are handles.max_are];
        c.Ticks = [handles.min_are, (handles.min_are + handles.mean_are)/2, handles.mean_are, (handles.max_are + handles.mean_are)/2, handles.max_are];
        c.TickLabels = {num2str(handles.min_are,'%0.2f'), num2str((handles.min_are + handles.mean_are)/2,'%0.2f'), num2str(handles.mean_are,'%0.2f'), num2str((handles.max_are + handles.mean_are)/2,'%0.2f'), num2str(handles.max_are,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Area [cm^{2}]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_are handles.max_are]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        set(handles.slider1,'visible','off')
        set(handles.text1,'visible','off')
        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_area = 1; % new data Julio Sotelo
        handles.id_axial_circulation = 0; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;   

    case 26 % axial circulation

        set(handles.popupmenu1,'Value',1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.axial_circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
        colormap(handles.axes1,'autumn');
        c = colorbar(handles.axes1);
        handles.min_acir = min(handles.axial_circulation(:));
        handles.max_acir = max(handles.axial_circulation(:));
        handles.mean_acir = (handles.min_acir + handles.max_acir)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_acir handles.max_acir];
        c.Ticks = [handles.min_acir, (handles.min_acir + handles.mean_acir)/2, handles.mean_acir, (handles.max_acir + handles.mean_acir)/2, handles.max_acir];
        c.TickLabels = {num2str(handles.min_acir,'%0.2f'), num2str((handles.min_acir + handles.mean_acir)/2,'%0.2f'), num2str(handles.mean_acir,'%0.2f'), num2str((handles.max_acir + handles.mean_acir)/2,'%0.2f'), num2str(handles.max_acir,'%0.2f')};
        c.Color = [1 1 1]; % color
        c.Location = 'manual';
        c.Position = [0.2 0.2 0.02 0.3];
        c.FontWeight = 'bold';
        c.Label.String = 'Axial Circulation [cm^{2}/s]';
        windows_screen_size = get(0,'ScreenSize');
        if windows_screen_size(4)<=900
            c.FontSize = 9;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.FontSize = 10;
        else
            c.FontSize = 11;
        end
        if windows_screen_size(4)<=900
            c.Label.FontSize = 10;
        elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
            c.Label.FontSize = 12;
        else
            c.Label.FontSize = 14;
        end
        c.Label.FontWeight = 'bold';
        c.Label.Color = [1 1 1];% color
        caxis(handles.axes1, [handles.min_acir handles.max_acir]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        % slider adjustment
        slider_step(1) = 1/size(handles.axial_circulation,2);
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/size(handles.axial_circulation,2),'sliderstep',slider_step,'max',1,'min',0,'visible','on')
        set(handles.text1,'String',['# CP: ',num2str(handles.peak_flow)],'visible','on')

        handles.id_mesh = 0;
        handles.id_mesh_vel = 0;
        handles.id_mesh_inlet = 0;
        handles.id_mesh_outlet = 0;
        handles.id_mesh_laplace = 0;
        handles.id_centerline = 0;
        handles.id_diameter = 0;
        handles.id_radius = 0;
        handles.id_axial_unit_vectors = 0;
        handles.id_circumferential_unit_vectors = 0;
        handles.id_wss_a = 0;
        handles.id_wss_c = 0;
        handles.id_axial_angle = 0;
        handles.id_forward_vel = 0;
        handles.id_backward_vel = 0;
        handles.id_regurgitant_flow = 0;
        handles.id_centerline_flow = 0;
        handles.id_eccentricity = 0;
        handles.id_curvature = 0; % new data Julio Sotelo
        handles.id_flattening = 0; % new data Julio Sotelo
        handles.id_ellipticity = 0; % new data Julio Sotelo
        handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
        handles.id_forward_vortex = 0; % new data Julio Sotelo
        handles.id_area = 0; % new data Julio Sotelo
        handles.id_axial_circulation = 1; % new data Julio Sotelo
        handles.id_ipcmra = 0;
        handles.id_mag = 0;   

        set(handles.popupmenu1,'value',1)
        set(handles.popupmenu3,'value',1)
        set(handles.pushbutton1,'visible','off')
        set(handles.pushbutton2,'visible','off')
end
handles.output = hObject;  
guidata(hObject, handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


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


% --- Executes when user attempts to close GUIDE_LAPLACE.
function GUIDE_LAPLACE_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to GUIDE_LAPLACE (see GCBO)
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
% Hint: delete(hObject) closes the figure
% delete(hObject);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.id_vwerp  == 1
    
    addpath(genpath('vWERP/'));
    
    inlet_seg = handles.PLANE_INLET;
    inlet = handles.PLANE_INLET;
    n_inlet = handles.normal_inlet/norm(handles.normal_inlet);

    outlet_seg = handles.PLANE_OUTLET;
    outlet = handles.PLANE_OUTLET;
    n_outlet = handles.normal_outlet/norm(handles.normal_outlet);
    
    mask = handles.SEG_for_vwerp;
    
    % adjusting the velocities
%     v = [];
%     if max(abs(handles.MR_PCA_FH(:)))>10
%         v{1}.im = handles.MR_PCA_AP.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_AP,4))/100;
%         v{2}.im = (handles.MR_PCA_FH.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_FH,4))/100);
%         v{3}.im = (handles.MR_PCA_RL.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_RL,4))/100);
%     else
%         v{1}.im = handles.MR_PCA_AP.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_AP,4));
%         v{2}.im = (handles.MR_PCA_FH.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_FH,4)));
%         v{3}.im = (handles.MR_PCA_RL.*repmat(handles.SEG,1,1,1,size(handles.MR_PCA_RL,4)));
%     end

    % adjusting the velocities
    v = [];
    if max(abs(handles.MR_PCA_FH(:)))>10
        v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4))/100)*-1;
        v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4))/100)*-1;
        v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4))/100);
    else
        v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4)))*-1;
        v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4)))*-1;
        v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4)));
    end
    
    % Check visualization - hardcoded by David 
    % Modified Julio Sotelo 21-09-2022
    handles.v_for_plot = v;
    handles.mask_for_plot = mask;
    handles.inlet_for_plot = inlet;
    handles.outlet_for_plot = outlet;
    handles.n_inlet_for_plot = n_inlet;
    handles.n_outlet_for_plot = n_outlet;
    
%     figure,
%     visu_mask_inlet_outlet(v, mask,inlet,outlet,n_inlet,n_outlet)

    % adjusting the dt
    t = linspace(0,(60/handles.heart_rate),size(handles.MR_PCA_AP,4));
    v{1}.dt = t(2);
    v{2}.dt = t(2);
    v{3}.dt = t(2);
    
    if handles.voxel_MR(1)>0.1
        v{1}.dx = handles.voxel_MR/1000;
        v{2}.dx = handles.voxel_MR/1000;
        v{3}.dx = handles.voxel_MR/1000;
    else
        v{1}.dx = handles.voxel_MR;
        v{2}.dx = handles.voxel_MR;
        v{3}.dx = handles.voxel_MR;
    end
    
    % adjusting the resolution
    v{1}.res = size(handles.MR_PCA_AP);
    v{2}.res = size(handles.MR_PCA_AP);
    v{3}.res = size(handles.MR_PCA_AP);
    
    v_old = v;
   
           
    werp = handles.werp; %werp --> 1, vwerp --> 2
    volInt = handles.volInt; %true --> volInt, false --> surfInt, for advective energy
    percentage_noise = handles.percentage_noise;
    resamp_factor = handles.resamp_factor; %resampling factor
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    turbulence = false;
    T = 1;
    omega = 2*pi/T;
    period = 2*pi/omega;
    rho = 1060;
    mu = 0.004;
    pa2mmhg = 0.00750061683;

%     %% Load data  (v,p,inlet plane, outlet plane, segmentation mask)
%     if ~turbulence
%         load('test/v.mat');
%         load('test/mask.mat');
%         load('test/inlet.mat'); inlet = plane; n_inlet = n;
%         load('test/outlet.mat'); outlet = plane; n_outlet = -n;
% 
%     %     load('temp_output_vwerp/v.mat');
%     %     load('temp_output_vwerp/mask.mat');
%     %     load('temp_output_vwerp/inlet.mat'); n_inlet = -n_inlet;
%     %     load('temp_output_vwerp/outlet.mat'); n_outlet = n_outlet;
% 
%     else % message with a questions Julio Sotelo
%         load('test_turb/v.mat');
%         load('test_turb/mask.mat');
%         load('test_turb/inlet.mat'); inlet = plane; n_inlet = n;
%         load('test_turb/outlet.mat'); outlet = plane; n_outlet = -n;
%         load('test_turb/cov.mat');
%     end

    %% Make compatible with previous notation
    for i = 1:3
        v{i}.PixDim = v{i}.dx;
        if size(v{1}.im,4) > 1
            dt = v{i}.dt;
        end
        if turbulence
            for j = 1:3
                Cov{i,j}.PixDim = Cov{i,j}.dx;
            end
        end
    end

    %% Add empty slices on all sides of masks and v field
    inlet = pad_zeros(inlet);
    outlet = pad_zeros(outlet);
    mask = pad_zeros(mask);
    % v field
    for i = 1:3
        v{i}.im = pad_zeros(v{i}.im);
        if turbulence
            for j = 1:3
                Cov{i,j}.im = pad_zeros(Cov{i,j}.im);
            end
        end
    end

    maskin = inlet;
    outlet_orig = outlet;
    inlet_mask = inlet;
    outlet_mask = outlet;
    clear inlet outlet

    if size(v{1}.im,4) > 1
        number_of_time_slices = size(v{1}.im,4);
        times = 0:period/number_of_time_slices:period;
    else %static
        number_of_time_slices = size(v{1}.im,4) + 1; %since evaluation is set between time steps
    end

    %% Clean up v - spatial
    v_orig = v;
    if percentage_noise ~= 0
        for i = 1:3 
            %Find ideal filter parameters
            if i == 1
                Venc = 2;
                D = 10*2*10^-3; 
                discretization_points = ceil(D/v{1}.dx(1));
                % Find ideal filter parameter
                [ksize, smooth_interp_order, type] = find_vel_filter(Venc, D, v{1}.PixDim, discretization_points, percentage_noise);
                if ksize == [1 1 1]
                    ksize = [3 3 3];
                end
            end

            %Apply filter
            for T = 1:size(v{1}.im,4);
                fprintf('Filtering dimension %i, time step %i\n',i,T)
                u = v{i}.im(:,:,:,T);
                h = v{i}.PixDim;
                ncalc = 2;
                global exterior
                % Least squares smoothing
                [mask_smoothing, u, ux, uy, uz, ~, ~, ~, ~, ~, ~] = computeLeastSquaresDerivatives_mask( u, h, smooth_interp_order, ksize, mask, ncalc);

                v_smooth{i}.im = u .* mask;
                v_smooth{i}.PixDim = v{1}.PixDim;
                Grv_smooth{i,1}.im = ux; Grv_smooth{i,2}.im = uy; Grv_smooth{i,3}.im = uz;
                for j = 1:3
                    Grv_smooth{i,j}.PixDim = v{1}.PixDim;
                end
                v{i}.im(:,:,:,T) = v_smooth{i}.im;
            end
        end
    end

    %% Clean Covariance if turbulent data
    if turbulence
        %Mask Cov such that diagonal indices are strictly positive Cov{i,i}
        cov_mask = mask;
        for i = 1:3
           cov_mask = cov_mask .* (Cov{i,i}.im > 0);
        end

        for i = 1:3
           for j = 1:3
                Cov{i,j}.im = cov_mask .* Cov{i,j}.im;
            end
        end
    end

    %%
    v_temporal = v;
    mask_orig = mask;
    clear v

    for T = 1:number_of_time_slices-1; %Compute solution in-between time steps    
        if size(v_temporal{1}.im,4) > 1
            %Compute solution in-between time steps if temporal data...
            %
            % That means that we compute for v_{i+1/2} = 0.5*(v_{i}+v_{i+1})
            for i = 1:length(v_temporal{1}.PixDim);
                v1{i}.im = v_temporal{i}.im(:,:,:,T);
                v1{i}.PixDim = v_temporal{i}.PixDim;
                v2{i}.im = v_temporal{i}.im(:,:,:,T+1);
                v2{i}.PixDim = v_temporal{i}.PixDim;
            end

            %% Get v_{i+1/2}
            for i = 1:length(v_temporal{1}.PixDim)
                v{i}.im = (1/2)*(v1{i}.im + v2{i}.im);
                v{i}.PixDim = v1{i}.PixDim;
            end

            % Check velocity shape
            if isrow(v)
                v = v';
            end

            %% Get dvdt at t_{i+1/2}
            dvdt = ddt_image(v1,v2,dt);
        else %static data
            v = v_temporal;
            % Check velocity shape
            if isrow(v)
                v = v';
            end
            dvdt = v; % Dummy dvdt... (only used to initialize kinetic, which then = 0)
        end

        %% Compute gradients \nabla v
        Grv = gradient_vimage_mask(v, mask_orig);

        %% Data interpolation
        if resamp_factor ~= 1
            v = resample_vel_data_v2(v, resamp_factor);
            dvdt = resample_vel_data_v2(dvdt, resamp_factor);
            Grv = resample_Grv_data_v2(Grv, resamp_factor);
            if turbulence
                Cov = resample_Grv_data_v2(Cov, resamp_factor);
            end
            if T == 1 %Only resample mask, inlet and outlet once
                mask_resampled = resample_mask_v2(mask, resamp_factor);
                mask = mask_resampled;
                maskin = resample_mask_v2(maskin, resamp_factor);
                outlet_orig = resample_mask_v2(outlet_orig, resamp_factor);
            else
                mask = mask_resampled;
            end
        end

        % Check velocity shape
        if isrow(v)
            v = v';
        end

        v_orig = v;

        outid = 1; %You can change this when you add more outlet planes
        maskout = (outlet_orig == outid);

        %% Get label map (exterior, interior, inlet, outlet)...
        if T == 1;

            %%%%%%%%% incorporate labels to the visualization Julio Sotelo
            Labels = werp2_make_labels(mask, maskin, maskout);
        end
        
        %% Solve w using Stoke's equation
        % Get the chosen velocity field w from solving Stokes equation (FD meth)...
        % This is probably the trickiest part!
        if werp == 1 %werp1        w = v;
        elseif werp == 2 %werp2
            if T == 1 %Only need to solve w for one time step or 1 noise sample

                % the w need to be incorpotated to the visualization Julio
                % Sotelo
               [w, div_err] = werp2_get_w_field_normal_vector(v, mask, maskin, maskout, Labels, n_inlet);
               %[w, div_err] = werp2_get_w_field(v, mask, maskin, maskout, Labels);
            end
        end

        %% Calculate gradients, and WERP
        % Declare global variables...
        global interior inlet exterior outlet

        % Compute gradients \nabla v and \nabla w...
        Grw = gradient_vimage_mask(w, (Labels ~= exterior));

        % Get vel magnitude (for fun)...
        vmag = magnitude_vimage(v); %vel_mag
        wmag = magnitude_vimage(w); %w_mag

        % Get individual energy contributors...
        visc   = werp2_get_viscous_normal_vector(Grv, Grw, mu, w, Labels, n_inlet, n_outlet);
        lambda = werp2_get_lambda_normal_vector(w, Labels, inlet, n_inlet); %Lambda = Q = flow
        if volInt
            advec  = werp2_get_advective_volInt(v, w, Grv, Labels, rho); %For volInt
        else
            advec  = werp2_get_advective_surfInt_normal_vector(v, w, Grw, Labels, rho, n_inlet, n_outlet);
        end
        kinetic = werp2_get_kinetic_energy(dvdt, w, rho) * (size(v_temporal{1}.im,4) > 1);
        if turbulence
            turb = werp2_get_turbulence(Cov,Grw,rho);
        else
            turb = 0;
        end

        % Derive relative pressure drop estimated by WERP...
        dP = - (visc + advec + kinetic + turb) / lambda; 

        % Print results...
        if size(v_temporal{1}.im,4) > 1
            if T == 1
                disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            end
            fprintf('Time: %f sec out of %f sec.\n', 0.5*(times(T)+times(T+1)), 0.5*(times(end-1)+times(end)))
        end
        fprintf('Estimated  dP: mmHg %f\n', dP.*pa2mmhg)
        fprintf('kinetic: %f, advec: %f, visc: %f, turb: %f mmHg',...
                -kinetic./lambda .* pa2mmhg,...
                -advec./lambda .* pa2mmhg,...
                -visc./lambda .* pa2mmhg,...
                -turb./lambda .* pa2mmhg);
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        LAMBDA(T) = lambda;
        KINETIC(T) = kinetic;
        ADVEC(T) = advec;
        VISC(T) = visc;  
        if turbulence
            TURB(T) = turbulence;
        end

    end

    %% Visualize
    if turbulence
        DP_VWERP = - (VISC + ADVEC + KINETIC + TURB) ./ LAMBDA; 
    else
        DP_VWERP = - (VISC + ADVEC + KINETIC) ./ LAMBDA; 
    end

    % mofied by Julio Sotelo 21-09-2022
    if ~turbulence %No visu for static turblence test data
%         figure
        time = linspace(dt+dt/2,(number_of_time_slices*dt)-dt/2,number_of_time_slices-1);
%         plot(time, [DP_VWERP(end), DP_VWERP(1:end-1)]*pa2mmhg,'b--*')
%         hold on
%         plot(time, zeros(size(time)),'k-')
%         xlabel('time [s]')
%         ylabel('Relative pressure drop [mmHg]')
%         grid on
%         title(['Estimated relative pressure drops, WERP']);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.inlet_seg   = inlet_seg;
    handles.n_inlet     = n_inlet;
    handles.inlet       = inlet;
    handles.outlet_seg  = outlet_seg;
    handles.n_outlet    = n_outlet;
    handles.outlet      = outlet;
    handles.mask        = mask;
    handles.v           = v;
    handles.v_old       = v_old;
    handles.time        = time;
    handles.DP_VWERP    = DP_VWERP;
    handles.pa2mmhg     = pa2mmhg;
    handles.VISC        = VISC;
    handles.ADVEC       = ADVEC;
    handles.KINETIC     = KINETIC;
    handles.LAMBDA      = LAMBDA;
    handles.w           = w;
    handles.Labels      = Labels;

    set(handles.pushbutton7,'Visible','on')

    f = msgbox("vWERP quantification finished","Success");
    pause(1)
    close(f)


    % we add two more rows in the popupmenu
    handles.list_string = {'Select Results ...','Mesh','Mesh + Velocity','Segmentation + Velocity','Segmentation + Inlet & Outlet','Mesh + Inlet','Mesh + Outlet','Relative pressure drop [mmHg]'};
    set(handles.popupmenu2,'String',handles.list_string,'Value',8);


    set(handles.uipanel2,'BackgroundColor',[1,1,1])

    axes(handles.axes1);
    plot(time, [DP_VWERP(end), DP_VWERP(1:end-1)]*pa2mmhg,'b--*','LineWidth',2, 'MarkerSize', 10)
    hold on
    plot(time, zeros(size(time)),'k-')
    hold off
    xlabel('time [s]', 'Fontsize', 14,'FontWeight','bold')
    ylabel('Relative pressure drop [mmHg]', 'Fontsize', 14,'FontWeight','bold')
    grid on
    title(['Estimated relative pressure drops, WERP'], 'Fontsize', 16,'FontWeight','bold');


    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    set(handles.popupmenu1,'value',1)
    set(handles.popupmenu3,'value',1)
    set(handles.pushbutton1,'visible','off')
    set(handles.pushbutton2,'visible','off')
    set(handles.pushbutton3,'Visible','off')
    set(handles.pushbutton5,'Visible','off')
    set(handles.pushbutton6,'Visible','off')
    set(handles.text3,'Visible','off')
    set(handles.text4,'Visible','off')
    set(handles.text5,'Visible','off')
    set(handles.text6,'Visible','off')
    set(handles.radiobutton1,'Visible','off')
    set(handles.radiobutton2,'Visible','off')
    set(handles.radiobutton3,'Visible','off')
    set(handles.radiobutton4,'Visible','off')
    set(handles.edit2,'Visible','off')
    set(handles.edit3,'Visible','off')



%     mkdir('Output_vWERP/')
%     save('Output_vWERP/inlet.mat','inlet_seg','n_inlet', 'inlet')
%     save('Output_vWERP/outlet.mat','outlet_seg','n_outlet', 'outlet')
%     save('Output_vWERP/mask.mat','mask')
%     save('Output_vWERP/v.mat','v','v_old')
%     save('Output_vWERP/results.mat','time','DP_VWERP','pa2mmhg','VISC', 'ADVEC', 'KINETIC', 'LAMBDA')
%     save('Output_vWERP/w.mat','w')
%     save('Output_vWERP/labels.mat','Labels')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else

    [out]                           = FAST_LAPLACE(handles.elem,handles.nodes,handles.id_inlet,handles.id_outlet);
    
    handles.Laplace                 = out.laplace; % Laplace
    h                               = out.h;
    close(h)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.id_wall = unique(handles.faces(:));
    id_n = ismember(handles.id_wall,handles.id_inlet);
    handles.id_wall(id_n)=[];
    id_n = ismember(handles.id_wall,handles.id_outlet);
    handles.id_wall(id_n)=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Centerline 
    input = [];
    input.faces                     = handles.faces;
    input.nodes                     = handles.nodes;
    input.id_wall                   = handles.id_wall;
    input.Laplace                   = handles.Laplace;
    
    [out]                           = CENTERLINE(input);
    
    handles.centerline              = out.centerline;
    handles.centerline_lapid        = out.centerline_lapid;
    handles.centerline_n_id_wall    = out.centerline_n_id_wall;
    h                               = out.h;
    close(h)
    
    if handles.id_section_loaded == 1
       handles.normal_inlet = (handles.centerline(2,:) - handles.centerline(1,:))/norm(handles.centerline(2,:) - handles.centerline(1,:));
       handles.normal_outlet = (handles.centerline(end,:) - handles.centerline(end-1,:))/norm(handles.centerline(end,:) - handles.centerline(end-1,:));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % length vessel and curvature
    input = [];    
    input.centerline                = handles.centerline;
    input.centerline_n_id_wall      = handles.centerline_n_id_wall;
    input.id_wall                   = handles.id_wall;
    input.Laplace                   = handles.Laplace;
    
    [out]                           = LENGTH_AND_CURVATURE(input);
    
    handles.length_vessel           = out.length_vessel/10; % length_vessel in cm
    handles.curvature               = out.curvature;
    h                               = out.h;
    close(h)
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Diameter and Radius quantification
    input = [];
    input.Laplace                   = handles.Laplace;
    input.id_wall                   = handles.id_wall;
    input.centerline_n_id_wall      = handles.centerline_n_id_wall;
    input.centerline                = handles.centerline;
    input.nodes                     = handles.nodes;
    input.faces                     = handles.faces;
    
    [out]                           = RADIUS_DIAMETER_FLATTENING_ELLIPTICITY(input);
        
    handles.diameter                = out.diameter/10; % diameter in cm
    handles.radius                  = out.radius/10;% radius in cm
    handles.flattening              = out.flattening;
    handles.ellipticity             = out.ellipticity;
    h                               = out.h;
    close(h)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fluid zone
    handles.id_fluidzone            = unique(handles.elem(:));
    id_n                            = ismember(handles.id_fluidzone,handles.id_inlet);
    handles.id_fluidzone(id_n)      = [];
    id_n                            = ismember(handles.id_fluidzone,handles.id_outlet);
    handles.id_fluidzone(id_n)      = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input = [];
    input.elem                      = handles.elem;
    input.faces                     = handles.faces;
    input.nodes                     = handles.nodes;
    input.Laplace                   = handles.Laplace;
    input.id_inlet                  = handles.id_inlet;
    input.id_outlet                 = handles.id_outlet;
    input.centerline                = handles.centerline;
    input.centerline_lapid          = handles.centerline_lapid;
    input.id_section_loaded         = handles.id_section_loaded;
    input.normal_inlet              = handles.normal_inlet;
    input.normal_outlet             = handles.normal_outlet;
    input.id_fluidzone              = handles.id_fluidzone;

    [out]                           = AXIAL_CIRCUMFERENTIAL_UNIT_VECTORS(input);
    
    handles.axial_unit_vectors      = out.axial_unit_vectors;
    handles.circumferential_unit_vectors = out.circumferential_unit_vectors;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% WSS decomposition axial and circumferential direction
    handles.WSS_A = repmat(sum(handles.WSS.*repmat(handles.axial_unit_vectors,1,1,size(handles.WSS,3)),2),1,3,1).* handles.axial_unit_vectors;
    handles.WSS_C = repmat(sum(handles.WSS.*repmat(handles.circumferential_unit_vectors,1,1,size(handles.WSS,3)),2),1,3,1).* handles.circumferential_unit_vectors;
    
    handles.mag_WSS_A = sqrt(sum(handles.WSS_A.^2,2));
    handles.mag_WSS_C = sqrt(sum(handles.WSS_C.^2,2));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% axial angle
    handles.angle_axial_direction = acos(sum(handles.veset.*repmat(handles.axial_unit_vectors,1,1,size(handles.veset,3)),2)./(sqrt(sum(handles.veset.^2,2)).*sqrt(sum(repmat(handles.axial_unit_vectors,1,1,size(handles.veset,3)).^2,2))))*180/pi;
    handles.angle_axial_direction(isnan(handles.angle_axial_direction))=0;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Fordward and Backward Velicity
    axial_projection_velocity = repmat(sum(handles.veset.*repmat(handles.axial_unit_vectors,1,1,size(handles.veset,3)),2),1,3,1).*handles.axial_unit_vectors;
    
    handles.forward_velocity = axial_projection_velocity.*repmat(handles.angle_axial_direction<=90,1,3,1);
    handles.backward_velocity = axial_projection_velocity.*repmat(handles.angle_axial_direction>90,1,3,1);

    % generation of diferent variables
    handles.mag_forward_velocity = sqrt(sum(handles.forward_velocity.^2,2));
    handles.mag_backward_velocity = sqrt(sum(handles.backward_velocity.^2,2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DIRECTION OF FORDWARD VORTEX
    handles.forward_vortex = squeeze(sum( repmat(handles.axial_unit_vectors,1,1,size(handles.vorticity,3)).*handles.vorticity,2));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Regurgitant Flow
    input = [];
    input.elem                      = handles.elem;
    input.faces                     = handles.faces;
    input.nodes                     = handles.nodes;
    input.heart_rate                = handles.heart_rate;
    input.Laplace                   = handles.Laplace;
    input.radius                    = handles.radius*10; % radius in mm
    input.diameter                  = handles.diameter*10; % diameter in mm
    input.veset                     = handles.veset;
    input.centerline                = handles.centerline;
    input.id_wall                   = handles.id_wall;
    input.centerline_lapid          = handles.centerline_lapid;
    input.centerline_n_id_wall      = handles.centerline_n_id_wall;
    input.id_fluidzone              = handles.id_fluidzone;
    input.mag_forward_velocity      = handles.mag_forward_velocity;
    input.mag_backward_velocity     = handles.mag_backward_velocity;
    input.forward_vortex            = handles.forward_vortex;
    input.list_n                    = handles.list_n;
    input.peak_flow_ori             = handles.peak_flow_ori;
        
    % Incluir area aca
    [out]                           = REGURGITANT_FLOW_AND_ECCENTRICITY(input);

    handles.regurgitant_flow        = out.regurgitant_flow;
    handles.centerline_flow         = out.centerline_flow;
    handles.eccentricity            = out.eccentricity;
    handles.area                    = out.area*0.01; % area un cm^2
    handles.axial_circulation       = out.axial_circulation*0.01; % axial_circulation un cm^2/s
    
    h = out.h;

    close(h)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Circulation quantification
%     input = [];
%     input.Laplace = handles.Laplace;
%     input.SEG = handles.SEG;
%     input.faces = handles.faces;
%     input.nodes = handles.nodes;
%     input.reduced_centerline = handles.reduced_centerline;
%     input.voxel_MR = handles.voxel_MR;
%     input.MR_PCA_AP = handles.MR_PCA_AP;
%     input.MR_PCA_FH = handles.MR_PCA_FH;
%     input.MR_PCA_RL = handles.MR_PCA_RL;
%     input.centerline = handles.centerline;
%     input.id_wall = handles.id_wall;
%     input.centerline_n_id_wall = handles.centerline_n_id_wall;
%         
%     [out] = CIRCULATION_MAP(input);
% 
%     handles.circulation = out.circulation;
%     h=out.h;
%     close(h)
% 
%     handles.circulation = zeros(size(handles.axial_circulation));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     handles.list_string = {'Select Results ...','Mesh','Velocity Vectors','Mesh + Inlet IDs','Mesh + Outlet IDs','Laplace',...
%                                                 'Centerline','Diameter','Radius','Axial Unit Vectors',...
%                                                 'Circum. Unit Vectors','WSSa','WSSc','Axial Angle','Forward Velocity','Backward Velocity',...
%                                                 'Regurgitant Flow','Centerline Flow','Eccentricity',...
%                                                 'Curvature','Flattening','Ellipticity','Length of Vessel','Circulation','Forward Vortex','Area','Axial Circulation'};

    handles.list_string = {'Select Results ...','Mesh','Velocity Vectors','Mesh + Inlet IDs','Mesh + Outlet IDs','Laplace',...
                                                'Centerline','Diameter','Radius','Axial Unit Vectors',...
                                                'Circum. Unit Vectors','WSSa','WSSc','Axial Angle','Forward Velocity','Backward Velocity',...
                                                'Regurgitant Flow','Centerline Flow','Eccentricity',...
                                                'Curvature','Flattening','Ellipticity','Length of Vessel','Forward Vortex','Area','Axial Circulation'};
                                            
    set(handles.popupmenu2,'String',handles.list_string,'Value',6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.id_mesh = 0;
    handles.id_mesh_vel = 0;
    handles.id_mesh_inlet = 0;
    handles.id_mesh_outlet = 0;
    handles.id_mesh_laplace = 1;
    handles.id_centerline = 0;
    handles.id_diameter = 0;
    handles.id_radius = 0;
    handles.id_axial_unit_vectors = 0;
    handles.id_circumferential_unit_vectors = 0;
    handles.id_wss_a = 0;
    handles.id_wss_c = 0;
    handles.id_axial_angle = 0;
    handles.id_forward_vel = 0;
    handles.id_backward_vel = 0;
    handles.id_regurgitant_flow = 0;
    handles.id_centerline_flow = 0;
    handles.id_eccentricity = 0;
    handles.id_curvature = 0; % new data Julio Sotelo
    handles.id_flattening = 0; % new data Julio Sotelo
    handles.id_ellipticity = 0; % new data Julio Sotelo
    handles.id_length_vessel = 0; % new data Julio Sotelo
%     handles.id_circulation = 0; % new data Julio Sotelo
    handles.id_forward_vortex = 0; % new data Julio Sotelo
    handles.id_area = 0; % new data Julio Sotelo
    handles.id_axial_circulation = 0; % new data Julio Sotelo
    handles.id_ipcmra = 0;
    handles.id_mag = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.popupmenu1,'Value',1)
        
    axes(handles.axes1);
    plot(0.0)
    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace,'CDataMapping','Scaled')
    colormap(handles.axes1,'cool');
    c = colorbar(handles.axes1);
    handles.min_lap = min(handles.Laplace(:));
    handles.max_lap = max(handles.Laplace(:));
    handles.mean_lap = (handles.min_lap + handles.max_lap)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_lap handles.max_lap];
    c.Ticks = [handles.min_lap, (handles.min_lap + handles.mean_lap)/2, handles.mean_lap, (handles.max_lap + handles.mean_lap)/2, handles.max_lap];
    c.TickLabels = {num2str(handles.min_lap,'%0.2f'), num2str((handles.min_lap + handles.mean_lap)/2,'%0.2f'), num2str(handles.mean_lap,'%0.2f'), num2str((handles.max_lap + handles.mean_lap)/2,'%0.2f'), num2str(handles.max_lap,'%0.2f')};
    c.Color = [1 1 1]; % color
    c.Location = 'manual';
    c.Position = [0.2 0.2 0.02 0.3];
    c.FontWeight = 'bold';
    c.Label.String = 'Laplace Value [-]';
    windows_screen_size = get(0,'ScreenSize');
    if windows_screen_size(4)<=900
        c.FontSize = 9;
    elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
        c.FontSize = 10;
    else
        c.FontSize = 11;
    end
    if windows_screen_size(4)<=900
        c.Label.FontSize = 10;
    elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
        c.Label.FontSize = 12;
    else
        c.Label.FontSize = 14;
    end
    c.Label.FontWeight = 'bold';
    c.Label.Color = [1 1 1];% color
    caxis(handles.axes1, [handles.min_lap handles.max_lap]);
    hold off
    axis vis3d
    daspect([1,1,1])
    axis off
    view([-34,-51])

    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    
end

handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

set(handles.uipanel2,'BackgroundColor',[0,0,0])
switch get(handles.popupmenu3,'Value')
    
    case 1
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        handles.id_view = 1;
        if handles.id_ipcmra ==1
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.slider_axes1 = round(handles.c/2);
            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            
        elseif handles.id_mag == 1

            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.slider_axes1 = round(handles.c/2);
            
            % slider adjustment
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,:,handles.slider_axes1)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%             himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            
        end

    case 2
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        handles.id_view = 2;
        if handles.id_ipcmra ==1
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');

            handles.slider_axes1 = round(handles.b/2);
            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        elseif handles.id_mag == 1
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.slider_axes1 = round(handles.b/2);
            % slider adjustment
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.MAG(:,handles.slider_axes1,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes1,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes1,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            
        end
        

    case 3
        set(handles.pushbutton5,'Visible','off')
        set(handles.pushbutton6,'Visible','off')
        handles.id_view = 3;
        if handles.id_ipcmra ==1

            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.slider_axes1 = round(handles.a/2);
            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
        
        elseif handles.id_mag == 1
            
            set(handles.pushbutton1,'Visible','on');
            set(handles.pushbutton2,'Visible','on');
            
            handles.slider_axes1 = round(handles.a/2);
            % slider adjustment
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])

            axes(handles.axes1);
            plot(0.0)
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.MAG(handles.slider_axes1,:,:)))
            hold on
            Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes1,:,:,handles.peak_flow,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            set(himage, 'AlphaData', cdata);
%             Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes1,:,:,:));
%             himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%             cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%             cdata = double(cdata)*0.5;
%             set(himage, 'AlphaData', cdata);
            hold off
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray
            
        end
end
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


% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.rot == 1
    rotate3d(handles.axes1,'on')
    handles.rot = handles.rot +  1;
else
    rotate3d(handles.axes1,'off')
    handles.rot = handles.rot -  1;
end

handles.output = hObject;  
guidata(hObject, handles);


% --- Executes when GUIDE_LAPLACE is resized.
function GUIDE_LAPLACE_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to GUIDE_LAPLACE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.GUIDE_LAPLACE, 'Units', 'pixels');
    FigPos = get(handles.GUIDE_LAPLACE, 'Position');
    set(handles.GUIDE_LAPLACE, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.GUIDE_LAPLACE, 'Units', 'normalized');

    set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.37)
    set(handles.popupmenu2,'FontUnits','Normalized','FontSize',0.37)
    set(handles.popupmenu3,'FontUnits','Normalized','FontSize',0.37)
    set(handles.text1,'FontUnits','Normalized','FontSize',0.52)
    set(handles.text3,'FontUnits','Normalized','FontSize',0.52)
    set(handles.text4,'FontUnits','Normalized','FontSize',0.52)
    set(handles.text5,'FontUnits','Normalized','FontSize',0.52)
    set(handles.text6,'FontUnits','Normalized','FontSize',0.52)
    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton3,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton5,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton6,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton7,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton10,'FontUnits','Normalized','FontSize',0.24)
    set(handles.pushbutton11,'FontUnits','Normalized','FontSize',0.24)

handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb_method = get(handles.radiobutton1,'Value');

if rb_method == 1
    handles.werp = 2;
    set(handles.radiobutton2,'Value',0);
    set(handles.pushbutton3,'String','vWERP','FontUnits','Normalized','FontSize',0.24)
elseif rb_method == 0
    
    handles.werp = 1;
    set(handles.radiobutton2,'Value',1);
    set(handles.pushbutton3,'String','WERP','FontUnits','Normalized','FontSize',0.24)
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb_method = get(handles.radiobutton2,'Value');

if rb_method == 1
    handles.werp = 1;
    set(handles.radiobutton1,'Value',0);
    set(handles.pushbutton3,'String','WERP','FontUnits','Normalized','FontSize',0.24)
elseif rb_method == 0
    handles.werp = 2;
    set(handles.radiobutton1,'Value',1);
    set(handles.pushbutton3,'String','vWERP','FontUnits','Normalized','FontSize',0.24)
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb_volint = get(handles.radiobutton3,'Value');

if rb_volint == 1
    handles.volInt = true;
    set(handles.radiobutton4,'Value',0);
elseif rb_volint == 0
    
    handles.volInt = false;
    set(handles.radiobutton4,'Value',1);
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rb_volint = get(handles.radiobutton4,'Value');

if rb_volint == 1
    handles.volInt = false;
    set(handles.radiobutton3,'Value',0);
elseif rb_volint == 0
    
    handles.volInt = true;
    set(handles.radiobutton3,'Value',1);
end

handles.output = hObject;
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton4



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.resamp_factor = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);
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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.percentage_noise = str2double(get(hObject,'String'));

handles.output = hObject;
guidata(hObject, handles);
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


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.normal_outlet = handles.normal_outlet*-1;
    
    axes(handles.axes1);
    plot(0.0)
    patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    hold on
    plot3(handles.nodes(handles.id_outlet,1),handles.nodes(handles.id_outlet,2),handles.nodes(handles.id_outlet,3),'*y')
    center_out = mean(handles.nodes(handles.id_outlet,1:3));
    quiver3(center_out(1),center_out(2),center_out(3),handles.normal_outlet(1),handles.normal_outlet(2),handles.normal_outlet(3),30,'r','Linewidth',3)
    hold off
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])

    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    set(handles.popupmenu1,'value',1)
    set(handles.popupmenu3,'value',1)
    set(handles.pushbutton1,'visible','off')
    set(handles.pushbutton2,'visible','off')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    handles.id_mesh = 0;
    handles.id_mesh_vel = 0;
    handles.id_mesh_inlet = 0;
    handles.id_mesh_outlet = 1;
    handles.id_mesh_laplace = 0;
    handles.id_centerline = 0;
    handles.id_diameter = 0;
    handles.id_radius = 0;
    handles.id_axial_unit_vectors = 0;
    handles.id_circumferential_unit_vectors = 0;
    handles.id_wss_a = 0;
    handles.id_wss_c = 0;
    handles.id_axial_angle = 0;
    handles.id_forward_vel = 0;
    handles.id_backward_vel = 0;
    handles.id_regurgitant_flow = 0;
    handles.id_centerline_flow = 0;
    handles.id_eccentricity = 0;
    handles.id_curvature = 0; % new data Julio Sotelo
    handles.id_flattening = 0; % new data Julio Sotelo
    handles.id_ellipticity = 0; % new data Julio Sotelo
    handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
    handles.id_forward_vortex = 0; % new data Julio Sotelo
    handles.id_ipcmra = 0;
    handles.id_mag = 0;
            
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    handles.normal_inlet = handles.normal_inlet*-1;

    axes(handles.axes1);
    plot(0.0)
    patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5],'edgealpha',0.5)
    hold on
    plot3(handles.nodes(handles.id_inlet,1),handles.nodes(handles.id_inlet,2),handles.nodes(handles.id_inlet,3),'*y')
    center_in = mean(handles.nodes(handles.id_inlet,1:3));
    quiver3(center_in(1),center_in(2),center_in(3),handles.normal_inlet(1),handles.normal_inlet(2),handles.normal_inlet(3),30,'r','Linewidth',3)
    hold off
    axis vis3d
    lighting gouraud
    daspect([1,1,1])
    axis off
    view([-34,-51])


    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    set(handles.popupmenu1,'value',1)
    set(handles.popupmenu3,'value',1)
    set(handles.pushbutton1,'visible','off')
    set(handles.pushbutton2,'visible','off')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    handles.id_mesh = 0;
    handles.id_mesh_vel = 0;
    handles.id_mesh_inlet = 1;
    handles.id_mesh_outlet = 0;
    handles.id_mesh_laplace = 0;
    handles.id_centerline = 0;
    handles.id_diameter = 0;
    handles.id_radius = 0;
    handles.id_axial_unit_vectors = 0;
    handles.id_circumferential_unit_vectors = 0;
    handles.id_wss_a = 0;
    handles.id_wss_c = 0;
    handles.id_axial_angle = 0;
    handles.id_forward_vel = 0;
    handles.id_backward_vel = 0;
    handles.id_regurgitant_flow = 0;
    handles.id_centerline_flow = 0;
    handles.id_eccentricity = 0;
    handles.id_curvature = 0; % new data Julio Sotelo
    handles.id_flattening = 0; % new data Julio Sotelo
    handles.id_ellipticity = 0; % new data Julio Sotelo
    handles.id_length_vessel = 0; % new data Julio Sotelo
%         handles.id_circulation = 0; % new data Julio Sotelo
    handles.id_forward_vortex = 0; % new data Julio Sotelo
    handles.id_ipcmra = 0;
    handles.id_mag = 0;
            
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    inlet_seg   = handles.inlet_seg;
    n_inlet     = handles.n_inlet;
    inlet       = handles.inlet;
    outlet_seg  = handles.outlet_seg;
    n_outlet    = handles.n_outlet;
    outlet      = handles.outlet;
    mask        = handles.mask;
    v           = handles.v;
    v_old       = handles.v_old;
    time        = handles.time;
    DP_VWERP    = handles.DP_VWERP;
    pa2mmhg     = handles.pa2mmhg;
    VISC        = handles.VISC;
    ADVEC       = handles.ADVEC;
    KINETIC     = handles.KINETIC;
    LAMBDA      = handles.LAMBDA;
    w           = handles.w;
    Labels      = handles.Labels;

    directory = uigetdir(pwd, 'Select Directory');
    if directory~=0

        save([directory,'\inlet.mat'],'inlet_seg','n_inlet', 'inlet')
        save([directory,'\outlet.mat'],'outlet_seg','n_outlet', 'outlet')
        save([directory,'\mask.mat'],'mask')
        save([directory,'\v.mat'],'v','v_old')
        save([directory,'\results.mat'],'time','DP_VWERP','pa2mmhg','VISC', 'ADVEC', 'KINETIC', 'LAMBDA')
        save([directory,'\w.mat'],'w')
        save([directory,'\labels.mat'],'Labels')
        
        uiwait(msgbox("The data has been successfully saved","Success"));
        
    else

        uiwait(msgbox("The data was not saved","Error","error"));

    end


handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dilate

    % adjusting values for plot
    mask = handles.SEG_for_vwerp;
    
    % we erode the mask in one voxel
    se = strel('sphere',1);
    mask = imdilate(mask, se);

    % adjusting the velocities
        v = [];
        if max(abs(handles.MR_PCA_FH(:)))>10
            v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4))/100)*-1;
            v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4))/100)*-1;
            v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4))/100);
        else
            v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4)))*-1;
            v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4)))*-1;
            v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4)));
        end

    handles.v_for_plot = v;
    handles.mask_for_plot = mask;

    % we extract the voxelization using binsurface iso2mesh
    [node_bin,elem_bin] = binsurface(handles.mask_for_plot);


    axes(handles.axes1);
    plot(0.0)
    set(handles.uipanel2,'BackgroundColor',[1,1,1])


    %visu_mask_inlet_outlet(handles.v_for_plot, handles.mask_for_plot, handles.inlet_for_plot, handles.outlet_for_plot, handles.n_inlet_for_plot, handles.n_outlet_for_plot)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization of the velocities vector field %%%%%%%%%%%%%%%%
    % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % code
    time = 5;
    n = sum(sum(sum(handles.mask_for_plot)));
    pt = zeros(n,3);
    V = zeros(n,3);

    n = 0;
    for i = 1:size(handles.mask_for_plot,1)
       for j = 1:size(handles.mask_for_plot,2)
           for k = 1:size(handles.mask_for_plot,3)
               if(handles.mask_for_plot(i,j,k) ~= 1)
                   continue
               end

               n = n + 1;
               pt(n,:) = [i, j, k];
               V(n,:) = [handles.v_for_plot{1}.im(i,j,k,time), handles.v_for_plot{2}.im(i,j,k,time), handles.v_for_plot{3}.im(i,j,k,time)];
           end
       end
    end

    vs = 10;
    jump = 1;
    q = quiver3(pt(1:jump:end,1),pt(1:jump:end,2),pt(1:jump:end,3),vs*V(1:jump:end,1),vs*V(1:jump:end,2),vs*V(1:jump:end,3), 'Linewidth',1);
    mag_vel = sqrt(sum(V(1:jump:end,:).^2,2));
    currentColormap = colormap(handles.axes1,'jet');
    [~, ~, ind] = histcounts(mag_vel(:), size(currentColormap, 1));
    handles.cmap = uint8(ind2rgb(ind, currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(mag_vel(:));
    handles.max_vel = max(mag_vel(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [0 0 0]; % color
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 11;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [0 0 0];% color
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    axis equal
    q.AutoScale='off';

    hold on

    % binsurface plot
    trisurf(elem_bin,node_bin(:,1)+0.5,node_bin(:,2)+0.5,node_bin(:,3)+0.5, 'facecolor','none');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization of the velocities vector field %%%%%%%%%%%%%%%%
    % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel('A-P', 'Fontsize', 14,'FontWeight','bold')
    ylabel('F-H', 'Fontsize', 14,'FontWeight','bold')
    zlabel('R-L', 'Fontsize', 14,'FontWeight','bold')
    title('Vector Field','Fontsize', 16,'FontWeight','bold')
    hold off
    grid on
    daspect([1,1,1])
    view(3)

    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    set(handles.popupmenu1,'value',1)
    set(handles.popupmenu3,'value',1)
    

    handles.SEG_for_vwerp = mask;
    
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Erode

    % adjusting values for plot
    mask = handles.SEG_for_vwerp;
    
    % we erode the mask in one voxel
    se = strel('sphere',1);
    mask = imerode(mask, se);

    % adjusting the velocities
        v = [];
        if max(abs(handles.MR_PCA_FH(:)))>10
            v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4))/100)*-1;
            v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4))/100)*-1;
            v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4))/100);
        else
            v{3}.im = (handles.MR_PCA_RL.*repmat(mask,1,1,1,size(handles.MR_PCA_RL,4)))*-1;
            v{2}.im = (handles.MR_PCA_FH.*repmat(mask,1,1,1,size(handles.MR_PCA_FH,4)))*-1;
            v{1}.im = (handles.MR_PCA_AP.*repmat(mask,1,1,1,size(handles.MR_PCA_AP,4)));
        end

    handles.v_for_plot = v;
    handles.mask_for_plot = mask;

    % we extract the voxelization using binsurface iso2mesh
    [node_bin,elem_bin] = binsurface(handles.mask_for_plot);


    axes(handles.axes1);
    plot(0.0)
    set(handles.uipanel2,'BackgroundColor',[1,1,1])


    %visu_mask_inlet_outlet(handles.v_for_plot, handles.mask_for_plot, handles.inlet_for_plot, handles.outlet_for_plot, handles.n_inlet_for_plot, handles.n_outlet_for_plot)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization of the velocities vector field %%%%%%%%%%%%%%%%
    % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % code
    time = 5;
    n = sum(sum(sum(handles.mask_for_plot)));
    pt = zeros(n,3);
    V = zeros(n,3);

    n = 0;
    for i = 1:size(handles.mask_for_plot,1)
       for j = 1:size(handles.mask_for_plot,2)
           for k = 1:size(handles.mask_for_plot,3)
               if(handles.mask_for_plot(i,j,k) ~= 1)
                   continue
               end

               n = n + 1;
               pt(n,:) = [i, j, k];
               V(n,:) = [handles.v_for_plot{1}.im(i,j,k,time), handles.v_for_plot{2}.im(i,j,k,time), handles.v_for_plot{3}.im(i,j,k,time)];
           end
       end
    end

    vs = 10;
    jump = 1;
    q = quiver3(pt(1:jump:end,1),pt(1:jump:end,2),pt(1:jump:end,3),vs*V(1:jump:end,1),vs*V(1:jump:end,2),vs*V(1:jump:end,3), 'Linewidth',1);
    mag_vel = sqrt(sum(V(1:jump:end,:).^2,2));
    currentColormap = colormap(handles.axes1,'jet');
    [~, ~, ind] = histcounts(mag_vel(:), size(currentColormap, 1));
    handles.cmap = uint8(ind2rgb(ind, currentColormap) * 255);
    handles.cmap(:,:,4) = 255;
    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
    c = colorbar(handles.axes1);
    handles.min_vel = min(mag_vel(:));
    handles.max_vel = max(mag_vel(:));
    handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
    c.LimitsMode = 'manual';
    c.Limits = [handles.min_vel handles.max_vel];
    c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
    c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
    c.Color = [0 0 0]; % color
    c.FontWeight = 'bold';
    c.Label.String = 'Velocity [m/s]';
    c.FontSize = 11;
    c.Label.FontSize = 14;
    c.Label.FontWeight = 'bold';
    c.Label.Color = [0 0 0];% color
    caxis(handles.axes1, [handles.min_vel handles.max_vel]);
    axis equal
    q.AutoScale='off';

    hold on

    % binsurface plot
    trisurf(elem_bin,node_bin(:,1)+0.5,node_bin(:,2)+0.5,node_bin(:,3)+0.5, 'facecolor','none');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visualization of the velocities vector field %%%%%%%%%%%%%%%%
    % David Marlevi code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xlabel('A-P', 'Fontsize', 14,'FontWeight','bold')
    ylabel('F-H', 'Fontsize', 14,'FontWeight','bold')
    zlabel('R-L', 'Fontsize', 14,'FontWeight','bold')
    title('Vector Field','Fontsize', 16,'FontWeight','bold')
    hold off
    grid on
    daspect([1,1,1])
    view(3)

    set(handles.slider1,'visible','off')
    set(handles.text1,'visible','off')
    set(handles.popupmenu1,'value',1)
    set(handles.popupmenu3,'value',1)
    

    handles.SEG_for_vwerp = mask;
    
handles.output = hObject;
guidata(hObject, handles);
