function varargout = GUIDE_FLOW(varargin)
% GUIDE_FLOW MATLAB code for GUIDE_FLOW.fig
%      GUIDE_FLOW, by itself, creates a new GUIDE_FLOW or raises the existing
%      singleton*.
%
%      H = GUIDE_FLOW returns the handle to a new GUIDE_FLOW or the handle to
%      the existing singleton*.
%
%      GUIDE_FLOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_FLOW.M with the given input arguments.
%
%      GUIDE_FLOW('Property','Value',...) creates a new GUIDE_FLOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_FLOW_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_FLOW_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_FLOW

% Last Modified by GUIDE v2.5 19-Jun-2022 18:34:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_FLOW_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_FLOW_OutputFcn, ...
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


% --- Executes just before GUIDE_FLOW is made visible.
function GUIDE_FLOW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE_FLOW (see VARARGIN)

% Choose default command line output for GUIDE_FLOW
handles.output = hObject;


    path(path,'iso2mesh/')

    handles.SEG                 = varargin{1}.SEG;
    handles.IPCMRA              = varargin{1}.IPCMRA;
    handles.voxel_MR            = varargin{1}.voxel_MR;
    handles.L                   = varargin{1}.L;
    handles.Lrgb                = varargin{1}.Lrgb;
    handles.Lrgb_vel            = varargin{1}.Lrgb_vel;
    handles.NUM                 = varargin{1}.NUM;
    handles.xd                  = varargin{1}.xd;
    handles.yd                  = varargin{1}.yd;
    handles.zd                  = varargin{1}.zd;
    handles.a                   = varargin{1}.a;
    handles.b                   = varargin{1}.b;
    handles.c                   = varargin{1}.c;
    handles.d                   = varargin{1}.d;
    handles.slider_axes1        = varargin{1}.slider_id_axes1;
    handles.MR_FFE_FH           = varargin{1}.MR_FFE_FH;
    handles.MAG                 = mean(handles.MR_FFE_FH,4); 
    handles.MR_PCA_FH           = varargin{1}.MR_PCA_FH;
    handles.MR_PCA_AP           = varargin{1}.MR_PCA_AP;
    handles.MR_PCA_RL           = varargin{1}.MR_PCA_RL;
    handles.MR_PCA_FH_smooth    = varargin{1}.MR_PCA_FH_smooth;
    handles.MR_PCA_AP_smooth    = varargin{1}.MR_PCA_AP_smooth;
    handles.MR_PCA_RL_smooth    = varargin{1}.MR_PCA_RL_smooth;
    handles.heart_rate          = varargin{1}.heart_rate;
    handles.Lrgb_vel            = varargin{1}.Lrgb_vel;
    handles.mags_vel            = varargin{1}.mags_vel;
    handles.Lrgb                = varargin{1}.Lrgb;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % debo corregir el valor de slider axes, para que abra siempre el mismo
    % valor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read mesh
    handles.faces               = varargin{1}.faces;
    handles.nodes               = varargin{1}.nodes;
    handles.elem                = varargin{1}.elem;
    handles.veset               = varargin{1}.veset;
    handles.veset_out           = varargin{1}.veset_out;
    
    handles.veset(unique(handles.faces(:)),:,:) = handles.veset(unique(handles.faces(:)),:,:)*0;
    handles.nodes_id            = (1:size(handles.nodes,1))'; 
    handles.elem_id             = (1:size(handles.elem,1))';
    handles.faces_id            = (1:size(handles.faces,1))';
    handles.id_mesh_inlet_flow  = varargin{1}.id_mesh_inlet_flow;
    handles.id_image            = varargin{1}.id_image;    
    handles.id_view             = varargin{1}.id_view;
    handles.id_result           = varargin{1}.id_result;
    
    handles.id_selected_faceid  = varargin{1}.id_selected_faceid;    
    handles.faceid              = varargin{1}.faceid;    
    handles.cutpos              = varargin{1}.cutpos;    
    handles.time                = varargin{1}.time;    
    handles.peak_flow           = varargin{1}.peak_flow;    
    handles.flow                = varargin{1}.flow;    
    handles.net_flow            = varargin{1}.net_flow;    
    handles.max_velocity        = varargin{1}.max_velocity;    
    handles.min_velocity        = varargin{1}.min_velocity; 
    handles.velocity_proj       = varargin{1}.velocity_proj;
    handles.id_section          = varargin{1}.id_section;
    
    % id loop
    handles.id_while = 0;
    handles.id_move_peak_flow = 0;
    
    handles.list_string = {'Select View ...','Sagital View','Axial View','Coronal View'};
    set(handles.popupmenu1,'String',handles.list_string,'Value',handles.id_view+1);
    handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
    set(handles.popupmenu2,'String',handles.list_string,'Value',handles.id_image+1);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % adjust the image pointer
    if handles.id_mesh_inlet_flow == 1 && handles.id_result==1
        
        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','on');
        set(handles.slider1,'Visible','on');
        
        pond = squeeze(max(abs(handles.velocity_proj)))/max(squeeze(max(abs(handles.velocity_proj))))*5;
        
        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
        hold on
        patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
        q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                handles.veset(handles.id_section,1,handles.peak_flow),handles.veset(handles.id_section,2,handles.peak_flow),handles.veset(handles.id_section,3,handles.peak_flow),pond(handles.peak_flow),'LineWidth',1);
        hold off
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.velocity_proj(:));
        handles.max_vel = max(handles.velocity_proj(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.3];
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
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        % slider adjustment
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.peak_flow)])
        
        
    elseif handles.id_mesh_inlet_flow == 1 && handles.id_result==2
        
        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.flow,'-r','LineWidth',2)
        hold on 
        plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
        hold off
        text_1 = '\leftarrow Peak flow';
        text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
        txt = {text_1,text_2};
        text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.flow) max(handles.flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.XLabel.String = 'Time [s]';
        ax.YLabel.String = 'Flow [cm^{3}/s]';
        
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.flow) max(handles.flow)])
%         xlabel('Time [s]')
%         ylabel('Flow [cm^{3}/s]')
            
    elseif handles.id_mesh_inlet_flow == 1 && handles.id_result==3
        
        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.net_flow,'-b','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.net_flow) max(handles.net_flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.XLabel.String = 'Time [s]';
        ax.YLabel.String = 'Net Flow [ml]';
        
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.net_flow) max(handles.net_flow)])
%         xlabel('Time [s]')
%         ylabel('Net Flow [ml]')
            
    elseif handles.id_mesh_inlet_flow == 1 && handles.id_result==4
        
        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.max_velocity,'-g','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.max_velocity) max(handles.max_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.XLabel.String = 'Time [s]';
        ax.YLabel.String = 'Maximum Velocity [cm/s]';
        
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.max_velocity) max(handles.max_velocity)])
%         xlabel('Time [s]')
%         ylabel('Maximum Velocity [cm/s]')
        
    elseif handles.id_mesh_inlet_flow == 1 && handles.id_result==5
        
        
        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.min_velocity,'-k','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.min_velocity) max(handles.min_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.XLabel.String = 'Time [s]';
        ax.YLabel.String = 'Minimum Velocity [cm/s]';
        
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.min_velocity) max(handles.min_velocity)])
%         xlabel('Time [s]')
%         ylabel('Minimum Velocity [cm/s]')

    else
        
        if handles.id_image == 1

            if handles.id_view == 1

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

        elseif handles.id_image == 2

            if handles.id_view == 1

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
    end

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes GUIDE_FLOW wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_FLOW_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.id_while == 0
    
    varargout{1} = handles.output;

    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
    id_mesh_inlet_flow = handles.id_mesh_inlet_flow;
    setappdata(0,'id_mesh_inlet_flow',id_mesh_inlet_flow);
    id_view = handles.id_view;
    setappdata(0,'id_view',id_view);
    id_image = handles.id_image;
    setappdata(0,'id_image',id_image);
    id_result = handles.id_result;
    setappdata(0,'id_result',id_result);
    slider_id_axes1 = handles.slider_axes1;
    setappdata(0,'slider_id_axes1',slider_id_axes1);
    id_selected_faceid = handles.id_selected_faceid;
    setappdata(0,'id_selected_faceid',id_selected_faceid);
    faceid = handles.faceid;
    setappdata(0,'faceid',faceid);
    cutpos = handles.cutpos;
    setappdata(0,'cutpos',cutpos);
    time = handles.time;
    setappdata(0,'time',time);
    peak_flow = handles.peak_flow;
    setappdata(0,'peak_flow',peak_flow);
    flow = handles.flow;
    setappdata(0,'flow',flow);
    net_flow = handles.net_flow;
    setappdata(0,'net_flow',net_flow);
    max_velocity = handles.max_velocity;
    setappdata(0,'max_velocity',max_velocity);
    min_velocity = handles.min_velocity;
    setappdata(0,'min_velocity',min_velocity);
    velocity_proj = handles.velocity_proj;
    setappdata(0,'velocity_proj',velocity_proj);
    id_section = handles.id_section;
    setappdata(0,'id_section',id_section);
    veset_out = handles.veset_out;
    setappdata(0,'veset_out',veset_out);
    Lrgb_vel = handles.Lrgb_vel;
    setappdata(0,'Lrgb_vel',Lrgb_vel);
    mags_vel = handles.mags_vel;
    setappdata(0,'mags_vel',mags_vel);
    MR_PCA_FH = handles.MR_PCA_FH;
    setappdata(0,'MR_PCA_FH',MR_PCA_FH);
    MR_PCA_AP = handles.MR_PCA_AP;
    setappdata(0,'MR_PCA_AP',MR_PCA_AP);
    MR_PCA_RL = handles.MR_PCA_RL;
    setappdata(0,'MR_PCA_RL',MR_PCA_RL);
    MR_PCA_FH_smooth = handles.MR_PCA_FH_smooth;
    setappdata(0,'MR_PCA_FH_smooth',MR_PCA_FH_smooth);
    MR_PCA_AP_smooth = handles.MR_PCA_AP_smooth;
    setappdata(0,'MR_PCA_AP_smooth',MR_PCA_AP_smooth);
    MR_PCA_RL_smooth = handles.MR_PCA_RL_smooth;
    setappdata(0,'MR_PCA_RL_smooth',MR_PCA_RL_smooth);


elseif handles.id_while == 1
    
    varargout{1} = handles.output;

    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
    id_mesh_inlet_flow = handles.id_mesh_inlet_flow;
    setappdata(0,'id_mesh_inlet_flow',id_mesh_inlet_flow);
    id_view = handles.id_view;
    setappdata(0,'id_view',id_view);
    id_image = handles.id_image;
    setappdata(0,'id_image',id_image);
    id_result = handles.id_result;
    setappdata(0,'id_result',id_result);
    slider_id_axes1 = handles.slider_axes1;
    setappdata(0,'slider_id_axes1',slider_id_axes1);
    id_selected_faceid = handles.id_selected_faceid;
    setappdata(0,'id_selected_faceid',id_selected_faceid);
    faceid = handles.faceid;
    setappdata(0,'faceid',faceid);
    cutpos = handles.cutpos;
    setappdata(0,'cutpos',cutpos);
    time = handles.time;
    setappdata(0,'time',time);
    peak_flow = handles.peak_flow;
    setappdata(0,'peak_flow',peak_flow);
    flow = handles.flow;
    setappdata(0,'flow',flow);
    net_flow = handles.net_flow;
    setappdata(0,'net_flow',net_flow);
    max_velocity = handles.max_velocity;
    setappdata(0,'max_velocity',max_velocity);
    min_velocity = handles.min_velocity;
    setappdata(0,'min_velocity',min_velocity);
    velocity_proj = handles.velocity_proj;
    setappdata(0,'velocity_proj',velocity_proj);
    id_section = handles.id_section;
    setappdata(0,'id_section',id_section);
    veset_out = handles.veset_out;
    setappdata(0,'veset_out',veset_out);
    Lrgb_vel = handles.Lrgb_vel;
    setappdata(0,'Lrgb_vel',Lrgb_vel);
    mags_vel = handles.mags_vel;
    setappdata(0,'mags_vel',mags_vel);
    MR_PCA_FH = handles.MR_PCA_FH;
    setappdata(0,'MR_PCA_FH',MR_PCA_FH);
    MR_PCA_AP = handles.MR_PCA_AP;
    setappdata(0,'MR_PCA_AP',MR_PCA_AP);
    MR_PCA_RL = handles.MR_PCA_RL;
    setappdata(0,'MR_PCA_RL',MR_PCA_RL);
    MR_PCA_FH_smooth = handles.MR_PCA_FH_smooth;
    setappdata(0,'MR_PCA_FH_smooth',MR_PCA_FH_smooth);
    MR_PCA_AP_smooth = handles.MR_PCA_AP_smooth;
    setappdata(0,'MR_PCA_AP_smooth',MR_PCA_AP_smooth);
    MR_PCA_RL_smooth = handles.MR_PCA_RL_smooth;
    setappdata(0,'MR_PCA_RL_smooth',MR_PCA_RL_smooth);
    delete(handles.figure1);
end

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.popupmenu1,'Value')
    
    case 1
        
        set(handles.popupmenu1,'Value',1)
        set(handles.popupmenu2,'Value',1)
        set(handles.popupmenu3,'Value',1)
        set(handles.slider1,'Visible','off')
        axes(handles.axes1);
        plot(0.0)
        axis off
        set(handles.pushbutton1,'Visible','off')
        set(handles.text1,'Visible','off')
         
    case 2
        
        set(handles.popupmenu3,'Value',1)
        handles.id_view = 1;
        handles.id_image = get(handles.popupmenu2,'Value')-1;
        
        
        if handles.id_image == 0

            axes(handles.axes1);
            plot(0.0)
            axis off
            
            set(handles.slider1,'Visible','off')
            set(handles.pushbutton1,'Visible','off')
            set(handles.text1,'Visible','off')
            
        else
            
            set(handles.slider1,'Visible','on')
            set(handles.pushbutton1,'Visible','on')
            set(handles.text1,'Visible','on')

        end
        
        handles.slider_axes1 = round(handles.c/2);
        
        if handles.id_image == 1

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


        elseif handles.id_image == 2


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

        end

    case 3
        
        set(handles.popupmenu3,'Value',1)
        handles.id_view = 2;
        handles.id_image = get(handles.popupmenu2,'Value')-1;
        
        if handles.id_image == 0

            axes(handles.axes1);
            plot(0.0)
            axis off
            
            set(handles.slider1,'Visible','off')
            set(handles.pushbutton1,'Visible','off')
            set(handles.text1,'Visible','off')
            
        else
            
            set(handles.slider1,'Visible','on')
            set(handles.pushbutton1,'Visible','on')
            set(handles.text1,'Visible','on')

        end

        handles.slider_axes1 = round(handles.b/2);
        
        if handles.id_image == 1

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


        elseif handles.id_image == 2

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

        end

    case 4
        set(handles.popupmenu3,'Value',1)
        handles.id_view = 3;
        handles.id_image = get(handles.popupmenu2,'Value')-1;
        
        if handles.id_image == 0

            axes(handles.axes1);
            plot(0.0)
            axis off
            
            set(handles.slider1,'Visible','off')
            set(handles.pushbutton1,'Visible','off')
            set(handles.text1,'Visible','off')
            
        else
            
            set(handles.slider1,'Visible','on')
            set(handles.pushbutton1,'Visible','on')
            set(handles.text1,'Visible','on')

        end
        
        handles.slider_axes1 = round(handles.a/2);
        
        if handles.id_image == 1

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


        elseif handles.id_image == 2

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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.list_string = {'Select Image ...','IPCMRA','MAGNITUDE'};
% set(handles.popupmenu2,'String',handles.list_string,'Value',handles.id_image);
% handles.list_string = {'Select Image ...','Sagital View','Axial View','Coronal View'};
% set(handles.popupmenu1,'String',handles.list_string,'Value',handles.id_view);

switch get(handles.popupmenu2,'Value')
    
    case 1
        
        set(handles.popupmenu1,'Value',1)
        set(handles.popupmenu2,'Value',1)
        set(handles.popupmenu3,'Value',1)
        set(handles.slider1,'Visible','off')
        axes(handles.axes1);
        plot(0.0)
        axis off
        set(handles.pushbutton1,'Visible','off')
        set(handles.text1,'Visible','off')
         
    case 2
        set(handles.popupmenu3,'Value',1)
        handles.id_image = 1;
        handles.id_view = get(handles.popupmenu1,'Value')-1;
        
        
        if handles.id_view == 0

            axes(handles.axes1);
            plot(0.0)
            axis off
            
            set(handles.slider1,'Visible','off')
            set(handles.pushbutton1,'Visible','off')
            set(handles.text1,'Visible','off')
            
        else
            
            set(handles.slider1,'Visible','on')
            set(handles.pushbutton1,'Visible','on')
            set(handles.text1,'Visible','on')

        end
        


        if handles.id_view == 1

            set(handles.pushbutton1,'Visible','on')
            
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

            set(handles.pushbutton1,'Visible','on')
            
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

            set(handles.pushbutton1,'Visible','on')
            
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

    case 3
        
        set(handles.popupmenu3,'Value',1)
        handles.id_image = 2;
        handles.id_view = get(handles.popupmenu1,'Value')-1;
        
         if handles.id_view == 0

            axes(handles.axes1);
            plot(0.0)
            axis off
            
            set(handles.slider1,'Visible','off')
            set(handles.pushbutton1,'Visible','off')
            set(handles.text1,'Visible','off')
            
        else
            
            set(handles.slider1,'Visible','on')
            set(handles.pushbutton1,'Visible','on')
            set(handles.text1,'Visible','on')

         end

        if handles.id_view == 1
            
            set(handles.pushbutton1,'Visible','on')
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

        elseif handles.id_view == 2

            set(handles.pushbutton1,'Visible','on')
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

        elseif handles.id_view == 3

            set(handles.pushbutton1,'Visible','on')
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

        end

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


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
switch get(handles.popupmenu3,'Value')
    
    case 1
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        axes(handles.axes1);
        plot(0.0)
        axis off
        
    case 2
        
        handles.id_result = 1;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','on');
        set(handles.slider1,'Visible','on');

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
        hold on
        patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
        q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                handles.veset(handles.id_section,1,handles.peak_flow),handles.veset(handles.id_section,2,handles.peak_flow),handles.veset(handles.id_section,3,handles.peak_flow),5,'LineWidth',1);
        hold off
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.velocity_proj(:));
        handles.max_vel = max(handles.velocity_proj(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.3];
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
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        

        % slider adjustment
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.peak_flow)])
        
            
    case 3

        handles.id_result = 2;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');

        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.flow,'-r','LineWidth',2)
        hold on 
        plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
        hold off
        text_1 = '\leftarrow Peak flow';
        text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
        txt = {text_1,text_2};
        text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.flow) max(handles.flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Flow [cm^{3}/s]';
        ax.XLabel.String = 'Time [s]';
        
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.flow) max(handles.flow)])
%         ax = gca;
%         ax.XColor = 'w';
%         ax.YColor = 'w';
%         ax.FontWeight = 'bold';
%         ax.FontSize = 10;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridColor = 'k';
%         ax.YLabel.String = 'Flow [cm^{3}/s]';
%         ax.XLabel.String = 'Time [s]';


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        answer1 = questdlg('¿Do you want to change the position of peak systole?','Question','Yes','No','No');
        switch answer1
            case 'Yes'
                
                set(handles.slider1,'visible','on')
%                 set(handles.pushbutton1,'visible','on','string','SAVE NEW PEAK FLOW')
                set(handles.text1,'visible','on','string',['# Cardiac Phase: ',num2str(handles.peak_flow)])
                handles.id_move_peak_flow = 1;
                
            case 'No'
                
%                 answer = questdlg(  '¿Do you want to average the peak systolic cardiac phases?','Question','Yes','No','No');
%                 switch answer
%                     case 'Yes'
%                         handles.id_average = 1;
%                         prompt = {'How many cardiac phases before peak systole:','How many cardiac phases after peak systole:'};
%                         dlgtitle = 'How many cardiac phases';
%                         definput = {'2','1'};
%                         dims = [1 35];
%                         answer = inputdlg(prompt,dlgtitle,dims,definput);
% 
%                         if handles.peak_flow - str2num(answer{1}) <1
% 
%                             f = warndlg('The number of cardiac phases before peak systole exceeds the limits','Warning');
% 
%                         elseif handles.peak_flow + str2num(answer{2})>length(handles.time)
% 
%                             f = warndlg('The number of cardiac phases after peak systole exceeds the limits','Warning');
% 
%                         else
%                             shad = [handles.time(handles.peak_flow - str2num(answer{1})),min(handles.flow);...
%                                 handles.time(handles.peak_flow + str2num(answer{2})),min(handles.flow);...
%                                 handles.time(handles.peak_flow + str2num(answer{2})),max(handles.flow);...
%                                 handles.time(handles.peak_flow - str2num(answer{1})),max(handles.flow)];
% 
%                             axes(handles.axes1);
%                             plot(0.0)
%                             plot(handles.time,handles.flow,'-r','LineWidth',2)
%                             hold on 
%                             plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
%                             fill(shad(:,1),shad(:,2),'g','edgecolor','none','facealpha',0.4)
%                             hold off
%                             text_1 = '\leftarrow Peak flow';
%                             text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
%                             txt = {text_1,text_2};
%                             text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
%                             axis square
%                             ax = gca;
%                             ax.XColor = 'w';
%                             ax.YColor = 'w';
%                             ax.XLim = [handles.time(1) handles.time(end)];
%                             ax.YLim = [min(handles.flow) max(handles.flow)];
%                             ax.FontWeight = 'bold';
%                             ax.FontSize = FontS;
%                             ax.XGrid = 'on';
%                             ax.YGrid = 'on';
%                             ax.GridColor = 'k';
%                             ax.YLabel.String = 'Flow [cm^{3}/s]';
%                             ax.XLabel.String = 'Time [s]';
% 
%                             answer1 = questdlg(  '¿Do you agree with the cardiac phases selected?','Question','Yes','No','No');
%                             switch answer1
%                                 case 'Yes'
% 
%                                     h = waitbar(0,['Averaging velocities...']);
%                                     pause(1)
% 
%                                     cp = handles.peak_flow - str2num(answer{1}):handles.peak_flow + str2num(answer{2});
%                                     veset_av = mean(handles.veset(:,:,cp),3);
%                                     handles.veset_out(:,:,cp) = repmat(veset_av,1,1,length(cp));
% 
%                                     handles.MR_PCA_FH(:,:,:,cp) = repmat(mean(handles.MR_PCA_FH(:,:,:,cp),4),1,1,1,length(cp));
%                                     handles.MR_PCA_AP(:,:,:,cp) = repmat(mean(handles.MR_PCA_AP(:,:,:,cp),4),1,1,1,length(cp));
%                                     handles.MR_PCA_RL(:,:,:,cp) = repmat(mean(handles.MR_PCA_RL(:,:,:,cp),4),1,1,1,length(cp));
%                                     handles.MR_PCA_FH_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_FH_smooth(:,:,:,cp),4),1,1,1,length(cp));
%                                     handles.MR_PCA_AP_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_AP_smooth(:,:,:,cp),4),1,1,1,length(cp));
%                                     handles.MR_PCA_RL_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_RL_smooth(:,:,:,cp),4),1,1,1,length(cp));
% 
% 
%                                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     % matrix color
%                                     handles.mags_vel = squeeze(sqrt(sum(handles.veset_out.^2,2)));
%                                     handles.Lrgb_vel = ones(size(handles.MR_FFE_FH,1)+4,size(handles.MR_FFE_FH,2)+4,size(handles.MR_FFE_FH,3)+4,size(handles.MR_FFE_FH,4),3);
% 
% 
%                                     MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[2,1,3]);
%                                     MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
%                                     xd_seg = MASK.*handles.xd;
%                                     yd_seg = MASK.*handles.yd;
%                                     zd_seg = MASK.*handles.zd;
%                                     xd_seg(MASK==0) = [];
%                                     yd_seg(MASK==0) = [];
%                                     zd_seg(MASK==0) = [];
%                                     pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
%                                     [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,3));
%                                     X_seg = MASK.*X;
%                                     Y_seg = MASK.*Y;
%                                     Z_seg = MASK.*Z;
%                                     X_seg(MASK==0) = [];
%                                     Y_seg(MASK==0) = [];
%                                     Z_seg(MASK==0) = [];
% 
%                                     pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
%                                     [~, ~, ind] = histcounts(handles.mags_vel(:), size(jet(128), 1));
%                                     ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
%                                     cmap_vol = uint8(ind2rgb(ind, jet(128)) * 255);
%                                     rgb_vel = cmap_vol/255;
% 
%                                     for n=1:length(pos_voxel(:,1))
%                                         d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
%                                         handles.Lrgb_vel(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
%                                     end
% 
%                                     handles.Lrgb_vel = handles.Lrgb_vel(3:end-2,3:end-2,3:end-2,:,:);
% 
%                                     close(h)
% 
%                                 case 'No'
% 
%                             end
%                         end

%                     case 'No'
%                 end
                
        end

    case 4
            
        handles.id_result = 3;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.net_flow,'-b','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.net_flow) max(handles.net_flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Net Flow [ml]';
        ax.XLabel.String = 'Time [s]';
%         grid on         
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.net_flow) max(handles.net_flow)])        
%         ax = gca;
%         ax.XColor = 'w';
%         ax.YColor = 'w';
%         ax.FontWeight = 'bold';
%         ax.FontSize = 10;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridColor = 'k';
%         ax.YLabel.String = 'Net Flow [ml]';
%         ax.XLabel.String = 'Time [s]';
            
    case 5
            
        handles.id_result = 4;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.max_velocity,'-g','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.max_velocity) max(handles.max_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Maximum Velocity [cm/s]';
        ax.XLabel.String = 'Time [s]';
        
        
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.max_velocity) max(handles.max_velocity)])
%         ax = gca;
%         ax.XColor = 'w';
%         ax.YColor = 'w';
%         ax.FontWeight = 'bold';
%         ax.FontSize = 10;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridColor = 'k';
%         ax.YLabel.String = 'Maximum Velocity [cm/s]';
%         ax.XLabel.String = 'Time [s]';
            
    case 6
            
        handles.id_result = 5;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.min_velocity,'-k','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.min_velocity) max(handles.min_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Minimum Velocity [cm/s]';
        ax.XLabel.String = 'Time [s]';
%         axis square
%         grid on 
%         xlim([handles.time(1) handles.time(end)])
%         ylim([min(handles.min_velocity) max(handles.min_velocity)])
%         ax = gca;
%         ax.XColor = 'w';
%         ax.YColor = 'w';
%         ax.FontWeight = 'bold';
%         ax.FontSize = 10;
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.GridColor = 'k';
%         ax.YLabel.String = 'Minimum Velocity [cm/s]';
%         ax.XLabel.String = 'Time [s]';
            
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


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.id_move_peak_flow == 1
    
    
        pp=1/handles.d;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = handles.d;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.slider_axes_peak = handles.slider_value;
        handles.peak_flow = handles.slider_axes_peak;
%         handles.id_result = 2;
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);

        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26,'string','SAVE PEAK FLOW');
        set(handles.text1,'visible','on','string',['# Cardiac Phase: ',num2str(handles.peak_flow)])
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.flow,'-r','LineWidth',2)
        hold on 
        plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
        hold off
        text_1 = '\leftarrow Peak flow';
        text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
        txt = {text_1,text_2};
        text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.flow) max(handles.flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Flow [cm^{3}/s]';
        ax.XLabel.String = 'Time [s]';
 
else
    
    id = get(handles.popupmenu3,'Value');
    if id == 2
        if handles.id_mesh_inlet_flow == 1 && handles.id_result==1

            pp=1/handles.d;
            slider_step(1) = pp;
            slider_step(2) = 0.1;
            set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
            maxslice = handles.d;
            handles.slider_value = get(hObject,'Value');
            if handles.slider_value<=pp
                    handles.slider_value = 1;
            else
                    handles.slider_value = round(handles.slider_value*maxslice);
            end
            handles.slider_axes_peak = handles.slider_value;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handles.id_view = 0;
            handles.id_image = 0;
            set(handles.popupmenu1,'Value',1);
            set(handles.popupmenu2,'Value',1);
            set(handles.pushbutton1,'Visible','off');
            set(handles.text1,'Visible','on');
            set(handles.slider1,'Visible','on');


            pond = squeeze(max(abs(handles.velocity_proj)))/max(squeeze(max(abs(handles.velocity_proj))))*5;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
            hold on
            patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
            q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                    handles.veset(handles.id_section,1,handles.slider_axes_peak),handles.veset(handles.id_section,2,handles.slider_axes_peak),handles.veset(handles.id_section,3,handles.slider_axes_peak),pond(handles.slider_axes_peak),'LineWidth',1);
            hold off
            currentColormap = colormap(handles.axes1,'jet');
            [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,handles.slider_axes_peak), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vel = min(handles.velocity_proj(:));
            handles.max_vel = max(handles.velocity_proj(:));
            handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vel handles.max_vel];
            c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
            c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.3];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_vel handles.max_vel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])

            % slider adjustment
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes_peak/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.slider_axes_peak)])

        end
    end

    if handles.id_image == 1
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


        end


    elseif handles.id_image == 2

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

    end

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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% if handles.id_move_peak_flow == 1
%     
%     handles.id_move_peak_flow = 0;
%     set(handles.popupmenu1,'Value',1);
%     set(handles.popupmenu2,'Value',1);
% 
%     set(handles.pushbutton1, 'FontUnits', 'points');
%     FontS = get(handles.pushbutton1, 'FontSize');
%     set(handles.pushbutton1, 'Visible','off','FontUnits', 'normalized','FontSize',0.26,'string','SELECT PLANE');
%     set(handles.slider1,'visible','off')
%     set(handles.text1,'visible','off','string',['# Slice: ',num2str(handles.slider_axes1)])
%     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     answer = questdlg(  '¿Do you want to average the peak systolic cardiac phases?','Question','Yes','No','No');
% %     switch answer
% %         case 'Yes'
% %             handles.id_average = 1;
% %             prompt = {'How many cardiac phases before peak systole:','How many cardiac phases after peak systole:'};
% %             dlgtitle = 'How many cardiac phases';
% %             definput = {'2','1'};
% %             dims = [1 35];
% %             answer = inputdlg(prompt,dlgtitle,dims,definput);
% % 
% %             if handles.peak_flow - str2num(answer{1}) <1
% % 
% %                 f = warndlg('The number of cardiac phases before peak systole exceeds the limits','Warning');
% % 
% %             elseif handles.peak_flow + str2num(answer{2})>length(handles.time)
% % 
% %                 f = warndlg('The number of cardiac phases after peak systole exceeds the limits','Warning');
% % 
% %             else
% %                 shad = [handles.time(handles.peak_flow - str2num(answer{1})),min(handles.flow);...
% %                     handles.time(handles.peak_flow + str2num(answer{2})),min(handles.flow);...
% %                     handles.time(handles.peak_flow + str2num(answer{2})),max(handles.flow);...
% %                     handles.time(handles.peak_flow - str2num(answer{1})),max(handles.flow)];
% % 
% %                 axes(handles.axes1);
% %                 plot(0.0)
% %                 plot(handles.time,handles.flow,'-r','LineWidth',2)
% %                 hold on 
% %                 plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
% %                 fill(shad(:,1),shad(:,2),'g','edgecolor','none','facealpha',0.4)
% %                 hold off
% %                 text_1 = '\leftarrow Peak flow';
% %                 text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
% %                 txt = {text_1,text_2};
% %                 text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
% %                 axis square
% %                 ax = gca;
% %                 ax.XColor = 'w';
% %                 ax.YColor = 'w';
% %                 ax.XLim = [handles.time(1) handles.time(end)];
% %                 ax.YLim = [min(handles.flow) max(handles.flow)];
% %                 ax.FontWeight = 'bold';
% %                 ax.FontSize = FontS;
% %                 ax.XGrid = 'on';
% %                 ax.YGrid = 'on';
% %                 ax.GridColor = 'k';
% %                 ax.YLabel.String = 'Flow [cm^{3}/s]';
% %                 ax.XLabel.String = 'Time [s]';
% % 
% %                 answer1 = questdlg(  '¿Do you agree with the cardiac phases selected?','Question','Yes','No','No');
% %                 switch answer1
% %                     case 'Yes'
% % 
% %                         h = waitbar(0,['Averaging velocities...']);
% %                         pause(1)
% % 
% %                         cp = handles.peak_flow - str2num(answer{1}):handles.peak_flow + str2num(answer{2});
% %                         veset_av = mean(handles.veset(:,:,cp),3);
% %                         handles.veset_out(:,:,cp) = repmat(veset_av,1,1,length(cp));
% % 
% %                         handles.MR_PCA_FH(:,:,:,cp) = repmat(mean(handles.MR_PCA_FH(:,:,:,cp),4),1,1,1,length(cp));
% %                         handles.MR_PCA_AP(:,:,:,cp) = repmat(mean(handles.MR_PCA_AP(:,:,:,cp),4),1,1,1,length(cp));
% %                         handles.MR_PCA_RL(:,:,:,cp) = repmat(mean(handles.MR_PCA_RL(:,:,:,cp),4),1,1,1,length(cp));
% %                         handles.MR_PCA_FH_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_FH_smooth(:,:,:,cp),4),1,1,1,length(cp));
% %                         handles.MR_PCA_AP_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_AP_smooth(:,:,:,cp),4),1,1,1,length(cp));
% %                         handles.MR_PCA_RL_smooth(:,:,:,cp) = repmat(mean(handles.MR_PCA_RL_smooth(:,:,:,cp),4),1,1,1,length(cp));
% % 
% % 
% %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                         % matrix color
% %                         handles.mags_vel = squeeze(sqrt(sum(handles.veset_out.^2,2)));
% %                         handles.Lrgb_vel = ones(size(handles.MR_FFE_FH,1)+4,size(handles.MR_FFE_FH,2)+4,size(handles.MR_FFE_FH,3)+4,size(handles.MR_FFE_FH,4),3);
% % 
% % 
% %                         MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[2,1,3]);
% %                         MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
% %                         xd_seg = MASK.*handles.xd;
% %                         yd_seg = MASK.*handles.yd;
% %                         zd_seg = MASK.*handles.zd;
% %                         xd_seg(MASK==0) = [];
% %                         yd_seg(MASK==0) = [];
% %                         zd_seg(MASK==0) = [];
% %                         pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
% %                         [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,3));
% %                         X_seg = MASK.*X;
% %                         Y_seg = MASK.*Y;
% %                         Z_seg = MASK.*Z;
% %                         X_seg(MASK==0) = [];
% %                         Y_seg(MASK==0) = [];
% %                         Z_seg(MASK==0) = [];
% % 
% %                         pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
% %                         [~, ~, ind] = histcounts(handles.mags_vel(:), size(jet(128), 1));
% %                         ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
% %                         cmap_vol = uint8(ind2rgb(ind, jet(128)) * 255);
% %                         rgb_vel = cmap_vol/255;
% % 
% %                         for n=1:length(pos_voxel(:,1))
% %                             d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
% %                             handles.Lrgb_vel(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_vel(d<mean(handles.voxel_MR)*2,:,:),1);
% %                         end
% % 
% %                         handles.Lrgb_vel = handles.Lrgb_vel(3:end-2,3:end-2,3:end-2,:,:);
% % 
% %                         close(h)
% % 
% %                     case 'No'
% % 
% %                 end
% %             end
% % 
% %         case 'No'
% %     end
% 
% else

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

%         center_ori = [center(2),center(1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
%         point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),mean(mean(handles.zd(:,:,handles.slider_axes1)))];
%         point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

        center_ori = [center(2),center(1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point1 = [PUNTOS_C(1,2),PUNTOS_C(1,1),(handles.slider_axes1-1)*handles.voxel_MR(3)];
        point2 = [center_ori(1),center_ori(2), center_ori(3) + sqrt(sum((center_ori-point1).^2))];

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure, % julio
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
%         hold on
%         plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
%         plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
%         plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
%         hold off
%         daspect([1 1 1])
%         axis vis3d
%         xlabel('ejex')
%         ylabel('ejey')
%         zlabel('ejez')
%         
        
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

        handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        id_selected_element = zeros(size(Selected_element(:)));
        for n=1:length(r)
            [rn,~,~] = find(Selected_element(:)==r(n)); 
            id_selected_element(rn) = 1;
        end
        id_selected_element = reshape(id_selected_element,[size(Selected_element,1),size(Selected_element,2)]);
        [rn,~,~] = find(sum(id_selected_element,2)==4);
        handles.id_selected_faceid = rn;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % node volume
        nodevol = nodevolume(handles.nodes,handles.elem);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % positive direction of the velocity verification.
        v1 = negative_vector/norm(negative_vector);
        v2 = VR_n/norm(VR_n);

        angle_dir = acos(sum(v1.*v2))*180/pi;
        if angle_dir<90
            VR_n = VR_n*-1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyection of the velocity
        normal = VR_n/norm(VR_n);
        handles.velocity_proj = sum(handles.veset(r,:,:).*repmat(normal,size(handles.veset(r,:,:),1),1,size(handles.veset(r,:,:),3)),2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow quantification    
        nod_volume_selected = nodevol(r)/sum(nodevol(r));
        handles.flow = squeeze(sum((handles.velocity_proj*100).*nod_volume_selected*(area/100)));
        handles.time = linspace(0,60/handles.heart_rate,size(handles.veset,3))';
        [~,I] = max(abs(handles.flow));
        handles.peak_flow = I;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % maximum and minimum velocity
        handles.max_velocity = squeeze(max(handles.velocity_proj*100));
        handles.min_velocity = squeeze(min(handles.velocity_proj*100));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % netflow
        handles.net_flow = cumsum(handles.flow.*handles.time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % show results in the popupmenu3
        handles.list_string = {'Select Results ...','Mesh + Slice','Flow','Net Flow','Maximum velocity','Minimum velocity'};
        set(handles.popupmenu3,'Visible','on','String',handles.list_string,'Value',2);

        pond = squeeze(max(abs(handles.velocity_proj)))/max(squeeze(max(abs(handles.velocity_proj))))*5;

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
        hold on
        patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
        q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                handles.veset(handles.id_section,1,handles.peak_flow),handles.veset(handles.id_section,2,handles.peak_flow),handles.veset(handles.id_section,3,handles.peak_flow),pond(handles.peak_flow),'LineWidth',1);
        hold off
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.velocity_proj(:));
        handles.max_vel = max(handles.velocity_proj(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.3];
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
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])


        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.slider1,'visible','on')
        set(handles.text1,'visible','on')
        set(handles.pushbutton1,'visible','off')

        % slider adjustment
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.peak_flow)])

        handles.id_mesh_inlet_flow = 1;
        handles.id_result = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(98/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        close(h)

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

       
%         center_ori = [center(3),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),center(2)]
%         point1 = [PUNTOS_C(1,2),mean(mean(squeeze(handles.xd(:,handles.slider_axes1,:)))),PUNTOS_C(1,1)];
%         point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];

        center_ori = [center(3),(handles.slider_axes1-1)*handles.voxel_MR(2),center(2)];
        point1 = [PUNTOS_C(1,2),(handles.slider_axes1-1)*handles.voxel_MR(2),PUNTOS_C(1,1)];
        point2 = [center_ori(1),center_ori(2) + sqrt(sum((center_ori-point1).^2)), center_ori(3)];
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure, % julio
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
%         hold on
%         plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
%         plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
%         plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
%         hold off
%         daspect([1 1 1])
%         axis vis3d
%         xlabel('ejex')
%         ylabel('ejey')
%         zlabel('ejez')

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

        handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        id_selected_element = zeros(size(Selected_element(:)));
        for n=1:length(r)
            [rn,~,~] = find(Selected_element(:)==r(n)); 
            id_selected_element(rn) = 1;
        end
        id_selected_element = reshape(id_selected_element,[size(Selected_element,1),size(Selected_element,2)]);
        [rn,~,~] = find(sum(id_selected_element,2)==4);
        handles.id_selected_faceid = rn;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % node volume
        nodevol = nodevolume(handles.nodes,handles.elem);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % positive direction of the velocity verification.
        v1 = negative_vector/norm(negative_vector);
        v2 = VR_n/norm(VR_n);

        angle_dir = acos(sum(v1.*v2))*180/pi;
        if angle_dir<90
            VR_n = VR_n*-1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyection of the velocity
        normal = VR_n/norm(VR_n);
        handles.velocity_proj = sum(handles.veset(r,:,:).*repmat(normal,size(handles.veset(r,:,:),1),1,size(handles.veset(r,:,:),3)),2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow quantification
        nod_volume_selected = nodevol(r)/sum(nodevol(r));
        handles.flow = squeeze(sum((handles.velocity_proj*100).*nod_volume_selected*(area/100)));
        handles.time = linspace(0,60/handles.heart_rate,size(handles.veset,3))';
        [~,I] = max(abs(handles.flow));
        handles.peak_flow = I;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % maximum and minimum velocity
        handles.max_velocity = squeeze(max(handles.velocity_proj*100));
        handles.min_velocity = squeeze(min(handles.velocity_proj*100));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % netflow
        handles.net_flow = cumsum(handles.flow.*handles.time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % show results in the popupmenu3
        handles.list_string = {'Select Results ...','Mesh + Slice','Flow','Net Flow','Maximum velocity','Minimum velocity'};
        set(handles.popupmenu3,'Visible','on','String',handles.list_string,'Value',2);

        pond = squeeze(max(abs(handles.velocity_proj)))/max(squeeze(max(abs(handles.velocity_proj))))*5;

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
        hold on
        patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
        q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                handles.veset(handles.id_section,1,handles.peak_flow),handles.veset(handles.id_section,2,handles.peak_flow),handles.veset(handles.id_section,3,handles.peak_flow),pond(handles.peak_flow),'LineWidth',1);
        hold off
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.velocity_proj(:));
        handles.max_vel = max(handles.velocity_proj(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.3];
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
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.slider1,'visible','on')
        set(handles.text1,'visible','on')
        set(handles.pushbutton1,'visible','off')

        % slider adjustment
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.peak_flow)])

        handles.id_mesh_inlet_flow = 1;
        handles.id_result = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(98/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        close( h)

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


%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure, % julio
%         patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor','k')
%         hold on
%         plot3(center_ori(1),center_ori(2),center_ori(3),'*g','Linewidth',10)
%         plot3(point1(1),point1(2),point1(3),'*r','Linewidth',10)
%         plot3(point2(1),point2(2),point2(3),'*b','Linewidth',10)
%         hold off
%         daspect([1 1 1])
%         axis([min(handles.xd(:)) max(handles.xd(:)) min(handles.yd(:)) max(handles.yd(:)) min(handles.zd(:)) max(handles.zd(:))])
%         grid on 
%         axis vis3d
%         xlabel('ejex')
%         ylabel('ejey')
%         zlabel('ejez')

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

        handles.id_section = r;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        id_selected_element = zeros(size(Selected_element(:)));
        for n=1:length(r)
            [rn,~,~] = find(Selected_element(:)==r(n)); 
            id_selected_element(rn) = 1;
        end
        id_selected_element = reshape(id_selected_element,[size(Selected_element,1),size(Selected_element,2)]);
        [rn,~,~] = find(sum(id_selected_element,2)==4);
        handles.id_selected_faceid = rn;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % node volume
        nodevol = nodevolume(handles.nodes,handles.elem);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % positive direction of the velocity verification.
        v1 = negative_vector/norm(negative_vector);
        v2 = VR_n/norm(VR_n);

        angle_dir = acos(sum(v1.*v2))*180/pi;
        if angle_dir<90
            VR_n = VR_n*-1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proyection of the velocity
        normal = VR_n/norm(VR_n);
        handles.velocity_proj = sum(handles.veset(r,:,:).*repmat(normal,size(handles.veset(r,:,:),1),1,size(handles.veset(r,:,:),3)),2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flow quantification
        nod_volume_selected = nodevol(r)/sum(nodevol(r));
        handles.flow = squeeze(sum((handles.velocity_proj*100).*nod_volume_selected*(area/100)));
        handles.time = linspace(0,60/handles.heart_rate,size(handles.veset,3))';
        [~,I] = max(abs(handles.flow));
        handles.peak_flow = I;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % maximum and minimum velocity
        handles.max_velocity = squeeze(max(handles.velocity_proj*100));
        handles.min_velocity = squeeze(min(handles.velocity_proj*100));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % netflow
        handles.net_flow = cumsum(handles.flow.*handles.time);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % show results in the popupmenu3
        handles.list_string = {'Select Results ...','Mesh + Slice','Flow','Net Flow','Maximum velocity','Minimum velocity'};
        set(handles.popupmenu3,'Visible','on','String',handles.list_string,'Value',2);

        pond = squeeze(max(abs(handles.velocity_proj)))/max(squeeze(max(abs(handles.velocity_proj))))*5;

        axes(handles.axes1);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
        hold on
        patch('faces',handles.faceid(handles.id_selected_faceid,:),'vertices',handles.cutpos,'facecolor','r','edgecolor','r')
        q = quiver3(handles.nodes(handles.id_section,1),handles.nodes(handles.id_section,2),handles.nodes(handles.id_section,3),...
                handles.veset(handles.id_section,1,handles.peak_flow),handles.veset(handles.id_section,2,handles.peak_flow),handles.veset(handles.id_section,3,handles.peak_flow),pond(handles.peak_flow),'LineWidth',1);
        hold off
        currentColormap = colormap(handles.axes1,'jet');
        [~, ~, ind] = histcounts(handles.velocity_proj(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.velocity_proj,1) size(handles.velocity_proj,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes1);
        handles.min_vel = min(handles.velocity_proj(:));
        handles.max_vel = max(handles.velocity_proj(:));
        handles.mean_vel = (handles.min_vel + handles.max_vel)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_vel handles.max_vel];
        c.Ticks = [handles.min_vel, (handles.min_vel + handles.mean_vel)/2, handles.mean_vel, (handles.max_vel + handles.mean_vel)/2, handles.max_vel];
        c.TickLabels = {num2str(handles.min_vel,'%0.2f'), num2str((handles.min_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.mean_vel,'%0.2f'), num2str((handles.max_vel + handles.mean_vel)/2,'%0.2f'), num2str(handles.max_vel,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.3];
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
        c.Label.Color = [1 1 1];
        caxis(handles.axes1, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])

        handles.id_view = 0;
        handles.id_image = 0;
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.slider1,'visible','on')
        set(handles.text1,'visible','on')
        set(handles.pushbutton1,'visible','off')

        % slider adjustment
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text1,'String',['# Cardiac Phase: ',num2str(handles.peak_flow)])

        handles.id_mesh_inlet_flow = 1;
        handles.id_result = 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % waitbar percentaje
        waitbar(98/100)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        close(h)
    end
% end
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
% Hint: delete(hObject) closes the figure


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(handles.popupmenu3,'Value');

if val == 3
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');

        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.flow,'-r','LineWidth',2)
        hold on 
        plot([handles.time(handles.peak_flow),handles.time(handles.peak_flow)],[min(handles.flow) max(handles.flow)],'-k','LineWidth',2)
        hold off
        text_1 = '\leftarrow Peak flow';
        text_2 = ['\leftarrow CP#', num2str(handles.peak_flow)];
        txt = {text_1,text_2};
        text(handles.time(handles.peak_flow),mean([min(handles.flow) max(handles.flow)]),txt,'FontWeight','bold','Color','b')
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.flow) max(handles.flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Flow [cm^{3}/s]';
        ax.XLabel.String = 'Time [s]';

elseif val == 4
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.net_flow,'-b','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.net_flow) max(handles.net_flow)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Net Flow [ml]';
        ax.XLabel.String = 'Time [s]';

elseif val == 5
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.max_velocity,'-g','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.max_velocity) max(handles.max_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Maximum Velocity [cm/s]';
        ax.XLabel.String = 'Time [s]';
        
elseif val == 6
        
        set(handles.popupmenu1,'Value',1);
        set(handles.popupmenu2,'Value',1);
        set(handles.pushbutton1,'Visible','off');
        set(handles.text1,'Visible','off');
        set(handles.slider1,'Visible','off');
        
        set(handles.pushbutton1, 'FontUnits', 'points');
        FontS = get(handles.pushbutton1, 'FontSize');
        set(handles.pushbutton1, 'FontUnits', 'normalized','FontSize',0.26);
        
        axes(handles.axes1);
        plot(0.0)
        plot(handles.time,handles.min_velocity,'-k','LineWidth',2)
        axis square
        ax = gca;
        ax.XColor = 'w';
        ax.YColor = 'w';
        ax.XLim = [handles.time(1) handles.time(end)];
        ax.YLim = [min(handles.min_velocity) max(handles.min_velocity)];
        ax.FontWeight = 'bold';
        ax.FontSize = FontS;
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        ax.GridColor = 'k';
        ax.YLabel.String = 'Minimum Velocity [cm/s]';
        ax.XLabel.String = 'Time [s]';

end

    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');

    set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.37)
    set(handles.popupmenu2,'FontUnits','Normalized','FontSize',0.37)
    set(handles.popupmenu3,'FontUnits','Normalized','FontSize',0.37)
    set(handles.text1,'FontUnits','Normalized','FontSize',0.78)
    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.26)

handles.output = hObject;
guidata(hObject, handles);  


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filter = {'*.xls'};
[file,path] = uiputfile(filter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,['saving flow values...']);

mean_velocity = squeeze(mean(handles.velocity_proj*100));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VV = {'Time [s]';...
      'Flow [ml/s]';...
      'Net Flow [ml]';...
      'Mean Velocity [cm/s]';
      'Maximum Velocity [cm/s]';...
      'Minimum Velocity [cm/s]'};
MM = [handles.time';handles.flow';handles.net_flow';mean_velocity';handles.max_velocity';handles.min_velocity'];
TT = table(VV,MM);
TT.Properties.VariableNames = {'Parameter','Cardiac_Phase'};

filename = [path,file];
writetable(TT,fullfile(filename),'Sheet','2D Flow','Range','A1')

close(h)

handles.output = hObject;
guidata(hObject, handles);
