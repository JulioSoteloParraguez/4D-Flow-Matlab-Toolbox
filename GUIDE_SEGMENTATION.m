function varargout = GUIDE_SEGMENTATION(varargin)
% GUIDE_SEGMENTATION MATLAB code for GUIDE_SEGMENTATION.fig
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
                   'gui_OpeningFcn', @GUIDE_SEGMENTATION_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_SEGMENTATION_OutputFcn, ...
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
function GUIDE_SEGMENTATION_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    set(handles.pushbutton6,'String','RESET');
    set(handles.text1,'String','# Slice: 0.0');
    set(handles.pushbutton1, 'Units', 'pixels');
    handles.pushbutton1_size = get(handles.pushbutton1, 'Position');
    set(handles.pushbutton1, 'Units', 'normalized');
    handles.id = varargin{1}.id;
    
    if handles.id ==1
        idx = mod(min(handles.pushbutton1_size(3:4)),2)>1;
        w = floor(min(handles.pushbutton1_size(3:4)));
        w(idx) = w(idx)+1;
        IM = imread('Symbols/P_Cut.tiff');
        IM_Cut_rgb1 = double(IM(:,:,1:3)<=1);
        IM_Cut_rgb1(:,:,2) = IM_Cut_rgb1(:,:,2)*0;
        IM_Cut_rgb1(:,:,3) = IM_Cut_rgb1(:,:,3)*0;
        f1=double(imresize(IM_Cut_rgb1,[w-20 w-20],'method','nearest')>0);
        IM = imread('Symbols/P_Magic.tiff');
        IM_Magic_rgb1 = double(IM(:,:,1:3)<=1);
        IM_Magic_rgb1(:,:,2) = IM_Magic_rgb1(:,:,2)*0;
        IM_Magic_rgb1(:,:,3) = IM_Magic_rgb1(:,:,3)*0;
        f2=double(imresize(IM_Magic_rgb1,[w-20 w-20],'method','nearest')>0);
        IM = imread('Symbols/P_Ellipse.tiff');
        IM_Ellipse_rgb1 = double(IM(:,:,1:3)<=1);
        IM_Ellipse_rgb1(:,:,2) = IM_Ellipse_rgb1(:,:,2)*0;
        IM_Ellipse_rgb1(:,:,3) = IM_Ellipse_rgb1(:,:,3)*0;
        f3=double(imresize(IM_Ellipse_rgb1,[w-20 w-20],'method','nearest')>0);
        IM = imread('Symbols/P_Cut.tiff');
        IM_Cut_rgb2 = double(IM(:,:,1:3)<=1);
        IM_Cut_rgb2(:,:,1) = IM_Cut_rgb2(:,:,1)*0;
        IM_Cut_rgb2(:,:,3) = IM_Cut_rgb2(:,:,3)*0;
        f4=double(imresize(IM_Cut_rgb2,[w-20 w-20],'method','nearest')>0);
        IM = imread('Symbols/P_Fill.tiff');
        IM_Ellipse_rgb2 = double(IM(:,:,1:3)<=1);
        IM_Ellipse_rgb2(:,:,1) = IM_Ellipse_rgb2(:,:,1)*0;
        IM_Ellipse_rgb2(:,:,3) = IM_Ellipse_rgb2(:,:,3)*0;
        f5=double(imresize(IM_Ellipse_rgb2,[w-20 w-20],'method','nearest')>0);
        set(handles.pushbutton1,'CData',f1,'visible','on');
        set(handles.pushbutton2,'CData',f2,'visible','on');
        set(handles.pushbutton3,'CData',f3,'visible','on');
        set(handles.pushbutton4,'CData',f4,'visible','on');
        set(handles.pushbutton5,'CData',f5,'visible','on');
        
        handles.IPCMRA          = varargin{1}.IPCMRA;
        handles.ANG             = varargin{1}.ANG; % Julio Sotelo
        handles.MAG             = varargin{1}.MAG; % Julio Sotelo
        handles.id_ang          = varargin{1}.id_ang;% Julio Sotelo
        handles.id_mag          = varargin{1}.id_mag;% Julio Sotelo
        handles.idangmag        = varargin{1}.idangmag; % Julio Sotelo
        handles.SEG             = varargin{1}.SEG;
        handles.L               = varargin{1}.L;
        handles.NUM             = varargin{1}.NUM;
        handles.Lrgb            = varargin{1}.Lrgb;
        handles.min_value_th    = varargin{1}.min_value_th;
        handles.max_value_th    = varargin{1}.max_value_th;
        handles.xd              = varargin{1}.xd;
        handles.yd              = varargin{1}.yd;
        handles.zd              = varargin{1}.zd;
        handles.a               = varargin{1}.a;
        handles.b               = varargin{1}.b;
        handles.c               = varargin{1}.c;
        handles.slider_axes1    = varargin{1}.slider_axes1;
        handles.slider_axes2    = varargin{1}.slider_axes2;
        handles.slider_axes3    = varargin{1}.slider_axes3;
        handles.id_seg          = 1;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.view_sac        = varargin{1}.view_sac;
        handles.SEG_old         = varargin{1}.SEG_old;
        handles.Lrgb_old        = varargin{1}.Lrgb_old;
        handles.L_old           = varargin{1}.L_old;
        handles.NUM_old         = varargin{1}.NUM_old;
        handles.id_seg        = varargin{1}.id_seg;
        handles.id_vel        = varargin{1}.id_vel;
        handles.id_vor        = varargin{1}.id_vor;
        handles.id_hd         = varargin{1}.id_hd;
        handles.id_rhd        = varargin{1}.id_rhd;
        handles.id_vd         = varargin{1}.id_vd;
        handles.id_el         = varargin{1}.id_el;
        handles.id_ke         = varargin{1}.id_ke;
        handles.Lrgb_vel      = varargin{1}.Lrgb_vel;
        handles.Lrgb_vor      = varargin{1}.Lrgb_vor;
        handles.Lrgb_hd       = varargin{1}.Lrgb_hd;
        handles.Lrgb_rhd      = varargin{1}.Lrgb_rhd;
        handles.Lrgb_vd       = varargin{1}.Lrgb_vd;
        handles.Lrgb_el       = varargin{1}.Lrgb_el;
        handles.Lrgb_ke       = varargin{1}.Lrgb_ke;
        handles.peak_flow      = varargin{1}.peak_flow;
        
        handles.id_fve        = varargin{1}.id_fve;
        handles.Lrgb_fve      = varargin{1}.Lrgb_fve;
        handles.id_bve        = varargin{1}.id_bve;
        handles.Lrgb_bve      = varargin{1}.Lrgb_bve;
        handles.id_aan        = varargin{1}.id_aan;
        handles.Lrgb_aan      = varargin{1}.Lrgb_aan;
        handles.id_fov        = varargin{1}.id_fov;
        handles.Lrgb_fov      = varargin{1}.Lrgb_fov;
            
        if handles.idangmag == 1% Julio Sotelo
            set(handles.popupmenu2,'Value',1)% Julio Sotelo
            handles.IPCMRA = handles.ANG;
        elseif handles.idangmag == 2% Julio Sotelo
            set(handles.popupmenu2,'Value',2)% Julio Sotelo
            handles.IPCMRA = handles.MAG;
        end% Julio Sotelo
        

        handles.id_while = 0;
        if handles.view_sac == 1
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 1;
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])
            
        elseif handles.view_sac == 2
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 2;
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes2)])
            
        elseif handles.view_sac == 3
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 3;
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes3)])
        end

        set(handles.popupmenu1,'Value',handles.view_sac)

    elseif handles.id ==2
        handles.idangmag        = varargin{1}.idangmag; % Julio Sotelo
        handles.IPCMRA          = varargin{1}.IPCMRA;
        handles.SEG             = varargin{1}.SEG;
        handles.L               = varargin{1}.L;
        handles.NUM             = varargin{1}.NUM;
        handles.Lrgb            = varargin{1}.Lrgb;
        handles.xd              = varargin{1}.xd;
        handles.yd              = varargin{1}.yd;
        handles.zd              = varargin{1}.zd;
        handles.a               = varargin{1}.a;
        handles.b               = varargin{1}.b;
        handles.c               = varargin{1}.c;
        handles.slider_axes1    = varargin{1}.slider_axes1;
        handles.slider_axes2    = varargin{1}.slider_axes2;
        handles.slider_axes3    = varargin{1}.slider_axes3;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.view_sac        = varargin{1}.view_sac;
        handles.id_seg        = varargin{1}.id_seg;
        handles.id_vel        = varargin{1}.id_vel;
        handles.id_vor        = varargin{1}.id_vor;
        handles.id_hd         = varargin{1}.id_hd;
        handles.id_rhd        = varargin{1}.id_rhd;
        handles.id_vd         = varargin{1}.id_vd;
        handles.id_el         = varargin{1}.id_el;
        handles.id_ke         = varargin{1}.id_ke;
        handles.Lrgb_vel      = varargin{1}.Lrgb_vel;
        handles.Lrgb_vor      = varargin{1}.Lrgb_vor;
        handles.Lrgb_hd       = varargin{1}.Lrgb_hd;
        handles.Lrgb_rhd      = varargin{1}.Lrgb_rhd;
        handles.Lrgb_vd       = varargin{1}.Lrgb_vd;
        handles.Lrgb_el       = varargin{1}.Lrgb_el;
        handles.Lrgb_ke       = varargin{1}.Lrgb_ke;
        
        handles.id_fve        = varargin{1}.id_fve;
        handles.Lrgb_fve      = varargin{1}.Lrgb_fve;
        handles.id_bve        = varargin{1}.id_bve;
        handles.Lrgb_bve      = varargin{1}.Lrgb_bve;
        handles.id_aan        = varargin{1}.id_aan;
        handles.Lrgb_aan      = varargin{1}.Lrgb_aan;
        handles.id_fov        = varargin{1}.id_fov;
        handles.Lrgb_fov      = varargin{1}.Lrgb_fov;
        
        
        handles.peak_flow      = varargin{1}.peak_flow;
        set(handles.popupmenu1,'visible','off');
        set(handles.popupmenu2,'visible','off');
        set(handles.pushbutton1,'visible','off');
        set(handles.pushbutton2,'visible','off');
        set(handles.pushbutton3,'visible','off');
        set(handles.pushbutton4,'visible','off');
        set(handles.pushbutton5,'visible','off');
        set(handles.pushbutton6,'visible','off');
        set(handles.pushbutton7,'visible','off');
        set(handles.pushbutton8,'visible','off');
        if handles.view_sac == 1
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 1;
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])
        elseif handles.view_sac == 2
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 2;
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes2)])
        elseif handles.view_sac == 3
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 3;
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes3)])
        end
    elseif handles.id ==3
        
        handles.idangmag        = varargin{1}.idangmag; % Julio Sotelo
        handles.ref             = varargin{1}.ref;
        handles.IPCMRA          = varargin{1}.IPCMRA;
        handles.voxel_MR        = varargin{1}.voxel_MR;
        handles.L               = varargin{1}.L;
        handles.Lrgb            = varargin{1}.Lrgb;
        handles.SEG             = varargin{1}.SEG;
        handles.faces           = varargin{1}.faces;
        handles.nodes           = varargin{1}.nodes;
        handles.d               = varargin{1}.d;
        handles.peak_flow       = varargin{1}.peak_flow;
        handles.veset           = varargin{1}.veset;
        handles.mags_vel        = varargin{1}.mags_vel;
        handles.WSS             = varargin{1}.WSS;
        handles.VOR             = varargin{1}.VOR;
        handles.mags_vor        = varargin{1}.mags_vor;
        handles.mags_wss        = varargin{1}.mags_wss;
        handles.mags_osi        = varargin{1}.mags_osi;
        handles.mags_hd         = varargin{1}.mags_hd;
        handles.mags_rhd        = varargin{1}.mags_rhd;
        handles.mags_vd         = varargin{1}.mags_vd;
        handles.mags_el         = varargin{1}.mags_el;
        handles.mags_ke         = varargin{1}.mags_ke;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.Laplace                       = varargin{1}.Laplace;                         
        handles.centerline                    = varargin{1}.centerline;                      
        handles.radius                        = varargin{1}.radius;                        
        handles.diameter                      = varargin{1}.diameter;                        
        handles.axial_unit_vectors            = varargin{1}.axial_unit_vectors;             
        handles.circumferential_unit_vectors  = varargin{1}.circumferential_unit_vectors;   
        handles.WSS_A                         = varargin{1}.WSS_A;                          
        handles.WSS_C                         = varargin{1}.WSS_C;                           
        handles.mag_WSS_A                     = varargin{1}.mag_WSS_A;                       
        handles.mag_WSS_C                     = varargin{1}.mag_WSS_C;                       
        handles.angle_axial_direction         = varargin{1}.angle_axial_direction;           
        handles.forward_velocity              = varargin{1}.forward_velocity;                
        handles.backward_velocity             = varargin{1}.backward_velocity;               
        handles.mag_forward_velocity          = varargin{1}.mag_forward_velocity;           
        handles.mag_backward_velocity         = varargin{1}.mag_backward_velocity;           
        handles.regurgitant_flow              = varargin{1}.regurgitant_flow;                
        handles.centerline_flow               = varargin{1}.centerline_flow;                 
        handles.eccentricity                  = varargin{1}.eccentricity;    
        handles.curvature                     = varargin{1}.curvature; % Julio Sotelo 05062019
        handles.ellipticity                   = varargin{1}.ellipticity; % Julio Sotelo 05062019
        handles.length_vessel                 = varargin{1}.length_vessel; % Julio Sotelo 05062019
%         handles.circulation                   = varargin{1}.circulation; % Julio Sotelo 05062019
        handles.forward_vortex                = varargin{1}.forward_vortex; % Julio Sotelo 05062019
        handles.flattening                    = varargin{1}.flattening; % Julio Sotelo 05062019
        handles.area                          = varargin{1}.area; % Julio Sotelo 05062019
        handles.axial_circulation             = varargin{1}.axial_circulation; % Julio Sotelo 05062019
        
        handles.cont_show = 1; % Julio Sotelo
        
        set(handles.popupmenu1,'visible','off');
        set(handles.popupmenu2,'visible','off')
        set(handles.pushbutton1,'visible','off');
        set(handles.pushbutton2,'visible','off');
        set(handles.pushbutton3,'visible','off');
        set(handles.pushbutton4,'visible','off');
        set(handles.pushbutton5,'visible','off');
        set(handles.pushbutton6,'visible','off');
        set(handles.pushbutton7,'visible','off');
        set(handles.pushbutton8,'visible','off');
        set(handles.slider1,'visible','off');
        set(handles.text1,'visible','off');
        set(handles.uitoolbar1,'visible','on'); % Julio Sotelo
        if handles.ref==2
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
            xd = X*handles.voxel_MR(1);
            yd = Y*handles.voxel_MR(2);
            zd = Z*handles.voxel_MR(3);
            xd = permute(xd,[2 1 3]);
            yd = permute(yd,[2 1 3]);
            zd = permute(zd,[2 1 3]);
            id_L = unique(sort(handles.L(:)));
            id_L(id_L==0)=[];
            axes(handles.axes1);
            plot(0,0)
            axis off
            axes(handles.axes1);
            for n=1:length(id_L)
                S_SEG = double(handles.L==id_L(n));
                S = S_SEG(:);
                R = squeeze(handles.Lrgb(:,:,:,1));
                G = squeeze(handles.Lrgb(:,:,:,2));
                B = squeeze(handles.Lrgb(:,:,:,3));
                R = R(:);
                G = G(:);
                B = B(:);
                data = smooth3(S_SEG,'box',3);
                fv = isosurface(xd,yd,zd,data,.5);
                p1 = patch(fv,'FaceColor',[mean(R(S==1)), mean(G(S==1)), mean(B(S==1))],'EdgeColor','k');
                hold on
                axis vis3d
                lighting gouraud
                daspect([1,1,1])
                axis off
                view([-34,-51])
            end
            hold off
        elseif handles.ref==3
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
            xd = X*handles.voxel_MR(1);
            yd = Y*handles.voxel_MR(2);
            zd = Z*handles.voxel_MR(3);
            xd = permute(xd,[2 1 3]);
            yd = permute(yd,[2 1 3]);
            zd = permute(zd,[2 1 3]);
            IMG = handles.SEG;
            [node_bin_aorta,elem_bin_aorta] = binsurface(IMG,4);
            node_bin_aorta(:,1) = node_bin_aorta(:,1)*handles.voxel_MR(1)-0.5*handles.voxel_MR(1);
            node_bin_aorta(:,2) = node_bin_aorta(:,2)*handles.voxel_MR(2)-0.5*handles.voxel_MR(2);
            node_bin_aorta(:,3) = node_bin_aorta(:,3)*handles.voxel_MR(3)-0.5*handles.voxel_MR(3);
            axes(handles.axes1);
            plot(0,0)
            axis off
            axes(handles.axes1);
            patch('Vertices', node_bin_aorta, 'Faces', elem_bin_aorta,'FaceColor',[0.85 0.85 0.85],'EdgeColor','k');
            lighting gouraud
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        elseif handles.ref==4
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
            xd = X*handles.voxel_MR(1);
            yd = Y*handles.voxel_MR(2);
            zd = Z*handles.voxel_MR(3);
            xd = permute(xd,[2 1 3]);
            yd = permute(yd,[2 1 3]);
            zd = permute(zd,[2 1 3]);
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])
        elseif handles.ref==5
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==6
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            c.LimitsMode = 'manual';
            handles.min_wss = min(handles.mags_wss(:));
            handles.max_wss = max(handles.mags_wss(:));
            handles.mean_wss = (handles.min_wss + handles.max_wss)/2;
            c.Limits = [handles.min_wss handles.max_wss];
            c.Ticks = [handles.min_wss, (handles.min_wss + handles.mean_wss)/2, handles.mean_wss, (handles.max_wss + handles.mean_wss)/2, handles.max_wss];
            c.TickLabels = {num2str(handles.min_wss,'%0.2f'), num2str((handles.min_wss + handles.mean_wss)/2,'%0.2f'), num2str(handles.mean_wss,'%0.2f'), num2str((handles.mean_wss + handles.max_wss)/2,'%0.2f'), num2str(handles.max_wss,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wss handles.max_wss]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==7
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_osi,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            c.LimitsMode = 'manual';
            handles.min_osi = min(handles.mags_osi(:));
            handles.max_osi = max(handles.mags_osi(:));
            handles.mean_osi = (handles.min_osi + handles.max_osi)/2;
            c.Limits = [handles.min_osi handles.max_osi];
            c.Ticks = [handles.min_osi, (handles.min_osi + handles.mean_osi)/2, handles.mean_osi, (handles.max_osi + handles.mean_osi)/2, handles.max_osi];
            c.TickLabels = {num2str(handles.min_osi,'%0.2f'), num2str((handles.min_osi + handles.mean_osi)/2,'%0.2f'), num2str(handles.mean_osi,'%0.2f'), num2str((handles.mean_osi + handles.max_osi)/2,'%0.2f'), num2str(handles.max_osi,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'OSI [-]';
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
            caxis(handles.axes1, [handles.min_osi handles.max_osi]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        elseif handles.ref==8
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'cool');
            [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vor = min(handles.mags_vor(:));
            handles.max_vor = max(handles.mags_vor(:));
            handles.mean_vor = (handles.min_vor + handles.max_vor)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vor handles.max_vor];
            c.Ticks = [handles.min_vor, (handles.min_vor + handles.mean_vor)/2, handles.mean_vor, (handles.max_vor + handles.mean_vor)/2, handles.max_vor];
            c.TickLabels = {num2str(handles.min_vor,'%0.2f'), num2str((handles.min_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.mean_vor,'%0.2f'), num2str((handles.max_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.max_vor,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Vorticity [1/s]';
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
            caxis(handles.axes1, [handles.min_vor handles.max_vor]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==9
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_hd = min(handles.mags_hd(:));
            handles.max_hd = max(handles.mags_hd(:));
            handles.mean_hd = (handles.min_hd + handles.max_hd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_hd handles.max_hd];
            c.Ticks = [handles.min_hd, (handles.min_hd + handles.mean_hd)/2, handles.mean_hd, (handles.max_hd + handles.mean_hd)/2, handles.max_hd];
            c.TickLabels = {num2str(handles.min_hd,'%0.2f'), num2str((handles.min_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.mean_hd,'%0.2f'), num2str((handles.max_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.max_hd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Helicity Density [m/s^{2}]';
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
            caxis(handles.axes1, [handles.min_hd handles.max_hd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==10
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'summer');
            c = colorbar(handles.axes1);
            handles.min_rhd = min(handles.mags_rhd(:));
            handles.max_rhd = max(handles.mags_rhd(:));
            handles.mean_rhd = (handles.min_rhd + handles.max_rhd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_rhd handles.max_rhd];
            c.Ticks = [handles.min_rhd, (handles.min_rhd + handles.mean_rhd)/2, handles.mean_rhd, (handles.max_rhd + handles.mean_rhd)/2, handles.max_rhd];
            c.TickLabels = {num2str(handles.min_rhd,'%0.2f'), num2str((handles.min_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.mean_rhd,'%0.2f'), num2str((handles.max_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.max_rhd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Relative Helicity Density [-]';
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
            caxis(handles.axes1, [handles.min_rhd handles.max_rhd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==11
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'winter');
            c = colorbar(handles.axes1);
            handles.min_vd = min(handles.mags_vd(:));
            handles.max_vd = max(handles.mags_vd(:));
            handles.mean_vd = (handles.min_vd + handles.max_vd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vd handles.max_vd];
            c.Ticks = [handles.min_vd, (handles.min_vd + handles.mean_vd)/2, handles.mean_vd, (handles.max_vd + handles.mean_vd)/2, handles.max_vd];
            c.TickLabels = {num2str(handles.min_vd,'%0.2f'), num2str((handles.min_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.mean_vd,'%0.2f'), num2str((handles.max_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.max_vd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Viscous Dissipation [1e^{3}/s^{2}]'; 
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
            caxis(handles.axes1, [handles.min_vd handles.max_vd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==12
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'autumn');
            c = colorbar(handles.axes1);
            handles.min_el = min(handles.mags_el(:));
            handles.max_el = max(handles.mags_el(:));
            handles.mean_el = (handles.min_el + handles.max_el)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_el handles.max_el];
            c.Ticks = [handles.min_el, (handles.min_el + handles.mean_el)/2, handles.mean_el, (handles.max_el + handles.mean_el)/2, handles.max_el];
            c.TickLabels = {num2str(handles.min_el,'%0.2f'), num2str((handles.min_el + handles.mean_el)/2,'%0.2f'), num2str(handles.mean_el,'%0.2f'), num2str((handles.max_el + handles.mean_el)/2,'%0.2f'), num2str(handles.max_el,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Energy Loss [\muW]';
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
            caxis(handles.axes1, [handles.min_el handles.max_el]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref==13
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'spring');
            c = colorbar(handles.axes1);
            handles.min_ke = min(handles.mags_ke(:));
            handles.max_ke = max(handles.mags_ke(:));
            handles.mean_ke = (handles.min_ke + handles.max_ke)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_ke handles.max_ke];
            c.Ticks = [handles.min_ke, (handles.min_ke + handles.mean_ke)/2, handles.mean_ke, (handles.max_ke + handles.mean_ke)/2, handles.max_ke];
            c.TickLabels = {num2str(handles.min_ke,'%0.2f'), num2str((handles.min_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.mean_ke,'%0.2f'), num2str((handles.max_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.max_ke,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Kinetic Energy [\muJ]';
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
            caxis(handles.axes1, [handles.min_ke handles.max_ke]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        
        elseif handles.ref==14 % Laplace
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            c.LimitsMode = 'manual';
            handles.min_lap = min(handles.Laplace(:));
            handles.max_lap = max(handles.Laplace(:));
            handles.mean_lap = (handles.min_lap + handles.max_lap)/2;
            c.Limits = [handles.min_lap handles.max_lap];
            c.Ticks = [handles.min_lap, (handles.min_lap + handles.mean_lap)/2, handles.mean_lap, (handles.max_lap + handles.mean_lap)/2, handles.max_lap];
            c.TickLabels = {num2str(handles.min_lap,'%0.2f'), num2str((handles.min_lap + handles.mean_lap)/2,'%0.2f'), num2str(handles.mean_lap,'%0.2f'), num2str((handles.mean_lap + handles.max_lap)/2,'%0.2f'), num2str(handles.max_lap,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_lap handles.max_lap]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        
        elseif handles.ref==15 % centerline
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            plot3(handles.centerline(:,1),handles.centerline(:,2),handles.centerline(:,3),'-c','LineWidth',3)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])
 
        elseif handles.ref==16 % diameter
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.diameter,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            c.LimitsMode = 'manual';
            handles.min_dia = min(handles.diameter(:));
            handles.max_dia = max(handles.diameter(:));
            handles.mean_dia = (handles.min_dia + handles.max_dia)/2;
            c.Limits = [handles.min_dia handles.max_dia];
            c.Ticks = [handles.min_dia, (handles.min_dia + handles.mean_dia)/2, handles.mean_dia, (handles.max_dia + handles.mean_dia)/2, handles.max_dia];
            c.TickLabels = {num2str(handles.min_dia,'%0.2f'), num2str((handles.min_dia + handles.mean_dia)/2,'%0.2f'), num2str(handles.mean_dia,'%0.2f'), num2str((handles.max_dia + handles.mean_dia)/2,'%0.2f'), num2str(handles.max_dia,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_dia handles.max_dia]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        elseif handles.ref==17 % radius
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.radius,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            handles.min_rad = min(handles.radius(:));
            handles.max_rad = max(handles.radius(:));
            handles.mean_rad = (handles.min_rad + handles.max_rad)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_rad handles.max_rad];
            c.Ticks = [handles.min_rad, (handles.min_rad + handles.mean_rad)/2, handles.mean_rad, (handles.max_rad + handles.mean_rad)/2, handles.max_rad];
            c.TickLabels = {num2str(handles.min_rad,'%0.2f'), num2str((handles.min_rad + handles.mean_rad)/2,'%0.2f'), num2str(handles.mean_rad,'%0.2f'), num2str((handles.max_rad + handles.mean_rad)/2,'%0.2f'), num2str(handles.max_rad,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_rad handles.max_rad]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        elseif handles.ref==18 % axial unit vector 
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.axial_unit_vectors(:,1),handles.axial_unit_vectors(:,2),handles.axial_unit_vectors(:,3),1,'g','LineWidth',1)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])
 
        elseif handles.ref==19 % circumferential unit vector  
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.circumferential_unit_vectors(:,1),handles.circumferential_unit_vectors(:,2),handles.circumferential_unit_vectors(:,3),1,'r','LineWidth',1)
            hold off
            axis vis3d
            lighting gouraud
            daspect([1,1,1])
            axis off
            view([-34,-51])
 
        elseif handles.ref==20 % WSSA   
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS-Axial [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
            
        elseif handles.ref==21 % WSSC
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS-Circum. [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
 
        elseif handles.ref==22 % Axial Angle
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_axa = min(handles.angle_axial_direction(:));
            handles.max_axa = max(handles.angle_axial_direction(:));
            handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_axa handles.max_axa];
            c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
            c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_axa handles.max_axa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
 
        elseif handles.ref==23 % Forward Velocity
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
 
        elseif handles.ref==24 % Backward Velocity    
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
 
        elseif handles.ref==25 % Regurgitant Flow
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.regurgitant_flow,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'summer');
            c = colorbar(handles.axes1);
            handles.min_rfl = min(handles.regurgitant_flow(:));
            handles.max_rfl = max(handles.regurgitant_flow(:));
            handles.mean_rfl = (handles.min_rfl + handles.max_rfl)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_rfl handles.max_rfl];
            c.Ticks = [handles.min_rfl, (handles.min_rfl + handles.mean_rfl)/2, handles.mean_rfl, (handles.max_rfl + handles.mean_rfl)/2, handles.max_rfl];
            c.TickLabels = {num2str(handles.min_rfl,'%0.2f'), num2str((handles.min_rfl + handles.mean_rfl)/2,'%0.2f'), num2str(handles.mean_rfl,'%0.2f'), num2str((handles.max_rfl + handles.mean_rfl)/2,'%0.2f'), num2str(handles.max_rfl,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_rfl handles.max_rfl]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
 
        elseif handles.ref==26 % Centerline Flow   
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
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
 
        elseif handles.ref==27 % Eccentricity
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.eccentricity,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'winter');
            c = colorbar(handles.axes1);
            handles.min_ecc = min(handles.eccentricity(:));
            handles.max_ecc = max(handles.eccentricity(:));
            handles.mean_ecc = (handles.min_ecc + handles.max_ecc)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_ecc handles.max_ecc];
            c.Ticks = [handles.min_ecc, (handles.min_ecc + handles.mean_ecc)/2, handles.mean_ecc, (handles.max_ecc + handles.mean_ecc)/2, handles.max_ecc];
            c.TickLabels = {num2str(handles.min_ecc,'%0.2f'), num2str((handles.min_ecc + handles.mean_ecc)/2,'%0.2f'), num2str(handles.mean_ecc,'%0.2f'), num2str((handles.max_ecc + handles.mean_ecc)/2,'%0.2f'), num2str(handles.max_ecc,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_ecc handles.max_ecc]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        
        elseif handles.ref==28 % curvature
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            if min(handles.curvature(:)) == 0 && max(handles.curvature(:)) ==0 % julio
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none')
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
            else    
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.curvature,'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'hot');
                c = colorbar(handles.axes1);
                handles.min_cur = min(handles.curvature(:));
                handles.max_cur = max(handles.curvature(:));
                handles.mean_cur = (handles.min_cur + handles.max_cur)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_cur handles.max_cur];
                c.Ticks = [handles.min_cur, (handles.min_cur + handles.mean_cur)/2, handles.mean_cur, (handles.max_cur + handles.mean_cur)/2, handles.max_cur];
                c.TickLabels = {num2str(handles.min_cur,'%0.2f'), num2str((handles.min_cur + handles.mean_cur)/2,'%0.2f'), num2str(handles.mean_cur,'%0.2f'), num2str((handles.max_cur + handles.mean_cur)/2,'%0.2f'), num2str(handles.max_cur,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_cur handles.max_cur]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
            end
            
        elseif handles.ref==29 % ellipticity
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.ellipticity,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            handles.min_ell = min(handles.ellipticity(:));
            handles.max_ell = max(handles.ellipticity(:));
            handles.mean_ell = (handles.min_ell + handles.max_ell)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_ell handles.max_ell];
            c.Ticks = [handles.min_ell, (handles.min_ell + handles.mean_ell)/2, handles.mean_ell, (handles.max_ell + handles.mean_ell)/2, handles.max_ell];
            c.TickLabels = {num2str(handles.min_ell,'%0.2f'), num2str((handles.min_ell + handles.mean_ell)/2,'%0.2f'), num2str(handles.mean_ell,'%0.2f'), num2str((handles.max_ell + handles.mean_ell)/2,'%0.2f'), num2str(handles.max_ell,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_ell handles.max_ell]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        elseif handles.ref==30 % length_vessel
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.length_vessel,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            handles.min_len = min(handles.length_vessel(:));
            handles.max_len = max(handles.length_vessel(:));
            handles.mean_len = (handles.min_len + handles.max_len)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_len handles.max_len];
            c.Ticks = [handles.min_len, (handles.min_len + handles.mean_len)/2, handles.mean_len, (handles.max_len + handles.mean_len)/2, handles.max_len];
            c.TickLabels = {num2str(handles.min_len,'%0.2f'), num2str((handles.min_len + handles.mean_len)/2,'%0.2f'), num2str(handles.mean_len,'%0.2f'), num2str((handles.max_len + handles.mean_len)/2,'%0.2f'), num2str(handles.max_len,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_len handles.max_len]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
%         elseif handles.ref==31 % circulation
%             
%             axes(handles.axes1);
%             plot(0.0)
%             axis off
%             axes(handles.axes1);
%             patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,handles.peak_flow),'CDataMapping','Scaled')
%             hold on
%             colormap(handles.axes1,'winter');
%             c = colorbar(handles.axes1);
%             handles.min_cir = min(handles.circulation(:));
%             handles.max_cir = max(handles.circulation(:));
%             handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%             c.LimitsMode = 'manual';
%             c.Limits = [handles.min_cir handles.max_cir];
%             c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%             c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%             c.Color = [1 1 1];
%             c.Location = 'manual';
%             c.Position = [0.2 0.1 0.02 0.5];
%             c.FontWeight = 'bold';
%             c.Label.String = 'Circulation [mm^{2}/s]';
%             windows_screen_size = get(0,'ScreenSize');
%             if windows_screen_size(4)<=900
%                 c.FontSize = 9;
%             elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                 c.FontSize = 10;
%             else
%                 c.FontSize = 11;
%             end
%             if windows_screen_size(4)<=900
%                 c.Label.FontSize = 10;
%             elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                 c.Label.FontSize = 12;
%             else
%                 c.Label.FontSize = 14;
%             end
%             c.Label.FontWeight = 'bold';
%             c.Label.Color = [1 1 1];
%             caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%             hold off
%             axis vis3d
%             daspect([1,1,1])
%             axis off
%             view([-34,-51])
%             slider_step(1) = 1/handles.d;
%             slider_step(2) = 0.1;
%             set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
%             set(handles.pushbutton6,'visible','on','String','PLAY');
%             set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
            
        elseif handles.ref==31 % forward_vortex
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_fov = min(handles.forward_vortex(:));
            handles.max_fov = max(handles.forward_vortex(:));
            handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_fov handles.max_fov];
            c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
            c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_fov handles.max_fov]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
            
        elseif handles.ref==32 % flattening
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.flattening,'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'parula');
            c = colorbar(handles.axes1);
            handles.min_fla = min(handles.flattening(:));
            handles.max_fla = max(handles.flattening(:));
            handles.mean_fla = (handles.min_fla + handles.max_fla)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_fla handles.max_fla];
            c.Ticks = [handles.min_fla, (handles.min_fla + handles.mean_fla)/2, handles.mean_fla, (handles.max_fla + handles.mean_fla)/2, handles.max_fla];
            c.TickLabels = {num2str(handles.min_fla,'%0.2f'), num2str((handles.min_fla + handles.mean_fla)/2,'%0.2f'), num2str(handles.mean_fla,'%0.2f'), num2str((handles.max_fla + handles.mean_fla)/2,'%0.2f'), num2str(handles.max_fla,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_fla handles.max_fla]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        elseif handles.ref==33 % area
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.area,'CDataMapping','Scaled')
            hold on
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
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_are handles.max_are]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            
        elseif handles.ref==34 % axial circulation
            
            axes(handles.axes1);
            plot(0.0)
            axis off
            axes(handles.axes1);
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
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
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_acir handles.max_acir]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'visible','on','Value', handles.peak_flow/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.pushbutton6,'visible','on','String','PLAY');
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);


        end
        
        handles.SEG = [];
        handles.L = [];
        handles.Lrgb = [];
        handles.NUM = [];
        handles.id_while = [];
        handles.slider_axes1 = [];
        handles.slider_axes2 = [];
        handles.slider_axes3 = [];
        handles.view_sac = [];
    end
    set(handles.uipanel1, 'Visible', 'on');
guidata(hObject, handles);
uiwait(handles.figure1);
function varargout = GUIDE_SEGMENTATION_OutputFcn(hObject, eventdata, handles)
if handles.id_while == 0
    varargout{1} = handles.output;
    IPCMRA = handles.IPCMRA;
    setappdata(0,'IPCMRA',IPCMRA);
    idangmag = handles.idangmag;
    setappdata(0,'idangmag',idangmag);
    SEG = handles.SEG;
    setappdata(0,'SEG',SEG);
    L = handles.L;
    setappdata(0,'L',L);
    Lrgb = handles.Lrgb;
    setappdata(0,'Lrgb',Lrgb);
    NUM = handles.NUM;
    setappdata(0,'NUM',NUM);
    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
    slider_axes1 = handles.slider_axes1;
    setappdata(0,'slider_axes1',slider_axes1);
    slider_axes2 = handles.slider_axes2;
    setappdata(0,'slider_axes2',slider_axes2);
    slider_axes3 = handles.slider_axes3;
    setappdata(0,'slider_axes3',slider_axes3);
    view_sac = handles.view_sac;
    setappdata(0,'view_sac',view_sac);
    
elseif handles.id_while == 1
    varargout{1} = handles.output;
    IPCMRA = handles.IPCMRA;
    setappdata(0,'IPCMRA',IPCMRA);
    idangmag = handles.idangmag;
    setappdata(0,'idangmag',idangmag);
    SEG = handles.SEG;
    setappdata(0,'SEG',SEG);
    L = handles.L;
    setappdata(0,'L',L);
    Lrgb = handles.Lrgb;
    setappdata(0,'Lrgb',Lrgb);
    NUM = handles.NUM;
    setappdata(0,'NUM',NUM);
    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
    slider_axes1 = handles.slider_axes1;
    setappdata(0,'slider_axes1',slider_axes1);
    slider_axes2 = handles.slider_axes2;
    setappdata(0,'slider_axes2',slider_axes2);
    slider_axes3 = handles.slider_axes3;
    setappdata(0,'slider_axes3',slider_axes3);
    view_sac = handles.view_sac;
    setappdata(0,'view_sac',view_sac);
    delete(handles.figure1);
end
function slider1_Callback(hObject, eventdata, handles)
    if handles.id == 3
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
        handles.peak_flow = handles.slider_value;
        if handles.ref == 5
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            set(handles.text1,'visible','on','String',['Cardiac Phase #:', num2str(handles.peak_flow)]);
        elseif handles.ref == 6
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            c.LimitsMode = 'manual';
            handles.min_wss = min(handles.mags_wss(:));
            handles.max_wss = max(handles.mags_wss(:));
            handles.mean_wss = (handles.min_wss + handles.max_wss)/2;
            c.Limits = [handles.min_wss handles.max_wss];
            c.Ticks = [handles.min_wss, (handles.min_wss + handles.mean_wss)/2, handles.mean_wss, (handles.max_wss + handles.mean_wss)/2, handles.max_wss];
            c.TickLabels = {num2str(handles.min_wss,'%0.2f'), num2str((handles.min_wss + handles.mean_wss)/2,'%0.2f'), num2str(handles.mean_wss,'%0.2f'), num2str((handles.mean_wss + handles.max_wss)/2,'%0.2f'), num2str(handles.max_wss,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wss handles.max_wss]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==8
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'cool');
            [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes1);
            handles.min_vor = min(handles.mags_vor(:));
            handles.max_vor = max(handles.mags_vor(:));
            handles.mean_vor = (handles.min_vor + handles.max_vor)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vor handles.max_vor];
            c.Ticks = [handles.min_vor, (handles.min_vor + handles.mean_vor)/2, handles.mean_vor, (handles.max_vor + handles.mean_vor)/2, handles.max_vor];
            c.TickLabels = {num2str(handles.min_vor,'%0.2f'), num2str((handles.min_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.mean_vor,'%0.2f'), num2str((handles.max_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.max_vor,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Vorticity [1/s]';
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
            caxis(handles.axes1, [handles.min_vor handles.max_vor]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==9
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_hd = min(handles.mags_hd(:));
            handles.max_hd = max(handles.mags_hd(:));
            handles.mean_hd = (handles.min_hd + handles.max_hd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_hd handles.max_hd];
            c.Ticks = [handles.min_hd, (handles.min_hd + handles.mean_hd)/2, handles.mean_hd, (handles.max_hd + handles.mean_hd)/2, handles.max_hd];
            c.TickLabels = {num2str(handles.min_hd,'%0.2f'), num2str((handles.min_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.mean_hd,'%0.2f'), num2str((handles.max_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.max_hd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Helicity Density [m/s^{2}]';
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
            caxis(handles.axes1, [handles.min_hd handles.max_hd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==10
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'summer');
            c = colorbar(handles.axes1);
            handles.min_rhd = min(handles.mags_rhd(:));
            handles.max_rhd = max(handles.mags_rhd(:));
            handles.mean_rhd = (handles.min_rhd + handles.max_rhd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_rhd handles.max_rhd];
            c.Ticks = [handles.min_rhd, (handles.min_rhd + handles.mean_rhd)/2, handles.mean_rhd, (handles.max_rhd + handles.mean_rhd)/2, handles.max_rhd];
            c.TickLabels = {num2str(handles.min_rhd,'%0.2f'), num2str((handles.min_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.mean_rhd,'%0.2f'), num2str((handles.max_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.max_rhd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Relative Helicity Density [-]';
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
            caxis(handles.axes1, [handles.min_rhd handles.max_rhd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==11
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'winter');
            c = colorbar(handles.axes1);
            handles.min_vd = min(handles.mags_vd(:));
            handles.max_vd = max(handles.mags_vd(:));
            handles.mean_vd = (handles.min_vd + handles.max_vd)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_vd handles.max_vd];
            c.Ticks = [handles.min_vd, (handles.min_vd + handles.mean_vd)/2, handles.mean_vd, (handles.max_vd + handles.mean_vd)/2, handles.max_vd];
            c.TickLabels = {num2str(handles.min_vd,'%0.2f'), num2str((handles.min_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.mean_vd,'%0.2f'), num2str((handles.max_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.max_vd,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Viscous Dissipation [1e^{3}/s^{2}]'; 
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
            caxis(handles.axes1, [handles.min_vd handles.max_vd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==12
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'autumn');
            c = colorbar(handles.axes1);
            handles.min_el = min(handles.mags_el(:));
            handles.max_el = max(handles.mags_el(:));
            handles.mean_el = (handles.min_el + handles.max_el)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_el handles.max_el];
            c.Ticks = [handles.min_el, (handles.min_el + handles.mean_el)/2, handles.mean_el, (handles.max_el + handles.mean_el)/2, handles.max_el];
            c.TickLabels = {num2str(handles.min_el,'%0.2f'), num2str((handles.min_el + handles.mean_el)/2,'%0.2f'), num2str(handles.mean_el,'%0.2f'), num2str((handles.max_el + handles.mean_el)/2,'%0.2f'), num2str(handles.max_el,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Energy Loss [\muW]';
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
            caxis(handles.axes1, [handles.min_el handles.max_el]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
        elseif handles.ref==13
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'spring');
            c = colorbar(handles.axes1);
            handles.min_ke = min(handles.mags_ke(:));
            handles.max_ke = max(handles.mags_ke(:));
            handles.mean_ke = (handles.min_ke + handles.max_ke)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_ke handles.max_ke];
            c.Ticks = [handles.min_ke, (handles.min_ke + handles.mean_ke)/2, handles.mean_ke, (handles.max_ke + handles.mean_ke)/2, handles.max_ke];
            c.TickLabels = {num2str(handles.min_ke,'%0.2f'), num2str((handles.min_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.mean_ke,'%0.2f'), num2str((handles.max_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.max_ke,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Kinetic Energy [\muJ]';
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
            caxis(handles.axes1, [handles.min_ke handles.max_ke]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
            
        elseif handles.ref==20 % WSSA   
        
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS-Axial [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
            
        elseif handles.ref==21 % WSSC
        
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
            currentColormap = colormap(handles.axes1,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'WSS-Circum. [N/m^{2}]';
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
            caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])

        elseif handles.ref==22 % Axial Angle
        
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_axa = min(handles.angle_axial_direction(:));
            handles.max_axa = max(handles.angle_axial_direction(:));
            handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_axa handles.max_axa];
            c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
            c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_axa handles.max_axa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])

        elseif handles.ref==23 % Forward Velocity
        
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])

        elseif handles.ref==24 % Backward Velocity    
        
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
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])   
%          elseif handles.ref==31 % circulation
%             
%             axes(handles.axes1);
%             plot(0.0)
%             patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,handles.peak_flow),'CDataMapping','Scaled')
%             hold on
%             colormap(handles.axes1,'winter');
%             c = colorbar(handles.axes1);
%             handles.min_cir = min(handles.circulation(:));
%             handles.max_cir = max(handles.circulation(:));
%             handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%             c.LimitsMode = 'manual';
%             c.Limits = [handles.min_cir handles.max_cir];
%             c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%             c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%             c.Color = [1 1 1];
%             c.Location = 'manual';
%             c.Position = [0.2 0.1 0.02 0.5];
%             c.FontWeight = 'bold';
%             c.Label.String = 'Circulation [mm^{2}/s]';
%             windows_screen_size = get(0,'ScreenSize');
%             if windows_screen_size(4)<=900
%                 c.FontSize = 9;
%             elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                 c.FontSize = 10;
%             else
%                 c.FontSize = 11;
%             end
%             if windows_screen_size(4)<=900
%                 c.Label.FontSize = 10;
%             elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                 c.Label.FontSize = 12;
%             else
%                 c.Label.FontSize = 14;
%             end
%             c.Label.FontWeight = 'bold';
%             c.Label.Color = [1 1 1];
%             caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%             hold off
%             axis vis3d
%             daspect([1,1,1])
%             axis off
%             view([-34,-51])
%             set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])  
            
        elseif handles.ref==31 % forward_vortex
            
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes1,'cool');
            c = colorbar(handles.axes1);
            handles.min_fov = min(handles.forward_vortex(:));
            handles.max_fov = max(handles.forward_vortex(:));
            handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_fov handles.max_fov];
            c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
            c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_fov handles.max_fov]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])  
            
        elseif handles.ref==34 % axial circulation
            
            axes(handles.axes1);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,handles.peak_flow),'CDataMapping','Scaled')
            hold on
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
            c.Position = [0.2 0.1 0.02 0.5];
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
            c.Label.Color = [1 1 1];
            caxis(handles.axes1, [handles.min_acir handles.max_acir]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])  
            
        end
    else
        if handles.view_sac == 1
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
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % Backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle 
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])
        elseif handles.view_sac == 2
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
            handles.slider_axes2 = handles.slider_value;
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % Backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle 
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes2)])
        elseif handles.view_sac == 3
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
            handles.slider_axes3 = handles.slider_value;
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward velocity
            if handles.id_fve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fve(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % Backward velocity
            if handles.id_bve == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_bve(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % axial angle 
            if handles.id_aan == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_aan(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            % forward vortex
            if handles.id_fov == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_fov(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes3)])
        end
    end
handles.output = hObject;
guidata(hObject, handles);
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function pushbutton1_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1
        axes(handles.axes1);
        k = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand(handles.axes1);
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(:,:,handles.slider_axes1,:)=squeeze(handles.Lrgb(:,:,handles.slider_axes1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 2
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand(handles.axes1);
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(:,handles.slider_axes2,:,:)=squeeze(handles.Lrgb(:,handles.slider_axes2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        plot(0.0)
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 3
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand(handles.axes1);
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(handles.slider_axes3,:,:,:)=squeeze(handles.Lrgb(handles.slider_axes3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    end
guidata(hObject, handles);
function pushbutton2_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        [x,y] = getpts(handles.axes1);
        x = x/handles.voxel_MR(1);
        y = y/handles.voxel_MR(2);
        L = bwlabel(squeeze(handles.SEG(:,:,handles.slider_axes1)));
        p = zeros(size(L));
        pn = zeros(size(L));
        for n=1:length(x)
            pn = (L==L(round(y(n)),round(x(n))));
            p = p + pn;
        end
        BW1 = p;
        BW2 = abs(BW1-1);
        handles.Lrgb(:,:,handles.slider_axes1,:)=squeeze(handles.Lrgb(:,:,handles.slider_axes1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 2
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        [x,y] = getpts(handles.axes1);
        x = x/handles.voxel_MR(3);
        y = y/handles.voxel_MR(2);
        L = bwlabel(squeeze(handles.SEG(:,handles.slider_axes2,:)));
        p = zeros(size(L));
        pn = zeros(size(L));
        for n=1:length(x)
            pn = (L==L(round(y(n)),round(x(n))));
            p = p + pn;
        end
        BW1 = p;
        BW2 = abs(BW1-1);
        handles.Lrgb(:,handles.slider_axes2,:,:)=squeeze(handles.Lrgb(:,handles.slider_axes2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 3
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        [x,y] = getpts(handles.axes1);
        x = x/handles.voxel_MR(3);
        y = y/handles.voxel_MR(1);
        L = bwlabel(squeeze(handles.SEG(handles.slider_axes3,:,:)));
        p = zeros(size(L));
        pn = zeros(size(L));
        for n=1:length(x)
            pn = (L==L(round(y(n)),round(x(n))));
            p = p + pn;
        end
        BW1 = p;
        BW2 = abs(BW1-1);
        handles.Lrgb(handles.slider_axes3,:,:,:)=squeeze(handles.Lrgb(handles.slider_axes3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    end
guidata(hObject, handles);
function pushbutton3_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1
        axes(handles.axes1);
        k = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imellipse;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(:,:,handles.slider_axes1,:)=squeeze(handles.Lrgb(:,:,handles.slider_axes1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 2
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imellipse;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(:,handles.slider_axes2,:,:)=squeeze(handles.Lrgb(:,handles.slider_axes2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 3
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imellipse;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        handles.Lrgb(handles.slider_axes3,:,:,:)=squeeze(handles.Lrgb(handles.slider_axes3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    end
guidata(hObject, handles);
function pushbutton4_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1
        axes(handles.axes1);
        k = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        BW3 = imdilate(BW1,[0 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;0 1 1 1 0]);
        H = handles.Lrgb(:,:,handles.slider_axes1-1:handles.slider_axes1+1,:).*repmat(BW3,1,1,3,3);
        MASK = double(sum(H,4)>0 & sum(H,4)<3);
        MASK_Lrgb_R = squeeze(handles.Lrgb(:,:,handles.slider_axes1-1:handles.slider_axes1+1,1));
        MASK_Lrgb_G = squeeze(handles.Lrgb(:,:,handles.slider_axes1-1:handles.slider_axes1+1,2));
        MASK_Lrgb_B = squeeze(handles.Lrgb(:,:,handles.slider_axes1-1:handles.slider_axes1+1,3));
        MASK = MASK(:);
        MASK_Lrgb_R = MASK_Lrgb_R(:);
        MASK_Lrgb_G = MASK_Lrgb_G(:);
        MASK_Lrgb_B = MASK_Lrgb_B(:);
        val = [mean(MASK_Lrgb_R(MASK==1)) mean(MASK_Lrgb_G(MASK==1)) mean(MASK_Lrgb_B(MASK==1))];
        MASK2 = double(repmat(BW1,1,1,3));
        MASK2(:,:,1) = MASK2(:,:,1)*val(1);
        MASK2(:,:,2) = MASK2(:,:,2)*val(2);
        MASK2(:,:,3) = MASK2(:,:,3)*val(3);
        handles.Lrgb(:,:,handles.slider_axes1,:) = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:)).*repmat(BW2,1,1,3) + MASK2;
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 2
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        BW3 = imdilate(BW1,[0 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;0 1 1 1 0]);
        H = handles.Lrgb(:,handles.slider_axes2-1:handles.slider_axes2+1,:,:).*permute(repmat(BW3,1,1,3,3),[1 3 2 4]);
        MASK = double(sum(H,4)>0 & sum(H,4)<3);
        MASK_Lrgb_R = squeeze(handles.Lrgb(:,handles.slider_axes2-1:handles.slider_axes2+1,:,1));
        MASK_Lrgb_G = squeeze(handles.Lrgb(:,handles.slider_axes2-1:handles.slider_axes2+1,:,2));
        MASK_Lrgb_B = squeeze(handles.Lrgb(:,handles.slider_axes2-1:handles.slider_axes2+1,:,3));
        MASK = MASK(:);
        MASK_Lrgb_R = MASK_Lrgb_R(:);
        MASK_Lrgb_G = MASK_Lrgb_G(:);
        MASK_Lrgb_B = MASK_Lrgb_B(:);
        val = [mean(MASK_Lrgb_R(MASK==1)) mean(MASK_Lrgb_G(MASK==1)) mean(MASK_Lrgb_B(MASK==1))];
        MASK = double(repmat(BW1,1,1,3));
        MASK(:,:,1) = MASK(:,:,1)*val(1);
        MASK(:,:,2) = MASK(:,:,2)*val(2);
        MASK(:,:,3) = MASK(:,:,3)*val(3);
        handles.Lrgb(:,handles.slider_axes2,:,:) = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:)).*repmat(BW2,1,1,3) + MASK;
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 3
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = imfreehand;
        BW1 = createMask(p,k);
        BW2 = abs(BW1-1);
        BW3 = imdilate(BW1,[0 1 1 1 0;1 1 1 1 1;1 1 1 1 1;1 1 1 1 1;0 1 1 1 0]);
        H = handles.Lrgb(handles.slider_axes3-1:handles.slider_axes3+1,:,:,:).*permute(repmat(BW3,1,1,3,3),[3,1,2,4]);
        MASK = double(sum(H,4)>0 & sum(H,4)<3);
        MASK_Lrgb_R = squeeze(handles.Lrgb(handles.slider_axes3-1:handles.slider_axes3+1,:,:,1));
        MASK_Lrgb_G = squeeze(handles.Lrgb(handles.slider_axes3-1:handles.slider_axes3+1,:,:,2));
        MASK_Lrgb_B = squeeze(handles.Lrgb(handles.slider_axes3-1:handles.slider_axes3+1,:,:,3));
        MASK = MASK(:);
        MASK_Lrgb_R = MASK_Lrgb_R(:);
        MASK_Lrgb_G = MASK_Lrgb_G(:);
        MASK_Lrgb_B = MASK_Lrgb_B(:);
        val = [mean(MASK_Lrgb_R(MASK==1)) mean(MASK_Lrgb_G(MASK==1)) mean(MASK_Lrgb_B(MASK==1))];
        MASK = double(repmat(BW1,1,1,3));
        MASK(:,:,1) = MASK(:,:,1)*val(1);
        MASK(:,:,2) = MASK(:,:,2)*val(2);
        MASK(:,:,3) = MASK(:,:,3)*val(3);
        handles.Lrgb(handles.slider_axes3,:,:,:) = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:)).*repmat(BW2,1,1,3) + MASK;
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    end
guidata(hObject, handles);

function pushbutton5_Callback(hObject, eventdata, handles)

    if handles.view_sac == 1
        axes(handles.axes1);
        k = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = impoint;
        BW1 = createMask(p,k);
        [r,c,~] = find(BW1==1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        val_Lrgb_R = handles.Lrgb(r(1),c(1),handles.slider_axes1,1);
        val_Lrgb_G = handles.Lrgb(r(1),c(1),handles.slider_axes1,2);
        val_Lrgb_B = handles.Lrgb(r(1),c(1),handles.slider_axes1,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.MM = zeros(size(handles.IPCMRA));% Julio Sotelo
        handles.MM(:,:,handles.slider_axes1) = double(BW1);% Julio Sotelo
        L = double(bwlabeln(handles.SEG,6));% Julio Sotelo
        LM = handles.MM.*double(bwlabeln(handles.SEG,6));
        val = max(LM(:));% Julio Sotelo
        handles.SEG = imfill(double(L==val),8,'holes');% Julio Sotelo
        
        
        MASK1 = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
        MASK1(:,:,:,1) = MASK1(:,:,:,1) - handles.SEG;
        MASK1(:,:,:,2) = MASK1(:,:,:,2) - handles.SEG;
        MASK1(:,:,:,3) = MASK1(:,:,:,3) - handles.SEG;
        
        MASK2 = repmat(handles.SEG,1,1,1,3);
        MASK2(:,:,:,1) = MASK2(:,:,:,1).*val_Lrgb_R;
        MASK2(:,:,:,2) = MASK2(:,:,:,2).*val_Lrgb_G;
        MASK2(:,:,:,3) = MASK2(:,:,:,3).*val_Lrgb_B;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.Lrgb = MASK1 + MASK2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        

    elseif handles.view_sac == 2
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        p = impoint;
        BW1 = createMask(p,k);
        [r,c,~] = find(BW1==1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        val_Lrgb_R = handles.Lrgb(r(1),handles.slider_axes2,c(1),1);
        val_Lrgb_G = handles.Lrgb(r(1),handles.slider_axes2,c(1),2);
        val_Lrgb_B = handles.Lrgb(r(1),handles.slider_axes2,c(1),3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.MM = zeros(size(handles.IPCMRA));% Julio Sotelo
        handles.MM(:,handles.slider_axes2,:) = double(BW1);% Julio Sotelo
        L = double(bwlabeln(handles.SEG,6));% Julio Sotelo
        LM = handles.MM.*double(bwlabeln(handles.SEG,6));
        val = max(LM(:));% Julio Sotelo
        handles.SEG = imfill(double(L==val),8,'holes');% Julio Sotelo
        
        
        MASK1 = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
        MASK1(:,:,:,1) = MASK1(:,:,:,1) - handles.SEG;
        MASK1(:,:,:,2) = MASK1(:,:,:,2) - handles.SEG;
        MASK1(:,:,:,3) = MASK1(:,:,:,3) - handles.SEG;
        
        MASK2 = repmat(handles.SEG,1,1,1,3);
        MASK2(:,:,:,1) = MASK2(:,:,:,1).*val_Lrgb_R;
        MASK2(:,:,:,2) = MASK2(:,:,:,2).*val_Lrgb_G;
        MASK2(:,:,:,3) = MASK2(:,:,:,3).*val_Lrgb_B;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.Lrgb = MASK1 + MASK2;

        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    elseif handles.view_sac == 3
        axes(handles.axes1);
        k = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        
        p = impoint;
        BW1 = createMask(p,k);
        [r,c,~] = find(BW1==1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        val_Lrgb_R = handles.Lrgb(handles.slider_axes3,r(1),c(1),1);
        val_Lrgb_G = handles.Lrgb(handles.slider_axes3,r(1),c(1),2);
        val_Lrgb_B = handles.Lrgb(handles.slider_axes3,r(1),c(1),3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.MM = zeros(size(handles.IPCMRA));% Julio Sotelo
        handles.MM(handles.slider_axes3,:,:) = double(BW1);% Julio Sotelo
        L = double(bwlabeln(handles.SEG,6));% Julio Sotelo
        LM = handles.MM.*double(bwlabeln(handles.SEG,6));
        val = max(LM(:));% Julio Sotelo
        handles.SEG = imfill(double(L==val),8,'holes');% Julio Sotelo
        
        
        MASK1 = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
        MASK1(:,:,:,1) = MASK1(:,:,:,1) - handles.SEG;
        MASK1(:,:,:,2) = MASK1(:,:,:,2) - handles.SEG;
        MASK1(:,:,:,3) = MASK1(:,:,:,3) - handles.SEG;
        
        MASK2 = repmat(handles.SEG,1,1,1,3);
        MASK2(:,:,:,1) = MASK2(:,:,:,1).*val_Lrgb_R;
        MASK2(:,:,:,2) = MASK2(:,:,:,2).*val_Lrgb_G;
        MASK2(:,:,:,3) = MASK2(:,:,:,3).*val_Lrgb_B;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.Lrgb = MASK1 + MASK2;
        
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.SEG = double(sum(handles.Lrgb,4)>0 & sum(handles.Lrgb,4)<3);
    end
guidata(hObject, handles);
function pushbutton6_Callback(hObject, eventdata, handles)
    if handles.id ==3
        for n=1:size(handles.veset,3)
            slider_step(1) = 1/handles.d;
            slider_step(2) = 0.1;
            set(handles.slider1,'Visible','on','Value', n/handles.d,'sliderstep',slider_step,'max',1,'min',0)
            if handles.ref == 5
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,n),handles.veset(1:end,2,n),handles.veset(1:end,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'jet');
                [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
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
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
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
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow)])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 6
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,n),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,n),handles.WSS(:,2,n),handles.WSS(:,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'hot');
                [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes1);
                c.LimitsMode = 'manual';
                handles.min_wss = min(handles.mags_wss(:));
                handles.max_wss = max(handles.mags_wss(:));
                handles.mean_wss = (handles.min_wss + handles.max_wss)/2;
                c.Limits = [handles.min_wss handles.max_wss];
                c.Ticks = [handles.min_wss, (handles.min_wss + handles.mean_wss)/2, handles.mean_wss, (handles.max_wss + handles.mean_wss)/2, handles.max_wss];
                c.TickLabels = {num2str(handles.min_wss,'%0.2f'), num2str((handles.min_wss + handles.mean_wss)/2,'%0.2f'), num2str(handles.mean_wss,'%0.2f'), num2str((handles.mean_wss + handles.max_wss)/2,'%0.2f'), num2str(handles.max_wss,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'WSS [N/m^{2}]';
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
                caxis(handles.axes1, [handles.min_wss handles.max_wss]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
                    currentColormap = colormap(handles.axes1,'hot');
                    [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
                    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
                    handles.cmap(:,:,4) = 255;
                    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                    c = colorbar(handles.axes1);
                    c.LimitsMode = 'manual';
                    handles.min_wss = min(handles.mags_wss(:));
                    handles.max_wss = max(handles.mags_wss(:));
                    handles.mean_wss = (handles.min_wss + handles.max_wss)/2;
                    c.Limits = [handles.min_wss handles.max_wss];
                    c.Ticks = [handles.min_wss, (handles.min_wss + handles.mean_wss)/2, handles.mean_wss, (handles.max_wss + handles.mean_wss)/2, handles.max_wss];
                    c.TickLabels = {num2str(handles.min_wss,'%0.2f'), num2str((handles.min_wss + handles.mean_wss)/2,'%0.2f'), num2str(handles.mean_wss,'%0.2f'), num2str((handles.mean_wss + handles.max_wss)/2,'%0.2f'), num2str(handles.max_wss,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'WSS [N/m^{2}]';
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
                    caxis(handles.axes1, [handles.min_wss handles.max_wss]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 8
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,n),handles.VOR(1:end,2,n),handles.VOR(1:end,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'cool');
                [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes1);
                handles.min_vor = min(handles.mags_vor(:));
                handles.max_vor = max(handles.mags_vor(:));
                handles.mean_vor = (handles.min_vor + handles.max_vor)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_vor handles.max_vor];
                c.Ticks = [handles.min_vor, (handles.min_vor + handles.mean_vor)/2, handles.mean_vor, (handles.max_vor + handles.mean_vor)/2, handles.max_vor];
                c.TickLabels = {num2str(handles.min_vor,'%0.2f'), num2str((handles.min_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.mean_vor,'%0.2f'), num2str((handles.max_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.max_vor,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Vorticity [1/s]';
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
                caxis(handles.axes1, [handles.min_vor handles.max_vor]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
                    hold on
                    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
                    currentColormap = colormap(handles.axes1,'cool');
                    [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
                    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                    handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
                    handles.cmap(:,:,4) = 255;
                    handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                    set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                    set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                    c = colorbar(handles.axes1);
                    handles.min_vor = min(handles.mags_vor(:));
                    handles.max_vor = max(handles.mags_vor(:));
                    handles.mean_vor = (handles.min_vor + handles.max_vor)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_vor handles.max_vor];
                    c.Ticks = [handles.min_vor, (handles.min_vor + handles.mean_vor)/2, handles.mean_vor, (handles.max_vor + handles.mean_vor)/2, handles.max_vor];
                    c.TickLabels = {num2str(handles.min_vor,'%0.2f'), num2str((handles.min_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.mean_vor,'%0.2f'), num2str((handles.max_vor + handles.mean_vor)/2,'%0.2f'), num2str(handles.max_vor,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Vorticity [1/s]';
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
                    caxis(handles.axes1, [handles.min_vor handles.max_vor]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 9
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,n),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'cool');
                c = colorbar(handles.axes1);
                handles.min_hd = min(handles.mags_hd(:));
                handles.max_hd = max(handles.mags_hd(:));
                handles.mean_hd = (handles.min_hd + handles.max_hd)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_hd handles.max_hd];
                c.Ticks = [handles.min_hd, (handles.min_hd + handles.mean_hd)/2, handles.mean_hd, (handles.max_hd + handles.mean_hd)/2, handles.max_hd];
                c.TickLabels = {num2str(handles.min_hd,'%0.2f'), num2str((handles.min_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.mean_hd,'%0.2f'), num2str((handles.max_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.max_hd,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Helicity Density [m/s^{2}]';
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
                caxis(handles.axes1, [handles.min_hd handles.max_hd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'cool');
                    c = colorbar(handles.axes1);
                    handles.min_hd = min(handles.mags_hd(:));
                    handles.max_hd = max(handles.mags_hd(:));
                    handles.mean_hd = (handles.min_hd + handles.max_hd)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_hd handles.max_hd];
                    c.Ticks = [handles.min_hd, (handles.min_hd + handles.mean_hd)/2, handles.mean_hd, (handles.max_hd + handles.mean_hd)/2, handles.max_hd];
                    c.TickLabels = {num2str(handles.min_hd,'%0.2f'), num2str((handles.min_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.mean_hd,'%0.2f'), num2str((handles.max_hd + handles.mean_hd)/2,'%0.2f'), num2str(handles.max_hd,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Helicity Density [m/s^{2}]';
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
                    caxis(handles.axes1, [handles.min_hd handles.max_hd]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 10
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,n),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'summer');
                c = colorbar(handles.axes1);
                handles.min_rhd = min(handles.mags_rhd(:));
                handles.max_rhd = max(handles.mags_rhd(:));
                handles.mean_rhd = (handles.min_rhd + handles.max_rhd)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_rhd handles.max_rhd];
                c.Ticks = [handles.min_rhd, (handles.min_rhd + handles.mean_rhd)/2, handles.mean_rhd, (handles.max_rhd + handles.mean_rhd)/2, handles.max_rhd];
                c.TickLabels = {num2str(handles.min_rhd,'%0.2f'), num2str((handles.min_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.mean_rhd,'%0.2f'), num2str((handles.max_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.max_rhd,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Relative Helicity Density [-]';
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
                caxis(handles.axes1, [handles.min_rhd handles.max_rhd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'summer');
                    c = colorbar(handles.axes1);
                    handles.min_rhd = min(handles.mags_rhd(:));
                    handles.max_rhd = max(handles.mags_rhd(:));
                    handles.mean_rhd = (handles.min_rhd + handles.max_rhd)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_rhd handles.max_rhd];
                    c.Ticks = [handles.min_rhd, (handles.min_rhd + handles.mean_rhd)/2, handles.mean_rhd, (handles.max_rhd + handles.mean_rhd)/2, handles.max_rhd];
                    c.TickLabels = {num2str(handles.min_rhd,'%0.2f'), num2str((handles.min_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.mean_rhd,'%0.2f'), num2str((handles.max_rhd + handles.mean_rhd)/2,'%0.2f'), num2str(handles.max_rhd,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Relative Helicity Density [-]';
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
                    caxis(handles.axes1, [handles.min_rhd handles.max_rhd]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 11
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,n),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'winter');
                c = colorbar(handles.axes1);
                handles.min_vd = min(handles.mags_vd(:));
                handles.max_vd = max(handles.mags_vd(:));
                handles.mean_vd = (handles.min_vd + handles.max_vd)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_vd handles.max_vd];
                c.Ticks = [handles.min_vd, (handles.min_vd + handles.mean_vd)/2, handles.mean_vd, (handles.max_vd + handles.mean_vd)/2, handles.max_vd];
                c.TickLabels = {num2str(handles.min_vd,'%0.2f'), num2str((handles.min_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.mean_vd,'%0.2f'), num2str((handles.max_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.max_vd,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Viscous Dissipation [1e^{3}/s^{2}]'; 
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
                caxis(handles.axes1, [handles.min_vd handles.max_vd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'winter');
                    c = colorbar(handles.axes1);
                    handles.min_vd = min(handles.mags_vd(:));
                    handles.max_vd = max(handles.mags_vd(:));
                    handles.mean_vd = (handles.min_vd + handles.max_vd)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_vd handles.max_vd];
                    c.Ticks = [handles.min_vd, (handles.min_vd + handles.mean_vd)/2, handles.mean_vd, (handles.max_vd + handles.mean_vd)/2, handles.max_vd];
                    c.TickLabels = {num2str(handles.min_vd,'%0.2f'), num2str((handles.min_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.mean_vd,'%0.2f'), num2str((handles.max_vd + handles.mean_vd)/2,'%0.2f'), num2str(handles.max_vd,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Viscous Dissipation [1e^{3}/s^{2}]'; 
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
                    caxis(handles.axes1, [handles.min_vd handles.max_vd]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 12
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,n),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'autumn');
                c = colorbar(handles.axes1);
                handles.min_el = min(handles.mags_el(:));
                handles.max_el = max(handles.mags_el(:));
                handles.mean_el = (handles.min_el + handles.max_el)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_el handles.max_el];
                c.Ticks = [handles.min_el, (handles.min_el + handles.mean_el)/2, handles.mean_el, (handles.max_el + handles.mean_el)/2, handles.max_el];
                c.TickLabels = {num2str(handles.min_el,'%0.2f'), num2str((handles.min_el + handles.mean_el)/2,'%0.2f'), num2str(handles.mean_el,'%0.2f'), num2str((handles.max_el + handles.mean_el)/2,'%0.2f'), num2str(handles.max_el,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Energy Loss [\muW]';
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
                caxis(handles.axes1, [handles.min_el handles.max_el]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'autumn');
                    c = colorbar(handles.axes1);
                    handles.min_el = min(handles.mags_el(:));
                    handles.max_el = max(handles.mags_el(:));
                    handles.mean_el = (handles.min_el + handles.max_el)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_el handles.max_el];
                    c.Ticks = [handles.min_el, (handles.min_el + handles.mean_el)/2, handles.mean_el, (handles.max_el + handles.mean_el)/2, handles.max_el];
                    c.TickLabels = {num2str(handles.min_el,'%0.2f'), num2str((handles.min_el + handles.mean_el)/2,'%0.2f'), num2str(handles.mean_el,'%0.2f'), num2str((handles.max_el + handles.mean_el)/2,'%0.2f'), num2str(handles.max_el,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Energy Loss [\muW]';
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
                    caxis(handles.axes1, [handles.min_el handles.max_el]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            elseif handles.ref == 13
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,n),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'spring');
                c = colorbar(handles.axes1);
                handles.min_ke = min(handles.mags_ke(:));
                handles.max_ke = max(handles.mags_ke(:));
                handles.mean_ke = (handles.min_ke + handles.max_ke)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_ke handles.max_ke];
                c.Ticks = [handles.min_ke, (handles.min_ke + handles.mean_ke)/2, handles.mean_ke, (handles.max_ke + handles.mean_ke)/2, handles.max_ke];
                c.TickLabels = {num2str(handles.min_ke,'%0.2f'), num2str((handles.min_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.mean_ke,'%0.2f'), num2str((handles.max_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.max_ke,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Kinetic Energy [\muJ]';
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
                caxis(handles.axes1, [handles.min_ke handles.max_ke]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                if n==size(handles.veset,3)
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'spring');
                    c = colorbar(handles.axes1);
                    handles.min_ke = min(handles.mags_ke(:));
                    handles.max_ke = max(handles.mags_ke(:));
                    handles.mean_ke = (handles.min_ke + handles.max_ke)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_ke handles.max_ke];
                    c.Ticks = [handles.min_ke, (handles.min_ke + handles.mean_ke)/2, handles.mean_ke, (handles.max_ke + handles.mean_ke)/2, handles.max_ke];
                    c.TickLabels = {num2str(handles.min_ke,'%0.2f'), num2str((handles.min_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.mean_ke,'%0.2f'), num2str((handles.max_ke + handles.mean_ke)/2,'%0.2f'), num2str(handles.max_ke,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'Kinetic Energy [\muJ]';
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
                    caxis(handles.axes1, [handles.min_ke handles.max_ke]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
               
            elseif handles.ref==20 % WSSA   
        
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,n)),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,n),handles.WSS_A(:,2,n),handles.WSS_A(:,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'hot');
                [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
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
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'WSS-Axial [N/m^{2}]';
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
                caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                
                if n==size(handles.veset,3)
                    
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
                    hold on
                    q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
                    currentColormap = colormap(handles.axes1,'hot');
                    [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
                    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'WSS-Axial [N/m^{2}]';
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
                    caxis(handles.axes1, [handles.min_wssa handles.max_wssa]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end

            elseif handles.ref==21 % WSSC

                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,n)),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,n),handles.WSS_C(1:end,2,n),handles.WSS_C(1:end,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'hot');
                [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
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
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'WSS-Circum. [N/m^{2}]';
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
                caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);
                
                if n==size(handles.veset,3)
                    
                    axes(handles.axes1);
                    plot(0.0)   
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
                    hold on
                    q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
                    currentColormap = colormap(handles.axes1,'hot');
                    [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
                    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                    handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
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
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
                    c.FontWeight = 'bold';
                    c.Label.String = 'WSS-Circum. [N/m^{2}]';
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
                    caxis(handles.axes1, [handles.min_wssc handles.max_wssc]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end

            elseif handles.ref==22 % Axial Angle

                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,n)),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes1,'cool');
                c = colorbar(handles.axes1);
                handles.min_axa = min(handles.angle_axial_direction(:));
                handles.max_axa = max(handles.angle_axial_direction(:));
                handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_axa handles.max_axa];
                c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
                c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_axa handles.max_axa]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);

                if n==size(handles.veset,3)
                    
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
                    hold on
                    colormap(handles.axes1,'cool');
                    c = colorbar(handles.axes1);
                    handles.min_axa = min(handles.angle_axial_direction(:));
                    handles.max_axa = max(handles.angle_axial_direction(:));
                    handles.mean_axa = (handles.min_axa + handles.max_axa)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_axa handles.max_axa];
                    c.Ticks = [handles.min_axa, (handles.min_axa + handles.mean_axa)/2, handles.mean_axa, (handles.max_axa + handles.mean_axa)/2, handles.max_axa];
                    c.TickLabels = {num2str(handles.min_axa,'%0.2f'), num2str((handles.min_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.mean_axa,'%0.2f'), num2str((handles.max_axa + handles.mean_axa)/2,'%0.2f'), num2str(handles.max_axa,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    c.Label.Color = [1 1 1];
                    caxis(handles.axes1, [handles.min_axa handles.max_axa]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
                

            elseif handles.ref==23 % Forward Velocity


                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,n),handles.forward_velocity(1:end,2,n),handles.forward_velocity(1:end,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'jet');
                [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
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
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);

                if n==size(handles.veset,3)
                    
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
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    c.Label.Color = [1 1 1];
                    caxis(handles.axes1, [handles.min_fvel handles.max_fvel]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
                
            elseif handles.ref==24 % Backward Velocity    

                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,n),handles.backward_velocity(1:end,2,n),handles.backward_velocity(1:end,3,n),5,'Linewidth',1);
                currentColormap = colormap(handles.axes1,'jet');
                [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
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
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);

                if n==size(handles.veset,3)
                    
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
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    c.Label.Color = [1 1 1];
                    caxis(handles.axes1, [handles.min_bvel handles.max_bvel]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
                
%             elseif handles.ref==31 % circulation
%             
%                 axes(handles.axes1);
%                 plot(0.0)
%                 patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,n),'CDataMapping','Scaled')                
%                 hold on
%                 colormap(handles.axes1,'winter');
%                 c = colorbar(handles.axes1);
%                 handles.min_cir = min(handles.circulation(:));
%                 handles.max_cir = max(handles.circulation(:));
%                 handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%                 c.LimitsMode = 'manual';
%                 c.Limits = [handles.min_cir handles.max_cir];
%                 c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%                 c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%                 c.Color = [1 1 1];
%                 c.Location = 'manual';
%                 c.Position = [0.2 0.1 0.02 0.5];
%                 c.FontWeight = 'bold';
%                 c.Label.String = 'Circulation [mm^{2}/s]';
%                 windows_screen_size = get(0,'ScreenSize');
%                 if windows_screen_size(4)<=900
%                     c.FontSize = 9;
%                 elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                     c.FontSize = 10;
%                 else
%                     c.FontSize = 11;
%                 end
%                 if windows_screen_size(4)<=900
%                     c.Label.FontSize = 10;
%                 elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                     c.Label.FontSize = 12;
%                 else
%                     c.Label.FontSize = 14;
%                 end
%                 c.Label.FontWeight = 'bold';
%                 c.Label.Color = [1 1 1];
%                 caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%                 hold off
%                 axis vis3d
%                 daspect([1,1,1])
%                 axis off
%                 view([-34,-51])
%                 set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
%                 pause(0.05);
% 
%                 if n==size(handles.veset,3)
%                     
%                     axes(handles.axes1);
%                     plot(0.0)
%                     patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,handles.peak_flow),'CDataMapping','Scaled')
%                     hold on
%                     colormap(handles.axes1,'winter');
%                     c = colorbar(handles.axes1);
%                     handles.min_cir = min(handles.circulation(:));
%                     handles.max_cir = max(handles.circulation(:));
%                     handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%                     c.LimitsMode = 'manual';
%                     c.Limits = [handles.min_cir handles.max_cir];
%                     c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%                     c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%                     c.Color = [1 1 1];
%                     c.Location = 'manual';
%                     c.Position = [0.2 0.1 0.02 0.5];
%                     c.FontWeight = 'bold';
%                     c.Label.String = 'Circulation [mm^{2}/s]';
%                     windows_screen_size = get(0,'ScreenSize');
%                     if windows_screen_size(4)<=900
%                         c.FontSize = 9;
%                     elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                         c.FontSize = 10;
%                     else
%                         c.FontSize = 11;
%                     end
%                     if windows_screen_size(4)<=900
%                         c.Label.FontSize = 10;
%                     elseif windows_screen_size(4)>900 && windows_screen_size(4)<=1100
%                         c.Label.FontSize = 12;
%                     else
%                         c.Label.FontSize = 14;
%                     end
%                     c.Label.FontWeight = 'bold';
%                     c.Label.Color = [1 1 1];
%                     caxis(handles.axes1, [handles.min_cir handles.max_cir]);
%                     hold off
%                     axis vis3d
%                     daspect([1,1,1])
%                     axis off
%                     view([-34,-51])
%                     set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
%                     set(handles.slider1,'Value', handles.peak_flow/handles.d)
%                 end

            elseif handles.ref==31 % forward_vortex

                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,n),'CDataMapping','Scaled')                
                hold on
                colormap(handles.axes1,'cool');
                c = colorbar(handles.axes1);
                handles.min_fov = min(handles.forward_vortex(:));
                handles.max_fov = max(handles.forward_vortex(:));
                handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_fov handles.max_fov];
                c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
                c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_fov handles.max_fov]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);

                if n==size(handles.veset,3)
                    
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')                    
                    hold on
                    colormap(handles.axes1,'cool');
                    c = colorbar(handles.axes1);
                    handles.min_fov = min(handles.forward_vortex(:));
                    handles.max_fov = max(handles.forward_vortex(:));
                    handles.mean_fov = (handles.min_fov + handles.max_fov)/2;
                    c.LimitsMode = 'manual';
                    c.Limits = [handles.min_fov handles.max_fov];
                    c.Ticks = [handles.min_fov, (handles.min_fov + handles.mean_fov)/2, handles.mean_fov, (handles.max_fov + handles.mean_fov)/2, handles.max_fov];
                    c.TickLabels = {num2str(handles.min_fov,'%0.2f'), num2str((handles.min_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.mean_fov,'%0.2f'), num2str((handles.max_fov + handles.mean_fov)/2,'%0.2f'), num2str(handles.max_fov,'%0.2f')};
                    c.Color = [1 1 1];
                    c.Location = 'manual';
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    c.Label.Color = [1 1 1];
                    caxis(handles.axes1, [handles.min_fov handles.max_fov]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
                
            elseif handles.ref==34 % axialcirculation
            
                axes(handles.axes1);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,n),'CDataMapping','Scaled')                
                hold on
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
                c.Position = [0.2 0.1 0.02 0.5];
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
                c.Label.Color = [1 1 1];
                caxis(handles.axes1, [handles.min_acir handles.max_acir]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(n)])
                pause(0.05);

                if n==size(handles.veset,3)
                    
                    axes(handles.axes1);
                    plot(0.0)
                    patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,handles.peak_flow),'CDataMapping','Scaled')
                    hold on
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
                    c.Position = [0.2 0.1 0.02 0.5];
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
                    c.Label.Color = [1 1 1];
                    caxis(handles.axes1, [handles.min_acir handles.max_acir]);
                    hold off
                    axis vis3d
                    daspect([1,1,1])
                    axis off
                    view([-34,-51])
                    set(handles.text1,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                    set(handles.slider1,'Value', handles.peak_flow/handles.d)
                end
            end
        end
    else
        handles.SEG = handles.SEG_old;
        handles.Lrgb = handles.Lrgb_old;
        handles.L = handles.L_old;
        handles.NUM = handles.NUM_old;
        if handles.view_sac == 1
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        elseif handles.view_sac == 2
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        elseif handles.view_sac == 3
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
            hold on
            Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        end
    end
guidata(hObject, handles);
function pushbutton7_Callback(hObject, eventdata, handles)
    input.id            = 2;
    input.IPCMRA        = handles.IPCMRA;
    input.SEG           = handles.SEG;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.NUM           = handles.NUM;
    input.xd            = handles.xd;
    input.yd            = handles.yd;
    input.zd            = handles.zd;
    input.a             = handles.a;
    input.b             = handles.b;
    input.c             = handles.c;
    input.slider_axes1  = handles.slider_axes1;
    input.slider_axes2  = handles.slider_axes2;
    input.slider_axes3  = handles.slider_axes3;
    input.voxel_MR      = handles.voxel_MR;
    GUIDE_THRESHOLDING(input)
    handles.SEG = getappdata(0,'SEG');
    handles.L = getappdata(0,'L');
    handles.NUM = getappdata(0,'NUM');
    handles.Lrgb = getappdata(0,'Lrgb');
    if handles.view_sac == 1
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    elseif handles.view_sac == 2
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    elseif handles.view_sac == 3
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
guidata(hObject, handles);
function pushbutton8_Callback(hObject, eventdata, handles)
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.Lrgb = label2rgb3d(L,'jet');
    if handles.view_sac == 1
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    elseif handles.view_sac == 2
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    elseif handles.view_sac == 3
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)));
        hold on
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
guidata(hObject, handles);
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    handles.id_while = 1;
    handles.output = hObject;
    guidata(hObject, handles);
    uiresume(hObject);
else
    delete(hObject);
end
function popupmenu1_Callback(hObject, eventdata, handles)
switch get(handles.popupmenu1,'Value')
    case 1
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.view_sac = 1;
        slider_step(1) = 1/handles.c;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
    case 2
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.view_sac = 2;
        slider_step(1) = 1/handles.b;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
    case 3
        axes(handles.axes1);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        handles.view_sac = 3;
        slider_step(1) = 1/handles.a;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
end
guidata(hObject, handles);

function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_SizeChangedFcn(hObject, eventdata, handles)

    % change size of figure buttons 

    set(handles.pushbutton1, 'Units', 'pixels');
    handles.pushbutton1_size = get(handles.pushbutton1, 'Position');
    set(handles.pushbutton1, 'Units', 'normalized');
    idx = mod(min(handles.pushbutton1_size(3:4)),2)>1;
    w = floor(min(handles.pushbutton1_size(3:4)));
    w(idx) = w(idx)+1;
    IM = imread('Symbols/P_Cut.tiff');
    IM_Cut_rgb1 = double(IM(:,:,1:3)<=1);
    IM_Cut_rgb1(:,:,2) = IM_Cut_rgb1(:,:,2)*0;
    IM_Cut_rgb1(:,:,3) = IM_Cut_rgb1(:,:,3)*0;
    f1=double(imresize(IM_Cut_rgb1,[w-20 w-20],'method','nearest')>0);
    IM = imread('Symbols/P_Magic.tiff');
    IM_Magic_rgb1 = double(IM(:,:,1:3)<=1);
    IM_Magic_rgb1(:,:,2) = IM_Magic_rgb1(:,:,2)*0;
    IM_Magic_rgb1(:,:,3) = IM_Magic_rgb1(:,:,3)*0;
    f2=double(imresize(IM_Magic_rgb1,[w-20 w-20],'method','nearest')>0);
    IM = imread('Symbols/P_Ellipse.tiff');
    IM_Ellipse_rgb1 = double(IM(:,:,1:3)<=1);
    IM_Ellipse_rgb1(:,:,2) = IM_Ellipse_rgb1(:,:,2)*0;
    IM_Ellipse_rgb1(:,:,3) = IM_Ellipse_rgb1(:,:,3)*0;
    f3=double(imresize(IM_Ellipse_rgb1,[w-20 w-20],'method','nearest')>0);
    IM = imread('Symbols/P_Cut.tiff');
    IM_Cut_rgb2 = double(IM(:,:,1:3)<=1);
    IM_Cut_rgb2(:,:,1) = IM_Cut_rgb2(:,:,1)*0;
    IM_Cut_rgb2(:,:,3) = IM_Cut_rgb2(:,:,3)*0;
    f4=double(imresize(IM_Cut_rgb2,[w-20 w-20],'method','nearest')>0);
    IM = imread('Symbols/P_Fill.tiff');
    IM_Ellipse_rgb2 = double(IM(:,:,1:3)<=1);
    IM_Ellipse_rgb2(:,:,1) = IM_Ellipse_rgb2(:,:,1)*0;
    IM_Ellipse_rgb2(:,:,3) = IM_Ellipse_rgb2(:,:,3)*0;
    f5=double(imresize(IM_Ellipse_rgb2,[w-20 w-20],'method','nearest')>0);
    set(handles.pushbutton1,'CData',f1);
    set(handles.pushbutton2,'CData',f2);
    set(handles.pushbutton3,'CData',f3);
    set(handles.pushbutton4,'CData',f4);
    set(handles.pushbutton5,'CData',f5);

    % change size windows
    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');
    
    set(handles.text1,'FontUnits','Normalized','FontSize',0.58)
    set(handles.pushbutton6,'FontUnits','Normalized','FontSize',0.14)
    set(handles.pushbutton7,'FontUnits','Normalized','FontSize',0.14)
    set(handles.pushbutton8,'FontUnits','Normalized','FontSize',0.14)
    set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.47)
    set(handles.popupmenu2,'FontUnits','Normalized','FontSize',0.47)

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.popupmenu2,'Value')
    case 1
        handles.IPCMRA = handles.ANG; % Julio Sotelo
        handles.idangmag = 1;
        if handles.view_sac == 1
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 1;
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])
        elseif handles.view_sac == 2
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 2;
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes2)])
        elseif handles.view_sac == 3
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 3;
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes3)])
        end
    case 2
        handles.IPCMRA = handles.MAG; % Julio Sotelo
        handles.idangmag = 2;
        if handles.view_sac == 1
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 1;
            slider_step(1) = 1/handles.c;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes1)])
        elseif handles.view_sac == 2
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 2;
            slider_step(1) = 1/handles.b;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes2)])
        elseif handles.view_sac == 3
            axes(handles.axes1);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vel == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vor == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_hd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_rhd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_vd == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_el == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_ke == 1 && handles.id_seg == 0
                Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            handles.view_sac = 3;
            slider_step(1) = 1/handles.a;
            slider_step(2) = 0.1;
            set(handles.slider1,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
            set(handles.text1,'String',['# Slice: ',num2str(handles.slider_axes3)])
        end
        
end
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


% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
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
