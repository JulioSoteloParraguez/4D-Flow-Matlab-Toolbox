function varargout = GUIDE_4D_FLOW(varargin)
% GUIDE_4D_FLOW MATLAB code for GUIDE_4D_FLOW.fig
% This code is develop by Dr. Julio Sotelo, that together with 
% Dr. Sergio Uribe and Dr. Daniel Hurtado, have been working in Cardiac MR 
% and particularly in the quantification of 4D flow data for about 7 years 
% iaddpath(genpath('vWERP/'));n the Biomedical Imaging center at Universidad Catolica of Chile
% (www.mri.cl).
% We have developed a methodology for the non-invasive quantification of 
% hemodynamics from 4D flow data sets based on Finite Elements methods. 
% This technique is unique and is possible obtain several hemodynamic 
% parameters in 3D as WSS, OSI, vorticity, helicity density, relative helicity
% density, viscouss dissipation, energy loss and kinetic energy. 
%   Author:      Dr. Julio Sotelo
%   Time-stamp:  2019-06-06
%   E-mail:      jasotelo@uc.cl
gui_Singleton = 1; 
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_4D_FLOW_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_4D_FLOW_OutputFcn, ...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GUIDE_4D_FLOW_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Showing the logo image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    txt = 'Symbols/LOGO.tiff';
    imlogo = imread(txt);
    [av,bv,~] = size(imlogo);
    windows_screen_size = get(0,'ScreenSize');
    imlogo = imresize(imlogo,[round(windows_screen_size(4)*av/(av+bv)) round(windows_screen_size(4)*bv/(av+bv))]);
    Sz= size(imlogo);
    flogo = figure('Position',[10000 10000 Sz(2) + 4 Sz(1) + 4],'name','4D FLOW APP','numbertitle','off','menubar','none');
    movegui(flogo,'center');
    set(flogo,'Units', 'pixels');
    image(imlogo(:,:,1:3))
    set(gca,'Visible','off','Units','pixels','Position', [2 2 Sz(2) Sz(1)]);
    pause(5)
    delete(flogo)
    clear imlogo
    warning off
    
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = GUIDE_4D_FLOW_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton1_Callback(hObject, eventdata, handles)
    
    set(handles.uipanel1,'Visible','on')
    set(handles.uipanel6,'Visible','off')
    set(handles.pushbutton1,'BackgroundColor',[0 0 0])
    set(handles.pushbutton15,'BackgroundColor',[0.2 0.2 0.2])

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton15_Callback(hObject, eventdata, handles)

    % change the color of the button and change the position of the panel
    pan1pos=get(handles.uipanel1,'Position');
    set(handles.pushbutton1,'BackgroundColor',[0.2 0.2 0.2])
    set(handles.pushbutton15,'BackgroundColor',[0 0 0])
    set(handles.uipanel1,'Visible','off')
    set(handles.uipanel6,'Position',pan1pos,'Visible','off')
    
    
    % change the velocity values between +pi and -pi
    handles.id_unwrappping = 1;
    handles.PHASE_FH = (handles.MR_PCA_FH*pi)/handles.VENC;
    handles.PHASE_AP = (handles.MR_PCA_AP*pi)/handles.VENC;
    handles.PHASE_RL = (handles.MR_PCA_RL*pi)/handles.VENC;
    
    handles.dt = (60/handles.heart_rate)/(size(handles.MR_PCA_FH,4)-1);
    
    % we use the segmentation as mask 
    if sum(handles.SEG(:))==0
        waitfor(warndlg('Segmentation of the vessel is recommended ...','Warning'))
        h = waitbar(0,'Loading files ...');
        [X, Y, Z, T] = ndgrid((0:size(handles.PHASE_FH,1)-1)-size(handles.PHASE_FH,1)/2, (0:size(handles.PHASE_FH,2)-1)-size(handles.PHASE_FH,2)/2, (0:size(handles.PHASE_FH,3)-1)-size(handles.PHASE_FH,3)/2, (0:size(handles.PHASE_FH,4)-1)-size(handles.PHASE_FH,4)/2);
        mod1 = 2.*cos(pi*X./size(handles.PHASE_FH,1)) + 2.*cos(pi*Y./size(handles.PHASE_FH,2)) + 2.*cos(pi*Z./size(handles.PHASE_FH,3)) + 2.*cos(pi*T./size(handles.PHASE_FH,4))  - 8;
        mod2 = mod1;
        mod2(floor(size(handles.PHASE_FH,1)/2)+1,floor(size(handles.PHASE_FH,2)/2)+1,floor(size(handles.PHASE_FH,3)/2)+1,floor(size(handles.PHASE_FH,4)/2)+1) = 1;
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_FH))).*mod1))));
        LAP2 = cos(handles.PHASE_FH).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_FH)))).*mod1)))) - sin(handles.PHASE_FH).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_FH)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_FH = LAP3;
        waitbar(0.33,h,'Loading files ...');
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_AP))).*mod1))));
        LAP2 = cos(handles.PHASE_AP).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_AP)))).*mod1)))) - sin(handles.PHASE_AP).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_AP)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_AP = LAP3;
        waitbar(0.66,h,'Loading files ...');
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_RL))).*mod1))));
        LAP2 = cos(handles.PHASE_RL).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_RL)))).*mod1)))) - sin(handles.PHASE_RL).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_RL)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_RL = LAP3;
        waitbar(0.99,h,'Loading files ...');
        pause(0.0001)
        close(h)
    else
        h = waitbar(0,'Loading files ...');
        IM_process = permute(handles.SEG,[1,3,2]);
        IM_p = zeros(size(IM_process,1)+2,size(IM_process,2)+2,size(IM_process,3)+2);
        IM_p(2:end-1,2:end-1,2:end-1) = IM_process;
        var = squeeze(double(sum(sum(IM_p))>0));
        [r,c,v] = find(var ==1);
        slices_s = [r(1),r(end)];
        kernel = zeros(3,3,3);
        for ii=1:size(kernel,1)
            for jj=1:size(kernel,2)
                for kk=1:size(kernel,3)
                    if sqrt(sum(([ii,jj,kk]-[round(size(kernel,1)/2),round(size(kernel,2)/2),round(size(kernel,3)/2)]).^2))<=floor(size(kernel,3)/2)
                        kernel(ii,jj,kk) = 1;
                    end
                end
            end
        end
        waitbar(0.25,h,'Loading files ...');
        IM_OUT = imdilate(IM_p,kernel);
        IM_OUT(:,:,slices_s(1)) = imdilate(IM_p(:,:,slices_s(1)),kernel(:,:,round(size(kernel,3)/2)));
        IM_OUT(:,:,slices_s(2)) = imdilate(IM_p(:,:,slices_s(2)),kernel(:,:,round(size(kernel,3)/2)));
        IM_OUT = IM_OUT(2:end-1,2:end-1,2:end-1);
        IM_OUT = permute(IM_OUT,[1,3,2]);
        handles.SEG_Aorta = IM_OUT;
        [X, Y, Z, T] = ndgrid((0:size(handles.PHASE_FH,1)-1)-size(handles.PHASE_FH,1)/2, (0:size(handles.PHASE_FH,2)-1)-size(handles.PHASE_FH,2)/2, (0:size(handles.PHASE_FH,3)-1)-size(handles.PHASE_FH,3)/2, (0:size(handles.PHASE_FH,4)-1)-size(handles.PHASE_FH,4)/2);
        mod1 = 2.*cos(pi*X./size(handles.PHASE_FH,1)) + 2.*cos(pi*Y./size(handles.PHASE_FH,2)) + 2.*cos(pi*Z./size(handles.PHASE_FH,3)) + 2.*cos(pi*T./size(handles.PHASE_FH,4))  - 8;
        mod2 = mod1;
        mod2(floor(size(handles.PHASE_FH,1)/2)+1,floor(size(handles.PHASE_FH,2)/2)+1,floor(size(handles.PHASE_FH,3)/2)+1,floor(size(handles.PHASE_FH,4)/2)+1) = 1;
       
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_FH))).*mod1))));
        LAP2 = cos(handles.PHASE_FH).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_FH)))).*mod1)))) - sin(handles.PHASE_FH).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_FH)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_FH = LAP3.*repmat(handles.SEG,[1 1 1 size(LAP3,4)]);
        waitbar(0.5,h,'Loading files ...');
        
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_AP))).*mod1))));
        LAP2 = cos(handles.PHASE_AP).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_AP)))).*mod1)))) - sin(handles.PHASE_AP).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_AP)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_AP = LAP3.*repmat(handles.SEG,[1 1 1 size(LAP3,4)]);
        waitbar(0.75,h,'Loading files ...');
        
        LAP1 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(handles.PHASE_RL))).*mod1))));
        LAP2 = cos(handles.PHASE_RL).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(sin(handles.PHASE_RL)))).*mod1)))) - sin(handles.PHASE_RL).*real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(cos(handles.PHASE_RL)))).*mod1))));
        LAP3 = real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(LAP2 - LAP1)))./mod2))))./2./pi;
        handles.LAP_RL = LAP3.*repmat(handles.SEG,[1 1 1 size(LAP3,4)]);
        waitbar(0.99,h,'Loading files ...');

        pause(0.0001)
        close(h)
    end
    
    LAPFH = handles.LAP_FH(:);
    LAPFH(LAPFH>0) = LAPFH(LAPFH>0)/max(LAPFH);
    LAPFH(LAPFH<0) = (LAPFH(LAPFH<0)/min(LAPFH))*-1;
    LAPFH = reshape(LAPFH,[size(handles.PHASE_FH,1),size(handles.PHASE_FH,2),size(handles.PHASE_FH,3),size(handles.PHASE_FH,4)]);
    
    LAPAP = handles.LAP_AP(:);
    LAPAP(LAPAP>0) = LAPAP(LAPAP>0)/max(LAPAP);
    LAPAP(LAPAP<0) = (LAPAP(LAPAP<0)/min(LAPAP))*-1;
    LAPAP = reshape(LAPAP,[size(handles.PHASE_FH,1),size(handles.PHASE_FH,2),size(handles.PHASE_FH,3),size(handles.PHASE_FH,4)]);

    LAPRL = handles.LAP_RL(:);
    LAPRL(LAPRL>0) = LAPRL(LAPRL>0)/max(LAPRL);
    LAPRL(LAPRL<0) = (LAPRL(LAPRL<0)/min(LAPRL))*-1;
    LAPRL = reshape(LAPRL,[size(handles.PHASE_FH,1),size(handles.PHASE_FH,2),size(handles.PHASE_FH,3),size(handles.PHASE_FH,4)]);
    
    handles.LAP_FH = LAPFH;
    handles.LAP_AP = LAPAP;
    handles.LAP_RL = LAPRL;
    
    
    % we set both popupmenu to the default values
    handles.view_sac = 1;
    set(handles.popupmenu2,'Value',1)
    handles.flow_enc = 1;
    set(handles.popupmenu3,'Value',1)
    
    handles.FLOW_IM = handles.PHASE_FH;
    handles.LAP = handles.LAP_FH;
    handles.slider_axes_uwrap_1 = round(handles.c/2);
    handles.slider_axes_uwrap_2 = round(handles.d/2);
    
    % we show the FH direction as default image, in the center of cardiac
    % phases and the center of the slice, to adjust the sliders. 
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)))
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/handles.c;
    slider_step(2) = 0.1;
    set(handles.slider6,'Value', handles.slider_axes_uwrap_1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
    slider_step(1) = 1/handles.d;
    slider_step(2) = 0.1;
    set(handles.slider5,'Value', handles.slider_axes_uwrap_2/handles.d,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.text5,'String',['Cardiac Phase # ',num2str(handles.slider_axes_uwrap_2),' of ',num2str(handles.d)])
    set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.c)])

    slider_step(1) = 0.05;
    slider_step(2) = 0.1;
    set(handles.slider7,'Value', 0,'sliderstep',slider_step,'max',1,'min',0)
    slider_step(1) = 0.05;
    slider_step(2) = 0.1;
    set(handles.slider8,'Value', 0,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.text8,'String',['Value : ',num2str(0)])
    set(handles.text10,'String',['Value : ',num2str(0)])
    

    handles.slider_value_positive = 0; % defauls values for the slider to adjust the threshold
    handles.slider_value_negative = 0; % defauls values for the slider to adjust the threshold
    
    handles.SEG_positive = zeros(size(handles.LAP));
    handles.SEG_negative = zeros(size(handles.LAP));

    handles.LAP_SEG_positive = handles.LAP.*double(handles.LAP>0); % defauls values for the slider to adjust the threshold
    handles.LAP_SEG_negative = handles.LAP.*double(handles.LAP<0); % defauls values for the slider to adjust the threshold
     
    handles.LRGB_SEG_positive = repmat(ones(size(handles.SEG_positive)),[1 1 1 1 3]); % defauls colors for positive threshold
    handles.LRGB_SEG_negative = repmat(ones(size(handles.SEG_negative)),[1 1 1 1 3]); % defauls colors for negative threshold
    
    handles.PI_W_Positive = 1;  % defauls values weithed for positive unwrapping
    handles.PI_W_Negative = -1; % defauls values weithed for negative unwrapping
    
    set(handles.uipanel6,'Visible','on')
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton2_Callback(hObject, eventdata, handles)
    input.id            = 2;
    input.idangmag      = 1;
    input.IPCMRA        = handles.IPCMRA;
    input.SEG           = handles.SEG;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.NUM           = handles.NUM;
    input.min_value_th  = handles.min_value_th;
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
    input.view_sac      = 1;
    input.id_seg        = handles.id_seg;
    input.id_vel        = handles.id_vel;
    input.id_vor        = handles.id_vor;
    input.id_hd         = handles.id_hd;
    input.id_rhd        = handles.id_rhd;
    input.id_vd         = handles.id_vd;
    input.id_el         = handles.id_el;
    input.id_ke         = handles.id_ke;
    input.Lrgb_vel      = handles.Lrgb_vel;
    input.Lrgb_vor      = handles.Lrgb_vor;
    input.Lrgb_hd       = handles.Lrgb_hd;
    input.Lrgb_rhd      = handles.Lrgb_rhd;
    input.Lrgb_vd       = handles.Lrgb_vd;
    input.Lrgb_el       = handles.Lrgb_el;
    input.Lrgb_ke       = handles.Lrgb_ke;
    input.id_fve        = handles.id_fve;
    input.Lrgb_fve      = handles.Lrgb_fve;
    input.id_bve        = handles.id_bve;
    input.Lrgb_bve      = handles.Lrgb_bve;
    input.id_aan        = handles.id_aan;
    input.Lrgb_aan      = handles.Lrgb_aan;
    input.id_fov        = handles.id_fov;
    input.Lrgb_fov      = handles.Lrgb_fov;
    
    input.peak_flow      = handles.peak_flow;

    GUIDE_SEGMENTATION(input)
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton6_Callback(hObject, eventdata, handles)
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    h = imrect(gca,handles.pos1);
    setColor(h,'c');
    position = wait(h);
    
    %%% Julio Sotelo 23-11-2018
    handles.position_sag = position;
    handles.position_axi = [handles.position_axi(1) position(2) handles.position_axi(3) position(4)];
    handles.position_cor = [handles.position_cor(1) position(1) handles.position_cor(3) position(3)];

    handles.c1 = position(2);
    handles.c2 = position(1);
    handles.c3 = handles.pos2(1);
    handles.c4 = position(4);
    handles.c5 = position(3);
    handles.c6 = handles.pos2(3);
    handles.pos1 = [handles.c2 handles.c1 handles.c5 handles.c4];
    handles.pos2 = [handles.c3 handles.c1 handles.c6 handles.c4];
    handles.pos3 = [handles.c3 handles.c2 handles.c6 handles.c5];
    
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton9_Callback(hObject, eventdata, handles)

    
    handles.position_sag_cor = round([(handles.position_sag(1)/handles.voxel_MR(2))+1 (handles.position_sag(2)/handles.voxel_MR(1))+1 (handles.position_sag(3)/handles.voxel_MR(2))+1 (handles.position_sag(4)/handles.voxel_MR(1))+1]);
    handles.position_axi_cor = round([(handles.position_axi(1)/handles.voxel_MR(3))+1 (handles.position_axi(2)/handles.voxel_MR(1))+1 (handles.position_axi(3)/handles.voxel_MR(3))+1 (handles.position_axi(4)/handles.voxel_MR(1))+1]);
    handles.position_cor_cor = round([(handles.position_cor(1)/handles.voxel_MR(3))+1 (handles.position_cor(2)/handles.voxel_MR(2))+1 (handles.position_cor(3)/handles.voxel_MR(3))+1 (handles.position_cor(4)/handles.voxel_MR(2))+1]);
    
    f1 = handles.position_sag_cor(1);
    f2 = handles.position_sag_cor(1) + handles.position_sag_cor(3)-1;
    f3 = handles.position_sag_cor(2);
    f4 = handles.position_sag_cor(2) + handles.position_sag_cor(4)-1;
    f5 = handles.position_axi_cor(1);
    f6 = handles.position_axi_cor(1) + handles.position_axi_cor(3)-1;
    
    
    handles.position_sag = [0 0 handles.position_sag(3) handles.position_sag(4)];
    handles.position_axi = [0 0 handles.position_axi(3) handles.position_axi(4)];
    handles.position_cor = [0 0 handles.position_cor(3) handles.position_cor(4)];

    handles.IPCMRA      = handles.IPCMRA(f3:f4,f1:f2,f5:f6);
    handles.Lrgb        = handles.Lrgb(f3:f4,f1:f2,f5:f6,:);
    handles.SEG         = handles.SEG(f3:f4,f1:f2,f5:f6);
    handles.L           = handles.L(f3:f4,f1:f2,f5:f6);
    handles.MR_FFE_FH   = handles.MR_FFE_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_AP   = handles.MR_FFE_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_RL   = handles.MR_FFE_RL(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_FH   = handles.MR_PCA_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_AP   = handles.MR_PCA_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_RL   = handles.MR_PCA_RL(f3:f4,f1:f2,f5:f6,:);
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.ANG(f3:f4,f1:f2,f5:f6);
    handles.MAG = handles.MAG(f3:f4,f1:f2,f5:f6);

    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;

    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    handles.id_resizing = 0;
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(1))-handles.voxel_MR(1); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(2))-handles.voxel_MR(2); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    set(handles.pushbutton6,'visible','off')
    set(handles.pushbutton7,'visible','off')
    set(handles.pushbutton8,'visible','off')
    set(handles.pushbutton9,'visible','off')
    set(handles.pushbutton10,'visible','off')
    set(handles.pushbutton11,'visible','off')
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton3_Callback(hObject, eventdata, handles)
    input.id            = 2;
    input.idangmag      = 1;
    input.IPCMRA        = handles.IPCMRA;
    input.SEG           = handles.SEG;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.NUM           = handles.NUM;
    input.min_value_th  = handles.min_value_th;
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
    input.view_sac      = 2;
    input.id_seg        = handles.id_seg;
    input.id_vel        = handles.id_vel;
    input.id_vor        = handles.id_vor;
    input.id_hd         = handles.id_hd;
    input.id_rhd        = handles.id_rhd;
    input.id_vd         = handles.id_vd;
    input.id_el         = handles.id_el;
    input.id_ke         = handles.id_ke;
    input.Lrgb_vel      = handles.Lrgb_vel;
    input.Lrgb_vor      = handles.Lrgb_vor;
    input.Lrgb_hd       = handles.Lrgb_hd;
    input.Lrgb_rhd      = handles.Lrgb_rhd;
    input.Lrgb_vd       = handles.Lrgb_vd;
    input.Lrgb_el       = handles.Lrgb_el;
    input.Lrgb_ke       = handles.Lrgb_ke;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input.id_fve        = handles.id_fve;
    input.Lrgb_fve      = handles.Lrgb_fve;
    input.id_bve        = handles.id_bve;
    input.Lrgb_bve      = handles.Lrgb_bve;
    input.id_aan        = handles.id_aan;
    input.Lrgb_aan      = handles.Lrgb_aan;
    input.id_fov        = handles.id_fov;
    input.Lrgb_fov      = handles.Lrgb_fov;
    input.peak_flow      = handles.peak_flow;
    GUIDE_SEGMENTATION(input)
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton7_Callback(hObject, eventdata, handles)
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    h = imrect(gca,handles.pos2);
    setColor(h,'c');
    position = wait(h);
    
    %%% Julio Sotelo 23-11-2018
    handles.position_axi = position;
    handles.position_sag = [handles.position_sag(1) position(2) handles.position_sag(3) position(4)];
    handles.position_cor = [position(1) handles.position_cor(2) position(3) handles.position_cor(4)];

    handles.c1 = position(2);
    handles.c2 = handles.pos1(1);
    handles.c3 = position(1);
    handles.c4 = position(4);
    handles.c5 = handles.pos1(3);
    handles.c6 = position(3);
    handles.pos1 = [handles.c2 handles.c1 handles.c5 handles.c4];
    handles.pos2 = [handles.c3 handles.c1 handles.c6 handles.c4];
    handles.pos3 = [handles.c3 handles.c2 handles.c6 handles.c5];
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(1))-handles.voxel_MR(1); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(2))-handles.voxel_MR(2); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton10_Callback(hObject, eventdata, handles)
    handles.position_sag_cor = round([(handles.position_sag(1)/handles.voxel_MR(2))+1 (handles.position_sag(2)/handles.voxel_MR(1))+1 (handles.position_sag(3)/handles.voxel_MR(2))+1 (handles.position_sag(4)/handles.voxel_MR(1))+1]);
    handles.position_axi_cor = round([(handles.position_axi(1)/handles.voxel_MR(3))+1 (handles.position_axi(2)/handles.voxel_MR(1))+1 (handles.position_axi(3)/handles.voxel_MR(3))+1 (handles.position_axi(4)/handles.voxel_MR(1))+1]);
    handles.position_cor_cor = round([(handles.position_cor(1)/handles.voxel_MR(3))+1 (handles.position_cor(2)/handles.voxel_MR(2))+1 (handles.position_cor(3)/handles.voxel_MR(3))+1 (handles.position_cor(4)/handles.voxel_MR(2))+1]);
    
    f1 = handles.position_sag_cor(1);
    f2 = handles.position_sag_cor(1) + handles.position_sag_cor(3)-1;
    f3 = handles.position_sag_cor(2);
    f4 = handles.position_sag_cor(2) + handles.position_sag_cor(4)-1;
    f5 = handles.position_axi_cor(1);
    f6 = handles.position_axi_cor(1) + handles.position_axi_cor(3)-1;
    
    
    handles.position_sag = [0 0 handles.position_sag(3) handles.position_sag(4)];
    handles.position_axi = [0 0 handles.position_axi(3) handles.position_axi(4)];
    handles.position_cor = [0 0 handles.position_cor(3) handles.position_cor(4)];
    
    
    handles.IPCMRA      = handles.IPCMRA(f3:f4,f1:f2,f5:f6);
    handles.Lrgb        = handles.Lrgb(f3:f4,f1:f2,f5:f6,:);
    handles.SEG         = handles.SEG(f3:f4,f1:f2,f5:f6);
    handles.L           = handles.L(f3:f4,f1:f2,f5:f6);
    handles.MR_FFE_FH   = handles.MR_FFE_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_AP   = handles.MR_FFE_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_RL   = handles.MR_FFE_RL(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_FH   = handles.MR_PCA_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_AP   = handles.MR_PCA_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_RL   = handles.MR_PCA_RL(f3:f4,f1:f2,f5:f6,:);
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.ANG(f3:f4,f1:f2,f5:f6);
    handles.MAG = handles.MAG(f3:f4,f1:f2,f5:f6);
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    handles.id_resizing = 0;
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(1))-handles.voxel_MR(1); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(2))-handles.voxel_MR(2); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    set(handles.pushbutton6,'visible','off')
    set(handles.pushbutton7,'visible','off')
    set(handles.pushbutton8,'visible','off')
    set(handles.pushbutton9,'visible','off')
    set(handles.pushbutton10,'visible','off')
    set(handles.pushbutton11,'visible','off')
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton4_Callback(hObject, eventdata, handles)
    input.id            = 2;
    input.idangmag      = 1;
    input.IPCMRA        = handles.IPCMRA;
    input.SEG           = handles.SEG;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.NUM           = handles.NUM;
    input.min_value_th  = handles.min_value_th;
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
    input.view_sac      = 3;
    input.id_seg        = handles.id_seg;
    input.id_vel        = handles.id_vel;
    input.id_vor        = handles.id_vor;
    input.id_hd         = handles.id_hd;
    input.id_rhd        = handles.id_rhd;
    input.id_vd         = handles.id_vd;
    input.id_el         = handles.id_el;
    input.id_ke         = handles.id_ke;
    input.Lrgb_vel      = handles.Lrgb_vel;
    input.Lrgb_vor      = handles.Lrgb_vor;
    input.Lrgb_hd       = handles.Lrgb_hd;
    input.Lrgb_rhd      = handles.Lrgb_rhd;
    input.Lrgb_vd       = handles.Lrgb_vd;
    input.Lrgb_el       = handles.Lrgb_el;
    input.Lrgb_ke       = handles.Lrgb_ke;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input.id_fve        = handles.id_fve;
    input.Lrgb_fve      = handles.Lrgb_fve;
    input.id_bve        = handles.id_bve;
    input.Lrgb_bve      = handles.Lrgb_bve;
    input.id_aan        = handles.id_aan;
    input.Lrgb_aan      = handles.Lrgb_aan;
    input.id_fov        = handles.id_fov;
    input.Lrgb_fov      = handles.Lrgb_fov;
    input.peak_flow      = handles.peak_flow;
    GUIDE_SEGMENTATION(input)
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton8_Callback(hObject, eventdata, handles)
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    h = imrect(gca,handles.pos3);
    setColor(h,'c');
    position = wait(h);
    
    %%% Julio Sotelo 23-11-2018
    handles.position_cor = position;
    handles.position_axi = [position(1) handles.position_axi(2) position(3) handles.position_axi(4)];
    handles.position_sag = [position(2) handles.position_sag(2) position(4) handles.position_sag(4)];

    handles.c1 = handles.pos2(2);
    handles.c2 = position(2);
    handles.c3 = position(1);
    handles.c4 = handles.pos2(4);
    handles.c5 = position(4);
    handles.c6 = position(3);
    handles.pos1 = [handles.c2 handles.c1 handles.c5 handles.c4];
    handles.pos2 = [handles.c3 handles.c1 handles.c6 handles.c4];
    handles.pos3 = [handles.c3 handles.c2 handles.c6 handles.c5];
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(1))-handles.voxel_MR(1); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(2))-handles.voxel_MR(2); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton11_Callback(hObject, eventdata, handles)
    handles.position_sag_cor = round([(handles.position_sag(1)/handles.voxel_MR(2))+1 (handles.position_sag(2)/handles.voxel_MR(1))+1 (handles.position_sag(3)/handles.voxel_MR(2))+1 (handles.position_sag(4)/handles.voxel_MR(1))+1]);
    handles.position_axi_cor = round([(handles.position_axi(1)/handles.voxel_MR(3))+1 (handles.position_axi(2)/handles.voxel_MR(1))+1 (handles.position_axi(3)/handles.voxel_MR(3))+1 (handles.position_axi(4)/handles.voxel_MR(1))+1]);
    handles.position_cor_cor = round([(handles.position_cor(1)/handles.voxel_MR(3))+1 (handles.position_cor(2)/handles.voxel_MR(2))+1 (handles.position_cor(3)/handles.voxel_MR(3))+1 (handles.position_cor(4)/handles.voxel_MR(2))+1]);
    
    f1 = handles.position_sag_cor(1);
    f2 = handles.position_sag_cor(1) + handles.position_sag_cor(3)-1;
    f3 = handles.position_sag_cor(2);
    f4 = handles.position_sag_cor(2) + handles.position_sag_cor(4)-1;
    f5 = handles.position_axi_cor(1);
    f6 = handles.position_axi_cor(1) + handles.position_axi_cor(3)-1;

    
    handles.position_sag = [0 0 handles.position_sag(3) handles.position_sag(4)];
    handles.position_axi = [0 0 handles.position_axi(3) handles.position_axi(4)];
    handles.position_cor = [0 0 handles.position_cor(3) handles.position_cor(4)];

    handles.IPCMRA      = handles.IPCMRA(f3:f4,f1:f2,f5:f6);
    handles.Lrgb        = handles.Lrgb(f3:f4,f1:f2,f5:f6,:);
    handles.SEG         = handles.SEG(f3:f4,f1:f2,f5:f6);
    handles.L           = handles.L(f3:f4,f1:f2,f5:f6);
    handles.MR_FFE_FH   = handles.MR_FFE_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_AP   = handles.MR_FFE_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_FFE_RL   = handles.MR_FFE_RL(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_FH   = handles.MR_PCA_FH(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_AP   = handles.MR_PCA_AP(f3:f4,f1:f2,f5:f6,:);
    handles.MR_PCA_RL   = handles.MR_PCA_RL(f3:f4,f1:f2,f5:f6,:);
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.ANG(f3:f4,f1:f2,f5:f6);
    handles.MAG = handles.MAG(f3:f4,f1:f2,f5:f6);
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    handles.id_resizing = 0;
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    set(handles.pushbutton6,'visible','off')
    set(handles.pushbutton7,'visible','off')
    set(handles.pushbutton8,'visible','off')
    set(handles.pushbutton9,'visible','off')
    set(handles.pushbutton10,'visible','off')
    set(handles.pushbutton11,'visible','off')
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton12_Callback(hObject, eventdata, handles)
    
    input               = [];
    input.SEG           = handles.SEG;
    input.IPCMRA        = handles.IPCMRA;
    input.voxel_MR      = handles.voxel_MR;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.NUM           = handles.NUM;
    input.MR_PCA_FH     = handles.MR_PCA_FH;
    input.MR_PCA_AP     = handles.MR_PCA_AP;
    input.MR_PCA_RL     = handles.MR_PCA_RL;
    input.xd            = handles.xd;
    input.yd            = handles.yd;
    input.zd            = handles.zd;
    input.a             = handles.a;
    input.b             = handles.b;
    input.c             = handles.c;
    input.d             = handles.d;
    input.VENC          = handles.VENC;
    input.heart_rate    = handles.heart_rate;
        
    GUIDE_FE_MESH(input)
    
    handles.faces = getappdata(0,'faces');
    handles.nodes = getappdata(0,'nodes');
    handles.elem = getappdata(0,'elem');
    handles.veset = getappdata(0,'veset');
    handles.cmap = getappdata(0,'cmap');
    handles.mags_vol = getappdata(0,'mags_vol');
    handles.Lrgb_vel = getappdata(0,'Lrgb_vel');
    handles.peak_flow = getappdata(0,'peak_flow');

    
    FH = getappdata(0,'MR_PCA_FH');
    AP = getappdata(0,'MR_PCA_AP');
    RL = getappdata(0,'MR_PCA_RL');
    
    FH_ori = getappdata(0,'MR_PCA_FH_ori');
    AP_ori = getappdata(0,'MR_PCA_AP_ori');
    RL_ori = getappdata(0,'MR_PCA_RL_ori');
    
    if isempty(handles.veset)==0
        handles.mags_vel = handles.mags_vol;
        handles.Lrgb_vel = handles.Lrgb_vel(3:end-2,3:end-2,3:end-2,:,:);
        handles.id_mesh = 1;
        handles.id_vel = 1;
        handles.save_id_mesh_mat = 1;
        handles.save_id_vel_mat = 1;
        handles.save_id_mesh_vtu = 1;
        handles.save_id_vel_vtu = 1;
        handles.save_id_vel_csv=1;
        
        handles.MR_PCA_FH_smooth = FH(3:end-2,3:end-2,3:end-2,:);
        handles.MR_PCA_AP_smooth = AP(3:end-2,3:end-2,3:end-2,:);
        handles.MR_PCA_RL_smooth = RL(3:end-2,3:end-2,3:end-2,:);
        
        handles.MR_PCA_FH = FH_ori(3:end-2,3:end-2,3:end-2,:);
        handles.MR_PCA_AP = AP_ori(3:end-2,3:end-2,3:end-2,:);
        handles.MR_PCA_RL = RL_ori(3:end-2,3:end-2,3:end-2,:);
    

        list_string = {'...','Surface','Voxel','Mesh','Velocity'};
        set(handles.popupmenu1,'visible','on','String',list_string);
        set(handles.pushbutton81,'visible','on');
        
        if isfolder('vWERP')
            set(handles.pushbutton83,'visible','on');
        end
    else
        msgbox('The tetrahedral mesh has not been generated ...','Warning','warn')
    end
    
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider1_Callback(hObject, eventdata, handles)
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
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider2_Callback(hObject, eventdata, handles)
    pp=1/handles.b;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider2,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.b;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.slider_axes2 = handles.slider_value;
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider3_Callback(hObject, eventdata, handles)
    pp=1/handles.a;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider3,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.a;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.slider_axes3 = handles.slider_value;
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton13_Callback(hObject, eventdata, handles)
    input.elem = handles.elem;
    input.faces = handles.faces;
    input.nodes = handles.nodes/1000;
    input.veset = handles.veset;    
    input.veset(unique(sort(handles.faces(:))),:,:) = handles.veset(unique(sort(handles.faces(:))),:,:).*0;    
    input.heart_rate = handles.heart_rate;
    input.mu = 4.5e-3;
    input.density = 1060;
    
    [out] = FAST_FEM(input); % Change
    
    handles.Ten = out.Ten;
    handles.normal = out.normal;
    handles.IWSS = out.IWSS;
    handles.WSS = out.WSS;
    handles.mags_wss = out.mags_wss;
    handles.OSI = out.OSI;
    handles.mags_osi = out.mags_osi;
    handles.VOR = out.VOR;
    handles.mags_vor = out.mags_vor;
    handles.HD = out.HD;    
    handles.mags_hd = out.mags_hd;
    handles.RHD = out.RHD;
    handles.mags_rhd = out.mags_rhd;
    handles.nodevol = out.nodevol;
    handles.VD = out.VD;
    handles.mags_vd = out.mags_vd;
    handles.EL = out.EL;
    handles.mags_el = out.mags_el;
    handles.KE = out.KE;
    handles.mags_ke = out.mags_ke;
    h = out.h;

   
    handles.VORc = handles.VOR;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % eliminate higher values of:
    % vorticity, viscouss dissipation, EL. 
    
    
    [conn,~,~] = meshconn(handles.elem,size(handles.nodes,1));
    sid = unique(sort(handles.faces(:))); 
    handles.list_n = [];
    for n=1:length(sid)
        handles.list_n = [handles.list_n,conn{sid(n)}];
    end
    handles.list_n = sort(unique(handles.list_n'));
    
    handles.VOR(handles.list_n,:,:) = handles.VOR(handles.list_n,:,:).*0;
    handles.VD(handles.list_n,:) = handles.VD(handles.list_n,:).*0;
    handles.EL(handles.list_n,:) = handles.EL(handles.list_n,:).*0;
    
    handles.mags_vor(handles.list_n,:) = handles.mags_vor(handles.list_n,:).*0;
    handles.mags_vd(handles.list_n,:) = handles.mags_vd(handles.list_n,:).*0;
    handles.mags_el(handles.list_n,:) = handles.mags_el(handles.list_n,:).*0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.Lrgb_vor = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    handles.Lrgb_hd  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    handles.Lrgb_rhd = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    handles.Lrgb_vd  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    handles.Lrgb_el  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    handles.Lrgb_ke  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
    
    MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[2,1,3]);
    MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
    xd_seg = MASK.*handles.xd;
    yd_seg = MASK.*handles.yd;
    zd_seg = MASK.*handles.zd;
    xd_seg(MASK==0) = [];
    yd_seg(MASK==0) = [];
    zd_seg(MASK==0) = [];
    pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
    [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,3));
    X_seg = MASK.*X;
    Y_seg = MASK.*Y;
    Z_seg = MASK.*Z;
    X_seg(MASK==0) = [];
    Y_seg(MASK==0) = [];
    Z_seg(MASK==0) = [];
    
    pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];
    [~, ~, ind] = histcounts(handles.mags_vor(:), size(cool(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, cool(128)) * 255);
    rgb_vor = cmap_vol/255;
    [~, ~, ind] = histcounts(handles.mags_hd(:), size(cool(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, cool(128)) * 255);
    rgb_hd = cmap_vol/255;
    [~, ~, ind] = histcounts(handles.mags_rhd(:), size(summer(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, summer(128)) * 255);
    rgb_rhd = cmap_vol/255;
    [~, ~, ind] = histcounts(handles.mags_vd(:), size(winter(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, winter(128)) * 255);
    rgb_vd = cmap_vol/255;
    [~, ~, ind] = histcounts(handles.mags_el(:), size(autumn(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, autumn(128)) * 255);
    rgb_el = cmap_vol/255;
    [~, ~, ind] = histcounts(handles.mags_ke(:), size(spring(128), 1));
    ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
    cmap_vol = uint8(ind2rgb(ind, spring(128)) * 255);
    rgb_ke = cmap_vol/255;
    
    for n=1:length(pos_voxel(:,1))
        d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));
        handles.Lrgb_vor(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_vor(d<mean(handles.voxel_MR)*2,:,:),1);
        handles.Lrgb_hd(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_hd(d<mean(handles.voxel_MR)*2,:,:),1);
        handles.Lrgb_rhd(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_rhd(d<mean(handles.voxel_MR)*2,:,:),1);
        handles.Lrgb_vd(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_vd(d<mean(handles.voxel_MR)*2,:,:),1);
        handles.Lrgb_el(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_el(d<mean(handles.voxel_MR)*2,:,:),1);
        handles.Lrgb_ke(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_ke(d<mean(handles.voxel_MR)*2,:,:),1);
    end
    
    handles.Lrgb_vor = handles.Lrgb_vor(3:end-2,3:end-2,3:end-2,:,:);
    handles.Lrgb_hd = handles.Lrgb_hd(3:end-2,3:end-2,3:end-2,:,:);
    handles.Lrgb_rhd = handles.Lrgb_rhd(3:end-2,3:end-2,3:end-2,:,:);
    handles.Lrgb_vd = handles.Lrgb_vd(3:end-2,3:end-2,3:end-2,:,:);
    handles.Lrgb_el = handles.Lrgb_el(3:end-2,3:end-2,3:end-2,:,:);
    handles.Lrgb_ke = handles.Lrgb_ke(3:end-2,3:end-2,3:end-2,:,:);
    
    list_string = {'...','Surface','Voxel','Mesh','Velocity','WSS','OSI',...
                   'Vorticity','Helicity Density','R. Helicity Density',...
                   'Viscous Dissipation','Energy Loss','Kinetic Energy'};
    set(handles.popupmenu1,'visible','on','String',list_string);
    handles.save_id_wss_mat = 1;
    handles.save_id_osi_mat = 1;
    handles.save_id_vor_mat = 1;
    handles.save_id_hd_mat = 1;
    handles.save_id_rhd_mat = 1;
    handles.save_id_vd_mat = 1;
    handles.save_id_el_mat = 1;
    handles.save_id_ke_mat = 1;
    handles.save_id_wss_vtu = 1;
    handles.save_id_osi_vtu = 1;
    handles.save_id_vor_vtu = 1;
    handles.save_id_hd_vtu = 1;
    handles.save_id_rhd_vtu = 1;
    handles.save_id_vd_vtu = 1;
    handles.save_id_el_vtu = 1;
    handles.save_id_ke_vtu = 1;
    close(h)
    disp('The finite element quantification is finish ...')
    
    set(handles.pushbutton62,'visible','on');% include
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider4_Callback(hObject, eventdata, handles)
    pp=1/handles.d;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider4,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.d;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.peak_flow = handles.slider_value;
    
    if handles.id_vel == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
        %%% Forward Velocity
    elseif handles.id_fve == 1 && handles.id_seg == 0
            
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,handles.peak_flow),handles.forward_velocity(1:end,2,handles.peak_flow),handles.forward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_fvel handles.max_fvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
            
    %%% backward Velocity
    elseif handles.id_bve == 1 && handles.id_seg == 0
            
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,handles.peak_flow),handles.backward_velocity(1:end,2,handles.peak_flow),handles.backward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_bvel handles.max_bvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    %%% wss
    elseif handles.id_wss == 1 && handles.id_seg == 1
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wss handles.max_wss]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    %%% WSS_A
    elseif handles.id_wssa == 1 && handles.id_seg == 1
            
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wssa handles.max_wssa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
            
    %%%%%%%% WSS_C
    elseif handles.id_wssc == 1 && handles.id_seg == 1
            
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wssc handles.max_wssc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    %%% vorticity
    elseif handles.id_vor == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'cool');
        [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vor handles.max_vor]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    %%%%%% Axial Angle
    elseif handles.id_aan == 1 && handles.id_seg == 0
        
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'cool');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_axa handles.max_axa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])

            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    %%%%% Helicity Density
    elseif handles.id_hd == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_hd handles.max_hd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    %%% Relative Helicity Density
    elseif handles.id_rhd == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'summer');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_rhd handles.max_rhd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    %%% Viscous Dissipation
    elseif handles.id_vd == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'winter');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vd handles.max_vd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    %%% Energy Loss
    elseif handles.id_el == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'autumn');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_el handles.max_el]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    %%% Kinetic Energy
    elseif handles.id_ke == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'spring');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_ke handles.max_ke]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    % Circulation
    elseif handles.id_cir == 1 && handles.id_seg == 1
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'winter');
        c = colorbar(handles.axes4);
        handles.min_cir = min(handles.circulation(:));
        handles.max_cir = max(handles.circulation(:));
        handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_cir handles.max_cir];
        c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
        c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
        c.Color = [1 1 1];
        c.Location = 'manual';
        c.Position = [0.2 0.1 0.02 0.5];
        c.FontWeight = 'bold';
        c.Label.String = 'Circulation [mm^{2}/s]';
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
        caxis(handles.axes4, [handles.min_cir handles.max_cir]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))]);
        
    % Forward Vortex
    elseif handles.id_fov == 1 && handles.id_seg == 0
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_fov handles.max_fov]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
   elseif handles.id_aci == 1 && handles.id_seg == 1  
        
        axes(handles.axes4);
        plot(0.0)
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.axial_circulation(:,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'autumn');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_acir handles.max_acir]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))]);

    end
        
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider4_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton14_Callback(hObject, eventdata, handles)
    for n=1:size(handles.veset,3)
        slider_step(1) = 1/handles.d;
        slider_step(2) = 0.1;
        set(handles.slider4,'Visible','on','Value', n/handles.d,'sliderstep',slider_step,'max',1,'min',0)
        %%%% Velocity
        if handles.id_vel == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,n),handles.veset(1:end,2,n),handles.veset(1:end,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'jet');
            [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_vel handles.max_vel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'jet');
                [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_vel handles.max_vel]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%% Forward Velocity
        elseif handles.id_fve == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,n),handles.forward_velocity(1:end,2,n),handles.forward_velocity(1:end,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'jet');
            [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_fvel handles.max_fvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,handles.peak_flow),handles.forward_velocity(1:end,2,handles.peak_flow),handles.forward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'jet');
                [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_fvel handles.max_fvel]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
            %%% backward Velocity
        elseif handles.id_bve == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,n),handles.backward_velocity(1:end,2,n),handles.backward_velocity(1:end,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'jet');
            [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_bvel handles.max_bvel]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,handles.peak_flow),handles.backward_velocity(1:end,2,handles.peak_flow),handles.backward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'jet');
                [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_bvel handles.max_bvel]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
            %%%% WSS
        elseif handles.id_wss == 1 && handles.id_seg == 1
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,n),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,n),handles.WSS(:,2,n),handles.WSS(:,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'hot');
            [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_wss handles.max_wss]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'hot');
                [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_wss handles.max_wss]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%%%%%% WSS_A
        elseif handles.id_wssa == 1 && handles.id_seg == 1
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,n)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,n),handles.WSS_A(:,2,n),handles.WSS_A(:,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_wssa handles.max_wssa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'hot');
                [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_wssa handles.max_wssa]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end  
        %%%%%%%% WSS_C
        elseif handles.id_wssc == 1 && handles.id_seg == 1
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,n)),'CDataMapping','Scaled')
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,n),handles.WSS_C(1:end,2,n),handles.WSS_C(1:end,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'hot');
            [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(:,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_wssc handles.max_wssc]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'hot');
                [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_wssc handles.max_wssc]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%%%%%% Vorticity
        elseif handles.id_vor == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
            hold on
            q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,n),handles.VOR(1:end,2,n),handles.VOR(1:end,3,n),5,'Linewidth',1);
            currentColormap = colormap(handles.axes4,'cool');
            [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
            ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
            handles.cmap = uint8(ind2rgb(ind(1:end,n), currentColormap) * 255);
            handles.cmap(:,:,4) = 255;
            handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
            set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
            set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_vor handles.max_vor]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
                hold on
                q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
                currentColormap = colormap(handles.axes4,'cool');
                [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
                handles.cmap(:,:,4) = 255;
                handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
                set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
                set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_vor handles.max_vor]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%%%% Axial Angle
        elseif handles.id_aan == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,n)),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'cool');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_axa handles.max_axa]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'cool');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_axa handles.max_axa]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
        
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        
        %%%%%% Helidicty Density
        elseif handles.id_hd == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'cool');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_hd handles.max_hd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'cool');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_hd handles.max_hd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% relative helicity density
        elseif handles.id_rhd == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'summer');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_rhd handles.max_rhd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'summer');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_rhd handles.max_rhd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% viscouss dissipation
        elseif handles.id_vd == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'winter');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_vd handles.max_vd]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'winter');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_vd handles.max_vd]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% energy loss
        elseif handles.id_el == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'autumn');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_el handles.max_el]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'autumn');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_el handles.max_el]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% kynetic energy
        elseif handles.id_ke == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'spring');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_ke handles.max_ke]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'spring');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_ke handles.max_ke]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% circulation
        elseif handles.id_cir == 1 && handles.id_seg == 1
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'winter');
            c = colorbar(handles.axes4);
            handles.min_cir = min(handles.circulation(:));
            handles.max_cir = max(handles.circulation(:));
            handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
            c.LimitsMode = 'manual';
            c.Limits = [handles.min_cir handles.max_cir];
            c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
            c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
            c.Color = [1 1 1];
            c.Location = 'manual';
            c.Position = [0.2 0.1 0.02 0.5];
            c.FontWeight = 'bold';
            c.Label.String = 'Circulation [mm^{2}/s]';
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
            caxis(handles.axes4, [handles.min_cir handles.max_cir]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'winter');
                c = colorbar(handles.axes4);
                handles.min_cir = min(handles.circulation(:));
                handles.max_cir = max(handles.circulation(:));
                handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
                c.LimitsMode = 'manual';
                c.Limits = [handles.min_cir handles.max_cir];
                c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
                c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
                c.Color = [1 1 1];
                c.Location = 'manual';
                c.Position = [0.2 0.1 0.02 0.5];
                c.FontWeight = 'bold';
                c.Label.String = 'Circulation [mm^{2}/s]';
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
                caxis(handles.axes4, [handles.min_cir handles.max_cir]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        %%%% forward vortex
        elseif handles.id_fov == 1 && handles.id_seg == 0
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,n),'CDataMapping','Scaled')                
            hold on
            colormap(handles.axes4,'cool');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_fov handles.max_fov]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')                
                hold on
                colormap(handles.axes4,'cool');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_fov handles.max_fov]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end 
        %%%% axial circulation
        elseif handles.id_aci == 1 && handles.id_seg == 1
            axes(handles.axes4);
            plot(0.0)
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,n),'CDataMapping','Scaled')
            hold on
            colormap(handles.axes4,'autumn');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_acir handles.max_acir]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(n),' of ',num2str(size(handles.veset,3))])
            pause(0.05);
            if n==size(handles.veset,3)
                axes(handles.axes4);
                plot(0.0)
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,handles.peak_flow),'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'autumn');
                c = colorbar(handles.axes4);
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
                c.Label.String = 'Axial Circulation [mm^{2}/s]';
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
                caxis(handles.axes4, [handles.min_acir handles.max_acir]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
                set(handles.slider4,'Value', handles.peak_flow/handles.d)
            end
        end
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton5_Callback(hObject, eventdata, handles)
    
    input.id            = 3;
    input.ref           = get(handles.popupmenu1,'Value');
    input.IPCMRA        = handles.IPCMRA;
    input.idangmag      = 1;
    input.voxel_MR      = handles.voxel_MR;
    input.L             = handles.L;
    input.Lrgb          = handles.Lrgb;
    input.SEG           = handles.SEG;
    input.peak_flow     = handles.peak_flow;
    input.d             = handles.d;
    
    if get(handles.popupmenu1,'Value')>=4
        input.faces         = handles.faces;
        input.nodes         = handles.nodes;
    else
        input.faces         = [];
        input.nodes         = [];
    end
    
    if isempty(handles.veset)==0
        input.veset         = handles.veset;
        input.mags_vel      = handles.mags_vel;
    else
        input.veset         = [];
        input.mags_vel      = [];
    end
    
    if isempty(handles.mags_wss)==0
        
        input.WSS           = handles.WSS;
        input.VOR           = handles.VOR;
        input.mags_vor      = handles.mags_vor;
        input.mags_wss      = handles.mags_wss;
        input.mags_osi      = handles.mags_osi;
        input.mags_hd       = handles.mags_hd;
        input.mags_rhd      = handles.mags_rhd;
        input.mags_vd       = handles.mags_vd;
        input.mags_el       = handles.mags_el;
        input.mags_ke       = handles.mags_ke;

    else
        input.WSS           = [];
        input.WSS           = [];
        input.VOR           = [];
        input.mags_wss      = [];
        input.mags_vor      = [];
        input.mags_osi      = [];
        input.mags_hd       = [];
        input.mags_rhd      = [];
        input.mags_vd       = [];
        input.mags_el       = [];
        input.mags_ke       = [];

        
        
        
    end
        
    if isempty(handles.Laplace)==0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        input.Laplace = handles.Laplace;                         
        input.centerline = handles.centerline;                      
        input.radius = handles.radius;                        
        input.diameter = handles.diameter;                        
        input.axial_unit_vectors = handles.axial_unit_vectors;             
        input.circumferential_unit_vectors = handles.circumferential_unit_vectors;   
        input.WSS_A = handles.WSS_A;                          
        input.WSS_C = handles.WSS_C;                           
        input.mag_WSS_A = handles.mag_WSS_A;                       
        input.mag_WSS_C = handles.mag_WSS_C;                       
        input.angle_axial_direction = handles.angle_axial_direction;           
        input.forward_velocity = handles.forward_velocity;                
        input.backward_velocity = handles.backward_velocity;               
        input.mag_forward_velocity = handles.mag_forward_velocity;           
        input.mag_backward_velocity = handles.mag_backward_velocity;           
        input.regurgitant_flow = handles.regurgitant_flow;                
        input.centerline_flow = handles.centerline_flow;                 
        input.eccentricity = handles.eccentricity;
        
        input.curvature = handles.curvature; % Julio Sotelo 05062019
        input.ellipticity = handles.ellipticity; % Julio Sotelo 05062019
        input.length_vessel = handles.length_vessel; % Julio Sotelo 05062019
%         input.circulation = handles.circulation; % Julio Sotelo 05062019
        input.forward_vortex = handles.forward_vortex; % Julio Sotelo 05062019
        input.flattening = handles.flattening; % Julio Sotelo 05062019
        input.area = handles.area; % Julio Sotelo 05062019
        input.axial_circulation = handles.axial_circulation; % Julio Sotelo 05062019
        
    else 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        input.Laplace = [];                         
        input.centerline = [];                      
        input.radius = [];                        
        input.diameter = [];                        
        input.axial_unit_vectors = [];           
        input.circumferential_unit_vectors = [];   
        input.WSS_A = [];                         
        input.WSS_C = [];                          
        input.mag_WSS_A = [];                      
        input.mag_WSS_C = [];                     
        input.angle_axial_direction = [];           
        input.forward_velocity = [];               
        input.backward_velocity = [];               
        input.mag_forward_velocity = [];           
        input.mag_backward_velocity = [];           
        input.regurgitant_flow = [];               
        input.centerline_flow = [];                
        input.eccentricity = []; 
        input.curvature = []; % Julio Sotelo 05062019
        input.ellipticity = []; % Julio Sotelo 05062019
        input.length_vessel = []; % Julio Sotelo 05062019
%         input.circulation = []; % Julio Sotelo 05062019
        input.forward_vortex = []; % Julio Sotelo 05062019
        input.flattening = []; % Julio Sotelo 05062019
        input.area = []; % Julio Sotelo 05062019
        input.axial_circulation = []; % Julio Sotelo 05062019
        
    end
    GUIDE_SEGMENTATION(input)
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu1_Callback(hObject, eventdata, handles)
switch get(handles.popupmenu1,'Value')
    case 1
        axes(handles.axes4);
        plot(0,0)
        axis off
        set(handles.pushbutton5,'visible','off')
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    case 2
        set(handles.pushbutton5,'visible','on')
        [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
        xd = X*handles.voxel_MR(1);
        yd = Y*handles.voxel_MR(2);
        zd = Z*handles.voxel_MR(3);
        xd = permute(xd,[2 1 3]);
        yd = permute(yd,[2 1 3]);
        zd = permute(zd,[2 1 3]);
        id_L = unique(sort(handles.L(:)));
        id_L(id_L==0)=[];
        axes(handles.axes4);
        plot(0,0)
        axis off
        axes(handles.axes4);
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
        set(handles.pushbutton14,'visible','off')
        set(handles.slider4,'visible','off')
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    case 3
        set(handles.pushbutton5,'visible','on')
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
        axes(handles.axes4);
        plot(0,0)
        axis off
        axes(handles.axes4);
        patch('Vertices', node_bin_aorta, 'Faces', elem_bin_aorta,'FaceColor',[0.85 0.85 0.85],'EdgeColor','k');
        lighting gouraud
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.pushbutton14,'visible','off')
        set(handles.slider4,'Visible','off')
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    case 4
        set(handles.pushbutton5,'visible','on')
        [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
        xd = X*handles.voxel_MR(1);
        yd = Y*handles.voxel_MR(2);
        zd = Z*handles.voxel_MR(3);
        xd = permute(xd,[2 1 3]);
        yd = permute(yd,[2 1 3]);
        zd = permute(zd,[2 1 3]);
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'Vertices',handles.nodes,'FaceColor','r','EdgeColor','k');
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        set(handles.pushbutton14,'visible','off')
        set(handles.slider4,'Visible','off')
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    case 5
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.veset(1:end,1,handles.peak_flow),handles.veset(1:end,2,handles.peak_flow),handles.veset(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mags_vel(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vel handles.max_vel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vel(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vel(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vel(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 1;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 6
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_wss(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS(:,1,handles.peak_flow),handles.WSS(:,2,handles.peak_flow),handles.WSS(:,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mags_wss(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wss handles.max_wss]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 1;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 7
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_osi,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_osi handles.max_osi]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 1;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 8
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor','k','facecolor',[0.5 0.5 0.5],'facealpha',0.1)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.VOR(1:end,1,handles.peak_flow),handles.VOR(1:end,2,handles.peak_flow),handles.VOR(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'cool');
        [~, ~, ind] = histcounts(handles.mags_vor(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vor handles.max_vor]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vor(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vor(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vor(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 1;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 9
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_hd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_hd handles.max_hd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_hd(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_hd(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_hd(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 1;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 10
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_rhd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'summer');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_rhd handles.max_rhd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_rhd(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_rhd(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_rhd(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 1;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 11
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_vd(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'winter');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_vd handles.max_vd]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vd(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vd(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_vd(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 1;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 12
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_el(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'autumn');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_el handles.max_el]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
	    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
	    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
	    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
	    axes(handles.axes1);
	    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
	    hold on
	    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
	    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_el(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes1,'gray')
	    axis off
	    daspect([1 1 1])
	    axes(handles.axes2);
	    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
	    hold on
	    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
	    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_el(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes2,'gray')
	    axis off
	    daspect([1 1 1])
	    axes(handles.axes3);
	    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
	    hold on
	    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
	    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_el(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes3,'gray')
	    axis off
	    daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 1;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    case 13
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.mags_ke(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'spring');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_ke handles.max_ke]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
	    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
	    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
	    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
	    axes(handles.axes1);
	    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
	    hold on
	    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
	    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_ke(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes1,'gray')
	    axis off
	    daspect([1 1 1])
	    axes(handles.axes2);
	    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
	    hold on
	    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
	    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_ke(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes2,'gray')
	    axis off
	    daspect([1 1 1])
	    axes(handles.axes3);
	    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
	    hold on
	    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
	    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_ke(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
	    hold off
	    axis image
	    colormap(handles.axes3,'gray')
	    axis off
	    daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 1;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 14 % Laplace
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_lap handles.max_lap]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 1;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
        
    case 15 % centerline
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        plot3(handles.centerline(:,1),handles.centerline(:,2),handles.centerline(:,3),'-c','LineWidth',3)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 1;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 16 % diameter
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.diameter,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_dia handles.max_dia]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 1;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 17 % radius
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.radius,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_rad handles.max_rad]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 1;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 18 % axial unit vector 
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.axial_unit_vectors(:,1),handles.axial_unit_vectors(:,2),handles.axial_unit_vectors(:,3),1,'g','LineWidth',1)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 1;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 19 % circumferential unit vector  
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.circumferential_unit_vectors(:,1),handles.circumferential_unit_vectors(:,2),handles.circumferential_unit_vectors(:,3),1,'r','LineWidth',1)
        hold off
        axis vis3d
        lighting gouraud
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 1;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 20 % WSSA   
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_A(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(:,1),handles.nodes(:,2),handles.nodes(:,3),handles.WSS_A(:,1,handles.peak_flow),handles.WSS_A(:,2,handles.peak_flow),handles.WSS_A(:,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_A(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wssa handles.max_wssa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 1;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 21 % WSSC
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.mag_WSS_C(:,1,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.WSS_C(1:end,1,handles.peak_flow),handles.WSS_C(1:end,2,handles.peak_flow),handles.WSS_C(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'hot');
        [~, ~, ind] = histcounts(handles.mag_WSS_C(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
        handles.cmap = uint8(ind2rgb(ind(:,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_wssc handles.max_wssc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 1;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 22 % Axial Angle
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',squeeze(handles.angle_axial_direction(:,:,handles.peak_flow)),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_axa handles.max_axa]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_aan(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_aan(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_aan(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 1;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 23 % Forward Velocity
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.forward_velocity(1:end,1,handles.peak_flow),handles.forward_velocity(1:end,2,handles.peak_flow),handles.forward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.forward_velocity,1) size(handles.forward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_fvel handles.max_fvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fve(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fve(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fve(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 1;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 24 % Backward Velocity    
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'edgecolor',[0.5 0.5 0.5],'facecolor','none','edgealpha',0.5)
        hold on
        q = quiver3(handles.nodes(1:end,1),handles.nodes(1:end,2),handles.nodes(1:end,3),handles.backward_velocity(1:end,1,handles.peak_flow),handles.backward_velocity(1:end,2,handles.peak_flow),handles.backward_velocity(1:end,3,handles.peak_flow),5,'Linewidth',1);
        currentColormap = colormap(handles.axes4,'jet');
        [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(currentColormap, 1));
        ind = reshape(ind,[size(handles.backward_velocity,1) size(handles.backward_velocity,3)]);
        handles.cmap = uint8(ind2rgb(ind(1:end,handles.peak_flow), currentColormap) * 255);
        handles.cmap(:,:,4) = 255;
        handles.cmap = permute(repmat(handles.cmap, [1 3 1]), [2 1 3]);
        set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:3,:,:), [], 4).');
        set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', reshape(handles.cmap(1:2,:,:), [], 4).');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_bvel handles.max_bvel]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_bve(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_bve(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_bve(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 1;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 25 % Regurgitant Flow
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.regurgitant_flow,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'summer');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_rfl handles.max_rfl]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 1;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 26 % Centerline Flow   
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
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
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 1;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 27 % Eccentricity
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.eccentricity,'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'winter');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_ecc handles.max_ecc]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 1;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

   case 28 % Curvature
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        
        if min(handles.curvature(:)) == 0 && max(handles.curvature(:)) ==0 % julio
            
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none')
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
            
        else
            
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.curvature,'CDataMapping','Scaled')
            colormap(handles.axes4,'hot');
            c = colorbar(handles.axes4);
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
            caxis(handles.axes4, [handles.min_cur handles.max_cur]);
            hold off
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        end

        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 1; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

    case 29 % Ellipticity
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.ellipticity,'CDataMapping','Scaled')
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_ell handles.max_ell]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 1; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

   case 30 % Length of the Vessel
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.length_vessel,'CDataMapping','Scaled')
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_len handles.max_len]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 1; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
%     case 31 % circulation
%         
%         set(handles.pushbutton5,'visible','on')
%         axes(handles.axes4);
%         plot(0.0)
%         axis off
%         axes(handles.axes4);
%         patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.circulation(:,handles.peak_flow),'CDataMapping','Scaled')
%         hold on
%         colormap(handles.axes4,'winter');
%         c = colorbar(handles.axes4);
%         handles.min_cir = min(handles.circulation(:));
%         handles.max_cir = max(handles.circulation(:));
%         handles.mean_cir = (handles.min_cir + handles.max_cir)/2;
%         c.LimitsMode = 'manual';
%         c.Limits = [handles.min_cir handles.max_cir];
%         c.Ticks = [handles.min_cir, (handles.min_cir + handles.mean_cir)/2, handles.mean_cir, (handles.max_cir + handles.mean_cir)/2, handles.max_cir];
%         c.TickLabels = {num2str(handles.min_cir,'%0.2f'), num2str((handles.min_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.mean_cir,'%0.2f'), num2str((handles.max_cir + handles.mean_cir)/2,'%0.2f'), num2str(handles.max_cir,'%0.2f')};
%         c.Color = [1 1 1];
%         c.Location = 'manual';
%         c.Position = [0.2 0.1 0.02 0.5];
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
%         c.Label.Color = [1 1 1];
%         caxis(handles.axes4, [handles.min_cir handles.max_cir]);
%         hold off
%         axis vis3d
%         daspect([1,1,1])
%         axis off
%         view([-34,-51])
%         
% 	    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
%         if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
%         if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
%         axes(handles.axes1);
%         imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
%         hold on
%         plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
%         plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
%         Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
%         himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%         cdata = double(cdata)*0.5;
%         set(himage, 'AlphaData', cdata);
%         hold off
%         axis image
%         colormap(handles.axes1,'gray')
%         axis off
%         daspect([1 1 1])
%         axes(handles.axes2);
%         imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
%         hold on
%         plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
%         plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
%         Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
%         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
%         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%         cdata = double(cdata)*0.5;
%         set(himage, 'AlphaData', cdata);
%         hold off
%         axis image
%         colormap(handles.axes2,'gray')
%         axis off
%         daspect([1 1 1])
%         axes(handles.axes3);
%         imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
%         hold on
%         plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
%         plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
%         Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
%         himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
%         cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
%         cdata = double(cdata)*0.5;
%         set(himage, 'AlphaData', cdata);
%         hold off
%         axis image
%         colormap(handles.axes3,'gray')
%         axis off
%         daspect([1 1 1])
%         handles.id_seg = 1;
%         handles.id_vel = 0;
%         handles.id_wss = 0;
%         handles.id_osi = 0;
%         handles.id_vor = 0;
%         handles.id_hd = 0;
%         handles.id_rhd = 0;
%         handles.id_vd = 0;
%         handles.id_el = 0;
%         handles.id_ke = 0;
%         handles.id_lap  = 0;
%         handles.id_cen  = 0;
%         handles.id_rad  = 0;
%         handles.id_dia  = 0;
%         handles.id_auv  = 0;
%         handles.id_cuv  = 0;
%         handles.id_wssa = 0;
%         handles.id_wssc = 0;
%         handles.id_aan  = 0;
%         handles.id_fve  = 0;
%         handles.id_bve  = 0;
%         handles.id_ref  = 0;
%         handles.id_cenf = 0;
%         handles.id_ecc  = 0;
%         handles.id_cur  = 0; % Julio Sotelo 28-05-2019
%         handles.id_ell  = 0; % Julio Sotelo 28-05-2019
%         handles.id_len  = 0; % Julio Sotelo 28-05-2019
% %         handles.id_cir  = 1; % Julio Sotelo 28-05-2019
%         handles.id_fov  = 0; % Julio Sotelo 28-05-2019
%         handles.id_fla  = 0; % Julio Sotelo 28-05-2019
%         handles.id_are  = 0; % Julio Sotelo 28-05-2019
%         handles.id_aci  = 0; % Julio Sotelo 28-05-2019
%         slider_step(1) = 1/size(handles.veset,3);
%         slider_step(2) = 0.1;
%         set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
%         set(handles.slider4,'visible','on')
%         set(handles.pushbutton14,'visible','on')
%         set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    case 31 % forward vortex
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.forward_vortex(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'cool');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_fov handles.max_fov]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fov(:,:,handles.slider_axes1,handles.peak_flow,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fov(:,handles.slider_axes2,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb_fov(handles.slider_axes3,:,:,handles.peak_flow,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 0;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 1; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
    case 32 % Flattening
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.flattening,'CDataMapping','Scaled')
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_fla handles.max_fla]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 1; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
        
     case 33 % Area
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.area,'CDataMapping','Scaled')
        colormap(handles.axes4,'parula');
        c = colorbar(handles.axes4);
        handles.min_are = min(handles.area(:));
        handles.max_are = max(handles.area(:));
        handles.mean_are = (handles.min_are + handles.max_are)/2;
        c.LimitsMode = 'manual';
        c.Limits = [handles.min_are handles.max_are];
        c.Ticks = [handles.min_are, (handles.min_are + handles.mean_are)/2, handles.mean_are, (handles.max_are + handles.mean_are)/2, handles.max_are];
        c.TickLabels = {num2str(handles.min_are,'%0.2f'), num2str((handles.min_are + handles.mean_are)/2,'%0.2f'), num2str(handles.mean_are,'%0.2f'), num2str((handles.max_are + handles.mean_are)/2,'%0.2f'), num2str(handles.max_are,'%0.2f')};
        c.Color = [1 1 1];
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
        caxis(handles.axes4, [handles.min_are handles.max_are]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
        
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 1; % Julio Sotelo 28-05-2019
        handles.id_aci  = 0; % Julio Sotelo 28-05-2019
        set(handles.slider4,'visible','off')
        set(handles.pushbutton14,'visible','off')
        set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
    
    case 34 % axial circulation    
        
        set(handles.pushbutton5,'visible','on')
        axes(handles.axes4);
        plot(0.0)
        axis off
        axes(handles.axes4);
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.axial_circulation(:,handles.peak_flow),'CDataMapping','Scaled')
        hold on
        colormap(handles.axes4,'autumn');
        c = colorbar(handles.axes4);
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
        caxis(handles.axes4, [handles.min_acir handles.max_acir]);
        hold off
        axis vis3d
        daspect([1,1,1])
        axis off
        view([-34,-51])
        
	    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
        cdata = double(cdata)*0.5;
        set(himage, 'AlphaData', cdata);
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        handles.id_seg = 1;
        handles.id_vel = 0;
        handles.id_wss = 0;
        handles.id_osi = 0;
        handles.id_vor = 0;
        handles.id_hd = 0;
        handles.id_rhd = 0;
        handles.id_vd = 0;
        handles.id_el = 0;
        handles.id_ke = 0;
        handles.id_lap  = 0;
        handles.id_cen  = 0;
        handles.id_rad  = 0;
        handles.id_dia  = 0;
        handles.id_auv  = 0;
        handles.id_cuv  = 0;
        handles.id_wssa = 0;
        handles.id_wssc = 0;
        handles.id_aan  = 0;
        handles.id_fve  = 0;
        handles.id_bve  = 0;
        handles.id_ref  = 0;
        handles.id_cenf = 0;
        handles.id_ecc  = 0;
        handles.id_cur  = 0; % Julio Sotelo 28-05-2019
        handles.id_ell  = 0; % Julio Sotelo 28-05-2019
        handles.id_len  = 0; % Julio Sotelo 28-05-2019
%         handles.id_cir  = 0; % Julio Sotelo 28-05-2019
        handles.id_fov  = 0; % Julio Sotelo 28-05-2019
        handles.id_fla  = 0; % Julio Sotelo 28-05-2019
        handles.id_are  = 0; % Julio Sotelo 28-05-2019
        handles.id_aci  = 1; % Julio Sotelo 28-05-2019
        slider_step(1) = 1/size(handles.veset,3);
        slider_step(2) = 0.1;
        set(handles.slider4,'Value', handles.peak_flow/size(handles.veset,3),'sliderstep',slider_step,'max',1,'min',0)
        set(handles.slider4,'visible','on')
        set(handles.pushbutton14,'visible','on')
        set(handles.text4,'Visible','on','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Folder_Callback(hObject, eventdata, handles)

    
    % THIS FUNCTION LOAD THE INFORMATION TO BE PROCESSED 
    handles.id_unwrappping = 0;
    set(handles.uipanel1,'Visible','off')
    cla(handles.uipanel1,'reset');
    path(path,['IO_CODES',filesep]) % cambiar
    path(path,'iso2mesh/')
    folder_name = uigetdir([],'Load Folder...');
    
    if folder_name==0
        return;  
    end
    
    % reading folder and subfolders
    allSubFolders = genpath(folder_name);
    remain = allSubFolders;
    
    if ismac || isunix 
      
      listOfFolderNames = {};
      while true
        [singleSubFolder, remain] = strtok(remain, ':');
        if isempty(singleSubFolder)
          break;
        end
        listOfFolderNames = [listOfFolderNames singleSubFolder];
      end
      numberOfFolders = length(listOfFolderNames);
      
    elseif ispc
      listOfFolderNames = {};
      while true
        [singleSubFolder, remain] = strtok(remain, ';');
        if isempty(singleSubFolder)
          break;
        end
        listOfFolderNames = [listOfFolderNames singleSubFolder];
      end
      numberOfFolders = length(listOfFolderNames);
    end
   
    % reading files names matlab, par-rec, dcm, dat
    files_names_mat = [];
    files_names_par = [];
    files_names_rec = [];
    files_names_dcm = [];
    
    cont_mat = 1;
    cont_par = 1;
    cont_rec = 1;
    cont_dcm = 1;
    
    cont_dat_vd_1 = 1;
    cont_dat_vd_2 = 1;
    cont_dat_vd_3 = 1;
    
    files_names_dat_vd_1 = [];
    files_names_dat_vd_2 = [];
    files_names_dat_vd_3 = [];
    files_names_txt_hd = [];
    files_names_dat_CD = [];
    
    % reading matlab files
    for k = 1 : numberOfFolders
        
        thisFolder = listOfFolderNames{k};
        filePattern_mat = sprintf('%s/*.mat', thisFolder);
        baseFileNames_mat = dir(filePattern_mat);
        numberOfFiles_mat = length(baseFileNames_mat);

        if numberOfFiles_mat>=1
            for n = 1:numberOfFiles_mat
                files_names_mat{cont_mat} = [thisFolder,filesep,baseFileNames_mat(n).name];
                cont_mat = cont_mat + 1;
            end
        end
        
    end
    
    % reading other files
    for k = 1 : numberOfFolders
        
        if isempty(files_names_mat)==1
            
            
            thisFolder = listOfFolderNames{k};
            
            filePattern_par = sprintf('%s/*.PAR', thisFolder);
            baseFileNames_par = dir(filePattern_par);
            numberOfFiles_par = length(baseFileNames_par);
            
            filePattern_rec = sprintf('%s/*.REC', thisFolder);
            baseFileNames_rec = dir(filePattern_rec);
            numberOfFiles_rec = length(baseFileNames_rec);
            
            filePattern_dcm = sprintf('%s/*.dcm', thisFolder);
            baseFileNames_dcm = dir(filePattern_dcm);
            numberOfFiles_dcm = length(baseFileNames_dcm);
            
            filePattern_dat = sprintf('%s/*.dat', thisFolder);
            baseFileNames_dat = dir(filePattern_dat);
            numberOfFiles_dat = length(baseFileNames_dat);

            if numberOfFiles_par>=1 
                for n = 1:numberOfFiles_par
                    files_names_par{cont_par} = [thisFolder,filesep,baseFileNames_par(n).name];
                    cont_par = cont_par + 1;
                end
            end

            if numberOfFiles_rec>=1 
                for n = 1:numberOfFiles_rec
                    files_names_rec{cont_rec} = [thisFolder,filesep,baseFileNames_rec(n).name];
                    cont_rec = cont_rec + 1;
                end
            end

            if numberOfFiles_dcm>=1 
                for n = 1:numberOfFiles_dcm
                    files_names_dcm{cont_dcm} = [thisFolder,filesep,baseFileNames_dcm(n).name];
                    cont_dcm = cont_dcm + 1;
                end
            end

            if numberOfFiles_dat>=1 
                files_names_txt_hd = [thisFolder,filesep,'pcvipr_header.txt'];
                files_names_dat_CD = [thisFolder,filesep,'CD.dat'];
                for n = 1:numberOfFiles_dat
                    if strncmp(flip(baseFileNames_dat(n).name),'tad.1_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
                        files_names_dat_vd_1{cont_dat_vd_1} = [thisFolder,filesep,baseFileNames_dat(n).name];
                        cont_dat_vd_1 = cont_dat_vd_1 + 1;
                    elseif  strncmp(flip(baseFileNames_dat(n).name),'tad.2_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
                        files_names_dat_vd_2{cont_dat_vd_2} = [thisFolder,filesep,baseFileNames_dat(n).name];
                        cont_dat_vd_2 = cont_dat_vd_2 + 1;
                    elseif  strncmp(flip(baseFileNames_dat(n).name),'tad.3_dv_',9)==1 && strncmp(baseFileNames_dat(n).name,'ph_',3)==1
                        files_names_dat_vd_3{cont_dat_vd_3} = [thisFolder,filesep,baseFileNames_dat(n).name];
                        cont_dat_vd_3 = cont_dat_vd_3 + 1;
                    end
                end
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % READING MATLAB FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(files_names_mat)==0
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove hidden files
        mat_temp = [];
        cont=1;
        for ff = 1:length(files_names_mat)
            if isempty(strfind(files_names_mat{ff},'._'))==1
                mat_temp{cont}=files_names_mat{ff};
                cont = cont + 1;
            end
        end
        files_names_mat = mat_temp;
        
        
        % loading structure data
        h = msgbox({'Please wait ...','Loading structure data ...'});
        
        fullFileName = files_names_mat{1};
        load(fullFileName);
        
        handles.VENC = data.VENC;
        handles.voxel_MR = data.voxel_MR;
        handles.heart_rate = data.heart_rate;
        handles.type = data.type;
        handles.MR_FFE_FH = data.MR_FFE_FH;
        handles.MR_FFE_AP = data.MR_FFE_AP;
        handles.MR_FFE_RL = data.MR_FFE_RL;
        handles.MR_PCA_FH = data.MR_PCA_FH;
        handles.MR_PCA_AP = data.MR_PCA_AP;
        handles.MR_PCA_RL = data.MR_PCA_RL;
        close(h)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data offset error JSOTELO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.type = 'MAT';
        input.VENC = handles.VENC;
        input.voxel_MR = handles.voxel_MR;
        input.heart_rate = handles.heart_rate;
        input.type = handles.type;
        input.MR_FFE_FH = handles.MR_FFE_FH;
        input.MR_FFE_AP = handles.MR_FFE_AP;
        input.MR_FFE_RL = handles.MR_FFE_RL;
        input.MR_PCA_FH = handles.MR_PCA_FH;
        input.MR_PCA_AP = handles.MR_PCA_AP;
        input.MR_PCA_RL = handles.MR_PCA_RL;
        input.id = 0;
        input.id_while = 0;        
        id_while = 0;
        while(1)
            while(id_while == 0)
                OFFSET_ERR_AND_NOISE_MASKING(input)
                input.VENC = getappdata(0,'VENC');
                input.voxel_MR = getappdata(0,'voxel_MR');
                input.heart_rate = getappdata(0,'heart_rate');
                input.type = getappdata(0,'type');
                input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
                input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
                input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
                input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
                input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
                input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
                input.id = getappdata(0,'id');
                id_while = getappdata(0,'id_while');
            end
            handles.VENC = getappdata(0,'VENC');
            handles.voxel_MR = getappdata(0,'voxel_MR');
            handles.heart_rate = getappdata(0,'heart_rate');
            handles.type = getappdata(0,'type');
            handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
            handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
            handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
            handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
            handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
            handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
            break
        end
        % load data offset error JSOTELO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        MR_FFE_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_FH;
        MR_FFE_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_AP;
        MR_FFE_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_RL;
        MR_PCA_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_FH;
        MR_PCA_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_AP;
        MR_PCA_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_RL;
        handles.MR_FFE_FH   = MR_FFE_FH_n;
        handles.MR_FFE_AP   = MR_FFE_AP_n;
        handles.MR_FFE_RL   = MR_FFE_RL_n;
        handles.MR_PCA_FH   = MR_PCA_FH_n;
        handles.MR_PCA_AP   = MR_PCA_AP_n;
        handles.MR_PCA_RL   = MR_PCA_RL_n;
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        handles.IPCMRA = (1/d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
        handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
        
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % READING PAR-REC FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(files_names_par)==0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove hidden files
        par_temp = [];
        rec_temp = [];
        cont=1;
        for ff = 1:length(files_names_par)
            if isempty(strfind(files_names_par{ff},'._'))==1
                par_temp{cont}=files_names_par{ff};
                cont = cont + 1;
            end
        end
        cont=1;
        for ff = 1:length(files_names_rec)
            if isempty(strfind(files_names_rec{ff},'._'))==1
                rec_temp{cont}=files_names_rec{ff};
                cont = cont + 1;
            end
        end
        
        files_names_par = par_temp;
        files_names_rec = rec_temp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        version = 0;
        h = waitbar(0,['Reading ',num2str(0),' files of ',num2str(length(files_names_rec)),' ...']);
        steps = round(length(files_names_rec)/3);
        st = steps;
        steps = 1;
        st = steps;
        handles.MR_FFE_FH = [];
        handles.MR_FFE_RL = [];
        handles.MR_FFE_AP = [];
        handles.MR_PCA_FH = [];
        handles.MR_PCA_RL = [];
        handles.MR_PCA_AP = [];
        for f = 1:length(files_names_par)
            fullFileName = files_names_par{f};
            fid_r1 = fopen(fullFileName);
            tline = fgetl(fid_r1);
            cont = 1;
            M = [];
            ID_M_1 = 0;
            ID_M_2 = 0;
            ID_M_3 = 0;
            while ischar(tline)
                  tline = fgetl(fid_r1);
                  if isempty(strfind(tline,'.    Acquisition nr                     :'))==0
                    p = textscan(tline,'%s','delimiter', ':');
                    g = p{1};
                    acquisition_number = str2double(g(2));
                  end
                  if isempty(strfind(tline,'.    Reconstruction nr                  :'))==0
                    p = textscan(tline,'%s','delimiter', ':');
                    g = p{1};
                    reconstruction_number = str2double(g(2));
                  end
                  if isempty(strfind(tline,'.    Phase encoding velocity [cm/sec]   :'))==0
                    p = textscan(tline,'%s','delimiter', ':');
                    g = textscan(p{1}{2},'%f');
                    phase_encoding = g{1};
                    [r,c,v] = find(phase_encoding>0);
                    if r == 1; id = 'RL'; end
                    if r == 2; id = 'AP'; end
                    if r == 3; id = 'FH'; end
                  end
                  if strcmp(tline,'#  sl ec  dyn ph ty    idx pix scan% rec size                (re)scale              window        angulation              offcentre        thick   gap   info      spacing     echo     dtime   ttime    diff  avg  flip    freq   RR-int  turbo delay b grad cont anis         diffusion       L.ty')==1  || ID_M_1==1
                      ID_M_1 = 1;
                      version = 1;
                      if cont>=100
                          if isempty(tline)==1
                              break
                          end
                          g = textscan(tline,'%f');
                          M = [M;g{1}'];
                      end
                  end
                  if strcmp(tline,'#sl ec dyn ph ty  idx pix % rec size (re)scale     window       angulation      offcentre         thick  gap   info   spacing   echo  dtime ttime    diff avg  flip  freq RR_int  turbo  delay b grad cont anis diffusion')==1  || ID_M_2==1
                      ID_M_2 = 1;
                      version = 2;
                      if cont>=100
                          if isempty(tline)==1
                              break
                          end
                          g = textscan(tline,'%f');
                          M = [M;g{1}'];
                      end
                  end
                  if strcmp(tline,'#  sl ec  dyn ph ty    idx pix scan% rec size                (re)scale              window        angulation              offcentre        thick   gap   info      spacing     echo     dtime   ttime    diff  avg  flip    freq   RR-int  turbo delay b grad cont anis         diffusion       L.ty  contagent   controute  contvolume  conttime  contdose  contingr  contingrconcen')==1  || ID_M_3==1
                      ID_M_3 = 1;
                      version = 3;
                      if cont>=107
                          if isempty(tline)==1
                              break
                          end
                          g = textscan(tline,'%f');
                          M = [M;g{1}'];
                      end
                  end
                  cont = cont + 1;
            end
            slice_number = M(:,1);
            cardiac_phase_number = M(:,4);
            Image_type = M(:,5);
            Reconstruction_resolution = M(:,10:11);
            Rescale_intercept = M(:,12);
            Rescale_slope = M(:,13);
            Scale_slope = M(:,14);
            if version == 1 || version == 2 || version == 3 
                Slice_gap = M(:,23) + M(:,24); 
            end
            Pixel_spacing = M(:,29:30);
            if version == 1 || version == 3
                Cardiac_frequency = abs(M(:,37));
            elseif version == 2
                Cardiac_frequency = 60;
            end
            fid2 = fopen(files_names_rec{f});
            [data] = fread(fid2,'int16');
            fclose(fid2);
            if version == 1 || version == 3 
                IM = reshape(data, Reconstruction_resolution(1,1), Reconstruction_resolution(1,2), length(slice_number));
                for n=1:length(slice_number)
                    if Image_type(n)==0
                        name_file = ['handles.MR_FFE_',id];
                    else
                        name_file = ['handles.MR_PCA_',id];
                    end
                    eval([name_file,'(:,:,slice_number(n),cardiac_phase_number(n))= double(IM(:,:,n))*Rescale_slope(n)+Rescale_intercept(n);'])
                end
            end
            if version ==2
                IM = reshape(data, Reconstruction_resolution(1,1), Reconstruction_resolution(1,2), length(slice_number));
                for n=1:length(slice_number)
                    if Image_type(n)==0
                        name_file = ['handles.MR_FFE_',id];
                        eval([name_file,'(:,:,slice_number(n),cardiac_phase_number(n))= double(IM(:,:,n))*Rescale_slope(n)+Rescale_intercept(n);'])
                    else
                        name_file = ['handles.MR_PCA_',id];
                        eval([name_file,'(:,:,slice_number(n),cardiac_phase_number(n))= ((double(IM(:,:,n))*Rescale_slope(n)+Rescale_intercept(n))/abs(Rescale_intercept(n)))*max(abs(phase_encoding));'])
                    end
                end
            end
            if st==f
                waitbar(st/ length(files_names_par),h,['Reading ',num2str(st),' files of ',num2str(length(files_names_par)),' ...']);
                st = st+steps;
            end
        end
        close(h)
        eval('handles.voxel_MR = [Pixel_spacing(1,:),Slice_gap(1)];')
        eval('handles.heart_rate = Cardiac_frequency(1,1);')
        handles.VENC = max(abs(phase_encoding));
        handles.MR_FFE_FH = double(handles.MR_FFE_FH);
        handles.MR_FFE_AP = double(handles.MR_FFE_AP);
        handles.MR_FFE_RL = double(handles.MR_FFE_RL);
        handles.MR_PCA_FH = double(handles.MR_PCA_FH);
        handles.MR_PCA_AP = double(handles.MR_PCA_AP);
        handles.MR_PCA_RL = double(handles.MR_PCA_RL);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data offset error JSOTELO
        handles.type = 'PAR-REC';
        input.VENC = handles.VENC;
        input.voxel_MR = handles.voxel_MR;
        input.heart_rate = handles.heart_rate;
        input.type = handles.type;
        input.MR_FFE_FH = handles.MR_FFE_FH;
        input.MR_FFE_AP = handles.MR_FFE_AP;
        input.MR_FFE_RL = handles.MR_FFE_RL;
        input.MR_PCA_FH = handles.MR_PCA_FH;
        input.MR_PCA_AP = handles.MR_PCA_AP;
        input.MR_PCA_RL = handles.MR_PCA_RL;
        input.id = 0;
        input.id_while = 0;        
        id_while = 0;
        while(1)
            while(id_while == 0)
                OFFSET_ERR_AND_NOISE_MASKING(input)
                input.VENC = getappdata(0,'VENC');
                input.voxel_MR = getappdata(0,'voxel_MR');
                input.heart_rate = getappdata(0,'heart_rate');
                input.type = getappdata(0,'type');
                input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
                input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
                input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
                input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
                input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
                input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
                input.id = getappdata(0,'id');
                id_while = getappdata(0,'id_while');
            end
            handles.VENC = getappdata(0,'VENC');
            handles.voxel_MR = getappdata(0,'voxel_MR');
            handles.heart_rate = getappdata(0,'heart_rate');
            handles.type = getappdata(0,'type');
            handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
            handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
            handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
            handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
            handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
            handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
            break
        end
        % load data offset error JSOTELO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        id_empty = double([isempty(handles.MR_FFE_FH),isempty(handles.MR_FFE_AP),isempty(handles.MR_FFE_RL)]);
        if sum(id_empty)>0
            [r,c,v] = find(id_empty==0);
            if c(1)==1
                handles.MR_FFE_FH = handles.MR_FFE_FH;
                handles.MR_FFE_AP = handles.MR_FFE_FH;
                handles.MR_FFE_RL = handles.MR_FFE_FH;
            elseif c(1)==2
                handles.MR_FFE_FH = handles.MR_FFE_AP;
                handles.MR_FFE_AP = handles.MR_FFE_AP;
                handles.MR_FFE_RL = handles.MR_FFE_AP;
            elseif c(1)==3
                handles.MR_FFE_FH = handles.MR_FFE_RL;
                handles.MR_FFE_AP = handles.MR_FFE_RL;
                handles.MR_FFE_RL = handles.MR_FFE_RL;
                    
            end
        end
        
        MASK = double(handles.MR_PCA_FH~=max(phase_encoding)*-1);
        handles.MR_PCA_FH = handles.MR_PCA_FH.*MASK;
        handles.MR_PCA_AP = handles.MR_PCA_AP.*MASK;
        handles.MR_PCA_RL = handles.MR_PCA_RL.*MASK;
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        MR_FFE_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_PCA_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
        MR_FFE_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_FH;
        MR_FFE_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_AP;
        MR_FFE_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_RL;
        MR_PCA_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_FH;
        MR_PCA_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_AP;
        MR_PCA_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_RL;
        handles.MR_FFE_FH   = MR_FFE_FH_n;
        handles.MR_FFE_AP   = MR_FFE_AP_n;
        handles.MR_FFE_RL   = MR_FFE_RL_n;
        handles.MR_PCA_FH   = MR_PCA_FH_n;
        handles.MR_PCA_AP   = MR_PCA_AP_n;
        handles.MR_PCA_RL   = MR_PCA_RL_n;
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        handles.IPCMRA = (1/d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
        handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
        
    elseif isempty(files_names_dcm)==0
        
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % remove hidden files
            dcm_temp = [];
            cont=1;
            for ff = 1:length(files_names_dcm)
                if isempty(strfind(files_names_dcm{ff},'._'))==1
                    dcm_temp{cont}=files_names_dcm{ff};
                    cont = cont + 1;
                end
            end
            files_names_dcm = dcm_temp;
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%             h = msgbox({'Please wait ...','Reading Manufacturer Information'});
%             info_dicom = dicominfo(files_names_dcm{1});
%             close(h)
%             if strcmp(info_dicom.Manufacturer,'SIEMENS')==1
%                 images = zeros(info_dicom.Rows,info_dicom.Columns,length(files_names_dcm));
%                 infoStr = cell(length(files_names_dcm),1);
%                 info = cell(length(files_names_dcm),1);
%                 NumberofCardiacPhases = info_dicom.CardiacNumberOfImages;
%                 NumberofSlices = (length(files_names_dcm)/NumberofCardiacPhases)/4;
%                 NumberofSeries = [];
%                 NumberinSeries = [];
%                 h = waitbar(0,['Reading ',num2str(0),' files of ',num2str(length(files_names_dcm)),' ...']);
%                 steps = round(length(files_names_dcm)/50);
%                 st = steps;
%                 for k=1:length(files_names_dcm)
%                     [infoStruct, dicInfoStruct, msg] = dicom_scan_singlefile(files_names_dcm{k});
%                     infoStr{k} = infoStruct;
%                     info{k} = dicInfoStruct;
%                     NumberofSeries(k) = info{k}.SeriesNumber;
%                     NumberinSeries(k) = infoStr{k}.image;
%                     images(:,:,k) = dicom_read_singlefile(files_names_dcm{k});
%                     if st==k 
%                         waitbar(st/ length(files_names_dcm),h,['Reading ',num2str(st),' files of ',num2str(length(files_names_dcm)),' ...']);
%                         st = st+steps;
%                     end
%                 end
%                 serie1 = min(NumberofSeries);
%                 serie2 = max(NumberofSeries);
%                 NumberofImagesSerie1 = sum(NumberofSeries==serie1);
%                 NumberofImagesSerie2 = sum(NumberofSeries==serie2);
%                 if NumberofImagesSerie1>NumberofImagesSerie2
%                     PhaseImagesID = serie1;
%                     MagnitudeImagesID = serie2;
%                 else
%                     PhaseImagesID = serie2;
%                     MagnitudeImagesID = serie1;
%                 end
%                 close(h)          
%                 voxelsize = infoStr{1}.voxVc(1:3);
%                 heart_rate = round(60/(info{1}.NominalInterval/1000));
%                 phaseRange  = 4096;
%                 if isfield(info{1},'ImageOrientationPatient')
%                     read_vector  = info{1}.ImageOrientationPatient(1:3);
%                     phase_vector = info{1}.ImageOrientationPatient(4:6);
%                     slice_vector = cross(read_vector,phase_vector);
%                     mainOrient   = find(abs(slice_vector) == max(abs(slice_vector)));
%                     if mainOrient == 1
%                         orientation = 'sag';
%                     elseif mainOrient == 2
%                         orientation = 'cor';
%                     else
%                         orientation = 'tra';
%                     end
%                 else
%                     errStr =  sprintf('%s\n%s',errStr,'Field ''ImageOrientationPatient'' does not exist in dicomInfoStruct.');
%                 end
%                 if isfield(info{1},'Columns')&& isfield(info{end},'Rows')
%                     imageSize = [info{1}.Rows,info{1}.Columns];
%                 else
%                     errStr =  sprintf('%s\n%s',errStr,'Field ''Columns & Rows'' does not exist in dicomInfoStruct.');
%                 end
%                 if isfield(info{1},'PhaseEncodingDirection')
%                     if strcmp(info{1}.PhaseEncodingDirection,'ROW')
%                         peDir = 'i';
%                     else
%                         peDir = 'j';
%                     end
%                 else
%                     errStr =  sprintf('%s\n%s',errStr,'Field ''PhaseEncodingDirection'' does not exist in dicomInfoStruct.');
%                 end
%                 if ~isempty(strfind(info{1}.SoftwareVersions,'VB13')) || ~isempty(strfind(info{1}.SoftwareVersions,'B15')) || ~isempty(strfind(info{1}.SoftwareVersions,'B17'))
%                    swVersion = 'VB13';
%                 elseif ~isempty(strfind(info{1}.SoftwareVersions,'VA25'))
%                    swVersion = 'VA25';
%                 elseif ~isempty(strfind(info{1}.SoftwareVersions,'VB12'))
%                    swVersion = 'VB12';
%                 else
%                    swVersion = 'VB13';
%                 end
%                 if (~strcmp(swVersion,'VB12') && ~strcmp(swVersion,'VA25'))
%                    if strcmp(orientation,'sag')
%                        if strcmp(peDir,'i')
%                            signVijk(1)= -1;
%                            signVijk(2)= -1;
%                            signVijk(3)= -1;
%                        else
%                            signVijk(1)= -1;
%                            signVijk(2)=  1;
%                            signVijk(3)= -1;
%                        end
%                    elseif strcmp(orientation,'cor')
%                        if strcmp(peDir,'i')
%                            signVijk(1)= -1;
%                            signVijk(2)=  1;
%                            signVijk(3)= -1;
%                        else
%                            signVijk(1)= -1;
%                            signVijk(2)= -1;
%                            signVijk(3)= -1;
%                        end
%                    elseif strcmp(orientation,'tra')
%                        if strcmp(peDir,'i')
%                            signVijk(1)= -1;
%                            signVijk(2)=  1;
%                            signVijk(3)= -1;
%                        else
%                            signVijk(1)=  1;
%                            signVijk(2)=  1;
%                            signVijk(3)= -1;
%                        end
%                    end
%                 else
%                    signVijk(1)= 1;
%                    signVijk(2)= 1;
%                    signVijk(3)= 1;
%                 end
%                 if isfield (info{1},'Private_0029_1020')
%                     tempStr = lower(info{1}.Private_0029_1020);
%                     if strcmp(swVersion,'VB13')
%                         posVencInplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[8]'));
%                         posVencThplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[9]'));
%                         posVencNextField  = strfind(tempStr,lower('sWiPMemBlock.adFree[10]'));
%                     elseif (strcmp(swVersion,'VB12')||strcmp(swVersion,'VA25'))
%                         posVencInplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[9]'));
%                         posVencThplane    = strfind(tempStr,lower('sWiPMemBlock.adFree[10]'));
%                         posVencNextField  = strfind(tempStr,lower('sWiPMemBlock.adFree[11]'));
%                     else
%                         posVencInplane    = [];
%                         posVencThplane    = [];
%                         posVencNextField  = [];
%                     end
%                     if ~isempty(posVencInplane) && ~isempty(posVencThplane) && ~isempty(posVencNextField)
%                         tempStr = tempStr(posVencInplane:posVencNextField-1);
%                         pos_eq = strfind(tempStr,'=');
%                         [~, idx] = regexp(tempStr, '\n', 'match', 'start');
%                         vencInPlane = str2double(tempStr(pos_eq(1)+1:idx(1)-1));
%                         vencThPlane = str2double(tempStr(pos_eq(2)+1:idx(2)-1));
%                     else
%                         vencInPlane = 1.5;
%                         vencThPlane = 1.5;
%                     end
%                 end
%                 if vencInPlane == vencThPlane
%                     VENC = vencInPlane;
%                 else
%                     msgbox('The VENC need to be the same in in-plane and th-plane ...','Warning','warn')
%                 end
%                 images_n = double(images);
%                 for k=1:size(images_n,3)
%                     if NumberofSeries(k)==PhaseImagesID
%                         images_n(:,:,k) = ((images(:,:,k)*info{k}.RescaleSlope + info{k}.RescaleIntercept)/phaseRange)*VENC*100;
%                     else
%                         images_n(:,:,k) = images(:,:,k);
%                     end
%                 end
%                 phase_sort = NumberinSeries(NumberofSeries==PhaseImagesID);
%                 magnitude_sort = NumberinSeries(NumberofSeries==MagnitudeImagesID);
%                 PHASE = images_n(:,:,NumberofSeries==PhaseImagesID);
%                 MAGNITUDE = images_n(:,:,NumberofSeries==MagnitudeImagesID);
%                 PHASE = PHASE(:,:,phase_sort);
%                 MAGNITUDE = MAGNITUDE(:,:,magnitude_sort);
%                 handles.MR_PCA_FH = flip(permute(reshape(PHASE(:,:,1:NumberofSlices*NumberofCardiacPhases),[size(PHASE,1) size(PHASE,2) NumberofSlices NumberofCardiacPhases]),[2,1,3,4])*-1,3);
%                 handles.MR_PCA_AP = flip(permute(reshape(PHASE(:,:,(NumberofSlices*NumberofCardiacPhases)+1:(NumberofSlices*NumberofCardiacPhases)*2),[size(PHASE,1) size(PHASE,2) NumberofSlices NumberofCardiacPhases]),[2,1,3,4])*-1,3);
%                 handles.MR_PCA_RL = flip(permute(reshape(PHASE(:,:,((NumberofSlices*NumberofCardiacPhases)*2)+1:(NumberofSlices*NumberofCardiacPhases)*3),[size(PHASE,1) size(PHASE,2) NumberofSlices NumberofCardiacPhases]),[2,1,3,4]),3);
%                 handles.MR_FFE = flip(permute(reshape(MAGNITUDE,[size(PHASE,1) size(PHASE,2) NumberofSlices NumberofCardiacPhases]),[2,1,3,4]),3);
%                 handles.VENC = VENC*100;
%                 handles.voxel_MR = voxelsize;
%                 handles.heart_rate = heart_rate;
% 
%                 id_empty = double([isempty(handles.MR_FFE_FH),isempty(handles.MR_FFE_AP),isempty(handles.MR_FFE_RL)]);
%                 if sum(id_empty)>0
%                     [r,c,v] = find(id_empty==0);
%                     if c(1)==1
%                         handles.MR_FFE_FH = handles.MR_FFE_FH;
%                         handles.MR_FFE_AP = handles.MR_FFE_FH;
%                         handles.MR_FFE_RL = handles.MR_FFE_FH;
%                     elseif c(1)==2
%                         handles.MR_FFE_FH = handles.MR_FFE_AP;
%                         handles.MR_FFE_AP = handles.MR_FFE_AP;
%                         handles.MR_FFE_RL = handles.MR_FFE_AP;
%                     elseif c(1)==3
%                         handles.MR_FFE_FH = handles.MR_FFE_RL;
%                         handles.MR_FFE_AP = handles.MR_FFE_RL;
%                         handles.MR_FFE_RL = handles.MR_FFE_RL;
% 
%                     end
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % load data offset error JSOTELO
%                 handles.type = 'DCM';
%                 input.VENC = handles.VENC;
%                 input.voxel_MR = handles.voxel_MR;
%                 input.heart_rate = handles.heart_rate;
%                 input.type = handles.type;
%                 input.MR_FFE_FH = handles.MR_FFE_FH;
%                 input.MR_FFE_AP = handles.MR_FFE_AP;
%                 input.MR_FFE_RL = handles.MR_FFE_RL;
%                 input.MR_PCA_FH = handles.MR_PCA_FH;
%                 input.MR_PCA_AP = handles.MR_PCA_AP;
%                 input.MR_PCA_RL = handles.MR_PCA_RL;
%                 input.id = 0;
%                 input.id_while = 0;        
%                 id_while = 0;
%                 while(1)
%                     while(id_while == 0)
%                         OFFSET_ERR_AND_NOISE_MASKING(input)
%                         input.VENC = getappdata(0,'VENC');
%                         input.voxel_MR = getappdata(0,'voxel_MR');
%                         input.heart_rate = getappdata(0,'heart_rate');
%                         input.type = getappdata(0,'type');
%                         input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
%                         input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
%                         input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
%                         input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
%                         input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
%                         input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
%                         input.id = getappdata(0,'id');
%                         id_while = getappdata(0,'id_while');
%                     end
%                     handles.VENC = getappdata(0,'VENC');
%                     handles.voxel_MR = getappdata(0,'voxel_MR');
%                     handles.heart_rate = getappdata(0,'heart_rate');
%                     handles.type = getappdata(0,'type');
%                     handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
%                     handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
%                     handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
%                     handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
%                     handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
%                     handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
%                     break
%                 end
%                 % load data offset error JSOTELO
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 [a,b,c,d] = size(handles.MR_FFE);
%                 handles.a = a;
%                 handles.b = b;
%                 handles.c = c;
%                 handles.d = d;
%                 MR_FFE_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_FFE_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_FFE_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_PCA_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_PCA_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_PCA_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
%                 MR_FFE_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE;
%                 MR_FFE_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE;
%                 MR_FFE_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE;
%                 MR_PCA_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_FH;
%                 MR_PCA_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_AP;
%                 MR_PCA_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_RL;
%                 handles.MR_FFE_FH   = MR_FFE_FH_n;
%                 handles.MR_FFE_AP   = MR_FFE_AP_n;
%                 handles.MR_FFE_RL   = MR_FFE_RL_n;
%                 handles.MR_PCA_FH   = MR_PCA_FH_n;
%                 handles.MR_PCA_AP   = MR_PCA_AP_n;
%                 handles.MR_PCA_RL   = MR_PCA_RL_n;
%                 [a,b,c,d] = size(handles.MR_FFE_FH);
%                 handles.a = a;
%                 handles.b = b;
%                 handles.c = c;
%                 handles.d = d;
%                 handles.IPCMRA = (1/d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
%                 handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
%                 
% % %             elseif strcmp(info_dicom.Manufacturer,'Philips Medical Systems')==1

            if strcmp(info_dicom.Manufacturer,'Philips')==1       
                h = waitbar(0,['Reading ',num2str(0),' files of ',num2str(length(files_names_dcm)),' ...']);
                steps = round(length(files_names_dcm)/50);
                st = steps;
                handles.MR_FFE = [];
                handles.MR_PCA_FH = [];
                handles.MR_PCA_AP = [];
                handles.MR_PCA_RL = [];

                if length(files_names_dcm)>4
                    for k = 1:length(files_names_dcm)
                        info_temp = dicominfo(files_names_dcm{k});
                        RescaleSlope = info_temp.RescaleSlope;
                        RescaleIntercept = info_temp.RescaleIntercept;
%                         RescaleSlope = info_temp.Private_2005_140a;
%                         RescaleIntercept = info_temp.Private_2005_1409;
                        im_temp = double(dicomread(files_names_dcm{k}))*RescaleSlope + RescaleIntercept;
                        [r,c,~] = find(info_temp.Private_2001_101a>0);
                        if r == 1; id = 'RL'; end
                        if r == 2; id = 'AP'; end
                        if r == 3; id = 'FH'; end
                        if strcmp(info_temp.Private_2005_106e,'FFE')==1
                            name_file = 'handles.MR_FFE';
                        else
                            name_file = ['handles.MR_',info_temp.Private_2005_106e,'_',id];
                        end
                        eval([name_file,'(:,:,info_temp.Private_2001_100a,info_temp.Private_2001_1008)= im_temp;'])
                        if st==k;
                            waitbar(st/ length(files_names_dcm),h,['Reading ',num2str(st),' files of ',num2str(length(files_names_dcm)),' ...']);
                            st = st+steps;
                        end
                    end
                    close(h)
                    eval(['handles.voxel_MR = [info_temp.PixelSpacing;info_temp.SpacingBetweenSlices]',char(39),';'])
                    eval('handles.heart_rate = info_temp.HeartRate;')
                else
                    steps = round(length(files_names_dcm)/3);
                    st = steps;
                    for k = 1:length(files_names_dcm)
                        info_temp = dicominfo(files_names_dcm{k});
                        im_temp = squeeze(double(dicomread(files_names_dcm{k})));
                        for f = 1:info_temp.NumberOfFrames
                            PixelSpacing = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.PixelSpacing']);
                            SpacingBetweenSlices = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.SpacingBetweenSlices']);
                            RescaleSlope = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.RescaleSlope']);
                            RescaleIntercept = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.RescaleIntercept']);
                            Image_Id = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.Private_2005_106e']);
                            CPhase = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.Private_2001_1008']);
                            Slice = eval(['info_temp.PerFrameFunctionalGroupsSequence.Item_',num2str(f),'.Private_2005_140f.Item_1.Private_2001_100a']);
                            [r,c,v] = find(info_temp.Private_2001_101a>0);
                            if r == 1; id = 'RL'; end
                            if r == 2; id = 'AP'; end
                            if r == 3; id = 'FH'; end
                            if strcmp(Image_Id,'FFE')==1
                                name_file = 'handles.MR_FFE';
                            else
                                name_file = ['handles.MR_',Image_Id,'_',id];
                            end
                            eval([name_file,'(:,:,Slice,CPhase)= squeeze(im_temp(:,:,f))*RescaleSlope + RescaleIntercept;'])
                        end
                        if st==k;
                            waitbar(st/ length(files_names_dcm),h,['Reading ',num2str(st),' files of ',num2str(length(files_names_dcm)),' ...']);
                            st = st+steps;
                        end
                    end
                    close(h)
                    eval(['handles.voxel_MR = [PixelSpacing;SpacingBetweenSlices]',char(39),';'])
                    eval('handles.heart_rate = info_temp.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.HeartRate;')
                end
                handles.VENC = max(abs(info_temp.Private_2001_101a));
                handles.MR_FFE_FH = permute(handles.MR_FFE,[2,1,3,4]);
                handles.MR_FFE_AP = permute(handles.MR_FFE,[2,1,3,4]);
                handles.MR_FFE_RL = permute(handles.MR_FFE,[2,1,3,4]);
                handles.MR_PCA_FH = permute(handles.MR_PCA_FH,[2,1,3,4]);
                handles.MR_PCA_AP = permute(handles.MR_PCA_AP,[2,1,3,4]);
                handles.MR_PCA_RL = permute(handles.MR_PCA_RL,[2,1,3,4]);
                
                id_empty = double([isempty(handles.MR_FFE_FH),isempty(handles.MR_FFE_AP),isempty(handles.MR_FFE_RL)]);
                if sum(id_empty)>0
                    [r,c,v] = find(id_empty==0);
                    if c(1)==1
                        handles.MR_FFE_FH = handles.MR_FFE_FH;
                        handles.MR_FFE_AP = handles.MR_FFE_FH;
                        handles.MR_FFE_RL = handles.MR_FFE_FH;
                    elseif c(1)==2
                        handles.MR_FFE_FH = handles.MR_FFE_AP;
                        handles.MR_FFE_AP = handles.MR_FFE_AP;
                        handles.MR_FFE_RL = handles.MR_FFE_AP;
                    elseif c(1)==3
                        handles.MR_FFE_FH = handles.MR_FFE_RL;
                        handles.MR_FFE_AP = handles.MR_FFE_RL;
                        handles.MR_FFE_RL = handles.MR_FFE_RL;

                    end
                    
                end
                
                
                MASK = double(handles.MR_PCA_FH~=max(info_temp.Private_2001_101a)*-1);
                handles.MR_PCA_FH = handles.MR_PCA_FH.*MASK;
                handles.MR_PCA_AP = handles.MR_PCA_AP.*MASK;
                handles.MR_PCA_RL = handles.MR_PCA_RL.*MASK;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % load data offset error JSOTELO
                handles.type = 'DCM';
                input.VENC = handles.VENC;
                input.voxel_MR = handles.voxel_MR;
                input.heart_rate = handles.heart_rate;
                input.type = handles.type;
                input.MR_FFE_FH = handles.MR_FFE_FH;
                input.MR_FFE_AP = handles.MR_FFE_AP;
                input.MR_FFE_RL = handles.MR_FFE_RL;
                input.MR_PCA_FH = handles.MR_PCA_FH;
                input.MR_PCA_AP = handles.MR_PCA_AP;
                input.MR_PCA_RL = handles.MR_PCA_RL;
                input.id = 0;
                input.id_while = 0;        
                id_while = 0;
                while(1)
                    while(id_while == 0)
                        OFFSET_ERR_AND_NOISE_MASKING(input)
                        input.VENC = getappdata(0,'VENC');
                        input.voxel_MR = getappdata(0,'voxel_MR');
                        input.heart_rate = getappdata(0,'heart_rate');
                        input.type = getappdata(0,'type');
                        input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
                        input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
                        input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
                        input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
                        input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
                        input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
                        input.id = getappdata(0,'id');
                        id_while = getappdata(0,'id_while');
                    end
                    handles.VENC = getappdata(0,'VENC');
                    handles.voxel_MR = getappdata(0,'voxel_MR');
                    handles.heart_rate = getappdata(0,'heart_rate');
                    handles.type = getappdata(0,'type');
                    handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
                    handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
                    handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
                    handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
                    handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
                    handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
                    break
                end
                % load data offset error JSOTELO
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [a,b,c,d] = size(handles.MR_FFE_FH);
                handles.a = a;
                handles.b = b;
                handles.c = c;
                handles.d = d;
                MR_FFE_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_FFE_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_FFE_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_PCA_FH_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_PCA_AP_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_PCA_RL_n = zeros(handles.a+2,handles.b+2,handles.c+2,handles.d);
                MR_FFE_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_FH;
                MR_FFE_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_AP;
                MR_FFE_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_FFE_RL;
                MR_PCA_FH_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_FH;
                MR_PCA_AP_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_AP;
                MR_PCA_RL_n(2:end-1,2:end-1,2:end-1,:)  = handles.MR_PCA_RL;
                handles.MR_FFE_FH   = MR_FFE_FH_n;
                handles.MR_FFE_AP   = MR_FFE_AP_n;
                handles.MR_FFE_RL   = MR_FFE_RL_n;
                handles.MR_PCA_FH   = MR_PCA_FH_n;
                handles.MR_PCA_AP   = MR_PCA_AP_n;
                handles.MR_PCA_RL   = MR_PCA_RL_n;
                [a,b,c,d] = size(handles.MR_FFE_FH);
                handles.a = a;
                handles.b = b;
                handles.c = c;
                handles.d = d;
                handles.IPCMRA = (1/d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
                handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
                
            end
    elseif isempty(files_names_dat_vd_1)==0
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove hidden files
        dat_temp = [];
        cont=1;
        for ff = 1:length(files_names_dat_vd_1)
            if isempty(strfind(files_names_dat_vd_1{ff},'._'))==1
                dat_temp{cont}=files_names_dat_vd_1{ff};
                cont = cont + 1;
            end
        end
        files_names_dat_vd_1 = dat_temp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        handles.heart_rate = GUIDE_HEART_RATE;
        [names, number] = textread(files_names_txt_hd,'%s %s');
        val_header = str2double(number);
        val_ap = val_header(13:15);
        val_fh = val_header(16:18);
        val_rl = val_header(19:21);
        matrixx = val_header(4);
        matrixy = val_header(5);
        matrixz = val_header(6);
        cardiac_phases = val_header(10);
        h = waitbar(0,['Reading ',num2str(0),' files of ',num2str(cardiac_phases),' ...']);
        steps = round(length(cardiac_phases)/length(cardiac_phases));
        st = steps;
        fid1 = fopen(files_names_dat_CD);
        [data1] = fread(fid1,'int16');
        fclose(fid1);
        CD = double(reshape(data1,[matrixx matrixy matrixz]));
        handles.MR_FFE_AP = repmat(CD,1,1,1,cardiac_phases);
        handles.MR_FFE_FH = repmat(CD,1,1,1,cardiac_phases);
        handles.MR_FFE_RL = repmat(CD,1,1,1,cardiac_phases);
        handles.MR_PCA_AP = zeros(matrixx, matrixy, matrixz, cardiac_phases);
        handles.MR_PCA_FH = zeros(matrixx, matrixy, matrixz, cardiac_phases);
        handles.MR_PCA_RL = zeros(matrixx, matrixy, matrixz, cardiac_phases);
        for f = 1:cardiac_phases
            fid = fopen(files_names_dat_vd_1{f});
            [data] = fread(fid,'int16');
            fclose(fid);
            handles.MR_PCA_RL(:,:,:,f) = reshape(data,[matrixx matrixy matrixz]);
            fid = fopen(files_names_dat_vd_2{f});
            [data] = fread(fid,'int16');
            fclose(fid);
            handles.MR_PCA_AP(:,:,:,f) = reshape(data,[matrixx matrixy matrixz]);
            fid = fopen(files_names_dat_vd_3{f});
            [data] = fread(fid,'int16');
            fclose(fid);
            handles.MR_PCA_FH(:,:,:,f) = reshape(data,[matrixx matrixy matrixz]);
            if st==f;
                waitbar(st/ cardiac_phases,h,['Reading ',num2str(st),' files of ',num2str(cardiac_phases),' ...']);
                st = st+steps;
            end
        end
        close(h)
        handles.voxel_MR = [max(abs(val_ap)),max(abs(val_fh)),max(abs(val_rl))];
        handles.VENC = 200;
        handles.MR_FFE_FH = flip(permute(handles.MR_FFE_FH,[2,3,1,4]),3);
        handles.MR_FFE_AP = flip(permute(handles.MR_FFE_AP,[2,3,1,4]),3);
        handles.MR_FFE_RL = flip(permute(handles.MR_FFE_RL,[2,3,1,4]),3);
        handles.MR_PCA_FH = (flip(permute(handles.MR_PCA_FH,[2,3,1,4]),3))/10;
        handles.MR_PCA_AP = (flip(permute(handles.MR_PCA_AP,[2,3,1,4]),3)*-1)/10;
        handles.MR_PCA_RL = (flip(permute(handles.MR_PCA_RL,[2,3,1,4]),3)*-1)/10;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load data offset error JSOTELO
        handles.type = 'DAT';
        input.VENC = handles.VENC;
        input.voxel_MR = handles.voxel_MR;
        input.heart_rate = handles.heart_rate;
        input.type = handles.type;
        input.MR_FFE_FH = handles.MR_FFE_FH;
        input.MR_FFE_AP = handles.MR_FFE_AP;
        input.MR_FFE_RL = handles.MR_FFE_RL;
        input.MR_PCA_FH = handles.MR_PCA_FH;
        input.MR_PCA_AP = handles.MR_PCA_AP;
        input.MR_PCA_RL = handles.MR_PCA_RL;
        input.id = 0;
        input.id_while = 0;        
        id_while = 0;
        while(1)
            while(id_while == 0)
                OFFSET_ERR_AND_NOISE_MASKING(input)
                input.VENC = getappdata(0,'VENC');
                input.voxel_MR = getappdata(0,'voxel_MR');
                input.heart_rate = getappdata(0,'heart_rate');
                input.type = getappdata(0,'type');
                input.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
                input.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
                input.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
                input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
                input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
                input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
                input.id = getappdata(0,'id');
                id_while = getappdata(0,'id_while');
            end
            handles.VENC = getappdata(0,'VENC');
            handles.voxel_MR = getappdata(0,'voxel_MR');
            handles.heart_rate = getappdata(0,'heart_rate');
            handles.type = getappdata(0,'type');
            handles.MR_FFE_FH = getappdata(0,'MR_FFE_FH');
            handles.MR_FFE_AP = getappdata(0,'MR_FFE_AP');
            handles.MR_FFE_RL = getappdata(0,'MR_FFE_RL');
            handles.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
            handles.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
            handles.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
            break
        end
        % load data offset error JSOTELO
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        [a,b,c,d] = size(handles.MR_FFE_FH);
        handles.a = a;
        handles.b = b;
        handles.c = c;
        handles.d = d;
        handles.IPCMRA = squeeze(handles.MR_FFE_FH(:,:,:,1));
        handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
        
    end
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
    
    set(handles.pushbutton2, 'Units', 'pixels');
    handles.pushbutton2_size = get(handles.pushbutton2, 'Position');
    set(handles.pushbutton2, 'Units', 'normalized');
    idx = mod(min(handles.pushbutton2_size(3:4)),2)>1;
    w = floor(min(handles.pushbutton2_size(3:4)));
    w(idx) = w(idx)+1;
    im = imread('Symbols/P3.png');
    g = double(imresize(im,[w-4 w-4],'method','nearest')>0);
    set(handles.pushbutton2,'CData',g,'visible','on')
    set(handles.pushbutton3,'CData',g,'visible','on')
    set(handles.pushbutton4,'CData',g,'visible','on')
    set(handles.pushbutton5,'CData',g)
    set(handles.uipanel1,'Visible','on')
    set(handles.uipanel6,'Visible','on')
    set(handles.pushbutton1,'visible','on','BackgroundColor',[0 0 0])
    set(handles.pushbutton15,'visible','on','BackgroundColor',[0.2 0.2 0.2])
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
        
    handles.id_mesh = 0;
    handles.id_unwrappping = 0;
    handles.id_vel = 0;
    handles.id_wss = 0;
    handles.id_osi = 0;
    handles.id_vor = 0;
    handles.id_hd   = 0;
    handles.id_rhd  = 0;
    handles.id_vd   = 0;
    handles.id_el   = 0;
    handles.id_ke   = 0;
    handles.id_lap  = 0;
    handles.id_cen  = 0;
    handles.id_rad  = 0;
    handles.id_dia  = 0;
    handles.id_auv  = 0;
    handles.id_cuv  = 0;
    handles.id_wssa = 0;
    handles.id_wssc = 0;
    handles.id_aan  = 0;
    handles.id_fve  = 0;
    handles.id_bve  = 0;
    handles.id_ref  = 0;
    handles.id_cenf = 0;
    handles.id_ecc  = 0;
    handles.id_cur  = 0; % Julio Sotelo 28-05-2019
    handles.id_ell  = 0; % Julio Sotelo 28-05-2019
    handles.id_len  = 0; % Julio Sotelo 28-05-2019
    handles.id_cir  = 0; % Julio Sotelo 28-05-2019
    handles.id_fov  = 0; % Julio Sotelo 28-05-2019
    handles.id_fla  = 0; % Julio Sotelo 28-05-2019
    handles.id_are  = 0; % Julio Sotelo 28-05-2019
    handles.id_aci  = 0; % Julio Sotelo 28-05-2019
    
    handles.Lrgb_vel      = [];
    handles.Lrgb_fve      = [];
    handles.Lrgb_bve      = [];
    handles.Lrgb_aan      = [];
    handles.Lrgb_vor      = [];
    handles.Lrgb_hd       = [];
    handles.Lrgb_rhd      = [];
    handles.Lrgb_vd       = [];
    handles.Lrgb_el       = [];
    handles.Lrgb_ke       = [];
    handles.Lrgb_fov       = [];% Julio Sotelo 28-05-2019
        
    handles.time            = [];
    handles.peak_flow       = [];
    handles.peak_flow_ori   = [];
    handles.flow            = [];
    handles.net_flow        = [];
    handles.max_velocity    = [];
    handles.min_velocity    = [];
    
    handles.veset = [];
    handles.VOR = [];
    handles.WSS = [];
    handles.mags_vel = [];
    handles.mags_wss = [];
    handles.mags_osi = [];
    handles.mags_vor = [];
    handles.mags_hd = [];
    handles.mags_rhd = [];
    handles.mags_vd = [];
    handles.mags_el = [];
    handles.mags_ke = [];
    
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
    
    handles.curvature  = []; % Julio Sotelo 28-05-2019
    handles.ellipticity  = []; % Julio Sotelo 28-05-2019
    handles.length_vessel  = []; % Julio Sotelo 28-05-2019
    handles.circulation  = []; % Julio Sotelo 28-05-2019
    handles.forward_vortex  = []; % Julio Sotelo 28-05-2019
    handles.flattening  = []; % Julio Sotelo 28-05-2019
    handles.area  = []; % Julio Sotelo 28-05-2019
    handles.axial_circulation  = []; % Julio Sotelo 28-05-2019
    
    handles.save_id_SEG_mat = 0;
    handles.save_id_SEG_vti = 0;
    handles.save_id_IPCMRA_mat = 1;
    handles.save_id_IPCMRA_vti = 1;
    handles.save_id_MR_PCA_mat = 1;
    handles.save_id_MR_PCA_vti = 1;
    handles.save_id_MR_FFE_mat = 1;
    handles.save_id_MR_FFE_vti = 1;
    handles.save_id_mesh_mat = 0;
    handles.save_id_vel_mat = 0;
    handles.save_id_wss_mat = 0;
    handles.save_id_osi_mat = 0;
    handles.save_id_vor_mat = 0;
    handles.save_id_hd_mat = 0;
    handles.save_id_rhd_mat = 0;
    handles.save_id_vd_mat = 0;
    handles.save_id_el_mat = 0;
    handles.save_id_ke_mat = 0;
    handles.save_id_mesh_vtu = 0;
    handles.save_id_vel_vtu = 0;
    handles.save_id_wss_vtu = 0;
    handles.save_id_osi_vtu = 0;
    handles.save_id_vor_vtu = 0;
    handles.save_id_hd_vtu = 0;
    handles.save_id_rhd_vtu = 0;
    handles.save_id_vd_vtu = 0;
    handles.save_id_el_vtu = 0;
    handles.save_id_ke_vtu = 0;
    
    handles.save_id_lap_mat     = 0;
    handles.save_id_cen_mat     = 0;
    handles.save_id_rad_mat     = 0;
    handles.save_id_dia_mat     = 0;
    handles.save_id_auv_mat     = 0;
    handles.save_id_cuv_mat     = 0;
    handles.save_id_wssa_mat    = 0;
    handles.save_id_wssc_mat    = 0;
    handles.save_id_aan_mat     = 0;
    handles.save_id_fve_mat     = 0;
    handles.save_id_bve_mat     = 0;
    handles.save_id_ref_mat     = 0;
    handles.save_id_cebf_mat    = 0;
    handles.save_id_ecc_mat     = 0;

    handles.save_id_cur_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_mat     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_mat     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_tim_mat     = 0; % Julio Sotelo 28-05-2019 time
    handles.save_id_flo_mat     = 0; % Julio Sotelo 28-05-2019 flow
    handles.save_id_nfl_mat     = 0; % Julio Sotelo 28-05-2019 net_flow
    handles.save_id_mav_mat     = 0; % Julio Sotelo 28-05-2019 max_velocity
    handles.save_id_miv_mat     = 0; % Julio Sotelo 28-05-2019 min_velocity
    
    
    handles.save_id_lap_vtu     = 0;
    handles.save_id_cen_vtu     = 0;
    handles.save_id_rad_vtu     = 0;
    handles.save_id_dia_vtu     = 0;
    handles.save_id_auv_vtu     = 0;
    handles.save_id_cuv_vtu     = 0;
    handles.save_id_wssa_vtu    = 0;
    handles.save_id_wssc_vtu    = 0;
    handles.save_id_aan_vtu     = 0;
    handles.save_id_fve_vtu     = 0;
    handles.save_id_bve_vtu     = 0;
    handles.save_id_ref_vtu     = 0;
    handles.save_id_cebf_vtu    = 0;
    handles.save_id_ecc_vtu     = 0;
    
    handles.save_id_cur_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_vtu     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_vtu     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_vel_csv   = 0;
    handles.save_id_wss_csv   = 0;
    handles.save_id_osi_csv   = 0;
    handles.save_id_vor_csv   = 0;
    handles.save_id_hd_csv    = 0;
    handles.save_id_rhd_csv   = 0;
    handles.save_id_vd_csv    = 0;
    handles.save_id_el_csv    = 0;
    handles.save_id_ke_csv    = 0;
    handles.save_id_rad_csv   = 0;
    handles.save_id_dia_csv   = 0;
    handles.save_id_wssa_csv  = 0;
    handles.save_id_wssc_csv  = 0;
    handles.save_id_aan_csv   = 0;
    handles.save_id_fve_csv   = 0;
    handles.save_id_bve_csv   = 0;
    handles.save_id_ref_csv   = 0;
    handles.save_id_ecc_csv   = 0;
    
    handles.save_id_cur_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_ell_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_len_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_cir_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fov_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_fla_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_are_csv     = 0;% Julio Sotelo 28-05-2019
    handles.save_id_aci_csv     = 0;% Julio Sotelo 28-05-2019
    
    handles.save_id_tim_csv     = 0; % Julio Sotelo 28-05-2019 time
    handles.save_id_flo_csv     = 0; % Julio Sotelo 28-05-2019 flow
    handles.save_id_nfl_csv     = 0; % Julio Sotelo 28-05-2019 net_flow
    handles.save_id_mav_csv     = 0; % Julio Sotelo 28-05-2019 max_velocity
    handles.save_id_miv_csv     = 0; % Julio Sotelo 28-05-2019 min_velocity
    
    handles.faces = [];
    handles.elem = [];
    handles.nodes = [];
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [r,c,v] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    handles.id_ang = 1;
    handles.id_mag = 0;
    
    handles.id_csv_save = 0;
    handles.SECTIONS_S = [];
    handles.SECTIONS_V = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Showing the logo image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    txt = 'Symbols/VIEWS.png';
    imlogo = imread(txt);
    [av,bv,~]=size(imlogo);
    windows_screen_size = get(0,'ScreenSize');
    imlogo = imresize(imlogo,[round(windows_screen_size(4)*av/(av+bv)) round(windows_screen_size(4)*bv/(av+bv))]);
    Sz= size(imlogo);
    flogo = figure('Position',[10000 10000 Sz(2) + 4 Sz(1) + 4],'name','VIEWS','numbertitle','off','menubar','none');
    movegui(flogo,'center');
    set(flogo,'Units', 'pixels');
    image(imlogo(:,:,1:3))
    set(gca,'Visible','off','Units','pixels','Position', [2 2 Sz(2) Sz(1)]);
    waitfor(flogo)
        

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Data_Callback(hObject, eventdata, handles)
    
    input.save_id_SEG_mat       = handles.save_id_SEG_mat;
    input.save_id_SEG_vti       = handles.save_id_SEG_vti;
    input.save_id_IPCMRA_mat    = handles.save_id_IPCMRA_mat;
    input.save_id_IPCMRA_vti    = handles.save_id_IPCMRA_vti;
    input.save_id_MR_PCA_mat    = handles.save_id_MR_PCA_mat;
    input.save_id_MR_PCA_vti    = handles.save_id_MR_PCA_vti;
    input.save_id_MR_FFE_mat    = handles.save_id_MR_FFE_mat;
    input.save_id_MR_FFE_vti    = handles.save_id_MR_FFE_vti;
    input.save_id_mesh_mat      = handles.save_id_mesh_mat;
    input.save_id_vel_mat       = handles.save_id_vel_mat;
    input.save_id_wss_mat       = handles.save_id_wss_mat;
    input.save_id_osi_mat       = handles.save_id_osi_mat;
    input.save_id_vor_mat       = handles.save_id_vor_mat;
    input.save_id_hd_mat        = handles.save_id_hd_mat;
    input.save_id_rhd_mat       = handles.save_id_rhd_mat;
    input.save_id_vd_mat        = handles.save_id_vd_mat;
    input.save_id_el_mat        = handles.save_id_el_mat;
    input.save_id_ke_mat        = handles.save_id_ke_mat;
    input.save_id_mesh_vtu      = handles.save_id_mesh_vtu;
    input.save_id_vel_vtu       = handles.save_id_vel_vtu;
    input.save_id_wss_vtu       = handles.save_id_wss_vtu;
    input.save_id_osi_vtu       = handles.save_id_osi_vtu;
    input.save_id_vor_vtu       = handles.save_id_vor_vtu;
    input.save_id_hd_vtu        = handles.save_id_hd_vtu;
    input.save_id_rhd_vtu       = handles.save_id_rhd_vtu;
    input.save_id_vd_vtu        = handles.save_id_vd_vtu;
    input.save_id_el_vtu        = handles.save_id_el_vtu;
    input.save_id_ke_vtu        = handles.save_id_ke_vtu;
    
    input.save_id_lap_mat       = handles.save_id_lap_mat;
    input.save_id_cen_mat       = handles.save_id_cen_mat;
    input.save_id_rad_mat       = handles.save_id_rad_mat;
    input.save_id_dia_mat       = handles.save_id_dia_mat;
    input.save_id_auv_mat       = handles.save_id_auv_mat;
    input.save_id_cuv_mat       = handles.save_id_cuv_mat;
    input.save_id_wssa_mat      = handles.save_id_wssa_mat;
    input.save_id_wssc_mat      = handles.save_id_wssc_mat;
    input.save_id_aan_mat       = handles.save_id_aan_mat;
    input.save_id_fve_mat       = handles.save_id_fve_mat;
    input.save_id_bve_mat       = handles.save_id_bve_mat;
    input.save_id_ref_mat       = handles.save_id_ref_mat;
    input.save_id_cebf_mat      = handles.save_id_cebf_mat;
    input.save_id_ecc_mat       = handles.save_id_ecc_mat;

    input.save_id_cur_mat       = handles.save_id_cur_mat; % Julio Sotelo 28-05-2018
    input.save_id_ell_mat       = handles.save_id_ell_mat; % Julio Sotelo 28-05-2018
    input.save_id_len_mat       = handles.save_id_len_mat; % Julio Sotelo 28-05-2018
%     input.save_id_cir_mat       = handles.save_id_cir_mat; % Julio Sotelo 28-05-2018
    input.save_id_fov_mat       = handles.save_id_fov_mat; % Julio Sotelo 28-05-2018
    input.save_id_fla_mat       = handles.save_id_fla_mat; % Julio Sotelo 28-05-2018
    input.save_id_are_mat       = handles.save_id_are_mat; % Julio Sotelo 28-05-2018
    input.save_id_aci_mat       = handles.save_id_aci_mat; % Julio Sotelo 28-05-2018
    
    input.save_id_tim_mat       = handles.save_id_tim_mat; % Julio Sotelo 28-05-2019 time
    input.save_id_flo_mat       = handles.save_id_flo_mat; % Julio Sotelo 28-05-2019 flow
    input.save_id_nfl_mat       = handles.save_id_nfl_mat; % Julio Sotelo 28-05-2019 net_flow
    input.save_id_mav_mat       = handles.save_id_mav_mat; % Julio Sotelo 28-05-2019 max_velocity
    input.save_id_miv_mat       = handles.save_id_miv_mat; % Julio Sotelo 28-05-2019 min_velocity

    input.save_id_lap_vtu       = handles.save_id_lap_vtu;
    input.save_id_cen_vtu       = handles.save_id_cen_vtu;
    input.save_id_rad_vtu       = handles.save_id_rad_vtu;
    input.save_id_dia_vtu       = handles.save_id_dia_vtu;
    input.save_id_auv_vtu       = handles.save_id_auv_vtu;
    input.save_id_cuv_vtu       = handles.save_id_cuv_vtu;
    input.save_id_wssa_vtu      = handles.save_id_wssa_vtu;
    input.save_id_wssc_vtu      = handles.save_id_wssc_vtu;
    input.save_id_aan_vtu       = handles.save_id_aan_vtu;
    input.save_id_fve_vtu       = handles.save_id_fve_vtu;
    input.save_id_bve_vtu       = handles.save_id_bve_vtu;
    input.save_id_ref_vtu       = handles.save_id_ref_vtu;
    input.save_id_cebf_vtu      = handles.save_id_cebf_vtu;
    input.save_id_ecc_vtu       = handles.save_id_ecc_vtu;

    input.save_id_cur_vtu       = handles.save_id_cur_vtu; % Julio Sotelo 28-05-2018
    input.save_id_ell_vtu       = handles.save_id_ell_vtu; % Julio Sotelo 28-05-2018
    input.save_id_len_vtu       = handles.save_id_len_vtu; % Julio Sotelo 28-05-2018
%     input.save_id_cir_vtu       = handles.save_id_cir_vtu; % Julio Sotelo 28-05-2018
    input.save_id_fov_vtu       = handles.save_id_fov_vtu; % Julio Sotelo 28-05-2018
    input.save_id_fla_vtu       = handles.save_id_fla_vtu; % Julio Sotelo 28-05-2018
    input.save_id_are_vtu       = handles.save_id_are_vtu; % Julio Sotelo 28-05-2018
    input.save_id_aci_vtu       = handles.save_id_aci_vtu; % Julio Sotelo 28-05-2018
    
    input.save_id_vel_csv       = 0;
    input.save_id_wss_csv       = 0;
    input.save_id_osi_csv       = 0;
    input.save_id_vor_csv       = 0;
    input.save_id_hd_csv        = 0;
    input.save_id_rhd_csv       = 0;
    input.save_id_vd_csv        = 0;
    input.save_id_el_csv        = 0;
    input.save_id_ke_csv        = 0;
    input.save_id_rad_csv       = 0;
    input.save_id_dia_csv       = 0;
    input.save_id_wssa_csv      = 0;
    input.save_id_wssc_csv      = 0;
    input.save_id_aan_csv       = 0;
    input.save_id_fve_csv       = 0;
    input.save_id_bve_csv       = 0;
    input.save_id_ref_csv       = 0;
    input.save_id_ecc_csv       = 0;
 
    input.save_id_cur_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_ell_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_len_csv       = 0; % Julio Sotelo 28-05-2018
%     input.save_id_cir_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_fov_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_fla_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_are_csv       = 0; % Julio Sotelo 28-05-2018
    input.save_id_aci_csv       = 0; % Julio Sotelo 28-05-2018
    
    input.save_id_tim_csv       = 0; % Julio Sotelo 28-05-2019 time
    input.save_id_flo_csv       = 0; % Julio Sotelo 28-05-2019 flow
    input.save_id_nfl_csv       = 0; % Julio Sotelo 28-05-2019 net_flow
    input.save_id_mav_csv       = 0; % Julio Sotelo 28-05-2019 max_velocity
    input.save_id_miv_csv       = 0; % Julio Sotelo 28-05-2019 min_velocity
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Include message is they want to save excel data
    
    if handles.save_id_lap_mat == 1

        answer1 = questdlg('Would you like save an excel file?', 'Question','Yes','No','No');
        switch answer1
            case 'Yes'
                
                input.save_id_vel_csv   = handles.save_id_vel_mat;
                input.save_id_wss_csv   = handles.save_id_wss_mat;
                input.save_id_osi_csv   = handles.save_id_osi_mat;
                input.save_id_vor_csv   = handles.save_id_vor_mat;
                input.save_id_hd_csv    = handles.save_id_hd_mat;
                input.save_id_rhd_csv   = handles.save_id_rhd_mat;
                input.save_id_vd_csv    = handles.save_id_vd_mat;
                input.save_id_el_csv    = handles.save_id_el_mat;
                input.save_id_ke_csv    = handles.save_id_ke_mat;
                input.save_id_rad_csv   = handles.save_id_rad_mat;
                input.save_id_dia_csv   = handles.save_id_dia_mat;
                input.save_id_wssa_csv  = handles.save_id_wssa_mat;
                input.save_id_wssc_csv  = handles.save_id_wssc_mat;
                input.save_id_aan_csv   = handles.save_id_aan_mat;
                input.save_id_fve_csv   = handles.save_id_fve_mat;
                input.save_id_bve_csv   = handles.save_id_bve_mat;
                input.save_id_ref_csv   = handles.save_id_ref_mat;
                input.save_id_ecc_csv   = handles.save_id_ecc_mat;
                
                input.save_id_cur_csv     = handles.save_id_cur_mat; % Julio Sotelo 28-05-2018
                input.save_id_ell_csv     = handles.save_id_ell_mat; % Julio Sotelo 28-05-2018
                input.save_id_len_csv     = handles.save_id_len_mat; % Julio Sotelo 28-05-2018
%                 input.save_id_cir_csv     = handles.save_id_cir_mat; % Julio Sotelo 28-05-2018
                input.save_id_fov_csv     = handles.save_id_fov_mat; % Julio Sotelo 28-05-2018
                input.save_id_fla_csv     = handles.save_id_fla_mat; % Julio Sotelo 28-05-2018
                input.save_id_are_csv     = handles.save_id_are_mat; % Julio Sotelo 28-05-2018
                input.save_id_aci_csv     = handles.save_id_aci_mat; % Julio Sotelo 28-05-2018
                
                input.save_id_tim_csv       = handles.save_id_tim_mat; % Julio Sotelo 28-05-2019 time
                input.save_id_flo_csv       = handles.save_id_flo_mat; % Julio Sotelo 28-05-2019 flow
                input.save_id_nfl_csv       = handles.save_id_nfl_mat; % Julio Sotelo 28-05-2019 net_flow
                input.save_id_mav_csv       = handles.save_id_mav_mat; % Julio Sotelo 28-05-2019 max_velocity
                input.save_id_miv_csv       = handles.save_id_miv_mat; % Julio Sotelo 28-05-2019 min_velocity
        
                while(1)
                    
                    handles.id_csv_save = 1;

                    prompt = {'Enter the number of sections:'};
                    dlgtitle = 'Input';
                    dims = [1 35];
                    definput = {'1'};
                    answer = inputdlg(prompt,dlgtitle,dims,definput);
                    if isempty(answer)==1
                        input.save_id_vel_csv       = 0;
                        input.save_id_wss_csv       = 0;
                        input.save_id_osi_csv       = 0;
                        input.save_id_vor_csv       = 0;
                        input.save_id_hd_csv        = 0;
                        input.save_id_rhd_csv       = 0;
                        input.save_id_vd_csv        = 0;
                        input.save_id_el_csv        = 0;
                        input.save_id_ke_csv        = 0;
                        input.save_id_rad_csv       = 0;
                        input.save_id_dia_csv       = 0;
                        input.save_id_wssa_csv      = 0;
                        input.save_id_wssc_csv      = 0;
                        input.save_id_aan_csv       = 0;
                        input.save_id_fve_csv       = 0;
                        input.save_id_bve_csv       = 0;
                        input.save_id_ref_csv       = 0;
                        input.save_id_ecc_csv       = 0;

                        input.save_id_cur_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_ell_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_len_csv       = 0; % Julio Sotelo 28-05-2018
%                         input.save_id_cir_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_fov_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_fla_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_are_csv       = 0; % Julio Sotelo 28-05-2018
                        input.save_id_aci_csv       = 0; % Julio Sotelo 28-05-2018
                        
                        input.save_id_tim_csv       = 0; % Julio Sotelo 28-05-2019 time
                        input.save_id_flo_csv       = 0; % Julio Sotelo 28-05-2019 flow
                        input.save_id_nfl_csv       = 0; % Julio Sotelo 28-05-2019 net_flow
                        input.save_id_mav_csv       = 0; % Julio Sotelo 28-05-2019 max_velocity
                        input.save_id_miv_csv       = 0; % Julio Sotelo 28-05-2019 min_velocity
                        
                        waitfor(msgbox('Sections will not be generated ...','Warning','warn'))
                    else
                        handles.number_of_sections = str2num(answer{1});
                        %%% include the sections here %%%%%%%%%%%%%%%
                        % handles.centerline
                        % handles.nodes
                        % handles.faces
                        % handles.elem

                        if handles.number_of_sections > 1

                            dd = cumsum(sqrt(sum((handles.centerline(1:end-1,:)-handles.centerline(2:end,:)).^2,2)));
                            total_length = dd(end);
                            step_sec = total_length/handles.number_of_sections;
                            vect_sec = 1:handles.number_of_sections-1;
                            sec_val = vect_sec.*step_sec;

                            lap_selection = zeros(handles.number_of_sections-1,1);
                            for n=1:length(sec_val)
                                di = abs(dd-sec_val(n));
                                [~,I] = min(di);
                                lap_selection(n)=handles.centerline_lapid(I);
                            end


                            % volume sections
                            ID_V = [1:size(handles.nodes,1)]';
                            handles.SECTIONS_V{1} = ID_V(handles.Laplace>0 & handles.Laplace<=lap_selection(1));
                            for n=2:length(lap_selection)
                                handles.SECTIONS_V{n} = ID_V(handles.Laplace>lap_selection(n-1) & handles.Laplace<=lap_selection(n));
                            end
                            handles.SECTIONS_V{handles.number_of_sections} = ID_V(handles.Laplace>lap_selection(end) & handles.Laplace<1);


                            % surface sections
                            ID_S = sort(unique(handles.faces(:)));
                            Laplace_S = handles.Laplace(ID_S);
                            handles.SECTIONS_S{1} = ID_S(Laplace_S>0 & Laplace_S<=lap_selection(1));
                            for n=2:length(lap_selection)
                                handles.SECTIONS_S{n} = ID_S(Laplace_S>lap_selection(n-1) & Laplace_S<=lap_selection(n));
                            end
                            handles.SECTIONS_S{handles.number_of_sections} = ID_S(Laplace_S>lap_selection(end) & Laplace_S<1);


                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            g = figure(1);
                            set(g,'Color','k','Name','Sections')
                            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
                            hold on
                            for n=1:length(lap_selection)+1
                                if mod(n,2)==0
                                    plot3(handles.nodes(handles.SECTIONS_V{n},1),handles.nodes(handles.SECTIONS_V{n},2),handles.nodes(handles.SECTIONS_V{n},3),'*r','LineWidth',2)
                                else
                                    plot3(handles.nodes(handles.SECTIONS_V{n},1),handles.nodes(handles.SECTIONS_V{n},2),handles.nodes(handles.SECTIONS_V{n},3),'*b','LineWidth',2)
                                end
                            end
                            hold off
                            axis vis3d
                            daspect([1,1,1])
                            axis off
                            view([-34,-51])
                            title('\color{white}Suface Sections')
                            waitfor(g)

                            answer2 = questdlg('Do you agree with the sections generated? ', 'Question','Yes','No','No');

                            switch answer2
                                case 'Yes'
                                    handles.segments = zeros(size(handles.nodes(:,1)));
                                    for n=1:length(handles.SECTIONS_V)
                                        handles.segments(handles.SECTIONS_V{n}) = handles.segments(handles.SECTIONS_V{n})+n;
                                    end
%                                     close(g)
                                    break
                                case 'No'

                            end

                        else
                            
                            ID_V = [1:size(handles.nodes,1)]';
                            handles.SECTIONS_V{1} = ID_V(handles.Laplace>0 & handles.Laplace<1);
                            
                            % surface sections
                            ID_S = sort(unique(handles.faces(:)));
                            Laplace_S = handles.Laplace(ID_S);
                            handles.SECTIONS_S{1} = ID_S(Laplace_S>0 & Laplace_S<1);
                            
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            g = figure(1);
                            set(g,'Color','k','Name','Sections')
                            patch('faces',handles.faces,'vertices',handles.nodes,'facecolor','none','edgecolor',[0.5 0.5 0.5])
                            hold on
                            plot3(handles.nodes(handles.SECTIONS_V{1},1),handles.nodes(handles.SECTIONS_V{1},2),handles.nodes(handles.SECTIONS_V{1},3),'*r','LineWidth',2)
                            hold off
                            axis vis3d
                            daspect([1,1,1])
                            axis off
                            view([-34,-51])
                            title('\color{white}Suface Sections')
                            waitfor(g)
                            
                            answer2 = questdlg('Do you agree with the sections generated? ', 'Question','Yes','No','No');
                            
                            switch answer2
                                case 'Yes'
                                    handles.segments = zeros(size(handles.nodes(:,1)));
                                    for n=1:length(handles.SECTIONS_V)
                                        handles.segments(handles.SECTIONS_V{n}) = handles.segments(handles.SECTIONS_V{n})+n;
                                    end
                                    break
                                case 'No'
                            end
                        end
                    end
                end
            case 'No'
                
                handles.id_csv_save = 0;

        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    GUIDE_SAVE(input)
    
    if handles.id_mesh == 1
        handles.veset(unique(sort(handles.faces(:))),:,:) = handles.veset(unique(sort(handles.faces(:))),:,:).*0;    
    end

    
    s_selection = getappdata(0,'save_selection'); % id of the selection saved
    handles.save_selection = [s_selection(:,2),s_selection(:,3),s_selection(:,4),s_selection(:,1)];
    
    if sum(handles.save_selection(:))==0
        msgbox('The data was not saved ...','Warning','warn')
    else
        directory = uigetdir(pwd, 'Select Directory');

        if directory~=0
            if sum(handles.save_selection(:,1))>0
                h = waitbar(0,['Saving MAT Files ',num2str(0),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                steps = 1;
                st = steps;

                if handles.save_selection(1,1)==1
                    mkdir([directory,'/MATLAB FILES/Segmentation ROI/'])
                    SEG = handles.SEG;

                    % Julio Sotelo
                    if isempty(handles.position_sag_cor)==0
                        position_sag_cor = handles.position_sag_cor;
                        position_axi_cor = handles.position_axi_cor;
                        position_cor_cor = handles.position_cor_cor;

                        save([directory,'/MATLAB FILES/Segmentation ROI/Coor.mat'],'position_sag_cor','position_axi_cor','position_cor_cor')
                    end

                    VoxelSize = handles.voxel_MR;

                    save([directory,'/MATLAB FILES/Segmentation ROI/SEG.mat'],'SEG')
                    save([directory,'/MATLAB FILES/Segmentation ROI/VoxelSize.mat'],'VoxelSize')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(2,1)==1
                    mkdir([directory,'/MATLAB FILES/Angiography ROI/'])
                    IPCMRA = handles.IPCMRA;
                    save([directory,'/MATLAB FILES/Angiography ROI/IPCMRA.mat'],'IPCMRA')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(3,1)==1
                    mkdir([directory,'/MATLAB FILES/Magnitude ROI/'])
                    MR_FFE_FH = handles.MR_FFE_FH;
                    MR_FFE_AP = handles.MR_FFE_AP;
                    MR_FFE_RL = handles.MR_FFE_RL;
                    save([directory,'/MATLAB FILES/Magnitude ROI/MR_FFE_FH.mat'],'MR_FFE_FH','-v7.3')
                    save([directory,'/MATLAB FILES/Magnitude ROI/MR_FFE_AP.mat'],'MR_FFE_AP','-v7.3')
                    save([directory,'/MATLAB FILES/Magnitude ROI/MR_FFE_RL.mat'],'MR_FFE_RL','-v7.3')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(4,1)==1
                    mkdir([directory,'/MATLAB FILES/Velocity ROI/'])
                    MR_PCA_FH = handles.MR_PCA_FH;
                    MR_PCA_AP = handles.MR_PCA_AP;
                    MR_PCA_RL = handles.MR_PCA_RL;
                    save([directory,'/MATLAB FILES/Velocity ROI/MR_PCA_FH.mat'],'MR_PCA_FH','-v7.3')
                    save([directory,'/MATLAB FILES/Velocity ROI/MR_PCA_AP.mat'],'MR_PCA_AP','-v7.3')
                    save([directory,'/MATLAB FILES/Velocity ROI/MR_PCA_RL.mat'],'MR_PCA_RL','-v7.3')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(5,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Mesh/'])
                    nodes = handles.nodes/1000;
                    faces = handles.faces;
                    elem = handles.elem;
                    save([directory,'/MATLAB FILES/FE Mesh/nodes.mat'],'nodes')
                    save([directory,'/MATLAB FILES/FE Mesh/faces.mat'],'faces')
                    save([directory,'/MATLAB FILES/FE Mesh/elem.mat'],'elem')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(6,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Velocity/'])
                    VEL = handles.veset;
                    save([directory,'/MATLAB FILES/FE Velocity/VEL.mat'],'VEL')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(7,1)==1
                    mkdir([directory,'/MATLAB FILES/FE WSS/'])
                    WSS = handles.WSS;
                    save([directory,'/MATLAB FILES/FE WSS/WSS.mat'],'WSS')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(8,1)==1
                    mkdir([directory,'/MATLAB FILES/FE OSI/'])
                    OSI = handles.OSI;
                    save([directory,'/MATLAB FILES/FE OSI/OSI.mat'],'OSI')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(9,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Vorticity/'])
                    VOR = handles.VOR;
                    save([directory,'/MATLAB FILES/FE Vorticity/VOR.mat'],'VOR')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(10,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Helicity Density/'])
                    HD = handles.HD;
                    save([directory,'/MATLAB FILES/FE Helicity Density/HD.mat'],'HD')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(11,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Relative Helicity Density/'])
                    RHD = handles.RHD;
                    save([directory,'/MATLAB FILES/FE Relative Helicity Density/RHD.mat'],'RHD')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(12,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Viscous Dissipation/'])
                    VD = handles.mags_vd;
                    save([directory,'/MATLAB FILES/FE Viscous Dissipation/VD.mat'],'VD')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(13,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Energy Loss/'])
                    EL = handles.mags_el;
                    save([directory,'/MATLAB FILES/FE Energy Loss/EL.mat'],'EL')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(14,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Kinetic Energy/'])
                    KE = handles.mags_ke;
                    save([directory,'/MATLAB FILES/FE Kinetic Energy/KE.mat'],'KE')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(15,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Laplace/'])
                    Laplace = handles.Laplace;
                    save([directory,'/MATLAB FILES/FE Laplace/Laplace.mat'],'Laplace')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(16,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Centerline/'])
                    Centerline = handles.centerline/1000;
                    save([directory,'/MATLAB FILES/FE Centerline/Centerline.mat'],'Centerline')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(17,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Radius/'])
                    Radius = handles.radius;
                    save([directory,'/MATLAB FILES/FE Radius/Radius.mat'],'Radius')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(18,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Diameter/'])
                    Diameter = handles.diameter;
                    save([directory,'/MATLAB FILES/FE Diameter/Diameter.mat'],'Diameter')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(19,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Axial Unit Vectors/'])
                    Axial_unit_vectors = handles.axial_unit_vectors;
                    save([directory,'/MATLAB FILES/FE Axial Unit Vectors/Axial_unit_vectors.mat'],'Axial_unit_vectors')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(20,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Circumferential Unit Vectors/'])
                    Circumferential_unit_vectors = handles.circumferential_unit_vectors;
                    save([directory,'/MATLAB FILES/FE Circumferential Unit Vectors/Circumferential_unit_vectors.mat'],'Circumferential_unit_vectors')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(21,1)==1
                    mkdir([directory,'/MATLAB FILES/FE WSS Axial/'])
                    WSS_A = handles.WSS_A;
                    save([directory,'/MATLAB FILES/FE WSS Axial/WSS_A.mat'],'WSS_A')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(22,1)==1
                    mkdir([directory,'/MATLAB FILES/FE WSS Circumferential/'])
                    WSS_C = handles.WSS_C;
                    save([directory,'/MATLAB FILES/FE WSS Circumferential/WSS_C.mat'],'WSS_C')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(23,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Axial Angle/'])
                    Angle_axial_direction = handles.angle_axial_direction;
                    save([directory,'/MATLAB FILES/FE Axial Angle/Angle_axial_direction.mat'],'Angle_axial_direction')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(24,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Forward Velocity/'])
                    Forward_velocity = handles.forward_velocity;
                    save([directory,'/MATLAB FILES/FE Forward Velocity/Forward_velocity.mat'],'Forward_velocity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(25,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Backward Velocity/'])
                    Backward_velocity = handles.backward_velocity;
                    save([directory,'/MATLAB FILES/FE Backward Velocity/Backward_velocity.mat'],'Backward_velocity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(26,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Regurgitant Flow/'])
                    Regurgitant_flow = handles.regurgitant_flow;
                    save([directory,'/MATLAB FILES/FE Regurgitant Flow/Regurgitant_flow.mat'],'Regurgitant_flow')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(27,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Centerline Flow/'])
                    Centerline_flow = handles.centerline_flow/1000;
                    save([directory,'/MATLAB FILES/FE Centerline Flow/Centerline_flow.mat'],'Centerline_flow')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(28,1)==1
                    mkdir([directory,'/MATLAB FILES/FE Eccentricity/'])
                    Eccentricity = handles.eccentricity;
                    save([directory,'/MATLAB FILES/FE Eccentricity/Eccentricity.mat'],'Eccentricity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end

                if handles.save_selection(29,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Curvature/'])
                    Curvature = handles.curvature;
                    save([directory,'/MATLAB FILES/FE Curvature/Curvature.mat'],'Curvature')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(30,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Ellipticity/'])
                    Ellipticity = handles.ellipticity;
                    save([directory,'/MATLAB FILES/FE Ellipticity/Ellipticity.mat'],'Ellipticity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(31,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Length Vessel/'])
                    Length_vessel = handles.length_vessel;
                    save([directory,'/MATLAB FILES/FE Length Vessel/Length_vessel.mat'],'Length_vessel')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
    %             if handles.save_selection(32,1)==1 % Julio Sotelo 28-05-2019
    %                 mkdir([directory,'/MATLAB FILES/FE Circulation/'])
    %                 Circulation = handles.circulation;
    %                 save([directory,'/MATLAB FILES/FE Circulation/Circulation.mat'],'Circulation')
    %                 waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
    %                 st = st+steps;
    %             end
                if handles.save_selection(32,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Forward Vortex/'])
                    Forward_vortex = handles.forward_vortex;
                    save([directory,'/MATLAB FILES/FE Forward Vortex/Forward_vortex.mat'],'Forward_vortex')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(33,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Flattening/'])
                    Flattening = handles.flattening;
                    save([directory,'/MATLAB FILES/FE Flattening/Flattening.mat'],'Flattening')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(34,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Area/'])
                    Area = handles.area;
                    save([directory,'/MATLAB FILES/FE Area/Area.mat'],'Area')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(35,1)==1 % Julio Sotelo 28-05-2019
                    mkdir([directory,'/MATLAB FILES/FE Axial Circulation/'])
                    Axial_circulation = handles.axial_circulation;
                    save([directory,'/MATLAB FILES/FE Axial Circulation/Axial_circulation.mat'],'Axial_circulation')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if handles.save_selection(36,1)==1 % Julio Sotelo 04-07-2019
                    mkdir([directory,'/MATLAB FILES/2D Flow/'])
                    Time = handles.time;
                    save([directory,'/MATLAB FILES/2D Flow/Time.mat'],'Time')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(37,1)==1 % Julio Sotelo 04-07-2019
                    mkdir([directory,'/MATLAB FILES/2D Flow/'])
                    Flow = handles.flow;
                    save([directory,'/MATLAB FILES/2D Flow/Flow.mat'],'Flow')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(38,1)==1 % Julio Sotelo 04-07-2019
                    mkdir([directory,'/MATLAB FILES/2D Flow/'])
                    Net_Flow = handles.net_flow;
                    save([directory,'/MATLAB FILES/2D Flow/Net_Flow.mat'],'Net_Flow')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(39,1)==1 % Julio Sotelo 04-07-2019
                    mkdir([directory,'/MATLAB FILES/2D Flow/'])
                    Maximum_Velocity = handles.max_velocity;
                    save([directory,'/MATLAB FILES/2D Flow/Maximum_Velocity.mat'],'Maximum_Velocity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                if handles.save_selection(40,1)==1 % Julio Sotelo 04-07-2019
                    mkdir([directory,'/MATLAB FILES/2D Flow/'])
                    Minimum_Velocity = handles.min_velocity;
                    save([directory,'/MATLAB FILES/2D Flow/Minimum_Velocity.mat'],'Minimum_Velocity')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if handles.id_csv_save == 1
                    mkdir([directory,'/MATLAB FILES/FE Segments/'])
                    Segments = handles.segments;
                    save([directory,'/MATLAB FILES/FE Segments/Segments.mat'],'Segments')
                    waitbar(st/ sum(handles.save_selection(:,1)),h,['Saving MAT Files ',num2str(st),' files of ',num2str(sum(handles.save_selection(:,1))),' ...']);
                    st = st+steps;
                end 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                close(h)
            end

            %%%% VTI

            if sum(handles.save_selection(:,2))>0
                h = waitbar(0,['Saving VTI Files ',num2str(0),' files of ',num2str(size(handles.MR_FFE_FH,4)),' ...']);
                steps = 1;
                st = steps;
                mkdir([directory,'/VTI FILES/'])

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for n=1:size(handles.MR_FFE_FH,4)
    %             for n=1
    %                 n = handles.peak_flow_ori;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    name = sprintf(['Volume_','%04d','.vti'],n);                
                    fileID = fopen([directory,'/VTI FILES/',name],'w');
                    fprintf(fileID,'%s\n','<?xml version="1.0"?>');
                    fprintf(fileID,'%s\n','<VTKFile type="ImageData"  version="0.1"  >');
                    fprintf(fileID,'%s\n',['<ImageData  WholeExtent="1 ',num2str(size(handles.SEG,1)),' 1 ',num2str(size(handles.SEG,2)),' 1 ',num2str(size(handles.SEG,3)),'" Origin="',num2str(-handles.voxel_MR(1)/1000),' ',num2str(-handles.voxel_MR(2)/1000),' ',num2str(-handles.voxel_MR(3)/1000),'" Spacing="',num2str(handles.voxel_MR(1)/1000),' ',num2str(handles.voxel_MR(2)/1000),' ',num2str(handles.voxel_MR(3)/1000),'" >']);
                    fprintf(fileID,'%s\n',['<Piece  Extent="1 ',num2str(size(handles.SEG,1)),' 1 ',num2str(size(handles.SEG,2)),' 1 ',num2str(size(handles.SEG,3)),'" >']);

                    scalar_variables = [1 2 3];
                    vector_variables = 4;
                    tensor_variables = [];
                    names_scalar_variables = {'Segmentation', 'IPCMRA', 'Magnitude Image' };
                    names_vector_variables = {'Velocity [m/s]'};
                    names_tensor_variables = {};

                    if sum(handles.save_selection(scalar_variables,2))==0
                        name1 = 'Scalars=""';
                    else
                        [r,~,~]=find(handles.save_selection(scalar_variables,2)==1);
                        name1 = ['Scalars="',names_scalar_variables{r(1)},'"'];
                    end
                    if sum(handles.save_selection(vector_variables,2))==0
                        name2 = 'Vectors=""';
                    else
                        [r,~,~]=find(handles.save_selection(vector_variables,2)==1);
                        name2 = ['Vectors="',names_vector_variables{r(1)},'"'];
                    end
                    if sum(handles.save_selection(tensor_variables,2))==0
                        name3 = 'Tensor=""';
                    else
                        [r,~,~]=find(handles.save_selection(tensor_variables,2)==1);
                        name3 = ['Tensor="',names_tensor_variables{r(1)},'"'];
                    end

                    fprintf(fileID,'%s\n','<PointData ',name1,' ',name2,' ',name3,'>');

                    if handles.save_selection(1,2)==1
                        fprintf(fileID,'%s\n','<DataArray type="UInt8" Name="Segmentation" format="ascii">');
                        fprintf(fileID,'%d\n',uint8(handles.SEG(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(2,2)==1
                        fprintf(fileID,'%s\n','<DataArray type="Float32" Name="IPCMRA" format="ascii">');
                        fprintf(fileID,'%f\n',handles.IPCMRA(:)');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(3,2)==1
                        fprintf(fileID,'%s\n','<DataArray type="UInt8" Name="Magnitude Image" format="ascii">');
                        var = squeeze(handles.MR_FFE_FH(:,:,:,n));
                        var = (var/max(var(:)))*255;
                        fprintf(fileID,'%d\n',uint8(var(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(4,2)==1
                        fprintf(fileID,'%s\n','<DataArray type="Float32" Name="Velocity [cm/s]" NumberOfComponents="3" format="ascii">');
                        var1 = squeeze(handles.MR_PCA_FH(:,:,:,n))*-1/100;
                        var2 = squeeze(handles.MR_PCA_RL(:,:,:,n))*-1/100;
                        var3 = squeeze(handles.MR_PCA_AP(:,:,:,n))/100;
                        var = [var3(:),var1(:),var2(:)];
                        fprintf(fileID,'%f %f %f\n',var');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end

                    fprintf(fileID,'%s\n','</PointData>');
                    fprintf(fileID,'%s\n','<CellData>');
                    fprintf(fileID,'%s\n','</CellData>');
                    fprintf(fileID,'%s\n','</Piece>');
                    fprintf(fileID,'%s\n','</ImageData>');
                    fprintf(fileID,'%s\n','</VTKFile>');

                    waitbar(st/ size(handles.MR_FFE_FH,4),h,['Saving VTI Files ',num2str(st),' files of ',num2str(size(handles.MR_FFE_FH,4)),' ...']);
                    st = st+steps;

                end

                name_all_vti = 'all_vti.pvd';
                fileID_2 = fopen([directory,'/VTI FILES/',name_all_vti],'w');
                fprintf(fileID_2,'%s\n','<?xml version="1.0"?>');
                fprintf(fileID_2,'%s\n','<VTKFile type="Collection" version="0.1">');
                fprintf(fileID_2,'%s\n','<Collection>');
                dt = (60/handles.heart_rate)/(size(handles.MR_FFE_FH,4)-1);
                for n = 1:size(handles.MR_FFE_FH,4)
                    name1 = sprintf(['Volume_','%04d','.vti'],n);
                    delta_time = squeeze(dt)*(n-1);
                    fprintf(fileID_2,'%s\n',['<DataSet timestep="',num2str(delta_time,'%0.8f'),'" part="0" file="',name1,'" />']);
                end
                fprintf(fileID_2,'%s\n','</Collection>');
                fprintf(fileID_2,'%s\n','</VTKFile>');

                close(h)
            end
            if sum(handles.save_selection(:,3))>0
                mkdir([directory,'/VTU FILES/'])
                elem_n = handles.elem-1;
                faces_n = handles.faces-1;
                nodes_n = handles.nodes/1000;
                ID_nodes = [1:size(nodes_n,1)]';
                tetra = elem_n;
                cell_type1 = zeros(1,size(elem_n,1))'+10;
                cell_ID1 = ([1:length(tetra(:,1))]*4)';
                tria = faces_n;
                cell_type2 = zeros(1,size(faces_n,1))'+5;
                cell_ID2 = ([1:length(tria(:,1))]*3)'+ cell_ID1(end);
                h = waitbar(0,['Saving VTU Files ',num2str(0),' files of ',num2str(size(handles.veset,3)),' ...']);
                steps = 1;
                st = steps;
                nodevol = squeeze(nodevolume(nodes_n,handles.elem))*1e+9;% variable is saved in mm3
                for n =1:size(handles.veset,3)
                    name = sprintf(['vessel_','%04d','.vtu'],n);
                    fileID = fopen([directory,'/VTU FILES/',name],'w');
                    fprintf(fileID,'%s\n','<?xml version="1.0"?>');
                    fprintf(fileID,'%s\n','<VTKFile type="UnstructuredGrid"  version="0.1"  >');
                    fprintf(fileID,'%s\n','<UnstructuredGrid>');
                    fprintf(fileID,'%s\n',['<Piece  NumberOfPoints="',num2str(length(nodes_n(:,1))),'" NumberOfCells="',num2str(length(elem_n(:,1)) + length(faces_n(:,1))),'">']);
                    fprintf(fileID,'%s\n','<Points>');
                    fprintf(fileID,'%s\n','<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii">');
                    fprintf(fileID,'%f %f %f\n',nodes_n');
                    fprintf(fileID,'%s\n','</DataArray>');
                    fprintf(fileID,'%s\n','</Points>');
                    fprintf(fileID,'%s\n','<Cells>');
                    fprintf(fileID,'%s\n','<DataArray  type="UInt32"  Name="connectivity"  format="ascii">');
                    fprintf(fileID,'%d %d %d %d\n',tetra');
                    fprintf(fileID,'%d %d %d\n',tria');
                    fprintf(fileID,'%s\n','</DataArray>');
                    fprintf(fileID,'%s\n','<DataArray  type="UInt32"  Name="offsets"  format="ascii">');
                    fprintf(fileID,'%d\n',cell_ID1');
                    fprintf(fileID,'%d\n',cell_ID2');
                    fprintf(fileID,'%s\n','</DataArray>');
                    fprintf(fileID,'%s\n','<DataArray  type="UInt8"  Name="types"  format="ascii">');
                    fprintf(fileID,'%d\n',cell_type1');
                    fprintf(fileID,'%d\n',cell_type2');
                    fprintf(fileID,'%s\n','</DataArray>');
                    fprintf(fileID,'%s\n','</Cells>');
                    scalar_variables = [8 10 11 12 13 14 15 17 18 23 26 28 29 30 31 32 33 34 35];
                    vector_variables = [6 7 9 19 20 21 22 24 25];
                    tensor_variables = [];
                    names_scalar_variables = {'OSI [-]', 'Helicity Density [m/s$^2$]', 'Relative Helicity Density [-]','Viscouss Dissipation [1e$^3$/s$^2$]','Energy Loss [$/mu$ W]','Kinetic Energy [$/mu$ J]',...
                                              'Laplace [-]','Radius [cm]','Diameter [cm]','Axial Angle [$^o$]','Regurgitant Flow [%]','Eccentricity [%]','Curvature [1/m]','Ellipticity [-]','Length of Vessel [cm]','Forward Vortex [1/s]','Flattening [-]','Area [cm$^2$]','Axial Circulation [cm$^2$/s]'};
                    names_vector_variables = {'Velocity [m/s]', 'WSS [N/m$^2$]', 'Vorticity [1/s]',...
                                              'Axial Unit Vector [-]','Circum. Unit Vector [-]','WSS-A [N/m$^2$]','WSS-C [N/m$^2$]','Forward Velocity [m/s]','Backward Velocity [m/s]'};

                    names_tensor_variables = {};

                    if sum(handles.save_selection(scalar_variables,3))==0
                        name1 = 'Scalars="Node Volume [m$^3$]"';
                    else
                        [r,c,v]=find(handles.save_selection(scalar_variables,3)==1);
                        name1 = ['Scalars="',names_scalar_variables{r(1)},'"'];
                    end
                    if sum(handles.save_selection(vector_variables,3))==0
                        name2 = 'Vectors=""';
                    else
                        [r,c,v]=find(handles.save_selection(vector_variables,3)==1);
                        name2 = ['Vectors="',names_vector_variables{r(1)},'"'];
                    end
                    if sum(handles.save_selection(tensor_variables,3))==0
                        name3 = 'Tensor=""';
                    else
                        [r,c,v]=find(handles.save_selection(tensor_variables,3)==1);
                        name3 = ['Tensor="',names_tensor_variables{r(1)},'"'];
                    end


                    fprintf(fileID,'%s\n','<PointData ',name1,' ',name2,' ',name3,'>');
                    if handles.save_selection(scalar_variables(1),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="OSI [-]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.OSI(:,1))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(2),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Helicity Density [m/s$^2$]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.HD(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(3),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Relative Helicity Density [-]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.RHD(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(4),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Viscouss Dissipation [1e$^3$/s$^2$]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.mags_vd(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(5),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Energy Loss [$\mu$ W]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.mags_el(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(6),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Kinetic Energy [$\mu$ J]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.mags_ke(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(7),3)==1 % Laplace [-]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Laplace [-]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.Laplace(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(8),3)==1 % Radius [mm]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Radius [cm]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.radius(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(9),3)==1 % 'Diameter [mm]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Diameter [cm]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.diameter(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(10),3)==1 % 'Axial Angle [^{o}]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Axial Angle [$^o$]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.angle_axial_direction(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(11),3)==1 % 'Regurgitant Flow [%]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Regurgitant Flow [%]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.regurgitant_flow(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(12),3)==1 % 'Eccentricity [%]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Eccentricity [%]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.eccentricity(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(13),3)==1 % 'Curvature [1/m]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Curvature [1/m]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.curvature(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(14),3)==1 % 'Ellipticity [-]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Ellipticity [-]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.ellipticity(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(15),3)==1 % 'Length Vessel [mm]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Length Vessel [cm]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.length_vessel(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
    %                 if handles.save_selection(scalar_variables(16),3)==1 % 'Circulation [mm$^2$/s]'
    %                     fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Circulation [mm$^2$/s]"  NumberOfComponents="1" format="ascii">');
    %                     fprintf(fileID,'%f\n',squeeze(handles.circulation(:,n))');
    %                     fprintf(fileID,'%s\n','</DataArray>');
    %                 end
                    if handles.save_selection(scalar_variables(16),3)==1 % 'Forward Vortex [1/s]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Forward Vortex [1/s]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.forward_vortex(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(17),3)==1 % 'Flattening [-]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Flattening [-]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.flattening(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(18),3)==1 % 'Area [mm$^2$]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Area [cm$^2$]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.area(:))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(scalar_variables(19),3)==1 % 'Axial Circulation [mm$^2$/s]'
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Axial Circulation [cm$^2$/s]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.axial_circulation(:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end

                    fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Node Volume [mm$^3$]"  NumberOfComponents="1" format="ascii">');
                    fprintf(fileID,'%f\n',squeeze(nodevol(:,1))');
                    fprintf(fileID,'%s\n','</DataArray>');
                    fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="ID nodes[#]"  NumberOfComponents="1" format="ascii">');
                    fprintf(fileID,'%f\n',squeeze(ID_nodes(:,1))');
                    fprintf(fileID,'%s\n','</DataArray>');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if handles.id_csv_save == 1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Segments[#]"  NumberOfComponents="1" format="ascii">');
                        fprintf(fileID,'%f\n',squeeze(handles.segments(:)));
                        fprintf(fileID,'%s\n','</DataArray>');
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if handles.save_selection(vector_variables(1),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Velocity [m/s]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.veset(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(2),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="WSS [N/m$^2$]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.WSS(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(3),3)==1
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Vorticity [1/s]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.VOR(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end


                    if handles.save_selection(vector_variables(4),3)==1 %Axial Unit Vector [-]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Axial Unit Vector [-]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.axial_unit_vectors)');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(5),3)==1 %Circum. Unit Vector [-]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Circum. Unit Vector [-]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.circumferential_unit_vectors)');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(6),3)==1 %WSS-A [N/m$^2$]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="WSS-A [N/m$^2$]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.WSS_A(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(7),3)==1 %WSS-C [N/m$^2$]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="WSS-C [N/m$^2$]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.WSS_C(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(8),3)==1 %Forward Velocity [m/s]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Forward Velocity [m/s]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.forward_velocity(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end
                    if handles.save_selection(vector_variables(9),3)==1 %Backward Velocity [m/s]
                        fprintf(fileID,'%s\n','<DataArray  type="Float64"  Name="Backward Velocity [m/s]"  NumberOfComponents="3" format="ascii">');
                        fprintf(fileID,'%f %f %f\n',squeeze(handles.backward_velocity(:,:,n))');
                        fprintf(fileID,'%s\n','</DataArray>');
                    end

                    fprintf(fileID,'%s\n','</PointData>');
                    fprintf(fileID,'%s\n','</Piece>');
                    fprintf(fileID,'%s\n','</UnstructuredGrid>');
                    fprintf(fileID,'%s\n','</VTKFile>');
                    if st==n
                        waitbar(st/ size(handles.veset,3),h,['Saving VTU Files ',num2str(st),' files of ',num2str(size(handles.veset,3)),' ...']);
                        st = st+steps;
                    end
                end
                name_all_vtu = 'all_vtu.pvd';
                fileID_2 = fopen([directory,'/VTU FILES/',name_all_vtu],'w');
                fprintf(fileID_2,'%s\n','<?xml version="1.0"?>');
                fprintf(fileID_2,'%s\n','<VTKFile type="Collection" version="0.1">');
                fprintf(fileID_2,'%s\n','<Collection>');
                dt = (60/handles.heart_rate)/(size(handles.veset,3)-1);
                for n = 1:size(handles.veset,3)
                    name1 = sprintf(['vessel_','%04d','.vtu'],n);
                    delta_time = squeeze(dt)*(n-1);
                    fprintf(fileID_2,'%s\n',['<DataSet timestep="',num2str(delta_time,'%0.8f'),'" part="0" file="',name1,'" />']);
                end
                fprintf(fileID_2,'%s\n','</Collection>');
                fprintf(fileID_2,'%s\n','</VTKFile>');

                if handles.save_selection(16,3)==1

                    ID = [1:length(handles.centerline(:,1))]-1;
                    M = [ID(1:end-1)',ID(2:end)'];
                    pp = zeros(1,size(handles.centerline,1)-1)'+2;
                    cell = [pp,M];

                    cell_type = zeros(1,size(handles.centerline,1)-1)'+3;

                    name_c = 'centerline.vtk';
                    fileID = fopen([directory,'/VTU FILES/',name_c],'w');
                    fprintf(fileID,'%s\n','# vtk DataFile Version 2.0');
                    fprintf(fileID,'%s\n',['Unstructured Grid example']);
                    fprintf(fileID,'%s\n','ASCII');
                    fprintf(fileID,'%s\n','DATASET UNSTRUCTURED_GRID');
                    fprintf(fileID,'%s\n',['POINTS ',num2str(length(handles.centerline(:,1))),' float']);
                    fprintf(fileID,'%f %f %f\n',(handles.centerline/1000)');
                    fprintf(fileID,'%s\n',['CELLS ',num2str(size(cell,1)),' ',num2str(size(cell,1)*size(cell,2))]);
                    fprintf(fileID,'%d %d %d\n',cell');
                    fprintf(fileID,'%s\n',['CELL_TYPES ',num2str(size(cell,1))]);
                    fprintf(fileID,'%d\n',cell_type');

                end

                if handles.save_selection(27,3)==1

                    ID = [1:length(handles.centerline_flow(:,1))]-1;
                    M = [ID(1:end-1)',ID(2:end)'];
                    pp = zeros(1,size(handles.centerline_flow,1)-1)'+2;
                    cell = [pp,M];

                    cell_type = zeros(1,size(handles.centerline_flow,1)-1)'+3;

                    name_c = 'centerline_flow.vtk';
                    fileID = fopen([directory,'/VTU FILES/',name_c],'w');
                    fprintf(fileID,'%s\n','# vtk DataFile Version 2.0');
                    fprintf(fileID,'%s\n',['Unstructured Grid example']);
                    fprintf(fileID,'%s\n','ASCII');
                    fprintf(fileID,'%s\n','DATASET UNSTRUCTURED_GRID');
                    fprintf(fileID,'%s\n',['POINTS ',num2str(length(handles.centerline_flow(:,1))),' float']);
                    fprintf(fileID,'%f %f %f\n',(handles.centerline_flow/1000)');
                    fprintf(fileID,'%s\n',['CELLS ',num2str(size(cell,1)),' ',num2str(size(cell,1)*size(cell,2))]);
                    fprintf(fileID,'%d %d %d\n',cell');
                    fprintf(fileID,'%s\n',['CELL_TYPES ',num2str(size(cell,1))]);
                    fprintf(fileID,'%d\n',cell_type');

                end

                close(h)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(handles.save_selection(:,4))>0

                mkdir([directory,'/XLS FILES/'])
                nodes_n = handles.nodes/1000;
                nodevol = squeeze(nodevolume(nodes_n,handles.elem));

                vol_V = [];
                vol_V_n = [];
                vol_S = [];
                id_vol_V_n =[];

                for n=1:handles.number_of_sections
                    vol_V{n} = nodevol(handles.SECTIONS_V{n})  /sum(nodevol(handles.SECTIONS_V{n}));
                    vol_S{n} = nodevol(handles.SECTIONS_S{n})  /sum(nodevol(handles.SECTIONS_S{n}));

                    temp = [];
                    [~,ia,~] = intersect(handles.SECTIONS_V{n},handles.list_n,'legacy');
                    temp = handles.SECTIONS_V{n};
                    temp(ia) = [];
                    id_vol_V_n{n} = temp;
                    vol_V_n{n} = nodevol(id_vol_V_n{n})  /sum(nodevol(id_vol_V_n{n}));

                end

                peak_systole = handles.peak_flow_ori;

                % MAGNUTUD OF VECTORS
                mag_velocity        = squeeze(sqrt(sum(handles.veset.^2,2)));% +1 time
                mag_wss             = squeeze(sqrt(sum(handles.WSS.^2,2)));% +1 time
                mag_wssa            = squeeze(sqrt(sum(handles.WSS_A.^2,2)));% +1 time
                mag_wssc            = squeeze(sqrt(sum(handles.WSS_C.^2,2)));% +1 time
                mag_vor             = squeeze(sqrt(sum(handles.VOR.^2,2)));% +1 time
                mag_forward_vel     = squeeze(sqrt(sum(handles.forward_velocity.^2,2)));% +1 time
                mag_backward_vel    = squeeze(sqrt(sum(handles.backward_velocity.^2,2)));% +1 time
                % SCALARS
                osi_values = repmat(handles.OSI,1,size(handles.veset,3)); % 1 time
                hd_values = handles.HD(:,:); % +1 time
                rhd_values = handles.RHD(:,:); % +1 time
                vd_values = handles.mags_vd(:,:); % +1 time
                el_values = handles.mags_el(:,:); % +1 time
                ke_values = handles.mags_ke(:,:); % +1 time
                aa_values = handles.angle_axial_direction(:,:); % 1 time
                rf_values = repmat(handles.regurgitant_flow,1,size(handles.veset,3)); % 1 time
                ecc_values = repmat(handles.eccentricity,1,size(handles.veset,3)); % 1 time
                dia_values = repmat(handles.diameter,1,size(handles.veset,3)); % 1 time
                rad_values = repmat(handles.radius,1,size(handles.veset,3)); % 1 time
                cur_values = repmat(handles.curvature,1,size(handles.veset,3)); % Julio Sotelo 28-05-2019
                ell_values = repmat(handles.ellipticity,1,size(handles.veset,3)); % Julio Sotelo 28-05-2019
                len_values = repmat(handles.length_vessel,1,size(handles.veset,3)); % Julio Sotelo 28-05-2019
    %             cir_values = handles.circulation(:,:); % Julio Sotelo 28-05-2019
                fov_values = handles.forward_vortex(:,:); % Julio Sotelo 28-05-2019
                fla_values = repmat(handles.flattening,1,size(handles.veset,3)); % Julio Sotelo 28-05-2019
                are_values = repmat(handles.area,1,size(handles.veset,3)); % Julio Sotelo 28-05-2019
                aci_values = handles.axial_circulation(:,:); % Julio Sotelo 28-05-2019

                %%% Parameter names
                Parameter = {   'Velocity [m/s]';...
                                'WSS [N/m^2]';...
                                'WSSA [N/m^2]';...
                                'WSSC [N/m^2]';...
                                'Vorticity [1/s]';...
                                'Forward Velocity [m/s]';...
                                'Backward Velocity [m/s]';...
                                'OSI [-]';...
                                'Helicity Density [m/s^2]';...
                                'Relative Helicity Density [-]';...
                                'Viscous Dissipation [1/s^2]';...
                                'Energy Loss [uW]';...
                                'Kinetic Energy [uJ]';...
                                'Angle in Axial Direction [°]';...
                                'Regurgitant Flow [%]';...
                                'Eccentricity [%]';...
                                'Diameter [cm]';...
                                'Radius [cm]';...
                                'Curvature [1/m]';...
                                'Ellipticity [-]';...
                                'Length Vessel [cm]';...
                                'Forward Vortex [1/s]';...
                                'Flattering [-]';...
                                'Area [cm^2]';...
                                'Axial Circulation [cm^2/s]';...
                                'Peak Systole CP#'};




                %%% Mean Values
                MAT = zeros(26,handles.number_of_sections,size(handles.veset,3));
                for m=1:size(handles.veset,3)
                    for n=1:handles.number_of_sections

                        MAT(1,n,m) = sum(squeeze(mag_velocity(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(2,n,m) = sum(squeeze(mag_wss(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(3,n,m) = sum(squeeze(mag_wssa(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(4,n,m) = sum(squeeze(mag_wssc(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(5,n,m) = sum(squeeze(mag_vor(id_vol_V_n{n},m)).*vol_V_n{n}); % vorticity
                        MAT(6,n,m) = sum(squeeze(mag_forward_vel(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(7,n,m) = sum(squeeze(mag_backward_vel(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(8,n,m) = sum(squeeze(osi_values(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(9,n,m) = sum(squeeze(hd_values(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(10,n,m) = sum(squeeze(rhd_values(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(11,n,m) = sum(squeeze(vd_values(id_vol_V_n{n},m)).*vol_V_n{n}); % viscous dissipation
                        MAT(12,n,m) = sum(squeeze(el_values(id_vol_V_n{n},m)).*vol_V_n{n}); % energy loss 
                        MAT(13,n,m) = sum(squeeze(ke_values(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(14,n,m) = sum(squeeze(aa_values(handles.SECTIONS_V{n},m)).*vol_V{n});
                        MAT(15,n,m) = sum(squeeze(rf_values(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(16,n,m) = sum(squeeze(ecc_values(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(17,n,m) = sum(squeeze(dia_values(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(18,n,m) = sum(squeeze(rad_values(handles.SECTIONS_S{n},m)).*vol_S{n});
                        MAT(19,n,m) = sum(squeeze(cur_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(20,n,m) = sum(squeeze(ell_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(21,n,m) = sum(squeeze(len_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
    %                     MAT(22,n,m) = sum(squeeze(cir_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(22,n,m) = sum(squeeze(fov_values(id_vol_V_n{n},m)).*vol_V_n{n}); % Julio Sotelo 28-05-2019
                        MAT(23,n,m) = sum(squeeze(fla_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(24,n,m) = sum(squeeze(are_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(25,n,m) = sum(squeeze(aci_values(handles.SECTIONS_S{n},m)).*vol_S{n}); % Julio Sotelo 28-05-2019
                        MAT(26,n,m) = peak_systole;

                    end
                end


                if handles.save_selection(6,4)==0; MAT(1,:,:) = NaN; end
                if handles.save_selection(7,4)==0; MAT(2,:,:) = NaN; end
                if handles.save_selection(8,4)==0; MAT(8,:,:) = NaN; end
                if handles.save_selection(9,4)==0; MAT(5,:,:) = NaN; end
                if handles.save_selection(10,4)==0; MAT(9,:,:) = NaN; end
                if handles.save_selection(11,4)==0; MAT(10,:,:) = NaN; end
                if handles.save_selection(12,4)==0; MAT(11,:,:) = NaN; end
                if handles.save_selection(13,4)==0; MAT(12,:,:) = NaN; end
                if handles.save_selection(14,4)==0; MAT(13,:,:) = NaN; end
                if handles.save_selection(17,4)==0; MAT(18,:,:) = NaN; end
                if handles.save_selection(18,4)==0; MAT(17,:,:) = NaN; end
                if handles.save_selection(21,4)==0; MAT(3,:,:) = NaN; end
                if handles.save_selection(22,4)==0; MAT(4,:,:) = NaN; end
                if handles.save_selection(23,4)==0; MAT(14,:,:) = NaN; end
                if handles.save_selection(24,4)==0; MAT(6,:,:) = NaN; end
                if handles.save_selection(25,4)==0; MAT(7,:,:) = NaN; end
                if handles.save_selection(26,4)==0; MAT(15,:,:) = NaN; end
                if handles.save_selection(28,4)==0; MAT(16,:,:) = NaN; end
                if handles.save_selection(29,4)==0; MAT(19,:,:) = NaN; end
                if handles.save_selection(30,4)==0; MAT(20,:,:) = NaN; end
                if handles.save_selection(31,4)==0; MAT(21,:,:) = NaN; end
                if handles.save_selection(32,4)==0; MAT(22,:,:) = NaN; end
                if handles.save_selection(33,4)==0; MAT(23,:,:) = NaN; end
                if handles.save_selection(34,4)==0; MAT(24,:,:) = NaN; end
                if handles.save_selection(35,4)==0; MAT(25,:,:) = NaN; end
    %             if handles.save_selection(36,4)==0; MAT(26,:,:) = NaN; end


                h = waitbar(0,['Saving mean values (XLS file) cardiac phase ',num2str(1),' of ',num2str(size(handles.veset,3)),' ...']);
                steps = 1;
                st = steps;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VV = {'Time [s]';...
                      'Flow [ml/s]';...
                      'Net Flow [ml]';...
                      'Maximum Velocity [cm/s]';...
                      'Minimum Velocity [cm/s]'};
                MM = [handles.time';handles.flow';handles.net_flow';handles.max_velocity';handles.min_velocity'];
                TT = table(VV,MM);
                TT.Properties.VariableNames = {'Parameter','Cardiac_Phase'};

                filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MEAN_VALUES.xls'];
                writetable(TT,fullfile(filename),'Sheet','2D Flow','Range','A1')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for n=1:size(handles.veset,3)
                    eval(['Sections_from_1_to_',num2str(handles.number_of_sections),' = MAT(:,:,n);']);
                    T = table(Parameter,eval(['Sections_from_1_to_',num2str(handles.number_of_sections)]));
                    T.Properties.VariableNames = {'Parameter','Section'};

                    filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MEAN_VALUES.xls'];
                    writetable(T,fullfile(filename),'Sheet',['Cardiac Phase #',num2str(n)],'Range','A1')

                    if st==n
                        waitbar(st/ size(handles.veset,3),h,['Saving mean values (XLS file) cardiac phase ',num2str(st),' of ',num2str(size(handles.veset,3)),' ...']);
                        st = st+steps;
                    end
                end
                close(h)

                %%% Maximun Values
                MAT = zeros(26,handles.number_of_sections,size(handles.veset,3));
                for m=1:size(handles.veset,3)
                    for n=1:handles.number_of_sections

                        MAT(1,n,m) = max(mag_velocity(handles.SECTIONS_V{n},m));
                        MAT(2,n,m) = max(mag_wss(handles.SECTIONS_S{n},m));
                        MAT(3,n,m) = max(mag_wssa(handles.SECTIONS_S{n},m));
                        MAT(4,n,m) = max(mag_wssc(handles.SECTIONS_S{n},m));
                        MAT(5,n,m) = max(mag_vor(id_vol_V_n{n},m));
                        MAT(6,n,m) = max(mag_forward_vel(handles.SECTIONS_V{n},m));
                        MAT(7,n,m) = max(mag_backward_vel(handles.SECTIONS_V{n},m));
                        MAT(8,n,m) = max(osi_values(handles.SECTIONS_S{n},m));
                        MAT(9,n,m) = max(hd_values(handles.SECTIONS_V{n},m));
                        MAT(10,n,m) = max(rhd_values(handles.SECTIONS_V{n},m));
                        MAT(11,n,m) = max(vd_values(id_vol_V_n{n},m));
                        MAT(12,n,m) = max(el_values(id_vol_V_n{n},m));
                        MAT(13,n,m) = max(ke_values(handles.SECTIONS_V{n},m));
                        MAT(14,n,m) = max(aa_values(handles.SECTIONS_V{n},m));
                        MAT(15,n,m) = max(rf_values(handles.SECTIONS_S{n},m));
                        MAT(16,n,m) = max(ecc_values(handles.SECTIONS_S{n},m));
                        MAT(17,n,m) = max(dia_values(handles.SECTIONS_S{n},m));
                        MAT(18,n,m) = max(rad_values(handles.SECTIONS_S{n},m));
                        MAT(19,n,m) = max(cur_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(20,n,m) = max(ell_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(21,n,m) = max(len_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
    %                     MAT(22,n,m) = max(cir_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(22,n,m) = max(fov_values(id_vol_V_n{n},m)); % Julio Sotelo 28-05-2019
                        MAT(23,n,m) = max(fla_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(24,n,m) = max(are_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(25,n,m) = max(aci_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(26,n,m) = peak_systole;

                    end
                end


                if handles.save_selection(6,4)==0; MAT(1,:,:) = NaN; end
                if handles.save_selection(7,4)==0; MAT(2,:,:) = NaN; end
                if handles.save_selection(8,4)==0; MAT(8,:,:) = NaN; end
                if handles.save_selection(9,4)==0; MAT(5,:,:) = NaN; end
                if handles.save_selection(10,4)==0; MAT(9,:,:) = NaN; end
                if handles.save_selection(11,4)==0; MAT(10,:,:) = NaN; end
                if handles.save_selection(12,4)==0; MAT(11,:,:) = NaN; end
                if handles.save_selection(13,4)==0; MAT(12,:,:) = NaN; end
                if handles.save_selection(14,4)==0; MAT(13,:,:) = NaN; end
                if handles.save_selection(17,4)==0; MAT(18,:,:) = NaN; end
                if handles.save_selection(18,4)==0; MAT(17,:,:) = NaN; end
                if handles.save_selection(21,4)==0; MAT(3,:,:) = NaN; end
                if handles.save_selection(22,4)==0; MAT(4,:,:) = NaN; end
                if handles.save_selection(23,4)==0; MAT(14,:,:) = NaN; end
                if handles.save_selection(24,4)==0; MAT(6,:,:) = NaN; end
                if handles.save_selection(25,4)==0; MAT(7,:,:) = NaN; end
                if handles.save_selection(26,4)==0; MAT(15,:,:) = NaN; end
                if handles.save_selection(28,4)==0; MAT(16,:,:) = NaN; end
                if handles.save_selection(29,4)==0; MAT(19,:,:) = NaN; end
                if handles.save_selection(30,4)==0; MAT(20,:,:) = NaN; end
                if handles.save_selection(31,4)==0; MAT(21,:,:) = NaN; end
                if handles.save_selection(32,4)==0; MAT(22,:,:) = NaN; end
                if handles.save_selection(33,4)==0; MAT(23,:,:) = NaN; end
                if handles.save_selection(34,4)==0; MAT(24,:,:) = NaN; end
                if handles.save_selection(35,4)==0; MAT(25,:,:) = NaN; end
    %             if handles.save_selection(36,4)==0; MAT(26,:,:) = NaN; end

                h = waitbar(0,['Saving maximun values (XLS file) cardiac phase ',num2str(1),' of ',num2str(size(handles.veset,3)),' ...']);
                steps = 1;
                st = steps;
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VV = {'Time [s]';...
                      'Flow [ml/s]';...
                      'Net Flow [ml]';...
                      'Maximum Velocity [cm/s]';...
                      'Minimum Velocity [cm/s]'};
                MM = [handles.time';handles.flow';handles.net_flow';handles.max_velocity';handles.min_velocity'];
                TT = table(VV,MM);
                TT.Properties.VariableNames = {'Parameter','Cardiac_Phase'};

                filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MAX_VALUES.xls'];
                writetable(TT,fullfile(filename),'Sheet','2D Flow','Range','A1')

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for n=1:size(handles.veset,3)
                    eval(['Sections_from_1_to_',num2str(handles.number_of_sections),' = MAT(:,:,n);']);
                    T = table(Parameter,eval(['Sections_from_1_to_',num2str(handles.number_of_sections)]));
                    T.Properties.VariableNames = {'Parameter','Section'};

                    filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MAX_VALUES.xls'];
                    writetable(T,fullfile(filename),'Sheet',['Cardiac Phase #',num2str(n)],'Range','A1')

                    if st==n
                        waitbar(st/ size(handles.veset,3),h,['Saving maximun values (XLS file) cardiac phase ',num2str(st),' of ',num2str(size(handles.veset,3)),' ...']);
                        st = st+steps;
                    end
                end
                close(h)

                %%% Minimun Values
                MAT = zeros(26,handles.number_of_sections,size(handles.veset,3));
                for m=1:size(handles.veset,3)
                    for n=1:handles.number_of_sections

                        MAT(1,n,m) = min(mag_velocity(handles.SECTIONS_V{n},m));
                        MAT(2,n,m) = min(mag_wss(handles.SECTIONS_S{n},m));
                        MAT(3,n,m) = min(mag_wssa(handles.SECTIONS_S{n},m));
                        MAT(4,n,m) = min(mag_wssc(handles.SECTIONS_S{n},m));
                        MAT(5,n,m) = min(mag_vor(id_vol_V_n{n},m));
                        MAT(6,n,m) = min(mag_forward_vel(handles.SECTIONS_V{n},m));
                        MAT(7,n,m) = min(mag_backward_vel(handles.SECTIONS_V{n},m));
                        MAT(8,n,m) = min(osi_values(handles.SECTIONS_S{n},m));
                        MAT(9,n,m) = min(hd_values(handles.SECTIONS_V{n},m));
                        MAT(10,n,m) = min(rhd_values(handles.SECTIONS_V{n},m));
                        MAT(11,n,m) = min(vd_values(id_vol_V_n{n},m));
                        MAT(12,n,m) = min(el_values(id_vol_V_n{n},m));
                        MAT(13,n,m) = min(ke_values(handles.SECTIONS_V{n},m));
                        MAT(14,n,m) = min(aa_values(handles.SECTIONS_V{n},m));
                        MAT(15,n,m) = min(rf_values(handles.SECTIONS_S{n},m));
                        MAT(16,n,m) = min(ecc_values(handles.SECTIONS_S{n},m));
                        MAT(17,n,m) = min(dia_values(handles.SECTIONS_S{n},m));
                        MAT(18,n,m) = min(rad_values(handles.SECTIONS_S{n},m));
                        MAT(19,n,m) = min(cur_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(20,n,m) = min(ell_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(21,n,m) = min(len_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
    %                     MAT(22,n,m) = min(cir_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(22,n,m) = min(fov_values(id_vol_V_n{n},m)); % Julio Sotelo 28-05-2019
                        MAT(23,n,m) = min(fla_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(24,n,m) = min(are_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(25,n,m) = min(aci_values(handles.SECTIONS_S{n},m)); % Julio Sotelo 28-05-2019
                        MAT(26,n,m) = peak_systole;

                    end
                end


                if handles.save_selection(6,4)==0; MAT(1,:,:) = NaN; end
                if handles.save_selection(7,4)==0; MAT(2,:,:) = NaN; end
                if handles.save_selection(8,4)==0; MAT(8,:,:) = NaN; end
                if handles.save_selection(9,4)==0; MAT(5,:,:) = NaN; end
                if handles.save_selection(10,4)==0; MAT(9,:,:) = NaN; end
                if handles.save_selection(11,4)==0; MAT(10,:,:) = NaN; end
                if handles.save_selection(12,4)==0; MAT(11,:,:) = NaN; end
                if handles.save_selection(13,4)==0; MAT(12,:,:) = NaN; end
                if handles.save_selection(14,4)==0; MAT(13,:,:) = NaN; end
                if handles.save_selection(17,4)==0; MAT(18,:,:) = NaN; end
                if handles.save_selection(18,4)==0; MAT(17,:,:) = NaN; end
                if handles.save_selection(21,4)==0; MAT(3,:,:) = NaN; end
                if handles.save_selection(22,4)==0; MAT(4,:,:) = NaN; end
                if handles.save_selection(23,4)==0; MAT(14,:,:) = NaN; end
                if handles.save_selection(24,4)==0; MAT(6,:,:) = NaN; end
                if handles.save_selection(25,4)==0; MAT(7,:,:) = NaN; end
                if handles.save_selection(26,4)==0; MAT(15,:,:) = NaN; end
                if handles.save_selection(28,4)==0; MAT(16,:,:) = NaN; end
                if handles.save_selection(29,4)==0; MAT(19,:,:) = NaN; end
                if handles.save_selection(30,4)==0; MAT(20,:,:) = NaN; end
                if handles.save_selection(31,4)==0; MAT(21,:,:) = NaN; end
                if handles.save_selection(32,4)==0; MAT(22,:,:) = NaN; end
                if handles.save_selection(33,4)==0; MAT(23,:,:) = NaN; end
                if handles.save_selection(34,4)==0; MAT(24,:,:) = NaN; end
                if handles.save_selection(35,4)==0; MAT(25,:,:) = NaN; end
    %             if handles.save_selection(36,4)==0; MAT(26,:,:) = NaN; end

                h = waitbar(0,['Saving minimun values (XLS file) cardiac phase ',num2str(1),' of ',num2str(size(handles.veset,3)),' ...']);
                steps = 1;
                st = steps;
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                VV = {'Time [s]';...
                      'Flow [ml/s]';...
                      'Net Flow [ml]';...
                      'Maximum Velocity [cm/s]';...
                      'Minimum Velocity [cm/s]'};
                MM = [handles.time';handles.flow';handles.net_flow';handles.max_velocity';handles.min_velocity'];
                TT = table(VV,MM);
                TT.Properties.VariableNames = {'Parameter','Cardiac_Phase'};

                filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MIN_VALUES.xls'];
                writetable(TT,fullfile(filename),'Sheet','2D Flow','Range','A1')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for n=1:size(handles.veset,3)
                    eval(['Sections_from_1_to_',num2str(handles.number_of_sections),' = MAT(:,:,n);']);
                    T = table(Parameter,eval(['Sections_from_1_to_',num2str(handles.number_of_sections)]));
                    T.Properties.VariableNames = {'Parameter','Section'};

                    filename = [directory,'/XLS FILES/RESULTS_SECTIONS_MIN_VALUES.xls'];
                    writetable(T,fullfile(filename),'Sheet',['Cardiac Phase #',num2str(n)],'Range','A1')

                    if st==n
                        waitbar(st/ size(handles.veset,3),h,['Saving minimun values (XLS file) cardiac phase ',num2str(st),' of ',num2str(size(handles.veset,3)),' ...']);
                        st = st+steps;
                    end
                end
                close(h)

            end
            disp('The data save process is finish ...')
        else
            msgbox('The data was not saved ...','Warning','warn')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure1_SizeChangedFcn(hObject, eventdata, handles)

% change size windows
set(handles.figure1, 'Units', 'pixels');
FigPos = get(handles.figure1, 'Position');
set(handles.figure1, 'Position', [FigPos(1:2),FigPos(4),FigPos(4) ]);
set(handles.figure1, 'Units', 'normalized');

% change size zoom buttons
set(handles.pushbutton2, 'Units', 'pixels');
handles.pushbutton2_size = get(handles.pushbutton2, 'Position');
set(handles.pushbutton2, 'Units', 'normalized');
idx = mod(min(handles.pushbutton2_size(3:4)),2)>1;
w = floor(min(handles.pushbutton2_size(3:4)));
w(idx) = w(idx)+1;
im = imread('Symbols/P3.png');
g = double(imresize(im,[w-4 w-4],'method','nearest')>0);
set(handles.pushbutton2,'CData',g,'visible','on')
set(handles.pushbutton3,'CData',g,'visible','on')
set(handles.pushbutton4,'CData',g,'visible','on')
set(handles.pushbutton5,'CData',g)

% 4D FLOW
set(handles.text1,'FontUnits','Normalized','FontSize',0.59)
set(handles.text2,'FontUnits','Normalized','FontSize',0.59)
set(handles.text3,'FontUnits','Normalized','FontSize',0.59)
set(handles.text4,'FontUnits','Normalized','FontSize',0.59)
set(handles.text99,'FontUnits','Normalized','FontSize',0.58)
set(handles.text101,'FontUnits','Normalized','FontSize',0.58)
set(handles.text102,'FontUnits','Normalized','FontSize',0.58)
set(handles.text111,'FontUnits','Normalized','FontSize',0.58)
set(handles.text112,'FontUnits','Normalized','FontSize',0.58)
set(handles.text100,'FontUnits','Normalized','FontSize',0.58)
set(handles.text103,'FontUnits','Normalized','FontSize',0.58)
set(handles.text104,'FontUnits','Normalized','FontSize',0.58)
set(handles.text109,'FontUnits','Normalized','FontSize',0.58)
set(handles.text110,'FontUnits','Normalized','FontSize',0.58)
set(handles.text98,'FontUnits','Normalized','FontSize',0.58)
set(handles.text105,'FontUnits','Normalized','FontSize',0.58)
set(handles.text106,'FontUnits','Normalized','FontSize',0.58)
set(handles.text107,'FontUnits','Normalized','FontSize',0.58)
set(handles.text108,'FontUnits','Normalized','FontSize',0.58)
set(handles.pushbutton9,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton10,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton11,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton12,'FontUnits','Normalized','FontSize',0.22)
set(handles.pushbutton13,'FontUnits','Normalized','FontSize',0.22)
set(handles.pushbutton14,'FontUnits','Normalized','FontSize',0.51)
set(handles.pushbutton62,'FontUnits','Normalized','FontSize',0.22)
set(handles.pushbutton69,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton70,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton71,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton78,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton79,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton80,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton81,'FontUnits','Normalized','FontSize',0.22)
set(handles.pushbutton83,'FontUnits','Normalized','FontSize',0.22)
set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.59)

% UNWRAPPING
set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.47)
set(handles.pushbutton15,'FontUnits','Normalized','FontSize',0.47)
set(handles.pushbutton16,'FontUnits','Normalized','FontSize',0.29)
set(handles.pushbutton17,'FontUnits','Normalized','FontSize',0.29)
set(handles.pushbutton18,'FontUnits','Normalized','FontSize',0.29)
set(handles.pushbutton19,'FontUnits','Normalized','FontSize',0.29)
set(handles.pushbutton20,'FontUnits','Normalized','FontSize',0.29)
set(handles.pushbutton21,'FontUnits','Normalized','FontSize',0.29)
set(handles.text5,'FontUnits','Normalized','FontSize',0.64)
set(handles.text6,'FontUnits','Normalized','FontSize',0.64)
set(handles.text7,'FontUnits','Normalized','FontSize',0.59)
set(handles.text8,'FontUnits','Normalized','FontSize',0.59)
set(handles.text9,'FontUnits','Normalized','FontSize',0.59)
set(handles.text10,'FontUnits','Normalized','FontSize',0.59)
set(handles.text11,'FontUnits','Normalized','FontSize',0.59)
set(handles.text12,'FontUnits','Normalized','FontSize',0.59)
set(handles.edit1,'FontUnits','Normalized','FontSize',0.59)
set(handles.edit2,'FontUnits','Normalized','FontSize',0.59)
set(handles.popupmenu2,'FontUnits','Normalized','FontSize',0.53)
set(handles.popupmenu3,'FontUnits','Normalized','FontSize',0.53)

handles.output = hObject;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
    handles.input1 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
    handles.input2 = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
    id_while = 0;
    while(id_while == 0)
        GUIDE_CONTRAST(handles.input1, handles.input2);
        handles.input1 = getappdata(0,'OUT');
        id_while = getappdata(0,'closed_loop');
        handles.IPCMRA = handles.input1;
        handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;

        % Julio Sotelo 23-11-2018
        if handles.id_ang == 1% Julio Sotelo 23-11-2018
            handles.ANG = handles.IPCMRA;% Julio Sotelo 23-11-2018
        elseif handles.id_mag == 1% Julio Sotelo 23-11-2018
            handles.MAG = handles.IPCMRA;% Julio Sotelo 23-11-2018
        end% Julio Sotelo 23-11-2018

        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool2_ClickedCallback(hObject, eventdata, handles)

if handles.id_seg == 1 
    
    answer = questdlg({'If you change the size of the image the segmentation will be removed';' Do you want to continue? ...'},'Warning','Yes','No','No');
    switch answer
        case 'Yes'
           set(handles.axes4,'Color',[0 0 0])
            axes(handles.axes4);
            colorbar('off')
            plot(0.0)
            axis off
            set(handles.text4, 'visible', 'off');
            list_string = {'...'};
            set(handles.popupmenu1,'visible','off','String',list_string,'value',1);
            set(handles.pushbutton5, 'visible', 'off');
            set(handles.pushbutton12, 'visible', 'off');
            set(handles.pushbutton13, 'visible', 'off');
            set(handles.pushbutton14, 'visible', 'off');
            set(handles.slider4, 'visible', 'off');
            handles.SEG = zeros(size(handles.SEG));
            handles.Lrgb = zeros(size(handles.Lrgb));
            handles.L = zeros(size(handles.L));
            handles.id_unwrappping = 0;
            set(handles.pushbutton6, 'Units', 'pixels');
            handles.pushbutton6_size = get(handles.pushbutton6, 'Position');
            set(handles.pushbutton6, 'Units', 'normalized');
            idx = mod(min(handles.pushbutton6_size(3:4)),2)>1;
            w = floor(min(handles.pushbutton6_size(3:4)));
            w(idx) = w(idx)+1;
            im = imread('Symbols/P2.png');
            g = double(imresize(im,[w-4 w-4],'method','nearest')>0);
            set(handles.pushbutton6,'CData',g,'visible','on')
            set(handles.pushbutton9,'visible','on')
            set(handles.pushbutton7,'CData',g,'visible','on')
            set(handles.pushbutton10,'visible','on')
            set(handles.pushbutton8,'CData',g,'visible','on')
            set(handles.pushbutton11,'visible','on')
            handles.id_resizing = 1;
            handles.c1 = handles.xd(1,1);
            handles.c2 = handles.yd(1,1);
            handles.c3 = handles.zd(1,1);
            handles.c4 = handles.xd(end,end);
            handles.c5 = handles.yd(end,end);
            handles.c6 = handles.zd(end,end);
            handles.pos1 = [handles.c2 handles.c1 handles.c5 handles.c4];
            handles.pos2 = [handles.c3 handles.c1 handles.c6 handles.c4];
            handles.pos3 = [handles.c3 handles.c2 handles.c6 handles.c5];
            [a,b,c,d] = size(handles.IPCMRA);
            set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
            set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
            set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
            if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
            if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
            if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
            if handles.id_resizing == 1
                rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axes(handles.axes2);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
            if handles.id_resizing == 1
                rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            axes(handles.axes3);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
            if handles.id_resizing == 1
                rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        case 'No'
            
    end
else 
    set(handles.axes4,'Color',[0 0 0])
    axes(handles.axes4);
    colorbar('off')
    plot(0.0)
    axis off
    set(handles.text4, 'visible', 'off');
    list_string = {'...'};
    set(handles.popupmenu1,'visible','off','String',list_string,'value',1);
    set(handles.pushbutton5, 'visible', 'off');
    set(handles.pushbutton12, 'visible', 'off');
    set(handles.pushbutton13, 'visible', 'off');
    set(handles.pushbutton14, 'visible', 'off');
    set(handles.slider4, 'visible', 'off');
    handles.SEG = zeros(size(handles.SEG));
    handles.Lrgb = zeros(size(handles.Lrgb));
    handles.L = zeros(size(handles.L));
    handles.id_unwrappping = 0;
    set(handles.pushbutton6, 'Units', 'pixels');
    handles.pushbutton6_size = get(handles.pushbutton6, 'Position');
    set(handles.pushbutton6, 'Units', 'normalized');
    idx = mod(min(handles.pushbutton6_size(3:4)),2)>1;
    w = floor(min(handles.pushbutton6_size(3:4)));
    w(idx) = w(idx)+1;
    im = imread('Symbols/P2.png');
    g = double(imresize(im,[w-4 w-4],'method','nearest')>0);
    set(handles.pushbutton6,'CData',g,'visible','on')
    set(handles.pushbutton9,'visible','on')
    set(handles.pushbutton7,'CData',g,'visible','on')
    set(handles.pushbutton10,'visible','on')
    set(handles.pushbutton8,'CData',g,'visible','on')
    set(handles.pushbutton11,'visible','on')
    handles.id_resizing = 1;
    handles.c1 = handles.xd(1,1);
    handles.c2 = handles.yd(1,1);
    handles.c3 = handles.zd(1,1);
    handles.c4 = handles.xd(end,end);
    handles.c5 = handles.yd(end,end);
    handles.c6 = handles.zd(end,end);
    handles.pos1 = [handles.c2 handles.c1 handles.c5 handles.c4];
    handles.pos2 = [handles.c3 handles.c1 handles.c6 handles.c4];
    handles.pos3 = [handles.c3 handles.c2 handles.c6 handles.c5];
    [a,b,c,d] = size(handles.IPCMRA);
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,1)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,3)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool3_ClickedCallback(hObject, eventdata, handles)
    set(handles.axes4,'Color',[0 0 0])
    axes(handles.axes4);
    colorbar('off')
    plot(0.0)
    axis off
    set(handles.text4, 'visible', 'off');
    list_string = {'...'};
    set(handles.popupmenu1,'visible','off','String',list_string,'value',1);
    set(handles.pushbutton5, 'visible', 'off');
    set(handles.pushbutton12, 'visible', 'off');
    set(handles.pushbutton13, 'visible', 'off');
    set(handles.pushbutton14, 'visible', 'off');
    set(handles.slider4, 'visible', 'off');
    handles.id_unwrappping = 0;
    input.id            = 1;
    input.IPCMRA        = (handles.IPCMRA./max(handles.IPCMRA(:)))*255;
    input.SEG           = handles.SEG;
    input.L             = handles.L;
    input.NUM           = handles.NUM;
    input.min_value_th  = min(input.IPCMRA(:));
    input.max_value_th  = max(input.IPCMRA(:));
    input.xd            = handles.xd;
    input.yd            = handles.yd;
    input.zd            = handles.zd;
    input.a             = handles.a;
    input.b             = handles.b;
    input.c             = handles.c;
    input.slider_axes1  = handles.slider_axes1;
    input.slider_axes2  = handles.slider_axes2;
    input.slider_axes3  = handles.slider_axes3;
    
    
    GUIDE_THRESHOLDING(input)
    handles.SEG = getappdata(0,'SEG');
    handles.L = getappdata(0,'L');
    handles.NUM = getappdata(0,'NUM');
    handles.Lrgb = getappdata(0,'Lrgb');

    handles.min_value_th = getappdata(0,'min_value_th');
    handles.max_value_th = getappdata(0,'max_value_th');
    if max(max(max(handles.NUM)))~=0
        handles.id_seg = 1;
        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap(handles.axes1,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap(handles.axes2,'gray')
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap(handles.axes3,'gray')
        axis off
        daspect([1 1 1])
        set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
        set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
        set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
        list_string = {'...','Surface','Voxel'};
        set(handles.popupmenu1,'visible','on','String',list_string);
        popupmnenu1_pos = get(handles.popupmenu1,'Value');
        if popupmnenu1_pos == 2
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
            xd = X*handles.voxel_MR(1);
            yd = Y*handles.voxel_MR(2);
            zd = Z*handles.voxel_MR(3);
            xd = permute(xd,[2 1 3]);
            yd = permute(yd,[2 1 3]);
            zd = permute(zd,[2 1 3]);
            id_L = unique(sort(handles.L(:)));
            id_L(id_L==0)=[];
            axes(handles.axes4);
            plot(0,0)
            axis off
            axes(handles.axes4);
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
        elseif popupmnenu1_pos == 3
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
            axes(handles.axes4);
            plot(0,0)
            axis off
            axes(handles.axes4);
            patch('Vertices', node_bin_aorta, 'Faces', elem_bin_aorta,'FaceColor',[0.85 0.85 0.85],'EdgeColor','k');
            lighting gouraud
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        end
    end
    handles.save_id_SEG_mat = 1;
    handles.save_id_SEG_vti = 1;
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool4_ClickedCallback(hObject, eventdata, handles)
% This function excecute the manual segmentation process


    if sum(handles.SEG(:))>0
        
        handles.id_unwrappping = 0;
        input.id            = 1;
        input.IPCMRA        = handles.IPCMRA;
        input.ANG           = handles.ANG; % Julio Sotelo
        input.MAG           = handles.MAG; % Julio Sotelo
        input.id_ang        = handles.id_ang;% Julio Sotelo
        input.id_mag        = handles.id_mag;% Julio Sotelo
        
        if handles.id_ang == 1
            input.idangmag      = 1;
        else 
            input.idangmag      = 2;
        end
        
        input.SEG           = handles.SEG;
        input.L             = handles.L;
        input.Lrgb          = handles.Lrgb;
        input.NUM           = handles.NUM;
        input.SEG_old       = handles.SEG;
        input.Lrgb_old      = handles.Lrgb;
        input.L_old         = handles.L ;
        input.NUM_old       = handles.NUM;
        input.min_value_th  = handles.min_value_th;
        input.max_value_th  = handles.max_value_th;
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
        input.view_sac      = 1;
        input.id_seg        = handles.id_seg;
        input.id_vel        = handles.id_vel;
        input.id_vor        = handles.id_vor;
        input.id_hd         = handles.id_hd;
        input.id_rhd        = handles.id_rhd;
        input.id_vd         = handles.id_vd;
        input.id_el         = handles.id_el;
        input.id_ke         = handles.id_ke;
        input.Lrgb_vel      = handles.Lrgb_vel;
        input.Lrgb_vor      = handles.Lrgb_vor;
        input.Lrgb_hd       = handles.Lrgb_hd;
        input.Lrgb_rhd      = handles.Lrgb_rhd;
        input.Lrgb_vd       = handles.Lrgb_vd;
        input.Lrgb_el       = handles.Lrgb_el;
        input.Lrgb_ke       = handles.Lrgb_ke;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        input.id_fve        = handles.id_fve;
        input.Lrgb_fve      = handles.Lrgb_fve;
        input.id_bve        = handles.id_bve;
        input.Lrgb_bve      = handles.Lrgb_bve;
        input.id_aan        = handles.id_aan;
        input.Lrgb_aan      = handles.Lrgb_aan;
        input.id_fov        = handles.id_fov;
        input.Lrgb_fov      = handles.Lrgb_fov;
        
        
        input.peak_flow      = handles.peak_flow;
        
        SEG_Original = handles.SEG;
        id_while = 0;
        while(1)
            while(id_while == 0)
                GUIDE_SEGMENTATION(input)
                input.SEG = getappdata(0,'SEG');
                input.idangmag = getappdata(0,'idangmag');
                input.IPCMRA = getappdata(0,'IPCMRA');
                input.L = getappdata(0,'L');
                input.NUM = getappdata(0,'NUM');
                input.Lrgb = getappdata(0,'Lrgb');
                input.slider_axes1  = getappdata(0,'slider_axes1');
                input.slider_axes2  = getappdata(0,'slider_axes2');
                input.slider_axes3  = getappdata(0,'slider_axes3');
                input.view_sac  = getappdata(0,'view_sac');
                id_while = getappdata(0,'id_while');
            end
            handles.SEG = input.SEG;
            handles.L = input.L;
            handles.NUM = input.NUM;
            handles.Lrgb = input.Lrgb;
            segdif = handles.SEG-SEG_Original;
            if sum(segdif(:))==0
                answer = questdlg(  'Do you sure that the segmentation is OK ...','Warning','Yes','No','No');
                switch answer
                    case 'Yes'
                        set(handles.pushbutton12,'visible','on');
                        break
                    case 'No'
                        id_while = 0;
                end
            else
                set(handles.pushbutton12,'visible','on');
                break
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if input.idangmag == 1
            handles.id_ang      = 1;
            handles.id_mag      = 0;
            handles.IPCMRA      = handles.ANG;
        else 
            handles.id_ang      = 0;
            handles.id_mag      = 1;
            handles.IPCMRA      = handles.MAG;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % showing the information in the 4D flow app

        if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
        if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
        if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        hold on
        plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axes(handles.axes2);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        axes(handles.axes3);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
        hold on
        plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
        plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
        if handles.id_seg == 1 && sum(handles.L(:))~=0
            Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
        end
        if handles.id_resizing == 1
            rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        popupmnenu1_pos = get(handles.popupmenu1,'Value');
        if popupmnenu1_pos == 2
            [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
            xd = X*handles.voxel_MR(1);
            yd = Y*handles.voxel_MR(2);
            zd = Z*handles.voxel_MR(3);
            xd = permute(xd,[2 1 3]);
            yd = permute(yd,[2 1 3]);
            zd = permute(zd,[2 1 3]);
            id_L = unique(sort(handles.L(:)));
            id_L(id_L==0)=[];
            axes(handles.axes4);
            plot(0,0)
            axis off
            axes(handles.axes4);
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
        elseif popupmnenu1_pos == 3
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
            axes(handles.axes4);
            plot(0,0)
            axis off
            axes(handles.axes4);
            patch('Vertices', node_bin_aorta, 'Faces', elem_bin_aorta,'FaceColor',[0.85 0.85 0.85],'EdgeColor','k');
            lighting gouraud
            axis vis3d
            daspect([1,1,1])
            axis off
            view([-34,-51])
        end
    else
        msgbox('The segmentation is necessary to run this option ...','Warning','warn')
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu2_Callback(hObject, eventdata, handles)

switch get(handles.popupmenu2,'Value')
    
    case 1
        handles.slider_axes_uwrap_1 = round(handles.c/2);
        
        axes(handles.axes5);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
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
        set(handles.slider6,'Value', handles.slider_axes_uwrap_1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.c)])
        
    case 2
        handles.slider_axes_uwrap_1 = round(handles.b/2);
        
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
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
        set(handles.slider6,'Value', handles.slider_axes_uwrap_1/handles.b,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.b)])
        
    case 3
        handles.slider_axes_uwrap_1 = round(handles.a/2);
        
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
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
        set(handles.slider6,'Value', handles.slider_axes_uwrap_1/handles.a,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.a)])
end

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu3_Callback(hObject, eventdata, handles)

switch get(handles.popupmenu3,'Value')
    case 1
        
        handles.FLOW_IM = handles.PHASE_FH;
        handles.LAP = handles.LAP_FH;
        handles.flow_enc = 1;
        handles.LAP_SEG_positive = handles.LAP.*double(handles.LAP>0); % defauls values for the slider to adjust the threshold
        handles.LAP_SEG_negative = handles.LAP.*double(handles.LAP<0); % defauls values for the slider to adjust the threshold
        
        if handles.view_sac == 1
            axes(handles.axes5);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
         elseif handles.view_sac == 2
             
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
        elseif handles.view_sac == 3
            
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        end
        
    case 2
        
        handles.FLOW_IM = handles.PHASE_AP;
        handles.LAP = handles.LAP_AP;
        handles.flow_enc = 2;
        handles.LAP_SEG_positive = handles.LAP.*double(handles.LAP>0); % defauls values for the slider to adjust the threshold
        handles.LAP_SEG_negative = handles.LAP.*double(handles.LAP<0); % defauls values for the slider to adjust the threshold
        
        if handles.view_sac == 1
            
            axes(handles.axes5);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
         elseif handles.view_sac == 2
             
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
        elseif handles.view_sac == 3
            
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        end
    case 3
        
        handles.FLOW_IM = handles.PHASE_RL;
        handles.LAP = handles.LAP_RL;
        handles.flow_enc = 3;
        handles.LAP_SEG_positive = handles.LAP.*double(handles.LAP>0); % defauls values for the slider to adjust the threshold
        handles.LAP_SEG_negative = handles.LAP.*double(handles.LAP<0); % defauls values for the slider to adjust the threshold
        
        if handles.view_sac == 1
            axes(handles.axes5);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
         elseif handles.view_sac == 2
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            
        elseif handles.view_sac == 3
            axes(handles.axes5);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
            hold on
            if sum(handles.SEG_positive(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            if sum(handles.SEG_negative(:)) > 0
                Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
                set(himage, 'AlphaData', cdata);
            end
            hold off
            axis image
            colormap gray
            axis off
            daspect([1 1 1])
        end
end

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider5_Callback(hObject, eventdata, handles)

    pp=1/handles.d;
    slider_step(1) = pp;
    slider_step(2) = 0.1;
    set(handles.slider5,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.d;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
            handles.slider_value = 1;
    else
            handles.slider_value = round(handles.slider_value*maxslice);
    end
    handles.slider_axes_uwrap_2 = handles.slider_value;
    set(handles.text5,'String',['Cardiac Phase # ',num2str(handles.slider_axes_uwrap_2),' of ',num2str(handles.d)])
    if handles.view_sac == 1
        axes(handles.axes5);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
     elseif handles.view_sac == 2
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    elseif handles.view_sac == 3
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider5_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider6_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1;
        pp=1/handles.c;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider5,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = handles.c;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.slider_axes_uwrap_1 = handles.slider_value;
        axes(handles.axes5);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.c)])
    elseif handles.view_sac == 2;
        pp=1/handles.b;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider6,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = handles.b;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.slider_axes_uwrap_1 = handles.slider_value;
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.c)])
    elseif handles.view_sac == 3;
        pp=1/handles.a;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider6,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = handles.a;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.slider_axes_uwrap_1 = handles.slider_value;
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        set(handles.text6,'String',['Slice # ',num2str(handles.slider_axes_uwrap_1),' of ',num2str(handles.c)])
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider6_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider7_Callback(hObject, eventdata, handles)
    
    % move the slider to set the SEG_positive
    slider_step(1) = 0.05;
    slider_step(2) = 0.1;
    set(handles.slider7,'sliderstep',slider_step,'max',1,'min',0)
    handles.slider_value = get(hObject,'Value');
    set(handles.text8,'String',['Value : ',num2str(handles.slider_value)])

    handles.slider_value_positive = handles.slider_value;
    handles.SEG_positive = double(handles.LAP_SEG_positive>handles.slider_value_positive);
    
    green = repmat(zeros(size(handles.SEG_positive)),[1,1,1,1,3]);
    green(:,:,:,:,2) = green(:,:,:,:,2)+1;
    handles.LRGB_SEG_positive = repmat(handles.SEG_positive,[1,1,1,1,3]).*green;
    
    if handles.view_sac == 1
        axes(handles.axes5);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        
     elseif handles.view_sac == 2
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        
    elseif handles.view_sac == 3
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider7_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider8_Callback(hObject, eventdata, handles)
    slider_step(1) = 0.05;
    slider_step(2) = 0.1;
    set(handles.slider8,'sliderstep',slider_step,'max',1,'min',0)
    handles.slider_value = get(hObject,'Value');
    set(handles.text10,'String',['Value : ',num2str(handles.slider_value)])
    
    handles.slider_value_negative = handles.slider_value*-1;
    handles.SEG_negative = double(handles.LAP_SEG_negative<handles.slider_value_negative);
    
    red = repmat(zeros(size(handles.SEG_negative)),[1,1,1,1,3]);
    red(:,:,:,:,1) = red(:,:,:,:,1)+1;
    handles.LRGB_SEG_negative = repmat(handles.SEG_negative,[1,1,1,1,3]).*red;
    
    if handles.view_sac == 1
        axes(handles.axes5);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        
     elseif handles.view_sac == 2
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        
    elseif handles.view_sac == 3
        axes(handles.axes5);
        imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
        hold on
        if sum(handles.SEG_positive(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        if sum(handles.SEG_negative(:)) > 0
            Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
            himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
            cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
            set(himage, 'AlphaData', cdata);
        end
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider8_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_Callback(hObject, eventdata, handles)
handles.PI_W_Positive = str2double(get(hObject,'String'));
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit2_Callback(hObject, eventdata, handles)
handles.PI_W_Negative = str2double(get(hObject,'String'));
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton16_Callback(hObject, eventdata, handles)
handles.FLOW_IM_temp = handles.FLOW_IM;
handles.FLOW_IM_temp(:,:,:,handles.slider_axes_uwrap_2) = handles.FLOW_IM_temp(:,:,:,handles.slider_axes_uwrap_2) + handles.SEG_positive(:,:,:,handles.slider_axes_uwrap_2)*handles.PI_W_Positive*pi;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM_temp(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton17_Callback(hObject, eventdata, handles)
handles.FLOW_IM_temp = handles.FLOW_IM + handles.SEG_positive*handles.PI_W_Positive*pi;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM_temp(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton18_Callback(hObject, eventdata, handles)
if handles.flow_enc == 1
    handles.PHASE_FH = handles.FLOW_IM_temp;
    handles.MR_PCA_FH = (handles.PHASE_FH*handles.VENC)/pi;
    handles.slider_value_positive = 0;
    handles.SEG_positive = zeros(size(handles.LAP));
    handles.LRGB_SEG_positive = repmat(ones(size(handles.SEG_positive)),[1 1 1 1 3]);
elseif handles.flow_enc == 2
    handles.PHASE_AP = handles.FLOW_IM_temp;
    handles.MR_PCA_AP = (handles.PHASE_AP*handles.VENC)/pi;
    handles.slider_value_positive = 0;
    handles.SEG_positive = zeros(size(handles.LAP));
    handles.LRGB_SEG_positive = repmat(ones(size(handles.SEG_positive)),[1 1 1 1 3]);
elseif handles.flow_enc == 3
    handles.PHASE_RL = handles.FLOW_IM_temp;
    handles.MR_PCA_RL = (handles.PHASE_RL*handles.VENC)/pi;
    handles.slider_value_positive = 0;
    handles.SEG_positive = zeros(size(handles.LAP));
    handles.LRGB_SEG_positive = repmat(ones(size(handles.SEG_positive)),[1 1 1 1 3]);
end
handles.FLOW_IM = handles.FLOW_IM_temp;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_negative(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_negative(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
slider_step(1) = 0.05;
slider_step(2) = 0.1;
set(handles.slider7,'Value', 0,'sliderstep',slider_step,'max',1,'min',0)
set(handles.text8,'String',['Value : ',num2str(0)])

handles.IPCMRA = (1/size(handles.MR_FFE_FH,4))*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
handles.ANG = handles.IPCMRA;

waitfor(msgbox('The data has been saved ... ...'))
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton19_Callback(hObject, eventdata, handles)
handles.FLOW_IM_temp = handles.FLOW_IM;
handles.FLOW_IM_temp(:,:,:,handles.slider_axes_uwrap_2) = handles.FLOW_IM_temp(:,:,:,handles.slider_axes_uwrap_2) - handles.SEG_negative(:,:,:,handles.slider_axes_uwrap_2)*abs(handles.PI_W_Negative)*pi;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM_temp(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton20_Callback(hObject, eventdata, handles)
handles.FLOW_IM_temp = handles.FLOW_IM - handles.SEG_negative*abs(handles.PI_W_Negative)*pi;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM_temp(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM_temp(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton21_Callback(hObject, eventdata, handles)
if handles.flow_enc == 1
    handles.PHASE_FH = handles.FLOW_IM_temp;
    handles.MR_PCA_FH = (handles.PHASE_FH*handles.VENC)/pi;
    handles.slider_value_negative = 0;
    handles.SEG_negative = zeros(size(handles.LAP));
    handles.LRGB_SEG_negative = repmat(ones(size(handles.SEG_negative)),[1 1 1 1 3]);
elseif handles.flow_enc == 2
    handles.PHASE_AP = handles.FLOW_IM_temp;
    handles.MR_PCA_AP = (handles.PHASE_AP*handles.VENC)/pi;
    handles.slider_value_negative = 0;
    handles.SEG_negative = zeros(size(handles.LAP));
    handles.LRGB_SEG_negative = repmat(ones(size(handles.SEG_negative)),[1 1 1 1 3]);
elseif handles.flow_enc == 3
    handles.PHASE_RL = handles.FLOW_IM_temp;
    handles.MR_PCA_RL = (handles.PHASE_RL*handles.VENC)/pi;
    handles.slider_value_negative = 0;
    handles.SEG_negative = zeros(size(handles.LAP));
    handles.LRGB_SEG_negative = repmat(ones(size(handles.SEG_negative)),[1 1 1 1 3]);
end
handles.FLOW_IM = handles.FLOW_IM_temp;
if handles.view_sac == 1
    axes(handles.axes5);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2)));
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,:,handles.slider_axes_uwrap_1,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
 elseif handles.view_sac == 2
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.FLOW_IM(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(:,handles.slider_axes_uwrap_1,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
elseif handles.view_sac == 3
    axes(handles.axes5);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.FLOW_IM(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2)))
    hold on
    if sum(handles.SEG_positive(:)) > 0
        Lrgb2d = squeeze(handles.LRGB_SEG_positive(handles.slider_axes_uwrap_1,:,:,handles.slider_axes_uwrap_2,:));
        himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
        cdata = ((double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))==1)*0.8;
        set(himage, 'AlphaData', cdata);
    end
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
end
slider_step(1) = 0.05;
slider_step(2) = 0.1;
set(handles.slider8,'Value', 0,'sliderstep',slider_step,'max',1,'min',0)
set(handles.text10,'String',['Value : ',num2str(0)])

handles.IPCMRA = (1/size(handles.MR_FFE_FH,4))*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;
handles.ANG = handles.IPCMRA;

waitfor(msgbox('The data has been saved ...'))
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Project_Callback(hObject, eventdata, handles)
    answer = questdlg(  {'This process can take a few second';'Do you want to continue? ...'},'Question','Yes','No','No');
    switch answer
        case 'Yes'
            [file,path] = uigetfile(pwd,'Select File Project');
            if file~=0
                waitfor(msgbox('This process can take a few seconds ...'));
                load([path,file])
                closereq
            end
        case 'No'
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Project_Callback(hObject, eventdata, handles)
    directory = uigetdir(pwd, 'Select Directory');
    if directory~=0
        c = msgbox('Saving Data...');
        save([directory,'/handles.mat'],'handles','-v7.3')
        close(c)
    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_SEG_Callback(hObject, eventdata, handles)
% hObject    handle to Load_SEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [file,path] = uigetfile('*.mat','Select Segmentation File');
    
    if file~=0
        load([path,file])
        name_file = file(1:end-4);        
        [as,bs,cs] = size(eval(name_file));
        
        
        if handles.a == as && handles.b == bs && handles.c == cs
            
            handles.SEG = eval(name_file);
            handles.Lrgb = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
            handles.Lrgb(:,:,:,1) =  handles.Lrgb(:,:,:,1)-handles.SEG + handles.SEG*0.1667;
            handles.Lrgb(:,:,:,3) =  handles.Lrgb(:,:,:,3)-handles.SEG + handles.SEG*0.1667;
            handles.NUM = 1;
            handles.L = handles.SEG;
            
            handles.id_seg = 1;
            handles.save_id_SEG_mat = 1;
            handles.save_id_SEG_vti = 1;
            
            if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
            if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
            if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes2);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes2,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes3);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes3,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
            set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
            set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
            list_string = {'...','Surface','Voxel'};
            set(handles.popupmenu1,'visible','on','String',list_string);
            
        elseif handles.a-2 == as && handles.b-2 == bs && handles.c-2 == cs % load original without modification
            
            handles.SEG(2:end-1,2:end-1,2:end-1) = eval(name_file);
            handles.Lrgb = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
            handles.Lrgb(:,:,:,1) =  handles.Lrgb(:,:,:,1)-handles.SEG + handles.SEG*0.1667;
            handles.Lrgb(:,:,:,3) =  handles.Lrgb(:,:,:,3)-handles.SEG + handles.SEG*0.1667;
            handles.NUM = 1;
            handles.L = handles.SEG;
            
            handles.id_seg = 1;
            handles.save_id_SEG_mat = 1;
            handles.save_id_SEG_vti = 1;
            
            if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
            if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
            if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes2);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes2,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes3);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes3,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
            set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
            set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
            list_string = {'...','Surface','Voxel'};
            set(handles.popupmenu1,'visible','on','String',list_string);
            
        else
            
            load([path,'Coor.mat'])
            
            f1 = position_sag_cor(1);
            f2 = position_sag_cor(1) + position_sag_cor(3)-1;
            f3 = position_sag_cor(2);
            f4 = position_sag_cor(2) + position_sag_cor(4)-1;
            f5 = position_axi_cor(1);
            f6 = position_axi_cor(1) + position_axi_cor(3)-1;
            
            handles.SEG(f3:f4,f1:f2,f5:f6) = eval(name_file);
            handles.Lrgb = ones(size(handles.SEG,1),size(handles.SEG,2),size(handles.SEG,3),3);
            handles.Lrgb(:,:,:,1) =  handles.Lrgb(:,:,:,1)-handles.SEG + handles.SEG*0.1667;
            handles.Lrgb(:,:,:,3) =  handles.Lrgb(:,:,:,3)-handles.SEG + handles.SEG*0.1667;
            handles.NUM = 1;
            handles.L = handles.SEG;
            
            handles.id_seg = 1;
            handles.save_id_SEG_mat = 1;
            handles.save_id_SEG_vti = 1;
            
            if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
            if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
            if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
            axes(handles.axes1);
            imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
            hold on
            plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes1,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes2);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes2,'gray')
            axis off
            daspect([1 1 1])
            axes(handles.axes3);
            imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
            hold on
            plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
            plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
            if handles.id_seg == 1 && sum(handles.L(:))~=0
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
            end
            if handles.id_resizing == 1
                rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
            end
            hold off
            axis image
            colormap(handles.axes3,'gray')
            axis off
            daspect([1 1 1])
            set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
            set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
            set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])
            list_string = {'...','Surface','Voxel'};
            set(handles.popupmenu1,'visible','on','String',list_string);

        end

    end
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool7_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.IPCMRA = (1/handles.d)*sum( (handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
    handles.IPCMRA = (handles.IPCMRA/max(handles.IPCMRA(:)))*255;

    handles.ANG = handles.IPCMRA;
    handles.id_ang = 1;
    handles.id_mag = 0;

    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MAG = mean(handles.MR_FFE_AP,4);
    handles.IPCMRA = handles.MAG;
    handles.id_ang = 0;
    handles.id_mag = 1;

    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos1,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes1,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos2,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes2,'gray')
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
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
    if handles.id_resizing == 1
        rectangle('Position',handles.pos3,'LineWidth',1,'EdgeColor','y','LineStyle','-')
    end
    hold off
    axis image
    colormap(handles.axes3,'gray')
    axis off
    daspect([1 1 1])
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(handles.slider_axes1)),' / ',num2str(handles.c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes2)),' / ',num2str(handles.b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(handles.slider_axes3)),' / ',num2str(handles.a),' Type: ',handles.type])

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton62_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    input.SEG = handles.SEG;
    input.IPCMRA = handles.IPCMRA;
    input.voxel_MR = handles.voxel_MR;
    input.L = handles.L;
    input.Lrgb = handles.Lrgb;
    input.Lrgb_vel = handles.Lrgb_vel;
    input.NUM = handles.NUM;
    input.xd = handles.xd;
    input.yd = handles.yd;
    input.zd = handles.zd;
    input.a = handles.a;
    input.b = handles.b;
    input.c = handles.c;
    input.d = handles.d;
    input.slider_id_axes1 = handles.slider_axes1;
    input.SEG = handles.SEG;
    input.MR_FFE_FH = handles.MR_FFE_FH;
    input.MR_PCA_AP = handles.MR_PCA_AP_smooth;
    input.MR_PCA_FH = handles.MR_PCA_FH_smooth;
    input.MR_PCA_RL = handles.MR_PCA_RL_smooth;
    input.voxel_MR = handles.voxel_MR;
    input.faces = handles.faces;
    input.nodes = handles.nodes;
    input.elem = handles.elem;
    input.veset = handles.veset;
    input.mags_vel = handles.mags_vel;
    input.peak_flow = handles.peak_flow;
    input.peak_flow_ori = handles.peak_flow_ori;
    input.WSS = handles.WSS;
    input.mags_wss = handles.mags_wss;    
    input.heart_rate = handles.heart_rate;
    input.vorticity = handles.VOR;
    input.id_mesh_inlet = 0;
    input.id_mesh_outlet = 0;
    input.id_inlet = 0;
    input.id_outlet = 0;
    
    input.inlet_exec = 0;
    input.outlet_exec = 0;
    
    input.id_view = 1;
    input.id_ipcmra = 1;
    input.id_mag = 0;
    
    input.list_n = handles.list_n;
    
    % identificador de vwerp
    input.id_vwerp = 0;
    
    id_while = 0;
    while(id_while == 0)
        
        GUIDE_LAPLACE(input)
        id_while = getappdata(0,'id_while');
        
        input.id_mesh_inlet = getappdata(0,'id_mesh_inlet');
        input.id_mesh_outlet = getappdata(0,'id_mesh_outlet');
        
        input.id_inlet = getappdata(0,'id_inlet');
        input.id_outlet = getappdata(0,'id_outlet');
        input.inlet_exec = getappdata(0,'inlet_exec');
        input.outlet_exec = getappdata(0,'outlet_exec');
        
        input.id_view = getappdata(0,'id_view');
        input.id_ipcmra = getappdata(0,'id_ipcmra');
        input.id_mag = getappdata(0,'id_mag');
        
        input.slider_id_axes1 = getappdata(0,'slider_id_axes1');
        
        handles.Laplace                         = getappdata(0,'Laplace');
        handles.centerline                      = getappdata(0,'centerline');
        handles.centerline_lapid                = getappdata(0,'centerline_lapid');
        handles.radius                          = getappdata(0,'radius');
        handles.diameter                        = getappdata(0,'diameter');
        handles.axial_unit_vectors              = getappdata(0,'axial_unit_vectors');
        handles.circumferential_unit_vectors    = getappdata(0,'circumferential_unit_vectors');
        handles.WSS_A                           = getappdata(0,'WSS_A');
        handles.WSS_C                           = getappdata(0,'WSS_C');
        handles.mag_WSS_A                       = getappdata(0,'mag_WSS_A');
        handles.mag_WSS_C                       = getappdata(0,'mag_WSS_C');
        handles.angle_axial_direction           = getappdata(0,'angle_axial_direction');
        handles.forward_velocity                = getappdata(0,'forward_velocity');
        handles.backward_velocity               = getappdata(0,'backward_velocity');
        handles.mag_forward_velocity            = getappdata(0,'mag_forward_velocity');
        handles.mag_backward_velocity           = getappdata(0,'mag_backward_velocity');
        handles.regurgitant_flow                = getappdata(0,'regurgitant_flow');
        handles.centerline_flow                 = getappdata(0,'centerline_flow');
        handles.eccentricity                    = getappdata(0,'eccentricity');
        
        handles.curvature                       = getappdata(0,'curvature'); % Julio Sotelo 28-05-2019
        handles.ellipticity                     = getappdata(0,'ellipticity'); % Julio Sotelo 28-05-2019
        handles.length_vessel                   = getappdata(0,'length_vessel'); % Julio Sotelo 28-05-2019
        handles.circulation                     = getappdata(0,'circulation'); % Julio Sotelo 28-05-2019
        handles.forward_vortex                  = getappdata(0,'forward_vortex'); % Julio Sotelo 28-05-2019
        handles.flattening                      = getappdata(0,'flattening'); % Julio Sotelo 28-05-2019
        handles.area                            = getappdata(0,'area'); % Julio Sotelo 28-05-2019
        handles.axial_circulation               = getappdata(0,'axial_circulation'); % Julio Sotelo 28-05-2019
        
        if id_while ==1
            
            if isempty(handles.Laplace)==0
                    
                h = msgbox('Please wait, transferring data ...','Wait','none');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                handles.Lrgb_aan  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
                handles.Lrgb_fve  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
                handles.Lrgb_bve  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3);
                handles.Lrgb_fov  = ones(size(handles.MR_PCA_FH,1)+4,size(handles.MR_PCA_FH,2)+4,size(handles.MR_PCA_FH,3)+4,size(handles.MR_PCA_FH,4),3); % Julio Sotelo 04-06-2019
                
                MASK = permute(double(abs(sum(handles.Lrgb,4)-3)>0),[2,1,3]);
                MASK(1:2,:,:) = 0; MASK(:,1:2,:) = 0; MASK(:,:,1:2) = 0; MASK(end-1:end,:,:) = 0; MASK(:,end-1:end,:) = 0; MASK(:,:,end-1:end) = 0;
                xd_seg = MASK.*handles.xd;
                yd_seg = MASK.*handles.yd;
                zd_seg = MASK.*handles.zd;
                xd_seg(MASK==0) = [];
                yd_seg(MASK==0) = [];
                zd_seg(MASK==0) = [];
                pos_voxel = [xd_seg(:),yd_seg(:),zd_seg(:)];
                
                [X,Y,Z] = meshgrid(1:size(handles.IPCMRA,1),1:size(handles.IPCMRA,2),1:size(handles.IPCMRA,3));
                X_seg = MASK.*X;
                Y_seg = MASK.*Y;
                Z_seg = MASK.*Z;
                X_seg(MASK==0) = [];
                Y_seg(MASK==0) = [];
                Z_seg(MASK==0) = [];
                pos_v = [X_seg(:),Y_seg(:),Z_seg(:)];

                % color asignation
                [~, ~, ind] = histcounts(handles.angle_axial_direction(:), size(cool(128), 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                cmap_vol = uint8(ind2rgb(ind, cool(128)) * 255);
                rgb_aan = cmap_vol/255;

                [~, ~, ind] = histcounts(handles.mag_forward_velocity(:), size(jet(128), 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                cmap_vol = uint8(ind2rgb(ind, jet(128)) * 255);
                rgb_fve = cmap_vol/255;

                [~, ~, ind] = histcounts(handles.mag_backward_velocity(:), size(jet(128), 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                cmap_vol = uint8(ind2rgb(ind, jet(128)) * 255);
                rgb_bve = cmap_vol/255;

                [~, ~, ind] = histcounts(handles.forward_vortex(:), size(cool(128), 1));
                ind = reshape(ind,[size(handles.veset,1) size(handles.veset,3)]);
                cmap_vol = uint8(ind2rgb(ind, cool(128)) * 255);
                rgb_fov = cmap_vol/255;

                for n=1:length(pos_voxel(:,1))
                    d = sqrt(sum((handles.nodes-repmat(pos_voxel(n,:),[length(handles.nodes(:,1)), 1])).^2,2));

                    handles.Lrgb_aan(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_aan(d<mean(handles.voxel_MR)*2,:,:),1);
                    handles.Lrgb_fve(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_fve(d<mean(handles.voxel_MR)*2,:,:),1);
                    handles.Lrgb_bve(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_bve(d<mean(handles.voxel_MR)*2,:,:),1);
                    handles.Lrgb_fov(pos_v(n,1)+2,pos_v(n,2)+2,pos_v(n,3)+2,:,:) = mean(rgb_fov(d<mean(handles.voxel_MR)*2,:,:),1);

                end
                handles.Lrgb_aan = handles.Lrgb_aan(3:end-2,3:end-2,3:end-2,:,:);
                handles.Lrgb_fve = handles.Lrgb_fve(3:end-2,3:end-2,3:end-2,:,:);
                handles.Lrgb_bve = handles.Lrgb_bve(3:end-2,3:end-2,3:end-2,:,:);
                handles.Lrgb_fov = handles.Lrgb_fov(3:end-2,3:end-2,3:end-2,:,:);

                list_string = {'...','Surface','Voxel','Mesh','Velocity','WSS','OSI',...
                               'Vorticity','Helicity Density','R. Helicity Density',...
                               'Viscous Dissipation','Energy Loss','Kinetic Energy',...
                               'Laplace','Centerline','Diameter','Radius','Axial Unit Vectors',...
                               'Circum. Unit Vectors','WSSa','WSSc','Axial Angle','Forward Velocity','Backward Velocity',...
                               'Regurgitant Flow','Centerline Flow','Eccentricity',...
                               'Curvature','Ellipticity','Length of Vessel','Forward Vortex','Flattening','Area','Axial Circulation'};
                set(handles.popupmenu1,'visible','on','String',list_string);


                handles.save_id_lap_mat     = 1;
                handles.save_id_cen_mat     = 1;
                handles.save_id_rad_mat     = 1;
                handles.save_id_dia_mat     = 1;
                handles.save_id_auv_mat     = 1;
                handles.save_id_cuv_mat     = 1;
                handles.save_id_wssa_mat    = 1;
                handles.save_id_wssc_mat    = 1;
                handles.save_id_aan_mat     = 1;
                handles.save_id_fve_mat     = 1;
                handles.save_id_bve_mat     = 1;
                handles.save_id_ref_mat     = 1;
                handles.save_id_cebf_mat    = 1;
                handles.save_id_ecc_mat     = 1;

                handles.save_id_cur_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_ell_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_len_mat     = 1; % Julio Sotelo 28-05-2019
%                 handles.save_id_cir_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_fov_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_fla_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_are_mat     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_aci_mat     = 1; % Julio Sotelo 28-05-2019
                
                handles.save_id_lap_vtu     = 1;
                handles.save_id_cen_vtu     = 1;
                handles.save_id_rad_vtu     = 1;
                handles.save_id_dia_vtu     = 1;
                handles.save_id_auv_vtu     = 1;
                handles.save_id_cuv_vtu     = 1;
                handles.save_id_wssa_vtu    = 1;
                handles.save_id_wssc_vtu    = 1;
                handles.save_id_aan_vtu     = 1;
                handles.save_id_fve_vtu     = 1;
                handles.save_id_bve_vtu     = 1;
                handles.save_id_ref_vtu     = 1;
                handles.save_id_cebf_vtu    = 1;
                handles.save_id_ecc_vtu     = 1;
                
                handles.save_id_cur_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_ell_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_len_vtu     = 1; % Julio Sotelo 28-05-2019
%                 handles.save_id_cir_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_fov_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_fla_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_are_vtu     = 1; % Julio Sotelo 28-05-2019
                handles.save_id_aci_vtu     = 1; % Julio Sotelo 28-05-2019
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                set(handles.popupmenu1,'Value',14)
                
                set(handles.pushbutton5,'visible','on')
                axes(handles.axes4);
                plot(0.0)
                axis off
                axes(handles.axes4);
                patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace,'CDataMapping','Scaled')
                hold on
                colormap(handles.axes4,'cool');
                c = colorbar(handles.axes4);
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
                caxis(handles.axes4, [handles.min_lap handles.max_lap]);
                hold off
                axis vis3d
                daspect([1,1,1])
                axis off
                view([-34,-51])
                if handles.slider_axes1 ==1, handles.slider_axes1_voxel = 0; else, handles.slider_axes1_voxel = (handles.slider_axes1*handles.voxel_MR(3))-handles.voxel_MR(3); end
                if handles.slider_axes2 ==1, handles.slider_axes2_voxel = 0; else, handles.slider_axes2_voxel = (handles.slider_axes2*handles.voxel_MR(2))-handles.voxel_MR(2); end
                if handles.slider_axes3 ==1, handles.slider_axes3_voxel = 0; else, handles.slider_axes3_voxel = (handles.slider_axes3*handles.voxel_MR(1))-handles.voxel_MR(1); end
                axes(handles.axes1);
                imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
                hold on
                plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
                plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
                Lrgb2d = squeeze(handles.Lrgb(:,:,handles.slider_axes1,:));
                himage = imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
                hold off
                axis image
                colormap(handles.axes1,'gray')
                axis off
                daspect([1 1 1])
                axes(handles.axes2);
                imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
                hold on
                plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
                plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
                Lrgb2d = squeeze(handles.Lrgb(:,handles.slider_axes2,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
                hold off
                axis image
                colormap(handles.axes2,'gray')
                axis off
                daspect([1 1 1])
                axes(handles.axes3);
                imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
                hold on
                plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
                plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
                Lrgb2d = squeeze(handles.Lrgb(handles.slider_axes3,:,:,:));
                himage = imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);
                hold off
                axis image
                colormap(handles.axes3,'gray')
                axis off
                daspect([1 1 1])
                handles.id_seg  = 1;
                handles.id_vel  = 0;
                handles.id_wss  = 0;
                handles.id_osi  = 0;
                handles.id_vor  = 0;
                handles.id_hd   = 0;
                handles.id_rhd  = 0;
                handles.id_vd   = 0;
                handles.id_el   = 0;
                handles.id_ke   = 0;
                handles.id_lap  = 1;
                handles.id_cen  = 0;
                handles.id_rad  = 0;
                handles.id_dia  = 0;
                handles.id_auv  = 0;
                handles.id_cuv  = 0;
                handles.id_wssa = 0;
                handles.id_wssc = 0;
                handles.id_aan  = 0;
                handles.id_fve  = 0;
                handles.id_bve  = 0;
                handles.id_ref  = 0;
                handles.id_cenf = 0;
                handles.id_ecc  = 0;
                
                handles.id_cur  = 0; % Julio Sotelo 28-05-2019
                handles.id_ell  = 0; % Julio Sotelo 28-05-2019
                handles.id_len  = 0; % Julio Sotelo 28-05-2019
%                 handles.id_cir  = 0; % Julio Sotelo 28-05-2019
                handles.id_fov  = 0; % Julio Sotelo 28-05-2019
                handles.id_fla  = 0; % Julio Sotelo 28-05-2019
                handles.id_are  = 0; % Julio Sotelo 28-05-2019
                handles.id_aci  = 0; % Julio Sotelo 28-05-2019
                
                set(handles.slider4,'visible','off')
                set(handles.pushbutton14,'visible','off')
                set(handles.text4,'Visible','off','String',['Cardiac Phase #: ',num2str(handles.peak_flow),' of ',num2str(size(handles.veset,3))])

                close(h)
                
                break

            else
                msgbox('The Laplace equation have not been calculated ...','Warning','warn')
            end
        end
    end

   
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton63_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[1,3,2,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[1,3,2,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[1,3,2,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[1,3,2,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[1,3,2,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[1,3,2,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[1,3,2]);
    
    handles.voxel_MR    = [handles.voxel_MR(1),handles.voxel_MR(3),handles.voxel_MR(2)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton64_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[3,2,1,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[3,2,1,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[3,2,1,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[3,2,1,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[3,2,1,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[3,2,1,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[3,2,1]);
    
    handles.voxel_MR    = [handles.voxel_MR(3),handles.voxel_MR(2),handles.voxel_MR(1)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;

    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton65_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton65 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[1,3,2,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[1,3,2,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[1,3,2,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[1,3,2,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[1,3,2,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[1,3,2,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[1,3,2]);
    
    handles.voxel_MR    = [handles.voxel_MR(1),handles.voxel_MR(3),handles.voxel_MR(2)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;

    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton66_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[2,1,3,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[2,1,3,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[2,1,3,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[2,1,3,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[2,1,3,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[2,1,3,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[2,1,3]);
    
    handles.voxel_MR    = [handles.voxel_MR(2),handles.voxel_MR(1),handles.voxel_MR(3)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;

    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton67_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton67 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[3,2,1,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[3,2,1,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[3,2,1,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[3,2,1,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[3,2,1,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[3,2,1,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[3,2,1]);
    
    handles.voxel_MR    = [handles.voxel_MR(3),handles.voxel_MR(2),handles.voxel_MR(1)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton68_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton68 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = permute(handles.MR_FFE_FH,[2,1,3,4]);
    handles.MR_FFE_AP   = permute(handles.MR_FFE_AP,[2,1,3,4]);
    handles.MR_FFE_RL   = permute(handles.MR_FFE_RL,[2,1,3,4]);
    handles.MR_PCA_FH   = permute(handles.MR_PCA_FH,[2,1,3,4]);
    handles.MR_PCA_AP   = permute(handles.MR_PCA_AP,[2,1,3,4]);
    handles.MR_PCA_RL   = permute(handles.MR_PCA_RL,[2,1,3,4]);
    handles.IPCMRA      = permute(handles.IPCMRA,[2,1,3]);
    
    handles.voxel_MR    = [handles.voxel_MR(2),handles.voxel_MR(1),handles.voxel_MR(3)];
    
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;

    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



    set(handles.pushbutton63, 'Units', 'pixels');
    set(handles.pushbutton64, 'Units', 'pixels');
    set(handles.pushbutton65, 'Units', 'pixels');
    set(handles.pushbutton66, 'Units', 'pixels');
    set(handles.pushbutton67, 'Units', 'pixels');
    set(handles.pushbutton68, 'Units', 'pixels');
    
    handles.pushbutton63_size = get(handles.pushbutton63, 'Position');
    handles.pushbutton64_size = get(handles.pushbutton64, 'Position');
    handles.pushbutton65_size = get(handles.pushbutton65, 'Position');
    handles.pushbutton66_size = get(handles.pushbutton66, 'Position');
    handles.pushbutton67_size = get(handles.pushbutton67, 'Position');
    handles.pushbutton68_size = get(handles.pushbutton68, 'Position');
    
    set(handles.pushbutton63, 'Units', 'normalized');
    set(handles.pushbutton64, 'Units', 'normalized');
    set(handles.pushbutton65, 'Units', 'normalized');
    set(handles.pushbutton66, 'Units', 'normalized');
    set(handles.pushbutton67, 'Units', 'normalized');
    set(handles.pushbutton68, 'Units', 'normalized');
    
    idx = mod(min(handles.pushbutton63_size(3:4)),2)>1;
    w = floor(min(handles.pushbutton64_size(3:4)));
    w(idx) = w(idx)+1;
    
    im_diagonal = imread('Symbols/DIAGONAL.png');
    im_left_right = imread('Symbols/LEFT-RIGHT.png');
    im_top_bottom = imread('Symbols/TOP-BOTTOM.png');
    
    set(handles.pushbutton63,'CData',double(imresize(im_left_right,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton64,'CData',double(imresize(im_top_bottom,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton65,'CData',double(imresize(im_left_right,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton66,'CData',double(imresize(im_diagonal,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton67,'CData',double(imresize(im_top_bottom,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton68,'CData',double(imresize(im_diagonal,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton69,'visible','on')
    set(handles.pushbutton70,'visible','on')
    set(handles.pushbutton71,'visible','on')

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton69_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton69 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.pushbutton63,'visible','off')
    set(handles.pushbutton64,'visible','off')
    set(handles.pushbutton65,'visible','off')
    set(handles.pushbutton66,'visible','off')
    set(handles.pushbutton67,'visible','off')
    set(handles.pushbutton68,'visible','off')
    set(handles.pushbutton69,'visible','off')
    set(handles.pushbutton70,'visible','off')
    set(handles.pushbutton71,'visible','off')

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton70_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.pushbutton63,'visible','off')
    set(handles.pushbutton64,'visible','off')
    set(handles.pushbutton65,'visible','off')
    set(handles.pushbutton66,'visible','off')
    set(handles.pushbutton67,'visible','off')
    set(handles.pushbutton68,'visible','off')
    set(handles.pushbutton69,'visible','off')
    set(handles.pushbutton70,'visible','off')
    set(handles.pushbutton71,'visible','off')

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton71_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton71 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.pushbutton63,'visible','off')
    set(handles.pushbutton64,'visible','off')
    set(handles.pushbutton65,'visible','off')
    set(handles.pushbutton66,'visible','off')
    set(handles.pushbutton67,'visible','off')
    set(handles.pushbutton68,'visible','off')
    set(handles.pushbutton69,'visible','off')
    set(handles.pushbutton70,'visible','off')
    set(handles.pushbutton71,'visible','off')

handles.output = hObject;
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton6_CreateFcn(hObject, eventdata, handles)
function File_Callback(hObject, eventdata, handles)
function Load_Callback(hObject, eventdata, handles)
function Save_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton72.
function pushbutton72_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton72 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,2);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,2);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,2);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,2);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,2);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,2);
    handles.IPCMRA      = flip(handles.IPCMRA,2);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton73.
function pushbutton73_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,1);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,1);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,1);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,1);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,1);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,1);
    handles.IPCMRA      = flip(handles.IPCMRA,1);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton74.
function pushbutton74_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,3);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,3);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,3);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,3);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,3);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,3);
    handles.IPCMRA      = flip(handles.IPCMRA,3);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton75.
function pushbutton75_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,1);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,1);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,1);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,1);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,1);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,1);
    handles.IPCMRA      = flip(handles.IPCMRA,1);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton76.
function pushbutton76_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,3);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,3);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,3);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,3);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,3);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,3);
    handles.IPCMRA      = flip(handles.IPCMRA,3);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton77.
function pushbutton77_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.MR_FFE_FH   = flip(handles.MR_FFE_FH,2);
    handles.MR_FFE_AP   = flip(handles.MR_FFE_AP,2);
    handles.MR_FFE_RL   = flip(handles.MR_FFE_RL,2);
    handles.MR_PCA_FH   = flip(handles.MR_PCA_FH,2);
    handles.MR_PCA_AP   = flip(handles.MR_PCA_AP,2);
    handles.MR_PCA_RL   = flip(handles.MR_PCA_RL,2);
    handles.IPCMRA      = flip(handles.IPCMRA,2);
        
    % Julio Sotelo 23-11-2018
    handles.ANG = handles.IPCMRA;
    
    var_mag = [min(handles.MR_FFE_AP(:)) min(handles.MR_FFE_FH(:)) min(handles.MR_FFE_RL(:))]; % Julio Sotelo 23-11-2018
    [~,c,~] = find(var_mag>=0); % Julio Sotelo 23-11-2018
    if sum(var_mag) ~= 3 % Julio Sotelo 23-11-2018
        if c(1)==1
            handles.MAG = mean(handles.MR_FFE_AP,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_FH = handles.MR_FFE_AP;
            handles.MR_FFE_RL = handles.MR_FFE_AP;
            
        elseif c(1)==2
            handles.MAG = mean(handles.MR_FFE_FH,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_FH;
            handles.MR_FFE_RL = handles.MR_FFE_FH;
            
        elseif c(1)==3
            handles.MAG = mean(handles.MR_FFE_RL,4); % Julio Sotelo 23-11-2018
            handles.MAG = (abs(handles.MAG)/max(handles.MAG(:)))*255; % Julio Sotelo 23-11-2018
            handles.MR_FFE_AP = handles.MR_FFE_RL;
            handles.MR_FFE_FH = handles.MR_FFE_RL;
            
        end % Julio Sotelo 23-11-2018
        
    end % Julio Sotelo 23-11-2018
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [a,b,c,d] = size(handles.MR_FFE_FH);
    handles.a = a;
    handles.b = b;
    handles.c = c;
    handles.d = d;
    
    % Julio Sotelo 23-11-2018
    handles.position_sag = [0 0 (b-1)*handles.voxel_MR(2) (a-1)*handles.voxel_MR(1)];
    handles.position_axi = [0 0 (c-1)*handles.voxel_MR(3) (a-1)*handles.voxel_MR(1)];
    handles.position_cor = [0 0 (c-1)*handles.voxel_MR(3) (b-1)*handles.voxel_MR(2)];
        
    handles.position_sag_cor = [];
    handles.position_axi_cor = [];
    handles.position_cor_cor = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    handles.input2 = handles.IPCMRA;
    handles.SEG = zeros(size(handles.IPCMRA));
    handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
    [L,NUM] = bwlabeln(handles.SEG,6);
    handles.L = L;
    handles.NUM = NUM;
    handles.id_seg = 0;
    handles.min_value_th = min(handles.IPCMRA(:));
    handles.max_value_th = max(handles.IPCMRA(:));
    handles.id_resizing = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y,Z] = meshgrid(0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,3)-1);
    handles.xd = X*handles.voxel_MR(1);
    handles.yd = Y*handles.voxel_MR(2);
    handles.zd = Z*handles.voxel_MR(3);
    handles.slider_axes1 = round(c/2);
    handles.slider_axes2 = round(b/2);
    handles.slider_axes3 = round(a/2);
    if handles.slider_axes1 ==1, handles.slider_axes1_voxel = handles.voxel_MR(3); else handles.slider_axes1_voxel = handles.slider_axes1*handles.voxel_MR(3); end
    if handles.slider_axes2 ==1, handles.slider_axes2_voxel = handles.voxel_MR(2); else handles.slider_axes2_voxel = handles.slider_axes2*handles.voxel_MR(2); end
    if handles.slider_axes3 ==1, handles.slider_axes3_voxel = handles.voxel_MR(1); else handles.slider_axes3_voxel = handles.slider_axes3*handles.voxel_MR(1); end
    axes(handles.axes1);
    imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
    hold on
    plot([handles.yd(1,1),handles.yd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes2_voxel,handles.slider_axes2_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--r','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes2);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,handles.slider_axes2,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes3_voxel,handles.slider_axes3_voxel]','--b','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.xd(1,1),handles.xd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    axes(handles.axes3);
    imagesc([handles.zd(1,1),handles.zd(end,end)]',[handles.yd(1,1),handles.yd(end,end)]',squeeze(handles.IPCMRA(handles.slider_axes3,:,:)))
    hold on
    plot([handles.zd(1,1),handles.zd(end,end)]',[handles.slider_axes2_voxel,handles.slider_axes2_voxel]','--r','Linewidth',1)
    plot([handles.slider_axes1_voxel,handles.slider_axes1_voxel]',[handles.yd(1,1),handles.yd(end,end)]','--g','Linewidth',1)
    hold off
    axis image
    colormap gray
    axis off
    daspect([1 1 1])
    slider_step(1) = 1/c;
    slider_step(2) = 0.1;
    set(handles.slider1,'Value', round(c/2)/c,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider1,'visible','on')
    slider_step(1) = 1/b;
    slider_step(2) = 0.1;
    set(handles.slider2,'Value', round(b/2)/b,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider2,'visible','on')
    slider_step(1) = 1/a;
    slider_step(2) = 0.1;
    set(handles.slider3,'Value', round(a/2)/a,'sliderstep',slider_step,'max',1,'min',0)
    set(handles.slider3,'visible','on')
    set(handles.text1,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,2)),' Im: ',num2str(round(c/2)),' / ',num2str(c),' Type: ',handles.type])
    set(handles.text2,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,1)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(b/2)),' / ',num2str(b),' Type: ',handles.type])
    set(handles.text3,'visible','on','String',['Image size: ',num2str(size(handles.IPCMRA,2)),' x ',num2str(size(handles.IPCMRA,3)),' Im: ',num2str(round(a/2)),' / ',num2str(a),' Type: ',handles.type])
    
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton78.
function pushbutton78_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    set(handles.pushbutton72,'visible','off')
    set(handles.pushbutton73,'visible','off')
    set(handles.pushbutton74,'visible','off')
    set(handles.pushbutton75,'visible','off')
    set(handles.pushbutton76,'visible','off')
    set(handles.pushbutton77,'visible','off')
    set(handles.pushbutton78,'visible','off')
    set(handles.pushbutton79,'visible','off')
    set(handles.pushbutton80,'visible','off')

handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton79.
function pushbutton79_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    set(handles.pushbutton72,'visible','off')
    set(handles.pushbutton73,'visible','off')
    set(handles.pushbutton74,'visible','off')
    set(handles.pushbutton75,'visible','off')
    set(handles.pushbutton76,'visible','off')
    set(handles.pushbutton77,'visible','off')
    set(handles.pushbutton78,'visible','off')
    set(handles.pushbutton79,'visible','off')
    set(handles.pushbutton80,'visible','off')

handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton80.
function pushbutton80_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton80 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    set(handles.pushbutton72,'visible','off')
    set(handles.pushbutton73,'visible','off')
    set(handles.pushbutton74,'visible','off')
    set(handles.pushbutton75,'visible','off')
    set(handles.pushbutton76,'visible','off')
    set(handles.pushbutton77,'visible','off')
    set(handles.pushbutton78,'visible','off')
    set(handles.pushbutton79,'visible','off')
    set(handles.pushbutton80,'visible','off')

handles.output = hObject;
guidata(hObject, handles);

% --------------------------------------------------------------------
function uipushtool9_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.pushbutton72, 'Units', 'pixels');
    set(handles.pushbutton73, 'Units', 'pixels');
    set(handles.pushbutton74, 'Units', 'pixels');
    set(handles.pushbutton75, 'Units', 'pixels');
    set(handles.pushbutton76, 'Units', 'pixels');
    set(handles.pushbutton77, 'Units', 'pixels');
    
    handles.pushbutton72_size = get(handles.pushbutton72, 'Position');
    handles.pushbutton73_size = get(handles.pushbutton73, 'Position');
    handles.pushbutton74_size = get(handles.pushbutton74, 'Position');
    handles.pushbutton75_size = get(handles.pushbutton75, 'Position');
    handles.pushbutton76_size = get(handles.pushbutton76, 'Position');
    handles.pushbutton77_size = get(handles.pushbutton77, 'Position');
    
    set(handles.pushbutton72, 'Units', 'normalized');
    set(handles.pushbutton73, 'Units', 'normalized');
    set(handles.pushbutton74, 'Units', 'normalized');
    set(handles.pushbutton75, 'Units', 'normalized');
    set(handles.pushbutton76, 'Units', 'normalized');
    set(handles.pushbutton77, 'Units', 'normalized');
    
    idx = mod(min(handles.pushbutton72_size(3:4)),2)>1;
    w = floor(min(handles.pushbutton73_size(3:4)));
    w(idx) = w(idx)+1;
    
    im_diagonal = imread('Symbols/DIAGONAL.png');
    im_left_right = imread('Symbols/LEFT-RIGHT.png');
    im_top_bottom = imread('Symbols/TOP-BOTTOM.png');
    
    set(handles.pushbutton72,'CData',double(imresize(im_left_right,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton73,'CData',double(imresize(im_top_bottom,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton74,'CData',double(imresize(im_left_right,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton75,'CData',double(imresize(im_top_bottom,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton76,'CData',double(imresize(im_left_right,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton77,'CData',double(imresize(im_top_bottom,[w-4 w-4],'method','nearest')>0),'visible','on')
    set(handles.pushbutton78,'visible','on')
    set(handles.pushbutton79,'visible','on')
    set(handles.pushbutton80,'visible','on')

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton81.
function pushbutton81_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton81 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Flow Quantification
    input                     = [];
    input.SEG                 = handles.SEG;
    input.IPCMRA              = handles.IPCMRA;
    input.voxel_MR            = handles.voxel_MR;
    input.L                   = handles.L;
    input.Lrgb                = handles.Lrgb;
    input.Lrgb_vel            = handles.Lrgb_vel;
    input.NUM                 = handles.NUM;
    input.xd                  = handles.xd;
    input.yd                  = handles.yd;
    input.zd                  = handles.zd;
    input.a                   = handles.a;
    input.b                   = handles.b;
    input.c                   = handles.c;
    input.d                   = handles.d;
    input.slider_id_axes1     = handles.slider_axes1;
    input.MR_FFE_FH           = handles.MR_FFE_FH;
    input.faces               = handles.faces;
    input.nodes               = handles.nodes;
    input.elem                = handles.elem;
    input.veset               = handles.veset;
    input.id_mesh_inlet_flow  = 0;
    input.id_image            = 1;    
    input.id_view             = 1;
    input.id_result           = 0;
    input.heart_rate          = handles.heart_rate;
    input.id_selected_faceid  = [];
    input.faceid              = [];
    input.cutpos              = [];
    input.time                = [];
    input.peak_flow           = handles.peak_flow;
    input.flow                = [];
    input.net_flow            = [];
    input.max_velocity        = [];
    input.min_velocity        = [];
    input.velocity_proj       = [];
    input.id_section          = [];
    input.veset_out           = handles.veset;
    input.Lrgb                = handles.Lrgb;
    input.Lrgb_vel            = handles.Lrgb_vel;
    input.mags_vel            = handles.mags_vel;
    input.MR_PCA_FH           = handles.MR_PCA_FH;
    input.MR_PCA_AP           = handles.MR_PCA_AP;
    input.MR_PCA_RL           = handles.MR_PCA_RL;
    input.MR_PCA_FH_smooth    = handles.MR_PCA_FH_smooth;
    input.MR_PCA_AP_smooth    = handles.MR_PCA_AP_smooth;
    input.MR_PCA_RL_smooth    = handles.MR_PCA_RL_smooth;
        
    id_while = 0;
    while(id_while == 0)

        GUIDE_FLOW(input);

        id_while = getappdata(0,'id_while');
        input.id_mesh_inlet_flow = getappdata(0,'id_mesh_inlet_flow');
        input.id_view = getappdata(0,'id_view');
        input.id_image = getappdata(0,'id_image');
        input.id_result = getappdata(0,'id_result');
        input.slider_id_axes1 = getappdata(0,'slider_id_axes1');
        input.id_selected_faceid = getappdata(0,'id_selected_faceid');
        input.faceid = getappdata(0,'faceid');
        input.cutpos = getappdata(0,'cutpos');
        input.time = getappdata(0,'time');
        input.peak_flow = getappdata(0,'peak_flow');
        input.flow = getappdata(0,'flow');
        input.net_flow = getappdata(0,'net_flow');
        input.max_velocity = getappdata(0,'max_velocity');
        input.min_velocity = getappdata(0,'min_velocity');
        input.velocity_proj = getappdata(0,'velocity_proj');
        input.id_section = getappdata(0,'id_section');
        input.veset_out = getappdata(0,'veset_out');
        input.Lrgb_vel = getappdata(0,'Lrgb_vel');
        input.mags_vel = getappdata(0,'mags_vel');
        input.MR_PCA_FH = getappdata(0,'MR_PCA_FH');
        input.MR_PCA_AP = getappdata(0,'MR_PCA_AP');
        input.MR_PCA_RL = getappdata(0,'MR_PCA_RL');
        input.MR_PCA_FH_smooth = getappdata(0,'MR_PCA_FH_smooth');
        input.MR_PCA_AP_smooth = getappdata(0,'MR_PCA_AP_smooth');
        input.MR_PCA_RL_smooth = getappdata(0,'MR_PCA_RL_smooth');


    end
    
    
    
    handles.MR_PCA_FH           = input.MR_PCA_FH;
    handles.MR_PCA_AP           = input.MR_PCA_AP;
    handles.MR_PCA_RL           = input.MR_PCA_RL;
    handles.MR_PCA_FH_smooth    = input.MR_PCA_FH_smooth;
    handles.MR_PCA_AP_smooth    = input.MR_PCA_AP_smooth;
    handles.MR_PCA_RL_smooth    = input.MR_PCA_RL_smooth;
    handles.mags_vel            = input.mags_vel;
    handles.Lrgb_vel            = input.Lrgb_vel;
    handles.veset               = input.veset_out;
    handles.time                = input.time;
    handles.peak_flow           = input.peak_flow;
    handles.peak_flow_ori       = input.peak_flow;
    handles.flow                = input.flow;
    handles.net_flow            = input.net_flow;
    handles.max_velocity        = input.max_velocity;
    handles.min_velocity        = input.min_velocity;
    
    if isempty(handles.flow)==0

        handles.save_id_tim_mat     = 1; % Julio Sotelo 28-05-2019 time
        handles.save_id_flo_mat     = 1; % Julio Sotelo 28-05-2019 flow
        handles.save_id_nfl_mat     = 1; % Julio Sotelo 28-05-2019 net_flow
        handles.save_id_mav_mat     = 1; % Julio Sotelo 28-05-2019 max_velocity
        handles.save_id_miv_mat     = 1; % Julio Sotelo 28-05-2019 min_velocity

        handles.save_id_tim_csv     = 1; % Julio Sotelo 28-05-2019 time
        handles.save_id_flo_csv     = 1; % Julio Sotelo 28-05-2019 flow
        handles.save_id_nfl_csv     = 1; % Julio Sotelo 28-05-2019 net_flow
        handles.save_id_mav_csv     = 1; % Julio Sotelo 28-05-2019 max_velocity
        handles.save_id_miv_csv     = 1; % Julio Sotelo 28-05-2019 min_velocity
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton13,'visible','on');
        
    end

handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton83.
function pushbutton83_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton83 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% vWERP Quantification
    input                       = [];
    input.SEG                   = handles.SEG;
    input.IPCMRA                = handles.IPCMRA;
    input.voxel_MR              = handles.voxel_MR;
    input.L                     = handles.L;
    input.Lrgb                  = handles.Lrgb;
    input.Lrgb_vel              = handles.Lrgb_vel;
    input.NUM                   = handles.NUM;
    input.xd                    = handles.xd;
    input.yd                    = handles.yd;
    input.zd                    = handles.zd;
    input.a                     = handles.a;
    input.b                     = handles.b;
    input.c                     = handles.c;
    input.d                     = handles.d;
    input.slider_id_axes1       = handles.slider_axes1;
    input.MR_FFE_FH             = handles.MR_FFE_FH;
    input.MR_PCA_AP             = handles.MR_PCA_AP;
    input.MR_PCA_FH             = handles.MR_PCA_FH;
    input.MR_PCA_RL             = handles.MR_PCA_RL;
    input.faces                 = handles.faces;
    input.nodes                 = handles.nodes;
    input.elem                  = handles.elem;
    input.veset                 = handles.veset;
    
    input.mags_vel              = handles.mags_vel;
    input.peak_flow             = handles.peak_flow;
    %input.peak_flow_ori         = handles.peak_flow_ori;
    %input.WSS                   = handles.WSS;
    %input.mags_wss              = handles.mags_wss;    
    input.heart_rate            = handles.heart_rate;
    %input.vorticity             = handles.VOR;
    input.id_mesh_inlet         = 0;
    input.id_mesh_outlet        = 0;
    input.id_inlet              = 0;
    input.id_outlet             = 0;
    
    input.inlet_exec            = 0;
    input.outlet_exec           = 0;
    
    input.id_view               = 1;
    input.id_ipcmra             = 1;
    input.id_mag                = 0;
    
    %input.list_n = handles.list_n;

    % identificador vwerp
    input.id_vwerp = 1;
    
    id_while = 0;
    while(id_while == 0)
        
        GUIDE_LAPLACE(input)
        id_while = getappdata(0,'id_while');
        
        input.id_mesh_inlet = getappdata(0,'id_mesh_inlet');
        input.id_mesh_outlet = getappdata(0,'id_mesh_outlet');
        
        input.id_inlet = getappdata(0,'id_inlet');
        input.id_outlet = getappdata(0,'id_outlet');
        input.inlet_exec = getappdata(0,'inlet_exec');
        input.outlet_exec = getappdata(0,'outlet_exec');
        
        input.id_view = getappdata(0,'id_view');
        input.id_ipcmra = getappdata(0,'id_ipcmra');
        input.id_mag = getappdata(0,'id_mag');
        
        input.slider_id_axes1 = getappdata(0,'slider_id_axes1');
        
%         handles.Laplace                         = getappdata(0,'Laplace');
%         handles.centerline                      = getappdata(0,'centerline');
%         handles.centerline_lapid                = getappdata(0,'centerline_lapid');
%         handles.radius                          = getappdata(0,'radius');
%         handles.diameter                        = getappdata(0,'diameter');
%         handles.axial_unit_vectors              = getappdata(0,'axial_unit_vectors');
%         handles.circumferential_unit_vectors    = getappdata(0,'circumferential_unit_vectors');
%         handles.WSS_A                           = getappdata(0,'WSS_A');
%         handles.WSS_C                           = getappdata(0,'WSS_C');
%         handles.mag_WSS_A                       = getappdata(0,'mag_WSS_A');
%         handles.mag_WSS_C                       = getappdata(0,'mag_WSS_C');
%         handles.angle_axial_direction           = getappdata(0,'angle_axial_direction');
%         handles.forward_velocity                = getappdata(0,'forward_velocity');
%         handles.backward_velocity               = getappdata(0,'backward_velocity');
%         handles.mag_forward_velocity            = getappdata(0,'mag_forward_velocity');
%         handles.mag_backward_velocity           = getappdata(0,'mag_backward_velocity');
%         handles.regurgitant_flow                = getappdata(0,'regurgitant_flow');
%         handles.centerline_flow                 = getappdata(0,'centerline_flow');
%         handles.eccentricity                    = getappdata(0,'eccentricity');
%         
%         handles.curvature                       = getappdata(0,'curvature'); % Julio Sotelo 28-05-2019
%         handles.ellipticity                     = getappdata(0,'ellipticity'); % Julio Sotelo 28-05-2019
%         handles.length_vessel                   = getappdata(0,'length_vessel'); % Julio Sotelo 28-05-2019
%         handles.circulation                     = getappdata(0,'circulation'); % Julio Sotelo 28-05-2019
%         handles.forward_vortex                  = getappdata(0,'forward_vortex'); % Julio Sotelo 28-05-2019
%         handles.flattening                      = getappdata(0,'flattening'); % Julio Sotelo 28-05-2019
%         handles.area                            = getappdata(0,'area'); % Julio Sotelo 28-05-2019
%         handles.axial_circulation               = getappdata(0,'axial_circulation'); % Julio Sotelo 28-05-2019
        
        if id_while ==1
            
           g = msgbox('The vWERP data was saved in the Output_vWERP folder ...');
           pause(2)
           close(g)
           
        end
    end

   
handles.output = hObject;
guidata(hObject, handles);
