function varargout = OFFSET_ERR_AND_NOISE_MASKING(varargin)
% OFFSET_ERR_AND_NOISE_MASKING MATLAB code for OFFSET_ERR_AND_NOISE_MASKING.fig
%      OFFSET_ERR_AND_NOISE_MASKING, by itself, creates a new OFFSET_ERR_AND_NOISE_MASKING or raises the existing
%      singleton*.
%
%      H = OFFSET_ERR_AND_NOISE_MASKING returns the handle to a new OFFSET_ERR_AND_NOISE_MASKING or the handle to
%      the existing singleton*.
%
%      OFFSET_ERR_AND_NOISE_MASKING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OFFSET_ERR_AND_NOISE_MASKING.M with the given input arguments.
%
%      OFFSET_ERR_AND_NOISE_MASKING('Property','Value',...) creates a new OFFSET_ERR_AND_NOISE_MASKING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OFFSET_ERR_AND_NOISE_MASKING_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OFFSET_ERR_AND_NOISE_MASKING_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OFFSET_ERR_AND_NOISE_MASKING

% Last Modified by GUIDE v2.5 03-Jun-2022 17:20:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OFFSET_ERR_AND_NOISE_MASKING_OpeningFcn, ...
                   'gui_OutputFcn',  @OFFSET_ERR_AND_NOISE_MASKING_OutputFcn, ...
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


% --- Executes just before OFFSET_ERR_AND_NOISE_MASKING is made visible.
function OFFSET_ERR_AND_NOISE_MASKING_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OFFSET_ERR_AND_NOISE_MASKING (see VARARGIN)

% Choose default command line output for OFFSET_ERR_AND_NOISE_MASKING
    handles.output = hObject;
    
    handles.id     = varargin{1}.id+1;
    handles.MR_FFE_FH    = varargin{1}.MR_FFE_FH;
    handles.MR_FFE_AP    = varargin{1}.MR_FFE_AP;
    handles.MR_FFE_RL    = varargin{1}.MR_FFE_RL;
    handles.MR_PCA_FH    = varargin{1}.MR_PCA_FH;
    handles.MR_PCA_AP    = varargin{1}.MR_PCA_AP;
    handles.MR_PCA_RL    = varargin{1}.MR_PCA_RL;
    handles.IPCMRA       = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);
    handles.Magnitude    = mean(handles.MR_FFE_FH,4);
    handles.voxel_MR     = varargin{1}.voxel_MR;
    handles.VENC         = varargin{1}.VENC;
    handles.heart_rate   = varargin{1}.heart_rate;
    handles.type         = varargin{1}.type;
    handles.id_while     = varargin{1}.id_while;
    
    if handles.id ==1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.list_string1 = {'Magnitud Image','IPCMRA'};
        set(handles.popupmenu1,'String',handles.list_string1,'Value',1);
        handles.list_string2 = {'Select the property to modify ...','Static Tissue','Noise','ROI'};
        set(handles.popupmenu2,'String',handles.list_string2,'Value',1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [handles.a, handles.b, handles.c, handles.d] = size(handles.MR_FFE_FH);
        [X,Y,Z] = meshgrid(0:handles.b-1,0:handles.a-1,0:handles.c-1);
        
        handles.xd = X.*handles.voxel_MR(2);
        handles.yd = Y.*handles.voxel_MR(1);
        handles.zd = Z.*handles.voxel_MR(3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        value1 = round(handles.c/2)/handles.c;
        value2 = round(handles.b/2)/handles.b;
        value3 = round(handles.a/2)/handles.a;

        pp = 1/(handles.c+1); slider_step(1) = pp; slider_step(2) = pp;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0,'Visible','on','value',value1)

        pp = 1/(handles.b+1); slider_step(1) = pp; slider_step(2) = pp;
        set(handles.slider2,'sliderstep',slider_step,'max',1,'min',0,'Visible','on','value',value2)

        pp = 1/(handles.a+1); slider_step(1) = pp; slider_step(2) = pp;
        set(handles.slider3,'sliderstep',slider_step,'max',1,'min',0,'Visible','on','value',value3)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.text2,'Visible','on','String',['Slice: ',num2str(round(handles.c/2))]); % text panel 1
        set(handles.text3,'Visible','on','String',['Slice: ',num2str(round(handles.b/2))]); % text panel 2
        set(handles.text4,'Visible','on','String',['Slice: ',num2str(round(handles.a/2))]); % text panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.slider_value1 = round(handles.c/2);
        handles.slider_value2 = round(handles.b/2);
        handles.slider_value3 = round(handles.a/2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of static tissue
        handles.std_mat = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.std_mat(:,:,:,1) = std(handles.MR_PCA_FH,1,4);
        handles.std_mat(:,:,:,2) = std(handles.MR_PCA_AP,1,4);
        handles.std_mat(:,:,:,3) = std(handles.MR_PCA_RL,1,4);

        handles.thres_static = 0.5;
        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);

        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of noise masking

        handles.thres_noise = 0.5;
        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);

        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of regions of interest (ROI)
        % initial points
        min_value = round(min([handles.a, handles.b, handles.c])*0.1);

        handles.points_view1 = [    min_value, min_value ;...
                                    min_value,handles.a-(min_value) ;... 
                                    handles.b-(min_value),min_value ;... 
                                    handles.b-(min_value),handles.a-(min_value)]; % modificado JSOTELO
        
        handles.points_view2 = [    min_value, min_value ; 
                                    min_value,handles.c-(min_value) ; 
                                    handles.a-(min_value),min_value ; 
                                    handles.a-(min_value),handles.c-(min_value)]; % modificado JSOTELO
                
        handles.points_view3 = [    min_value, min_value ; 
                                    min_value,handles.c-(min_value) ; 
                                    handles.b-(min_value),min_value ; 
                                    handles.b-(min_value),handles.c-(min_value)]; % modificado JSOTELO
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        handles.points_view1(:,1) = handles.points_view1(:,1).*handles.voxel_MR(2);
        handles.points_view1(:,2) = handles.points_view1(:,2).*handles.voxel_MR(1);

        handles.points_view2(:,1) = handles.points_view2(:,1).*handles.voxel_MR(2);
        handles.points_view2(:,2) = handles.points_view2(:,2).*handles.voxel_MR(3);

        handles.points_view3(:,1) = handles.points_view3(:,1).*handles.voxel_MR(1);
        handles.points_view3(:,2) = handles.points_view3(:,2).*handles.voxel_MR(3);

        handles.mask_roi = zeros(size(handles.MR_PCA_FH));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.id_image = get(handles.popupmenu1, 'Value');
        handles.id_variable = get(handles.popupmenu2,'value');

        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
                plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
                plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
                plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            end

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
                plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
                plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
                plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
                plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
                plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
                plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
                plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
                plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
                plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            end

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
                plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
                plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
                plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
                plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
                plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
                plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end
    else
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.id_image = get(handles.popupmenu1, 'Value');
        handles.id_variable = get(handles.popupmenu2,'value');

        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
                plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
                plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
                plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            end

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
                plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
                plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
                plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
                plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
                plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
                plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
                himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
                plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
                plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
                plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            end

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
                plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
                plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
                plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            if handles.id_variable==2
                Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==3
                Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
                himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
                cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
                cdata = double(cdata)*0.5;
                set(himage, 'AlphaData', cdata);

            elseif handles.id_variable==4
                plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
                plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
                plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
                plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            end

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end
        
    end

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes OFFSET_ERR_AND_NOISE_MASKING wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OFFSET_ERR_AND_NOISE_MASKING_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
 
if handles.id_while == 0
    varargout{1} = handles.output;
    MR_FFE_FH = handles.MR_FFE_FH;
    MR_FFE_AP = handles.MR_FFE_AP;
    MR_FFE_RL = handles.MR_FFE_RL;
    MR_PCA_FH = handles.MR_PCA_FH;
    MR_PCA_AP = handles.MR_PCA_AP;
    MR_PCA_RL = handles.MR_PCA_RL;
    voxel_MR = handles.voxel_MR;
    VENC = handles.VENC;
    heart_rate = handles.heart_rate;
    type = handles.type;
    setappdata(0,'MR_FFE_FH',MR_FFE_FH);
    setappdata(0,'MR_FFE_AP',MR_FFE_AP);
    setappdata(0,'MR_FFE_RL',MR_FFE_RL);
    setappdata(0,'MR_PCA_FH',MR_PCA_FH);
    setappdata(0,'MR_PCA_AP',MR_PCA_AP);
    setappdata(0,'MR_PCA_RL',MR_PCA_RL);
    setappdata(0,'voxel_MR',voxel_MR);
    setappdata(0,'VENC',VENC);
    setappdata(0,'heart_rate',heart_rate);
    setappdata(0,'type',type);
    id = handles.id;
    setappdata(0,'id',id);
    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
elseif handles.id_while == 1
    varargout{1} = handles.output;
    MR_FFE_FH = handles.MR_FFE_FH;
    MR_FFE_AP = handles.MR_FFE_AP;
    MR_FFE_RL = handles.MR_FFE_RL;
    MR_PCA_FH = handles.MR_PCA_FH;
    MR_PCA_AP = handles.MR_PCA_AP;
    MR_PCA_RL = handles.MR_PCA_RL;
    voxel_MR = handles.voxel_MR;
    VENC = handles.VENC;
    heart_rate = handles.heart_rate;
    type = handles.type;
    setappdata(0,'MR_FFE_FH',MR_FFE_FH);
    setappdata(0,'MR_FFE_AP',MR_FFE_AP);
    setappdata(0,'MR_FFE_RL',MR_FFE_RL);
    setappdata(0,'MR_PCA_FH',MR_PCA_FH);
    setappdata(0,'MR_PCA_AP',MR_PCA_AP);
    setappdata(0,'MR_PCA_RL',MR_PCA_RL);
    setappdata(0,'voxel_MR',voxel_MR);
    setappdata(0,'VENC',VENC);
    setappdata(0,'heart_rate',heart_rate);
    setappdata(0,'type',type);
    id = handles.id;
    setappdata(0,'id',id);
    id_while = handles.id_while;
    setappdata(0,'id_while',id_while);
    delete(handles.figure1);
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, ~, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.id_image = get(handles.popupmenu1, 'Value');

    
    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, ~, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, ~, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    handles.id_variable = get(handles.popupmenu2, 'Value');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if handles.id_variable==1

        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')


    elseif handles.id_variable==2

        set(handles.pushbutton1,'Visible','on'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','on'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','on'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','on'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','on'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','on'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','on'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','on'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','on'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','on')

        value = handles.thres_static;
        pp=1/100;
        slider_step(1) = pp;
        slider_step(2) = pp;
        set(handles.slider4,'sliderstep',slider_step,'max',1,'min',0,'Visible','on','value',value);    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of static tissue
        handles.std_mat = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.std_mat(:,:,:,1) = std(handles.MR_PCA_FH,1,4);
        handles.std_mat(:,:,:,2) = std(handles.MR_PCA_AP,1,4);
        handles.std_mat(:,:,:,3) = std(handles.MR_PCA_RL,1,4);

        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);
        
        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));
        

    elseif handles.id_variable==3

        set(handles.pushbutton1,'Visible','on'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','on'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','on'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','on'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','on'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','on'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','on'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','on'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','on'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','on')

        value = handles.thres_noise;
        pp=1/100;
        slider_step(1) = pp;
        slider_step(2) = pp;
        set(handles.slider4,'sliderstep',slider_step,'max',1,'min',0,'Visible','on','value',value);  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of noise masking

        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);

        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));


    elseif handles.id_variable==4

        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','on'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','on'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','on'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','on'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','on'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','on'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

    end


    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


guidata(hObject, handles); 


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, ~, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider3_Callback(hObject, ~, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    pp=1/handles.a;
    slider_step(1) = pp;
    slider_step(2) = pp;
    set(handles.slider3,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.a;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
        handles.slider_value3 = 1;
    else
        handles.slider_value3 = round(handles.slider_value*maxslice);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.text4,'Visible','on','String',['Slice: ',num2str(handles.slider_value3)]); % text panel 1

    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


handles.output = hObject;  
guidata(hObject, handles);  
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, ~, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, ~, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    if handles.id_image == 1

        axes(handles.axes3);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(handles.slider_value3,:,:,2) = double(squeeze(handles.static_Lrgb(handles.slider_value3,:,:,2)).*abs(1-BW1));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.noise_Lrgb(handles.slider_value3,:,:,1) = double(squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,1)).*abs(1-BW1));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif handles.id_image == 2

        axes(handles.axes3);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(handles.slider_value3,:,:,2) = double(squeeze(handles.static_Lrgb(handles.slider_value3,:,:,2)).*abs(1-BW1));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.noise_Lrgb(handles.slider_value3,:,:,1) = double(squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,1)).*abs(1-BW1));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        
    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, ~, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if handles.id_image == 1

        axes(handles.axes3);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(handles.slider_value3,:,:,:) = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(handles.slider_value3,:,:,:) = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif handles.id_image == 2

        axes(handles.axes3);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(handles.slider_value3,:,:,:) = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes3,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(handles.slider_value3,:,:,:) = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   

    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, ~, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.id_variable ==2


        h = msgbox({'Please wait ...','Adjusting Offset Error ...'});

        static = handles.static_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        
        handles.MR_PCA_FH = bface_correction(static,handles.MR_PCA_FH);
        handles.MR_PCA_AP = bface_correction(static,handles.MR_PCA_AP);
        handles.MR_PCA_RL = bface_correction(static,handles.MR_PCA_RL);

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of static tissue
        handles.std_mat = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.std_mat(:,:,:,1) = std(handles.MR_PCA_FH,1,4);
        handles.std_mat(:,:,:,2) = std(handles.MR_PCA_AP,1,4);
        handles.std_mat(:,:,:,3) = std(handles.MR_PCA_RL,1,4);

        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);

        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end


    elseif handles.id_variable ==3

        h = msgbox({'Please wait ...','Removing noise ...'});

        noise = handles.noise_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        handles.MR_FFE_FH = abs(repmat(noise,1,1,1,size(handles.MR_FFE_FH,4))-1).*handles.MR_FFE_FH;
        handles.MR_FFE_AP = abs(repmat(noise,1,1,1,size(handles.MR_FFE_AP,4))-1).*handles.MR_FFE_AP;
        handles.MR_FFE_RL = abs(repmat(noise,1,1,1,size(handles.MR_FFE_RL,4))-1).*handles.MR_FFE_RL;
        handles.MR_PCA_FH = abs(repmat(noise,1,1,1,size(handles.MR_PCA_FH,4))-1).*handles.MR_PCA_FH;
        handles.MR_PCA_AP = abs(repmat(noise,1,1,1,size(handles.MR_PCA_AP,4))-1).*handles.MR_PCA_AP;
        handles.MR_PCA_RL = abs(repmat(noise,1,1,1,size(handles.MR_PCA_RL,4))-1).*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of noise masking

        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);

        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end

    elseif handles.id_variable ==4

        h = msgbox({'Please wait ...','Masking Velocities ...'});

        col = round(handles.points_view1(:,1)./handles.voxel_MR(2));
        row = round(handles.points_view1(:,2)./handles.voxel_MR(2));
        slc = round(handles.points_view2(:,2)./handles.voxel_MR(3));

        mask_roi = handles.mask_roi;
        mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:) = handles.mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:)+1;

        handles.MR_PCA_FH = mask_roi.*handles.MR_PCA_FH;
        handles.MR_PCA_AP = mask_roi.*handles.MR_PCA_AP;
        handles.MR_PCA_RL = mask_roi.*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA       = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end
    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, ~, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view3(1,2),...
                handles.points_view3(1,1),...
                max(handles.points_view3(:,2)) - handles.points_view3(1,2),... 
                max(handles.points_view3(:,1)) - handles.points_view3(1,1)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);

        handles.points_view1 = [ position(2), handles.points_view1(1,2); 
                                 position(2), handles.points_view1(2,2); 
                                 position(2) + position(4), handles.points_view1(3,2); 
                                 position(2) + position(4), handles.points_view1(4,2)];

        handles.points_view2 = [ handles.points_view2(2,1), position(1); 
                                 handles.points_view2(2,1), position(1) + position(3); 
                                 handles.points_view2(3,1), position(1); 
                                 handles.points_view2(4,1), position(1) + position(3)];

        handles.points_view3 = [ position(2), position(1); 
                                 position(2), position(1) + position(3); 
                                 position(2) + position(4), position(1); 
                                 position(2) + position(4), position(1) + position(3)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray



    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view3(1,2),...
                handles.points_view3(1,1),...
                max(handles.points_view3(:,2)) - handles.points_view3(1,2),... 
                max(handles.points_view3(:,1)) - handles.points_view3(1,1)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);

        handles.points_view1 = [ position(2), handles.points_view1(1,2); 
                                 position(2), handles.points_view1(2,2); 
                                 position(2) + position(4), handles.points_view1(3,2); 
                                 position(2) + position(4), handles.points_view1(4,2)];

        handles.points_view2 = [ handles.points_view2(2,1), position(1); 
                                 handles.points_view2(2,1), position(1) + position(3); 
                                 handles.points_view2(3,1), position(1); 
                                 handles.points_view2(4,1), position(1) + position(3)];

        handles.points_view3 = [ position(2), position(1); 
                                 position(2), position(1) + position(3); 
                                 position(2) + position(4), position(1); 
                                 position(2) + position(4), position(1) + position(3)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray



    end


handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on slider movement.
function slider2_Callback(hObject, ~, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    pp=1/handles.b;
    slider_step(1) = pp;
    slider_step(2) = pp;
    set(handles.slider2,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.b;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
        handles.slider_value2 = 1;
    else
        handles.slider_value2 = round(handles.slider_value*maxslice);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.text3,'Visible','on','String',['Slice: ',num2str(handles.slider_value2)]); % text panel 1

    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


handles.output = hObject;  
guidata(hObject, handles);  
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, ~, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, ~, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.id_image == 1

        axes(handles.axes2);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(:,handles.slider_value2,:,2) = double(squeeze(handles.static_Lrgb(:,handles.slider_value2,:,2)).*abs(1-BW1));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.noise_Lrgb(:,handles.slider_value2,:,1) = double(squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,1)).*abs(1-BW1));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    elseif handles.id_image == 2

        axes(handles.axes2);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(:,handles.slider_value2,:,2) = double(squeeze(handles.static_Lrgb(:,handles.slider_value2,:,2)).*abs(1-BW1));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.noise_Lrgb(:,handles.slider_value2,:,1) = double(squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,1)).*abs(1-BW1));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        
    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, ~, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.id_image == 1

        axes(handles.axes2);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(:,handles.slider_value2,:,:) = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(:,handles.slider_value2,:,:) = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif handles.id_image == 2

        axes(handles.axes2);
        plot(0.0)
        k = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(:,handles.slider_value2,:,:) = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes2,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(:,handles.slider_value2,:,:) = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, ~, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
    if handles.id_variable ==2


        h = msgbox({'Please wait ...','Adjusting Offset Error ...'});

        static = handles.static_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        
        handles.MR_PCA_FH = bface_correction(static,handles.MR_PCA_FH);
        handles.MR_PCA_AP = bface_correction(static,handles.MR_PCA_AP);
        handles.MR_PCA_RL = bface_correction(static,handles.MR_PCA_RL);

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of static tissue
        handles.std_mat = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.std_mat(:,:,:,1) = std(handles.MR_PCA_FH,1,4);
        handles.std_mat(:,:,:,2) = std(handles.MR_PCA_AP,1,4);
        handles.std_mat(:,:,:,3) = std(handles.MR_PCA_RL,1,4);

        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);

        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end


    elseif handles.id_variable ==3

        h = msgbox({'Please wait ...','Removing noise ...'});

        noise = handles.noise_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        handles.MR_FFE_FH = abs(repmat(noise,1,1,1,size(handles.MR_FFE_FH,4))-1).*handles.MR_FFE_FH;
        handles.MR_FFE_AP = abs(repmat(noise,1,1,1,size(handles.MR_FFE_AP,4))-1).*handles.MR_FFE_AP;
        handles.MR_FFE_RL = abs(repmat(noise,1,1,1,size(handles.MR_FFE_RL,4))-1).*handles.MR_FFE_RL;
        handles.MR_PCA_FH = abs(repmat(noise,1,1,1,size(handles.MR_PCA_FH,4))-1).*handles.MR_PCA_FH;
        handles.MR_PCA_AP = abs(repmat(noise,1,1,1,size(handles.MR_PCA_AP,4))-1).*handles.MR_PCA_AP;
        handles.MR_PCA_RL = abs(repmat(noise,1,1,1,size(handles.MR_PCA_RL,4))-1).*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of noise masking

        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);

        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end

    elseif handles.id_variable ==4

        h = msgbox({'Please wait ...','Masking Velocities ...'});

        col = round(handles.points_view1(:,1)./handles.voxel_MR(2));
        row = round(handles.points_view1(:,2)./handles.voxel_MR(2));
        slc = round(handles.points_view2(:,2)./handles.voxel_MR(3));

        mask_roi = handles.mask_roi;
        mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:) = handles.mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:)+1;

        handles.MR_PCA_FH = mask_roi.*handles.MR_PCA_FH;
        handles.MR_PCA_AP = mask_roi.*handles.MR_PCA_AP;
        handles.MR_PCA_RL = mask_roi.*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA       = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end
    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, ~, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view2(1,2),...
                handles.points_view2(1,1),...
                max(handles.points_view2(:,2)) - handles.points_view2(1,2),... 
                max(handles.points_view2(:,1)) - handles.points_view2(1,1)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);


        handles.points_view1 = [ handles.points_view1(1,1), position(2); 
                                 handles.points_view1(2,1), position(2) + position(4); 
                                 handles.points_view1(3,1), position(2); 
                                 handles.points_view1(4,1), position(2) + position(4)];

        handles.points_view2 = [ position(2), position(1); 
                                 position(2), position(1) + position(3); 
                                 position(2) + position(4), position(1); 
                                 position(2) + position(4), position(1) + position(3)];

        handles.points_view3 = [ handles.points_view3(1,1), position(1); 
                                 handles.points_view3(2,1), position(1) + position(3); 
                                 handles.points_view3(3,1), position(1); 
                                 handles.points_view3(4,1), position(1) + position(3)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray



    elseif handles.id_image == 2

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view2(1,2),...
                handles.points_view2(1,1),...
                max(handles.points_view2(:,2)) - handles.points_view2(1,2),... 
                max(handles.points_view2(:,1)) - handles.points_view2(1,1)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);


        handles.points_view1 = [ handles.points_view1(1,1), position(2); 
                                 handles.points_view1(2,1), position(2) + position(4); 
                                 handles.points_view1(3,1), position(2); 
                                 handles.points_view1(4,1), position(2) + position(4)];

        handles.points_view2 = [ position(2), position(1); 
                                 position(2), position(1) + position(3); 
                                 position(2) + position(4), position(1); 
                                 position(2) + position(4), position(1) + position(3)];

        handles.points_view3 = [ handles.points_view3(1,1), position(1); 
                                 handles.points_view3(2,1), position(1) + position(3); 
                                 handles.points_view3(3,1), position(1); 
                                 handles.points_view3(4,1), position(1) + position(3)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    pp=1/(handles.c+1);
    slider_step(1) = pp;
    slider_step(2) = pp;
    set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
    maxslice = handles.c;
    handles.slider_value = get(hObject,'Value');
    if handles.slider_value<=pp
        handles.slider_value1 = 1;
    else
        handles.slider_value1 = round(handles.slider_value*maxslice);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.text2,'Visible','on','String',['Slice: ',num2str(handles.slider_value1)]); % text panel 1

    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end

handles.output = hObject;  
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, ~, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if handles.id_image == 1

        axes(handles.axes1);
        plot(0.0)
        k = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(:,:,handles.slider_value1,2) = double((handles.static_Lrgb(:,:,handles.slider_value1,2).*abs(1-BW1)));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            
            handles.noise_Lrgb(:,:,handles.slider_value1,1) = double((handles.noise_Lrgb(:,:,handles.slider_value1,1).*abs(1-BW1)));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif handles.id_image == 2

        axes(handles.axes1);
        plot(0.0)
        k = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            handles.static_Lrgb(:,:,handles.slider_value1,2) = double((handles.static_Lrgb(:,:,handles.slider_value1,2).*abs(1-BW1)));
            handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','g');
            BW1 = createMask(p,k);
            delete(p)

            
            handles.noise_Lrgb(:,:,handles.slider_value1,1) = double((handles.noise_Lrgb(:,:,handles.slider_value1,1).*abs(1-BW1)));
            handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray    

    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.id_image == 1

        axes(handles.axes1);
        plot(0.0)
        k = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(:,:,handles.slider_value1,:) = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(:,:,handles.slider_value1,:) = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
            
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif handles.id_image == 2

        axes(handles.axes1);
        plot(0.0)
        k = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)));
        hold on 
            
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.static_Lrgb(:,:,handles.slider_value1,:) = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

            axis image
            colormap gray
            axis off
            daspect([1 1 1])
            p = drawfreehand(handles.axes1,'color','r');
            BW1 = createMask(p,k);
            BW2 = abs(BW1-1);
            delete(p)

            handles.noise_Lrgb(:,:,handles.slider_value1,:) = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:)).*repmat(BW2,1,1,3)+repmat(BW1,1,1,3);
            handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 
        
        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);
            
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray  

    end

handles.output = hObject;  
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, ~, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if handles.id_variable ==2


        h = msgbox({'Please wait ...','Adjusting Offset Error ...'});

        static = handles.static_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        
        handles.MR_PCA_FH = bface_correction(static,handles.MR_PCA_FH);
        handles.MR_PCA_AP = bface_correction(static,handles.MR_PCA_AP);
        handles.MR_PCA_RL = bface_correction(static,handles.MR_PCA_RL);

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of static tissue
        handles.std_mat = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.std_mat(:,:,:,1) = std(handles.MR_PCA_FH,1,4);
        handles.std_mat(:,:,:,2) = std(handles.MR_PCA_AP,1,4);
        handles.std_mat(:,:,:,3) = std(handles.MR_PCA_RL,1,4);

        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);

        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end


    elseif handles.id_variable ==3

        h = msgbox({'Please wait ...','Removing noise ...'});

        noise = handles.noise_Lrgb_mask.*double(mean(handles.MR_FFE_FH,4)~=0);
        handles.MR_FFE_FH = abs(repmat(noise,1,1,1,size(handles.MR_FFE_FH,4))-1).*handles.MR_FFE_FH;
        handles.MR_FFE_AP = abs(repmat(noise,1,1,1,size(handles.MR_FFE_AP,4))-1).*handles.MR_FFE_AP;
        handles.MR_FFE_RL = abs(repmat(noise,1,1,1,size(handles.MR_FFE_RL,4))-1).*handles.MR_FFE_RL;
        handles.MR_PCA_FH = abs(repmat(noise,1,1,1,size(handles.MR_PCA_FH,4))-1).*handles.MR_PCA_FH;
        handles.MR_PCA_AP = abs(repmat(noise,1,1,1,size(handles.MR_PCA_AP,4))-1).*handles.MR_PCA_AP;
        handles.MR_PCA_RL = abs(repmat(noise,1,1,1,size(handles.MR_PCA_RL,4))-1).*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % masking magnitud
        handles.mask = double(handles.Magnitude>0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % cuantification of noise masking

        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);

        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end

    elseif handles.id_variable ==4

        h = msgbox({'Please wait ...','Masking Velocities ...'});

        col = round(handles.points_view1(:,1)./handles.voxel_MR(2));
        row = round(handles.points_view1(:,2)./handles.voxel_MR(2));
        slc = round(handles.points_view2(:,2)./handles.voxel_MR(3));

        mask_roi = handles.mask_roi;
        mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:) = handles.mask_roi(min(row):max(row), min(col):max(col),min(slc):max(slc),:)+1;

        handles.MR_PCA_FH = mask_roi.*handles.MR_PCA_FH;
        handles.MR_PCA_AP = mask_roi.*handles.MR_PCA_AP;
        handles.MR_PCA_RL = mask_roi.*handles.MR_PCA_RL;

        close(h)
        handles.IPCMRA       = (1./size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*(handles.MR_PCA_FH.^2 + handles.MR_PCA_AP.^2 + handles.MR_PCA_RL.^2),4);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.uipanel4,'visible','off')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.popupmenu2,'Value',1);
        handles.id_variable = get(handles.popupmenu2, 'Value');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(handles.pushbutton1,'Visible','off'); % pushbutton ADD panel 1
        set(handles.pushbutton2,'Visible','off'); % pushbutton DEL panel 1
        set(handles.pushbutton10,'Visible','off'); % pushbutton ROI panel 1
        set(handles.pushbutton3,'Visible','off'); % pushbutton SAVE panel 1

        set(handles.pushbutton4,'Visible','off'); % pushbutton ADD panel 2
        set(handles.pushbutton5,'Visible','off'); % pushbutton DEL panel 2
        set(handles.pushbutton11,'Visible','off'); % pushbutton ROI panel 2
        set(handles.pushbutton6,'Visible','off'); % pushbutton SAVE panel 2

        set(handles.pushbutton7,'Visible','off'); % pushbutton ADD panel 3
        set(handles.pushbutton8,'Visible','off'); % pushbutton DEL panel 3
        set(handles.pushbutton12,'Visible','off'); % pushbutton ROI panel 3
        set(handles.pushbutton9,'Visible','off'); % pushbutton SAVE panel 3

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if handles.id_image == 1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


        elseif handles.id_image == 2

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes1);
            plot(0.0)
            imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
            hold on 

            plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes1,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes2);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes2,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            axes(handles.axes3);
            plot(0.0)
            imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
            hold on 

            plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
            plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

            hold off
            daspect([1,1,1])
            axis image
            set(handles.axes3,'xticklabel',[],'yticklabel',[])
            axis off
            colormap gray

        end
    end

handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, ~, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view1(1,1),...
                handles.points_view1(1,2),...
                max(handles.points_view1(:,1)) - handles.points_view1(1,1),... 
                max(handles.points_view1(:,2)) - handles.points_view1(1,2)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);


        handles.points_view1 = [ position(1), position(2); 
                                 position(1), position(2) + position(4); 
                                 position(1) + position(3), position(2); 
                                 position(1) + position(3), position(2) + position(4)];

        handles.points_view2 = [ position(2), handles.points_view2(1,2); 
                                 position(2), handles.points_view2(2,2); 
                                 position(2) + position(4), handles.points_view2(3,2); 
                                 position(2) + position(4), handles.points_view2(4,2)];

        handles.points_view3 = [ position(1), handles.points_view3(1,2) ; 
                                 position(1), handles.points_view3(2,2); 
                                 position(1) + position(3), handles.points_view3(3,2); 
                                 position(1) + position(3), handles.points_view3(4,2)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray



    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        vect = [handles.points_view1(1,1),...
                handles.points_view1(1,2),...
                max(handles.points_view1(:,1)) - handles.points_view1(1,1),... 
                max(handles.points_view1(:,2)) - handles.points_view1(1,2)];

        h = imrect(gca,vect);
        setColor(h,'c');
        position = wait(h);


        handles.points_view1 = [ position(1), position(2); 
                                 position(1), position(2) + position(4); 
                                 position(1) + position(3), position(2); 
                                 position(1) + position(3), position(2) + position(4)];

        handles.points_view2 = [ position(2), handles.points_view2(1,2); 
                                 position(2), handles.points_view2(2,2); 
                                 position(2) + position(4), handles.points_view2(3,2); 
                                 position(2) + position(4), handles.points_view2(4,2)];

        handles.points_view3 = [ position(1), handles.points_view3(1,2) ; 
                                 position(1), handles.points_view3(2,2); 
                                 position(1) + position(3), handles.points_view3(3,2); 
                                 position(1) + position(3), handles.points_view3(4,2)];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
        plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
        plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
        plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
        plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)
        
        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 
        
        plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
        plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
        plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray



    end



handles.output = hObject;  
guidata(hObject, handles);


% --- Executes on slider movement.
function slider4_Callback(hObject, ~, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if handles.id_variable==2

        pp=1/1000;
        slider_step(1) = pp;
        slider_step(2) = pp;
        set(handles.slider4,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = 1000;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
            handles.slider_value4 = 1;
        else
            handles.slider_value4 = round(handles.slider_value*maxslice);
        end


        handles.thres_static = handles.slider_value4/1000;

        handles.static = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),4]);
        handles.static(:,:,:,1) = double(handles.std_mat(:,:,:,1) < mean(reshape(handles.std_mat(:,:,:,1),1,[]))*handles.thres_static);
        handles.static(:,:,:,2) = double(handles.std_mat(:,:,:,2) < mean(reshape(handles.std_mat(:,:,:,2),1,[]))*handles.thres_static);
        handles.static(:,:,:,3) = double(handles.std_mat(:,:,:,3) < mean(reshape(handles.std_mat(:,:,:,3),1,[]))*handles.thres_static);
        handles.static(:,:,:,4) = double(handles.Magnitude>mean(handles.Magnitude(:)));
        handles.static = double(sum(handles.static,4)==4);

        handles.static_Lrgb = ones(size(handles.static,1),size(handles.static,2),size(handles.static,3),3);
        handles.static_Lrgb(:,:,:,2) = handles.static_Lrgb(:,:,:,2)-(handles.static);
        handles.static_Lrgb = handles.static_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.static_Lrgb_mask = abs(1-handles.static_Lrgb(:,:,:,2));

    elseif handles.id_variable==3


        pp=1/1000;
        slider_step(1) = pp;
        slider_step(2) = pp;
        set(handles.slider4,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = 999;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
            handles.slider_value4 = 1;
        else
            handles.slider_value4 = round(handles.slider_value*maxslice);
        end

        handles.thres_noise = handles.slider_value4/1000;

        handles.noise = zeros([size(handles.MR_PCA_FH,1),size(handles.MR_PCA_FH,2),size(handles.MR_PCA_FH,3),3]);
        handles.noise(:,:,:,1) = double(handles.std_mat(:,:,:,1) > max(reshape(handles.std_mat(:,:,:,1),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,2) = double(handles.std_mat(:,:,:,2) > max(reshape(handles.std_mat(:,:,:,2),1,[]))*(1-handles.thres_noise));
        handles.noise(:,:,:,3) = double(handles.std_mat(:,:,:,3) > max(reshape(handles.std_mat(:,:,:,3),1,[]))*(1-handles.thres_noise));
        handles.noise = double(sum(handles.noise,4)==3);

        handles.noise_Lrgb = ones(size(handles.noise,1),size(handles.noise,2),size(handles.noise,3),3);
        handles.noise_Lrgb(:,:,:,1) = handles.noise_Lrgb(:,:,:,1)-(handles.noise);
        handles.noise_Lrgb = handles.noise_Lrgb.*repmat(handles.mask,1,1,1,3);
        handles.noise_Lrgb_mask = abs(1-handles.noise_Lrgb(:,:,:,1));
    end


    if handles.id_image == 1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,:,handles.slider_value1)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.Magnitude(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.Magnitude(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


    elseif handles.id_image == 2

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes1);
        plot(0.0)
        imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,:,handles.slider_value1)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,:,handles.slider_value1,:));
            himage = imagesc([min(handles.xd(:)),max(handles.xd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view1(1,1),handles.points_view1(2,1)]', [handles.points_view1(1,2),handles.points_view1(2,2)]','y', 'linewidth',2)
            plot([handles.points_view1(1,1),handles.points_view1(3,1)]', [handles.points_view1(1,2),handles.points_view1(3,2)]','y', 'linewidth',2)
            plot([handles.points_view1(2,1),handles.points_view1(4,1)]', [handles.points_view1(2,2),handles.points_view1(4,2)]','y', 'linewidth',2)
            plot([handles.points_view1(3,1),handles.points_view1(4,1)]', [handles.points_view1(3,2),handles.points_view1(4,2)]','y', 'linewidth',2)
        end

        plot([min(handles.xd(:)),max(handles.xd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),[min(handles.yd(:)),max(handles.yd(:))]','r', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes1,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes2);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',squeeze(handles.IPCMRA(:,handles.slider_value2,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(:,handles.slider_value2,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.yd(:)),max(handles.yd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view2(1,2),handles.points_view2(2,2)]', [handles.points_view2(1,1),handles.points_view2(2,1)]','y', 'linewidth',2)
            plot([handles.points_view2(1,2),handles.points_view2(3,2)]', [handles.points_view2(1,1),handles.points_view2(3,1)]','y', 'linewidth',2)
            plot([handles.points_view2(2,2),handles.points_view2(4,2)]', [handles.points_view2(2,1),handles.points_view2(4,1)]','y', 'linewidth',2)
            plot([handles.points_view2(3,2),handles.points_view2(4,2)]', [handles.points_view2(3,1),handles.points_view2(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value3,handles.slider_value3]'-1)*handles.voxel_MR(1),'g', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.yd(:)),max(handles.yd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes2,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        axes(handles.axes3);
        plot(0.0)
        imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',squeeze(handles.IPCMRA(handles.slider_value3,:,:)))
        hold on 

        if handles.id_variable==2
            Lrgb2d = squeeze(handles.static_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==3
            Lrgb2d = squeeze(handles.noise_Lrgb(handles.slider_value3,:,:,:));
            himage = imagesc([min(handles.zd(:)),max(handles.zd(:))]',[min(handles.xd(:)),max(handles.xd(:))]',Lrgb2d);
            cdata = (double(rgb2gray(Lrgb2d))/double(max(max(rgb2gray(Lrgb2d)))))~=1;
            cdata = double(cdata)*0.5;
            set(himage, 'AlphaData', cdata);

        elseif handles.id_variable==4
            plot([handles.points_view3(1,2),handles.points_view3(2,2)]', [handles.points_view3(1,1),handles.points_view3(2,1)]','y', 'linewidth',2)
            plot([handles.points_view3(1,2),handles.points_view3(3,2)]', [handles.points_view3(1,1),handles.points_view3(3,1)]','y', 'linewidth',2)
            plot([handles.points_view3(2,2),handles.points_view3(4,2)]', [handles.points_view3(2,1),handles.points_view3(4,1)]','y', 'linewidth',2)
            plot([handles.points_view3(3,2),handles.points_view3(4,2)]', [handles.points_view3(3,1),handles.points_view3(4,1)]','y', 'linewidth',2)
        end

        plot([min(handles.zd(:)),max(handles.zd(:))]',([handles.slider_value2,handles.slider_value2]'-1)*handles.voxel_MR(2),'r', 'linewidth',1)
        plot(([handles.slider_value1,handles.slider_value1]'-1)*handles.voxel_MR(3),[min(handles.xd(:)),max(handles.xd(:))]','b', 'linewidth',1)

        hold off
        daspect([1,1,1])
        axis image
        set(handles.axes3,'xticklabel',[],'yticklabel',[])
        axis off
        colormap gray

    end


guidata(hObject, handles); 

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, ~, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% change size windows
set(handles.figure1, 'Units', 'pixels');
FigPos = get(handles.figure1, 'Position');
windows_screen_size = get(0,'ScreenSize');

set(handles.figure1, 'Position', [FigPos(1:2),round(windows_screen_size(3)*0.9),round((windows_screen_size(3)*0.9)*0.4)]);
set(handles.figure1, 'Units', 'normalized');

set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.59)
set(handles.popupmenu2,'FontUnits','Normalized','FontSize',0.59)

set(handles.text2,'FontUnits','Normalized','FontSize',0.59)
set(handles.text3,'FontUnits','Normalized','FontSize',0.59)
set(handles.text4,'FontUnits','Normalized','FontSize',0.59)

set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton3,'FontUnits','Normalized','FontSize',0.25)
set(handles.pushbutton4,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton5,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton6,'FontUnits','Normalized','FontSize',0.25)
set(handles.pushbutton7,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton8,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton9,'FontUnits','Normalized','FontSize',0.25)
set(handles.pushbutton10,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton11,'FontUnits','Normalized','FontSize',0.38)
set(handles.pushbutton12,'FontUnits','Normalized','FontSize',0.38)

set(handles.text5,'FontUnits','Normalized','FontSize',0.59)
set(handles.text6,'FontUnits','Normalized','FontSize',0.59)
set(handles.text7,'FontUnits','Normalized','FontSize',0.59)

handles.output = hObject;
guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)
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
