function varargout = GUIDE_THRESHOLDING(varargin)
% GUIDE_THRESHOLDING MATLAB code for GUIDE_THRESHOLDING.fig
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
                   'gui_OpeningFcn', @GUIDE_THRESHOLDING_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_THRESHOLDING_OutputFcn, ...
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
function GUIDE_THRESHOLDING_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
    handles.id = varargin{1}.id;
    if handles.id == 1
        set(handles.uipanel2,'Visible','off')
        handles.IPCMRA = varargin{1}.IPCMRA;
        handles.SEG = varargin{1}.SEG;
        handles.L = varargin{1}.L;
        handles.NUM = varargin{1}.NUM;
        handles.Lrgb = zeros(size(handles.IPCMRA,1),size(handles.IPCMRA,2),size(handles.IPCMRA,3),3);
        handles.min_value_th = varargin{1}.min_value_th;
        handles.max_value_th = varargin{1}.max_value_th;
        handles.xd = varargin{1}.xd;
        handles.yd = varargin{1}.yd;
        handles.zd = varargin{1}.zd;
        handles.a = varargin{1}.a;
        handles.b = varargin{1}.b;
        handles.c = varargin{1}.c;
        handles.slider_axes1 = varargin{1}.slider_axes1;
        handles.slider_axes2 = varargin{1}.slider_axes2;
        handles.slider_axes3 = varargin{1}.slider_axes3;
        handles.id_seg = 0;
        slider_step(1) = 1/255;
        slider_step(2) = 0.1;
        set(handles.slider1,'Value', handles.min_value_th/handles.max_value_th,'sliderstep',slider_step,'max',1,'min',0)
        slider_step(1) = 1/255;
        slider_step(2) = 0.1;
        set(handles.slider2,'Value', handles.max_value_th/handles.max_value_th,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text5,'String',num2str(handles.min_value_th))
        set(handles.text6,'String',num2str(handles.max_value_th))
        set(handles.text8,'String',['from ',num2str(min(handles.IPCMRA(:))),' to ',num2str(max(handles.IPCMRA(:)))])
        axes(handles.axes1);
        imagesc([handles.yd(1,1),handles.yd(end,end)]',[handles.xd(1,1),handles.xd(end,end)]',squeeze(handles.IPCMRA(:,:,handles.slider_axes1)))
        axis image
        daspect([1 1 1])
        colormap gray
        axis off
        slider_step(1) = 1/handles.c;
        slider_step(2) = 0.1;
        set(handles.slider3,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text15,'String',num2str(handles.slider_axes1))
        set(handles.popupmenu1,'Value',1)
        handles.view_sac = 1;
    elseif handles.id == 2
        set(handles.text2,'Visible','off')
        set(handles.text3,'Visible','off')
        set(handles.text4,'Visible','off')
        set(handles.text5,'Visible','off')
        set(handles.text6,'Visible','off')
        set(handles.text7,'Visible','off')
        set(handles.text8,'Visible','off')
        set(handles.slider1,'Visible','off')
        set(handles.slider2,'Visible','off')
        set(handles.pushbutton1,'Visible','off')
        handles.id            = 2;
        handles.IPCMRA        = varargin{1}.IPCMRA;
        handles.SEG           = varargin{1}.SEG;
        handles.L             = varargin{1}.L;
        handles.Lrgb          = varargin{1}.Lrgb;
        handles.NUM           = varargin{1}.NUM;
        handles.xd            = varargin{1}.xd;
        handles.yd            = varargin{1}.yd;
        handles.zd            = varargin{1}.zd;
        handles.a             = varargin{1}.a;
        handles.b             = varargin{1}.b;
        handles.c             = varargin{1}.c;
        handles.slider_axes1  = varargin{1}.slider_axes1;
        handles.slider_axes2  = varargin{1}.slider_axes2;
        handles.slider_axes3  = varargin{1}.slider_axes3;
        handles.voxel_MR      = varargin{1}.voxel_MR;
        handles.id_seg        = 1;
        set(handles.text11,'visible','on')
        set(handles.text11,'String',num2str(handles.NUM))
        set(handles.pushbutton2,'visible','on')

        
        set(handles.uitable1, 'Units', 'pixels');
        AA = get(handles.uitable1,'Position');
        FigPos = get(handles.uitable1, 'Position');
        set(handles.uitable1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
        set(handles.uitable1, 'Units', 'normalized');
        
        set(handles.uitable1,'visible','on','ColumnName',{'Color','# voxels','Selection'},'ColumnEditable',true,'ColumnWidth',{round(AA(3)/3)*0.9 round(AA(3)/3)*0.9 round(AA(3)/3)*0.9})
        handles.data_colors = cell(handles.NUM,3);
        R = squeeze(handles.Lrgb(:,:,:,1));
        G = squeeze(handles.Lrgb(:,:,:,2));
        B = squeeze(handles.Lrgb(:,:,:,3));
        R = R(:);
        G = G(:);
        B = B(:);
        for n=1:handles.NUM
            handles.data_colors{n,1} = ['<html><table border=0 width=100 bgcolor=',rgb2hex([mean(R(handles.L(:)==n)), mean(G(handles.L(:)==n)), mean(B(handles.L(:)==n))]),'><TR><TD>',num2str(n),'</TD></TR> </table></html>'];
            handles.data_colors{n,2} = sum(handles.L(:)==n);
            handles.data_colors{n,3} = false;
            handles.data_colors{n,4} = n;
        end
        handles.data_colors = sortrows(handles.data_colors,2,'descend');
        set(handles.uitable1,'Data',handles.data_colors(:,1:3))
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
        daspect([1 1 1])
        axis image
        colormap gray
        axis off
        slider_step(1) = 1/handles.c;
        slider_step(2) = 0.1;
        set(handles.slider3,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
        set(handles.text15,'String',num2str(handles.slider_axes1))
        set(handles.popupmenu1,'Value',1)
        handles.view_sac = 1;
    end
guidata(hObject, handles);
uiwait(handles.figure1);
function varargout = GUIDE_THRESHOLDING_OutputFcn(hObject, eventdata, handles)
if handles.id == 1;
    varargout{1} = handles.output;
    SEG = handles.SEG;
    setappdata(0,'SEG',SEG);
    L = handles.L;
    setappdata(0,'L',L);
    NUM = handles.NUM;
    setappdata(0,'NUM',NUM);
    Lrgb = handles.Lrgb;
    setappdata(0,'Lrgb',Lrgb);
    min_value_th = handles.min_value_th;
    setappdata(0,'min_value_th',min_value_th);
    max_value_th = handles.max_value_th;
    setappdata(0,'max_value_th',max_value_th);
    delete(handles.figure1);
elseif handles.id == 2;
    varargout{1} = handles.output;
    SEG = handles.SEG;
    setappdata(0,'SEG',SEG);
    L = handles.L;
    setappdata(0,'L',L);
    NUM = handles.NUM;
    setappdata(0,'NUM',NUM);
    Lrgb = handles.Lrgb;
    setappdata(0,'Lrgb',Lrgb);
    delete(handles.figure1);
end
function slider1_Callback(hObject, eventdata, handles)
        pp=1/255;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider1,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = 255;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.min_value_th = handles.slider_value;
        set(handles.text5,'String',num2str(handles.slider_value));
handles.output = hObject;
guidata(hObject, handles);
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_Callback(hObject, eventdata, handles)
        pp=1/255;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider2,'sliderstep',slider_step,'max',1,'min',0)
        maxslice = 255;
        handles.slider_value = get(hObject,'Value');
        if handles.slider_value<=pp
                handles.slider_value = 1;
        else
                handles.slider_value = round(handles.slider_value*maxslice);
        end
        handles.max_value_th = handles.slider_value;
        set(handles.text6,'String',num2str(handles.slider_value));
handles.output = hObject;
guidata(hObject, handles);
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function pushbutton1_Callback(hObject, eventdata, handles)
        SEG1 = handles.IPCMRA > handles.min_value_th;
        if handles.max_value_th==255
            SEG2 = ones(size(handles.IPCMRA));
        else
            SEG2 = handles.IPCMRA < handles.max_value_th;
        end
        SEG = SEG1 + SEG2;
        handles.SEG = double(SEG==2);
        [L,NUM] = bwlabeln(handles.SEG,6);
        handles.L = L;
        handles.NUM = NUM;
        handles.Lrgb = label2rgb3d(handles.L,'jet');
        set(handles.text11,'visible','on')
        set(handles.text11,'String',num2str(handles.NUM))
        set(handles.pushbutton2,'visible','on')

        
        
        set(handles.uitable1, 'Units', 'pixels');
        AA = get(handles.uitable1,'Position');
        FigPos = get(handles.uitable1, 'Position');
        set(handles.uitable1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
        set(handles.uitable1, 'Units', 'normalized');
        
        
        set(handles.uitable1,'visible','on','ColumnName',{'Color','# voxels','Selection'},'ColumnEditable',true,'ColumnWidth',{round(AA(3)/3)*0.9 round(AA(3)/3)*0.9 round(AA(3)/3)*0.9})
        handles.data_colors = cell(NUM,3);
        R = squeeze(handles.Lrgb(:,:,:,1));
        G = squeeze(handles.Lrgb(:,:,:,2));
        B = squeeze(handles.Lrgb(:,:,:,3));
        R = R(:);
        G = G(:);
        B = B(:);
        for n=1:NUM
            handles.data_colors{n,1} = ['<html><table border=0 width=100 bgcolor=',rgb2hex([mean(R(handles.L(:)==n)), mean(G(handles.L(:)==n)), mean(B(handles.L(:)==n))]),'><TR><TD>',num2str(n),'</TD></TR> </table></html>'];
            handles.data_colors{n,2} = sum(handles.L(:)==n);
            handles.data_colors{n,3} = false;
            handles.data_colors{n,4} = n;
        end
        handles.data_colors = sortrows(handles.data_colors,2,'descend');
        set(handles.uitable1,'Data',handles.data_colors(:,1:3))
        handles.id_seg = 1;
        if handles.view_sac == 1;
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
            set(handles.text15,'String',num2str(handles.slider_axes1))
        elseif handles.view_sac == 2;
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
            set(handles.text15,'String',num2str(handles.slider_axes2))
        elseif handles.view_sac == 3;
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
            set(handles.text15,'String',num2str(handles.slider_axes3))
        end
        set(handles.uipanel2,'Visible','on')
handles.output = hObject;
guidata(hObject, handles);
function text5_CreateFcn(hObject, eventdata, handles)
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end
function pushbutton2_Callback(hObject, eventdata, handles)
    tableData = get(handles.uitable1, 'data');
    MAT = zeros(size(tableData,1),2);
    for n=1:size(tableData,1)
       MAT(n,:) = [tableData{n,2},tableData{n,3}];
    end
    [r,c,v] = find(MAT(:,2)==1);
    SEG_new = zeros(size(handles.IPCMRA));
    L_new = zeros(size(handles.IPCMRA));
    for n=1:length(r)
        r_id = handles.data_colors(r(n),4);
        SEG_new = SEG_new + double((handles.L == r_id{1}));
        L_new = L_new + double((handles.L == r_id{1}))*r_id{1};
    end
    handles.Lrgb = handles.Lrgb.*double(repmat(abs(SEG_new),1,1,1,3)) + double(repmat(abs(SEG_new-1),1,1,1,3));
    handles.SEG = SEG_new;
    handles.L = L_new;
    handles.NUM = length(r);
    if handles.view_sac == 1;
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
        set(handles.text15,'String',num2str(handles.slider_axes1))
    elseif handles.view_sac == 2;
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
        set(handles.text15,'String',num2str(handles.slider_axes2))
    elseif handles.view_sac == 3;
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
        set(handles.text15,'String',num2str(handles.slider_axes3))
    end
handles.output = hObject;
guidata(hObject, handles);
function slider3_Callback(hObject, eventdata, handles)
    if handles.view_sac == 1;
        pp=1/handles.c;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider3,'sliderstep',slider_step,'max',1,'min',0)
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
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        set(handles.text15,'String',num2str(handles.slider_axes1))
    elseif handles.view_sac == 2;
        pp=1/handles.b;
        slider_step(1) = pp;
        slider_step(2) = 0.1;
        set(handles.slider3,'sliderstep',slider_step,'max',1,'min',0)
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
        hold off
        axis image
        colormap gray
        axis off
        daspect([1 1 1])
        set(handles.text15,'String',num2str(handles.slider_axes2))
    elseif handles.view_sac == 3;
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
        set(handles.text15,'String',num2str(handles.slider_axes3))
    end
handles.output = hObject;
guidata(hObject, handles);
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
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
        set(handles.slider3,'Value', handles.slider_axes1/handles.c,'sliderstep',slider_step,'max',1,'min',0)
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
        set(handles.slider3,'Value', handles.slider_axes2/handles.b,'sliderstep',slider_step,'max',1,'min',0)
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
        set(handles.slider3,'Value', handles.slider_axes3/handles.a,'sliderstep',slider_step,'max',1,'min',0)
end
handles.output = hObject;
guidata(hObject, handles);
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    set(handles.uitable1, 'Units', 'pixels');
    AA = get(handles.uitable1,'Position');
    FigPos = get(handles.uitable1, 'Position');
    set(handles.uitable1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.uitable1, 'Units', 'normalized');

    set(handles.uitable1,'visible','on','ColumnName',{'Color','# voxels','Selection'},'ColumnEditable',true,'ColumnWidth',{round(AA(3)/3)*0.9 round(AA(3)/3)*0.9 round(AA(3)/3)*0.9})
        

    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');
    
    set(handles.text2,'FontUnits','Normalized','FontSize',0.56)
    set(handles.text3,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text4,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text5,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text6,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text7,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text8,'FontUnits','Normalized','FontSize',0.60)
    set(handles.text10,'FontUnits','Normalized','FontSize', 0.60)
    set(handles.text11,'FontUnits','Normalized','FontSize', 0.60)
    set(handles.text12,'FontUnits','Normalized','FontSize', 0.60)
    set(handles.text14,'FontUnits','Normalized','FontSize', 0.60)
    set(handles.text15,'FontUnits','Normalized','FontSize', 0.60)
    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.4)
    set(handles.uitable1,'FontUnits','Normalized','FontSize',0.03)
    set(handles.popupmenu1,'FontUnits','Normalized','FontSize',0.33)

handles.output = hObject;
guidata(hObject, handles);
