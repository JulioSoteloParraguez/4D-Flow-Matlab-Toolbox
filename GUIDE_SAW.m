function varargout = GUIDE_SAW(varargin)
% GUIDE_SAW MATLAB code for GUIDE_SAW.fig
%      GUIDE_SAW, by itself, creates a new GUIDE_SAW or raises the existing
%      singleton*.
%
%      H = GUIDE_SAW returns the handle to a new GUIDE_SAW or the handle to
%      the existing singleton*.
%
%      GUIDE_SAW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_SAW.M with the given input arguments.
%
%      GUIDE_SAW('Property','Value',...) creates a new GUIDE_SAW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_SAW_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_SAW_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_SAW

% Last Modified by GUIDE v2.5 24-Nov-2022 23:57:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUIDE_SAW_OpeningFcn, ...
    'gui_OutputFcn',  @GUIDE_SAW_OutputFcn, ...
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


% --- Executes just before GUIDE_SAW is made visible.
function GUIDE_SAW_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUIDE_SAW (see VARARGIN)

% Choose default command line output for GUIDE_SAW
handles.output = hObject;
path(path,'SAW/')

handles.VENC        = varargin{1}.VENC;
handles.voxel_MR    = varargin{1}.voxel_MR;
handles.heart_rate  = varargin{1}.heart_rate;
handles.MR_FFE_FH   = varargin{1}.MR_FFE_FH;
handles.MR_FFE_AP   = varargin{1}.MR_FFE_AP;
handles.MR_FFE_RL   = varargin{1}.MR_FFE_RL;
handles.MR_PCA_FH   = varargin{1}.MR_PCA_FH;
handles.MR_PCA_AP   = varargin{1}.MR_PCA_AP;
handles.MR_PCA_RL   = varargin{1}.MR_PCA_RL;
handles.SEG         = varargin{1}.SEG;
handles.Centerline  = varargin{1}.Centerline;

% for point location
handles.faces       = varargin{1}.faces;
handles.nodes       = varargin{1}.nodes;
handles.elem        = varargin{1}.elem;
handles.Laplace     = varargin{1}.Laplace;
handles.IPCMRA      = (1/size(handles.MR_FFE_FH,4)).*sum((handles.MR_FFE_FH.^2).*sqrt((handles.MR_PCA_FH.^2)+(handles.MR_PCA_AP.^2)+(handles.MR_PCA_RL.^2)),4);
% handles.IPCMRA      = histeq((handles.IPCMRA/max(handles.IPCMRA(:))));
handles.IPCMRA      = (handles.IPCMRA/max(handles.IPCMRA(:)));
for n=1:size(handles.IPCMRA,3)
    handles.IPCMRA(:,:,n) = imadjust(handles.IPCMRA(:,:,n),[0,0.3],[0,1]);
end
handles.IPCMRA      = handles.IPCMRA*255;

[handles.X_Pos, handles.Y_Pos, handles.Z_Pos] = meshgrid(0:size(handles.IPCMRA,2)-1,0:size(handles.IPCMRA,1)-1,0:size(handles.IPCMRA,3)-1);
handles.X_Pos = handles.X_Pos*handles.voxel_MR(2);
handles.Y_Pos = handles.Y_Pos*handles.voxel_MR(1);
handles.Z_Pos = handles.Z_Pos*handles.voxel_MR(3);

% if we execute the plane selected in the popupmenu 1 case 5, we replace
% the variables for handles.idx
handles.id_case_5 = 0;

%% working with the data information
handles.PM                  = PressureMappingClass;
handles.PM.Vel(1,:,:,:,:)   =  handles.MR_PCA_AP;
handles.PM.Vel(2,:,:,:,:)   = -handles.MR_PCA_FH;
handles.PM.Vel(3,:,:,:,:)   = -handles.MR_PCA_RL;
handles.PM.Vel              = permute( handles.PM.Vel, [ 5 1 2 3 4 ] );
handles.PM.Vel              = handles.PM.Vel*10; % From [ cm ] to [ mm ]
handles.PM.Venc             = handles.VENC*10; % From [ cm ] to [ mm ]
handles.PM.nFrames          = size( handles.PM.Vel, 1 );


%% Calculating time step 'dt' from heart rate 'data.heart_rate'
mpb             = 1/handles.heart_rate; % Time duration of a heart beat in minutes
spb             = 60*mpb;            % Time duration of a heart beat in seconds
% Time between frames, assuming that the total number of frames covers a full herat beat
handles.PM.dt   = spb / ( handles.PM.nFrames - 1 );
handles.PM.Time = handles.PM.dt*( 0 : handles.PM.nFrames-1 ); % Time vector of the acquisition

%% working with the segmentation information
handles.PM.Domain               = handles.SEG( 2:end-1, 2:end-1, 2:end-1 );
handles.PM.hd.spacing           = handles.voxel_MR;
handles.PM.hd.origin            = [ 0 0 0 ];
handles.PM.hd.TransformMatrix   = eye( 3 );
handles.PM.hd.Mv2w              = [ diag( handles.PM.hd.spacing ) [ 0 0 0 ]' ];
handles.PM.hd.Mw2v              = [ diag( 1./handles.PM.hd.spacing ) [ 0 0 0 ]' ];
handles.PM                      = handles.PM.reduceDomain( 3 );
handles.PM.Vel                  = handles.PM.Vel( :, :, handles.PM.hd.subDom.r, handles.PM.hd.subDom.c, handles.PM.hd.subDom.s );
handles.PM                      = handles.PM.calculatePeakFrame;


%% working with centerline
handles.clPoints            = handles.Centerline; % [ mm ]
% Centering the centerline with respect to the original binary domain
handles.clPoints            = handles.clPoints - handles.PM.hd.spacing;
% Registering the centerline with respect to the new reduced binary domain
handles.clPoints            = handles.clPoints - handles.PM.hd.spacing.*[ handles.PM.hd.subDom.r(1)-1 handles.PM.hd.subDom.c(1)-1 handles.PM.hd.subDom.s(1)-1 ];
% Calculating the slopes
handles.clSlopes            = diff( handles.clPoints, 1 )./vecnorm( diff( handles.clPoints, 1 ), 2, 2 );
handles.clSlopes(end+1,:)   = handles.clSlopes(end,:);

%%
% Length of the centerline at each centerline point
d                   = [ 0 cumsum( vecnorm( diff( handles.clPoints, 1 ), 2, 2 )' ) ];
% Number of cross sections based on the voxel dimensions
handles.PM.nPlanes          = ceil( d( end )/min( handles.PM.hd.spacing ) ) + 1;
handles.PM.PlanesEveryXmm   = d( end ) / ( handles.PM.nPlanes - 1 );

di                  = ( 0 : handles.PM.nPlanes - 1 )*handles.PM.PlanesEveryXmm;

disp([min(d),max(d)])
disp([min(di),max(di)])
disp(size(handles.clPoints))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%handles.clPointsi           = interpn( d, handles.clPointsi, di, 'linear' )';

handles.clPointsi                = zeros(length(di),3);
handles.clPointsi(:,1)           = interp1( d, handles.clPoints(:,1), di, 'linear' )';
handles.clPointsi(:,2)           = interp1( d, handles.clPoints(:,2), di, 'linear' )';
handles.clPointsi(:,3)           = interp1( d, handles.clPoints(:,3), di, 'linear' )';

%handles.clSlopesi           = interpn( d, handles.clSlopes, di, 'linear' )';

handles.clSlopesi                = zeros(length(di),3);
handles.clSlopesi(:,1)           = interp1( d, handles.clSlopes(:,1), di, 'linear' )';
handles.clSlopesi(:,2)           = interp1( d, handles.clSlopes(:,2), di, 'linear' )';
handles.clSlopesi(:,3)           = interp1( d, handles.clSlopes(:,3), di, 'linear' )';

handles.clSlopesi           = normr( handles.clSlopesi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% working with all planes
handles.PM.b_Clin = 1;
opts.bAllPlanes = 1;
opts.bRecalculateBplanes = 1;

handles.PM = handles.PM.buildAllPlanes( handles.clPointsi, handles.clSlopesi );
% Calculating the binary mask of the cross sections
handles.PM = handles.PM.interpolateMaskAtPlanes( opts );
% Calculating the velocity field for each cross section
opts.bRecalculateV = 1;
for iFrame = 1 : handles.PM.nFrames
    handles.PM = handles.PM.intersect_velocity_with_planes( iFrame, opts );
end

%% Calculating radius, curvature, flow Q, pressure difference SAW, Time-To-Foot, Pulse Wave Velocity, stiffness E
%  Metrics saved as handles.PM.r, handles.PM.c, handles.PM.Q, handles.PM.SAW, handles.PM.TTF, handles.PM.PWV, handles.PM.E
formulation             = 'all'; % formulation = 'all' ==> SAW, fSAW, backSAW, SB
handles.PM              = handles.PM.ComputeSAW( formulation, 1:handles.PM.nFrames); % SAW at each frame, also Q is calculated
handles.PM.c            = [ handles.PM.AllPlanes.curv ]'; % Curvature  at each point of the centerline in [ 1/m ]
handles.PM              = handles.PM.ComputeTTF( 'fromQ' ); % Time-To-Foot, needed to compute PWV
pwvOpts.dtPcent         = 1;
pwvOpts.pointsAfterDt   = '1st';
handles.PM              = handles.PM.ComputePWV( pwvOpts ); % PWV, needed to compute E (elasticity, stiffness)
handles.PM              = handles.PM.ComputeE;

%% Plotting the 3D velocity field at peak systole
axes(handles.axes1)
plot(0,0)
[~] = handles.PM.plotDomain;
hold on
iFrame = handles.PM.peakFrame;
options.frame = iFrame;
options.scale = 5;
handles.PM.plotVel( options );
hold off
grid on
axis on
ylabel('F-H')
xlabel('A-P')
zlabel('R-L')
view(142,39)
daspect([1 1 1])
axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');

list_string1 = {'Segmentation + 3D velocity field (Peak Systole)',...
    'Segmentation + 3D points of the centerline',...
    'Segmentation + 3D velocity field (Peak Systole) + 2D Planes',...
    'Segmentation + Streamlines of Velocity field (Peak Systole)',...
    'Segmentation + 2D Planes locations',...
    'Flow at peak frame along the vessel [l/min]',...
    'SAW at peak frame along the vessel [mmHg]',...
    'Pulse Wave Velocity along the vessel [m/s]',...
    'Elastic Modulus along the vessel [kPa]',...
    'Diameter along the vessel [mm]',...
    'Curvature along the vessel [1/m]'};
set(handles.popupmenu1,'visible','on','String',list_string1,'value',1);
set(handles.pushbutton1,'visible','on');
set(handles.pushbutton2,'visible','off');

handles.pressure_recovery_id = 0;
handles.arch_planes_selection_id = 0;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes GUIDE_SAW wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_SAW_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popuhandles.PMenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popuhandles.PMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popuhandles.PMenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popuhandles.PMenu1
switch get(handles.popupmenu1,'Value')
    case 1
        axes(handles.axes1)
        plot(0,0)
        [~] = handles.PM.plotDomain;
        hold on
        iFrame = handles.PM.peakFrame;
        options.frame = iFrame;
        options.scale = 5;
        handles.PM.plotVel( options );
        hold off
        grid on
        axis on
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');

        set(handles.pushbutton2,'visible','off');
        
    case 2
        axes(handles.axes1)
        plot(0,0)
        handles.PM.plotDomain();
        hold on
        plot3( handles.clPoints(:,1), handles.clPoints(:,2), handles.clPoints(:,3), 'r-' )
        plot3( handles.clPointsi(:,1), handles.clPointsi(:,2), handles.clPointsi(:,3), 'bo' )
        hold off
        grid on
        axis on
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');

        set(handles.pushbutton2,'visible','off');
        
    case 3
        axes(handles.axes1)
        plot(0,0)
        handles.PM.plotDomain();
        hold on
        iFrame = handles.PM.peakFrame;
        dFF = 6; % Across-planes decimation factor
        planes = 1 : dFF : handles.PM.nPlanes;
        handles.PM.plotPlanes( planes, iFrame );
        hold off
        grid on
        axis on
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');

        set(handles.pushbutton2,'visible','off');
        
    case 4
        iFrame = handles.PM.peakFrame; % Chosen time-frame
        dF  = 3; % In-plane decimation factor
        dFF = 6; % Across-planes decimation factor
        planes = 1 : dFF : handles.PM.nPlanes;
        
        axes(handles.axes1)
        plot(0,0)
        axis off
        handles.PM.plotStreamlines( planes, iFrame, dF );
        grid on
        axis on
        axis tight
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');

        set(handles.pushbutton2,'visible','off');
        
    case 5 

        set(handles.pushbutton2,'visible','on','string','ARCH PLANES RE-LOCATION');
        handles.id_case_5 = 1;

        if handles.arch_planes_selection_id == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
        else
            handles.idx = [1,handles.idx1, handles.idx2, handles.idx3, handles.PM.nPlanes];
        end
        
        axes(handles.axes1)
        plot(0,0)
        handles.PM.plotDomain();
        hold on
        plot3( handles.clPointsi(:,1), handles.clPointsi(:,2), handles.clPointsi(:,3), 'r-' )
        plot3( handles.clPointsi(handles.idx,1), handles.clPointsi(handles.idx,2), handles.clPointsi(handles.idx,3), '*g' )
        optPlotPlane.bPlotVelocityVectors = 0;
        optPlotPlane.bPlotKeyPoint = 0;
        optPlotPlane.bPlotKeySlope = 0;
        P2plot = handles.PM.AllPlanes( handles.idx );
        for iL = 1 : size( handles.idx, 2 )
            P = P2plot(iL);
            handles.PM.plotPlane( P, handles.PM.peakFrame, optPlotPlane )
        end
        hold off
        grid on
        axis on
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');
        
    case 6
        if handles.id_case_5 == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
            handles.idx(end) = handles.PM.nPlanes;
        end
        
        dist = [ handles.PM.AllPlanes.dist ];
        D = handles.PM.Q( :, handles.PM.peakFrame )*60;
        yUnit = '[l/m]';
        titleStr = 'Flow\nat peak frame';
        D = smooth( D, length( D )/6,'rlowess');
        
        axes(handles.axes1)
        plot(0,0)
        plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
        hold on
        ax = gca;
        axis( ax, 'tight' )
        for j = 1 : 2 : size( handles.idx, 2 )-1
            patch( 'Faces',[1 2 3 4], 'Vertices', ...
                [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
        end
        
        xlabel( 'Aortic Length [mm]' )
        ylabel( yUnit )
        ax.Title.String = sprintf( titleStr );
        
        ax.FontName = 'Calibri';
        ax.FontSize = 14;
        ax.FontWeight = 'normal';
        axis square
        
        ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
        hold off
        
        axtoolbar(handles.axes1,'Visible','off');
        set(handles.pushbutton2,'visible','off');
        
        
    case 7
        set(handles.pushbutton2,'visible','on','string','COMPUTE PRESSURE RECOVERY');

        if handles.pressure_recovery_id == 0

            if handles.id_case_5 == 0
                handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
                handles.idx(end) = handles.PM.nPlanes;
            end
            
            dist = [ handles.PM.AllPlanes.dist ];
            D = handles.PM.SAW( :, handles.PM.peakFrame );
            yUnit = '[mmHg]';
            titleStr = 'SAW\nat peak frame';
            D = smooth( D, length( D )/6,'rlowess');
            
            axes(handles.axes1)
            plot(0,0)
            plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
            hold on
            ax = gca;
            axis( ax, 'tight' )
            for j = 1 : 2 : size( handles.idx, 2 )-1
                patch( 'Faces',[1 2 3 4], 'Vertices', ...
                    [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                    [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
            end
            
            xlabel( 'Aortic Length [mm]' )
            ylabel( yUnit )
            ax.Title.String = sprintf( titleStr );
            
            ax.FontName = 'Calibri';
            ax.FontSize = 14;
            ax.FontWeight = 'normal';
            axis square
            
            ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
            hold off
            axtoolbar(handles.axes1,'Visible','off');
    
        else
            
            % update the figure when the windows is closed
            if handles.id_case_5 == 0
                handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
                handles.idx(end) = handles.PM.nPlanes;
            end
            
            dist = [ handles.PM.AllPlanes.dist ];
            D = handles.PM.SAW( :, handles.PM.peakFrame );
            yUnit = '[mmHg]';
            titleStr = 'SAW\nat peak frame';
            D = smooth( D, length( D )/6,'rlowess');
            
            axes(handles.axes1)
            plot(0,0)
            plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
            hold on
            ax = gca;
            axis( ax, 'tight' )
            for j = 1 : 2 : size( handles.idx, 2 )-1
                patch( 'Faces',[1 2 3 4], 'Vertices', ...
                    [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                    [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
            end
            
            xlabel( 'Aortic Length [mm]' )
            ylabel( yUnit )
            ax.Title.String = sprintf( titleStr );
            
            ax.FontName = 'Calibri';
            ax.FontSize = 14;
            ax.FontWeight = 'normal';
            axis square
            
            ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
            hold off
            axtoolbar(handles.axes1,'Visible','off');
            
            [~,EOAIdx] = max(D);   
            handles.EOAdist = dist(EOAIdx);
            xline(handles.EOAdist, 'g-.','LineWidth', 2 );
            xline(handles.EOAdist + handles.PM.PrecDist, 'b-.','LineWidth', 2 );

            text(handles.EOAdist + handles.PM.PrecDist, ...
                ((max(D( handles.idx(1):handles.idx(end) )) + min(D( handles.idx(1):handles.idx(end) )))/2),...
                ['\leftarrow PR Distance = ', num2str(handles.PM.PrecDist,4), ' [mm]'] ,'FontSize',14)

        end
        
    case 8
        
        if handles.id_case_5 == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
            handles.idx(end) = handles.PM.nPlanes;
        end
        
        dist = [ handles.PM.AllPlanes.dist ];
        D = handles.PM.PWV;
        yUnit = 'PWV [m/s]';
        titleStr = 'Pulse Wave Velocity';
        D = smooth( D, length( D )/6,'rlowess');
        
        axes(handles.axes1)
        plot(0,0)
        plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
        hold on
        ax = gca;
        axis( ax, 'tight' )
        for j = 1 : 2 : size( handles.idx, 2 )-1
            patch( 'Faces',[1 2 3 4], 'Vertices', ...
                [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
        end
        
        xlabel( 'Aortic Length [mm]' )
        ylabel( yUnit )
        ax.Title.String = sprintf( titleStr );
        
        ax.FontName = 'Calibri';
        ax.FontSize = 14;
        ax.FontWeight = 'normal';
        axis square
        
        ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
        hold off
        axtoolbar(handles.axes1,'Visible','off');

        set(handles.pushbutton2,'visible','off');
        
    case 9
        if handles.id_case_5 == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
            handles.idx(end) = handles.PM.nPlanes;
        end
        
        dist = [ handles.PM.AllPlanes.dist ];
        D = handles.PM.E;
        yUnit = 'E [kPa]';
        titleStr = 'Elastic Modulus';
        D = smooth( D, length( D )/6,'rlowess');
        
        axes(handles.axes1)
        plot(0,0)
        plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
        hold on
        ax = gca;
        axis( ax, 'tight' )
        for j = 1 : 2 : size( handles.idx, 2 )-1
            patch( 'Faces',[1 2 3 4], 'Vertices', ...
                [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
        end
        
        xlabel( 'Aortic Length [mm]' )
        ylabel( yUnit )
        ax.Title.String = sprintf( titleStr );
        
        ax.FontName = 'Calibri';
        ax.FontSize = 14;
        ax.FontWeight = 'normal';
        axis square
        
        ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
        hold off
        axtoolbar(handles.axes1,'Visible','off');

        set(handles.pushbutton2,'visible','off');
        
    case 10
        if handles.id_case_5 == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
            handles.idx(end) = handles.PM.nPlanes;
        end
        
        dist = [ handles.PM.AllPlanes.dist ];
        D = handles.PM.r * 2;
        yUnit = 'D [mm]';
        titleStr = 'Diameter';
        D = smooth( D, length( D )/6,'rlowess');
        
        axes(handles.axes1)
        plot(0,0)
        plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
        hold on
        ax = gca;
        axis( ax, 'tight' )
        for j = 1 : 2 : size( handles.idx, 2 )-1
            patch( 'Faces',[1 2 3 4], 'Vertices', ...
                [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
        end
        
        xlabel( 'Aortic Length [mm]' )
        ylabel( yUnit )
        ax.Title.String = sprintf( titleStr );
        
        ax.FontName = 'Calibri';
        ax.FontSize = 14;
        ax.FontWeight = 'normal';
        axis square
        
        ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
        hold off
        axtoolbar(handles.axes1,'Visible','off');

        set(handles.pushbutton2,'visible','off');
        
    case 11
        if handles.id_case_5 == 0
            handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
            handles.idx(end) = handles.PM.nPlanes;
        end
        
        dist = [ handles.PM.AllPlanes.dist ];
        D = [ handles.PM.AllPlanes.curv ];
        yUnit = '[1/m]';
        titleStr = 'Curvature';
        D = smooth( D, length( D )/6,'rlowess');
        
        axes(handles.axes1)
        plot(0,0)
        plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
        hold on
        ax = gca;
        axis( ax, 'tight' )
        for j = 1 : 2 : size( handles.idx, 2 )-1
            patch( 'Faces',[1 2 3 4], 'Vertices', ...
                [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
        end
        
        xlabel( 'Aortic Length [mm]' )
        ylabel( yUnit )
        ax.Title.String = sprintf( titleStr );
        
        ax.FontName = 'Calibri';
        ax.FontSize = 14;
        ax.FontWeight = 'normal';
        axis square
        
        ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
        hold off
        axtoolbar(handles.axes1,'Visible','off');

        set(handles.pushbutton2,'visible','off');
        
        
        
end

handles.output = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popuhandles.PMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popuhandles.PMenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

directory = uigetdir(pwd, 'Select Directory');

PM = handles.PM;

if ~isfield( handles, 'idx' )
    %  Landmarks
    idx = 1 : ceil(PM.nPlanes/4) : PM.nPlanes;
    idx(end+1) = PM.nPlanes;
else
    idx = handles.idx;
end

save([directory,'\PM.mat'],'PM')

ReportGen( PM, idx, directory, handles.pressure_recovery_id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist = [ PM.AllPlanes.dist ];
AL = [dist( idx(1):idx(end) )-dist( idx(1) )]';
FP = PM.Q( :, PM.peakFrame )*60;
FP = smooth( FP, length( FP )/6,'rlowess');

SP = PM.SAW( :, PM.peakFrame );
SP = smooth( SP, length( SP )/6,'rlowess');

PW = PM.PWV;
PW = smooth( PW, length( PW )/6,'rlowess');

EL = PM.E;
EL = smooth( EL, length( EL )/6,'rlowess');

DI = PM.r * 2;
DI = smooth( DI, length( DI )/6,'rlowess');

CU = [PM.AllPlanes.curv]';
CU = smooth( CU, length( CU )/6,'rlowess');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% table data
if handles.arch_planes_selection_id == 1
    handles.idx = [1,handles.idx1, handles.idx2, handles.idx3, length(AL)];
    sections = zeros(size(AL));
    for n =1:length(handles.idx)-1
        sections(handles.idx(n):handles.idx(n+1))= n;
    end
    if handles.pressure_recovery_id == 1
        EOAdist = handles.EOAdist;
        EOAdist_vect = zeros(size(AL)) + EOAdist;

        PR = ones(size(AL))*handles.PM.PrecDist;

        VV = {'Aortic Length [mm]';...
              'Flow at peak frame [l/min]';...
              'SAW at peak frame [mmHg]';...
              'Pulse Wave Velocity [m/s]';...
              'Elastic Modulus [KPa]';...
              'Diameter [mm]';...
              'Curvature [1/m]';...
              'EOAdist [mm]';...
              'Pressure Recovery Distance [mm]';...
              'Sections [-]'};
    
        MM = [AL';FP';SP';PW';EL';DI';CU';EOAdist_vect';PR';sections'];
        TT = table(VV,MM);
        TT.Properties.VariableNames = {' ','Point'};
    
    else 
    
        VV = {'Aortic Length [mm]';...
              'Flow at peak frame [l/min]';...
              'SAW at peak frame [mmHg]';...
              'Pulse Wave Velocity [m/s]';...
              'Elastic Modulus [KPa]';...
              'Diameter [mm]';...
              'Curvature [1/m]';...
              'Sections [-]'};
        
        MM = [AL';FP';SP';PW';EL';DI';CU';sections'];
        TT = table(VV,MM);
        TT.Properties.VariableNames = {' ','Point'};
    
    end

else
    handles.idx = 1 : floor(handles.PM.nPlanes/4) : length(AL);
    sections = zeros(size(AL));
    for n =1:length(handles.idx)-1
        sections(handles.idx(n):handles.idx(n+1))= n;
    end
    if handles.pressure_recovery_id == 1
        EOAdist = handles.EOAdist;
        EOAdist_vect = zeros(size(AL)) + EOAdist;

        PR = ones(size(AL))*handles.PM.PrecDist;

        VV = {'Aortic Length [mm]';...
              'Flow at peak frame [l/min]';...
              'SAW at peak frame [mmHg]';...
              'Pulse Wave Velocity [m/s]';...
              'Elastic Modulus [KPa]';...
              'Diameter [mm]';...
              'Curvature [1/m]';...
              'EOAdist [mm]';...
              'Pressure Recovery Distance [mm]';...
              'Sections [-]'};
    
        MM = [AL';FP';SP';PW';EL';DI';CU';EOAdist_vect';PR';sections'];
        TT = table(VV,MM);
        TT.Properties.VariableNames = {' ','Point'};
    
    else 
    
        VV = {'Aortic Length [mm]';...
              'Flow at peak frame [l/min]';...
              'SAW at peak frame [mmHg]';...
              'Pulse Wave Velocity [m/s]';...
              'Elastic Modulus [KPa]';...
              'Diameter [mm]';...
              'Curvature [1/m]';...
              'Sections [-]'};
        
        MM = [AL';FP';SP';PW';EL';DI';CU';sections'];
        TT = table(VV,MM);
        TT.Properties.VariableNames = {' ','Point'};
    
    end

end

filename = [directory,'/SAW_CURVES.xlsx'];
writetable(TT,fullfile(filename),'Sheet','SAW','Range','A1')


f = msgbox("The data has been successfully saved","Success");
pause(0.5)
close(f)

handles.output = hObject;
guidata(hObject, handles);




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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(handles.popupmenu1,'Value')
    case 7
        
            opts = [];
            opts.dist         = [ handles.PM.AllPlanes.dist ];
            opts.bManuPlateau = 1; %manual definition of the threshold
            opts.bSmooth      = 1; %binary to decide on smoothing (1) or not (0) the SAW trace
            opts.idx          = handles.idx; 
            D = handles.PM.SAW( :, handles.PM.peakFrame );
            D = smooth( D, length( D )/6,'rlowess');
            opts.D            = D;
            handles.PM = handles.PM.ComputePrec(opts);
            
            
            % update the figure when the windows is closed
            if handles.id_case_5 == 0
                handles.idx = 1 : floor(handles.PM.nPlanes/4) : handles.PM.nPlanes;
                handles.idx(end) = handles.PM.nPlanes;
            end
            
            dist = [ handles.PM.AllPlanes.dist ];
            D = handles.PM.SAW( :, handles.PM.peakFrame );
            yUnit = '[mmHg]';
            titleStr = 'SAW\nat peak frame';
            D = smooth( D, length( D )/6,'rlowess');
            
            axes(handles.axes1)
            plot(0,0)
            plot( dist( handles.idx(1):handles.idx(end) )-dist( handles.idx(1) ), D( handles.idx(1):handles.idx(end) ),  'r', 'LineWidth', 3 );
            hold on
            ax = gca;
            axis( ax, 'tight' )
            for j = 1 : 2 : size( handles.idx, 2 )-1
                patch( 'Faces',[1 2 3 4], 'Vertices', ...
                    [ [ dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(1)  ]; ...
                    [   dist(handles.idx(j+1)-handles.idx(1)+1) ax.YLim(2)  ]; ...
                    [   dist(handles.idx(j)  -handles.idx(1)+1) ax.YLim(2)  ] ], 'FaceAlpha', 0.1, 'EdgeColor', 'none' );
            end
            
            xlabel( 'Aortic Length [mm]' )
            ylabel( yUnit )
            ax.Title.String = sprintf( titleStr );
            
            ax.FontName = 'Calibri';
            ax.FontSize = 14;
            ax.FontWeight = 'normal';
            axis square
            
            ax.XTick = ax.XLim(1) : 25 : ax.XLim(2);
            hold off
            axtoolbar(handles.axes1,'Visible','off');
            
            [~,EOAIdx] = max(D);   
            handles.EOAdist = dist(EOAIdx);
            xline(handles.EOAdist, 'g-.','LineWidth', 2 );
            xline(handles.EOAdist + handles.PM.PrecDist, 'b-.','LineWidth', 2 );
            
            text(handles.EOAdist + handles.PM.PrecDist, ...
                ((max(D( handles.idx(1):handles.idx(end) )) + min(D( handles.idx(1):handles.idx(end) )))/2),...
                ['\leftarrow PR Distance = ', num2str(handles.PM.PrecDist,4), ' [mm]'] ,'FontSize',14)
            
            handles.pressure_recovery_id = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5

        set(handles.pushbutton2,'visible','on','string','ARCH PLANES RE-LOCATION');
    
        % we create the mip only with the slice that contain the aorta
        [~,I1] = min(abs(squeeze(handles.Z_Pos(1,1,:)) - min(handles.nodes(:,3))));
        [~,I2] = min(abs(squeeze(handles.Z_Pos(1,1,:)) - max(handles.nodes(:,3))));
        MIP = max(handles.IPCMRA(:,:,I1:I2),[],3);
        MIP = (MIP/max(MIP(:))*255);
        
        axes(handles.axes1)
        plot(0,0)
        surface(handles.Y_Pos(:,:,1),handles.X_Pos(:,:,1),handles.Z_Pos(:,:,end),MIP,'FaceColor','texturemap','EdgeColor','none','CData',uint8(repmat(MIP,1,1,3)));
        hold on
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace*255,'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        hold off
        axis off
        view(0,-90)
        daspect([1 1 1])
        
        position = zeros(2,3);
        
        while(1)
            
            uiwait(msgbox({'Please select one point Before Braquiocephalic Artery, and press enter','Close this message first...'}));
            dcmObject = datacursormode;
            pause
            datacursormode off
            cursor = getCursorInfo(dcmObject);
            pos = [cursor.Position(1),cursor.Position(2),cursor.Position(3)];
            
            % We select the node close to the selected position
            d = sqrt((handles.nodes(:,1)-pos(1)).^2 + (handles.nodes(:,2)-pos(2)).^2 + (handles.nodes(:,3)-pos(3)).^2);
            [~,I] = min(d);
            value = handles.Laplace(I);
            
            [cutpos,~,facedata,~] = qmeshcut(handles.elem,handles.nodes,handles.Laplace,value);
            
            % move the centerline location
            movement = handles.Centerline(1,:)-handles.clPointsi(1,:);
            
            axes(handles.axes1)
            plot(0,0)
            surface(handles.Y_Pos(:,:,1),handles.X_Pos(:,:,1),handles.Z_Pos(:,:,end),MIP,'FaceColor','texturemap','EdgeColor','none','CData',uint8(repmat(MIP,1,1,3)));
            hold on
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace*255,'CDataMapping','Scaled', 'facealpha',0.2)
            colormap(handles.axes1,'cool');
            patch('faces',facedata,'vertices',cutpos,'EdgeColor','y','FaceColor','y')
            plot3( handles.clPointsi(:,1) + movement(1), handles.clPointsi(:,2) + movement(2), handles.clPointsi(:,3) + movement(3), 'r-', 'linewidth',3)
            hold off
            axis off
            view(0,-90)
            daspect([1 1 1])
            
            
            answer = questdlg('Are you agree with the section selected?', ...
                'Question', ...
                'Yes','No','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    position(1,:) = mean(cutpos);
                    break
                case 'No'
                    disp('Select another location ...')
                    
            end
            
        end
        
        
        axes(handles.axes1)
        plot(0,0)
        surface(handles.Y_Pos(:,:,1),handles.X_Pos(:,:,1),handles.Z_Pos(:,:,end),MIP,'FaceColor','texturemap','EdgeColor','none','CData',uint8(repmat(MIP,1,1,3)));
        hold on
        patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace*255,'CDataMapping','Scaled')
        colormap(handles.axes1,'cool');
        hold off
        axis off
        view(0,-90)
        daspect([1 1 1])
        
        while(1)
            
            uiwait(msgbox({'Please select one point After Left Subclavian Artery, and press enter','Close this message first...'}));
            dcmObject = datacursormode;
            pause
            datacursormode off
            cursor = getCursorInfo(dcmObject);
            pos = [cursor.Position(1),cursor.Position(2),cursor.Position(3)];
            
            % We select the node close to the selected position
            d = sqrt((handles.nodes(:,1)-pos(1)).^2 + (handles.nodes(:,2)-pos(2)).^2 + (handles.nodes(:,3)-pos(3)).^2);
            [~,I] = min(d);
            value = handles.Laplace(I);
            
            [cutpos,~,facedata,~] = qmeshcut(handles.elem,handles.nodes,handles.Laplace,value);
            
            axes(handles.axes1)
            plot(0,0)
            surface(handles.Y_Pos(:,:,1),handles.X_Pos(:,:,1),handles.Z_Pos(:,:,end),MIP,'FaceColor','texturemap','EdgeColor','none','CData',uint8(repmat(MIP,1,1,3)));
            hold on
            patch('faces',handles.faces,'vertices',handles.nodes,'EdgeColor','k','FaceColor','interp','FaceVertexCData',handles.Laplace*255,'CDataMapping','Scaled', 'facealpha',0.2)
            colormap(handles.axes1,'cool');
            patch('faces',facedata,'vertices',cutpos,'EdgeColor','y','FaceColor','y')
            plot3( handles.clPointsi(:,1) + movement(1), handles.clPointsi(:,2) + movement(2), handles.clPointsi(:,3) + movement(3), 'r-', 'linewidth',3)
            hold off
            axis off
            view(0,-90)
            daspect([1 1 1])
            
            
            answer = questdlg( 'Are you agree with the section selected?', ...
                'Question', ...
                'Yes','No','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    position(2,:) = mean(cutpos);
                    break
                case 'No'
                    disp('Select another location ...')
                    
            end
            
        end
        
        disp(position)
        
        centerline_moved = [handles.clPointsi(:,1) + movement(1), handles.clPointsi(:,2) + movement(2), handles.clPointsi(:,3) + movement(3)];
        d1 =  sqrt(sum((centerline_moved - repmat(position(1,:),size(centerline_moved,1),1)).^2,2));
        d2 =  sqrt(sum((centerline_moved - repmat(position(2,:),size(centerline_moved,1),1)).^2,2));
        
        [~,handles.idx1] = min(d1); % plane before the arch
        [~,handles.idx2] = min(d2); % plane after the arch
        
        % we select point than in F-H direction has 6 mm of distance with
        % respect to the initial point in the centerline
        id_pos = 1:size(handles.clPointsi,1);
        id_pos_sel = id_pos((handles.clPointsi(:,2) > handles.clPointsi(1,2)-6) & (handles.clPointsi(:,2) < handles.clPointsi(1,2)+6));
        
        % Then we calculate the shortest distance between the first point and the selected points,
        % with the restriction that the distance need to be bigger than the
        % mean distance between all points selected
        
        mean_point_selected = mean(handles.clPointsi(id_pos_sel,:));
        distanca_mean_point = sqrt(sum((handles.clPointsi(id_pos_sel,:) - repmat(mean_point_selected,length(id_pos_sel),1)).^2,2));
        mean_distance = mean(distanca_mean_point);
        
        d =  sqrt(sum((handles.clPointsi(id_pos_sel,:) - repmat(handles.clPointsi(1,:),length(id_pos_sel),1)).^2,2));
        new_id_pos_sel = id_pos_sel(d>mean_distance);
        d =  sqrt(sum((handles.clPointsi(new_id_pos_sel,:) - repmat(handles.clPointsi(1,:),length(new_id_pos_sel),1)).^2,2));
        [~,selected_plane] = min(d);
        handles.idx3 = new_id_pos_sel(selected_plane);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.idx = [1,handles.idx1, handles.idx2, handles.idx3, handles.PM.nPlanes];
        
        axes(handles.axes1)
        plot(0,0)
        handles.PM.plotDomain();
        hold on
        plot3( handles.clPointsi(:,1), handles.clPointsi(:,2), handles.clPointsi(:,3), 'r-' )
        plot3( handles.clPointsi(handles.idx,1), handles.clPointsi(handles.idx,2), handles.clPointsi(handles.idx,3), '*g' )
        optPlotPlane.bPlotVelocityVectors = 0;
        optPlotPlane.bPlotKeyPoint = 0;
        optPlotPlane.bPlotKeySlope = 0;
        P2plot = handles.PM.AllPlanes( handles.idx );
        for iL = 1 : size( handles.idx, 2 )
            P = P2plot(iL);
            handles.PM.plotPlane( P, handles.PM.peakFrame, optPlotPlane )
        end
        hold off
        grid on
        axis on
        ylabel('F-H')
        xlabel('A-P')
        zlabel('R-L')
        view(142,39)
        daspect([1 1 1])
        axtoolbar(handles.axes1,{'rotate', 'restoreview'},'Visible','on');
        
        handles.arch_planes_selection_id = 1;
end

handles.output = hObject;
guidata(hObject, handles);
