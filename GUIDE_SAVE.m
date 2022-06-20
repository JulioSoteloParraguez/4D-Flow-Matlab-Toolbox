function varargout = GUIDE_SAVE(varargin)
% GUIDE_SAVE MATLAB code for GUIDE_SAVE.fig
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
                   'gui_OpeningFcn', @GUIDE_SAVE_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_SAVE_OutputFcn, ...
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
function GUIDE_SAVE_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

        % initialization of the handles to generate an automatic selecction
        % of the variables
        handles.con_all_data = 0;
        handles.con_xls_data = 0;
        handles.con_mat_data = 0;
        handles.con_vti_data = 0;
        handles.con_vtu_data = 0;


        col = zeros(40,4);
        
        %csv files
        col(6,1) = varargin{1}.save_id_vel_csv;
        col(7,1) = varargin{1}.save_id_wss_csv;
        col(8,1) = varargin{1}.save_id_osi_csv;
        col(9,1) = varargin{1}.save_id_vor_csv;
        col(10,1) = varargin{1}.save_id_hd_csv;
        col(11,1) = varargin{1}.save_id_rhd_csv;
        col(12,1) = varargin{1}.save_id_vd_csv;
        col(13,1) = varargin{1}.save_id_el_csv;
        col(14,1) = varargin{1}.save_id_ke_csv;
        col(17,1) = varargin{1}.save_id_rad_csv;
        col(18,1) = varargin{1}.save_id_dia_csv;
        col(21,1) = varargin{1}.save_id_wssa_csv;
        col(22,1) = varargin{1}.save_id_wssc_csv;
        col(23,1) = varargin{1}.save_id_aan_csv;
        col(24,1) = varargin{1}.save_id_fve_csv;
        col(25,1) = varargin{1}.save_id_bve_csv;
        col(26,1) = varargin{1}.save_id_ref_csv;
        col(28,1) = varargin{1}.save_id_ecc_csv;
        col(29,1) = varargin{1}.save_id_cur_csv; % Julio Sotelo 28-05-2019
        col(30,1) = varargin{1}.save_id_ell_csv; % Julio Sotelo 28-05-2019
        col(31,1) = varargin{1}.save_id_len_csv; % Julio Sotelo 28-05-2019
%         col(32,1) = varargin{1}.save_id_cir_csv; % Julio Sotelo 28-05-2019
        col(32,1) = varargin{1}.save_id_fov_csv; % Julio Sotelo 28-05-2019
        col(33,1) = varargin{1}.save_id_fla_csv; % Julio Sotelo 28-05-2019
        col(34,1) = varargin{1}.save_id_are_csv; % Julio Sotelo 28-05-2019
        col(35,1) = varargin{1}.save_id_aci_csv; % Julio Sotelo 28-05-2019
        col(36,1) = varargin{1}.save_id_tim_csv; % Julio Sotelo 28-05-2019 time
        col(37,1) = varargin{1}.save_id_flo_csv; % Julio Sotelo 28-05-2019 flow
        col(38,1) = varargin{1}.save_id_nfl_csv; % Julio Sotelo 28-05-2019 net_flow
        col(39,1) = varargin{1}.save_id_mav_csv; % Julio Sotelo 28-05-2019 max_velocity
        col(40,1) = varargin{1}.save_id_miv_csv; % Julio Sotelo 28-05-2019 min_velocity
        

        % matlab files
        col(1,2) = varargin{1}.save_id_SEG_mat;
        col(2,2) = varargin{1}.save_id_IPCMRA_mat;
        col(3,2) = varargin{1}.save_id_MR_FFE_mat;
        col(4,2) = varargin{1}.save_id_MR_PCA_mat;
        col(5,2) = varargin{1}.save_id_mesh_mat;
        col(6,2) = varargin{1}.save_id_vel_mat;
        col(7,2) = varargin{1}.save_id_wss_mat;
        col(8,2) = varargin{1}.save_id_osi_mat;
        col(9,2) = varargin{1}.save_id_vor_mat;
        col(10,2) = varargin{1}.save_id_hd_mat;
        col(11,2) = varargin{1}.save_id_rhd_mat;
        col(12,2) = varargin{1}.save_id_vd_mat;
        col(13,2) = varargin{1}.save_id_el_mat;
        col(14,2) = varargin{1}.save_id_ke_mat;
        col(15,2) = varargin{1}.save_id_lap_mat;
        col(16,2) = varargin{1}.save_id_cen_mat;
        col(17,2) = varargin{1}.save_id_rad_mat;
        col(18,2) = varargin{1}.save_id_dia_mat;
        col(19,2) = varargin{1}.save_id_auv_mat;
        col(20,2) = varargin{1}.save_id_cuv_mat;
        col(21,2) = varargin{1}.save_id_wssa_mat;
        col(22,2) = varargin{1}.save_id_wssc_mat;
        col(23,2) = varargin{1}.save_id_aan_mat;
        col(24,2) = varargin{1}.save_id_fve_mat;
        col(25,2) = varargin{1}.save_id_bve_mat;
        col(26,2) = varargin{1}.save_id_ref_mat;
        col(27,2) = varargin{1}.save_id_cebf_mat;
        col(28,2) = varargin{1}.save_id_ecc_mat;
        col(29,2) = varargin{1}.save_id_cur_mat; % Julio Sotelo 28-05-2019
        col(30,2) = varargin{1}.save_id_ell_mat; % Julio Sotelo 28-05-2019
        col(31,2) = varargin{1}.save_id_len_mat; % Julio Sotelo 28-05-2019
%         col(32,2) = varargin{1}.save_id_cir_mat; % Julio Sotelo 28-05-2019
        col(32,2) = varargin{1}.save_id_fov_mat; % Julio Sotelo 28-05-2019
        col(33,2) = varargin{1}.save_id_fla_mat; % Julio Sotelo 28-05-2019
        col(34,2) = varargin{1}.save_id_are_mat; % Julio Sotelo 28-05-2019
        col(35,2) = varargin{1}.save_id_aci_mat; % Julio Sotelo 28-05-2019
        col(36,2) = varargin{1}.save_id_tim_mat; % Julio Sotelo 28-05-2019 time
        col(37,2) = varargin{1}.save_id_flo_mat; % Julio Sotelo 28-05-2019 flow
        col(38,2) = varargin{1}.save_id_nfl_mat; % Julio Sotelo 28-05-2019 net_flow
        col(39,2) = varargin{1}.save_id_mav_mat; % Julio Sotelo 28-05-2019 max_velocity
        col(40,2) = varargin{1}.save_id_miv_mat; % Julio Sotelo 28-05-2019 min_velocity
        
        % vti files
        col(1,3) = varargin{1}.save_id_SEG_vti;
        col(2,3) = varargin{1}.save_id_IPCMRA_vti;
        col(3,3) = varargin{1}.save_id_MR_FFE_vti;
        col(4,3) = varargin{1}.save_id_MR_PCA_vti;
        
        % vtu files
        col(5,4) = varargin{1}.save_id_mesh_vtu;
        col(6,4) = varargin{1}.save_id_vel_vtu;
        col(7,4) = varargin{1}.save_id_wss_vtu;
        col(8,4) = varargin{1}.save_id_osi_vtu;
        col(9,4) = varargin{1}.save_id_vor_vtu;
        col(10,4) = varargin{1}.save_id_hd_vtu;
        col(11,4) = varargin{1}.save_id_rhd_vtu;
        col(12,4) = varargin{1}.save_id_vd_vtu;
        col(13,4) = varargin{1}.save_id_el_vtu;
        col(14,4) = varargin{1}.save_id_ke_vtu;
        col(15,4) = varargin{1}.save_id_lap_vtu;
        col(16,4) = varargin{1}.save_id_cen_vtu;
        col(17,4) = varargin{1}.save_id_rad_vtu;
        col(18,4) = varargin{1}.save_id_dia_vtu;
        col(19,4) = varargin{1}.save_id_auv_vtu;
        col(20,4) = varargin{1}.save_id_cuv_vtu;
        col(21,4) = varargin{1}.save_id_wssa_vtu;
        col(22,4) = varargin{1}.save_id_wssc_vtu;
        col(23,4) = varargin{1}.save_id_aan_vtu;
        col(24,4) = varargin{1}.save_id_fve_vtu;
        col(25,4) = varargin{1}.save_id_bve_vtu;
        col(26,4) = varargin{1}.save_id_ref_vtu;
        col(27,4) = varargin{1}.save_id_cebf_vtu;
        col(28,4) = varargin{1}.save_id_ecc_vtu;
        col(29,4) = varargin{1}.save_id_cur_vtu; % Julio Sotelo 28-05-2019
        col(30,4) = varargin{1}.save_id_ell_vtu; % Julio Sotelo 28-05-2019
        col(31,4) = varargin{1}.save_id_len_vtu; % Julio Sotelo 28-05-2019
%         col(32,4) = varargin{1}.save_id_cir_vtu; % Julio Sotelo 28-05-2019
        col(32,4) = varargin{1}.save_id_fov_vtu; % Julio Sotelo 28-05-2019
        col(33,4) = varargin{1}.save_id_fla_vtu; % Julio Sotelo 28-05-2019
        col(34,4) = varargin{1}.save_id_are_vtu; % Julio Sotelo 28-05-2019
        col(35,4) = varargin{1}.save_id_aci_vtu; % Julio Sotelo 28-05-2019


        if sum(col(:,1))==0
            set(handles.pushbutton3,'visible','off')
        end
        
        set(handles.uitable1,'ColumnName',{'Variable','.xls','.mat','.vti','.vtu'},'ColumnEditable',true)
        data = cell(40,5);
        data{1,1} = 'Segmentation ROI';
        data{2,1} = 'Angiography ROI';
        data{3,1} = 'Magnitude ROI';
        data{4,1} = 'Velocity ROI (m/s)';
        data{5,1} = 'FE - Mesh';
        data{6,1} = 'FE - Velocity (m/s)';
        data{7,1} = 'FE - WSS (N/m^{2}';
        data{8,1} = 'FE - OSI (-)';
        data{9,1} = 'FE - Vorticity (1/s)';
        data{10,1} = 'FE - Helicity Density (m/s^{2})';
        data{11,1} = 'FE - Relative Helicity Density (-)';
        data{12,1} = 'FE - Viscous Dissipation (1/s^{2})';
        data{13,1} = 'FE - Energy Loss (\muW)';
        data{14,1} = 'FE - Kinetic Energy (\muJ)';
        data{15,1} = 'FE - Laplace (-)';
        data{16,1} = 'FE - Centerline (-)';
        data{17,1} = 'FE - Radius (cm)';
        data{18,1} = 'FE - Diameter (cm)';
        data{19,1} = 'FE - Axial Unit Vector (-)';
        data{20,1} = 'FE - Circum. Unit Vector (-)';
        data{21,1} = 'FE - WSS-A (N/m^{2})';
        data{22,1} = 'FE - WSS-C (N/m^{2})';
        data{23,1} = 'FE - Axial Angle (^{o})';
        data{24,1} = 'FE - Forward Velocity (m/s)';
        data{25,1} = 'FE - Backward Velocity (m/s)';
        data{26,1} = 'FE - Regurgitant Flow (%)';
        data{27,1} = 'FE - Centerline Flow (-)';
        data{28,1} = 'FE - Eccentricity (%)';
        data{29,1} = 'FE - Curvature (1/m)'; % Julio Sotelo 28-05-2019
        data{30,1} = 'FE - Ellipticity (-)'; % Julio Sotelo 28-05-2019
        data{31,1} = 'FE - Length Vessel (mm)'; % Julio Sotelo 28-05-2019
%         data{32,1} = 'FE - Circulation (mm^{2}/m)'; % Julio Sotelo 28-05-2019
        data{32,1} = 'FE - Forward vortex (1/s)'; % Julio Sotelo 28-05-2019
        data{33,1} = 'FE - Flattening (-)'; % Julio Sotelo 28-05-2019
        data{34,1} = 'FE - Area (cm^{2})'; % Julio Sotelo 28-05-2019
        data{35,1} = 'FE - Axial Circulation (cm^{2}/m)'; % Julio Sotelo 28-05-2019
        data{36,1} = 'FE - Time (s)'; % Julio Sotelo 28-05-2019
        data{37,1} = 'FE - Flow (cm^{3}/s)'; % Julio Sotelo 28-05-2019
        data{38,1} = 'FE - Net Flow (cm^{3})'; % Julio Sotelo 28-05-2019
        data{39,1} = 'FE - Maximum Velocity (cm/m)'; % Julio Sotelo 28-05-2019
        data{40,1} = 'FE - Minimum Velocity (cm/m)'; % Julio Sotelo 28-05-2019

        if col(1,1)==1, data{1,2} = false; else, data{1,2} = '         NA'; end
        if col(2,1)==1, data{2,2} = false; else, data{2,2} = '         NA'; end
        if col(3,1)==1, data{3,2} = false; else, data{3,2} = '         NA'; end
        if col(4,1)==1, data{4,2} = false; else, data{4,2} = '         NA'; end
        if col(5,1)==1, data{5,2} = false; else, data{5,2} = '         NA'; end
        if col(6,1)==1, data{6,2} = false; else, data{6,2} = '         NA'; end
        if col(7,1)==1, data{7,2} = false; else, data{7,2} = '         NA'; end
        if col(8,1)==1, data{8,2} = false; else, data{8,2} = '         NA'; end
        if col(9,1)==1, data{9,2} = false; else, data{9,2} = '         NA'; end
        if col(10,1)==1, data{10,2} = false; else, data{10,2} = '         NA'; end
        if col(11,1)==1, data{11,2} = false; else, data{11,2} = '         NA'; end
        if col(12,1)==1, data{12,2} = false; else, data{12,2} = '         NA'; end
        if col(13,1)==1, data{13,2} = false; else, data{13,2} = '         NA'; end
        if col(14,1)==1, data{14,2} = false; else, data{14,2} = '         NA'; end
        if col(15,1)==1, data{15,2} = false; else, data{15,2} = '         NA'; end
        if col(16,1)==1, data{16,2} = false; else, data{16,2} = '         NA'; end
        if col(17,1)==1, data{17,2} = false; else, data{17,2} = '         NA'; end
        if col(18,1)==1, data{18,2} = false; else, data{18,2} = '         NA'; end
        if col(19,1)==1, data{19,2} = false; else, data{19,2} = '         NA'; end
        if col(20,1)==1, data{20,2} = false; else, data{20,2} = '         NA'; end
        if col(21,1)==1, data{21,2} = false; else, data{21,2} = '         NA'; end
        if col(22,1)==1, data{22,2} = false; else, data{22,2} = '         NA'; end
        if col(23,1)==1, data{23,2} = false; else, data{23,2} = '         NA'; end
        if col(24,1)==1, data{24,2} = false; else, data{24,2} = '         NA'; end
        if col(25,1)==1, data{25,2} = false; else, data{25,2} = '         NA'; end
        if col(26,1)==1, data{26,2} = false; else, data{26,2} = '         NA'; end
        if col(27,1)==1, data{27,2} = false; else, data{27,2} = '         NA'; end
        if col(28,1)==1, data{28,2} = false; else, data{28,2} = '         NA'; end
        if col(29,1)==1, data{29,2} = false; else, data{29,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(30,1)==1, data{30,2} = false; else, data{30,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(31,1)==1, data{31,2} = false; else, data{31,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(32,1)==1, data{32,2} = false; else, data{32,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(33,1)==1, data{33,2} = false; else, data{33,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(34,1)==1, data{34,2} = false; else, data{34,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(35,1)==1, data{35,2} = false; else, data{35,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(36,1)==1, data{36,2} = false; else, data{36,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(37,1)==1, data{37,2} = false; else, data{37,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(38,1)==1, data{38,2} = false; else, data{38,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(39,1)==1, data{39,2} = false; else, data{39,2} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(40,1)==1, data{40,2} = false; else, data{40,2} = '         NA'; end % Julio Sotelo 28-05-2019
%         if col(41,1)==1, data{41,2} = false; else, data{41,2} = '         NA'; end % Julio Sotelo 28-05-2019
        
        
        if col(1,2)==1, data{1,3} = false; else, data{1,3} = '         NA'; end
        if col(2,2)==1, data{2,3} = false; else, data{2,3} = '         NA'; end
        if col(3,2)==1, data{3,3} = false; else, data{3,3} = '         NA'; end
        if col(4,2)==1, data{4,3} = false; else, data{4,3} = '         NA'; end
        if col(5,2)==1, data{5,3} = false; else, data{5,3} = '         NA'; end
        if col(6,2)==1, data{6,3} = false; else, data{6,3} = '         NA'; end
        if col(7,2)==1, data{7,3} = false; else, data{7,3} = '         NA'; end
        if col(8,2)==1, data{8,3} = false; else, data{8,3} = '         NA'; end
        if col(9,2)==1, data{9,3} = false; else, data{9,3} = '         NA'; end
        if col(10,2)==1, data{10,3} = false; else, data{10,3} = '         NA'; end
        if col(11,2)==1, data{11,3} = false; else, data{11,3} = '         NA'; end
        if col(12,2)==1, data{12,3} = false; else, data{12,3} = '         NA'; end
        if col(13,2)==1, data{13,3} = false; else, data{13,3} = '         NA'; end
        if col(14,2)==1, data{14,3} = false; else, data{14,3} = '         NA'; end
        if col(15,2)==1, data{15,3} = false; else, data{15,3} = '         NA'; end
        if col(16,2)==1, data{16,3} = false; else, data{16,3} = '         NA'; end
        if col(17,2)==1, data{17,3} = false; else, data{17,3} = '         NA'; end
        if col(18,2)==1, data{18,3} = false; else, data{18,3} = '         NA'; end
        if col(19,2)==1, data{19,3} = false; else, data{19,3} = '         NA'; end
        if col(20,2)==1, data{20,3} = false; else, data{20,3} = '         NA'; end
        if col(21,2)==1, data{21,3} = false; else, data{21,3} = '         NA'; end
        if col(22,2)==1, data{22,3} = false; else, data{22,3} = '         NA'; end
        if col(23,2)==1, data{23,3} = false; else, data{23,3} = '         NA'; end
        if col(24,2)==1, data{24,3} = false; else, data{24,3} = '         NA'; end
        if col(25,2)==1, data{25,3} = false; else, data{25,3} = '         NA'; end
        if col(26,2)==1, data{26,3} = false; else, data{26,3} = '         NA'; end
        if col(27,2)==1, data{27,3} = false; else, data{27,3} = '         NA'; end
        if col(28,2)==1, data{28,3} = false; else, data{28,3} = '         NA'; end
        if col(29,2)==1, data{29,3} = false; else, data{29,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(30,2)==1, data{30,3} = false; else, data{30,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(31,2)==1, data{31,3} = false; else, data{31,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(32,2)==1, data{32,3} = false; else, data{32,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(33,2)==1, data{33,3} = false; else, data{33,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(34,2)==1, data{34,3} = false; else, data{34,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(35,2)==1, data{35,3} = false; else, data{35,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(36,2)==1, data{36,3} = false; else, data{36,3} = '         NA'; end % Julio Sotelo 28-05-2019
        
        if col(37,2)==1, data{37,3} = false; else, data{37,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(38,2)==1, data{38,3} = false; else, data{38,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(39,2)==1, data{39,3} = false; else, data{39,3} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(40,2)==1, data{40,3} = false; else, data{40,3} = '         NA'; end % Julio Sotelo 28-05-2019
%         if col(41,2)==1, data{41,3} = false; else, data{41,3} = '         NA'; end % Julio Sotelo 28-05-2019
        

        if col(1,3)==1, data{1,4} = false; else, data{1,4} = '         NA'; end
        if col(2,3)==1, data{2,4} = false; else, data{2,4} = '         NA'; end
        if col(3,3)==1, data{3,4} = false; else, data{3,4} = '         NA'; end
        if col(4,3)==1, data{4,4} = false; else, data{4,4} = '         NA'; end
        if col(5,3)==1, data{5,4} = false; else, data{5,4} = '         NA'; end
        if col(6,3)==1, data{6,4} = false; else, data{6,4} = '         NA'; end
        if col(7,3)==1, data{7,4} = false; else, data{7,4} = '         NA'; end
        if col(8,3)==1, data{8,4} = false; else, data{8,4} = '         NA'; end
        if col(9,3)==1, data{9,4} = false; else, data{9,4} = '         NA'; end
        if col(10,3)==1, data{10,4} = false; else, data{10,4} = '         NA'; end
        if col(11,3)==1, data{11,4} = false; else, data{11,4} = '         NA'; end
        if col(12,3)==1, data{12,4} = false; else, data{12,4} = '         NA'; end
        if col(13,3)==1, data{13,4} = false; else, data{13,4} = '         NA'; end
        if col(14,3)==1, data{14,4} = false; else, data{14,4} = '         NA'; end
        if col(15,3)==1, data{15,4} = false; else, data{15,4} = '         NA'; end
        if col(16,3)==1, data{16,4} = false; else, data{16,4} = '         NA'; end
        if col(17,3)==1, data{17,4} = false; else, data{17,4} = '         NA'; end
        if col(18,3)==1, data{18,4} = false; else, data{18,4} = '         NA'; end
        if col(19,3)==1, data{19,4} = false; else, data{19,4} = '         NA'; end
        if col(20,3)==1, data{20,4} = false; else, data{20,4} = '         NA'; end
        if col(21,3)==1, data{21,4} = false; else, data{21,4} = '         NA'; end
        if col(22,3)==1, data{22,4} = false; else, data{22,4} = '         NA'; end
        if col(23,3)==1, data{23,4} = false; else, data{23,4} = '         NA'; end
        if col(24,3)==1, data{24,4} = false; else, data{24,4} = '         NA'; end
        if col(25,3)==1, data{25,4} = false; else, data{25,4} = '         NA'; end
        if col(26,3)==1, data{26,4} = false; else, data{26,4} = '         NA'; end
        if col(27,3)==1, data{27,4} = false; else, data{27,4} = '         NA'; end
        if col(28,3)==1, data{28,4} = false; else, data{28,4} = '         NA'; end
        if col(29,3)==1, data{29,4} = false; else, data{29,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(30,3)==1, data{30,4} = false; else, data{30,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(31,3)==1, data{31,4} = false; else, data{31,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(32,3)==1, data{32,4} = false; else, data{32,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(33,3)==1, data{33,4} = false; else, data{33,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(34,3)==1, data{34,4} = false; else, data{34,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(35,3)==1, data{35,4} = false; else, data{35,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(36,3)==1, data{36,4} = false; else, data{36,4} = '         NA'; end % Julio Sotelo 28-05-2019
        
        if col(37,3)==1, data{37,4} = false; else, data{37,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(38,3)==1, data{38,4} = false; else, data{38,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(39,3)==1, data{39,4} = false; else, data{39,4} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(40,3)==1, data{40,4} = false; else, data{40,4} = '         NA'; end % Julio Sotelo 28-05-2019
%         if col(41,3)==1, data{41,4} = false; else, data{41,4} = '         NA'; end % Julio Sotelo 28-05-2019
        
        
        if col(1,4)==1, data{1,5} = false; else, data{1,5} = '         NA'; end
        if col(2,4)==1, data{2,5} = false; else, data{2,5} = '         NA'; end
        if col(3,4)==1, data{3,5} = false; else, data{3,5} = '         NA'; end
        if col(4,4)==1, data{4,5} = false; else, data{4,5} = '         NA'; end
        if col(5,4)==1, data{5,5} = false; else, data{5,5} = '         NA'; end
        if col(6,4)==1, data{6,5} = false; else, data{6,5} = '         NA'; end
        if col(7,4)==1, data{7,5} = false; else, data{7,5} = '         NA'; end
        if col(8,4)==1, data{8,5} = false; else, data{8,5} = '         NA'; end
        if col(9,4)==1, data{9,5} = false; else, data{9,5} = '         NA'; end
        if col(10,4)==1, data{10,5} = false; else, data{10,5} = '         NA'; end
        if col(11,4)==1, data{11,5} = false; else, data{11,5} = '         NA'; end
        if col(12,4)==1, data{12,5} = false; else, data{12,5} = '         NA'; end
        if col(13,4)==1, data{13,5} = false; else, data{13,5} = '         NA'; end
        if col(14,4)==1, data{14,5} = false; else, data{14,5} = '         NA'; end
        if col(15,4)==1, data{15,5} = false; else, data{15,5} = '         NA'; end
        if col(16,4)==1, data{16,5} = false; else, data{16,5} = '         NA'; end
        if col(17,4)==1, data{17,5} = false; else, data{17,5} = '         NA'; end
        if col(18,4)==1, data{18,5} = false; else, data{18,5} = '         NA'; end
        if col(19,4)==1, data{19,5} = false; else, data{19,5} = '         NA'; end
        if col(20,4)==1, data{20,5} = false; else, data{20,5} = '         NA'; end
        if col(21,4)==1, data{21,5} = false; else, data{21,5} = '         NA'; end
        if col(22,4)==1, data{22,5} = false; else, data{22,5} = '         NA'; end
        if col(23,4)==1, data{23,5} = false; else, data{23,5} = '         NA'; end
        if col(24,4)==1, data{24,5} = false; else, data{24,5} = '         NA'; end
        if col(25,4)==1, data{25,5} = false; else, data{25,5} = '         NA'; end
        if col(26,4)==1, data{26,5} = false; else, data{26,5} = '         NA'; end
        if col(27,4)==1, data{27,5} = false; else, data{27,5} = '         NA'; end
        if col(28,4)==1, data{28,5} = false; else, data{28,5} = '         NA'; end
        if col(29,4)==1, data{29,5} = false; else, data{29,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(30,4)==1, data{30,5} = false; else, data{30,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(31,4)==1, data{31,5} = false; else, data{31,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(32,4)==1, data{32,5} = false; else, data{32,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(33,4)==1, data{33,5} = false; else, data{33,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(34,4)==1, data{34,5} = false; else, data{34,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(35,4)==1, data{35,5} = false; else, data{35,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(36,4)==1, data{36,5} = false; else, data{36,5} = '         NA'; end % Julio Sotelo 28-05-2019
        
        if col(37,4)==1, data{37,5} = false; else, data{37,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(38,4)==1, data{38,5} = false; else, data{38,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(39,4)==1, data{39,5} = false; else, data{39,5} = '         NA'; end % Julio Sotelo 28-05-2019
        if col(40,4)==1, data{40,5} = false; else, data{40,5} = '         NA'; end % Julio Sotelo 28-05-2019
%         if col(41,4)==1, data{41,5} = false; else, data{41,5} = '         NA'; end % Julio Sotelo 28-05-2019
        
        
        
        set(handles.uitable1, 'Units', 'pixels');
        AA = get(handles.uitable1,'Position');
        FigPos = get(handles.uitable1, 'Position');
        set(handles.uitable1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
        set(handles.uitable1, 'Units', 'normalized');
        
        set(handles.uitable1,'data',data,'visible','on','ColumnWidth',{round(AA(3)/2) round(AA(3)/2)*0.2 round(AA(3)/2)*0.2 round(AA(3)/2)*0.2 round(AA(3)/2)*0.2})
        handles.save_selection = zeros(size(col,1),4);
        
guidata(hObject, handles);
uiwait(handles.figure1);
function varargout = GUIDE_SAVE_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;
    save_selection = handles.save_selection;
    setappdata(0,'save_selection',save_selection);
    delete(handles.figure1);
function pushbutton1_Callback(hObject, eventdata, handles)
    tableData = get(handles.uitable1, 'data');
    MAT = zeros(size(tableData,1),4);
    for n=1:size(tableData,1)
       if double(isempty(strfind(tableData{n,2},'NA')))==0, MAT(n,1)=0; else, MAT(n,1)=double(tableData{n,2}); end
       if double(isempty(strfind(tableData{n,3},'NA')))==0, MAT(n,2)=0; else, MAT(n,2)=double(tableData{n,3}); end
       if double(isempty(strfind(tableData{n,4},'NA')))==0, MAT(n,3)=0; else, MAT(n,3)=double(tableData{n,4}); end
       if double(isempty(strfind(tableData{n,5},'NA')))==0, MAT(n,4)=0; else, MAT(n,4)=double(tableData{n,5}); end
    end
    handles.save_selection = MAT;
handles.output = hObject;
guidata(hObject, handles);
close(handles.figure1);
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    tableData = get(handles.uitable1, 'data');

    if handles.con_all_data == 0
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,2},'NA')))==1 
               tableData{n,2} = true; 
           end
           if double(isempty(strfind(tableData{n,3},'NA')))==1 
               tableData{n,3} = true;
           end
           if double(isempty(strfind(tableData{n,4},'NA')))==1 
               tableData{n,4} = true;
           end
           if double(isempty(strfind(tableData{n,5},'NA')))==1
               tableData{n,5} = true;
           end
        end
        handles.con_all_data = 1;
    else
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,2},'NA')))==1 
               tableData{n,2} = false; 
           end
           if double(isempty(strfind(tableData{n,3},'NA')))==1 
               tableData{n,3} = false;
           end
           if double(isempty(strfind(tableData{n,4},'NA')))==1 
               tableData{n,4} = false;
           end
           if double(isempty(strfind(tableData{n,5},'NA')))==1
               tableData{n,5} = false;
           end
        end
        handles.con_all_data = 0;
    end
    
    
    set(handles.uitable1, 'data',tableData);
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tableData = get(handles.uitable1, 'data');
    
    
    if handles.con_xls_data == 0
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,2},'NA')))==1 
               tableData{n,2} = true; 
           end
        end
        handles.con_xls_data = 1;
    else
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,2},'NA')))==1 
               tableData{n,2} = false; 
           end
        end
        handles.con_xls_data = 0;
    end

    set(handles.uitable1, 'data',tableData);
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tableData = get(handles.uitable1, 'data');
    

    if handles.con_mat_data == 0
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,3},'NA')))==1 
               tableData{n,3} = true; 
           end
        end
        handles.con_mat_data = 1;
    else
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,3},'NA')))==1 
               tableData{n,3} = false; 
           end
        end
        handles.con_mat_data = 0;
    end

    set(handles.uitable1, 'data',tableData);
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tableData = get(handles.uitable1, 'data');

    if handles.con_vti_data == 0
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,4},'NA')))==1 
               tableData{n,4} = true; 
           end
        end
        handles.con_vti_data = 1;
    else
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,4},'NA')))==1 
               tableData{n,4} = false; 
           end
        end
        handles.con_vti_data = 0;
    end

    set(handles.uitable1, 'data',tableData);
handles.output = hObject;
guidata(hObject, handles);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tableData = get(handles.uitable1, 'data');

    if handles.con_vtu_data == 0
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,5},'NA')))==1 
               tableData{n,5} = true; 
           end
        end
        handles.con_vtu_data = 1;
    else
        for n=1:size(tableData,1)
           if double(isempty(strfind(tableData{n,5},'NA')))==1 
               tableData{n,5} = false; 
           end
        end
        handles.con_vtu_data = 0;
    end

    set(handles.uitable1, 'data',tableData);
handles.output = hObject;
guidata(hObject, handles);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set(handles.figure1, 'Units', 'pixels');
    FigPos = get(handles.figure1, 'Position');
    set(handles.figure1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.figure1, 'Units', 'normalized');

    set(handles.uitable1, 'Units', 'pixels');
    AA = get(handles.uitable1,'Position');
    FigPos = get(handles.uitable1, 'Position');
    set(handles.uitable1, 'Position', [FigPos(1:2),FigPos(3),FigPos(4) ]);
    set(handles.uitable1, 'Units', 'normalized');

    set(handles.uitable1,'visible','on','ColumnWidth',{round(AA(3)/2) round(AA(3)/2)*0.2 round(AA(3)/2)*0.2 round(AA(3)/2)*0.2 round(AA(3)/2)*0.2})
    
    
    set(handles.text1,'FontUnits','Normalized','FontSize',0.60)
    set(handles.uitable1,'FontUnits','Normalized','FontSize',0.033)
    set(handles.pushbutton1,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton2,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton3,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton4,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton5,'FontUnits','Normalized','FontSize',0.4)
    set(handles.pushbutton6,'FontUnits','Normalized','FontSize',0.4)
    
handles.output = hObject;  
guidata(hObject, handles);
