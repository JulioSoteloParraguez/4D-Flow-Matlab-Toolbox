# 4D-Flow-Matlab-Toolbox

![fcvm-09-885338-g001](https://user-images.githubusercontent.com/75703824/196751007-3f58efda-eefc-4b39-ac8a-b40c3217af2b.jpg)


This code was developed by Dr. Julio Sotelo, as part of his PhD thesis and also his postdoctoral research, guided by the professors Dr. Sergio Uribe and Dr. Daniel Hurtado, from the Pontificia Universidad Católica de Chile. It is also part of the contribution developed by the Center for Biomedical Imaging (www.mri.cl) and Nucleo Milenio CardioMR (https://cardiomr.cl/). We have developed a methodology for the non-invasive quantification of hemodynamics and geometrical parameters from 4D flow data sets based on Finite Elements methods. [Reference 1 to 6]

Also we include aditional Matlab-based toolkits to extract the pressure differences from 4D flow MRI images, developed by colaborators at Massachusetts Institute of Technology (MIT) and King's College London (KCL).

The application run on Windows and macOS, it can also be used in Linux. To generate the FEM mesh we make use of iso2mesh opensource toolbox (http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home). [Reference 7 and 8]

The 4D Flow MRI data set need to be loaded as MATLAB structure file (See the document USER-GUIDE-4D-FLOW-APP). 

The data structure must have the following format:

a) The file need to be named data.mat, and saved in a single folder.

b) The data.mat file contain the following:

  - data.MR_FFE_FH = magnitude image (4D matrix with: rows, columns, slices, cardiac phases)
  - data.MR_FFE_AP = magnitude image (4D matrix with: rows, columns, slices, cardiac phases)
  - data.MR_FFE_RL = magnitude image (4D matrix with: rows, columns, slices, cardiac phases)
  - data.MR_PCA_FH = velocity image in cm/s (4D matrix with: rows, columns, slices, cardiac phases)
  - data.MR_PCA_AP = velocity image in cm/s (4D matrix with: rows, columns, slices, cardiac phases)
  - data.MR_PCA_RL = velocity image in cm/s (4D matrix with: rows, columns, slices, cardiac phases)
  - data.voxel_MR = voxel size (row,columns,slices) in mm
  - data.VENC = velocity encoding in cm/s
  - data.heart_rate = cardiac frequency in bpm
  - data.type = you can write 'DCM' in this variable.
  - if you need more information to save you can add more data.XXX variables. 
  
c) unzip the zip file (iso2mesh), for the system version that you are using.

d) To excecute the app go to the app folder and write "run GUIDE_4D_FLOW.m" in the MATLAB command windows.

To include the vWERP module: You can contact to David Marlevi (marlevi@mit.edu), the author of this method and request the "vWERP/" folder that has to go inside the 4D flow APP folder.

To obtain more information about SAW method, please visit (http://cmib.website/resources/#PressMapTk).

If you have some problems to create this structure file, or if you need more assistance, please contact me to the email: julio.sotelo@usm.cl

# References APP and Hemodynamics Parameters

Matlab-based toolkit to extract the hemodynamics and geomtrical parameters from 4D-flow MRI acquisitions:

>1.- Sotelo J, Mura J, Hurtado DE, Uribe S. A NOVEL MATLAB TOOLBOX FOR PROCESSING 4D FLOW MRI DATA. Proc. Intl. Soc. Mag. Reson. Med. 27 (2019)

>2.- Sotelo J, Urbina J, Valverde I, Tejos C, Irarrazaval P, Andia ME, Uribe S, Hurtado DE. 3D Quantification of Wall Shear Stress and Oscillatory Shear Index Using a Finite-Element Method in 3D CINE PC-MRI Data of the Thoracic Aorta. IEEE Trans Med Imaging. 2016 Jun;35(6):1475-87. doi: 10.1109/TMI.2016.2517406. Epub 2016 Jan 14. PMID: 26780787.

>3.- Sotelo J, Urbina J, Valverde I, Mura J, Tejos C, Irarrazaval P, Andia ME, Hurtado DE, Uribe S. Three-dimensional quantification of vorticity and helicity from 3D cine PC-MRI using finite-element interpolations. Magn Reson Med. 2018 Jan;79(1):541-553. doi: 10.1002/mrm.26687. Epub 2017 Mar 31. PMID: 28370386.

>4.- Sotelo J, Dux-Santoy L, Guala A, Rodríguez-Palomares J, Evangelista A, Sing-Long C, Urbina J, Mura J, Hurtado DE, Uribe S. 3D axial and circumferential wall shear stress from 4D flow MRI data using a finite element method and a laplacian approach. Magn Reson Med. 2018 May;79(5):2816-2823. doi: 10.1002/mrm.26927. Epub 2017 Oct 4. PMID: 28980342.

>5.- Sotelo J, Valverde I, Martins D, Bonnet D, Boddaert N, Pushparajan K, Uribe S, Raimondi F. Impact of aortic arch curvature in flow haemodynamics in patients with transposition of the great arteries after arterial switch operation. Eur Heart J Cardiovasc Imaging. 2021 Jan 31:jeaa416. doi: 10.1093/ehjci/jeaa416. Epub ahead of print. PMID: 33517430.

>6.- Sotelo J, Bissell MM, Jiang Y, Mella H, Mura J, Uribe S. Three-Dimensional quantification of circulation using finite-element methods in 4D flow MR data of the thoracic aorta. Magn Reson Med. 2021 (In Press). doi: 10.1002/MRM.29004

# References for iso2mesh toolbox

Matlab-based toolkit to generate tetrahedral finite element mesh:

>7.-Anh Phong Tran, Shijie Yan and Qianqian Fang*, (2020) "Improving model-based fNIRS analysis using mesh-based anatomical and light-transport models," Neurophotonics, 7(1), 015008

>8.- Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary and gray-scale images," Proceedings of IEEE International Symposium on Biomedical Imaging 2009, pp. 1142-1145, 2009

# References for vWERP method (work-energy based method)

Matlab-based toolkit to extract the pressure differences from 4D-flow MRI acquisitions:

>9.- Marlevi, D., Ruijsink, B., Balmus, M. et al. Estimation of Cardiovascular Relative Pressure Using Virtual Work-Energy. Sci Rep 9, 1375 (2019). https://doi.org/10.1038/s41598-018-37714-0

# References for SAW method (simplified advective WERP formulation)

Matlab-based toolkit to extract the spatio-temporal maps of pressure differences from 4D-flow MRI acquisitions, and metrics such as:

>10.-	The conduit and reservoir function of arteries (e.g. the aorta), as described in Unlocking the Non-invasive Assessment of Conduit and Reservoir Function in the Aorta (https://doi.org/10.1007/s12265-022-10221-4)

>11.-	The distance to pressure recovery, as described in Non-invasive assessment of pressure recovery distance after aortic valve stenosis (in review)

>12.-	Pressure differences accordingly to the Work Energy Relative Pressure (WERP) method, as described in Non-invasive pressure difference estimation from PC-MRI using the work-energy equation. (https://doi.org/10.1016/j.media.2015.08.012)

>13.-	Pressure differences accordingly to the Virtual WERP (v-WERP) method, as described in Non-invasive estimation of relative pressure for intracardiac flows using virtual work-energy. (https://doi.org/10.1016/j.media.2020.101948)


# Networks

Google Scholar: https://scholar.google.com/citations?user=aqhjvOQAAAAJ&hl=es&authuser=3

Researchgate: https://www.researchgate.net/profile/Julio_Sotelo3
