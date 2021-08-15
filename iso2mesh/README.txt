----------------------------------------------------------------------
= Iso2Mesh: An Image-based 3D Surface and Volumetric Mesh Generator  =
----------------------------------------------------------------------

*Author: Qianqian Fang <q.fang at neu.edu>
      Department of Bioengineering
      Northeastern University
      360 Huntington Ave, Boston, MA 02115
*Version: 1.9.0 (Century Egg)
*License: GPL v2 or later (see COPYING) 
      (this license does not cover the binaries under the bin/ 
       directory, see Section III for more details)
*URL: http://iso2mesh.sf.net


== Table of Content ==
<toc>


== # Introduction ==

"Iso2Mesh" is a MATLAB/Octave-based mesh generation toolbox,
designed for easy creation of high quality surface and 
tetrahedral meshes from 3D volumetric images. It contains 
a rich set of mesh processing scripts/programs, working 
either independently or interacting with external free 
meshing utilities. Iso2Mesh toolbox can directly convert
a 3D image stack, including binary, segmented or gray-scale 
images such as MRI or CT scans, into quality volumetric 
meshes. This makes it particularly suitable for multi-modality 
medical imaging data analysis and multi-physics modeling.
Above all, iso2mesh is open-source. You can download it for 
free. You are also allowed to extend the toolbox for your
own research and share with other users. Iso2Mesh is 
cross-platform and is compatible with both MATLAB and GNU Octave 
(a free MATLAB clone).

The details of this toolbox can be found in the following
papers:

*Anh Phong Tran and Qianqian Fang, "Fast and high-quality tetrahedral \
 mesh generation from neuroanatomical scans," arXiv preprint arXiv:1708.08954, 2017
*Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary and \
 gray-scale images," Proceedings of IEEE International Symposium on Biomedical Imaging \
 (ISBI 2009), pp. 1142-1145, 2009


== # Overview of the functions ==

Creation of high-quality surface and tetrahedral meshes 
from volumetric images has been a challenging task. 
There are very limited software and resources available 
for this purpose. Commercial tools, such as Mimics 
and Simpleware, are both expensive and limited in flexibility. 
Iso2Mesh was developed as a free alternative to these 
expensive commercial tools and provides researchers a highly
flexible, modular and streamlined image-based mesh 
generation pipeline. Intuitive interfaces and rich
functionalities allow one to enjoy a wide range 
of mesh-based analyses, ranging from 3D volumetric 
image pre-processing (hole-filling, thinning and 
thickening), surface mesh modeling (extraction, 
remeshing, repairing, and smoothing) to volumetric 
mesh creation.

Converting 3D image stacks into quality surface and
tetrahedral meshes is one of the core features of iso2mesh.
We provide serveral automated functions to perform
the image->mesh and mesh->image conversion, including

* vol2mesh (v2m): convert a 3D volumetric image into a tetrahedral mesh
* vol2surf (v2s): extract triangular surfaces from a 3D image volume
* surf2mesh (s2m): create a tetrahedral mesh from a triangular surface mesh
* surf2vol (s2v): rasterize a close-surface to a volumetric image
* mesh2vol (m2v): rasterize a tetrahedral mesh to a volumetric image

Most of these function are associated with several meshing
options and parameters to give users full control to mesh 
density, adaptivity, region labeling and mesh quality.
The output data for some of these functions can be used
as the input for the others, giving endless combinations to
analyze your data. In addition to image-based mesh generation,
iso2mesh can also mesh geometry primitives such as spheres,
cubes and cylinders. This makes iso2mesh a CAD-capable software,
fully integrated in the MATLAB/Octave environments.

Another core feature of iso2mesh is surface mesh
processing. A surface mesh is the bridge between a voxelated
image and a tetrahedral mesh, and is the foundation for 
successful 3D mesh generation. In iso2mesh, we provide the
following key functions for surface mesh processing:

* smoothsurf (sms): smoothing a surface mesh
* surfboolean: boolean operations (join, intersect, diff) of two surfaces
* meshresample: downsample a surface mesh to a coarse surface
* meshrefine: refine a surface mesh
* remeshsurf: remesh a surface mesh, including up-sampling and down-sampling
* meshcheckrepair: test a surface mesh and remove defects and self-intersecting elements

A third area iso2mesh excels is the rich set of mesh analysis 
and inquiry functions. For both surface and tetrahedral meshes, 
iso2mesh can report the element volume (or area), mesh quality
metrics, node connectivity and neighbors, interior edges and 
boundaries, surface norms, element centroids, etc.

Last, but not the least, iso2mesh can import and export mesh
data from/to a dozen of file formats, including those used by the 
popular FEA software like COMSOL and ABAQUS. The STL format can
export an iso2mesh-generated mesh for 3D printing.
Users can not only export the data to a file, but also make
3D plots in MATLAB/Octave using the powerful "plotmesh" function.

These examples are only a small fraction in the over 100 optimized 
functions provided by iso2mesh. The modular design of iso2mesh 
makes the code easy to understand and easy to be intergrated in 
your data analysis workflow. Please check out the full iso2mesh 
function list and detailed help information in the following URL:

http://iso2mesh.sf.net/cgi-bin/index.cgi?Doc/FunctionList

== # Compiling Iso2Mesh ==

The default release of Iso2Mesh packages already contains pre-compiled
binaries for a wide range of platforms (32/64bit Windows, 32/64bit Linux
and Mac with 64bit Intel and 32 bit PowerPC CPUs). So, without needing 
to recompile, Iso2Mesh can be executed out-of-box on MATLAB or GNU Octave.

However, in the event that your operating system is not supported, or
due to license restrictions, such as creating a release for various
Linux distributions, you can recreate the mesh utility binaries under
iso2mesh/bin folder by following the below commands:

 git clone --recurse-submodules https://github.com/fangq/iso2mesh.git
 cd iso2mesh
 rm -rf bin/*.mex* bin/*.exe
 cd tools
 make clean
 make

This will download and recompile the below binaries in the bin folder:

 cgalmesh
 cgalsurf
 cgalsimp2
 jmeshlib
 meshfix
 tetgen1.5
 cork

Once these binary files are recreated, you can run all the major functionalities
of Iso2Mesh. The gtsset and gtrefine tools are depreciated and replaced by 
cork and tetgen.

To compile the above external tools, the below tools must be pre-installed
(tested on Ubuntu 14.04 LTS, if you use another Linux distribution, the package
names might be different)

 libcgal-dev
 clang
 cmake

you can install these by

 sudo apt-get install libcgal-dev clang cmake
 
on Ubuntu.

== # Acknowledgement ==

This toolbox interacts with a number external meshing tools 
to perform the essential functionalities. These tools are listed 
below:

=== bin/tetgen and bin/tetgen1.5 ===

*Summary:tetgen is a compact and fast 3D mesh generator
*License: GNU Affero General Public License version 3
*URL: http://tetgen.org/
*Author: Hang Si <si at wias-berlin.de>
::Research Group: Numerical Mathematics and Scientific Computing
::Weierstrass Institute for Applied Analysis and Stochastics
::Mohrenstr. 39, 10117 Berlin, Germany

=== bin/cgalsurf ===

*Summary: cgalsurf is a utility to extract a surface mesh from \
a gray-scale or a binary 3D image
*Source: it is a slightly modified version from Surface_mesher \
example from CGAL 3.4
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.3/doc_html/cgal_manual/Surface_mesher/Chapter_main.html

=== bin/cgalmesh and bin/cgalpoly ===

*Summary: cgalmesh and cgalpoly are utilities to produce surface \
and volumetric meshes from a multi-valued volumetric image
*Source: it is a slightly modified version from Mesh_3 \
example from CGAL 3.5
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.5/doc_html/cgal_manual/Mesh_3/Chapter_main.html

=== bin/cgalsimp2 ===

*Summary: cgalsimp2 performs surface mesh simplification in iso2mesh.
*Source: it is adapted from Surface_mesh_simplification example of CGAL library
*License: CGAL v3 core library is licensed under QPL (Q Public License) \
other modules are under the Lesser General Public License (LGPL)
*URL: http://www.cgal.org/Manual/3.4/doc_html/cgal_manual/Surface_mesh_simplification/Chapter_main.html

=== bin/jmeshlib ===

*Summary: meshfix is adapted from the sample code of JMeshLib
*License: GPL (GNU General Public License) v2 or later
*URL:http://jmeshlib.sourceforge.net/
*Author:Marco Attene <attene at ge.imati.cnr.it>
::Istituto di Matematica Applicata e Tecnologie Informatiche
::Consiglio Nazionale delle Ricerche
::Via De Marini, 6 (Torre di Francia)
::16149 Genoa - ITALY 

=== bin/meshfix ===

*Summary: meshfix is a mesh-repairing utility
*License: GPL (GNU General Public License) v2 or later
*URL:http://code.google.com/p/meshfix/
*Author: Marco Attene, Mirko Windhoff and Axel Thielscher.
::Istituto di Matematica Applicata e Tecnologie Informatiche
::Consiglio Nazionale delle Ricerche
::Via De Marini, 6 (Torre di Francia)
::16149 Genoa - ITALY 

=== bin/cork ===

*Summary: A robust surface mesh Boolean operation algorithm
*License: LGPL (GNU Lesser General Public License)
*URL:https://github.com/gilbo/cork
*Author: Gilbert Bernstein

=== bin/gtsset and bin/gtsrefine ===

*Summary: GTS is the GNU Triangulated Surface Library
*License: LGPL (GNU Lesser General Public License)
*URL:http://gts.sourceforge.net/
*Author: GTS developers


Note: Iso2Mesh and the above meshing utilities are considered 
as an "aggregate" rather than "derived work", based on the 
definitions in GPL FAQ (http://www.gnu.org/licenses/gpl-faq.html#MereAggregation)
Therefore, the license of Iso2Mesh and these utilities are independent.
The Iso2Mesh license only applies to the scripts and documentation/data
in this package and exclude those programs stored in the bin/ directory.
The source codes of the modified meshing utilities are available
separately at Iso2Mesh's website and retain their upstream licenses.

Your acknowledgement of Iso2Mesh in your publications or 
presentations would be greatly appreciated by the author of 
this toolbox. The citation information can be found in the
Introduction section.
