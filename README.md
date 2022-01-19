# FEMTIC
FEMTIC is a 3-D magnetotelluric inversion code based on the following studies. FEMTIC was made using object-oriented programming with C++. FEMTIC enables us to incorporate topography and bathymetry into an inversion model. FEMTIC is applicable to land magnetelluric survey data as well as ocean bottom magnetelluric survey data.

*Y. Usui, 3-D inversion of magnetotelluric data using unstructured tetrahedral elements: applicability to data affected by topography, Geophys. J. Int., 202 (2): 828-849, https://doi.org/10.1093/gji/ggv186, 2015.*

*Y. Usui, Y. Ogawa, K. Aizawa, W. Kanda, T. Hashimoto, T. Koyama, Y. Yamaya and T. Kagiyama, Three-dimensional resistivity structure of Asama Volcano revealed by data-space magnetotelluric inversion using unstructured tetrahedral elements, Geophys. J. Int., 208 (3): 1359-1372, https://doi.org/10.1093/gji/ggw459, 2017.*

*Y. Usui, T. Kasaya, Y. Ogawa and H. Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophys. J. Int., 214(2): 952-974, https://doi.org/10.1093/gji/ggy171, 2018.*

**The website of FEMTIC:**
https://sites.google.com/view/yoshiyausui/femtic

# Functional overview
FEMTIC gives a three-dimensional electrical resistivity structure from the response functions at observation points on the Earth's surface.

**Mesh type**: Tetrahedral mesh / Hexahedral brick mesh / Non-conforming deformed hexahedral mesh

**Data type**: Impedance tensor / Vertical magnetic transfer function / Inter-station horizontal magnetic transfer function / Phase tensor / Apparent resistivity / Phase.

**Inversion algorithm**: Model-space Gauss-Netwon method / Data-space Gauss-Netwon method

**Parallel computation**: Multiple processes parallel computation with MPI / Multiple threads parallel computation with OpenMP / MPI & OpenMP hybrid parallel computation

**Model parameter**: Subsurface electrical resistivity / Distortion matrix of galvanic distortion

# Release note
***v4.1*** Nov. 9, 2021: This new version supports difference filter. The error calculation of log10(apparent resistivity) is modified. Rotation angles of distortion matrix are limited to from -90 to 90 (deg.) when gains and rotations of galvanic distortion are estimated.

***v4.0*** Jun. 3, 2021: This new version supports non-conforming deformed hexahedral mesh.

***v3.5*** Jan. 11, 2021: This new version supports observed data of apparent resistivity and phase and unsymmetric roughening matrix.

***v3.4.7*** Sep. 4, 2020: The integer indices into the multiple right-hand-side vectors and solution vectors were changed from 32-bit to 64-bit.

***v3.4.6*** Sep. 2, 2020: This version allows us to make resistivity of every individual subsurface element to be a different model parameter, in analogy with other 3-D inversion code.

# Pre/post-processing tools for FEMTIC
Some pre/post-processing tools for FEMTIC, including meshing tools, and their manuals can be downloaded from GitHub. Results of FEMTIC can be visualized by [ParaView](https://www.paraview.org/).

[makeTetraMesh](https://github.com/yoshiya-usui/makeTetraMesh.git): By using this tool, you can make a surface mesh for creating a tetrahedral mesh.

[makeMtr](https://github.com/yoshiya-usui/makeMtr.git): This tool output .mtr file of TetGen by reading node and .ele files of TetGen.

[TetGen2Femtic](https://github.com/yoshiya-usui/TetGen2Femtic.git): This program converts output files of TetGen to FEMTIC.

[makeDHexaMesh](https://github.com/yoshiya-usui/makeDHexaMesh.git): Tool for making non-conforming deformed hexahedral mesh for FEMTIC

[makeHexaMesh](https://github.com/yoshiya-usui/makeHexaMesh.git): Tool for making hexahedral brick mesh for FEMTIC

[mergeResultOfFEMTIC](https://github.com/yoshiya-usui/mergeResultOfFEMTIC.git): By this program, you can merge result files (.csv) of FEMTIC.

[makeCutawayForGMT](https://github.com/yoshiya-usui/makeCutawayForGMT.git): By using this program, you can make a file needed to draw a cross-section of a resistivity structure by GMT.

[changeResistivityForFemtic](https://github.com/yoshiya-usui/changeResistivityForFemtic.git): By this program, you can change resistivity values of a specified area for sensitivity tests of FEMTIC inversion results
