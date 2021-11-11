# FEMTIC
FEMTIC is a 3-D magnetotelluric inversion code based on the following studies. 
FEMTIC was made using object-oriented programming with C++.

*Y. Usui, 3-D inversion of magnetotelluric data using unstructured tetrahedral elements: applicability to data affected by topography, Geophys. J. Int., 202 (2): 828-849, https://doi.org/10.1093/gji/ggv186, 2015.*

*Y. Usui, Y. Ogawa, K. Aizawa, W. Kanda, T. Hashimoto, T. Koyama, Y. Yamaya and T. Kagiyama, Three-dimensional resistivity structure of Asama Volcano revealed by data-space magnetotelluric inversion using unstructured tetrahedral elements, Geophys. J. Int., 208 (3): 1359-1372, https://doi.org/10.1093/gji/ggw459, 2017.*

*Y. Usui, T. Kasaya, Y. Ogawa and H. Iwamoto, Marine magnetotelluric inversion with an unstructured tetrahedral mesh, Geophys. J. Int., 214(2): 952-974, https://doi.org/10.1093/gji/ggy171, 2018.*

The website of FEMTIC:
https://sites.google.com/view/yoshiyausui/femtic

# Functional overview
FEMTIC gives a three-dimensional electrical resistivity structure from the response functions at observation points on the Earth's surface.

Mesh type: tetrahedral mesh, hexahedral brick mesh, and non-conforming deformed hexahedral mesh

Data type: impedance tensor, vertical magnetic transfer function, inter-station horizontal magnetic transfer function, phase tensor, apparent resistivity, and phase.

Parallel computation: Hybrid MPI/OpenMP parallel calculation can be used.

Model parameter: Subsurface electrical resistivity and distortion matrix of galvanic distortion.

# Release note
***v4.1*** Nov. 9, 2021 This new version supports difference filter. The error calculation of log10(apparent resistivity) is modified. Rotation angles of distortion matrix are limited to from -90 to 90 (deg.) when gains and rotations of galvanic distortion are estimated.

***v4.0*** Jun. 3, 2021 This new version supports non-conforming deformed hexahedral mesh.

***v3.5*** Jan. 11, 2021 This new version supports observed data of apparent resistivity and phase and unsymmetric roughening matrix.

***v3.4.7*** Sep. 4, 2020 The integer indices into the multiple right-hand-side vectors and solution vectors were changed from 32-bit to 64-bit.

***v3.4.6*** Sep. 2, 2020 This version allows us to make resistivity of every individual subsurface element to be a different model parameter, in analogy with other 3-D inversion code.

