//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 Yoshiya Usui
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//-------------------------------------------------------------------------------------------------------
#ifndef DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_EDGE_BASED
#define DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_EDGE_BASED

#include "Forward2DTriangleElement.h"

// Class of 2D forward calculation by using edge based element with triangle shape
class Forward2DTriangleElementEdgeBased : public Forward2DTriangleElement {

public:

	// Constructer
	explicit Forward2DTriangleElementEdgeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DTriangleElementEdgeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ) = 0;

protected:

	struct JacobianMatrix{
		double jacob11;
		double jacob12;
		double jacob21;
		double jacob22;
	};

	// Type of outer edge
	enum TypeOfOuterEdge{
		INNER_EDGE = -1,
		LOWER_EDGE = 0,
		UPPER_EDGE,
		LEFT_EDGE,
		RIGHT_EDGE,
	};

	//// Array converting global element IDs of 2D mesh to global edge IDs of 2D mesh
	//EdgeIDsOfTriangleElement* m_elemID2EdgeIDGlobal;

	// Array converting local IDs to the flag specifing whether direction of local vector basis function
	bool** m_signInversion;

	// Array edge IDs of 2D mesh to non-zero electric field values specified
	std::map< int, std::complex<double> > m_edgesIDGlobal2NonZeroValues;

	// Calculate Ex
	virtual std::complex<double> calcEx( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Ey
	virtual std::complex<double> calcEy( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Ez
	virtual std::complex<double> calcEz( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Hx
	virtual std::complex<double> calcHx( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Hy
	virtual std::complex<double> calcHy( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Hz
	virtual std::complex<double> calcHz( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

private:

	// Calculate horizontal electric field
	virtual std::complex<double> calcValueElectricFieldHorizontal( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate vertical electric field
	virtual std::complex<double> calcValueElectricFieldVertical( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate magnetic field perpendicular to the boundary plane
	virtual std::complex<double> calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Defailt constructer
	Forward2DTriangleElementEdgeBased();

	// Copy constructer
	Forward2DTriangleElementEdgeBased(const Forward2DTriangleElementEdgeBased& rhs);

	// Copy assignment operator
	Forward2DTriangleElementEdgeBased& operator=(const Forward2DTriangleElementEdgeBased& rhs);

};

#endif
