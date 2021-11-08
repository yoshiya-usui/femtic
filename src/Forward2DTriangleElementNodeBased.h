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
#ifndef DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_NODE_BASED
#define DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_NODE_BASED

#include "Forward2DTriangleElement.h"

// Class of 2D forward calculation by using node based element with triangle shape
class Forward2DTriangleElementNodeBased : public Forward2DTriangleElement {

public:

	// Constructer
	explicit Forward2DTriangleElementNodeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DTriangleElementNodeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ) = 0;

protected:

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

	// Defailt constructer
	Forward2DTriangleElementNodeBased();

	// Copy constructer
	Forward2DTriangleElementNodeBased(const Forward2DTriangleElementNodeBased& rhs);

	// Copy assignment operator
	Forward2DTriangleElementNodeBased& operator=(const Forward2DTriangleElementNodeBased& rhs);


};

#endif
