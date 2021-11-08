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
#ifndef DBLDEF_FORWARD_2D_SQUARE_ELEMENT_0TH_ORDER_EDGE_BASED
#define DBLDEF_FORWARD_2D_SQUARE_ELEMENT_0TH_ORDER_EDGE_BASED

#include "MeshDataBrickElement.h"
#include "Forward2DSquareElementEdgeBased.h"

// Class of 2D forward calculation by using 0th order edge element with square shape
class Forward2DSquareElement0thOrderEdgeBased : public Forward2DSquareElementEdgeBased {

public:

	// Constructer
	explicit Forward2DSquareElement0thOrderEdgeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DSquareElement0thOrderEdgeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement );

private:

	// Defailt constructer
	Forward2DSquareElement0thOrderEdgeBased();

	// Copy constructer
	Forward2DSquareElement0thOrderEdgeBased(const Forward2DSquareElement0thOrderEdgeBased& rhs);

	// Copy assignment operator
	Forward2DSquareElement0thOrderEdgeBased& operator=(const Forward2DSquareElement0thOrderEdgeBased& rhs);

	// Get shape functions of horizontal direction for 0th order edge-based elements
	inline double getShapeFuncHorizontal0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const;

	// Get shape functions of vertical direction for 0th order edge-based elements
	inline double getShapeFuncVertical0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const;

	// Get shape functions rotated for 0th order edge-based elements
	inline double getShapeFuncRotated0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const;

	// Calculate horizontal electric field values for 0th order edge-based elements
	virtual std::complex<double> calcValueElectricFieldHorizontal( const int iElem, const double wLocal, const double hLocal ) const;

	// Calculate vertical electric field values for 0th order edge-based elements
	virtual std::complex<double> calcValueElectricFieldVertical( const int iElem, const double wLocal, const double hLocal ) const;

	//// Calculate magnetic field values for 0th order edge-based elements
	//std::complex<double> calcValueMagneticField( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate magnetic field perpendicular to the boundary plane
	virtual std::complex<double> calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter V of Rodi(1976)
	virtual std::complex<double> calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter J of Rodi(1976)
	virtual std::complex<double> calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter I of Rodi(1976)
	virtual std::complex<double> calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;


};

#endif
