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
#ifndef DBLDEF_FORWARD_2D_SQUARE_ELEMENT_1ST_ORDER_NODE_BASED
#define DBLDEF_FORWARD_2D_SQUARE_ELEMENT_1ST_ORDER_NODE_BASED

#include "MeshDataBrickElement.h"
#include "Forward2DSquareElementNodeBased.h"

// Class of 2D forward calculation by using 1st order node element with square shape
class Forward2DSquareElement1stOrderNodeBased : public Forward2DSquareElementNodeBased {

public:

	// Constructer
	explicit Forward2DSquareElement1stOrderNodeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DSquareElement1stOrderNodeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement );

	//// Get electric field perpendicular to bondary plane from element ID and coordinate values
	//virtual std::complex<double> getElectricFieldPerpendicularToPlane( const int iElem, const double wLocal, const double hLocal ) const;

private:

	// Defailt constructer
	Forward2DSquareElement1stOrderNodeBased();

	// Copy constructer
	Forward2DSquareElement1stOrderNodeBased(const Forward2DSquareElement1stOrderNodeBased& rhs);

	// Copy assignment operator
	Forward2DSquareElement1stOrderNodeBased& operator=(const Forward2DSquareElement1stOrderNodeBased& rhs);

	// Get shape functions for 1st order node-based elements
	inline double getShapeFunc1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const;

	// Get shape functions defferentiated by local coordinate of w for 1st order node-based elements
	inline double getShapeFuncDiffByWLocal1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const;

	// Get shape functions defferentiated by local coordinate of h for 1st order node-based elements
	inline double getShapeFuncDiffByHLocal1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const;

	//// Calculate parameter V of Rodi(1976) for 1st order node-based elements
	//std::complex<double> calcValueV1stOrder( const int iElem, const double wLocal, const double hLocal ) const;

	//// Calculate parameter J of Rodi(1976) for 1st order node-based elements
	//std::complex<double> calcValueJ1stOrder( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	//// Calculate parameter I of Rodi(1976) for 1st order node-based elements
	//std::complex<double> calcValueI1stOrder( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter V of Rodi(1976)
	virtual std::complex<double> calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter J of Rodi(1976)
	virtual std::complex<double> calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate parameter I of Rodi(1976)
	virtual std::complex<double> calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate electric field perpendicular to the boundary plane
	virtual std::complex<double> calcValueElectricFieldPerpendicular( const int iElem, const double wCoord, const double hCoord ) const;

	// Calculate horizontal magnetic field
	virtual std::complex<double> calcValueMagneticFieldHorizontal( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate vertical magnetic field
	virtual std::complex<double> calcValueMagneticFieldVertical( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

};

#endif
