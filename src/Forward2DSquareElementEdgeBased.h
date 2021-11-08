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
#ifndef DBLDEF_FORWARD_2D_SQUARE_ELEMENT_EDGE_BASED
#define DBLDEF_FORWARD_2D_SQUARE_ELEMENT_EDGE_BASED

#include "MeshDataBrickElement.h"
#include "Forward2DSquareElement.h"

// Class of 2D forward calculation by using square edge element
class Forward2DSquareElementEdgeBased : public Forward2DSquareElement {

public:

	// Constructer
	explicit Forward2DSquareElementEdgeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DSquareElementEdgeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ) = 0;

protected:

	// Calculate Ex
	virtual std::complex<double> calcEx( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate Ey
	virtual std::complex<double> calcEy( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate Ez
	virtual std::complex<double> calcEz( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate Hx
	virtual std::complex<double> calcHx( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate Hy
	virtual std::complex<double> calcHy( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate Hz
	virtual std::complex<double> calcHz( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

private:

	// Defailt constructer
	Forward2DSquareElementEdgeBased();

	// Copy constructer
	Forward2DSquareElementEdgeBased(const Forward2DSquareElementEdgeBased& rhs);

	// Copy assignment operator
	Forward2DSquareElementEdgeBased& operator=(const Forward2DSquareElementEdgeBased& rhs);

	// Calculate parameter V of Rodi(1976)
	virtual std::complex<double> calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate parameter J of Rodi(1976)
	virtual std::complex<double> calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate parameter I of Rodi(1976)
	virtual std::complex<double> calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate horizontal electric field
	virtual std::complex<double> calcValueElectricFieldHorizontal( const int iElem, const double wCoord, const double hCoord ) const = 0;

	// Calculate vertical electric field
	virtual std::complex<double> calcValueElectricFieldVertical( const int iElem, const double wCoord, const double hCoord ) const = 0;

	// Calculate magnetic field perpendicular to the boundary plane
	virtual std::complex<double> calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

};

#endif
