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
#ifndef DBLDEF_FORWARD_2D_SQUARE_ELEMENT
#define DBLDEF_FORWARD_2D_SQUARE_ELEMENT

#include "MeshDataBrickElement.h"
#include "Forward2D.h"

// Class of 2D forward calculation by using square element
class Forward2DSquareElement : public Forward2D {

public:
	// Constructer
	explicit Forward2DSquareElement( const int planeID, const int iPol );

	// Destructer
	~Forward2DSquareElement();

	// Calculate EM fields of boundary planes by 2D forward calculcation
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ) = 0;

protected:
	
	// Calculate width of element
	double calcWidth( const int iElem, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate height of element
	double calcHeight( const int iElem, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate element division number of horizontal direction
	int calcNumElemHorizontal( const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Calculate element division number of vertical direction
	int calcNumElemVertical( const MeshDataBrickElement* const pMeshDataBrickElement ) const;

	// Output EM field and responses of boundary planes obtained by 2D analysis
	void output2DResult( const int type, const double freq, const int nElem, const int numElemW, const MeshDataBrickElement* const pMeshDataBrickElement ) const;

private:

	// Defailt constructer
	Forward2DSquareElement();

	// Copy constructer
	Forward2DSquareElement(const Forward2DSquareElement& rhs);

	// Copy assignment operator
	Forward2DSquareElement& operator=(const Forward2DSquareElement& rhs);

	// Calculate parameter V of Rodi(1976)
	virtual std::complex<double> calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate parameter J of Rodi(1976)
	virtual std::complex<double> calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate parameter I of Rodi(1976)
	virtual std::complex<double> calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Ex
	virtual std::complex<double> calcEx( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Ey
	virtual std::complex<double> calcEy( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Ez
	virtual std::complex<double> calcEz( const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Hx
	virtual std::complex<double> calcHx( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Hy
	virtual std::complex<double> calcHy( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

	// Calculate Hz
	virtual std::complex<double> calcHz( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const = 0;

};

#endif
