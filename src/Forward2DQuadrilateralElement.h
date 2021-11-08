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
#ifndef DBLDEF_FORWARD_2D_QUAD_ELEMENT
#define DBLDEF_FORWARD_2D_QUAD_ELEMENT

#include "Forward2D.h"
#include "MeshDataNonConformingHexaElement.h"
#include <utility>
#include <map>

// Class of 2D forward calculation by using quadrilateral element
class Forward2DQuadrilateralElement : public Forward2D {

public:

	// Constructer
	explicit Forward2DQuadrilateralElement( const int planeID, const int iPol );

	// Destructer
	~Forward2DQuadrilateralElement();

	// Calculate EM fields of boundary planes by 2D forward calculcation with 1st order nodal element
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataNonConformingHexaElement* const pMeshData ) = 0;

protected:
	
	// Total node number of equation
	int m_numEquations;

	// Total node number of equation after degenerated
	int m_numEquationsDegenerated;

	// Total node number of 2D mesh
	int m_numNodeTotal2D;

	// Array converting global node IDs of 3D mesh to the ones of 2D mesh
	std::map<int,int> m_nodeIDs3DTo2D;

	// Output EM field and responses of boundary planes obtained by 2D analysis
	void output2DResult( const double freq, const MeshDataNonConformingHexaElement* const pMeshData ) const;

private:

	// Defailt constructer
	Forward2DQuadrilateralElement();

	// Copy constructer
	Forward2DQuadrilateralElement(const Forward2DQuadrilateralElement& rhs);

	// Copy assignment operator
	Forward2DQuadrilateralElement& operator=(const Forward2DQuadrilateralElement& rhs);

	// Calculate Ex
	virtual std::complex<double> calcEx( const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;

	// Calculate Ey
	virtual std::complex<double> calcEy( const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;

	// Calculate Ez
	virtual std::complex<double> calcEz( const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;

	// Calculate Hx
	virtual std::complex<double> calcHx( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;

	// Calculate Hy
	virtual std::complex<double> calcHy( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;

	// Calculate Hz
	virtual std::complex<double> calcHz( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataNonConformingHexaElement* const pMeshData ) const = 0;
	
};

#endif
