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
#ifndef DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT
#define DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT

#include "Forward2D.h"
#include "MeshDataTetraElement.h"
#include <utility>
#include <map>

// Class of 2D forward calculation by using triangle element
class Forward2DTriangleElement : public Forward2D {

public:

	// Constructer
	explicit Forward2DTriangleElement( const int planeID, const int iPol );

	// Destructer
	~Forward2DTriangleElement();

	// Calculate EM fields of boundary planes by 2D forward calculcation with 1st order nodal element
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ) = 0;

	//// Count total number of node on plane
	//int countTotalNodeNumberOnPlane( const MeshDataTetraElement* const pMeshDataTetraElement );

	//// Count outer edges on plane
	//int countOuterEdgesOnPlane( const MeshDataTetraElement* const pMeshDataTetraElement ) const;

protected:

	//std::vector< std::pair< std::pair<int, int>, int> > m_globalEdgeID2GlobalNodePairAndElemID;
	
	//struct NodeIDsOfTriangleElement{
	//	int node0;
	//	int node1;
	//	int node2;
	//};

	// Total node number of equation
	int m_numEquations;

	// Total node number of equation after degenerated
	int m_numEquationsDegenerated;

	// Total node number of 2D mesh
	int m_numNodeTotal2D;

	// Array converting global node IDs of 3D mesh to the ones of 2D mesh
	std::map<int,int> m_nodeIDs3DTo2D;

	//// Calculate array converting global node IDs of 3D mesh to the ones of 2D mesh
	//void calcArrayConvertingNodeIDs3DTo2D( const MeshDataTetraElement* const pMeshDataTetraElement );

	// Output EM field and responses of boundary planes obtained by 2D analysis
	void output2DResult( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

private:

	// Defailt constructer
	Forward2DTriangleElement();

	// Copy constructer
	Forward2DTriangleElement(const Forward2DTriangleElement& rhs);

	// Copy assignment operator
	Forward2DTriangleElement& operator=(const Forward2DTriangleElement& rhs);

	//// Number global edge IDs
	//void numberGlobalEdgesIDs( const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate Ex
	virtual std::complex<double> calcEx( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate Ey
	virtual std::complex<double> calcEy( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate Ez
	virtual std::complex<double> calcEz( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate Hx
	virtual std::complex<double> calcHx( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate Hy
	virtual std::complex<double> calcHy( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;

	// Calculate Hz
	virtual std::complex<double> calcHz( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const = 0;
	
};

#endif
