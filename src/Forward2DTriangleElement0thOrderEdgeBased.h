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
#ifndef DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_0TH_ORDER_EDGE_BASED
#define DBLDEF_FORWARD_2D_TRIANGLE_ELEMENT_0TH_ORDER_EDGE_BASED

#include "Forward2DTriangleElementEdgeBased.h"

// Class of 2D forward calculation by using triangle element 0th order edge based 
class Forward2DTriangleElement0thOrderEdgeBased : public Forward2DTriangleElementEdgeBased {

public:

	// Constructer
	explicit Forward2DTriangleElement0thOrderEdgeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DTriangleElement0thOrderEdgeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation with 0tht order edge element
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshData );

private:

	static const int DIRICHLET_BOUNDARY_NONZERO_VALUE = -1;
	static const int DIRICHLET_BOUNDARY_ZERO_VALUE    = -2;

	// Number of integral points
	static const int m_numIntegralPoints = 3;

	// U coordinates of integral points
	double m_uCoord[m_numIntegralPoints];

	// V coordinates of integral points
	double m_vCoord[m_numIntegralPoints];

	// Weights of integral points
	double m_weight[m_numIntegralPoints];

	// Defailt constructer
	Forward2DTriangleElement0thOrderEdgeBased();

	// Copy constructer
	Forward2DTriangleElement0thOrderEdgeBased(const Forward2DTriangleElement0thOrderEdgeBased& rhs);

	// Copy assignment operator
	Forward2DTriangleElement0thOrderEdgeBased& operator=(const Forward2DTriangleElement0thOrderEdgeBased& rhs);

	// Calculate array converting local IDs to global ones
	void calcArrayConvertLocalID2Global( const MeshDataTetraElement* const pMeshData );

	// Get type of outer edge
	// [note] : You must confirm inputed edge is the outer edge. 
	int getTypeOfOuterEdgeOfBoundaryPlanes( const int elemIDLocal2D, const int edgeIDLocal2D,
		const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Get shape functions of the 1st direction with respect to the reference element coordinate system
	inline double getShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const int num ) const;

	// Get shape functions of the 2nd direction with respect to the reference element coordinate system
	inline double getShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const int num ) const;

	// Get shape functions rotated with respect to the reference element coordinate system
	inline double getShapeFuncRotated() const;

	// Calculate horizontal electric field
	virtual std::complex<double> calcValueElectricFieldHorizontal( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate vertical electric field
	virtual std::complex<double> calcValueElectricFieldVertical( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate magnetic field perpendicular to the boundary plane
	virtual std::complex<double> calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const;

	// Calculate jacobian matrix of the elements on the Z-X plane of boundary
	void calcJacobianMatrixOnZXPlaneOfBoundary( const MeshDataTetraElement* const pMeshDataTetraElement, const int elemID2D, Forward2DTriangleElementEdgeBased::JacobianMatrix& JacobMat, double& determinant ) const;

	// Calculate jacobian matrix of the elements on the Y-Z plane of boundary
	void calcJacobianMatrixOnYZPlaneOfBoundary( const MeshDataTetraElement* const pMeshDataTetraElement, const int elemID2D, Forward2DTriangleElementEdgeBased::JacobianMatrix& JacobMat, double& determinant ) const;

};

#endif
