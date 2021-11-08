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
#ifndef DBLDEF_FORWARD_2D_NONCONFORMING_QUAD_ELEMENT_0TH_ORDER_EDGE_BASED
#define DBLDEF_FORWARD_2D_NONCONFORMING_QUAD_ELEMENT_0TH_ORDER_EDGE_BASED

#include "Forward2DQuadrilateralElementEdgeBased.h"

// Class of 2D forward calculation by using quadrilateral element 0th order edge based 
class Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased : public Forward2DQuadrilateralElementEdgeBased {

public:

	// Constructer
	explicit Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( const int planeID, const int iPol );

	// Destructer
	~Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased();

	// Calculate EM fields of boundary planes by 2D forward calculcation with 0tht order edge element
	virtual void calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataNonConformingHexaElement* const pMeshData );

private:

	static const int DIRICHLET_BOUNDARY_NONZERO_VALUE = -1;
	static const int DIRICHLET_BOUNDARY_ZERO_VALUE    = -2;

	const static int SLAVE_DOFS = -3;

	const static int m_numGauss = 2;

	const static int m_numIntegralPoints = m_numGauss * m_numGauss;

	double m_integralPointXi[m_numIntegralPoints];

	double m_integralPointEta[m_numIntegralPoints];

	double m_weights[m_numIntegralPoints];

	// Array of reference coord xi values for each node
	double m_xiAtNode[4];

	// Array of reference coord eta values for each node
	double m_etaAtNode[4];

	// Array of reference coord xi values for each edge
	double m_xiAtEdge[4];

	// Array of reference coord eta values for each edge
	double m_etaAtEdge[4];

	// Flag specifing wether map converting master dofs after degeneration and MPC factors from slave dof after degeneration has been made
	bool m_hasMadeMapSlaveDofToMasterDofAndFactors;

	// Array converting global node IDs after degeneration to the ones after constraint
	int* m_IDsAfterDegenerated2AfterConstrained;

	// Total number of equations after degeneration and constraint
	int m_numEquationDegeneratedAndConstrained;

	// Map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
	// Even for master dofs, master dof and factors are inserted (in this case, master dof is equal to index and factor is one)
	std::vector< std::pair<int,double> >* m_slaveDofToMasterDofAndFactors;

	// Defailt constructer
	Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased();

	// Copy constructer
	Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased(const Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased& rhs);

	// Copy assignment operator
	Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased& operator=(const Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased& rhs);

	// Calculate array converting local IDs to global ones
	void calcArrayConvertLocalID2Global( const MeshDataNonConformingHexaElement* const pMeshData );

	// Make map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
	void makeMapSlaveDofToMasterDofAndFactors( const MeshDataNonConformingHexaElement* const pMeshData );

	// Get type of outer edge
	// [note] : You must confirm inputed edge is the outer edge. 
	int getTypeOfOuterEdgeOfBoundaryPlanes( const int edgeIDLocal2D ) const;

	// Get shape functions of the horizontal direction with respect to the reference element coordinate system
	double getShapeFuncH( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat ) const;

	// Get shape functions of the vertical direction with respect to the reference element coordinate system
	double getShapeFuncV( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat) const;

	// Get shape functions rotated with respect to the reference element coordinate system
	double getShapeFuncRotated( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat ) const;

	// Calculate horizontal electric field
	virtual std::complex<double> calcValueElectricFieldHorizontal( const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshDataElement ) const;

	// Calculate vertical electric field
	virtual std::complex<double> calcValueElectricFieldVertical( const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshDataElement ) const;

	// Calculate magnetic field perpendicular to the boundary plane
	virtual std::complex<double> calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshDataElement ) const;

	// Calculate jacobian matrix of the elements on the Z-X plane of boundary
	double calcJacobianMatrixOnZXPlaneOfBoundary( const MeshDataNonConformingHexaElement* const pMeshData, const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const;

	// Calculate jacobian matrix of the elements on the Y-Z plane of boundary
	double calcJacobianMatrixOnYZPlaneOfBoundary( const MeshDataNonConformingHexaElement* const pMeshData, const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const;

	// Calculate jacobian matrix
	double calcJacobianMatrix( const MeshDataNonConformingHexaElement* const pMeshData, const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const;

	// Calculate inverse of jacobian matrix  multiplied by determinant
	void calcInverseOfJacobianMatrix( const Forward2D::Matrix2x2& jacobMat, const double determinant, Forward2D::Matrix2x2& invJacobMat ) const;

	// Add master dof and factor pair to m_slaveDofToMasterDofAndFactors
	void addMasterDofAndFactorPair( const int slaveDof, const int masterDof, const double factor );

	// Set non-zero strucuture of matrix for forward calculation
	void setNonZeroStrucuture( const MeshDataNonConformingHexaElement* const pMeshData );

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	void setNonZeroValues( const double freq, const MeshDataNonConformingHexaElement* const pMeshData );

	// Get flag specifing whether an 2-D element faces slave element
	bool faceSlaveElements( const int iElem, const int iEdge, const MeshDataNonConformingHexaElement* const pMeshData ) const;

	// Get flag specifing whether the inputted element edge is an outer edge
	bool isOuterEdge( const int iElem, const int iEdge, const MeshDataNonConformingHexaElement* const pMeshData ) const;

	// Get neighbor face index from edge index
	int getNeighborFaceIndexFromEdgeIndex( const int iEdge ) const;

};

#endif
