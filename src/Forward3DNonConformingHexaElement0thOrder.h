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
#ifndef DBLDEF_FORWARD_3D_NONCONFORMING_HEXA_ELEMENT_0TH_ORDER
#define DBLDEF_FORWARD_3D_NONCONFORMING_HEXA_ELEMENT_0TH_ORDER

#include "Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased.h"
#include "Forward3D.h"
#include "MeshDataNonConformingHexaElement.h"

// Class of 3D forward calculation by using the 0th order nonconforming hexahedral element 
class Forward3DNonConformingHexaElement0thOrder : public Forward3D {

public:

	// Constructer
	Forward3DNonConformingHexaElement0thOrder();

	// Destructer
	virtual ~Forward3DNonConformingHexaElement0thOrder();

	// Run 3D forward calculation
	virtual void forwardCalculation( const double freq, const int iPol );

	// Calculate X component of electric field
	virtual std::complex<double> calcValueElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate Y component of electric field
	virtual std::complex<double> calcValueElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate Z component of electric field
	virtual std::complex<double> calcValueElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate Z component of rotated electric field
	virtual std::complex<double> calcValueRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate X component of electric field only from the edges on the Earth's surface
	std::complex<double> calcValueRotatedElectricFieldNormal( const int iElem, const double xLocal, const double yLocal ) const;

	// Calculate X component of electric field only from the edges on the Earth's surface
	virtual std::complex<double> calcValueElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const;

	// Calculate Y component of electric field only from the edges on the Earth's surface
	virtual std::complex<double> calcValueElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const;

	// Calculate tangential electric field directed to X from all edges of owner element
	virtual std::complex<double> calcValueElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate tangential electric field directed to Y from all edges of owner element
	virtual std::complex<double> calcValueElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate tangential electric field directed to X
	virtual std::complex<double> calcValueElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord ) const;

	// Calculate tangential electric field directed to Y
	virtual std::complex<double> calcValueElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord ) const;

	// Calculate X component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate Y component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate Z component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const;

	// Calculate interpolator vector of X component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Y component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of rotated electric field
	virtual void calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
	void calcInterpolatorVectorOfRotatedElectricFieldNormal( const int iElem, const double xLocal, const double yLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
	virtual void calcInterpolatorVectorOfElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Y component of electric field only from the edges on the Earth's surface
	virtual void calcInterpolatorVectorOfElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of tangential electric field directed to X from all edges
	virtual void calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of tangential electric field directed to Y from all edges
	virtual void calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of tangential electric field directed to X
	virtual void calcInterpolatorVectorOfElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of tangential electric field directed to Y
	virtual void calcInterpolatorVectorOfElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of X component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Y component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs );

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs );

	// Set non-zero strucuture of matrix for forward calculation
	virtual void setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix );

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	virtual void setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix );

	// Calculate vector x of the reciprocity algorithm of Rodi (1976)
	virtual void calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows );

	// Copy solution vector degenerated
	virtual void copySolutionVectorDegenerated( const int iPol, std::complex<double>* solutionVector ) const;

	// Call function inputMeshData of the class MeshData
	virtual void callInputMeshData();

	// Get pointer to the class MeshData
	virtual const MeshData* getPointerToMeshData() const;

	// Get pointer to the class MeshDataNonConformingHexaElement
	const MeshDataNonConformingHexaElement* getPointerToMeshDataNonConformingHexaElement() const;

	// Calculate difference of voltage for brick element
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
		const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const;

	// Calculate difference of voltage for tetra element
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const;

	// Get total number of equations finally solved
	virtual int getNumOfEquationFinallySolved() const;


private:

	const static int DIRICHLET_BOUNDARY_NONZERO_VALUE = -1;// This must be the same as the ones of other functions !!

	const static int DIRICHLET_BOUNDARY_ZERO_VALUE = -2;// This must be the same as the ones of other functions !!

	const static int SLAVE_ON_OUTER_EDGES = -3;

	const static int SLAVE_ON_INTERIOR_EDGES = -4;

	const static int m_numGauss = 2;

	const static int m_numIntegralPoints = m_numGauss * m_numGauss * m_numGauss;

	double m_integralPointXi[m_numIntegralPoints];

	double m_integralPointEta[m_numIntegralPoints];

	double m_integralPointZeta[m_numIntegralPoints];

	double m_weights[m_numIntegralPoints];

	// Array of reference coord xi values for each node
	double m_xiAtNode[8];

	// Array of reference coord eta values for each node
	double m_etaAtNode[8];

	// Array of reference coord zeta values for each node
	double m_zetaAtNode[8];

	// Array of reference coord xi values for each edge
	double m_xiAtEdge[12];

	// Array of reference coord eta values for each edge
	double m_etaAtEdge[12];

	// Array of reference coord zeta values for each edge
	double m_zetaAtEdge[12];

	const static double m_eps;

	// Flag specifing wether map converting master dofs after degeneration and MPC factors from slave dof after degeneration has been made
	bool m_hasMadeMapSlaveDofToMasterDofAndFactors;

	// Array converting global node IDs after degeneration to the ones after constraint
	int* m_IDsAfterDegenerated2AfterConstrained;

	// Total number of equations after degeneration and constraint
	int m_numOfEquationDegeneratedAndConstrained;

	// Solution vector after degeneration and constraint
	std::complex<double>* m_solutionVectorDegeneratedAndConstrained;

	// Map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
	// Even for master dofs, master dof and factors are inserted (in this case, master dof is equal to index and factor is one)
	std::vector< std::pair<int,double> >* m_slaveDofToMasterDofAndFactors;

	// Vector of MPC constants;
	std::complex<double>* m_vectorMPCConstants;

	// Copy constructer
	Forward3DNonConformingHexaElement0thOrder(const Forward3DNonConformingHexaElement0thOrder& rhs);

	// Copy assignment operator
	Forward3DNonConformingHexaElement0thOrder& operator=(const Forward3DNonConformingHexaElement0thOrder& rhs);

	// Class of 2D forward calculation using quadrilateral elements
	Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased* m_Fwd2DQuadrilateralElement[4];

	MeshDataNonConformingHexaElement m_MeshDataNonConformingHexaElement;

	// Calculate array converting local IDs to global ones
	void calcArrayConvertLocalID2Global();

	// Calculate array converting global IDs to the ones after degeneration
	void calcArrayConvertIDsGlobal2AfterDegenerated();

	// Calculate array converting global edge IDs non-zero electric field values specified to the edges
	void calcArrayConvertIDGlobal2NonZeroValues();

	// Make map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
	double calc2DJacobianMatrixForEarthSurface( const int elemID, const double xi, const double eta, 
		Forward3D::Matrix2x2& JacobMat ) const;

	// Make map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
	void makeMapSlaveDofToMasterDofAndFactors();

	// Calculate MPC constants
	void calcMPCConstants();

	// Add master dof and factor pair to m_slaveDofToMasterDofAndFactors
	bool doesIntegralXCompFirst( const CommonParameters::locationXY& startPoint, const CommonParameters::locationXY& endPoint,
		bool& rotationDirectionPlus, CommonParameters::locationXY& sharedPoint ) const;

	// Add master dof and factor pair to m_slaveDofToMasterDofAndFactors
	void addMasterDofAndFactorPair( const int slaveDof, const int masterDof, const double factor );
	
	// Get shape functions of the x direction with respect to the reference element coordinate system
	double getShapeFuncX( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Get shape functions of the y direction with respect to the reference element coordinate system
	double getShapeFuncY( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Get shape functions of the z direction with respect to the reference element coordinate system
	double getShapeFuncZ( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Get x component of shape function rotated for 0th order edge-based elements
	double getShapeFuncRotatedX( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Get y component of shape function rotated for 0th order edge-based elements
	double getShapeFuncRotatedY( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Get z component of shape function rotated for 0th order edge-based elements
	double getShapeFuncRotatedZ( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const;

	// Calculate jacobian matrix of the elements
	double get2DShapeFuncRotatedForEarthSurface( const double xi, const double eta, const int num, const Forward3D::Matrix2x2& invJacobMat ) const;

	// Calculate jacobian matrix of the elements
	double calcJacobianMatrix( const int elemID,  const double xi, const double eta, const double zeta, Forward3D::Matrix3x3& JacobMat ) const;

	// Calculate inverse of jacobian matrix  multiplied by determinant
	void calcInverseOfJacobianMatrix( const Forward3D::Matrix3x3& jacobMat, const double determinant, Forward3D::Matrix3x3& invJacobMat ) const;
	
	// Output results of forward calculation to VTK file
	virtual void outputResultToVTK() const;

	// Output results of forward calculation to binary file
	virtual void outputResultToBinary( const int iFreq, const int iPol ) const;
	
	// Add values to right-hand sides matrix consisting of interpolator vectors by taking into account MPC
	void addValuesToRhsVectorsByConsideringMPC( const int irow, const int irhs, const std::complex<double>& val );

	//// Make pair of master dof and factor for a slave dof
	//void makePairOfMasterDofAndFactorForASlaveDof( const int slaveDofAfterDegenerated, std::vector< std::pair<int,double> >& masterDofAndFactorsAfterDegenerated ) const;

	//void makePairOfMasterDofAndFactorForASlaveDofAux( const int slaveDofBeforeDegenerated, std::vector< std::pair<int,double> >& masterDofAndFactorsBeforeDegenerated ) const;

	//// Get dof before degeneration from the after degeneration
	//int getDofBeforeDegenerationFromDofAfterDegeneration( const int dofAfterDegeneration ) const;

};

#endif
