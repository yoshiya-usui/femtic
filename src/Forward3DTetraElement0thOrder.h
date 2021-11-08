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
#ifndef DBLDEF_FORWARD_3D_TETRA_ELEMENT_0TH_ORDER
#define DBLDEF_FORWARD_3D_TETRA_ELEMENT_0TH_ORDER

#include "Forward3D.h"
#include "MeshDataTetraElement.h"
#include "Forward2DTriangleElement0thOrderEdgeBased.h"
#include <set>

// Class of 3D forward calculation by using 0th order tetra element 
class Forward3DTetraElement0thOrder : public Forward3D {

public:

	// Constructer
	Forward3DTetraElement0thOrder();

	// Destructer
	virtual ~Forward3DTetraElement0thOrder();

	//Run 3D forward calculation
	virtual void forwardCalculation( const double freq, const int iPol );

	// Calculate X component of electric field
	virtual std::complex<double> calcValueElectricFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate Y component of electric field
	virtual std::complex<double> calcValueElectricFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate Z component of electric field
	virtual std::complex<double> calcValueElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

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

	// Calculate Z component of rotated electric field
	virtual std::complex<double> calcValueRotatedElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate normal component of rotated electric field
	virtual std::complex<double> calcValueRotatedElectricFieldNormal( const int iElem, const int iFace, const double uCoord, const double vCoord ) const;

	// Calculate X component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate Y component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate Z component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const;

	// Calculate difference of voltage for brick element
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
		const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const;

	// Calculate difference of voltage
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const;

	// Calculate interpolator vector of X component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Y component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of rotated electric field
	virtual void calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

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

	// Calculate interpolator vector of normal component of rotated electric field
	virtual void calcInterpolatorVectorOfRotatedElectricFieldNormal( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of X component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Y component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of Z component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) );

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs );

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs );

	// Set non-zero strucuture of matrix for forward calculation
	virtual void setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix );

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	virtual void setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix );

	//----- DO NOT DELETE FOR FUTURE USE >>>>>
	//// Set non-zero strucuture of matrix for calculating derivatives
	//virtual void setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID, std::set<int>& nonZeroRowsAndCols );

	//// Set non-zero values of matrix and right-hande side vector for calculating derivatives
	//virtual void setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID );
	//----- DO NOT DELETE FOR FUTURE USE <<<<<

	// Calculate vector x of the reciprocity algorithm of Rodi (1976)
	virtual void calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows );

	// Call function inputMeshData of the class MeshData
	virtual void callInputMeshData();

	// Get pointer to the class MeshData
	virtual const MeshData* getPointerToMeshData() const;

	// Get pointer to the class MeshDataTetraElement
	const MeshDataTetraElement* getPointerToMeshDataTetraElement() const;
	
private:

	const static int DIRICHLET_BOUNDARY_NONZERO_VALUE = -1;// This must be the same as the ones of other functions !!

	const static int DIRICHLET_BOUNDARY_ZERO_VALUE    = -2;// This must be the same as the ones of other functions !!

	const static int m_numIntegralPoints = 4;

	double m_uCoord[m_numIntegralPoints];

	double m_vCoord[m_numIntegralPoints];

	double m_wCoord[m_numIntegralPoints];

	double m_weights[m_numIntegralPoints];

	const static double m_eps;

	// Copy constructer
	Forward3DTetraElement0thOrder(const Forward3DTetraElement0thOrder& rhs);

	// Copy assignment operator
	Forward3DTetraElement0thOrder& operator=(const Forward3DTetraElement0thOrder& rhs);

	// Class of 2D forward calculation using triangle elements
	//Forward2DTriangleElement* m_Fwd2DTriangleElement[4][2];
	Forward2DTriangleElement0thOrderEdgeBased* m_Fwd2DTriangleElement[4];

	MeshDataTetraElement m_MeshDataTetraElement;

	bool** m_signInversion;

#ifdef _ANISOTOROPY
	// Set non-zero values of matrix and right-hande side vector for isotropic medium
	void setNonZeroValuesIsotropy( ComplexSparseSquareSymmetricMatrix& matrix );

	// Set non-zero values of matrix and right-hande side vector for anisotropic medium
	void setNonZeroValuesAnisotropy( ComplexSparseSquareSymmetricMatrix& matrix );
#endif

	// Calculate array converting local IDs to global ones
	void calcArrayConvertLocalID2Global();

	// Calculate array converting global IDs to the ones after degeneration
	void calcArrayConvertIDsGlobal2AfterDegenerated();

	// Calculate array converting global edge IDs non-zero electric field values specified to the edges
	void calcArrayConvertIDGlobal2NonZeroValues();

	// Get shape functions of the 1st direction with respect to the reference element coordinate system
	inline double getShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const double wLocal, const int num ) const;

	// Get shape functions of the 2nd direction with respect to the reference element coordinate system
	inline double getShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const double wLocal, const int num ) const;

	// Get shape functions of the 3rd direction with respect to the reference element coordinate system
	inline double getShapeFuncReferenceCoordW( const double uLocal, const double vLocal, const double wLocal, const int num ) const;

	// Get 2D shape functions of the 1st direction with respect to the reference element coordinate system
	inline double get2DShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const int num ) const;

	// Get 2D shape functions of the 2nd direction with respect to the reference element coordinate system
	inline double get2DShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const int num ) const;

	// Get 2D shape functions rotated with respect to the reference element coordinate system
	inline double get2DShapeFuncRotated() const;

	// Get shape functions rotated with respect to the reference element coordinate system ( The 1st direction )
	inline double getShapeFuncRotatedReferenceCoordU( const int num ) const;

	// Get shape functions rotated with respect to the reference element coordinate system ( The 2nd direction )
	inline double getShapeFuncRotatedReferenceCoordV( const int num ) const;

	// Get shape functions rotated with respect to the reference element coordinate system ( The 3rd direction )
	inline double getShapeFuncRotatedReferenceCoordW( const int num ) const;

	// Calculate jacobian matrix of the elements
	void calcJacobianMatrix( const int elemID, Forward3D::Matrix3x3& JacobMat, double& determinant ) const;

	// Calculate 2D jacobian matrix of the faces
	void calc2DJacobianMatrix( const int elemID, const int faceID, Forward3D::Matrix2x2& JacobMat, double& determinant ) const;

	// Calculate the inclinations of the element's face
	void calcInclinationsOfElementFace( const int elemID, const int faceID, double& dLengdX, double& dLengdY ) const;

	// Calculate inverse of jacobian matrix  multiplied by determinant
	void calcInverseOfJacobianMatrix( const Forward3D::Matrix3x3& jacobMat, Forward3D::Matrix3x3& invJacobMat ) const;

	// Calculate flag of sign inversion. This function used by calcVoltageDifference.
	bool calcRatioAndReverseFlag( const int faceID, const int edgeIDLocal2D, const CommonParameters::AreaCoords& startPoint, const CommonParameters::AreaCoords& endPoint, double& ratio ) const;

	// Output results of forward calculation to VTK file
	virtual void outputResultToVTK() const;

	// Output results of forward calculation to binary file
	virtual void outputResultToBinary( const int iFreq, const int iPol ) const;

	// Calculate integrals of which element matrix consisits
	void calcIntegrals( const int elemID, double* eMat, double* fMat ) const;
	
	// Calculate determinant of 3x3 matrix
	double calcDeterminant( const double* rowVec0, const double* rowVec1, const double* rowVec2, const int icol ) const;

};

#endif
