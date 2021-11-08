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
#ifndef DBLDEF_FORWARD_3D
#define DBLDEF_FORWARD_3D

//#include "ComplexSparseSquareMatrix.h"
#include "ComplexSparseSquareSymmetricMatrix.h"
#include "CommonParameters.h"
#include <map>
#include <set>
#include <iostream>
#include <stdlib.h>
#include "MeshData.h"

// Class of 3D forward analysis
class Forward3D{

public:
	// Constructer
	Forward3D();

	// Destructer
	virtual ~Forward3D();

	//Run 3D forward calculation
	virtual void forwardCalculation( const double freq, const int iPol ) = 0;

	// Calculate X component of electric field
	virtual std::complex<double> calcValueElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate Y component of electric field
	virtual std::complex<double> calcValueElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate Z component of electric field
	virtual std::complex<double> calcValueElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate Z component of rotated electric field
	virtual std::complex<double> calcValueRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate X component of electric field only from the edges on the Earth's surface
	virtual std::complex<double> calcValueElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const = 0;

	// Calculate Y component of electric field only from the edges on the Earth's surface
	virtual std::complex<double> calcValueElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const = 0;

	// Calculate tangential electric field directed to X from all edges of owner element
	virtual std::complex<double> calcValueElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate tangential electric field directed to Y from all edges of owner element
	virtual std::complex<double> calcValueElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate tangential electric field directed to X
	virtual std::complex<double> calcValueElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord ) const = 0;

	// Calculate tangential electric field directed to Y
	virtual std::complex<double> calcValueElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord ) const = 0;

	// Calculate X component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate Y component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate Z component of magnetic field
	virtual std::complex<double> calcValueMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const = 0;

	// Calculate interpolator vector of X component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Y component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Z component of electric field
	virtual void calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Z component of rotated electric field
	virtual void calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
	virtual void calcInterpolatorVectorOfElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Y component of electric field only from the edges on the Earth's surface
	virtual void calcInterpolatorVectorOfElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of tangential electric field directed to X from all edges
	virtual void calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of tangential electric field directed to Y from all edges
	virtual void calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of tangential electric field directed to X
	virtual void calcInterpolatorVectorOfElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of tangential electric field directed to Y
	virtual void calcInterpolatorVectorOfElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of X component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Y component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of Z component of magnetic field
	virtual void calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor = std::complex<double>(1.0,0.0) ) = 0;

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs ) = 0;

	// Calculate interpolator vector of difference of voltage
	virtual void calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs ) = 0;

	// Set non-zero strucuture of matrix for forward calculation
	virtual void setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix ) = 0;

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	virtual void setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ) = 0;

	//----- DO NOT DELETE FOR FUTURE USE >>>>>
	//// Set non-zero strucuture of matrix for calculating derivatives
	//virtual void setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID, std::set<int>& nonZeroRowsAndCols ) = 0;

	//// Set non-zero values of matrix and right-hande side vector for calculating derivatives
	//virtual void setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID ) = 0;
	//----- DO NOT DELETE FOR FUTURE USE <<<<<

	// Calculate vector x of the reciprocity algorithm of Rodi (1976)
	virtual void calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows ) = 0;

	// Copy solution vector degenerated
	virtual void copySolutionVectorDegenerated( const int iPol, std::complex<double>* solutionVector ) const;

	// Call function inputMeshData of the class MeshData
	virtual void callInputMeshData() = 0;

	// Get pointer to the class MeshData
	virtual const MeshData* getPointerToMeshData() const = 0;

	// Get polarization at present for which forward analysis is executed 
	int getPolarizationCurrent() const;

	// Get frequency at present for which forward analysis is executed 
	double getFrequencyCurrent() const;

	// Get order of finite element
	int getOrderOfFiniteElement() const;

	// Get total number of equations after degeneration
	int getNumOfEquationDegenerated() const;

	// Get total number of equations finally solved
	virtual int getNumOfEquationFinallySolved() const;

	// Release memory of total matrix and sparse solver
	void releaseMemoryOfMatrixAndSolver();

	// Initialize right-hand side vectors
	void initializeRhsVectors( const int nrhs );

	// Perform solve phase for right-hand sides consisting of interpolator vectors
	void solvePhaseForRhsConsistingInterpolatorVectors( const int numInterpolatorVectors, std::complex<double>* solutionForInterpolatorVectors );

	// Calculate derivative of EM field
	void calculateDerivativesOfEMField( const int numInterpolatorVectors, const std::complex<double>* const solutionForInterpolatorVectors, std::complex<double>* const derivatives );

	// Allocate memory for derivatives of interpolator vectors
	void allcateMemoryForDerivativeOfInterpolatorVectors( const int numInterpolatorVectors );

	// Calculate difference of voltage for brick element
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
		const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const = 0;

	// Calculate difference of voltage for tetra element
	virtual std::complex<double> calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const = 0;

protected:

	struct Matrix2x2{
		double mat11;
		double mat12;
		double mat21;
		double mat22;
	};

	struct Matrix3x3{
		double mat11;
		double mat12;
		double mat13;
		double mat21;
		double mat22;
		double mat23;
		double mat31;
		double mat32;
		double mat33;
	};

	// Total number of the equation
	int m_numOfEquation;

	// Total number of equations after degeneration
	int m_numOfEquationDegenerated;

	// Array converting local edge IDs to global ones
	int** m_IDsLocal2Global;

	// Whether array converting local edge IDs to global ones has already been set or not
	bool m_hasSetIDsLocal2Global;

	// Array converting global node IDs to the ones after degeneration
	int* m_IDsGlobal2AfterDegenerated[2];

	// Array converting global node IDs non-zero electric field values specified to the nodes
	std::map< int, std::complex<double> > m_globalID2NonZeroValues;

	// Whether array converting global node IDs to the ones after degeneration has already been set or not
	bool m_hasIDsGlobal2AfterDegenerated[2];

	// Array converting global node IDs of slaves to the ones of its master
	std::map<int, int> m_globalIDSlave2Master[2];

	// Array of the matrix of 3D anaysis
	ComplexSparseSquareSymmetricMatrix m_matrix3DAnalysis;

	// Whether matrix structure has already been set or not
	bool m_hasMatrixStructureSetAndAnalyzed;

	// Solution vector of 3D analysis
	std::complex<double>* m_solution;
	
	// Set polarization at present for which forward analysis is executed 
	void setPolarizationCurrent( const int iPol );

	// Set frequency at present for which forward analysis is executed 
	void setFrequencyCurrent( const double freq );

	// Set order of finite element
	void setOrderOfFiniteElement( const int order );

	// Add values to right-hand sides matrix consisting of interpolator vectors
	void addValuesToRhsVectors( const int irow, const int irhs, const std::complex<double>& val );

	// Output results of forward calculation to VTK file
	virtual void outputResultToVTK() const = 0;

	// Initialize sparse solver
	void initializeSparseSolver();

private:
	// Copy constructer
	Forward3D(const Forward3D& rhs){
		std::cerr << "Error : Copy constructer of the class Forward3D is not implemented." << std::endl;
		exit(1);
	}

	// Copy assignment operator
	Forward3D& operator=(const Forward3D& rhs){
		std::cerr << "Error : Copy constructer of the class Forward3D is not implemented." << std::endl;
		exit(1);
	}

	// Polarization at present for which forward analysis is executed 
	int m_polarizationCurrent;

	// Frequency at present for which forward analysis is executed 
	double m_frequencyCurrent;

	// Order of finite element
	int m_orderOfFiniteElement;

};

#endif
