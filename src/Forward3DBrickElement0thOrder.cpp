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
#include "Forward3DBrickElement0thOrder.h"
#include "MeshDataBrickElement.h"
#include "ResistivityBlock.h"
#include "Util.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "Forward2DSquareElement0thOrderEdgeBased.h"
#include "Forward2DSquareElement1stOrderNodeBased.h"
#include "ObservedData.h"
#include <stdio.h>
#include <string.h>

#include <map>
#include <new>
#include <vector>
#include <assert.h>

#ifdef _USE_OMP
#include <omp.h>
#endif

Forward3DBrickElement0thOrder::Forward3DBrickElement0thOrder():
	Forward3D()
{

	for( int i = 0; i < 4; ++i ){
		m_Fwd2DSquareElement[i] = NULL;
	}
	
	m_Fwd2DSquareElement[MeshData::ZXMinus] = new Forward2DSquareElement0thOrderEdgeBased( MeshData::ZXMinus, CommonParameters::EX_POLARIZATION );
	m_Fwd2DSquareElement[MeshData::ZXPlus ] = new Forward2DSquareElement0thOrderEdgeBased( MeshData::ZXPlus,  CommonParameters::EX_POLARIZATION );

	m_Fwd2DSquareElement[MeshData::YZMinus] = new Forward2DSquareElement0thOrderEdgeBased( MeshData::YZMinus, CommonParameters::EY_POLARIZATION );
	m_Fwd2DSquareElement[MeshData::YZPlus ] = new Forward2DSquareElement0thOrderEdgeBased( MeshData::YZPlus,  CommonParameters::EY_POLARIZATION );

	//---------------------------------------------------------------------------
	//--- Calculate integral points and weights of two point Gauss quadrature ---
	//---------------------------------------------------------------------------
	int ip(0);
	for( int i = 0; i < m_numGauss; ++i ){
		for( int j = 0; j < m_numGauss; ++j ){
			for( int k = 0; k < m_numGauss; ++k ){
				m_xLocal[ip] = CommonParameters::abscissas2Point[i];
				m_yLocal[ip] = CommonParameters::abscissas2Point[j];
				m_zLocal[ip] = CommonParameters::abscissas2Point[k];
				m_weights3D[ip] = CommonParameters::weights2Point[i] * CommonParameters::weights2Point[j] * CommonParameters::weights2Point[k];
				++ip;
			}
		}
	}

}

//Destructer
Forward3DBrickElement0thOrder::~Forward3DBrickElement0thOrder(){

	for( int i = 0; i < 4; ++i ){
		if( m_Fwd2DSquareElement[i] != NULL ){
			delete m_Fwd2DSquareElement[i];
		}
	}

}

//Run 3D forward calculation by using brick element
void Forward3DBrickElement0thOrder::forwardCalculation( const double freq, const int iPol ){

	setFrequencyCurrent(freq);
	setPolarizationCurrent(iPol);
	setOrderOfFiniteElement( 0 );	

	if( m_matrix3DAnalysis.getNumRightHandSideVectors() != 1 ){
		m_matrix3DAnalysis.reallocateMemoryForRightHandSideVectors(1);
	}

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	if( pAnalysisControl->getBoundaryConditionBottom() != AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		OutputFiles::m_logFile << "Error : When 0th order 3D edge-based element is used, electric field of the bottom must be zero !" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Start 3D forward calculation with 0th order 3D edge-based element." << pAnalysisControl->outputElapsedTime() << std::endl;

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

#ifdef _DEBUG_WRITE
	std::cout << "nElem " << nElem << std::endl;// For debug
#endif

	//----------------------------------------------------------------
	//--- Calculate array converting local edge IDs to global ones ---
	//----------------------------------------------------------------
	if( !m_hasSetIDsLocal2Global ){
		calcArrayConvertLocalID2Global();
	}

	//---------------------------------------------------------------------------------
	//--- Calculate array converting global node IDs to the ones after degeneration ---
	//---------------------------------------------------------------------------------
	if( !m_hasIDsGlobal2AfterDegenerated[iPol] ){
		calcArrayConvertIDsGlobal2AfterDegenerated();
	}

	OutputFiles::m_logFile << "# Number of equation = " << m_numOfEquation
		<< ", Number of equation after degeneration = " << m_numOfEquationDegenerated << std::endl;

	//--------------------------------------------------------------------------------------------------------
	//--- Calculate array converting global node IDs non-zero electric field values specified to the nodes ---
	//--------------------------------------------------------------------------------------------------------
	calcArrayConvertIDGlobal2NonZeroValues();

	//-----------------------------------------------------------------------
	//--- Set matrix structure and analyze the structure by matrix solver ---
	//-----------------------------------------------------------------------
	if( !m_hasSetIDsLocal2Global ){
		OutputFiles::m_logFile << "Error : m_hasSetIDsLocal2Global is NULL !!" << std::endl;
		exit(1);
	}
	if( !m_hasIDsGlobal2AfterDegenerated[iPol] ){
		OutputFiles::m_logFile << "Error : m_hasIDsGlobal2AfterDegenerated[iPol] is NULL !!" << std::endl;
		exit(1);
	}

	if( !m_hasMatrixStructureSetAndAnalyzed ){ 
		// If matrix structure has not been set yet, set matrix structure.
		OutputFiles::m_logFile << "# Set matrix structure. " << pAnalysisControl->outputElapsedTime() << std::endl;

		initializeSparseSolver();

		m_matrix3DAnalysis.setDegreeOfEquation( m_numOfEquationDegenerated );

		//setNonZeroStrucuture( m_matrix3DAnalysis, nElem, false, -1 );
		setNonZeroStrucuture( m_matrix3DAnalysis );

		//--- Convert the matrix from the triplet format to the CRS format
		if( m_matrix3DAnalysis.hasConvertedToCRSFormat() == false ){
			m_matrix3DAnalysis.convertToCRSFormat();
		}

		//--- Anaysis phase of matrix solver
		OutputFiles::m_logFile << "# Anaysis phase of matrix solver for 3D forward calculation."
			<<  " Polarization : " << iPol << ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
		m_matrix3DAnalysis.analysisPhaseMatrixSolver();

		//---- Output memory required to caluculate the 2D analysis
		m_matrix3DAnalysis.writeMemoryRequiredByMatrixSolver(); 
		
		m_hasMatrixStructureSetAndAnalyzed = true;

	}

	//-------------------------------------------------------
	//--- Set values of matrix and right hand side vector ---
	//-------------------------------------------------------
	OutputFiles::m_logFile << "# Set values of matrix and right hand side vector. " << pAnalysisControl->outputElapsedTime() << std::endl;

	m_matrix3DAnalysis.zeroClearNonZeroValues();// Zero clear matrix values
	m_matrix3DAnalysis.zeroClearRightHandSideVector();// Zero clear right hand side vector

	//setNonZeroValues( m_matrix3DAnalysis, nElem, false, -1 );
	setNonZeroValues( m_matrix3DAnalysis );

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	m_matrix3DAnalysis.debugWriteMatrix();
	m_matrix3DAnalysis.debugWriteRightHandSide();
#endif
	//----- debug <<<<<

	//-----------------------------------------------------
	//--- Numrical factorization phase of matrix solver ---
	//-----------------------------------------------------
	OutputFiles::m_logFile << "# Numerical factorization phase of matrix solver for 3D forward calculation."
			<<  " Polarization : " << iPol << "." 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix3DAnalysis.factorizationPhaseMatrixSolver();//Numrical factorization phase of matrix solver

	//------------------------------------
	//--- Solve phase of matrix solver ---
	//------------------------------------
	std::complex<double>* solutionDegenerated = new std::complex<double>[ m_numOfEquationDegenerated ];
	OutputFiles::m_logFile << "# Solve phase of matrix solver for 3D forward calculation."
			<<  " Polarization : " << iPol << "." 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix3DAnalysis.solvePhaseMatrixSolver( solutionDegenerated );//Solve phase of matrix solver

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		std::cout << "i solutionDegenerated : " << i << " " << solutionDegenerated[i] << std::endl;
	}
#endif
	//----- debug <<<<<

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[ m_numOfEquation ];
				
	bool* alreadyFound = new bool[ m_numOfEquation ];
	for( int i = 0; i < m_numOfEquation; ++i ){
		alreadyFound[i] = false;
	}

	for( int iElem = 0; iElem < nElem; ++iElem ){

		for( int iEdge = 0; iEdge < 12; ++iEdge ){

			const int iNum = m_IDsLocal2Global[iElem][iEdge];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ];

				if( iNumDegenerated == DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					m_solution[iNum] = m_globalID2NonZeroValues[ iNum ];
				}else{
					m_solution[iNum] = solutionDegenerated[ iNumDegenerated ]; 
				}

				alreadyFound[ iNum ] = true;
			}
			
		}// iEdge
	}// iElem
	
	delete[] solutionDegenerated;
	delete[] alreadyFound;

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquation; ++i ){
		std::cout << "i m_solution : " << i << " " << m_solution[i] << std::endl;
	}
#endif
	//----- debug <<<<<
	
	// Output EM field vectors
	if( pAnalysisControl->writeBinaryFormat() ){// BINARY
		outputResultToBinary( (ObservedData::getInstance())->getFreqIDs( freq ), iPol );
	}
	else{// ASCII
		outputResultToVTK();
	}

	//----------------------
	//--- Release memory ---
	//----------------------
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] IDsLocal2Global[iElem];
	//}
	//delete [] IDsLocal2Global;

	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] IDsLocal2GlobalDegenerated[iElem];
	//}
	//delete [] IDsLocal2GlobalDegenerated;
	//m_globalID2NonZeroValues.clear();
	//delete [] IDsGlobal2AfterDegenerated;
}

// Calculate X component electric field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	return calcValueElectricFieldXDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else if( order == 1 ){
	//	return calcValueElectricFieldXDirection1stOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//return std::complex<double>(0.0,0.0);

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueElectricFieldXDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueElectricFieldXDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncX( xLocal, yLocal, zLocal, i );
	}
	return val;

}

// Calculate Y component electric field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	return calcValueElectricFieldYDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else if( order == 1 ){
	//	return calcValueElectricFieldYDirection1stOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//return std::complex<double>(0.0,0.0);

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueElectricFieldYDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueElectricFieldYDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncY( xLocal, yLocal, zLocal, i );
	}
	return val;

}

// Calculate Z component electric field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	return calcValueElectricFieldZDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else if( order == 1 ){
	//	return calcValueElectricFieldZDirection1stOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//return std::complex<double>(0.0,0.0);

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueElectricFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueElectricFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncZ( xLocal, yLocal, zLocal, i );
	}
	return val;

}

// Calculate Z component of rotated electric field
std::complex<double> Forward3DBrickElement0thOrder::calcValueRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueMagneticFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueMagneticFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotatedZ( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i );
	}

	return val;

}

// Calculate X component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate Y component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate tangential electric field directed to X from all edges of owner element
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate tangential electric field directed to Y from all edges of owner element
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate tangential electric field directed to X
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate tangential electric field directed to Y
std::complex<double> Forward3DBrickElement0thOrder::calcValueElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate X component magnetic field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	return calcValueMagneticFieldXDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else if( order == 1 ){
	//	return calcValueMagneticFieldXDirection1stOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//return std::complex<double>(0.0,0.0);

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueMagneticFieldXDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueMagneticFieldXDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotatedX( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i );
	}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Y component magnetic field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	return calcValueMagneticFieldYDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else if( order == 1 ){
	//	return calcValueMagneticFieldYDirection1stOrderEdgeBased( iElem, xLocal, yLocal, zLocal );
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//return std::complex<double>(0.0,0.0);

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueMagneticFieldYDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueMagneticFieldYDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotatedY( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i );
	}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Z component magnetic field values for 0th order edge-based elements
std::complex<double> Forward3DBrickElement0thOrder::calcValueMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL in calcValueMagneticFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL in calcValueMagneticFieldZDirection0thOrderEdgeBased." << std::endl;
	//	exit(1);
	//}

	////const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	//const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	//const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	//const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	//std::complex<double> val(0.0, 0.0);
	//for( int i = 0; i < 12; ++i ){
	//	val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotatedZ( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i );
	//}

	//const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	//const double factor = omega * CommonParameters::mu;

	//val /= std::complex<double>(0.0, factor);

	//return val;

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	return calcValueRotatedElectricFieldZDirection( iElem, xLocal, yLocal, zLocal ) / std::complex<double>(0.0, omega * CommonParameters::mu);

}

// Calculate difference of voltage
std::complex<double> Forward3DBrickElement0thOrder::calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
	const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const{

	std::complex<double> voltageDifference = std::complex<double>(0.0, 0.0);

	for( int ielem = 0; ielem < nElem; ++ielem ){

		//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "elememtsIncludingDipole[ielem] : " << elememtsIncludingDipole[ielem] << std::endl;
#endif
//----- debug <<<<<

		bool rotationDirectionPlus = true;
		const bool integralXCompFirst = doesIntegralXCompFirst( localCoordinateValuesStartPoint[ielem], localCoordinateValuesEndPoint[ielem], rotationDirectionPlus );

		const double lengthX = m_MeshDataBrickElement.getEdgeLengthX( elememtsIncludingDipole[ielem] );
		const double lengthY = m_MeshDataBrickElement.getEdgeLengthY( elememtsIncludingDipole[ielem] );

		const double XCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].X - localCoordinateValuesStartPoint[ielem].X ) * lengthX * 0.5;
		const double YCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].Y - localCoordinateValuesStartPoint[ielem].Y ) * lengthY * 0.5;

//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "integralXCompFirst : " << integralXCompFirst << std::endl;
		std::cout << "lengthX, lengthY : " << lengthX << " " << lengthY << std::endl;
		std::cout << "XCoordDifferenceOfSegment, YCoordDifferenceOfSegment : " << XCoordDifferenceOfSegment << " " << YCoordDifferenceOfSegment << std::endl;
#endif
//----- debug <<<<<

		if( integralXCompFirst ){
			voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 );// Since electric field is constant on edges
			voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0 );// Since electric field is constant on edges
		}else{
			voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 );// Since electric field is constant on edges
			voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0 );// Since electric field is constant on edges
		}

//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "voltageDifference : " << voltageDifference << std::endl;
#endif
//----- debug <<<<<

		double areaWithSign = 0.5 * fabs( XCoordDifferenceOfSegment * YCoordDifferenceOfSegment );
		if( !rotationDirectionPlus ){
			areaWithSign *= -1.0;
		}

		voltageDifference += areaWithSign * calcValueRotatedElectricFieldZDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 );// Since magnetic field is constant in the area

//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "areaWithSign calcValueRotatedElectricFieldZDirection : " << areaWithSign << " " << calcValueRotatedElectricFieldZDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 ) << std::endl;
		std::cout << "voltageDifference  total : " << voltageDifference << std::endl;
#endif
//----- debug <<<<<

	}

	return voltageDifference;

//	const double eps = 1.0e-12;
//
//	std::complex<double> voltageDifference = std::complex<double>(0.0, 0.0);
//
//	for( int ielem = 0; ielem < nElem; ++ielem ){
//
//		const double lengthX = m_MeshDataBrickElement.getEdgeLengthX( elememtsIncludingDipole[ielem] );
//		const double lengthY = m_MeshDataBrickElement.getEdgeLengthY( elememtsIncludingDipole[ielem] );
//
//		const double XCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].X - localCoordinateValuesStartPoint[ielem].X ) * lengthX * 0.5;
//		const double YCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].Y - localCoordinateValuesStartPoint[ielem].Y ) * lengthY * 0.5;
//
//#ifdef _DEBUG_WRITE
//		std::cout << "lengthX : " << lengthX << std::endl;
//		std::cout << "lengthY : " << lengthY << std::endl;
//		std::cout << "XCoordDifferenceOfSegment : " << XCoordDifferenceOfSegment << std::endl;
//		std::cout << "YCoordDifferenceOfSegment : " << YCoordDifferenceOfSegment << std::endl;
//#endif
//
//		if( fabs( localCoordinateValuesStartPoint[ielem].X - 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].X - 1.0 ) < eps ){// Segment locate on edge #6
//#ifdef _DEBUG_WRITE
//			std::cout << "A" << std::endl;
//#endif
//			voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][6] ];
//			continue;// Go Next segment
//		}//-----------------------------------------------------------------------------------------------------------------------------------
//		else if( fabs( localCoordinateValuesStartPoint[ielem].X + 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].X + 1.0 ) < eps ){// Segment locate on edge #4
//#ifdef _DEBUG_WRITE
//			std::cout << "B" << std::endl;
//#endif
//			voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][4] ];
//			continue;// Go Next segment
//		}//-----------------------------------------------------------------------------------------------------------------------------------
//		else if( fabs( localCoordinateValuesStartPoint[ielem].Y - 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].Y - 1.0 ) < eps ){// Segment locate on edge #1
//#ifdef _DEBUG_WRITE
//			std::cout << "C" << std::endl;
//#endif
//			voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][1] ];
//			continue;// Go Next segment
//		}//-----------------------------------------------------------------------------------------------------------------------------------
//		else if( fabs( localCoordinateValuesStartPoint[ielem].Y + 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].Y + 1.0 ) < eps ){// Segment locate on edge #0
//#ifdef _DEBUG_WRITE
//			std::cout << "D" << std::endl;
//#endif
//			voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][0] ];
//			continue;// Go Next segment
//		}//-----------------------------------------------------------------------------------------------------------------------------------
//		else{
//
//			double areaWithSign(0.0);
//
//			if( fabs( localCoordinateValuesStartPoint[ielem].Y + 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].Y - 1.0 ) < eps ){// Trapezoidal area including edge #6
//				
//#ifdef _DEBUG_WRITE
//				std::cout << "E" << std::endl;
//#endif
//
//				const double diff1 = ( 1.0 - localCoordinateValuesStartPoint[ielem].X ) * lengthX * 0.5;
//				const double diff2 = ( localCoordinateValuesEndPoint[ielem].X - 1.0   ) * lengthX * 0.5;
//				voltageDifference -= std::complex<double>( diff1 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][0] ];
//				voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][6] ];
//				voltageDifference -= std::complex<double>( diff2 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][1] ];
//
//				areaWithSign = 0.5 * fabs( ( diff1 + diff2 ) * YCoordDifferenceOfSegment );
//
//#ifdef _DEBUG_WRITE
//				std::cout << "diff1 diff2 : " << diff1 << " " << diff2 << std::endl;
//				std::cout << "areaWithSign : " << areaWithSign << std::endl;
//#endif
//
//			}//-----------------------------------------------------------------------------------------------------------------------------------
//			else if( fabs( localCoordinateValuesStartPoint[ielem].Y - 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].Y + 1.0 ) < eps ){// Trapezoidal area including edge #4
//
//#ifdef _DEBUG_WRITE
//				std::cout << "F" << std::endl;
//#endif
//
//				const double diff1 = ( -1.0 - localCoordinateValuesStartPoint[ielem].X ) * lengthX * 0.5;
//				const double diff2 = ( localCoordinateValuesEndPoint[ielem].X + 1.0    ) * lengthX * 0.5;
//				voltageDifference -= std::complex<double>( diff1 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][1] ];
//				voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][4] ];
//				voltageDifference -= std::complex<double>( diff2 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][0] ];
//
//				areaWithSign = 0.5 * fabs( ( diff1 + diff2 ) * YCoordDifferenceOfSegment );
//
//#ifdef _DEBUG_WRITE
//				std::cout << "diff1 diff2 : " << diff1 << " " << diff2 << std::endl;
//				std::cout << "areaWithSign : " << areaWithSign << std::endl;
//#endif
//
//			}//-----------------------------------------------------------------------------------------------------------------------------------
//			else if( fabs( localCoordinateValuesStartPoint[ielem].X - 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].X + 1.0 ) < eps ){// Trapezoidal area including edge #1
//
//#ifdef _DEBUG_WRITE
//				std::cout << "G" << std::endl;
//#endif
//
//				const double diff1 = ( 1.0 - localCoordinateValuesStartPoint[ielem].Y ) * lengthY * 0.5;
//				const double diff2 = ( localCoordinateValuesEndPoint[ielem].Y - 1.0   ) * lengthY * 0.5;
//				voltageDifference -= std::complex<double>( diff1 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][6] ];
//				voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][1] ];
//				voltageDifference -= std::complex<double>( diff2 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][4] ];
//
//				areaWithSign = 0.5 * fabs( ( diff1 + diff2 ) * XCoordDifferenceOfSegment );
//
//#ifdef _DEBUG_WRITE
//				std::cout << "diff1 diff2 : " << diff1 << " " << diff2 << std::endl;
//				std::cout << "areaWithSign : " << areaWithSign << std::endl;
//#endif
//
//			}//-----------------------------------------------------------------------------------------------------------------------------------
//			else if( fabs( localCoordinateValuesStartPoint[ielem].X + 1.0 ) < eps && fabs( localCoordinateValuesEndPoint[ielem].X - 1.0 ) < eps ){// Trapezoidal area including edge #0
//
//#ifdef _DEBUG_WRITE
//				std::cout << "H" << std::endl;
//#endif
//
//				const double diff1 = ( -1.0 - localCoordinateValuesStartPoint[ielem].Y ) * lengthY * 0.5;
//				const double diff2 = ( localCoordinateValuesEndPoint[ielem].Y + 1.0    ) * lengthY * 0.5;
//				voltageDifference -= std::complex<double>( diff1 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][4] ];
//				voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][0] ];
//				voltageDifference -= std::complex<double>( diff2 , 0.0 ) * m_solution[ m_IDsLocal2Global[ielem][6] ];
//
//				areaWithSign = 0.5 * fabs( ( diff1 + diff2 ) * XCoordDifferenceOfSegment );
//
//#ifdef _DEBUG_WRITE
//				std::cout << "diff1 diff2 : " << diff1 << " " << diff2 << std::endl;
//				std::cout << "areaWithSign : " << areaWithSign << std::endl;
//#endif
//
//			}//-----------------------------------------------------------------------------------------------------------------------------------
//			else{// Triangle area
//
//				if( XCoordDifferenceOfSegment * YCoordDifferenceOfSegment > 0 ){// Rotation direction is +Z
//#ifdef _DEBUG_WRITE
//					std::cout << "I-1" << std::endl;
//#endif
//					voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 );// Since electric field is constant on edges
//					voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0 );// Since electric field is constant on edges
//				}else{// Rotation direction is -Z
//#ifdef _DEBUG_WRITE
//					std::cout << "I-2" << std::endl;
//#endif
//					voltageDifference -= std::complex<double>( YCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0 );// Since electric field is constant on edges
//					voltageDifference -= std::complex<double>( XCoordDifferenceOfSegment , 0.0 ) * calcValueElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0 );// Since electric field is constant on edges
//				}
//
//				areaWithSign = 0.5 * XCoordDifferenceOfSegment * YCoordDifferenceOfSegment;
//
//#ifdef _DEBUG_WRITE
//				std::cout << "areaWithSign : " << areaWithSign << std::endl;
//#endif
//
//			}
//
//			voltageDifference += areaWithSign * calcValueRotatedElectricFieldZDirection( elememtsIncludingDipole[ielem], 0.0, 0.0, -1.0 );// Since magnetic field is constant in the area
//	
//			continue;// Go Next segment
//
//		}
//		
//	}
//
//	return voltageDifference;

}

// Calculate difference of voltage for tetra element
std::complex<double> Forward3DBrickElement0thOrder::calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
	const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const{

	OutputFiles::m_logFile << "Error : calcVoltageDifference for tetrahedral element can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of X component of electric field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	calcInterpolatorVectorOfElectricFieldXDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal, irhs );
	//}else if( order == 1 ){
	//	OutputFiles::m_logFile << "Error : Interpolator vectors of X component of electric field cannot be calculated when 1st order element is used." << std::endl;
	//	exit(1);
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncX( xLocal, yLocal, zLocal, i ) , 0.0 ) * factor );
	}

}

// Calculate interpolator vector of Y component of electric field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	calcInterpolatorVectorOfElectricFieldYDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal, irhs );
	//}else if( order == 1 ){
	//	OutputFiles::m_logFile << "Error : Interpolator vectors of Y component of electric field cannot be calculated when 1st order element is used." << std::endl;
	//	exit(1);
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncY( xLocal, yLocal, zLocal, i ), 0.0 ) * factor );
	}

}

// Calculate interpolator vector of Z component of electric field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	calcInterpolatorVectorOfElectricFieldZDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal, irhs );
	//}else if( order == 1 ){
	//	OutputFiles::m_logFile << "Error : Interpolator vectors of Z component of electric field cannot be calculated when 1st order element is used." << std::endl;
	//	exit(1);
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncZ( xLocal, yLocal, zLocal, i ), 0.0 ) * factor );
	}

}

// Calculate interpolator vector of Z component of rotated electric field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncRotatedZ( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i ), 0.0 ) * factor );
	}

}

// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, 
	const int irhs, const std::complex<double>& factor ){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of Y component of electric field only from the edges on the Earth's surface
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){
		
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of tangential electric field directed to X from all edges
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, 
	const int irhs, const std::complex<double>& factor ){
	
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of tangential electric field directed to Y from all edges
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal,
	const int irhs, const std::complex<double>& factor ){
	
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of tangential electric field directed to X
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of tangential electric field directed to Y
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " can't be used for brick element ." << std::endl;
	exit(1);

}

// Calculate interpolator vector of X component of magnetic field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	calcInterpolatorVectorOfMagneticFieldXDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal, irhs );
	//}else if( order == 1 ){
	//	OutputFiles::m_logFile << "Error : Interpolator vectors of Z component of magnetic field cannot be calculated when 1st order element is used." << std::endl;
	//	exit(1);
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncRotatedX( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i ), 0.0 ) * constant );
	}

}

// Calculate interpolator vector of Y component of magnetic field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//const int order = getOrderOfFiniteElement();
	//
	//if( order == 0 ){
	//	calcInterpolatorVectorOfMagneticFieldYDirection0thOrderEdgeBased( iElem, xLocal, yLocal, zLocal, irhs );
	//}else if( order == 1 ){
	//	OutputFiles::m_logFile << "Error : Interpolator vectors of Z component of magnetic field cannot be calculated when 1st order element is used." << std::endl;
	//	exit(1);
	//}else{
	//	OutputFiles::m_logFile << "Error : Order of finite element must be 1 or 2 !! : order = " << order << "." << std::endl;
	//	exit(1);
	//}

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncRotatedY( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i ), 0.0 ) * constant );
	}

}

// Calculate interpolator vector of Z component of magnetic field
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	//const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	////const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	//const double lengX = m_MeshDataBrickElement.getEdgeLengthX( iElem );
	//const double lengY = m_MeshDataBrickElement.getEdgeLengthY( iElem );
	//const double lengZ = m_MeshDataBrickElement.getEdgeLengthZ( iElem );

	//const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	//const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	//for( int i = 0; i < 12; ++i ){
	//	const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
	//	if( irow < 0 ){
	//		continue;
	//	}
	//	addValuesToRhsVectors( irow, irhs, std::complex<double>( getShapeFuncRotatedZ( xLocal, yLocal, zLocal, lengX, lengY, lengZ, i ), 0.0 ) * constant );
	//}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	calcInterpolatorVectorOfRotatedElectricFieldZDirection( iElem, xLocal, yLocal, zLocal, irhs, constant );

}

// Calculate interpolator vector of difference of voltage
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
	const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	for( int ielem = 0; ielem < nElem; ++ielem ){

		bool rotationDirectionPlus = true;
		const bool integralXCompFirst = doesIntegralXCompFirst( localCoordinateValuesStartPoint[ielem], localCoordinateValuesEndPoint[ielem], rotationDirectionPlus );

		const double lengthX = m_MeshDataBrickElement.getEdgeLengthX( elememtsIncludingDipole[ielem] );
		const double lengthY = m_MeshDataBrickElement.getEdgeLengthY( elememtsIncludingDipole[ielem] );
		const double lengthZ = m_MeshDataBrickElement.getEdgeLengthZ( elememtsIncludingDipole[ielem] );

		const double XCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].X - localCoordinateValuesStartPoint[ielem].X ) * lengthX * 0.5;
		const double YCoordDifferenceOfSegment = ( localCoordinateValuesEndPoint[ielem].Y - localCoordinateValuesStartPoint[ielem].Y ) * lengthY * 0.5;

		if( integralXCompFirst ){
			calcInterpolatorVectorOfElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0, irhs, std::complex<double>( -XCoordDifferenceOfSegment , 0.0 ) );
			calcInterpolatorVectorOfElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0, irhs, std::complex<double>( -YCoordDifferenceOfSegment , 0.0 ) );
		}else{
			calcInterpolatorVectorOfElectricFieldYDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0, irhs, std::complex<double>( -YCoordDifferenceOfSegment , 0.0 ) );
			calcInterpolatorVectorOfElectricFieldXDirection( elememtsIncludingDipole[ielem], localCoordinateValuesEndPoint[ielem].X,   localCoordinateValuesEndPoint[ielem].Y,   -1.0, irhs, std::complex<double>( -XCoordDifferenceOfSegment , 0.0 ) );
		}

		double areaWithSign = 0.5 * fabs( XCoordDifferenceOfSegment * YCoordDifferenceOfSegment );
		if( !rotationDirectionPlus ){
			areaWithSign *= -1.0;
		}

		calcInterpolatorVectorOfRotatedElectricFieldZDirection( elememtsIncludingDipole[ielem], localCoordinateValuesStartPoint[ielem].X, localCoordinateValuesStartPoint[ielem].Y, -1.0, irhs, areaWithSign );

	}

}

// Calculate interpolator vector of difference of voltage
void Forward3DBrickElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
	const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs ){

	OutputFiles::m_logFile << "Error : calcInterpolatorVectorOfVoltageDifference for tetrahedral element can't be used for brick element ." << std::endl;
	exit(1);

}

// Set non-zero strucuture of matrix for forward calculation
void Forward3DBrickElement0thOrder::setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix ){

	const ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();

	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	int iElem(0);
	int elemID(0);
	int iEdge1(0);
	int iEdge2(0);
	int row(0);
	int col(0);
#ifdef _USE_OMP
	#pragma omp parallel for default(shared) \
		private( iElem, elemID, iEdge1, iEdge2, row, col )
#endif
	for( iElem = 0; iElem < nElem; ++iElem ){

		// [Attention] : You must use elemID instead of iElem from this line
		elemID = iElem;

		for( iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){

				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col < 0 ){
					continue;
				}

				if( col >= row ){// Store only upper triangle part
#ifdef _USE_OMP
					#pragma omp critical (addToMatrix)
#endif
					{
						matrix.setStructureByTripletFormat( row, col );
					}

				}

			}// iEdge2

		}// iEdge1

	}// iElem

}

// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward3DBrickElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ){

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	int iElem(0);
	int elemID(0);
	double lengX(0.0);
	double lengY(0.0);
	double lengZ(0.0);
	double jacobian(0.0);
	double sigma(0.0);
	double factor1(0.0);
	std::complex<double> factor2(0.0,0.0);
	int iEdge1(0);
	int iEdge2(0);
	int row(0);
	int col(0);
	int ip(0);
	double integral1(0.0);
	double integral2(0.0);
	std::complex<double> val(0.0,0.0);
	int loc(0);
	//------------------------------------------------------------ Start of parallel region >>>>>
//#ifdef _USE_OMP
//	#pragma omp parallel for default(shared) \
//		private( iElem, elemID, lengX, lengY, lengZ, jacobian, sigma, factor1, factor2, \
//			iEdge1, iEdge2, row, col, integral1, integral2, ip, val, loc )
//#endif
	for( iElem = 0; iElem < nElem; ++iElem ){

		// [Attention] : You must use elemID instead of iElem from this line
		elemID = iElem;

		//--- Calculate Jacobian
		lengX = m_MeshDataBrickElement.getEdgeLengthX( elemID );
		lengY = m_MeshDataBrickElement.getEdgeLengthY( elemID );
		lengZ = m_MeshDataBrickElement.getEdgeLengthZ( elemID );
		jacobian = 0.125 * lengX * lengY * lengZ;

		//--- Calculate omega * mu * sigma
		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

		factor1 = 1.0;
		//factor2 = std::complex<double>( omega * omega * CommonParameters::mu * CommonParameters::epsilon, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form
		factor2 = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form

		for( iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){

				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}

				integral1 = 0.0;
				integral2 = 0.0;
				for( ip = 0; ip < m_numIntegralPoints; ++ip ){
					integral1 += ( getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 )
						         * getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
								 + getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
						         * getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
								 + getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
						         * getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 ) )
							     * m_weights3D[ip];
					integral2 += ( getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
						         * getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
						         + getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
								 * getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
						         + getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
								 * getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 ) )
								 * m_weights3D[ip];
				}
				integral1 *= jacobian;
				integral2 *= jacobian;

				val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form

				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
//#ifdef _USE_OMP
//					#pragma omp critical (addToRhs)
//#endif
					{
						matrix.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
					}

				}else if( col >= row ){// Store only upper triangle part
					loc = matrix.checkAndGetLocationNonZeroValue( row, col );
//#ifdef _USE_OMP
//					#pragma omp critical (addToMatrix)
//#endif
					{
						//matrix.addNonZeroValues( row, col, val );// Add to matrix
						matrix.addNonZeroValuesWithoutSearchingLocation( loc, val );// Add to matrix
					}

				}

			}// iEdge2

		}// iEdge1		

	}
	//------------------------------------------------------------ End of parallel region <<<<<

}

//----- DO NOT DELETE FOR FUTURE USE >>>>>
//// Set non-zero strucuture of matrix for calculating derivatives
//void Forward3DBrickElement0thOrder::setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID, std::set<int>& nonZeroRowsAndCols ){
//
//	const ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();
//	if( ptrResistivityBlock->isFixedResistivityValue(blkID) ){
//		return;
//	}
//
//	const int iPol = getPolarizationCurrent();
//
//	const std::vector< std::pair<int,double> >& blk2Elem = ptrResistivityBlock->getBlockID2Elements(blkID);
//	const int nElem = static_cast<int>( blk2Elem.size() );
//
//	//------------------------------------------
//	//--- Components due to stiffness matrix ---
//	//------------------------------------------
//	int iElem(0);
//	int elemID(0);
//	int iEdge1(0);
//	int iEdge2(0);
//	int row(0);
//	int col(0);
////#ifdef _USE_OMP
////	#pragma omp parallel for default(shared) \
////		private( iElem, elemID, iEdge1, iEdge2, row, col )
////#endif
//	for( iElem = 0; iElem < nElem; ++iElem ){
//
//		// [Attention] : You must use elemID instead of iElem from this line
//		elemID = blk2Elem[iElem].first;
//
//		for( iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){
//
//			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
//			if( row < 0 ){
//				continue;
//			}
//
//			for( iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){
//
//				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
//				if( col < 0 ){
//					continue;
//				}
//
//				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
////#ifdef _USE_OMP
////					#pragma omp critical (addToRhs)
////#endif
//					{
//						nonZeroRowsAndCols.insert( row );
//					}
//
//				}else if( col >= row ){// Store only upper triangle part
////#ifdef _USE_OMP
////					#pragma omp critical (addToMatrix)
////#endif
//					{
//						matrix.setStructureByTripletFormat( row, col );
//						nonZeroRowsAndCols.insert( row );
//						nonZeroRowsAndCols.insert( col );
//					}
//
//				}
//
//			}// iEdge2
//
//		}// iEdge1
//
//	}// iElem
//
//}
//
//// Set non-zero values of matrix and right-hande side vector for calculating derivatives
//void Forward3DBrickElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID ){
//
//	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
//
//	if( pResistivityBlock->isFixedResistivityValue(blkID) ){
//		return;
//	}
//
//	const int iPol = getPolarizationCurrent();
//	const double freq = getFrequencyCurrent();
//	const double ln10 = 2.30258509299405;
//	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
//
//	const std::vector< std::pair<int,double> >& mdl2Elem = pResistivityBlock->getBlockID2Elements(blkID);
//	const int nElem = static_cast<int>( mdl2Elem.size() );
//
//	//------------------------------------------
//	//--- Components due to stiffness matrix ---
//	//------------------------------------------
//	int iElem(0);
//	int elemID(0);
//	double lengX(0.0);
//	double lengY(0.0);
//	double lengZ(0.0);
//	double jacobian(0.0);
//	double sigma(0.0);
//	double factor1(0.0);
//	std::complex<double> factor2(0.0,0.0);
//	int iEdge1(0);
//	int iEdge2(0);
//	int row(0);
//	int col(0);
//	int ip(0);
//	double integral1(0.0);
//	double integral2(0.0);
//	std::complex<double> val(0.0,0.0);
//	int loc(0);
//	//------------------------------------------------------------ Start of parallel region >>>>>
////#ifdef _USE_OMP
////	#pragma omp parallel for default(shared) \
////		private( iElem, elemID, lengX, lengY, lengZ, jacobian, sigma, factor1, factor2, \
////			iEdge1, iEdge2, row, col, integral1, integral2, ip, val, loc )
////#endif
//	for( iElem = 0; iElem < nElem; ++iElem ){
//
//		// [Attention] : You must use elemID instead of iElem from this line
//		elemID = mdl2Elem[iElem].first;
//
////----- debug >>>>>
//#ifdef _DEBUG_WRITE
//		std::cout << "blkID iElem elemID = " << blkID << " " << iElem << " " << elemID << std::endl;
//#endif
////----- debug <<<<<
//
//		//--- Calculate Jacobian
//		lengX = m_MeshDataBrickElement.getEdgeLengthX( elemID );
//		lengY = m_MeshDataBrickElement.getEdgeLengthY( elemID );
//		lengZ = m_MeshDataBrickElement.getEdgeLengthZ( elemID );
//		jacobian = 0.125 * lengX * lengY * lengZ;
//
//		//--- Calculate omega * mu * sigma
//		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
//
//		factor1 = 0.0;
//		factor2 = std::complex<double>( 0.0, - omega * CommonParameters::mu * sigma * ln10 * mdl2Elem[iElem].second );// exp(-i*omega*t) form
//
//		for( iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){
//
//			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
//			if( row < 0 ){
//				continue;
//			}
//
//			for( iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){
//
//				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
//				if( col <= DIRICHLET_BOUNDARY_ZERO_VALUE ){
//					continue;
//				}
//
//				integral1 = 0.0;
//				integral2 = 0.0;
//				for( ip = 0; ip < m_numIntegralPoints; ++ip ){
//					integral1 += ( getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 )
//						         * getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
//								 + getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
//						         * getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
//								 + getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
//						         * getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 ) )
//							     * m_weights3D[ip];
//					integral2 += ( getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
//						         * getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
//						         + getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
//								 * getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
//						         + getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
//								 * getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 ) )
//								 * m_weights3D[ip];
//				}
//				integral1 *= jacobian;
//				integral2 *= jacobian;
//
//				val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form
//
//				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
////#ifdef _USE_OMP
////					#pragma omp critical (addToRhs)
////#endif
//					{
//						matrix.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
//					}
//
//				}else if( col >= row ){// Store only upper triangle part
//					loc = matrix.checkAndGetLocationNonZeroValue( row, col );
////#ifdef _USE_OMP
////					#pragma omp critical (addToMatrix)
////#endif
//					{
//						matrix.addNonZeroValuesWithoutSearchingLocation( loc, val );// Add to matrix
//					}
//
//				}
//
//			}// iEdge2
//
//		}// iEdge1		
//
//	}
//	//------------------------------------------------------------ End of parallel region <<<<<
//
//}
//----- DO NOT DELETE FOR FUTURE USE <<<<<

// Calculate vector x of the reciprocity algorithm of Rodi (1976)
void Forward3DBrickElement0thOrder::calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows ){
	
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	if( pResistivityBlock->isFixedResistivityValue(blkID) ){
		return;
	}

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	const std::vector< std::pair<int,double> >& mdl2Elem = pResistivityBlock->getBlockID2Elements(blkID);
	const int nElem = static_cast<int>( mdl2Elem.size() );

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	int iElem(0);
	int elemID(0);
	double lengX(0.0);
	double lengY(0.0);
	double lengZ(0.0);
	double jacobian(0.0);
	double sigma(0.0);
	double factor1(0.0);
	std::complex<double> factor2(0.0,0.0);
	int iEdge1(0);
	int iEdge2(0);
	int row(0);
	int col(0);
	int ip(0);
	double integral1(0.0);
	double integral2(0.0);
	std::complex<double> val(0.0,0.0);
	int loc(0);
//#ifdef _USE_OMP
//	#pragma omp parallel for default(shared) \
//		private( iElem, elemID, lengX, lengY, lengZ, jacobian, sigma, factor1, factor2, \
//			iEdge1, iEdge2, row, col, integral1, integral2, ip, val, loc )
//#endif
	for( iElem = 0; iElem < nElem; ++iElem ){

		// [Attention] : You must use elemID instead of iElem from this line
		elemID = mdl2Elem[iElem].first;

//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "blkID iElem elemID = " << blkID << " " << iElem << " " << elemID << std::endl;
#endif
//----- debug <<<<<

		//--- Calculate Jacobian
		lengX = m_MeshDataBrickElement.getEdgeLengthX( elemID );
		lengY = m_MeshDataBrickElement.getEdgeLengthY( elemID );
		lengZ = m_MeshDataBrickElement.getEdgeLengthZ( elemID );
		jacobian = 0.125 * lengX * lengY * lengZ;

		//--- Calculate omega * mu * sigma
		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

		factor1 = 0.0;
		factor2 = std::complex<double>( 0.0, - omega * CommonParameters::mu * sigma * ln10 * mdl2Elem[iElem].second );// exp(-i*omega*t) form

		for( iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){

				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}

				integral1 = 0.0;
				integral2 = 0.0;
				for( ip = 0; ip < m_numIntegralPoints; ++ip ){
					integral1 += ( getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 )
						         * getShapeFuncRotatedX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
								 + getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
						         * getShapeFuncRotatedY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 )
								 + getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge1 ) 
						         * getShapeFuncRotatedZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], lengX, lengY, lengZ, iEdge2 ) )
							     * m_weights3D[ip];
					integral2 += ( getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
						         * getShapeFuncX( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
						         + getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
								 * getShapeFuncY( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 )
						         + getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge1 )
								 * getShapeFuncZ( m_xLocal[ip], m_yLocal[ip], m_zLocal[ip], iEdge2 ) )
								 * m_weights3D[ip];
				}
				integral1 *= jacobian;
				integral2 *= jacobian;

				val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form

				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					vecOut[row] -= val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ];
					nonZeroRows.push_back(row);

				}else{
					vecOut[row] -= val * vecIn[col];
					nonZeroRows.push_back(row);
				}

			}// iEdge2

		}// iEdge1		

	}

	std::sort(nonZeroRows.begin(), nonZeroRows.end());
	nonZeroRows.erase( std::unique( nonZeroRows.begin(), nonZeroRows.end() ), nonZeroRows.end() );

}

// Call function inputMeshData of the class MeshData
void Forward3DBrickElement0thOrder::callInputMeshData(){

	m_MeshDataBrickElement.inputMeshData();

}

// Get pointer to the class MeshData
const MeshData* Forward3DBrickElement0thOrder::getPointerToMeshData() const{

	return static_cast<const MeshData*>( &m_MeshDataBrickElement ) ;

}

// Get pointer to the class MeshDataBrickElement
const MeshDataBrickElement* Forward3DBrickElement0thOrder::getPointerToMeshDataBrickElement() const{

	return &m_MeshDataBrickElement;

}

// Get X component of shape function for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncX( const double xLocal, const double yLocal, const double zLocal, const int num ) const{

	switch( num ){
		case 0:
			return 0.25 * ( 1.0 - yLocal ) * ( 1.0 - zLocal );
			break;
		case 1:
			return 0.25 * ( 1.0 + yLocal ) * ( 1.0 - zLocal );
			break;
		case 2:
			return 0.25 * ( 1.0 - yLocal ) * ( 1.0 + zLocal );
			break;
		case 3:
			return 0.25 * ( 1.0 + yLocal ) * ( 1.0 + zLocal );
			break;
		case 4:
			return 0.0;
			break;
		case 5:
			return 0.0;
			break;
		case 6:
			return 0.0;
			break;
		case 7:
			return 0.0;
			break;
		case 8:
			return 0.0;
			break;
		case 9:
			return 0.0;
			break;
		case 10:
			return 0.0;
			break;
		case 11:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncX : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get Y component of shape function for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncY( const double xLocal, const double yLocal, const double zLocal, const int num ) const{

	switch( num ){
		case 0:
			return 0.0; 
			break;
		case 1:
			return 0.0; 
			break;
		case 2:
			return 0.0; 
			break;
		case 3:
			return 0.0; 
			break;
		case 4:
			return 0.25 * ( 1.0 - zLocal ) * ( 1.0 - xLocal );
			break;
		case 5:
			return 0.25 * ( 1.0 + zLocal ) * ( 1.0 - xLocal );
			break;
		case 6:
			return 0.25 * ( 1.0 - zLocal ) * ( 1.0 + xLocal );
			break;
		case 7:
			return 0.25 * ( 1.0 + zLocal ) * ( 1.0 + xLocal );
			break;
		case 8:
			return 0.0;
			break;
		case 9:
			return 0.0;
			break;
		case 10:
			return 0.0;
			break;
		case 11:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncY : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get Z component of shape function for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncZ( const double xLocal, const double yLocal, const double zLocal, const int num ) const{

	switch( num ){
		case 0:
			return 0.0; 
			break;
		case 1:
			return 0.0; 
			break;
		case 2:
			return 0.0; 
			break;
		case 3:
			return 0.0; 
			break;
		case 4:
			return 0.0;
			break;
		case 5:
			return 0.0;
			break;
		case 6:
			return 0.0;
			break;
		case 7:
			return 0.0;
			break;
		case 8:
			return 0.25 * ( 1.0 - xLocal ) * ( 1.0 - yLocal );
			break;
		case 9:
			return 0.25 * ( 1.0 + xLocal ) * ( 1.0 - yLocal );
			break;
		case 10:
			return 0.25 * ( 1.0 - xLocal ) * ( 1.0 + yLocal );
			break;
		case 11:
			return 0.25 * ( 1.0 + xLocal ) * ( 1.0 + yLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncZ : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get X component of shape function rotated for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncRotatedX( const double xLocal, const double yLocal, const double zLocal,
		const double lx, const double ly, const double lz, const int num ) const{

	switch( num ){
		case 0:
			return 0.0; 
			break;
		case 1:
			return 0.0; 
			break;
		case 2:
			return 0.0; 
			break;
		case 3:
			return 0.0; 
			break;
		case 4:
			return   0.5 / lz * ( 1.0 - xLocal );
			break;
		case 5:
			return - 0.5 / lz * ( 1.0 - xLocal );
			break;
		case 6:
			return   0.5 / lz * ( 1.0 + xLocal );
			break;
		case 7:
			return - 0.5 / lz * ( 1.0 + xLocal );
			break;
		case 8:
			return - 0.5 / ly * ( 1.0 - xLocal );
			break;
		case 9:
			return - 0.5 / ly * ( 1.0 + xLocal );
			break;
		case 10:
			return   0.5 / ly * ( 1.0 - xLocal );
			break;
		case 11:
			return   0.5 / ly * ( 1.0 + xLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedX : num = " << num << std::endl;
			exit(1);
			break;
	}
	
}

// Get Y component of shape function rotated for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncRotatedY( const double xLocal, const double yLocal, const double zLocal,
		const double lx, const double ly, const double lz, const int num ) const{

	switch( num ){
		case 0:
			return - 0.5 / lz * ( 1.0 - yLocal );
			break;
		case 1:
			return - 0.5 / lz * ( 1.0 + yLocal );
			break;
		case 2:
			return   0.5 / lz * ( 1.0 - yLocal );
			break;
		case 3:
			return   0.5 / lz * ( 1.0 + yLocal );
			break;
		case 4:
			return 0.0; 
			break;
		case 5:
			return 0.0; 
			break;
		case 6:
			return 0.0; 
			break;
		case 7:
			return 0.0; 
			break;
		case 8:
			return   0.5 / lx * ( 1.0 - yLocal );
			break;
		case 9:
			return - 0.5 / lx * ( 1.0 - yLocal );
			break;
		case 10:
			return   0.5 / lx * ( 1.0 + yLocal );
			break;
		case 11:
			return - 0.5 / lx * ( 1.0 + yLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedY : num = " << num << std::endl;
			exit(1);
			break;
	}
	
}

// Get Z component of shape function rotated for 0th order edge-based elements
inline double Forward3DBrickElement0thOrder::getShapeFuncRotatedZ( const double xLocal, const double yLocal, const double zLocal,
		const double lx, const double ly, const double lz, const int num ) const{

	switch( num ){
		case 0:
			return   0.5 / ly * ( 1.0 - zLocal );
			break;
		case 1:
			return - 0.5 / ly * ( 1.0 - zLocal );
			break;
		case 2:
			return   0.5 / ly * ( 1.0 + zLocal );
			break;
		case 3:
			return - 0.5 / ly * ( 1.0 + zLocal );
			break;
		case 4:
			return - 0.5 / lx * ( 1.0 - zLocal );
			break;
		case 5:
			return - 0.5 / lx * ( 1.0 + zLocal );
			break;
		case 6:
			return   0.5 / lx * ( 1.0 - zLocal );
			break;
		case 7:
			return   0.5 / lx * ( 1.0 + zLocal );
			break;
		case 8:
			return 0.0; 
			break;
		case 9:
			return 0.0; 
			break;
		case 10:
			return 0.0; 
			break;
		case 11:
			return 0.0; 
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedZ : num = " << num << std::endl;
			exit(1);
			break;
	}
	
}

// Calculate array converting local IDs to global ones
void Forward3DBrickElement0thOrder::calcArrayConvertLocalID2Global(){

	if( m_hasSetIDsLocal2Global ){// If array converting local edge IDs to global ones has already been set
		OutputFiles::m_logFile << "Warning : Array converting local edge IDs to global ones has already been set." << std::endl;
		return;
	}

	//--- Element number of each direction
	const int numElemX = m_MeshDataBrickElement.getNumElemX();
	const int numElemY = m_MeshDataBrickElement.getNumElemY();
	const int numElemZ = m_MeshDataBrickElement.getNumElemZ();

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	//--- Number of edges on boundary planes
	const int numEdgesXYPlane = m_MeshDataBrickElement.calcNumEdgesOnXYPlane();
	const int numEdgesYZPlane = m_MeshDataBrickElement.calcNumEdgesOnYZPlane();
	const int numEdgesZXPlane = m_MeshDataBrickElement.calcNumEdgesOnZXPlane();

	m_numOfEquation = numEdgesXYPlane * ( numElemZ + 1 )// Horizontal edges
	                + ( numElemX + 1 ) * ( numElemY + 1 ) * numElemZ;// Vertical edges

	//----------------------------------------------------------------
	//--- Calculate array converting local edge IDs to global ones ---
	//----------------------------------------------------------------
	if( m_IDsLocal2Global != NULL ){// Release memory
		const int num = sizeof( m_IDsLocal2Global ) / sizeof( m_IDsLocal2Global[0] );
		for( int i = 0; i < num; ++i ){
			delete [] m_IDsLocal2Global[i];
			m_IDsLocal2Global[i] = NULL;
		}
		delete [] m_IDsLocal2Global;
		m_IDsLocal2Global = NULL;
	}

	m_IDsLocal2Global = new int*[nElem];

	for( int iElem = 0; iElem < nElem; ++iElem ){

		const int offsetZ = numEdgesXYPlane + ( numElemX + 1 ) * ( numElemY + 1 );
		const int offsetY = 2 * numElemX + 1;

		const int iz =   iElem / ( numElemX * numElemY );
		const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
		const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

		m_IDsLocal2Global[iElem] = new int[12];

		m_IDsLocal2Global[iElem][0] = iz * offsetZ +   iy       * offsetY + ix;
		m_IDsLocal2Global[iElem][1] = iz * offsetZ + ( iy + 1 ) * offsetY + ix;

		m_IDsLocal2Global[iElem][2] = ( iz + 1 ) * offsetZ +   iy       * offsetY + ix;
		m_IDsLocal2Global[iElem][3] = ( iz + 1 ) * offsetZ + ( iy + 1 ) * offsetY + ix;

		m_IDsLocal2Global[iElem][4] = m_IDsLocal2Global[iElem][0] + numElemX;
		m_IDsLocal2Global[iElem][6] = m_IDsLocal2Global[iElem][0] + numElemX + 1;

		m_IDsLocal2Global[iElem][5] = m_IDsLocal2Global[iElem][2] + numElemX;
		m_IDsLocal2Global[iElem][7] = m_IDsLocal2Global[iElem][2] + numElemX + 1;

		m_IDsLocal2Global[iElem][8]  = iz * offsetZ + numEdgesXYPlane +   iy       * ( numElemX + 1 ) + ix;
		m_IDsLocal2Global[iElem][9]  = iz * offsetZ + numEdgesXYPlane +   iy       * ( numElemX + 1 ) + ix + 1;
		m_IDsLocal2Global[iElem][10] = iz * offsetZ + numEdgesXYPlane + ( iy + 1 ) * ( numElemX + 1 ) + ix;
		m_IDsLocal2Global[iElem][11] = iz * offsetZ + numEdgesXYPlane + ( iy + 1 ) * ( numElemX + 1 ) + ix + 1;

	}

	m_hasSetIDsLocal2Global = true;

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem IDsLocal2Global : " << iElem;
		for( int i = 0; i < 12; ++i ){
			std::cout << " " << m_IDsLocal2Global[iElem][i];
		}
		std::cout << std::endl;
	}
#endif
	//----- debug <<<<<

}

// Calculate array converting global IDs to the ones after degeneration
void Forward3DBrickElement0thOrder::calcArrayConvertIDsGlobal2AfterDegenerated(){

	const int iPol = getPolarizationCurrent();
	if( m_hasIDsGlobal2AfterDegenerated[iPol] ){
		OutputFiles::m_logFile << "Warning : Array converting global IDs to the ones after degeneration has already been set." << std::endl;
		return;
	}

	//--- Element number of each direction
	const int numElemX = m_MeshDataBrickElement.getNumElemX();
	const int numElemY = m_MeshDataBrickElement.getNumElemY();
	const int numElemZ = m_MeshDataBrickElement.getNumElemZ();

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	const int numEdgesXYPlane = m_MeshDataBrickElement.calcNumEdgesOnXYPlane();
	const int numEdgesYZPlane = m_MeshDataBrickElement.calcNumEdgesOnYZPlane();
	const int numEdgesZXPlane = m_MeshDataBrickElement.calcNumEdgesOnZXPlane();

	int nTmp = m_numOfEquation - 2 * numEdgesXYPlane;// Exclude edges on horizontal planes of the top and the bottom
	nTmp -= 2 * numEdgesYZPlane;// Exclude edges on the Y-Z planes of the boundary
	nTmp -= 2 * numEdgesZXPlane;// Exclude edges on the Z-X planes of the boundary
	nTmp += 4 * ( numElemX + numElemY + numElemZ );// Restore number of edges removed twice
	m_numOfEquationDegenerated = nTmp;// Number of equations after degeneration
	
	//-----------------------------------------------------------------
	//--- Initialize array converting local edge IDs to global ones ---
	//-----------------------------------------------------------------
	if( m_IDsGlobal2AfterDegenerated[iPol] != NULL ){// Release memory
		delete [] m_IDsGlobal2AfterDegenerated[iPol];
		m_IDsGlobal2AfterDegenerated[iPol] = NULL;
	}

	m_IDsGlobal2AfterDegenerated[iPol] = new int[m_numOfEquation];
	for( int i= 0; i < m_numOfEquation; ++i ){// Initialize
		m_IDsGlobal2AfterDegenerated[iPol][i] = i;
	}

	//----------------------------------------------------------- Ex Polarization >>>>>
	if( iPol == CommonParameters::EX_POLARIZATION ){//Ex Polarization

		//-----------------------------------------
		//--- Top of the model ( source field ) ---
		//-----------------------------------------
		for( int iElem = 0; iElem < numElemX * numElemY; ++iElem ){// Elements at the top
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][0] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][1] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][4] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][6] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
		}

		//---------------------------------------------------------------------------
		//--- Sides of the model ( EM fields obtained by 2D forward calculation ) ---
		//---------------------------------------------------------------------------
		{// Z-X plane ( Minus side )

			const int iPlane = MeshData::ZXMinus;
			for( int iElem = 0; iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][0] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][2] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][8] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][9] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;

				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}
			}
		}
		{// Z-X plane ( Plus side )

			const int iPlane = MeshData::ZXPlus;
			for( int iElem = numElemX * ( numElemY - 1 ); iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][1]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][3]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][10] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][11] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;

				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}

			}

		}
		{// Y-Z plane ( Minus side )

			const int iPlane = MeshData::YZMinus;
			for( int iElem = 0; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][4]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][5]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][8]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][10] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;

			}

		}
		{// Y-Z plane ( Plus side )

			const int iPlane = MeshData::YZPlus;
			for( int iElem = numElemX - 1; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][6]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][7]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][9]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][11] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			}

		}

		//---------------------------
		//--- Bottom of the model ---
		//---------------------------
		for( int iElem = numElemX * numElemY * ( numElemZ - 1 ); iElem < nElem; ++iElem ){// Elements at the bottom
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][2] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][3] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][5] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][7] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
		}

	//----------------------------------------------------------- Ex Polarization <<<<<
	}else{//Ey Polarization
	//----------------------------------------------------------- Ey Polarization >>>>>

		//-----------------------------------------
		//--- Top of the model ( source field ) ---
		//-----------------------------------------
		for( int iElem = 0; iElem < numElemX * numElemY; ++iElem ){// Elements at the top
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][0] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][1] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][4] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][6] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
		}

		//---------------------------------------------------------------------------
		//--- Sides of the model ( EM fields obtained by 2D forward calculation ) ---
		//---------------------------------------------------------------------------
		{// Y-Z plane ( Minus side )

			const int iPlane = MeshData::YZMinus;
			for( int iElem = 0; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Y-Z plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][4]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][5]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][8]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][10] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			}

		}
		{// Y-Z plane ( Plus side )

			const int iPlane = MeshData::YZPlus;
			for( int iElem = numElemX - 1; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Y-Z plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][6]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][7]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][9]  ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][11] ] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			}

		}
		{// Z-X plane ( Minus side )

			const int iPlane = MeshData::ZXMinus;
			for( int iElem = 0; iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][0] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][2] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][8] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][9] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;

				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}
			}
		}
		{// Z-X plane ( Plus side )

			const int iPlane = MeshData::ZXPlus;
			for( int iElem = numElemX * ( numElemY - 1 ); iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;

				// Edges on the Z-X plane
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][1]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][3]  ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][10] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][11] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				
				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}

			}

		}

		//---------------------------
		//--- Bottom of the model ---
		//---------------------------
		for( int iElem = numElemX * numElemY * ( numElemZ - 1 ); iElem < nElem; ++iElem ){// Elements at the bottom
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][2] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][3] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][5] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][7] ] = DIRICHLET_BOUNDARY_ZERO_VALUE;
		}

	}
	//----------------------------------------------------------- Ey Polarization <<<<<
		
	//------------------------------------------------------------------------------------
	//--- Converting global node IDs to the ones after degeneration in ascending order ---
	//--- excluding the nodes to be constrained                                        ---
	//------------------------------------------------------------------------------------
	int icount(0);
	for( int i = 0; i < m_numOfEquation; ++i ){
		if( m_IDsGlobal2AfterDegenerated[iPol][i] < 0 ){
			continue;
		}
		m_IDsGlobal2AfterDegenerated[iPol][i] = icount;
		++icount;
	}
	if( icount != m_numOfEquationDegenerated ){
		OutputFiles::m_logFile << "Error : Number of counter is not equal to m_numOfEquationDegenerated !! icount = "
			<< icount << ", m_numOfEquationDegenerated = " << m_numOfEquationDegenerated << std::endl;
		exit(1);
	}

//----- debug >>>>>
#ifdef _DEBUG_WRITE
for( int i = 0; i < m_numOfEquation; ++i ){
	std::cout << "i m_IDsGlobal2AfterDegenerated[iPol] : " << i << " " << m_IDsGlobal2AfterDegenerated[iPol][i] << std::endl;
}
#endif
//----- debug <<<<<

	//------------------------------------------------------------------------
	//--- Renumber global node IDs after degeneration by coordinate values ---
	//------------------------------------------------------------------------
	if( ( AnalysisControl::getInstance() )->getNumberingMethod() != AnalysisControl::NOT_ASSIGNED ){
		renumberNodes();
	}

	m_hasIDsGlobal2AfterDegenerated[iPol] = true;

}

// Calculate array converting global edge IDs non-zero electric field values specified to the edges
void Forward3DBrickElement0thOrder::calcArrayConvertIDGlobal2NonZeroValues(){

	const double sourceValueElectric = CommonParameters::sourceValueElectric;

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();

	//--- Element number of each direction
	const int numElemX = m_MeshDataBrickElement.getNumElemX();
	const int numElemY = m_MeshDataBrickElement.getNumElemY();
	const int numElemZ = m_MeshDataBrickElement.getNumElemZ();

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	//----------------------------------------------------------------------------
	//--- Calculating EM field of the boundary planes with 2D forward analysis ---
	//----------------------------------------------------------------------------
	OutputFiles::m_logFile << "# Calculating EM field of the boundary planes with 2D forward analysis.        " << std::endl;
	//for( int i = 0; i < 4; ++i ){
	//	OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << i << " ]-----------------------------" << std::endl;
	//	m_Fwd2DSquareElement[i][iPol]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataBrickElement );
	//}

	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex polarization

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DSquareElement[MeshData::ZXMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataBrickElement );

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DSquareElement[MeshData::ZXPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataBrickElement );

	}else{// Ey polarization

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DSquareElement[MeshData::YZMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataBrickElement );

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DSquareElement[MeshData::YZPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataBrickElement );

	}

	OutputFiles::m_logFile << "#------------------------------------------------------------------------------" << std::endl;
	//----------------------------------------------------------------------------

	// Initialize array converting global node IDs non-zero electric field values specified to the nodes
	m_globalID2NonZeroValues.clear();

	//----------------------------------------------------------- Ex Polarization >>>>>
	if( iPol == CommonParameters::EX_POLARIZATION ){//Ex Polarization

		//-----------------------------------------
		//--- Top of the model ( source field ) ---
		//-----------------------------------------
		for( int iElem = 0; iElem < numElemX * numElemY; ++iElem ){// Elements at the top
			m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][0], std::complex<double>(sourceValueElectric, 0.0) ) );
			m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][1], std::complex<double>(sourceValueElectric, 0.0) ) );
		}

		//---------------------------------------------------------------------------
		//--- Sides of the model ( EM fields obtained by 2D forward calculation ) ---
		//---------------------------------------------------------------------------
		{// Z-X plane ( Minus side )

			const int iPlane = MeshData::ZXMinus;
			for( int iElem = 0; iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;
				const int elemID2D = numElemX * iz + ix;

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][2], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 0 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][8], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 1 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][0], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 2 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][9], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 3 ) ) );

				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}
			}
		}
		{// Z-X plane ( Plus side )

			const int iPlane = MeshData::ZXPlus;
			for( int iElem = numElemX * ( numElemY - 1 ); iElem < nElem; ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;
				const int elemID2D = numElemX * iz + ix;

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][3],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 0 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][10], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 1 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][1],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 2 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][11], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 3 ) ) );

				if( ix == numElemX - 1 ){
					iElem += numElemX * ( numElemY - 1 ) + 1;
				}else{
					++iElem;
				}

			}

		}

	//----------------------------------------------------------- Ex Polarization <<<<<
	}else{//Ey Polarization
	//----------------------------------------------------------- Ey Polarization >>>>>

		//-----------------------------------------
		//--- Top of the model ( source field ) ---
		//-----------------------------------------
		for( int iElem = 0; iElem < numElemX * numElemY; ++iElem ){// Elements at the top
			m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][4], std::complex<double>(sourceValueElectric, 0.0) ) );
			m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][6], std::complex<double>(sourceValueElectric, 0.0) ) );
		}

		//---------------------------------------------------------------------------
		//--- Sides of the model ( EM fields obtained by 2D forward calculation ) ---
		//---------------------------------------------------------------------------
		{// Y-Z plane ( Minus side )

			const int iPlane = MeshData::YZMinus;
			for( int iElem = 0; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;
				const int elemID2D = numElemY * iz + iy;

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][5],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 0 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][8],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 1 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][4],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 2 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][10], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 3 ) ) );
			}

		}
		{// Y-Z plane ( Plus side )

			const int iPlane = MeshData::YZPlus;
			for( int iElem = numElemX - 1; iElem < nElem; iElem += numElemX ){
				const int iz =   iElem / ( numElemX * numElemY );
				const int iy = ( iElem % ( numElemX * numElemY ) ) / numElemX;
				const int ix = ( iElem % ( numElemX * numElemY ) ) % numElemX;
				const int elemID2D = numElemY * iz + iy;

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][7],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 0 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][9],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 1 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][6],  m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 2 ) ) );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][11], m_Fwd2DSquareElement[iPlane]->getSolutionFromLocalID( elemID2D, 3 ) ) );

			}

		}
		
	}
	//----------------------------------------------------------- Ey Polarization <<<<<

//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for ( std::map<int, std::complex<double> >::iterator itr = m_globalID2NonZeroValues.begin(); itr != m_globalID2NonZeroValues.end(); ++itr ){
		std::cout << "key value " << itr->first << " " << itr->second << std::endl;
	}
#endif
//----- debug <<<<<

}

// Renumber global node IDs after degeneration by coordinate values
void Forward3DBrickElement0thOrder::renumberNodes(){

	const int numberingMethod = ( AnalysisControl::getInstance() )->getNumberingMethod();

	if( numberingMethod == AnalysisControl::NOT_ASSIGNED ){
		return;
	}

	OutputFiles::m_logFile << "# Renumbering nodes. " << ( AnalysisControl::getInstance() )->outputElapsedTime() << std::endl;

	const int iPol = getPolarizationCurrent();

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	//--------------------------------------
	//--- Calculate coordinates of nodes ---
	//--------------------------------------
	double* IDsGlobalDegenerated2Xcoord = new double[ m_numOfEquationDegenerated ];
	double* IDsGlobalDegenerated2Ycoord = new double[ m_numOfEquationDegenerated ];
	double* IDsGlobalDegenerated2Zcoord = new double[ m_numOfEquationDegenerated ];
	{
		bool* alreadyFound = new bool[ m_numOfEquationDegenerated ];
		for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
			alreadyFound[i] = false;
		}

		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int node0 = m_MeshDataBrickElement.getNodesOfElements( iElem, 0 );
			const int node1 = m_MeshDataBrickElement.getNodesOfElements( iElem, 1 );
			const int node3 = m_MeshDataBrickElement.getNodesOfElements( iElem, 3 );
			const int node4 = m_MeshDataBrickElement.getNodesOfElements( iElem, 4 );

			const double originX = m_MeshDataBrickElement.getXCoordinatesOfNodes( node0 );
			const double lengthX = m_MeshDataBrickElement.getXCoordinatesOfNodes( node1 ) - originX;
			const double originY = m_MeshDataBrickElement.getYCoordinatesOfNodes( node0 );
			const double lengthY = m_MeshDataBrickElement.getYCoordinatesOfNodes( node3 ) - originY;
			const double originZ = m_MeshDataBrickElement.getZCoordinatesOfNodes( node0 );
			const double lengthZ = m_MeshDataBrickElement.getZCoordinatesOfNodes( node4 ) - originZ;
			
			int icount(0);
			{// Nodes on the edges parallel to X coordinate
				for( int iz = 0; iz < 2; ++iz ){
					for( int iy = 0; iy < 2; ++iy ){

						const int nodeIDGlobalDegenerated = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][icount] ];
						if( nodeIDGlobalDegenerated >= 0 && alreadyFound[ nodeIDGlobalDegenerated ] == false ){
							const double ratioY[2] = { 0.0, 1.0 };
							const double ratioZ[2] = { 0.0, 1.0 };
							IDsGlobalDegenerated2Xcoord[ nodeIDGlobalDegenerated ] = originX + lengthX * 0.5;
							IDsGlobalDegenerated2Ycoord[ nodeIDGlobalDegenerated ] = originY + lengthY * ratioY[iy];
							IDsGlobalDegenerated2Zcoord[ nodeIDGlobalDegenerated ] = originZ + lengthZ * ratioZ[iz];
							alreadyFound[ nodeIDGlobalDegenerated ] = true;
						}
						++icount;

					}
				}
			}

			{// Nodes on the edges parallel to Y coordinate
				for( int ix = 0; ix < 2; ++ix ){
					for( int iz = 0; iz < 2; ++iz ){

						const int nodeIDGlobalDegenerated = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][icount] ];
						if( nodeIDGlobalDegenerated >= 0 && alreadyFound[ nodeIDGlobalDegenerated ] == false ){
							const double ratioX[2] = { 0.0, 1.0 };
							const double ratioZ[2] = { 0.0, 1.0 };
							IDsGlobalDegenerated2Xcoord[ nodeIDGlobalDegenerated ] = originX + lengthX * ratioX[ix];
							IDsGlobalDegenerated2Ycoord[ nodeIDGlobalDegenerated ] = originY + lengthY * 0.5;
							IDsGlobalDegenerated2Zcoord[ nodeIDGlobalDegenerated ] = originZ + lengthZ * ratioZ[iz];
							alreadyFound[ nodeIDGlobalDegenerated ] = true;
						}
						++icount;

					}
				}
			}

			{// Nodes on the edges parallel to Z coordinate
				for( int iy = 0; iy < 2; ++iy ){
					for( int ix = 0; ix < 2; ++ix ){

						const int nodeIDGlobalDegenerated = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][icount] ];
						if( nodeIDGlobalDegenerated >= 0 && alreadyFound[ nodeIDGlobalDegenerated ] == false ){
							const double ratioX[2] = { 0.0, 1.0 };
							const double ratioY[2] = { 0.0, 1.0 };
							IDsGlobalDegenerated2Xcoord[ nodeIDGlobalDegenerated ] = originX + lengthX * ratioX[ix];
							IDsGlobalDegenerated2Ycoord[ nodeIDGlobalDegenerated ] = originY + lengthY * ratioY[iy];
							IDsGlobalDegenerated2Zcoord[ nodeIDGlobalDegenerated ] = originZ + lengthZ * 0.5;
							alreadyFound[ nodeIDGlobalDegenerated ] = true;
						}
						++icount;

					}
				}
			}

		}// iElem

		delete [] alreadyFound;

	}

	// Array containing global IDs after denegerated except minus values. Elements of this array are sorted by coordinate values.
	int* IDsNew2Old = new int[ m_numOfEquationDegenerated ];
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		IDsNew2Old[i] = i;
	}

//----- debug >>>>>
#ifdef _DEBUG_WRITE
for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
std::cout << "Before i IDsNew2Old X Y Z : "	<< i << " " << IDsNew2Old[i]
											<< " " << IDsGlobalDegenerated2Xcoord[IDsNew2Old[i]]
											<< " " << IDsGlobalDegenerated2Ycoord[IDsNew2Old[i]]
											<< " " << IDsGlobalDegenerated2Zcoord[IDsNew2Old[i]] << std::endl;
}
#endif
//----- debug <<<<<

	switch( numberingMethod ){
		case AnalysisControl::XYZ:
			quickSortThreeKeys( m_numOfEquationDegenerated, IDsNew2Old, IDsGlobalDegenerated2Xcoord, IDsGlobalDegenerated2Ycoord, IDsGlobalDegenerated2Zcoord );
			break;
		case AnalysisControl::YZX:
			quickSortThreeKeys( m_numOfEquationDegenerated, IDsNew2Old, IDsGlobalDegenerated2Ycoord, IDsGlobalDegenerated2Zcoord, IDsGlobalDegenerated2Xcoord );
			break;
		case AnalysisControl::ZXY:
			quickSortThreeKeys( m_numOfEquationDegenerated, IDsNew2Old, IDsGlobalDegenerated2Zcoord, IDsGlobalDegenerated2Xcoord, IDsGlobalDegenerated2Ycoord );
			break;
	}

//----- debug >>>>>
#ifdef _DEBUG_WRITE
for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
std::cout << "After  i IDsNew2Old X Y Z : "	<< i << " " << IDsNew2Old[i]
											<< " " << IDsGlobalDegenerated2Xcoord[IDsNew2Old[i]]
											<< " " << IDsGlobalDegenerated2Ycoord[IDsNew2Old[i]]
											<< " " << IDsGlobalDegenerated2Zcoord[IDsNew2Old[i]] << std::endl;
}
#endif
//----- debug <<<<<

	delete [] IDsGlobalDegenerated2Xcoord;
	delete [] IDsGlobalDegenerated2Ycoord;
	delete [] IDsGlobalDegenerated2Zcoord;

	int* IDsOld2New = new int[ m_numOfEquationDegenerated ];
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		IDsOld2New[ IDsNew2Old[i] ] = i;
	}
	delete [] IDsNew2Old;

//----- debug >>>>>
#ifdef _DEBUG_WRITE
for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
std::cout << " i IDsOld2New : "	<< i << " " << IDsOld2New[i] << std::endl;
}
#endif
//----- debug <<<<<

	int icount(0);
	for( int i = 0; i < m_numOfEquation; ++i ){
		if( m_IDsGlobal2AfterDegenerated[iPol][i] < 0 ){
			continue;
		}
		m_IDsGlobal2AfterDegenerated[iPol][i] = IDsOld2New[ icount ];
		++icount;
	}
	delete [] IDsOld2New;

}

//// Calculate flag specifing whether rotation direction of integral route is positive or not
//bool Forward3DBrickElement0thOrder::doesRotationDirectionPlus( const CommonParameters::locationXY& startPoint, const CommonParameters::locationXY& endPoint,
//	bool& integralXCompFirst ) const{
//
//	const double eps = 1.0e-12;
//
//	CommonParameters::locationXY sharedPoint = { 0.0, 0.0 };
//
//	if( ( fabs( startPoint.X - 1.0 ) < eps && fabs( endPoint.Y - 1.0 ) < eps ) ||
//			( fabs( startPoint.Y - 1.0 ) < eps && fabs( endPoint.X - 1.0 ) < eps ) ){
//		sharedPoint.X = 1.0;
//		sharedPoint.Y = 1.0;
//	}else if( ( fabs( startPoint.X - 1.0 ) < eps && fabs( endPoint.Y + 1.0 ) < eps ) ||
//				( fabs( startPoint.Y + 1.0 ) < eps && fabs( endPoint.X - 1.0 ) < eps ) ){
//		sharedPoint.X = 1.0;
//		sharedPoint.Y = -1.0;
//	}else if( ( fabs( startPoint.X + 1.0 ) < eps && fabs( endPoint.Y + 1.0 ) < eps ) ||
//				( fabs( startPoint.Y + 1.0 ) < eps && fabs( endPoint.X + 1.0 ) < eps ) ){
//		sharedPoint.X = -1.0;
//		sharedPoint.Y = -1.0;
//	}else if( ( fabs( startPoint.X + 1.0 ) < eps && fabs( endPoint.Y - 1.0 ) < eps ) ||
//				( fabs( startPoint.Y - 1.0 ) < eps && fabs( endPoint.X + 1.0 ) < eps ) ){
//		sharedPoint.X = -1.0;
//		sharedPoint.Y = 1.0;
//	}
//
//	const double outerProduct = ( sharedPoint.X - startPoint.X )*( endPoint.Y - startPoint.Y ) - ( sharedPoint.Y - startPoint.Y )*( endPoint.X - startPoint.X );
//
//	if( outerProduct ){
//		return true;
//	}else{
//		return false;
//	}
//
//}

// Calculate flag specifing whether integral X component first
bool Forward3DBrickElement0thOrder::doesIntegralXCompFirst( const CommonParameters::locationXY& startPoint, const CommonParameters::locationXY& endPoint,
	bool& rotationDirectionPlus ) const{

	const double eps = 1.0e-12;

	CommonParameters::locationXY sharedPoint = { 0.0, 0.0 };

	bool intersectTwoEdges = true;
	if( ( fabs( startPoint.X - 1.0 ) < eps && fabs( endPoint.Y - 1.0 ) < eps ) ||
			( fabs( startPoint.Y - 1.0 ) < eps && fabs( endPoint.X - 1.0 ) < eps ) ){
		sharedPoint.X = 1.0;
		sharedPoint.Y = 1.0;
	}else if( ( fabs( startPoint.X - 1.0 ) < eps && fabs( endPoint.Y + 1.0 ) < eps ) ||
				( fabs( startPoint.Y + 1.0 ) < eps && fabs( endPoint.X - 1.0 ) < eps ) ){
		sharedPoint.X = 1.0;
		sharedPoint.Y = -1.0;
	}else if( ( fabs( startPoint.X + 1.0 ) < eps && fabs( endPoint.Y + 1.0 ) < eps ) ||
				( fabs( startPoint.Y + 1.0 ) < eps && fabs( endPoint.X + 1.0 ) < eps ) ){
		sharedPoint.X = -1.0;
		sharedPoint.Y = -1.0;
	}else if( ( fabs( startPoint.X + 1.0 ) < eps && fabs( endPoint.Y - 1.0 ) < eps ) ||
				( fabs( startPoint.Y - 1.0 ) < eps && fabs( endPoint.X + 1.0 ) < eps ) ){
		sharedPoint.X = -1.0;
		sharedPoint.Y = 1.0;
	}else{// Not intersect two edges => integral X component first
		sharedPoint.X = startPoint.X;
		sharedPoint.Y = endPoint.Y;
		intersectTwoEdges = false;
	}

	const double outerProduct = ( sharedPoint.X - startPoint.X )*( endPoint.Y - startPoint.Y ) - ( sharedPoint.Y - startPoint.Y )*( endPoint.X - startPoint.X );

	if( outerProduct > 0 ){
		rotationDirectionPlus = true;
	}else{
		rotationDirectionPlus = false;
	}

//----- debug >>>>>
#ifdef _DEBUG_WRITE
	std::cout << "sharedPoint.X sharedPoint.Y : " << sharedPoint.X << " " << sharedPoint.Y << std::endl;
	std::cout << "outerProduct rotationDirectionPlus : " << outerProduct << " " << rotationDirectionPlus << std::endl;
#endif
//----- debug <<<<<

	if( intersectTwoEdges ){
		if( sharedPoint.X * sharedPoint.Y < 0.0 ){
			if( rotationDirectionPlus ){
				return true;
			}else{
				return false;
			}
		}else{
			if( rotationDirectionPlus ){
				return false;
			}else{
				return true;
			}
		}
	}else{// Not intersect two edges => integral X component first
		return false;
	}

}


//// Output results of forward calculation to VTK file
//void Forward3DBrickElement0thOrder::outputResultToVTK() const{
//
//	if( !OutputFiles::m_vtkFile.is_open() ){
//		return;
//	}
//	
//	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
//
//	std::string stringPolarization;
//	const int iPol = getPolarizationCurrent();
//	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex Polarization
//		stringPolarization = "Ex_polarization";
//	}else{// Ey Polarization
//		stringPolarization = "Ey_polarization";
//	}
//
//	//--- Total element number
//	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
//	//const int nElem = pMeshDataBrickElement->m_numElemTotal;
//	const int nElem = m_MeshDataBrickElement.getNumElemTotal();
//
//	const double freq = getFrequencyCurrent();
//	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK) ){// Output electric field vector
//
//		OutputFiles::m_vtkFile << "VECTORS real_part_of_electric_field_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const float Ex = static_cast<float>( real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Ey = static_cast<float>( real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Ez = static_cast<float>( real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
//		}
//
//		OutputFiles::m_vtkFile << "VECTORS imaginary_part_of_electric_field_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const float Ex = static_cast<float>( imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Ey = static_cast<float>( imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Ez = static_cast<float>( imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
//		}
//
//	}
//
//	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK) ){// Output magnetic field vector
//		OutputFiles::m_vtkFile << "VECTORS real_part_of_magnetic_field_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const float Hx = static_cast<float>( real( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Hy = static_cast<float>( real( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Hz = static_cast<float>( real( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
//		}
//
//		OutputFiles::m_vtkFile << "VECTORS imaginary_part_of_magnetic_field_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const float Hx = static_cast<float>( imag( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Hy = static_cast<float>( imag( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float Hz = static_cast<float>( imag( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
//		}
//	}
//
//	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK) ){// Output corrent density
//		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();
//		OutputFiles::m_vtkFile << "VECTORS real_part_of_current_density_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
//			const float jx = static_cast<float>( sigma * real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float jy = static_cast<float>( sigma * real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float jz = static_cast<float>( sigma * real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
//		}
//
//		OutputFiles::m_vtkFile << "VECTORS imaginary_part_of_current_density_frequency_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
//		for( int iElem = 0 ; iElem < nElem; ++iElem ){
//			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
//			const float jx = static_cast<float>( sigma * imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float jy = static_cast<float>( sigma * imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			const float jz = static_cast<float>( sigma * imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
//			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
//		}
//	}
//
//	
//}
//

// Output results of forward calculation to VTK file
void Forward3DBrickElement0thOrder::outputResultToVTK() const{

	if( !OutputFiles::m_vtkFile.is_open() ){
		return;
	}
	
	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();

	std::string stringPolarization;
	const int iPol = getPolarizationCurrent();
	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex Polarization
		//stringPolarization = "Ex_polarization";
		stringPolarization = "ExPol";
	}else{// Ey Polarization
		//stringPolarization = "Ey_polarization";
		stringPolarization = "EyPol";
	}

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	const double freq = getFrequencyCurrent();
	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK ) ){// Output electric field vector
		OutputFiles::m_vtkFile << "VECTORS Re(E)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Ex = static_cast<float>( real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Ey = static_cast<float>( real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Ez = static_cast<float>( real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(E)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Ex = static_cast<float>( imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Ey = static_cast<float>( imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Ez = static_cast<float>( imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
		}

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK ) ){// Output magnetic field vector
		OutputFiles::m_vtkFile << "VECTORS Re(H)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Hx = static_cast<float>( real( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Hy = static_cast<float>( real( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Hz = static_cast<float>( real( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(H)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Hx = static_cast<float>( imag( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Hy = static_cast<float>( imag( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float Hz = static_cast<float>( imag( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
		}
	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_CURRENT_DENSITY ) ){// Output corrent density
		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();
		OutputFiles::m_vtkFile << "VECTORS Re(j)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			const float jx = static_cast<float>( sigma * real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float jy = static_cast<float>( sigma * real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float jz = static_cast<float>( sigma * real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(j)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			const float jx = static_cast<float>( sigma * imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float jy = static_cast<float>( sigma * imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			const float jz = static_cast<float>( sigma * imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
		}
	}

	
}

// Output results of forward calculation to binary file
void Forward3DBrickElement0thOrder::outputResultToBinary( const int iFreq, const int iPol ) const{
	
	const std::string stringPolarization = iPol == CommonParameters::EX_POLARIZATION ? "ExPol" : "EyPol";

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();

	//--- Total element number
	const int nElem = m_MeshDataBrickElement.getNumElemTotal();

	const double freq = getFrequencyCurrent();

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK ) ){
		// Output real part of electric field vector

		std::ostringstream oss;
		oss << "ReE_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Real(E) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ex = static_cast<float>( real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ex, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ey = static_cast<float>( real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ey, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ez = static_cast<float>( real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ez, sizeof( float ) );
		}

		fout.close();

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK ) ){
		// Output imaginary part of electric field vector

		std::ostringstream oss;
		oss << "ImE_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Imag(E) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ex = static_cast<float>( imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ex, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ey = static_cast<float>( imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ey, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Ez = static_cast<float>( imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Ez, sizeof( float ) );
		}

		fout.close();

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK ) ){
		// Output real part of magnetic field vector

		std::ostringstream oss;
		oss << "ReH_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Real(H) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hx = static_cast<float>( real( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hx, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hy = static_cast<float>( real( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hy, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hz = static_cast<float>( real( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hz, sizeof( float ) );
		}

		fout.close();

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK ) ){
		// Output imaginary part of magnetic field vector

		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();

		std::ostringstream oss;
		oss << "ImH_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Imag(H) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hx = static_cast<float>( imag( calcValueMagneticFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hx, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hy = static_cast<float>( imag( calcValueMagneticFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hy, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			float Hz = static_cast<float>( imag( calcValueMagneticFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &Hz, sizeof( float ) );
		}

		fout.close();

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_CURRENT_DENSITY ) ){
		// Output real part of corrent density

		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();

		std::ostringstream oss;
		oss << "Rej_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Real(j) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jx = static_cast<float>( sigma * real( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jx, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jy = static_cast<float>( sigma * real( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jy, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jz = static_cast<float>( sigma * real( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jz, sizeof( float ) );
		}

		fout.close();

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_CURRENT_DENSITY ) ){
		// Output imaginary part of corrent density

		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();

		std::ostringstream oss;
		oss << "Imj_Freq" << iFreq << "_" << stringPolarization << ".iter" << pAnalysisControl->getIterationNumCurrent();
		std::ofstream fout;
		fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

		char line[80];
		std::ostringstream ossTitle;
		ossTitle << "Imag(j) of " << stringPolarization << " at " << freq << " [Hz]";
		strcpy( line, ossTitle.str().c_str() );
		fout.write( line, 80 );

		strcpy( line, "part" );
		fout.write( line, 80 );

		int ibuf(1);
		fout.write( (char*) &ibuf, sizeof( int ) );

		strcpy( line, "hexa8" );
		fout.write( line, 80 );

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jx = static_cast<float>( sigma * imag( calcValueElectricFieldXDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jx, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jy = static_cast<float>( sigma * imag( calcValueElectricFieldYDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jy, sizeof( float ) );
		}

		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			float jz = static_cast<float>( sigma * imag( calcValueElectricFieldZDirection( iElem, 0.0, 0.0, 0.0 ) ) );
			fout.write( (char*) &jz, sizeof( float ) );
		}

		fout.close();

	}

}

//// Get total number of element
//int Forward3DBrickElement0thOrder::getNumElemTotal() const{
//	
//	return m_MeshDataBrickElement.getNumElemTotal();
//
//}
