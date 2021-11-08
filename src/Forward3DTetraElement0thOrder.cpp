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
#include "Forward3DTetraElement0thOrder.h"
#include "MeshDataTetraElement.h"
#include "ResistivityBlock.h"
#include "CommonParameters.h"
#include "OutputFiles.h"
#include "Forward2DTriangleElement0thOrderEdgeBased.h"
#include "ObservedData.h"
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <assert.h>

const double Forward3DTetraElement0thOrder::m_eps = 1.0e-12;

Forward3DTetraElement0thOrder::Forward3DTetraElement0thOrder():
	Forward3D(),
	m_signInversion(NULL)
{

	m_uCoord[0] = 0.5854101966249685;
	m_uCoord[1] = 0.1381966011250105;
	m_uCoord[2] = 0.1381966011250105;
	m_uCoord[3] = 0.1381966011250105;

	m_vCoord[0] = 0.1381966011250105;
	m_vCoord[1] = 0.5854101966249685;
	m_vCoord[2] = 0.1381966011250105;
	m_vCoord[3] = 0.1381966011250105;

	m_wCoord[0] = 0.1381966011250105;
	m_wCoord[1] = 0.1381966011250105;
	m_wCoord[2] = 0.5854101966249685;
	m_wCoord[3] = 0.1381966011250105;

	m_weights[0] = 0.041666666666666666;
	m_weights[1] = 0.041666666666666666;
	m_weights[2] = 0.041666666666666666;
	m_weights[3] = 0.041666666666666666;


	for( int i = 0; i < 4; ++i ){
		//for( int j = 0; j < 2; ++j ){
		//	m_Fwd2DTriangleElement[i][j] = NULL;
		//}
		m_Fwd2DTriangleElement[i] = NULL;
	}

	//m_Fwd2DTriangleElement[MeshData::YZMinus][CommonParameters::EX_POLARIZATION] = new Forward2DTriangleElement1stOrderNodeBased( MeshData::YZMinus, CommonParameters::EX_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::YZPlus ][CommonParameters::EX_POLARIZATION] = new Forward2DTriangleElement1stOrderNodeBased( MeshData::YZPlus,  CommonParameters::EX_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::ZXMinus][CommonParameters::EY_POLARIZATION] = new Forward2DTriangleElement1stOrderNodeBased( MeshData::ZXMinus, CommonParameters::EY_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::ZXPlus ][CommonParameters::EY_POLARIZATION] = new Forward2DTriangleElement1stOrderNodeBased( MeshData::ZXPlus,  CommonParameters::EY_POLARIZATION );

	//m_Fwd2DTriangleElement[MeshData::YZMinus][CommonParameters::EY_POLARIZATION] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::YZMinus, CommonParameters::EY_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::YZPlus ][CommonParameters::EY_POLARIZATION] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::YZPlus,  CommonParameters::EY_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::ZXMinus][CommonParameters::EX_POLARIZATION] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::ZXMinus, CommonParameters::EX_POLARIZATION );
	//m_Fwd2DTriangleElement[MeshData::ZXPlus ][CommonParameters::EX_POLARIZATION] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::ZXPlus,  CommonParameters::EX_POLARIZATION );

	m_Fwd2DTriangleElement[MeshData::ZXMinus] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::ZXMinus, CommonParameters::EX_POLARIZATION );
	m_Fwd2DTriangleElement[MeshData::ZXPlus ] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::ZXPlus,  CommonParameters::EX_POLARIZATION );

	m_Fwd2DTriangleElement[MeshData::YZMinus] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::YZMinus, CommonParameters::EY_POLARIZATION );
	m_Fwd2DTriangleElement[MeshData::YZPlus ] = new Forward2DTriangleElement0thOrderEdgeBased( MeshData::YZPlus,  CommonParameters::EY_POLARIZATION );

}

//Destructer
Forward3DTetraElement0thOrder::~Forward3DTetraElement0thOrder(){

	for( int i = 0; i < 4; ++i ){
		//for( int j = 0; j < 2; ++j ){
		//	if( m_Fwd2DTriangleElement[i][j] != NULL ){
		//		delete m_Fwd2DTriangleElement[i][j];
		//	}
		//}
		if( m_Fwd2DTriangleElement[i] != NULL ){
			delete m_Fwd2DTriangleElement[i];
		}
	}

	if( m_signInversion != NULL ){
		const int num = sizeof( m_signInversion ) / sizeof( m_signInversion[0] );
		for( int i = 0; i < num; ++i ){
			delete [] m_signInversion[i];
			m_signInversion[i] = NULL;
		}
		delete [] m_signInversion;
		m_signInversion = NULL;
	}

}

//Run 3D forward calculation
void Forward3DTetraElement0thOrder::forwardCalculation( const double freq, const int iPol ){

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

	//-----------------------------------------------------------------------------
	//--- Set Number of equations and array converting local IDs to global ones ---
	//-----------------------------------------------------------------------------
	if( !m_hasSetIDsLocal2Global ){
		calcArrayConvertLocalID2Global();
	}
	
	if( !m_hasIDsGlobal2AfterDegenerated[iPol] ){
		calcArrayConvertIDsGlobal2AfterDegenerated();
	}

	OutputFiles::m_logFile << "# Number of equation = " << m_numOfEquation
		<< ", Number of equation after degeneration = " << m_numOfEquationDegenerated << std::endl;

	// Calculate array converting global edge IDs non-zero electric field values specified to the edges
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
	//----- debug <<<<<6

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[ m_numOfEquation ];

	bool* alreadyFound = new bool[ m_numOfEquation ];
	for( int i = 0; i < m_numOfEquation; ++i ){
		alreadyFound[i] = false;
	}

	const int nElem = m_MeshDataTetraElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){

		for( int iEdge = 0; iEdge < 6; ++iEdge ){

			const int iNum = m_IDsLocal2Global[iElem][iEdge];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ];

				if( iNumDegenerated == Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
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

}

// Calculate X component of electric field
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ]; 
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	return valU * std::complex<double>(invJacobMat.mat11/detJacob, 0.0) + valV * std::complex<double>(invJacobMat.mat12/detJacob, 0.0) + valW * std::complex<double>(invJacobMat.mat13/detJacob, 0.0);

}

// Calculate Y component of electric field
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ]; 
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	return valU * std::complex<double>(invJacobMat.mat21/detJacob, 0.0) + valV * std::complex<double>(invJacobMat.mat22/detJacob, 0.0) + valW * std::complex<double>(invJacobMat.mat23/detJacob, 0.0);

}

// Calculate Z component of electric field
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ]; 
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	return valU * std::complex<double>(invJacobMat.mat31/detJacob, 0.0) + valV * std::complex<double>(invJacobMat.mat32/detJacob, 0.0) + valW * std::complex<double>(invJacobMat.mat33/detJacob, 0.0);

}

// Calculate X component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Tangential electric field is multiply  by the inverse of the inclination in the X direction
	return calcValueElectricFieldTangentialX(iElem, iFace, uCoord, vCoord) / dLengdX;

}

// Calculate Y component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Tangential electric field is multiply  by the inverse of the inclination in the Y direction
	return calcValueElectricFieldTangentialY(iElem, iFace, uCoord, vCoord) / dLengdY;

}

// Calculate tangential electric field directed to X from all edges of owner element
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Calculate X component of electric field
	return calcValueElectricFieldXDirection(iElem, xLocal, yLocal, zLocal) * dLengdX;

}

// Calculate tangential electric field directed to Y from all edges of owner element
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Calculate Y component of electric field
	return calcValueElectricFieldYDirection(iElem, xLocal, yLocal, zLocal) * dLengdY;

}

// Calculate tangential electric field directed to X
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);

	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][iEdge] ]; 
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			valU -= val * std::complex<double>( get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV -= val * std::complex<double>( get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV += val * std::complex<double>( get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}
	}

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 

	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	return valU * std::complex<double>(jacobMat.mat22/detJacob, 0.0) - valV * std::complex<double>(jacobMat.mat12/detJacob, 0.0);

}

// Calculate tangential electric field directed to Y
std::complex<double> Forward3DTetraElement0thOrder::calcValueElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);

	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][iEdge] ]; 
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			valU -= val * std::complex<double>( get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV -= val * std::complex<double>( get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV += val * std::complex<double>( get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}
	}

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 

	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	return valV * std::complex<double>(jacobMat.mat11/detJacob, 0.0) - valU * std::complex<double>(jacobMat.mat21/detJacob, 0.0);

}

// Calculate Z component of rotated electric field
std::complex<double> Forward3DTetraElement0thOrder::calcValueRotatedElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ];
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	return valU * std::complex<double>(jacobMat.mat13/detJacob, 0.0) + valV * std::complex<double>(jacobMat.mat23/detJacob, 0.0) + valW * std::complex<double>(jacobMat.mat33/detJacob, 0.0);
	
}

// Calculate normal component of rotated electric field
std::complex<double> Forward3DTetraElement0thOrder::calcValueRotatedElectricFieldNormal( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> val(0.0, 0.0);

	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			val -= m_solution[ m_IDsLocal2Global[iElem][iEdge] ] * std::complex<double>( get2DShapeFuncRotated() * length, 0.0 );
		}else{
			val += m_solution[ m_IDsLocal2Global[iElem][iEdge] ] * std::complex<double>( get2DShapeFuncRotated() * length, 0.0 );
		}
	}

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 

	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	val /= std::complex<double>( detJacob, 0.0);
	return val;

}

// Calculate X component of magnetic field
std::complex<double> Forward3DTetraElement0thOrder::calcValueMagneticFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ];
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	std::complex<double> val = valU * std::complex<double>(jacobMat.mat11/detJacob, 0.0)
							 + valV * std::complex<double>(jacobMat.mat21/detJacob, 0.0)
							 + valW * std::complex<double>(jacobMat.mat31/detJacob, 0.0);

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Y component of magnetic field
std::complex<double> Forward3DTetraElement0thOrder::calcValueMagneticFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);
	std::complex<double> valW(0.0, 0.0);

	for( int i = 0; i < 6; ++i ){
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ];
		if( m_signInversion[iElem][i] ){
			valU -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}else{
			valU += val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
			valV += val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
			valW += val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
		}
	}

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	std::complex<double> val = valU * std::complex<double>(jacobMat.mat12/detJacob, 0.0) 
							 + valV * std::complex<double>(jacobMat.mat22/detJacob, 0.0) 
							 + valW * std::complex<double>(jacobMat.mat32/detJacob, 0.0);

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Z component of magnetic field
std::complex<double> Forward3DTetraElement0thOrder::calcValueMagneticFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL." << std::endl;
	//	exit(1);
	//}
	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	//std::complex<double> valU(0.0, 0.0);
	//std::complex<double> valV(0.0, 0.0);
	//std::complex<double> valW(0.0, 0.0);

	//for( int i = 0; i < 6; ++i ){
	//	const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );
	//	const std::complex<double> val = m_solution[ m_IDsLocal2Global[iElem][i] ];
	//	if( m_signInversion[iElem][i] ){
	//		valU -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
	//		valV -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
	//		valW -= val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
	//	}else{
	//		valU += val * std::complex<double>( getShapeFuncRotatedReferenceCoordU(i) * length, 0.0 );
	//		valV += val * std::complex<double>( getShapeFuncRotatedReferenceCoordV(i) * length, 0.0 );
	//		valW += val * std::complex<double>( getShapeFuncRotatedReferenceCoordW(i) * length, 0.0 );
	//	}
	//}

	//Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//double detJacob(0.0);
	//calcJacobianMatrix( iElem, jacobMat, detJacob );

	//std::complex<double> val = valU * std::complex<double>(jacobMat.mat13/detJacob, 0.0) + valV * std::complex<double>(jacobMat.mat23/detJacob, 0.0) + valW * std::complex<double>(jacobMat.mat33/detJacob, 0.0);

	//const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	//const double factor = omega * CommonParameters::mu;

	//val /= std::complex<double>(0.0, factor);

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	return calcValueRotatedElectricFieldZDirection( iElem, uCoord, vCoord, wCoord ) / std::complex<double>( 0.0, omega * CommonParameters::mu );

}

// Calculate difference of voltage for brick element
std::complex<double> Forward3DTetraElement0thOrder::calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
	const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const{

	OutputFiles::m_logFile << "Error : calcVoltageDifference for brick element can't be used for tetrahedral element ." << std::endl;
	exit(1);

}

// Calculate difference of voltage
std::complex<double> Forward3DTetraElement0thOrder::calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
		const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const{

	std::complex<double> voltageDifference = std::complex<double>(0.0, 0.0);
	const double eps = 1.0e-6;

	for( int ielem = 0; ielem < nElem; ++ielem ){

		const int elemID = elememtsIncludingDipole[ielem];
		const int faceID = facesIncludingDipole[ielem];

		const CommonParameters::AreaCoords* const areaCoord[2] = { &areaCoordValStartPoint[ielem], &areaCoordValEndPoint[ielem] };

		int edgeIDLocal2D[2] = { -1, -1 };
		if( areaCoord[0]->coord0 < eps && areaCoord[1]->coord0 < eps ){// Check whether both start and end point locate on line 1
			edgeIDLocal2D[0] = 1;
			edgeIDLocal2D[1] = 1;
		}else if( areaCoord[0]->coord1 < eps && areaCoord[1]->coord1 < eps ){// Check whether both start and end point locate on line 2
			edgeIDLocal2D[0] = 2;
			edgeIDLocal2D[1] = 2;
		}else if( areaCoord[0]->coord2 < eps && areaCoord[1]->coord2 < eps ){// Check whether both start and end point locate on line 0
			edgeIDLocal2D[0] = 0;
			edgeIDLocal2D[1] = 0;
		}else{// Start and end point locate on another line or inside of triangle

			for( int i = 0; i < 2; ++i ){
				if( areaCoord[i]->coord0 < eps ){
					edgeIDLocal2D[i] = 1;
				}else if( areaCoord[i]->coord1 < eps ){
					edgeIDLocal2D[i] = 2;
				}else if( areaCoord[i]->coord2 < eps ){
					edgeIDLocal2D[i] = 0;
				}
			}

		}

#ifdef _DEBUG_WRITE	
		std::cout << "ielem : " << ielem << std::endl;
		std::cout << "elemID : " << elemID << std::endl;
		std::cout << "faceID : " << faceID << std::endl;
		std::cout << "areaCoord[0] : " << areaCoord[0]->coord0 << " " << areaCoord[0]->coord1 << " " << areaCoord[0]->coord2 << std::endl;
		std::cout << "areaCoord[1] : " << areaCoord[1]->coord0 << " " << areaCoord[1]->coord1 << " " << areaCoord[1]->coord2 << std::endl;
		std::cout << "edgeIDLocal2D : " << edgeIDLocal2D[0] << " " << edgeIDLocal2D[1] << std::endl;
#endif

		if( edgeIDLocal2D[0] == -1 || edgeIDLocal2D[1] == -1 ){// Either start or end point locate inside of triangle

#ifdef _DEBUG_WRITE	
			std::cout << "A" << std::endl;
#endif

			const double diffX = m_MeshDataTetraElement.calcXCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcXCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );
			const double diffY = m_MeshDataTetraElement.calcYCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcYCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );
			const double diffZ = m_MeshDataTetraElement.calcZCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcZCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );

#ifdef _DEBUG_WRITE	
			std::cout << "diffX diffY diffZ : " << diffX << " " << diffY << " " << diffZ << std::endl;
#endif

			std::complex<double> val(0.0,0.0);
			for( int i = 0; i < 2; ++i ){
				CommonParameters::VolumeCoords volumCoord = { 0.0, 0.0, 0.0, 0.0 };
				m_MeshDataTetraElement.calcVolumeCoordFromAreaCoord( faceID, *areaCoord[i], volumCoord );
				const std::complex<double> Ex = calcValueElectricFieldXDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3 );
				const std::complex<double> Ey = calcValueElectricFieldYDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3 );
				const std::complex<double> Ez = calcValueElectricFieldZDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3 );
				val -= ( Ex * diffX + Ey * diffY + Ez * diffZ );// Reverse sign because V = - E * l    
			}
			voltageDifference += 0.5 * val;

#ifdef _DEBUG_WRITE	
			std::cout << "voltageDifference : " << voltageDifference << std::endl;
#endif

			continue;// Go Next segment

		}else if( edgeIDLocal2D[0] == edgeIDLocal2D[1] ){// Start poin and end point locate on the same edge

#ifdef _DEBUG_WRITE	
			std::cout << "B" << std::endl;
#endif
			
			double ratio(0.0);
			const bool reverse = calcRatioAndReverseFlag( faceID, edgeIDLocal2D[0], *areaCoord[0], *areaCoord[1], ratio );

			const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D[0] );
			const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, edgeIDLocal3D ) * ( ratio ) ;

			if( m_signInversion[elemID][edgeIDLocal3D] != reverse ){
				voltageDifference -= m_solution[ m_IDsLocal2Global[elemID][edgeIDLocal3D] ] * std::complex<double>( - length, 0.0 );// Reverse sign because V = - E * l
			}else{
				voltageDifference += m_solution[ m_IDsLocal2Global[elemID][edgeIDLocal3D] ] * std::complex<double>( - length, 0.0 );// Reverse sign because V = - E * l
			}

#ifdef _DEBUG_WRITE	
			std::cout << "reverse : " << reverse << std::endl;
			std::cout << "ratio : " << ratio << std::endl;
			std::cout << "edgeIDLocal3D : " << edgeIDLocal3D << std::endl;
			std::cout << "length : " << length << std::endl;
			std::cout << "voltageDifference : " << voltageDifference << std::endl;
#endif

			continue;// Go Next segment

		}else{// Start poin and end point locate on different edges

#ifdef _DEBUG_WRITE	
			std::cout << "C" << std::endl;
#endif

			////int sharedNodeIDLocal2D(-1);
			////if( ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 1 ) || ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 0 ) ){
			////	sharedNodeIDLocal2D = 1;
			////}else if( ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 2 ) || ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 1 ) ){
			////	sharedNodeIDLocal2D = 2;
			////}else if( ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 0 ) || ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 2 ) ){
			////	sharedNodeIDLocal2D = 0;
			////}
			CommonParameters::AreaCoords areaCoordOfSharePoint = { 0.0, 0.0, 0.0 };
			if( ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 1 ) || ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 0 ) ){
				areaCoordOfSharePoint.coord1 = 1.0;
			}else if( ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 2 ) || ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 1 ) ){
				areaCoordOfSharePoint.coord2 = 1.0;
			}else if( ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 0 ) || ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 2 ) ){
				areaCoordOfSharePoint.coord0 = 1.0;
			}

			double ratio[2] = { 0.0, 0.0 };
			const bool reverse[2] = { calcRatioAndReverseFlag( faceID, edgeIDLocal2D[0], *areaCoord[0], areaCoordOfSharePoint, ratio[0] ),
										calcRatioAndReverseFlag( faceID, edgeIDLocal2D[1], areaCoordOfSharePoint, *areaCoord[1], ratio[1] ) };

#ifdef _DEBUG_WRITE	
			std::cout << "areaCoordOfSharePoint : " << areaCoordOfSharePoint.coord0 << " " << areaCoordOfSharePoint.coord1 << " " << areaCoordOfSharePoint.coord2 << std::endl;
#endif

			for( int i = 0; i < 2; ++i ){// Start point, End point
				
#ifdef _DEBUG_WRITE	
				std::cout << "i = " << i << std::endl;
#endif

				const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D[i] );
				const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, edgeIDLocal3D ) * ratio[i];

				if( m_signInversion[elemID][edgeIDLocal3D] != reverse[i] ){
					voltageDifference -= m_solution[ m_IDsLocal2Global[elemID][edgeIDLocal3D] ] * std::complex<double>( - length, 0.0 );// Reverse sign because V = - E * l
				}else{
					voltageDifference += m_solution[ m_IDsLocal2Global[elemID][edgeIDLocal3D] ] * std::complex<double>( - length, 0.0 );// Reverse sign because V = - E * l
				}

#ifdef _DEBUG_WRITE	
				std::cout << "reverse[i] : " << reverse[i] << std::endl;
				std::cout << "ratio : " << ratio[i] << std::endl;
				std::cout << "edgeIDLocal3D : " << edgeIDLocal3D << std::endl;
				std::cout << "length : " << length << std::endl;
				std::cout << "voltageDifference : " << voltageDifference << std::endl;
#endif

			}

			const double areaWithSign = m_MeshDataTetraElement.calcAreaWithSignFromAreaCoords( elemID, faceID, *areaCoord[0], areaCoordOfSharePoint, *areaCoord[1] );

			voltageDifference += calcValueRotatedElectricFieldNormal( elemID, faceID, 0.333333333333, 0.333333333333 ) * areaWithSign;

#ifdef _DEBUG_WRITE	
			std::cout << "areaWithSign calcValueRotatedElectricFieldZDirection : " << areaWithSign << " " << calcValueRotatedElectricFieldZDirection( elemID, 0.333333333333, 0.333333333333, 0.333333333333 ) << std::endl;
			std::cout << "areaWithSign calcValueRotatedElectricFieldNormal : " << areaWithSign << " " << calcValueRotatedElectricFieldNormal( elemID, faceID, 0.333333333333, 0.333333333333 ) << std::endl;
			std::cout << "voltageDifference : " << voltageDifference << std::endl;
#endif

			continue;// Go Next segment

		}

	}

	return voltageDifference;

}

// Calculate interpolator vector of X component of electric field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat11/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat12/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat13/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat11/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat12/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat13/detJacob, 0.0 ) * factor );
		}
	}

}

// Calculate interpolator vector of Y component of electric field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();

	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}
	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat21/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat22/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat23/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat21/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat22/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat23/detJacob, 0.0 ) * factor );
		}
	}

}

// Calculate interpolator vector of Z component of electric field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	Forward3D::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);

	calcJacobianMatrix( iElem, jacobMat, detJacob );
	calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat31/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat32/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat33/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordU( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat31/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordV( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat32/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncReferenceCoordW( uCoord, vCoord, wCoord, i ) * length * invJacobMat.mat33/detJacob, 0.0 ) * factor );
		}
	}

}

// Calculate interpolator vector of Z component of rotated electric field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat13/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat23/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat33/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat13/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat23/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat33/detJacob, 0.0 ) * factor );
		}
	}

}

// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, 
	const int irhs, const std::complex<double>& factor ){

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Factor is multiplied by the inverse of the inclination in the X direction
	calcInterpolatorVectorOfElectricFieldTangentialX(iElem, iFace, uCoord, vCoord, irhs, factor * std::complex<double>(1.0/dLengdX, 0.0) );

}

// Calculate interpolator vector of Y component of electric field only from the edges on the Earth's surface
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Factor is multiplied by the inverse of the inclination in the Y direction
	calcInterpolatorVectorOfElectricFieldTangentialY(iElem, iFace, uCoord, vCoord, irhs, factor * std::complex<double>(1.0/dLengdY, 0.0) );

}

// Calculate interpolator vector of tangential electric field directed to X from all edges
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal,
	const int irhs, const std::complex<double>& factor ){

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Calculate interpolator vector of X component of electric field
	calcInterpolatorVectorOfElectricFieldXDirection( iElem, xLocal, yLocal, zLocal, irhs, factor * std::complex<double>(dLengdX, 0.0) );

}

// Calculate interpolator vector of tangential electric field directed to Y from all edges
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, 
	const int irhs, const std::complex<double>& factor ){

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(iElem, iFace, dLengdX, dLengdY);

	// Calculate interpolator vector of Y component of electric field
	calcInterpolatorVectorOfElectricFieldYDirection( iElem, xLocal, yLocal, zLocal, irhs, factor * std::complex<double>(dLengdY, 0.0) );
}

// Calculate interpolator vector of tangential electric field directed to X
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){

	const int iPol = getPolarizationCurrent();

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 
	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length * jacobMat.mat22/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length * jacobMat.mat12/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length * jacobMat.mat22/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length * jacobMat.mat12/detJacob, 0.0 ) * factor );
		}
	}
	
}

// Calculate interpolator vector of tangential electric field directed to Y
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord,
	const int irhs, const std::complex<double>& factor ){

	const int iPol = getPolarizationCurrent();

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 
	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length * jacobMat.mat11/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length * jacobMat.mat21/detJacob, 0.0 ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   get2DShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length * jacobMat.mat11/detJacob, 0.0 ) * factor );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - get2DShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length * jacobMat.mat21/detJacob, 0.0 ) * factor );
		}
	}

}

// Calculate interpolator vector of normal component of rotated electric field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfRotatedElectricFieldNormal( const int iElem, const int iFace, const double uCoord, const double vCoord,
																					   const int irhs,  const std::complex<double>& factor ){

	const int iPol = getPolarizationCurrent();

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 }; 
	double detJacob(0.0);
	calc2DJacobianMatrix( iElem, iFace, jacobMat, detJacob );
	
	for( int i = 0; i < 3; ++i ){
		const int iEdge = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( iFace, i );
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, iEdge );
		if( m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, i ) > m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( iElem, iFace, ( i + 1 ) % 3 ) ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - get2DShapeFuncRotated() * length / detJacob ) * factor );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   get2DShapeFuncRotated() * length / detJacob ) * factor );
		}
	}
	
}

// Calculate interpolator vector of X component of magnetic field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat11/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat21/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat31/detJacob, 0.0 ) * constant );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat11/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat21/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat31/detJacob, 0.0 ) * constant );
		}
	}

}

// Calculate interpolator vector of Y component of magnetic field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
	const int irhs,	const std::complex<double>& factor ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );

	Forward3D::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	calcJacobianMatrix( iElem, jacobMat, detJacob );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	for( int i = 0; i < 6; ++i ){

		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}

		const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

		if( m_signInversion[iElem][i] ){
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat12/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat22/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat32/detJacob, 0.0 ) * constant );
		}else{
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat12/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat22/detJacob, 0.0 ) * constant );
			addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat32/detJacob, 0.0 ) * constant );
		}
	}

}

// Calculate interpolator vector of Z component of magnetic field
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double uCoord, const double vCoord, const double wCoord,
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

	//Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//double detJacob(0.0);
	//calcJacobianMatrix( iElem, jacobMat, detJacob );

	//const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	//const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	//for( int i = 0; i < 6; ++i ){

	//	const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
	//	if( irow < 0 ){
	//		continue;
	//	}

	//	const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( iElem, i );

	//	if( m_signInversion[iElem][i] ){
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat13/detJacob, 0.0 ) * constant );
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat23/detJacob, 0.0 ) * constant );
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>( - getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat33/detJacob, 0.0 ) * constant );
	//	}else{
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordU(i) * length * jacobMat.mat13/detJacob, 0.0 ) * constant );
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordV(i) * length * jacobMat.mat23/detJacob, 0.0 ) * constant );
	//		addValuesToRhsVectors( irow, irhs, std::complex<double>(   getShapeFuncRotatedReferenceCoordW(i) * length * jacobMat.mat33/detJacob, 0.0 ) * constant );
	//	}
	//}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> constant = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	calcInterpolatorVectorOfRotatedElectricFieldZDirection( iElem, uCoord, vCoord, wCoord, irhs, constant );

}

// Calculate interpolator vector of difference of voltage
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs ){
	OutputFiles::m_logFile << "Error : calcInterpolatorVectorOfVoltageDifference for brick element can't be used for tetrahedral element ." << std::endl;
	exit(1);
}

// Calculate interpolator vector of difference of voltage
void Forward3DTetraElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
	const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs ){

	//if( m_IDsLocal2Global == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global is NULL." << std::endl;
	//	exit(1);
	//}

	const int iPol = getPolarizationCurrent();
	//if( m_IDsGlobal2AfterDegenerated[iPol] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsGlobal2AfterDegenerated[iPol] is NULL." << std::endl;
	//	exit(1);
	//}

	assert( m_IDsLocal2Global != NULL );
	assert( m_IDsGlobal2AfterDegenerated[iPol] != NULL );
	const double eps = 1.0e-6;

	for( int ielem = 0; ielem < nElem; ++ielem ){

		const int elemID = elememtsIncludingDipole[ielem];
		const int faceID = facesIncludingDipole[ielem];

		const CommonParameters::AreaCoords* const areaCoord[2] = { &areaCoordValStartPoint[ielem], &areaCoordValEndPoint[ielem] };

		int edgeIDLocal2D[2] = { -1, -1 };
		if( areaCoord[0]->coord0 < eps && areaCoord[1]->coord0 < eps ){// Check whether both start and end point locate on line 1
			edgeIDLocal2D[0] = 1;
			edgeIDLocal2D[1] = 1;
		}else if( areaCoord[0]->coord1 < eps && areaCoord[1]->coord1 < eps ){// Check whether both start and end point locate on line 2
			edgeIDLocal2D[0] = 2;
			edgeIDLocal2D[1] = 2;
		}else if( areaCoord[0]->coord2 < eps && areaCoord[1]->coord2 < eps ){// Check whether both start and end point locate on line 0
			edgeIDLocal2D[0] = 0;
			edgeIDLocal2D[1] = 0;
		}else{// Start and end point locate on another line or inside of triangle

			for( int i = 0; i < 2; ++i ){
				if( areaCoord[i]->coord0 < eps ){
					edgeIDLocal2D[i] = 1;
				}else if( areaCoord[i]->coord1 < eps ){
					edgeIDLocal2D[i] = 2;
				}else if( areaCoord[i]->coord2 < eps ){
					edgeIDLocal2D[i] = 0;
				}
			}

		}

		//-------------------------------------------------------------------------------------------------------------------------------------
		if( edgeIDLocal2D[0] == -1 || edgeIDLocal2D[1] == -1 ){// Either start or end point locate inside of triangle

			const double diffX = m_MeshDataTetraElement.calcXCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcXCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );
			const double diffY = m_MeshDataTetraElement.calcYCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcYCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );
			const double diffZ = m_MeshDataTetraElement.calcZCoordOfPointOnFace( elemID, faceID, *areaCoord[1] ) - m_MeshDataTetraElement.calcZCoordOfPointOnFace( elemID, faceID, *areaCoord[0] );

			for( int i = 0; i < 2; ++i ){
				CommonParameters::VolumeCoords volumCoord = { 0.0, 0.0, 0.0, 0.0 };
				m_MeshDataTetraElement.calcVolumeCoordFromAreaCoord( faceID, *areaCoord[i], volumCoord );
				calcInterpolatorVectorOfElectricFieldXDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3, irhs, std::complex<double>(-0.5 * diffX, 0.0) );
				calcInterpolatorVectorOfElectricFieldYDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3, irhs, std::complex<double>(-0.5 * diffY, 0.0) );
				calcInterpolatorVectorOfElectricFieldZDirection( elemID, volumCoord.coord1, volumCoord.coord2, volumCoord.coord3, irhs, std::complex<double>(-0.5 * diffZ, 0.0) );
			}

			continue;// Go Next segment

		//-------------------------------------------------------------------------------------------------------------------------------------
		}else if( edgeIDLocal2D[0] == edgeIDLocal2D[1] ){// Start poin and end point locate on the same edge
		
			double ratio(0.0);
			const bool reverse = calcRatioAndReverseFlag( faceID, edgeIDLocal2D[0], *areaCoord[0], *areaCoord[1], ratio );

			const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D[0] );
			const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, edgeIDLocal3D ) * ( ratio ) ;

			const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][edgeIDLocal3D] ];
			if( irow < 0 ){
				continue;
			}

			if( m_signInversion[elemID][edgeIDLocal3D] != reverse ){
				addValuesToRhsVectors( irow, irhs, std::complex<double>(  length, 0.0 ) );// Reverse sign because V = - E * l
			}else{
				addValuesToRhsVectors( irow, irhs, std::complex<double>( -length, 0.0 ) );// Reverse sign because V = - E * l
			}

			continue;// Go Next segment

		//-------------------------------------------------------------------------------------------------------------------------------------
		}else{// Start poin and end point locate on different edges

			CommonParameters::AreaCoords areaCoordOfSharePoint = { 0.0, 0.0, 0.0 };
			if( ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 1 ) || ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 0 ) ){
				areaCoordOfSharePoint.coord1 = 1.0;
			}else if( ( edgeIDLocal2D[0] == 1 && edgeIDLocal2D[1] == 2 ) || ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 1 ) ){
				areaCoordOfSharePoint.coord2 = 1.0;
			}else if( ( edgeIDLocal2D[0] == 2 && edgeIDLocal2D[1] == 0 ) || ( edgeIDLocal2D[0] == 0 && edgeIDLocal2D[1] == 2 ) ){
				areaCoordOfSharePoint.coord0 = 1.0;
			}

			double ratio[2] = { 0.0, 0.0 };
			const bool reverse[2] = { calcRatioAndReverseFlag( faceID, edgeIDLocal2D[0], *areaCoord[0], areaCoordOfSharePoint, ratio[0] ),
										calcRatioAndReverseFlag( faceID, edgeIDLocal2D[1], areaCoordOfSharePoint, *areaCoord[1], ratio[1] ) };

			for( int i = 0; i < 2; ++i ){// Start point, End point
				
				const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D[i] );
				const double length = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, edgeIDLocal3D ) * ratio[i];

				const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][edgeIDLocal3D] ];
				if( irow < 0 ){
					continue;
				}

				if( m_signInversion[elemID][edgeIDLocal3D] != reverse[i] ){
					addValuesToRhsVectors( irow, irhs, std::complex<double>(  length, 0.0 ) );// Reverse sign because V = - E * l
				}else{
					addValuesToRhsVectors( irow, irhs, std::complex<double>( -length, 0.0 ) );// Reverse sign because V = - E * l
				}

			}

			const double areaWithSign = m_MeshDataTetraElement.calcAreaWithSignFromAreaCoords( elemID, faceID, *areaCoord[0], areaCoordOfSharePoint, *areaCoord[1] );
			calcInterpolatorVectorOfRotatedElectricFieldNormal( elemID, faceID, 0.333333333333, 0.333333333333, irhs, areaWithSign );

			continue;// Go Next segment

		}

	}


}

// Set non-zero strucuture of matrix for forward calculation
void Forward3DTetraElement0thOrder::setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix ){

	const ResistivityBlock* const ptrResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();

	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

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

		for( iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){

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

#ifdef _ANISOTOROPY
// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward3DTetraElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ){

	if( (AnalysisControl::getInstance() )->isAnisotropyConsidered() ){
		setNonZeroValuesAnisotropy( matrix );
	}else{
		setNonZeroValuesIsotropy( matrix );
	}

}
	// Set non-zero values of matrix and right-hande side vector for isotropic medium
void Forward3DTetraElement0thOrder::setNonZeroValuesIsotropy( ComplexSparseSquareSymmetricMatrix& matrix ){

	assert( !(AnalysisControl::getInstance() )->isAnisotropyConsidered() );
#else
// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward3DTetraElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ){
#endif

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

	const int rowIndex[6] = { 0, 6, 11, 15, 18, 20 };

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	int iElem(0);
	int elemID(0);
	//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
	//Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//double detJacob(0.0);
	//double divDetJacob(0.0);
	//double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
	//int i(0);
	//int ip(0);
	//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<
	double sigma(0.0);
	double factor1(0.0);
	std::complex<double> factor2(0.0,0.0);
	double eMat[21];
	double fMat[21];
	int iEdge1(0);
	int iEdge2(0);
	int row(0);
	int col(0);
	double integral1(0.0);
	double integral2(0.0);
	std::complex<double> val(0.0,0.0);
	int loc(0);
	//------------------------------------------------------------ Start of parallel region >>>>>
#ifdef _USE_OMP
#pragma omp parallel for default(shared) \
		private( iElem, elemID, sigma, factor1, factor2, eMat, fMat, \
			iEdge1, iEdge2, row, col, integral1, integral2, val, loc )
#endif
	for( iElem = 0; iElem < nElem; ++iElem ){

		elemID = iElem;

		//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
		////--- Calculate Jacobian
		//calcJacobianMatrix( elemID, jacobMat, detJacob );
		//divDetJacob = 1.0 / detJacob;
		////--- Calculate inverse of Jacobian matrix multiplied by determinant
		//calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

		//for( i = 0; i < 6; ++i ){
		//	length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, i );
		//}
		//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<

		//--- Calculate omega * mu * sigma
		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

		factor1 = 1.0;
		factor2 = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form

		//--- Calculate integrals of which element matrix consisits
		calcIntegrals( elemID, eMat, fMat );

		for( iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){

				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}
				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use >>>>>
				//integral1 = 0.0;
				//integral2 = 0.0;
				//for( ip = 0; ip < m_numIntegralPoints; ++ip ){
				//	integral1 += (   ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31 )
				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32 )
				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33 ) ) * m_weights[ip];
				//	integral2 += ( (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat11
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat12
				//				    + getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat13 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat11
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat12
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat13 )
				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat21
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat22
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat23 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat21
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat22
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat23 )
				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat31
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat32
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat33 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat31
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat32
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat33 ) ) * m_weights[ip];
				//}

				//if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
				//	integral1 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//	integral2 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//}else{
				//	integral1 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//	integral2 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//}
				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use <<<<<

				if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
					integral1 = -eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
					integral2 = -fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
				}else{
					integral1 =  eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
					integral2 =  fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
				}

				val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form

				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
#ifdef _USE_OMP
					#pragma omp critical (addToRhs)
#endif
					{
						matrix.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
					}

				}else if( col >= row ){// Store only upper triangle part
					loc = matrix.checkAndGetLocationNonZeroValue( row, col );
#ifdef _USE_OMP
					#pragma omp critical (addToMatrix)
#endif
					{
						matrix.addNonZeroValuesWithoutSearchingLocation( loc, val );// Add to matrix
					}

				}

			}// iEdge2

		}// iEdge1		
	}

}

#ifdef _ANISOTOROPY
// Set non-zero values of matrix and right-hande side vector for forward calculation for anisotropic medium
void Forward3DTetraElement0thOrder::setNonZeroValuesAnisotropy( ComplexSparseSquareSymmetricMatrix& matrix ){

	assert( (AnalysisControl::getInstance() )->isAnisotropyConsidered() );

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	const int nElem = m_MeshDataTetraElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = iElem;

		//--- Calculate Jacobian
		Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double detJacob(0.0);
		calcJacobianMatrix( elemID, jacobMat, detJacob );
		const double divDetJacob = 1.0 / detJacob;
		//--- Calculate inverse of Jacobian matrix multiplied by determinant
		Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

		double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
		for( int i = 0; i < 6; ++i ){
			length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, i );
		}

		//--- Calculate omega * mu * sigma
		const double factor1 = 1.0;
		const std::complex<double> factor2 = std::complex<double>( 0.0, omega * CommonParameters::mu );// exp(-i*omega*t) form

		//--- Calculate coefficents of conductivity tensor
		//CommonParameters::Vector3D matX = { 0.0, 0.0, 0.0 };
		//CommonParameters::Vector3D matY = { 0.0, 0.0, 0.0 };
		//CommonParameters::Vector3D matZ = { 0.0, 0.0, 0.0 };
		//pResistivityBlock->calcAisotropicConductivityTensor( pResistivityBlock->getBlockIDFromElemID(iElem), matX, matY, matZ );
		double sigmaTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcAisotropicConductivityTensor( pResistivityBlock->getBlockIDFromElemID(iElem), sigmaTensor );

		for( int iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){
			const int row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){
				const int col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for( int ip = 0; ip < m_numIntegralPoints; ++ip ){
					integral1 += (  ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11
									+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21
									+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31 )
								  * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11
								    + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21  
									+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31 )
								  + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12 
								    + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22  
									+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32 )
								  * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12 
								    + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22  
									+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32 )
								  + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13 
								    + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23  
									+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33 )
								  * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13 
								    + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23  
									+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33 ) ) * m_weights[ip];
					const double Nx = getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat11
									+ getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat12
									+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat13;
					const double Ny = getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat21
									+ getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat22
									+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat23;
					const double Nz = getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat31
									+ getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat32
									+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat33; 
					integral2 += ( (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat11
								    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat12
								    + getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat13 )
								 * ( Nx * sigmaTensor[0][0] + Ny * sigmaTensor[0][1] + Nz * sigmaTensor[0][2] )
								 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat21
								    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat22
									+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat23 )
								 * ( Nx * sigmaTensor[1][0] + Ny * sigmaTensor[1][1] + Nz * sigmaTensor[1][2] )
								 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat31
								    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat32
									+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat33 )
								 * ( Nx * sigmaTensor[2][0] + Ny * sigmaTensor[2][1] + Nz * sigmaTensor[2][2] ) ) * m_weights[ip];
				}

				if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
					integral1 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
					integral2 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
				}else{
					integral1 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
					integral2 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
				}

				const std::complex<double> val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form

				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					matrix.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
				}else if( col >= row ){// Store only upper triangle part
					const int loc = matrix.checkAndGetLocationNonZeroValue( row, col );
					matrix.addNonZeroValuesWithoutSearchingLocation( loc, val );// Add to matrix
				}

			}// iEdge2
		}// iEdge1		
	}// iElem

}
#endif

//----- DO NOT DELETE FOR FUTURE USE >>>>>
//// Set non-zero strucuture of matrix for calculating derivatives
//void Forward3DTetraElement0thOrder::setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID, std::set<int>& nonZeroRowsAndCols ){
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
//		for( iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){
//
//			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
//			if( row < 0 ){
//				continue;
//			}
//
//			for( iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){
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
//void Forward3DTetraElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix, const int blkID ){
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
//	const int rowIndex[6] = { 0, 6, 11, 15, 18, 20 };
//
//	//------------------------------------------
//	//--- Components due to stiffness matrix ---
//	//------------------------------------------
//	int iElem(0);
//	int elemID(0);
//	//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
//	//Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//	//Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//	//double detJacob(0.0);
//	//double divDetJacob(0.0);
//	//double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
//	//int i(0);
//	//int ip(0);
//	//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<
//	double sigma(0.0);
//	double factor1(0.0);
//	std::complex<double> factor2(0.0,0.0);
//	double eMat[21];
//	double fMat[21];
//	int iEdge1(0);
//	int iEdge2(0);
//	int row(0);
//	int col(0);
//	double integral1(0.0);
//	double integral2(0.0);
//	std::complex<double> val(0.0,0.0);
//	int loc(0);
//	//------------------------------------------------------------ Start of parallel region >>>>>
////#ifdef _USE_OMP
////	#pragma omp parallel for default(shared) \
////		private( iElem, elemID, sigma, factor1, factor2, eMat, fMat, \
////			iEdge1, iEdge2, row, col, integral1, integral2, val, loc )
////#endif
//	for( iElem = 0; iElem < nElem; ++iElem ){
//
//		// [Attention] : You must use elemID instead of iElem from this line
//		elemID = mdl2Elem[iElem].first;
//
//#ifdef _DEBUG_WRITE
//		std::cout << "blkID iElem elemID = " << blkID << " " << iElem << " " << elemID << std::endl;
//#endif
//
//		//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
//		////--- Calculate Jacobian
//		//calcJacobianMatrix( elemID, jacobMat, detJacob );
//		//divDetJacob = 1.0 / detJacob;
//		////--- Calculate inverse of Jacobian matrix multiplied by determinant
//		//calcInverseOfJacobianMatrix( jacobMat, invJacobMat );
//
//		//for( i = 0; i < 6; ++i ){
//		//	length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, i );
//		//}
//		//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<
//
//		//--- Calculate omega * mu * sigma
//		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
//
//		factor1 = 0.0;
//		factor2 = std::complex<double>( 0.0, - omega * CommonParameters::mu * sigma * ln10 * mdl2Elem[iElem].second );// exp(-i*omega*t) form
//
//		//--- Calculate integrals of which element matrix consisits
//		calcIntegrals( elemID, eMat, fMat );
//
//		for( iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){
//
//			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
//			if( row < 0 ){
//				continue;
//			}
//
//			for( iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){
//
//				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
//				if( col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
//					continue;
//				}
//
//				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use >>>>>
//				//integral1 = 0.0;
//				//integral2 = 0.0;
//				//for( ip = 0; ip < m_numIntegralPoints; ++ip ){
//				//	integral1 += (   ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31 )
//				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31 )
//				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32 )
//				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32 )
//				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33 )
//				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33 ) ) * m_weights[ip];
//				//	integral2 += ( (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat11
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat12
//				//				    + getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat13 )
//				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat11
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat12
//				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat13 )
//				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat21
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat22
//				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat23 )
//				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat21
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat22
//				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat23 )
//				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat31
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat32
//				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat33 )
//				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat31
//				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat32
//				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat33 ) ) * m_weights[ip];
//				//}
//
//				//if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
//				//	integral1 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
//				//	integral2 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
//				//}else{
//				//	integral1 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
//				//	integral2 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
//				//}
//				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use <<<<<
//
//				if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
//					integral1 = -eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
//					integral2 = -fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
//				}else{
//					integral1 =  eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
//					integral2 =  fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
//				}
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
void Forward3DTetraElement0thOrder::calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows ){

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

	const int rowIndex[6] = { 0, 6, 11, 15, 18, 20 };

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	int iElem(0);
	int elemID(0);
	//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
	//Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//double detJacob(0.0);
	//double divDetJacob(0.0);
	//double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
	//int i(0);
	//int ip(0);
	//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<
	double sigma(0.0);
	double factor1(0.0);
	std::complex<double> factor2(0.0,0.0);
	double eMat[21];
	double fMat[21];
	int iEdge1(0);
	int iEdge2(0);
	int row(0);
	int col(0);
	double integral1(0.0);
	double integral2(0.0);
	std::complex<double> val(0.0,0.0);
//#ifdef _USE_OMP
//	#pragma omp parallel for default(shared) \
//		private( iElem, elemID, sigma, factor1, factor2, eMat, fMat, \
//			iEdge1, iEdge2, row, col, integral1, integral2, val, loc )
//#endif
	for( iElem = 0; iElem < nElem; ++iElem ){

		// [Attention] : You must use elemID instead of iElem from this line
		elemID = mdl2Elem[iElem].first;

#ifdef _DEBUG_WRITE
		std::cout << "blkID iElem elemID = " << blkID << " " << iElem << " " << elemID << std::endl;
#endif

		//----- modified by Y.Usui 2013.12.15 Do not delete for future use >>>>>
		////--- Calculate Jacobian
		//calcJacobianMatrix( elemID, jacobMat, detJacob );
		//divDetJacob = 1.0 / detJacob;
		////--- Calculate inverse of Jacobian matrix multiplied by determinant
		//calcInverseOfJacobianMatrix( jacobMat, invJacobMat );

		//for( i = 0; i < 6; ++i ){
		//	length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, i );
		//}
		//----- modified by Y.Usui 2013.12.15 Do not delete for future use <<<<<

		//--- Calculate omega * mu * sigma
		sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

		factor1 = 0.0;
		factor2 = std::complex<double>( 0.0, - omega * CommonParameters::mu * sigma * ln10 * mdl2Elem[iElem].second );// exp(-i*omega*t) form

		//--- Calculate integrals of which element matrix consisits
		calcIntegrals( elemID, eMat, fMat );

		for( iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){

			row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}

			for( iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){

				col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}

				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use >>>>>
				//integral1 = 0.0;
				//integral2 = 0.0;
				//for( ip = 0; ip < m_numIntegralPoints; ++ip ){
				//	integral1 += (   ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31 )
				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32 )
				//				   + ( getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33 )
				//				   * ( getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13 + getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23  + getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33 ) ) * m_weights[ip];
				//	integral2 += ( (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat11
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat12
				//				    + getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat13 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat11
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat12
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat13 )
				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat21
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat22
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat23 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat21
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat22
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat23 )
				//				 + (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat31
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat32
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1 ) * invJacobMat.mat33 )
				//				 * (  getShapeFuncReferenceCoordU( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat31
				//				    + getShapeFuncReferenceCoordV( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat32
				//					+ getShapeFuncReferenceCoordW( m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2 ) * invJacobMat.mat33 ) ) * m_weights[ip];
				//}

				//if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
				//	integral1 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//	integral2 *= - ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//}else{
				//	integral1 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//	integral2 *= ( divDetJacob * length[iEdge1] * length[iEdge2] );
				//}
				//----- modified by Y.Usui 2013.12.15 : Do not delete for future use <<<<<

				if( m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2] ){
					integral1 = -eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
					integral2 = -fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
				}else{
					integral1 =  eMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
					integral2 =  fMat[ iEdge1 < iEdge2 ? rowIndex[iEdge1] + iEdge2 - iEdge1 : rowIndex[iEdge2] + iEdge1 - iEdge2 ];
				}

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
void Forward3DTetraElement0thOrder::callInputMeshData(){

	m_MeshDataTetraElement.inputMeshData();

}

// Get pointer to the class MeshData
const MeshData* Forward3DTetraElement0thOrder::getPointerToMeshData() const{

	return static_cast<const MeshData*>( &m_MeshDataTetraElement ) ;

}

// Get pointer to the class MeshDataTetraElement
const MeshDataTetraElement* Forward3DTetraElement0thOrder::getPointerToMeshDataTetraElement() const{

	return &m_MeshDataTetraElement;

}

// Calculate array converting local IDs to global ones
void Forward3DTetraElement0thOrder::calcArrayConvertLocalID2Global(){
	
	typedef std::pair< int, int > NodeIDPair;
	typedef std::pair< int, int > ElemAndEdgeID;
	typedef std::pair< NodeIDPair, ElemAndEdgeID > InputedDataType;

	std::vector< InputedDataType > workVector;

	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 6; ++iEdge ){
			const int nodeID0 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 0 );
			const int nodeID1 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 1 ); 
			const NodeIDPair nodePairCur = nodeID1 > nodeID0 ? std::make_pair( nodeID0, nodeID1 ) : std::make_pair( nodeID1, nodeID0 );
			const ElemAndEdgeID elemAndEdgeIDCur = std::make_pair( iElem, iEdge );
			workVector.push_back( InputedDataType( std::pair< NodeIDPair, ElemAndEdgeID >(nodePairCur, elemAndEdgeIDCur) ) );
		}
	}

	std::sort( workVector.begin(), workVector.end() );

#ifdef _DEBUG_WRITE
	for( std::vector<InputedDataType>::iterator itr = workVector.begin(); itr != workVector.end(); ++itr ){
		std::cout << "workVector : " << (*itr).first.first << " " << (*itr).first.second << " " << (*itr).second.first << " " << (*itr).second.second << std::endl;
	}
#endif

	// Allocate memory to m_IDsLocal2Global ---
	if( m_IDsLocal2Global != NULL ){
		const int nElem = sizeof( m_IDsLocal2Global ) / sizeof( m_IDsLocal2Global[0] );
		for( int iElem = 0; iElem < nElem; ++iElem ){
			delete [] m_IDsLocal2Global[iElem];
		}
		delete [] m_IDsLocal2Global;
		m_IDsLocal2Global = NULL;
	}

	m_IDsLocal2Global = new int*[nElem]; 
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_IDsLocal2Global[iElem] = new int[6];
	}
	//-----------------------------------------

	// Allocate memory to m_signInversion ---
	if( m_signInversion != NULL ){
		const int num = sizeof( m_signInversion ) / sizeof( m_signInversion[0] );
		for( int i = 0; i < num; ++i ){
			delete [] m_signInversion[i];
			m_signInversion[i] = NULL;
		}
		delete [] m_signInversion;
		m_signInversion = NULL;
	}

	m_signInversion = new bool*[nElem];
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_signInversion[iElem] = new bool[6];
	}
	//----------------------------------------

	std::pair<int, int> nodePairPre;
	int edgeIDGlobal(-1);

	std::vector<InputedDataType>::iterator itrEnd = workVector.end();
	for( std::vector<InputedDataType>::iterator itr = workVector.begin(); itr != itrEnd; ++itr ){

		const int elemIDGlobal = (*itr).second.first;
		const int edgeIDLocal = (*itr).second.second;

		const int nodePairCur0 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 0 );
		const int nodePairCur1 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 1 );

//#ifdef _DEBUG_WRITE
//		std::cout << "nodePairCur0 nodePairCur1 nodePairPre.first nodePairPre.second :  " << nodePairCur0 << " " << nodePairCur1 << " " << nodePairPre.first << " " << nodePairPre.second << std::endl;
//#endif

		if( ( nodePairCur0 == nodePairPre.first && nodePairCur1 == nodePairPre.second ) ||
			( nodePairCur0 == nodePairPre.second && nodePairCur1 == nodePairPre.first ) ){// Same edge with privious one

		}else{// New edge

			++edgeIDGlobal;

			nodePairPre.first  = nodePairCur0;
			nodePairPre.second = nodePairCur1;

		}

		m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] = edgeIDGlobal;

		m_signInversion[elemIDGlobal][edgeIDLocal] = nodePairCur1 > nodePairCur0 ? false : true; 

#ifdef _DEBUG_WRITE
		std::cout << "elemIDGlobal edgeIDLocal nodePairCur0 nodePairCur1 m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] m_signInversion :  " << elemIDGlobal << " " << edgeIDLocal << " " << nodePairCur0 << " " << nodePairCur1 << " " << m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] << " " << m_signInversion[elemIDGlobal][edgeIDLocal] << std::endl;
#endif

	}

	m_numOfEquation = edgeIDGlobal + 1;

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	std::cout << "m_numOfEquation = " << m_numOfEquation << std::endl;
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	std::cout << "iElem m_IDsLocal2Global[iElem][iEdge] : " << iElem;
	//	for( int iEdge = 0; iEdge < 6; ++iEdge ){
	//		std::cout << " " << m_IDsLocal2Global[iElem][iEdge];
	//	}
	//	std::cout << std::endl;
	//}
#endif
	//----- debug <<<<<

	m_hasSetIDsLocal2Global = true;

}

// Calculate array converting global IDs to the ones after degeneration
void Forward3DTetraElement0thOrder::calcArrayConvertIDsGlobal2AfterDegenerated(){

	if( m_hasSetIDsLocal2Global == false ){
		OutputFiles::m_logFile << "Error : Array converting local edge IDs to global ones has not been set yet." << std::endl;
		exit(1);
	}

	const int iPol = getPolarizationCurrent();

	if( m_IDsGlobal2AfterDegenerated[iPol] != NULL ){
		delete [] m_IDsGlobal2AfterDegenerated[iPol];
		m_IDsGlobal2AfterDegenerated[iPol] = NULL;
	}
	m_IDsGlobal2AfterDegenerated[iPol] = new int[m_numOfEquation];
	for( int i = 0; i < m_numOfEquation; ++i ){
		m_IDsGlobal2AfterDegenerated[iPol][i] = 0;
	}

	int planesWithNonzeroValues[3] = {  -1, -1, -1 };
	int planesWithZeroValues[3] = {  -1, -1, -1 };
	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex Polarization
		planesWithNonzeroValues[0] = MeshData::ZXMinus;
		planesWithNonzeroValues[1] = MeshData::ZXPlus;
		planesWithNonzeroValues[2] = MeshData::XYMinus;
		planesWithZeroValues[0] = MeshData::YZMinus;
		planesWithZeroValues[1] = MeshData::YZPlus;
		planesWithZeroValues[2] = MeshData::XYPlus;
	}else{// Ey Polarization
		planesWithNonzeroValues[0] = MeshData::YZMinus;
		planesWithNonzeroValues[1] = MeshData::YZPlus;
		planesWithNonzeroValues[2] = MeshData::XYMinus;
		planesWithZeroValues[0] = MeshData::ZXMinus;
		planesWithZeroValues[1] = MeshData::ZXPlus;
		planesWithZeroValues[2] = MeshData::XYPlus;
	}

	// Planes on which dirichlet boundary condition with non-zero values is specified
	for( int iPlane = 0; iPlane < 3; ++iPlane ){

		const int planeID = planesWithNonzeroValues[iPlane];

		const int nElemOnPlane = m_MeshDataTetraElement.getNumElemOnBoundaryPlanes( planeID );

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){

			const int elemIDGlobal = m_MeshDataTetraElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataTetraElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );

			for( int iEdge = 0; iEdge < 3; ++iEdge ){
				const int edgeIDLocal = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] ] = Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE;
			}

		}

	}

	// Planes on which dirichlet boundary condition with zero values is specified
	for( int iPlane = 0; iPlane < 3; ++iPlane ){

		const int planeID = planesWithZeroValues[iPlane];

		const int nElemOnPlane = m_MeshDataTetraElement.getNumElemOnBoundaryPlanes( planeID );

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){

			const int elemIDGlobal = m_MeshDataTetraElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataTetraElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );

			for( int iEdge = 0; iEdge < 3; ++iEdge ){
				const int edgeIDLocal = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] ] = Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE;

//#ifdef _DEBUG_WRITE
//				std::cout << "[zero value] planeID elemID faceID edgeIDLocal edgeID : " << planeID << " " << elemIDGlobal << " " << faceIDLocal << " " << edgeIDLocal << " " << m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] << std::endl;
//#endif

			}

		}

	}

	int icount(-1);
	for( int i = 0; i < m_numOfEquation; ++i ){
		if( m_IDsGlobal2AfterDegenerated[iPol][i] < 0 ){
			continue;
		}
		m_IDsGlobal2AfterDegenerated[iPol][i] = ++icount;
	}

	m_numOfEquationDegenerated = icount + 1;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquation; ++i ){
		std::cout << "i m_IDsGlobal2AfterDegenerated[iPol][i] : " << i << " " << m_IDsGlobal2AfterDegenerated[iPol][i] << std::endl;
	}
	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 6; ++iEdge ){
			std::cout << "iElem iEdge m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ] : "  << iElem << " "  << iEdge << " " << m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ] << std::endl;
		}
	}

	std::cout << "m_numOfEquationDegenerated = " << m_numOfEquationDegenerated << std::endl;
#endif

	m_hasIDsGlobal2AfterDegenerated[iPol] = true;


}

// Calculate array converting global node IDs non-zero electric field values specified to the nodes
void Forward3DTetraElement0thOrder::calcArrayConvertIDGlobal2NonZeroValues(){

	if( m_hasSetIDsLocal2Global == false ){
		OutputFiles::m_logFile << "Error : Array converting local edge IDs to global ones has not been set yet." << std::endl;
		exit(1);
	}

	const int iPol = getPolarizationCurrent();

	if( m_hasIDsGlobal2AfterDegenerated[iPol] == false ){
		OutputFiles::m_logFile << "Error : Array converting global IDs to the ones after degeneration has not been set yet." << std::endl;
		exit(1);
	}

	//----------------------------------------------------------------------------
	//--- Calculating EM field of the boundary planes with 2D forward analysis ---
	//----------------------------------------------------------------------------
	const double freq = getFrequencyCurrent();
	OutputFiles::m_logFile << "# Calculating EM field of the boundary planes with 2D forward analysis.        " << std::endl;
	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex polarization

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DTriangleElement[MeshData::ZXMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataTetraElement );

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DTriangleElement[MeshData::ZXPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataTetraElement );

	}else{// Ey polarization

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DTriangleElement[MeshData::YZMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataTetraElement );

		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DTriangleElement[MeshData::YZPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataTetraElement );

	}
	OutputFiles::m_logFile << "#------------------------------------------------------------------------------" << std::endl;
	//----------------------------------------------------------------------------

	// Initialize array converting global edge IDs non-zero electric field values specified to the edges
	m_globalID2NonZeroValues.clear();

	int sidePlanes[2] = { -1, -1 };
	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex Polarization
		sidePlanes[0] = MeshData::ZXMinus;
		sidePlanes[1] = MeshData::ZXPlus;
	}else{// Ey Polarization
		sidePlanes[0] = MeshData::YZMinus;
		sidePlanes[1] = MeshData::YZPlus;
	}

	// Side planes on which dirichlet boundary condition with non-zero values is specified
	for( int iPlane = 0; iPlane < 2; ++iPlane ){

		const int planeID = sidePlanes[iPlane];

		int convertEdgeIDs[3] = { 0, 1, 2 };
		if( planeID == MeshData::ZXMinus || planeID == MeshData::YZMinus ){// If the plane is minus side
			convertEdgeIDs[0] = 2;
			convertEdgeIDs[2] = 0;
		}

		const int nElemOnPlane = m_MeshDataTetraElement.getNumElemOnBoundaryPlanes( planeID );

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){

			const int elemIDGlobal = m_MeshDataTetraElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataTetraElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );

			for( int iEdge = 0; iEdge < 3; ++iEdge ){
				const int edgeIDLocal = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );

				std::complex<double> val = m_Fwd2DTriangleElement[planeID]->getSolutionFromLocalID( iElem, convertEdgeIDs[iEdge] );

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[elemIDGlobal][edgeIDLocal], val ) );

//#ifdef _DEBUG_WRITE
//				std::cout << "planeID elemID faceID edgeIDLocal edgeID value : " << planeID << " " << elemIDGlobal << " " << faceIDLocal << " " << edgeIDLocal << " " << m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] << " " << val << std::endl;
//#endif

			}

		}

	}

	// Top planes on which dirichlet boundary condition with non-zero values ( source field ) is specified
	{
		const double sourceValueElectric = CommonParameters::sourceValueElectric;

		const int nElemOnPlane = m_MeshDataTetraElement.getNumElemOnBoundaryPlanes( MeshData::XYMinus );

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){

			const int elemIDGlobal = m_MeshDataTetraElement.getElemBoundaryPlanes( MeshData::XYMinus, iElem );
			const int faceIDLocal = m_MeshDataTetraElement.getFaceIDLocalFromElementBoundaryPlanes( MeshData::XYMinus, iElem );
		
			for( int iEdge = 0; iEdge < 3; ++iEdge ){
				const int edgeIDLocal = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				const int nodeID0 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 0 ); 
				const int nodeID1 = m_MeshDataTetraElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 1 ); 
			
				double val = ( iPol == CommonParameters::EX_POLARIZATION ) ? m_MeshDataTetraElement.caldDiffXOfTwoNodes( nodeID0, nodeID1 ) : m_MeshDataTetraElement.caldDiffYOfTwoNodes( nodeID0, nodeID1 );
				val /= m_MeshDataTetraElement.calcDistanceOfTwoNodes( nodeID0, nodeID1 );
				if( m_signInversion[elemIDGlobal][edgeIDLocal] ){
					val *= -1.0;
				}

				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[elemIDGlobal][edgeIDLocal], std::complex<double>(sourceValueElectric*val, 0.0) ) );
			}

		}
	}
	
//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for ( std::map<int, std::complex<double> >::iterator itr = m_globalID2NonZeroValues.begin(); itr != m_globalID2NonZeroValues.end(); ++itr ){
		std::cout << "key value " << itr->first << " " << itr->second << std::endl;
	}
#endif
//----- debug <<<<<

}

// Get shape functions of the 1st direction with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::getShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const double wLocal, const int num ) const{

	switch( num ){
		case 0:
			return 1.0 - vLocal - wLocal;
			break;
		case 1:
			return vLocal;
			break;
		case 2:
			return wLocal;
			break;
		case 3:
			return - vLocal;
			break;
		case 4:
			return wLocal;
			break;
		case 5:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncReferenceCoordU : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of the 2nd direction with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::getShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const double wLocal, const int num ) const{

	switch( num ){
		case 0:
			return uLocal;
			break;
		case 1:
			return 1.0 - uLocal - wLocal;
			break;
		case 2:
			return wLocal;
			break;
		case 3:
			return uLocal;
			break;
		case 4:
			return 0.0;
			break;
		case 5:
			return - wLocal;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncReferenceCoordV : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of the 3rd direction with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::getShapeFuncReferenceCoordW( const double uLocal, const double vLocal, const double wLocal, const int num ) const{

	switch( num ){
		case 0:
			return uLocal;
			break;
		case 1:
			return vLocal;
			break;
		case 2:
			return 1.0 - uLocal - vLocal;
			break;
		case 3:
			return 0.0;
			break;
		case 4:
			return - uLocal;
			break;
		case 5:
			return vLocal;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncReferenceCoordW : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get 2D shape functions of the 1st direction with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::get2DShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const int num ) const{

	switch( num ){
		case 0:
			return 1 - vLocal;
			break;
		case 1:
			return - vLocal;
			break;
		case 2:
			return - vLocal;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in get2DShapeFuncReferenceCoordU : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get 2D shape functions of the 2nd direction with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::get2DShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const int num ) const{

	switch( num ){
		case 0:
			return uLocal;
			break;
		case 1:
			return uLocal;
			break;
		case 2:
			return uLocal - 1;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in get2DShapeFuncReferenceCoordV : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get 2D shape functions rotated with respect to the reference element coordinate system
inline double Forward3DTetraElement0thOrder::get2DShapeFuncRotated() const{

	return 2.0;

}

// Get shape functions rotated with respect to the reference element coordinate system ( The 1st direction )
inline double Forward3DTetraElement0thOrder::getShapeFuncRotatedReferenceCoordU( const int num ) const{

	switch( num ){
		case 0:
			return 0.0;
			break;
		case 1:
			return 2.0;
			break;
		case 2:
			return - 2.0;
			break;
		case 3:
			return 0.0;
			break;
		case 4:
			return 0.0;
			break;
		case 5:
			return 2.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedReferenceCoordU : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions rotated with respect to the reference element coordinate system ( The 2nd direction )
inline double Forward3DTetraElement0thOrder::getShapeFuncRotatedReferenceCoordV( const int num ) const{

	switch( num ){
		case 0:
			return -2.0;
			break;
		case 1:
			return 0.0;
			break;
		case 2:
			return 2.0;
			break;
		case 3:
			return 0.0;
			break;
		case 4:
			return 2.0;
			break;
		case 5:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedReferenceCoordV : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions rotated with respect to the reference element coordinate system ( The 3rd direction )
inline double Forward3DTetraElement0thOrder::getShapeFuncRotatedReferenceCoordW( const int num ) const{

	switch( num ){
		case 0:
			return 2.0;
			break;
		case 1:
			return -2.0;
			break;
		case 2:
			return 0.0;
			break;
		case 3:
			return 2.0;
			break;
		case 4:
			return 0.0;
			break;
		case 5:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotatedReferenceCoordU : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Calculate jacobian matrix of the elements
void Forward3DTetraElement0thOrder::calcJacobianMatrix( const int elemID, Forward3D::Matrix3x3& JacobMat, double& determinant ) const{

	const int nodeID[4] = { m_MeshDataTetraElement.getNodesOfElements( elemID, 0 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 1 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 2 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 3 ) };

	const double x0 = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[0] );
	const double y0 = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[0] );
	const double z0 = m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[0] );

	JacobMat.mat11 = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[1] ) - x0;
	JacobMat.mat12 = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[1] ) - y0;
	JacobMat.mat13 = m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[1] ) - z0;

	JacobMat.mat21 = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[2] ) - x0;
	JacobMat.mat22 = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[2] ) - y0;
	JacobMat.mat23 = m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[2] ) - z0;

	JacobMat.mat31 = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[3] ) - x0;
	JacobMat.mat32 = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[3] ) - y0;
	JacobMat.mat33 = m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[3] ) - z0;

	determinant = JacobMat.mat11 * JacobMat.mat22 * JacobMat.mat33
				+ JacobMat.mat12 * JacobMat.mat23 * JacobMat.mat31
				+ JacobMat.mat13 * JacobMat.mat21 * JacobMat.mat32
				- JacobMat.mat13 * JacobMat.mat22 * JacobMat.mat31
				- JacobMat.mat12 * JacobMat.mat21 * JacobMat.mat33
				- JacobMat.mat11 * JacobMat.mat23 * JacobMat.mat32;

//----- debug >>>>>
//#ifdef _DEBUG_WRITE
//	std::cout << "jacob11 jacob12 jacob13 : " << JacobMat.mat11 << " " << JacobMat.mat12 << " " << JacobMat.mat13 << std::endl;
//	std::cout << "jacob21 jacob22 jacob23 : " << JacobMat.mat21 << " " << JacobMat.mat22 << " " << JacobMat.mat23 << std::endl;
//	std::cout << "jacob31 jacob32 jacob33 : " << JacobMat.mat31 << " " << JacobMat.mat32 << " " << JacobMat.mat33 << std::endl;
//	std::cout << "determinant : " << determinant << std::endl;
//#endif
//----- debug <<<<<

}

// Calculate 2D jacobian matrix of the faces
void Forward3DTetraElement0thOrder::calc2DJacobianMatrix( const int elemID, const int faceID, Forward3D::Matrix2x2& JacobMat, double& determinant ) const{

	const int nodeID[3] = { m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 0 ),
							m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 1 ),
							m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 2 ) };

	CommonParameters::locationXY coord[3];
	for( int iNode = 0 ; iNode < 3; ++iNode ){
		coord[iNode].X = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[iNode] );
		coord[iNode].Y = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[iNode] );
	}

	double dLengdX(0.0);
	double dLengdY(0.0);
	// Calculate the inclinations of the element's face
	calcInclinationsOfElementFace(elemID, faceID, dLengdX, dLengdY);

	JacobMat.mat11 = ( coord[1].X - coord[0].X ) * dLengdX;
	JacobMat.mat12 = ( coord[1].Y - coord[0].Y ) * dLengdY;
	JacobMat.mat21 = ( coord[2].X - coord[0].X ) * dLengdX;
	JacobMat.mat22 = ( coord[2].Y - coord[0].Y ) * dLengdY;

	determinant = JacobMat.mat11 * JacobMat.mat22 - JacobMat.mat12 * JacobMat.mat21;

}

// Calculate the inclinations of the element's face
void Forward3DTetraElement0thOrder::calcInclinationsOfElementFace( const int elemID, const int faceID, double& dLengdX, double& dLengdY ) const{

	const int nodeID[3] = { m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 0 ),
							m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 1 ),
							m_MeshDataTetraElement.getNodeIDGlobalFromElementAndFace( elemID, faceID, 2 ) };

	CommonParameters::locationXYZ coord[3];
	for( int iNode = 0 ; iNode < 3; ++iNode ){
		coord[iNode].X = m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[iNode] );
		coord[iNode].Y = m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[iNode] );
		coord[iNode].Z = m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[iNode] );
	}

	const double angle1 = atan2( coord[1].Y - coord[0].Y, coord[1].X - coord[0].X );
	const double angle2 = atan2( coord[2].Y - coord[0].Y, coord[2].X - coord[0].X );

	const double dzdl1 = ( coord[1].Z  - coord[0].Z ) / m_MeshDataTetraElement.calcHorizontalDistanceOfTwoNodes( nodeID[1], nodeID[0] );
	const double dzdl2 = ( coord[2].Z  - coord[0].Z ) / m_MeshDataTetraElement.calcHorizontalDistanceOfTwoNodes( nodeID[2], nodeID[0] );

	const double dzdx = ( dzdl1 * sin( angle2 ) - dzdl2 * sin( angle1 ) ) / sin( angle2 - angle1 );
	const double dzdy = ( dzdl1 * sin( angle2 - 0.5 * CommonParameters::PI ) - dzdl2 * sin( angle1 - 0.5 * CommonParameters::PI ) ) / sin( angle2 - angle1 );

	dLengdX = sqrt( 1.0 + dzdx * dzdx );
	dLengdY = sqrt( 1.0 + dzdy * dzdy );

}

// Calculate inverse of jacobian matrix  multiplied by determinant
void Forward3DTetraElement0thOrder::calcInverseOfJacobianMatrix( const Forward3D::Matrix3x3& jacobMat, Forward3D::Matrix3x3& invJacobMat ) const{

	invJacobMat.mat11 = jacobMat.mat22 * jacobMat.mat33 - jacobMat.mat23 * jacobMat.mat32; 
	invJacobMat.mat12 = jacobMat.mat13 * jacobMat.mat32 - jacobMat.mat12 * jacobMat.mat33; 
	invJacobMat.mat13 = jacobMat.mat12 * jacobMat.mat23 - jacobMat.mat13 * jacobMat.mat22; 

	invJacobMat.mat21 = jacobMat.mat23 * jacobMat.mat31 - jacobMat.mat21 * jacobMat.mat33; 
	invJacobMat.mat22 = jacobMat.mat11 * jacobMat.mat33 - jacobMat.mat13 * jacobMat.mat31; 
	invJacobMat.mat23 = jacobMat.mat13 * jacobMat.mat21 - jacobMat.mat11 * jacobMat.mat23; 

	invJacobMat.mat31 = jacobMat.mat21 * jacobMat.mat32 - jacobMat.mat22 * jacobMat.mat31; 
	invJacobMat.mat32 = jacobMat.mat12 * jacobMat.mat31 - jacobMat.mat11 * jacobMat.mat32; 
	invJacobMat.mat33 = jacobMat.mat11 * jacobMat.mat22 - jacobMat.mat12 * jacobMat.mat21; 
	
}

//// Calculate flag of sign inversion. This function used by calcVoltageDifference.
//bool Forward3DTetraElement0thOrder::calcReverseFlag( const int faceID, const int edgeIDLocal2D ) const{
//
//#ifdef _DEBUG_WRITE
//	std::cout << "faceID : " << faceID << std::endl;
//	std::cout << "edgeIDLocal2D : " << edgeIDLocal2D << std::endl;
//#endif
//
//	int ibuf1(0);
//	int ibuf2(0);
//	switch( edgeIDLocal2D ){
//		case 0:
//			ibuf1 = 0;
//			ibuf2 = 1;
//			break;
//		case 1:
//			ibuf1 = 1;
//			ibuf2 = 2;
//			break;
//		case 2:
//			ibuf1 = 2;
//			ibuf2 = 0;
//			break;
//		default:
//			OutputFiles::m_logFile << "Error : Wrong edge ID. edgeIDLocal2D =  " << edgeIDLocal2D << " !!" << std::endl;
//			exit(1);
//			break;
//	}
//
//	const int nodeIDLocal3DFromFace[2] = { m_MeshDataTetraElement.getNodeIDLocalFromFaceIDLocal( faceID, ibuf1 ),
//											m_MeshDataTetraElement.getNodeIDLocalFromFaceIDLocal( faceID, ibuf2 ) };
//
//	const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D );
//
//	const int nodeIDLocal3DFromEdge[2] = { m_MeshDataTetraElement.getNodeIDLocalFromEdgeIDLocal( edgeIDLocal3D, 0 ),
//											m_MeshDataTetraElement.getNodeIDLocalFromEdgeIDLocal( edgeIDLocal3D, 1 ) };
//	
//#ifdef _DEBUG_WRITE
//	std::cout << "edgeIDLocal3D : " << edgeIDLocal3D << std::endl;
//	std::cout << "nodeIDLocal2DFromEdge : " << nodeIDLocal3DFromEdge[0] << " " << nodeIDLocal3DFromEdge[1] << std::endl;
//	std::cout << "nodeIDLocal2DFromFace : " << nodeIDLocal3DFromFace[0] << " " << nodeIDLocal3DFromFace[1] << std::endl;
//#endif
//
//	if( nodeIDLocal3DFromFace[0] == nodeIDLocal3DFromEdge[0] && nodeIDLocal3DFromFace[1] == nodeIDLocal3DFromEdge[1] ){
//		return false;
//	}else if( nodeIDLocal3DFromFace[0] == nodeIDLocal3DFromEdge[1] && nodeIDLocal3DFromFace[1] == nodeIDLocal3DFromEdge[0] ){
//		return true;
//	}else{
//		OutputFiles::m_logFile << "Error : nodeIDLocal3DFromFace does not match nodeIDLocal3DFromEdge !! nodeIDLocal3DFromFace : " <<
//			nodeIDLocal3DFromFace[0] << " " << nodeIDLocal3DFromFace[1] << ", nodeIDLocal3DFromEdge : " << nodeIDLocal3DFromEdge[0] << " " << nodeIDLocal3DFromEdge[1] << std::endl;
//		exit(1);
//	}
//
//	return false;
//	
//}

// Calculate flag of sign inversion. This function used by calcVoltageDifference.
bool Forward3DTetraElement0thOrder::calcRatioAndReverseFlag( const int faceID, const int edgeIDLocal2D,
	const CommonParameters::AreaCoords& startPoint, const CommonParameters::AreaCoords& endPoint, double& ratio ) const{

#ifdef _DEBUG_WRITE
	std::cout << "faceID : " << faceID << std::endl;
	std::cout << "edgeIDLocal2D : " << edgeIDLocal2D << std::endl;
#endif

	int nodeIDLocal2D[2];
	switch( edgeIDLocal2D ){
		case 0:
			ratio = endPoint.coord1 - startPoint.coord1;
			if( ratio > 0){
				nodeIDLocal2D[0] = 0;
				nodeIDLocal2D[1] = 1;
			}else{
				nodeIDLocal2D[0] = 1;
				nodeIDLocal2D[1] = 0;
			}
			break;
		case 1:
			ratio = endPoint.coord2 - startPoint.coord2;
			if( ratio > 0 ){
				nodeIDLocal2D[0] = 1;
				nodeIDLocal2D[1] = 2;
			}else{
				nodeIDLocal2D[0] = 2;
				nodeIDLocal2D[1] = 1;
			}
			break;
		case 2:
			ratio = endPoint.coord0 - startPoint.coord0;
			if( ratio > 0 ){
				nodeIDLocal2D[0] = 2;
				nodeIDLocal2D[1] = 0;
			}else{
				nodeIDLocal2D[0] = 0;
				nodeIDLocal2D[1] = 2;
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong edge ID. edgeIDLocal2D =  " << edgeIDLocal2D << " !!" << std::endl;
			exit(1);
			break;
	}

	if( ratio < 0 ){
		ratio *= -1.0;
	}

	//std::cout << "ratio : " << ratio << std::endl;

	const int nodeIDLocal3D[2] = { m_MeshDataTetraElement.getNodeIDLocalFromFaceIDLocal( faceID, nodeIDLocal2D[0] ),
											m_MeshDataTetraElement.getNodeIDLocalFromFaceIDLocal( faceID, nodeIDLocal2D[1] ) };

	const int edgeIDLocal3D = m_MeshDataTetraElement.getEdgeIDLocalFromFaceIDLocal( faceID, edgeIDLocal2D );

	const int nodeIDLocal3DFromEdge[2] = { m_MeshDataTetraElement.getNodeIDLocalFromEdgeIDLocal( edgeIDLocal3D, 0 ),
											m_MeshDataTetraElement.getNodeIDLocalFromEdgeIDLocal( edgeIDLocal3D, 1 ) };
	
#ifdef _DEBUG_WRITE
	std::cout << "edgeIDLocal3D : " << edgeIDLocal3D << std::endl;
	std::cout << "nodeIDLocal3DFromEdge : " << nodeIDLocal3DFromEdge[0] << " " << nodeIDLocal3DFromEdge[1] << std::endl;
	std::cout << "nodeIDLocal3D : " << nodeIDLocal3D[0] << " " << nodeIDLocal3D[1] << std::endl;
#endif

	if( nodeIDLocal3D[0] == nodeIDLocal3DFromEdge[0] && nodeIDLocal3D[1] == nodeIDLocal3DFromEdge[1] ){
		return false;
	}else if( nodeIDLocal3D[0] == nodeIDLocal3DFromEdge[1] && nodeIDLocal3D[1] == nodeIDLocal3DFromEdge[0] ){
		return true;
	}else{
		OutputFiles::m_logFile << "Error : nodeIDLocal3D does not match nodeIDLocal3DFromEdge !! nodeIDLocal3D : " <<
			nodeIDLocal3D[0] << " " << nodeIDLocal3D[1] << ", nodeIDLocal3DFromEdge : " << nodeIDLocal3DFromEdge[0] << " " << nodeIDLocal3DFromEdge[1] << std::endl;
		exit(1);
	}

	return false;
}

//// Calculate area coordinates of shared point and rotation direction
//// Rotation direction is +Z : return true
//// Rotation direction is -Z : return false
//bool Forward3DTetraElement0thOrder::calcAreaCoordAndRotationDirection( const int edgeIDLocal2DStart, const int edgeIDLocal2DEnd, CommonParameters::AreaCoords& areaCoordOfSharePoint ) const{
//
//	if( edgeIDLocal2DStart == 0 && edgeIDLocal2DEnd == 1 ){
//		areaCoordOfSharePoint.coord1 = 1.0;
//		return false;
//	}else if( edgeIDLocal2DStart == 1 && edgeIDLocal2DEnd == 0 ){
//		areaCoordOfSharePoint.coord1 = 1.0;
//		return true;
//	}//-----------------------------------------------------------------
//	else if( edgeIDLocal2DStart == 1 && edgeIDLocal2DEnd == 2 ){
//		areaCoordOfSharePoint.coord2 = 1.0;
//		return false;
//	}else if( edgeIDLocal2DStart == 2 && edgeIDLocal2DEnd == 1 ){
//		areaCoordOfSharePoint.coord2 = 1.0;
//		return true;
//	}//-----------------------------------------------------------------
//	else if( edgeIDLocal2DStart == 2 && edgeIDLocal2DEnd == 0 ){
//		areaCoordOfSharePoint.coord0 = 1.0;
//		return false;
//	}else if( edgeIDLocal2DStart == 0 && edgeIDLocal2DEnd == 2 ){
//		areaCoordOfSharePoint.coord0 = 1.0;
//		return true;
//	}//-----------------------------------------------------------------
//
//	OutputFiles::m_logFile << "Error : Wrong edge IDs !! edgeIDLocal2DStart : " << edgeIDLocal2DStart << ", edgeIDLocal2DEnd : " << edgeIDLocal2DEnd << std::endl;
//	exit(1);
//
//	return false;
//
//}
//
//// Calculate area coordinates of shared point
//void Forward3DTetraElement0thOrder::calcAreaCoordOfSharedPoint( const int edgeIDLocal2DStart, const int edgeIDLocal2DEnd, CommonParameters::AreaCoords& areaCoordOfSharePoint ) const{
//
//	areaCoordOfSharePoint.coord0 = 0.0;
//	areaCoordOfSharePoint.coord1 = 0.0;
//	areaCoordOfSharePoint.coord2 = 0.0;
//	if( ( edgeIDLocal2DStart == 0 && edgeIDLocal2DEnd == 1 ) || ( edgeIDLocal2DStart == 1 && edgeIDLocal2DEnd == 0 ) ){
//		areaCoordOfSharePoint.coord1 = 1.0;
//	}else if( ( edgeIDLocal2DStart == 1 && edgeIDLocal2DEnd == 2 ) || ( edgeIDLocal2DStart == 2 && edgeIDLocal2DEnd == 1 ) ){
//		areaCoordOfSharePoint.coord2 = 1.0;
//	}else if( ( edgeIDLocal2DStart == 2 && edgeIDLocal2DEnd == 0 ) || ( edgeIDLocal2DStart == 0 && edgeIDLocal2DEnd == 2 ) ){
//		areaCoordOfSharePoint.coord0 = 1.0;
//	}
//
//	OutputFiles::m_logFile << "Error : Wrong edge IDs !! edgeIDLocal2DStart : " << edgeIDLocal2DStart << ", edgeIDLocal2DEnd : " << edgeIDLocal2DEnd << std::endl;
//	exit(1);
//
//}

// Output results of forward calculation to VTK file
void Forward3DTetraElement0thOrder::outputResultToVTK() const{

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
	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

	const double freq = getFrequencyCurrent();
	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_ELECTRIC_FIELD_VECTORS_TO_VTK ) ){// Output electric field vector
		OutputFiles::m_vtkFile << "VECTORS Re(E)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Ex = static_cast<float>( real( calcValueElectricFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Ey = static_cast<float>( real( calcValueElectricFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Ez = static_cast<float>( real( calcValueElectricFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(E)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Ex = static_cast<float>( imag( calcValueElectricFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Ey = static_cast<float>( imag( calcValueElectricFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Ez = static_cast<float>( imag( calcValueElectricFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << Ex << " " << Ey << " " << Ez << std::endl;
		}

	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_MAGNETIC_FIELD_VECTORS_TO_VTK ) ){// Output magnetic field vector
		OutputFiles::m_vtkFile << "VECTORS Re(H)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Hx = static_cast<float>( real( calcValueMagneticFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Hy = static_cast<float>( real( calcValueMagneticFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Hz = static_cast<float>( real( calcValueMagneticFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(H)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const float Hx = static_cast<float>( imag( calcValueMagneticFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Hy = static_cast<float>( imag( calcValueMagneticFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float Hz = static_cast<float>( imag( calcValueMagneticFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << Hx << " " << Hy << " " << Hz << std::endl;
		}
	}

	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_CURRENT_DENSITY ) ){// Output corrent density
		const ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();
		OutputFiles::m_vtkFile << "VECTORS Re(j)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			const float jx = static_cast<float>( sigma * real( calcValueElectricFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float jy = static_cast<float>( sigma * real( calcValueElectricFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float jz = static_cast<float>( sigma * real( calcValueElectricFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
		}

		OutputFiles::m_vtkFile << "VECTORS Im(j)_" << freq <<"(Hz)_" << stringPolarization << " float" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(iElem);
			const float jx = static_cast<float>( sigma * imag( calcValueElectricFieldXDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float jy = static_cast<float>( sigma * imag( calcValueElectricFieldYDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			const float jz = static_cast<float>( sigma * imag( calcValueElectricFieldZDirection( iElem, 0.25, 0.25, 0.25 ) ) );
			OutputFiles::m_vtkFile << jx << " " << jy << " " << jz << std::endl;
		}
	}
		
}

// Output results of forward calculation to binary file
void Forward3DTetraElement0thOrder::outputResultToBinary( const int iFreq, const int iPol ) const{
	
	const std::string stringPolarization = iPol == CommonParameters::EX_POLARIZATION ? "ExPol" : "EyPol";

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();

	//--- Total element number
	const int nElem = m_MeshDataTetraElement.getNumElemTotal();

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

		strcpy( line, "tetra4" );
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

		strcpy( line, "tetra4" );
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

		strcpy( line, "tetra4" );
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

		strcpy( line, "tetra4" );
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

		strcpy( line, "tetra4" );
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

		strcpy( line, "tetra4" );
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

// Calculate integrals of which element matrix consisits
void Forward3DTetraElement0thOrder::calcIntegrals( const int elemID, double* eMat, double* fMat ) const{

	const int nodeID[4] = { m_MeshDataTetraElement.getNodesOfElements( elemID, 0 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 1 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 2 ),
							m_MeshDataTetraElement.getNodesOfElements( elemID, 3 ) };

	double rowVecX[4] = { m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[0] ),
		                  m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[1] ),
						  m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[2] ),
						  m_MeshDataTetraElement.getXCoordinatesOfNodes( nodeID[3] ) };
	double rowVecY[4] = { m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[0] ),
		                  m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[1] ),
						  m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[2] ),
						  m_MeshDataTetraElement.getYCoordinatesOfNodes( nodeID[3] ) };
	double rowVecZ[4] = { m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[0] ),
		                  m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[1] ),
						  m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[2] ),
						  m_MeshDataTetraElement.getZCoordinatesOfNodes( nodeID[3] ) };
	const double rowOnes[4] = { 1.0, 1.0, 1.0, 1.0 };

	// irow = 0 -------------------------
	const double a[4] = {    calcDeterminant( rowVecX, rowVecY, rowVecZ, 0 ),
	                       - calcDeterminant( rowVecX, rowVecY, rowVecZ, 1 ),
	                         calcDeterminant( rowVecX, rowVecY, rowVecZ, 2 ),
  						   - calcDeterminant( rowVecX, rowVecY, rowVecZ, 3 ) };
	//-----------------------------------
	// irow = 1 -------------------------
	const double b[4] = {  - calcDeterminant( rowOnes, rowVecY, rowVecZ, 0 ),
	                         calcDeterminant( rowOnes, rowVecY, rowVecZ, 1 ),
	                       - calcDeterminant( rowOnes, rowVecY, rowVecZ, 2 ),
  						     calcDeterminant( rowOnes, rowVecY, rowVecZ, 3 ) };
	//-----------------------------------
	// irow = 2 -------------------------
	const double c[4] = {    calcDeterminant( rowVecX, rowOnes, rowVecZ, 0 ),
	                       - calcDeterminant( rowVecX, rowOnes, rowVecZ, 1 ),
	                         calcDeterminant( rowVecX, rowOnes, rowVecZ, 2 ),
  						   - calcDeterminant( rowVecX, rowOnes, rowVecZ, 3 ) };
	//-----------------------------------
	// irow = 3 -------------------------
	const double d[4] = {  - calcDeterminant( rowVecX, rowVecY, rowOnes, 0 ),
	                         calcDeterminant( rowVecX, rowVecY, rowOnes, 1 ),
	                       - calcDeterminant( rowVecX, rowVecY, rowOnes, 2 ),
  						     calcDeterminant( rowVecX, rowVecY, rowOnes, 3 ) };
	//-----------------------------------

	for( int i = 1; i < 4; ++i ){
		rowVecX[i] -= rowVecX[0];
		rowVecY[i] -= rowVecY[0];
		rowVecZ[i] -= rowVecZ[0];
	}
	const double detJacob = calcDeterminant( rowVecX, rowVecY, rowVecZ, 0 );
	const double vol = detJacob / 6.0;

	const int rowIndex[6] = { 0, 6, 11, 15, 18, 20 };

	//--------------------------
	//--- Calculate E Matrix ---
	//--------------------------
	const double factor = 2.0 / detJacob / detJacob; 
	const double length[6] = { m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 0 ),
	                           m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 1 ),
	                           m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 2 ),
	                           m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 3 ),
	                           m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 4 ),
							   m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge( elemID, 5 ) };

	const double coeff[6][3] = { { factor * length[0] * ( c[0]*d[1] - d[0]*c[1] ),
	                               factor * length[0] * ( d[0]*b[1] - b[0]*d[1] ),
				    			   factor * length[0] * ( b[0]*c[1] - c[0]*b[1] ) },
	                             { factor * length[1] * ( c[0]*d[2] - d[0]*c[2] ),
								   factor * length[1] * ( d[0]*b[2] - b[0]*d[2] ),
								   factor * length[1] * ( b[0]*c[2] - c[0]*b[2] ) },
	                             { factor * length[2] * ( c[0]*d[3] - d[0]*c[3] ),
	                               factor * length[2] * ( d[0]*b[3] - b[0]*d[3] ),
								   factor * length[2] * ( b[0]*c[3] - c[0]*b[3] ) },
 	                             { factor * length[3] * ( c[1]*d[2] - d[1]*c[2] ),
	                               factor * length[3] * ( d[1]*b[2] - b[1]*d[2] ),
								   factor * length[3] * ( b[1]*c[2] - c[1]*b[2] ) },
								 { factor * length[4] * ( c[3]*d[1] - d[3]*c[1] ),
	                               factor * length[4] * ( d[3]*b[1] - b[3]*d[1] ),
								   factor * length[4] * ( b[3]*c[1] - c[3]*b[1] ) },
                                 { factor * length[5] * ( c[2]*d[3] - d[2]*c[3] ),
	                               factor * length[5] * ( d[2]*b[3] - b[2]*d[3] ),
								   factor * length[5] * ( b[2]*c[3] - c[2]*b[3] ) } };

	for( int iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){
		for( int iEdge2 = iEdge1; iEdge2 < 6; ++iEdge2 ){
			eMat[ rowIndex[iEdge1] + iEdge2 - iEdge1 ] = vol * ( coeff[iEdge1][0]*coeff[iEdge2][0] + coeff[iEdge1][1]*coeff[iEdge2][1] + coeff[iEdge1][2]*coeff[iEdge2][2] );
		}
	}

	//--------------------------
	//--- Calculate F Matrix ---
	//--------------------------
	double f[16];
	for( int i = 0; i < 4; ++i ){
		for( int j = i; j < 4; ++j ){
			f[ rowIndex[i] + j - i ] = b[i]*b[j] + c[i]*c[j] + d[i]*d[j];
		}
	}

	fMat[0]  = length[0]*length[0]/(360.0*vol)*(     f[6] -     f[1]               +     f[0]  );// F00
	fMat[1]  = length[0]*length[1]/(720.0*vol)*( 2.0*f[7] -     f[1]   -     f[2]  +     f[0]  );// F01
	fMat[2]  = length[0]*length[2]/(720.0*vol)*( 2.0*f[8] -     f[1]   -     f[3]  +     f[0]  );// F02
	fMat[3]  = length[0]*length[3]/(720.0*vol)*(     f[7] -     f[6]   - 2.0*f[2]  +     f[1]  );// F03
	fMat[4]  = length[0]*length[4]/(720.0*vol)*(     f[6] -     f[8]   -     f[1]  + 2.0*f[3]  );// F04
	fMat[5]  = length[0]*length[5]/(720.0*vol)*(     f[8] -     f[7]   -     f[3]  +     f[2]  );// F05
 
	fMat[6]  = length[1]*length[1]/(360.0*vol)*(     f[11] -    f[2]               +     f[0]  );// F11
	fMat[7]  = length[1]*length[2]/(720.0*vol)*( 2.0*f[12] -    f[2]   -     f[3]  +     f[0]  );// F12
	fMat[8]  = length[1]*length[3]/(720.0*vol)*(     f[11] -    f[7]   -     f[2]  + 2.0*f[1]  );// F13
	fMat[9]  = length[1]*length[4]/(720.0*vol)*(     f[7]  -    f[12]  -     f[1]  +     f[3]  );// F14
	fMat[10] = length[1]*length[5]/(720.0*vol)*(     f[2]  -    f[11]  - 2.0*f[3]  +     f[12] );// F15

	fMat[11] = length[2]*length[2]/(360.0*vol)*(     f[15] -     f[3]              +     f[0]  );// F22
	fMat[12] = length[2]*length[3]/(720.0*vol)*(     f[12] -     f[8]  -     f[2]  +     f[1]  );// F23
	fMat[13] = length[2]*length[4]/(720.0*vol)*(     f[8]  -     f[15] - 2.0*f[1]  +     f[3]  );// F24
	fMat[14] = length[2]*length[5]/(720.0*vol)*(     f[15] -     f[12] -     f[3]  + 2.0*f[2]  );// F25

	fMat[15] = length[3]*length[3]/(360.0*vol)*(     f[11] -     f[7]              +     f[6]  );// F33
	fMat[16] = length[3]*length[4]/(720.0*vol)*(     f[7]  - 2.0*f[12] -     f[6]  +     f[8]  );// F34
	fMat[17] = length[3]*length[5]/(720.0*vol)*(     f[12] -     f[11] - 2.0*f[8]  +     f[7]  );// F35

	fMat[18] = length[4]*length[4]/(360.0*vol)*(     f[6]  -     f[8]              +     f[15] );// F44
	fMat[19] = length[4]*length[5]/(720.0*vol)*(     f[8]  - 2.0*f[7]  -     f[15] +     f[12] );// F45

	fMat[20] = length[5]*length[5]/(360.0*vol)*(     f[15] -     f[12]             +     f[11] );// F55

}

// Calculate determinant of 3x3 matrix
double Forward3DTetraElement0thOrder::calcDeterminant( const double* rowVec0, const double* rowVec1, const double* rowVec2, const int icol ) const{

	switch(icol){
		case 0:
			return   rowVec0[1] * rowVec1[2] * rowVec2[3] + rowVec0[2] * rowVec1[3] * rowVec2[1] + rowVec0[3] * rowVec1[1] * rowVec2[2]
				   - rowVec0[3] * rowVec1[2] * rowVec2[1] - rowVec0[2] * rowVec1[1] * rowVec2[3] - rowVec0[1] * rowVec1[3] * rowVec2[2];
			break;
		case 1:
			return   rowVec0[0] * rowVec1[2] * rowVec2[3] + rowVec0[2] * rowVec1[3] * rowVec2[0] + rowVec0[3] * rowVec1[0] * rowVec2[2]
				   - rowVec0[3] * rowVec1[2] * rowVec2[0] - rowVec0[2] * rowVec1[0] * rowVec2[3] - rowVec0[0] * rowVec1[3] * rowVec2[2];
			break;
		case 2:
			return   rowVec0[0] * rowVec1[1] * rowVec2[3] + rowVec0[1] * rowVec1[3] * rowVec2[0] + rowVec0[3] * rowVec1[0] * rowVec2[1]
				   - rowVec0[3] * rowVec1[1] * rowVec2[0] - rowVec0[1] * rowVec1[0] * rowVec2[3] - rowVec0[0] * rowVec1[3] * rowVec2[1];
			break;
		case 3:
			return   rowVec0[0] * rowVec1[1] * rowVec2[2] + rowVec0[1] * rowVec1[2] * rowVec2[0] + rowVec0[2] * rowVec1[0] * rowVec2[1]
				   - rowVec0[2] * rowVec1[1] * rowVec2[0] - rowVec0[1] * rowVec1[0] * rowVec2[2] - rowVec0[0] * rowVec1[2] * rowVec2[1];
			break;
		default:
			OutputFiles::m_logFile << "Error : icol must be 0, 1, 2, or 3 . icol = " << icol << std::endl;
			exit(1);
			break;
	}

	exit(1);
	return 0.0;

}
