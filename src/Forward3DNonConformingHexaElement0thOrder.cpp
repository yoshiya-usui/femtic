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
#include "Forward3DNonConformingHexaElement0thOrder.h"
#include "MeshDataNonConformingHexaElement.h"
#include "ResistivityBlock.h"
#include "CommonParameters.h"
#include "OutputFiles.h"
#include "ObservedData.h"
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <assert.h>

const double Forward3DNonConformingHexaElement0thOrder::m_eps = 1.0e-12;

Forward3DNonConformingHexaElement0thOrder::Forward3DNonConformingHexaElement0thOrder():
	Forward3D(),
	m_slaveDofToMasterDofAndFactors(NULL),
	m_hasMadeMapSlaveDofToMasterDofAndFactors(false),
	m_IDsAfterDegenerated2AfterConstrained(NULL),
	m_numOfEquationDegeneratedAndConstrained(0),
	m_solutionVectorDegeneratedAndConstrained(NULL),
	m_vectorMPCConstants(NULL)
{
	//---------------------------------------------------------------------------
	//--- Calculate integral points and weights of two point Gauss quadrature ---
	//---------------------------------------------------------------------------
	int ip(0);
	for( int i = 0; i < m_numGauss; ++i ){
		for( int j = 0; j < m_numGauss; ++j ){
			for( int k = 0; k < m_numGauss; ++k ){
				m_integralPointXi[ip]   = CommonParameters::abscissas2Point[i];
				m_integralPointEta[ip]  = CommonParameters::abscissas2Point[j];
				m_integralPointZeta[ip] = CommonParameters::abscissas2Point[k];
				m_weights[ip] = CommonParameters::weights2Point[i] * CommonParameters::weights2Point[j] * CommonParameters::weights2Point[k];
				++ip;
			}
		}
	}

	// Array of reference coord xi values for each node
	m_xiAtNode[0] = -1.0;
	m_xiAtNode[1] =  1.0;
	m_xiAtNode[2] =  1.0;
	m_xiAtNode[3] = -1.0;
	m_xiAtNode[4] = -1.0;
	m_xiAtNode[5] =  1.0;
	m_xiAtNode[6] =  1.0;
	m_xiAtNode[7] = -1.0;

	// Array of reference coord eta values for each node
	m_etaAtNode[0] = -1.0;
	m_etaAtNode[1] = -1.0;
	m_etaAtNode[2] =  1.0;
	m_etaAtNode[3] =  1.0;
	m_etaAtNode[4] = -1.0;
	m_etaAtNode[5] = -1.0;
	m_etaAtNode[6] =  1.0;
	m_etaAtNode[7] =  1.0;

	// Array of reference coord zeta values for each node
	m_zetaAtNode[0] = -1.0;
	m_zetaAtNode[1] = -1.0;
	m_zetaAtNode[2] = -1.0;
	m_zetaAtNode[3] = -1.0;
	m_zetaAtNode[4] =  1.0;
	m_zetaAtNode[5] =  1.0;
	m_zetaAtNode[6] =  1.0;
	m_zetaAtNode[7] =  1.0;

	// Array of reference coord xi values for each edge
	m_xiAtEdge[0]  = -9999.999;
	m_xiAtEdge[1]  = -9999.999;
	m_xiAtEdge[2]  = -9999.999;
	m_xiAtEdge[3]  = -9999.999;
	m_xiAtEdge[4]  = -1.0;
	m_xiAtEdge[5]  = -1.0;
	m_xiAtEdge[6]  =  1.0;
	m_xiAtEdge[7]  =  1.0;
	m_xiAtEdge[8]  = -1.0;
	m_xiAtEdge[9]  =  1.0;
	m_xiAtEdge[10] = -1.0;
	m_xiAtEdge[11] =  1.0;

	// Array of reference coord eta values for each edge
	m_etaAtEdge[0]  = -1.0;
	m_etaAtEdge[1]  =  1.0;
	m_etaAtEdge[2]  = -1.0;
	m_etaAtEdge[3]  =  1.0;
	m_etaAtEdge[4]  = -9999.999;
	m_etaAtEdge[5]  = -9999.999;
	m_etaAtEdge[6]  = -9999.999;
	m_etaAtEdge[7]  = -9999.999;
	m_etaAtEdge[8]  = -1.0;
	m_etaAtEdge[9]  = -1.0;
	m_etaAtEdge[10] =  1.0;
	m_etaAtEdge[11] =  1.0;

	// Array of reference coord zeta values for each edge
	m_zetaAtEdge[0]  = -1.0;
	m_zetaAtEdge[1]  = -1.0;
	m_zetaAtEdge[2]  =  1.0;
	m_zetaAtEdge[3]  =  1.0;
	m_zetaAtEdge[4]  = -1.0;
	m_zetaAtEdge[5]  =  1.0; 
	m_zetaAtEdge[6]  = -1.0;
	m_zetaAtEdge[7]  =  1.0;
	m_zetaAtEdge[8]  = -9999.999;
	m_zetaAtEdge[9]  = -9999.999;
	m_zetaAtEdge[10] = -9999.999;
	m_zetaAtEdge[11] = -9999.999;

	m_Fwd2DQuadrilateralElement[MeshData::ZXMinus] = new Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( MeshData::ZXMinus, CommonParameters::EX_POLARIZATION );
	m_Fwd2DQuadrilateralElement[MeshData::ZXPlus ] = new Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( MeshData::ZXPlus,  CommonParameters::EX_POLARIZATION );

	m_Fwd2DQuadrilateralElement[MeshData::YZMinus] = new Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( MeshData::YZMinus, CommonParameters::EY_POLARIZATION );
	m_Fwd2DQuadrilateralElement[MeshData::YZPlus ] = new Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( MeshData::YZPlus,  CommonParameters::EY_POLARIZATION );

}

//Destructer
Forward3DNonConformingHexaElement0thOrder::~Forward3DNonConformingHexaElement0thOrder(){

	if( m_slaveDofToMasterDofAndFactors != NULL ){
		delete [] m_slaveDofToMasterDofAndFactors;
		m_slaveDofToMasterDofAndFactors = NULL;
	}

	if( m_IDsAfterDegenerated2AfterConstrained != NULL ){
		delete [] m_IDsAfterDegenerated2AfterConstrained;
		m_IDsAfterDegenerated2AfterConstrained = NULL;
	}

	if( m_solutionVectorDegeneratedAndConstrained != NULL ){
		delete [] m_solutionVectorDegeneratedAndConstrained;
		m_solutionVectorDegeneratedAndConstrained = NULL;
	}

	if( m_vectorMPCConstants != NULL ){
		delete [] m_vectorMPCConstants;
		m_vectorMPCConstants = NULL;
	}
	
}

// Run 3D forward calculation
void Forward3DNonConformingHexaElement0thOrder::forwardCalculation( const double freq, const int iPol ){

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

	OutputFiles::m_logFile << "# Start 3D forward calculation with the 0th order 3D edge-based element." << pAnalysisControl->outputElapsedTime() << std::endl;

	//-----------------------------------------------------------------------------
	//--- Set Number of equations and array converting local IDs to global ones ---
	//-----------------------------------------------------------------------------
	if( !m_hasSetIDsLocal2Global ){
		calcArrayConvertLocalID2Global();
	}
	
	if( !m_hasIDsGlobal2AfterDegenerated[iPol] ){
		calcArrayConvertIDsGlobal2AfterDegenerated();
	}
	
	if( !m_hasMadeMapSlaveDofToMasterDofAndFactors ){
		makeMapSlaveDofToMasterDofAndFactors();
	}

	OutputFiles::m_logFile << "# Number of equation = " << m_numOfEquation
		<< ", Number of equation after degeneration = " << m_numOfEquationDegenerated
		<< ", Number of equation after constraint = " << m_numOfEquationDegeneratedAndConstrained << std::endl;

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

		m_matrix3DAnalysis.setDegreeOfEquation( m_numOfEquationDegeneratedAndConstrained );

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

	// Calculate MPC constants
	calcMPCConstants();

	// Set non-zero values of matrix and right-hande side vector for forward calculation
	setNonZeroValues( m_matrix3DAnalysis );

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
	//std::complex<double>* solutionConstrained = new std::complex<double>[m_numOfEquationDegeneratedAndConstrained];
	if( m_solutionVectorDegeneratedAndConstrained != NULL ){
		delete [] m_solutionVectorDegeneratedAndConstrained;
		m_solutionVectorDegeneratedAndConstrained = NULL;
	}
	m_solutionVectorDegeneratedAndConstrained = new std::complex<double>[m_numOfEquationDegeneratedAndConstrained];
	OutputFiles::m_logFile << "# Solve phase of matrix solver for 3D forward calculation."
			<<  " Polarization : " << iPol << "." 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix3DAnalysis.solvePhaseMatrixSolver( m_solutionVectorDegeneratedAndConstrained );//Solve phase of matrix solver

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquationDegeneratedAndConstrained; ++i ){
		std::cout << "i m_solutionVectorDegeneratedAndConstrained : " << i << " " << m_solutionVectorDegeneratedAndConstrained[i] << std::endl;
	}
#endif

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[ m_numOfEquation ];

	bool* alreadyFound = new bool[ m_numOfEquation ];
	for( int i = 0; i < m_numOfEquation; ++i ){
		alreadyFound[i] = false;
	}

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 12; ++iEdge ){
			const int iNum = m_IDsLocal2Global[iElem][iEdge];
			if( !alreadyFound[iNum] ){
				const int iNumDegenerated = m_IDsGlobal2AfterDegenerated[iPol][iNum];
				if( iNumDegenerated == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					m_solution[iNum] = m_globalID2NonZeroValues[iNum];
				}else{
					m_solution[iNum] = m_vectorMPCConstants[iNumDegenerated];
					const std::vector< std::pair<int,double> >& masterDofAndFactor = m_slaveDofToMasterDofAndFactors[iNumDegenerated];
					for( std::vector< std::pair<int,double> >::const_iterator itr = masterDofAndFactor.begin(); itr != masterDofAndFactor.end(); ++itr ){
						const int dof = m_IDsAfterDegenerated2AfterConstrained[itr->first];
						const std::complex<double> factor = std::complex<double>(itr->second, 0.0);
						m_solution[iNum] += m_solutionVectorDegeneratedAndConstrained[dof] * factor;
					}
				}
				alreadyFound[iNum] = true;
			}			
		}// iEdge
	}// iElem

	//delete[] solutionConstrained;
	delete[] alreadyFound;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquation; ++i ){
		std::cout << "i m_solution : " << i << " " << m_solution[i] << std::endl;
	}
#endif

	// Output EM field vectors
	if( pAnalysisControl->writeBinaryFormat() ){// BINARY
		outputResultToBinary( (ObservedData::getInstance())->getFreqIDs( freq ), iPol );
	}
	else{// ASCII
		outputResultToVTK();
	}

}

// Calculate X component of electric field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncX( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}
	return val;

}

// Calculate Y component of electric field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncY( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}
	return val;

}

// Calculate Z component of electric field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncZ( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}
	return val;

}

// Calculate Z component of rotated electric field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncRotatedZ( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}
	return val;

}

// Calculate X component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate Y component of electric field only from the edges on the Earth's surface
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate tangential electric field directed to X from all edges of owner element
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncX( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}

	const double dzdx = JacobMat.mat13 / JacobMat.mat11;
	const double factor = sqrt( 1.0 + dzdx * dzdx );
	return val * std::complex<double>(factor, 0.0);

}

// Calculate tangential electric field directed to Y from all edges of owner element
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncY( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}

	const double dzdy = JacobMat.mat23 / JacobMat.mat22;
	const double factor = sqrt( 1.0 + dzdy * dzdy );
	return val * std::complex<double>(factor, 0.0);

}

// Calculate tangential electric field directed to X
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate tangential electric field directed to Y
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate X component of magnetic field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncRotatedX( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Y component of magnetic field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncRotatedY( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 );
	}

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	val /= std::complex<double>(0.0, omega * CommonParameters::mu);

	return val;

}

// Calculate Z component of magnetic field
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcValueMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal ) const{

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	return calcValueRotatedElectricFieldZDirection( iElem, xLocal, yLocal, zLocal ) / std::complex<double>(0.0, omega * CommonParameters::mu);

}

// Calculate interpolator vector of X component of electric field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncX( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of Y component of electric field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncY( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of Z component of electric field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncZ( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of Z component of rotated electric field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfRotatedElectricFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncRotatedZ( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of X component of electric field only from the edges on the Earth's surface
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldXDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate interpolator vector of Y component of electric field only from the edges on the Earth's surface
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldYDirectionFromEdgesOnEarthSurface( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate interpolator vector of tangential electric field directed to X from all edges
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const double dzdx = JacobMat.mat13 / JacobMat.mat11;
	const double dLengdx = sqrt( 1.0 + dzdx * dzdx );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncX( xLocal, yLocal, zLocal, i, invJacobMat ) * length * dLengdx, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of tangential electric field directed to Y from all edges
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( const int iElem, const int iFace, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const double dzdy = JacobMat.mat23 / JacobMat.mat22;
	const double dLengdy = sqrt( 1.0 + dzdy * dzdy );

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncY( xLocal, yLocal, zLocal, i, invJacobMat ) * length * dLengdy, 0.0 ) * factor;
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of tangential electric field directed to X
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialX( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate interpolator vector of tangential electric field directed to Y
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfElectricFieldTangentialY( const int iElem, const int iFace, const double uCoord, const double vCoord, const int irhs, const std::complex<double>& factor ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate interpolator vector of X component of magnetic field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfMagneticFieldXDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncRotatedX( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor
									/ std::complex<double>(0.0, omega * CommonParameters::mu);
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of Y component of magnetic field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfMagneticFieldYDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global != NULL );

	Forward3D::Matrix3x3 JacobMat;
	const double detJacob = calcJacobianMatrix( iElem, xLocal, yLocal, zLocal, JacobMat );
	Forward3D::Matrix3x3 invJacobMat;
	calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency

	const int iPol = getPolarizationCurrent();
	for( int i = 0; i < 12; ++i ){
		const int irow = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][i] ];
		if( irow < 0 ){
			continue;
		}
		const double length = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( iElem, i );
		const std::complex<double> val = std::complex<double>( getShapeFuncRotatedY( xLocal, yLocal, zLocal, i, invJacobMat ) * length, 0.0 ) * factor
									/ std::complex<double>(0.0, omega * CommonParameters::mu);
		addValuesToRhsVectorsByConsideringMPC( irow, irhs, val );
	}

}

// Calculate interpolator vector of Z component of magnetic field
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfMagneticFieldZDirection( const int iElem, const double xLocal, const double yLocal, const double zLocal, const int irhs, const std::complex<double>& factor ){

	const double omega = 2.0 * CommonParameters::PI * getFrequencyCurrent();//Angular frequency
	const std::complex<double> factorMod = factor / std::complex<double>(0.0, omega * CommonParameters::mu);

	calcInterpolatorVectorOfRotatedElectricFieldZDirection( iElem, xLocal, yLocal, zLocal, irhs, factorMod );

}

// Calculate interpolator vector of difference of voltage
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint, const int irhs ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate interpolator vector of difference of voltage
void Forward3DNonConformingHexaElement0thOrder::calcInterpolatorVectorOfVoltageDifference( const int nElem, const int* elememtsIncludingDipole, const int* const facesIncludingDipole, 
	const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint, const int irhs ){
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Set non-zero strucuture of matrix for forward calculation
void Forward3DNonConformingHexaElement0thOrder::setNonZeroStrucuture( ComplexSparseSquareSymmetricMatrix& matrix ){

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = iElem;
		for( int iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){
			const int row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){
				const int col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col < 0 ){
					continue;
				}
				std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[row];
				std::vector< std::pair<int,double> >& colMasters= m_slaveDofToMasterDofAndFactors[col];
				for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					for( std::vector< std::pair<int,double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol ){
						const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
						if( colMod >= rowMod ){// Store only upper triangle part
							matrix.setStructureByTripletFormat( rowMod, colMod );
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}

}

// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward3DNonConformingHexaElement0thOrder::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ){

#ifdef _DEBUG_WRITE
	ComplexSparseSquareSymmetricMatrix matrixTemp(m_numOfEquationDegenerated);
#endif

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = iElem;
		//--- Calculate omega * mu * sigma
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		const double factor1 = 1.0;
		const std::complex<double> factor2 = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form
		double length[12];
		for( int i = 0; i < 12; ++i ){
			length[i] = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( elemID, i );
		}
		for( int iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){
			const int row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){
				const int col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for( int ip = 0; ip < m_numIntegralPoints; ++ip ){
					const double xi = m_integralPointXi[ip];
					const double eta = m_integralPointEta[ip];
					const double zeta = m_integralPointZeta[ip];
					Forward3D::Matrix3x3 JacobMat;
					const double detJacob = calcJacobianMatrix( elemID, xi, eta, zeta, JacobMat );
					Forward3D::Matrix3x3 invJacobMat;
					calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );
					integral1 += ( getShapeFuncRotatedX( xi, eta, zeta, iEdge1, invJacobMat )
						         * getShapeFuncRotatedX( xi, eta, zeta, iEdge2, invJacobMat )
								 + getShapeFuncRotatedY( xi, eta, zeta, iEdge1, invJacobMat ) 
						         * getShapeFuncRotatedY( xi, eta, zeta, iEdge2, invJacobMat )
								 + getShapeFuncRotatedZ( xi, eta, zeta, iEdge1, invJacobMat ) 
						         * getShapeFuncRotatedZ( xi, eta, zeta, iEdge2, invJacobMat ) )
							     * detJacob * m_weights[ip];
					integral2 += ( getShapeFuncX( xi, eta, zeta, iEdge1, invJacobMat )
						         * getShapeFuncX( xi, eta, zeta, iEdge2, invJacobMat )
						         + getShapeFuncY( xi, eta, zeta, iEdge1, invJacobMat )
								 * getShapeFuncY( xi, eta, zeta, iEdge2, invJacobMat )
						         + getShapeFuncZ( xi, eta, zeta, iEdge1, invJacobMat )
								 * getShapeFuncZ( xi, eta, zeta, iEdge2, invJacobMat ) )
								 * detJacob * m_weights[ip];
				}
				integral1 *= length[iEdge1] * length[iEdge2];
				integral2 *= length[iEdge1] * length[iEdge2];
				const std::complex<double> val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form
//#ifdef _DEBUG_WRITE
//				if( col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
//					matrixTemp.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
//				}else if( col >= row ){// Store only upper triangle part
//					matrixTemp.setStructureAndAddValueByTripletFormat( row, col, val );// Add to matrix
//				}
//#endif
				const std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[row];
				for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
					if( col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						// Add to right hand side vector
						matrix.addRightHandSideVector( rowMod, -valMod * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );
					}else{
						// Add to right hand side vector corresponding to MPC constants
						matrix.addRightHandSideVector( rowMod, valMod * m_vectorMPCConstants[col] );
						const std::vector< std::pair<int,double> >& colMasters= m_slaveDofToMasterDofAndFactors[col];
						for( std::vector< std::pair<int,double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol ){
							const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
							const std::complex<double> valModMod = valMod * std::complex<double>(itrCol->second, 0.0);
							if( colMod >= rowMod ){// Store only upper triangle part
								const int loc = matrix.checkAndGetLocationNonZeroValue( rowMod, colMod );
								matrix.addNonZeroValuesWithoutSearchingLocation( loc, valModMod );// Add to matrix
							}
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}

#ifdef _DEBUG_WRITE
	std::cout << "matrix" << std::endl;
	matrix.debugWriteMatrix();
	matrix.debugWriteRightHandSide();
#endif

//#ifdef _DEBUG_WRITE
//	matrixTemp.transformationByConstraintMatrix( m_numOfEquationDegeneratedAndConstrained, m_transposedConstraintMatrix );
//	matrixTemp.convertToCRSFormat();
//	std::cout << "matrixTemp" << std::endl;
//	matrixTemp.debugWriteMatrix();
//	matrixTemp.debugWriteRightHandSide();
//#endif

}

// Calculate vector x of the reciprocity algorithm of Rodi (1976)
void Forward3DNonConformingHexaElement0thOrder::calVectorXOfReciprocityAlgorithm( const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows ){

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

	for( int iElem = 0; iElem < nElem; ++iElem ){
		// [Attention] : You must use elemID instead of iElem from this line
		const int elemID = mdl2Elem[iElem].first;

//----- debug >>>>>
#ifdef _DEBUG_WRITE
		std::cout << "blkID iElem elemID = " << blkID << " " << iElem << " " << elemID << std::endl;
#endif
//----- debug <<<<<

		//--- Calculate omega * mu * sigma
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		const double factor1 = 0.0;
		const std::complex<double> factor2 = std::complex<double>(0.0, - omega * CommonParameters::mu * sigma * ln10 * mdl2Elem[iElem].second);// exp(-i*omega*t) form
		double length[12];
		for( int i = 0; i < 12; ++i ){
			length[i] = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge( elemID, i );
		}
		for( int iEdge1 = 0; iEdge1 < 12; ++iEdge1 ){
			const int row = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge1] ];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 12; ++iEdge2 ){
				const int col = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemID][iEdge2] ];
				if( col <= Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for( int ip = 0; ip < m_numIntegralPoints; ++ip ){
					const double xi = m_integralPointXi[ip];
					const double eta = m_integralPointEta[ip];
					const double zeta = m_integralPointZeta[ip];
					Forward3D::Matrix3x3 JacobMat;
					const double detJacob = calcJacobianMatrix( elemID, xi, eta, zeta, JacobMat );
					Forward3D::Matrix3x3 invJacobMat;
					calcInverseOfJacobianMatrix( JacobMat, detJacob, invJacobMat );
					integral1 += ( getShapeFuncRotatedX( xi, eta, zeta, iEdge1, invJacobMat )
						         * getShapeFuncRotatedX( xi, eta, zeta, iEdge2, invJacobMat )
								 + getShapeFuncRotatedY( xi, eta, zeta, iEdge1, invJacobMat ) 
						         * getShapeFuncRotatedY( xi, eta, zeta, iEdge2, invJacobMat )
								 + getShapeFuncRotatedZ( xi, eta, zeta, iEdge1, invJacobMat ) 
						         * getShapeFuncRotatedZ( xi, eta, zeta, iEdge2, invJacobMat ) )
							     * detJacob * m_weights[ip];
					integral2 += ( getShapeFuncX( xi, eta, zeta, iEdge1, invJacobMat )
						         * getShapeFuncX( xi, eta, zeta, iEdge2, invJacobMat )
						         + getShapeFuncY( xi, eta, zeta, iEdge1, invJacobMat )
								 * getShapeFuncY( xi, eta, zeta, iEdge2, invJacobMat )
						         + getShapeFuncZ( xi, eta, zeta, iEdge1, invJacobMat )
								 * getShapeFuncZ( xi, eta, zeta, iEdge2, invJacobMat ) )
								 * detJacob * m_weights[ip];
				}
				integral1 *= length[iEdge1] * length[iEdge2];
				integral2 *= length[iEdge1] * length[iEdge2];
				const std::complex<double> val = std::complex<double>( integral1 * factor1 , 0.0 ) - std::complex<double>( integral2, 0.0 ) * factor2;// exp(-i*omega*t) form
//#ifdef _DEBUG_WRITE
//				if( col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
//					matrixTemp.addRightHandSideVector( row, -val * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ] );// Add to right hand side vector
//				}else if( col >= row ){// Store only upper triangle part
//					matrixTemp.setStructureAndAddValueByTripletFormat( row, col, val );// Add to matrix
//				}
//#endif
				const std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[row];
				for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
					if( col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						nonZeroRows.push_back(rowMod);
						// Add to right hand side vector
						vecOut[rowMod] -= valMod * m_globalID2NonZeroValues[ m_IDsLocal2Global[elemID][iEdge2] ];
					}else{
						nonZeroRows.push_back(rowMod);
						// Add to right hand side vector corresponding to MPC constants
						vecOut[rowMod] += valMod * m_vectorMPCConstants[col];
						const std::vector< std::pair<int,double> >& colMasters= m_slaveDofToMasterDofAndFactors[col];
						for( std::vector< std::pair<int,double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol ){
							const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
							const std::complex<double> valModMod = valMod * std::complex<double>(itrCol->second, 0.0);
							vecOut[rowMod] -= valModMod * vecIn[colMod];
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}

	std::sort(nonZeroRows.begin(), nonZeroRows.end());
	nonZeroRows.erase( std::unique( nonZeroRows.begin(), nonZeroRows.end() ), nonZeroRows.end() );

}

// Copy solution vector degenerated
void Forward3DNonConformingHexaElement0thOrder::copySolutionVectorDegenerated( const int iPol, std::complex<double>* solutionVector ) const{

	for( int i = 0; i < m_numOfEquationDegeneratedAndConstrained; ++i ){
		solutionVector[i] = m_solutionVectorDegeneratedAndConstrained[i];
	}

}

// Call function inputMeshData of the class MeshData
void Forward3DNonConformingHexaElement0thOrder::callInputMeshData(){

	m_MeshDataNonConformingHexaElement.inputMeshData();

}

// Get pointer to the class MeshData
const MeshData* Forward3DNonConformingHexaElement0thOrder::getPointerToMeshData() const{

	return static_cast<const MeshData*>( &m_MeshDataNonConformingHexaElement ) ;

}

// Get pointer to the class MeshDataNonConformingHexaElement
const MeshDataNonConformingHexaElement* Forward3DNonConformingHexaElement0thOrder::getPointerToMeshDataNonConformingHexaElement() const{

	return &m_MeshDataNonConformingHexaElement;

}

// Calculate difference of voltage for brick element
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcVoltageDifference( const int nElem, const int* elememtsIncludingDipole,
	const CommonParameters::locationXY* localCoordinateValuesStartPoint, const CommonParameters::locationXY* localCoordinateValuesEndPoint ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Calculate difference of voltage for tetra element
std::complex<double> Forward3DNonConformingHexaElement0thOrder::calcVoltageDifference( const int nElem, const int* const elememtsIncludingDipole, const int* const facesIncludingDipole, 
	const CommonParameters::AreaCoords* const areaCoordValStartPoint, const CommonParameters::AreaCoords* const areaCoordValEndPoint ) const{
	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);
}

// Get total number of equations finally solved
int Forward3DNonConformingHexaElement0thOrder::getNumOfEquationFinallySolved() const{
	return m_numOfEquationDegeneratedAndConstrained;
}

// Calculate array converting local IDs to global ones
void Forward3DNonConformingHexaElement0thOrder::calcArrayConvertLocalID2Global(){

	typedef std::pair< int, int > NodeIDPair;
	typedef std::pair< int, int > ElemAndEdgeID;
	typedef std::pair< NodeIDPair, ElemAndEdgeID > InputedDataType;

	std::vector< InputedDataType > workVector;
	
	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 12; ++iEdge ){
			const int nodeID0 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 0 );
			const int nodeID1 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 1 ); 
			const NodeIDPair nodePairCur = nodeID1 > nodeID0 ? std::make_pair( nodeID0, nodeID1 ) : std::make_pair( nodeID1, nodeID0 );
			const ElemAndEdgeID elemAndEdgeIDCur = std::make_pair( iElem, iEdge );
			workVector.push_back( InputedDataType( std::pair< NodeIDPair, ElemAndEdgeID >(nodePairCur, elemAndEdgeIDCur) ) );
		}
	}
	std::sort( workVector.begin(), workVector.end() );

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
		m_IDsLocal2Global[iElem] = new int[12];
	}
	//-----------------------------------------

	std::pair<int, int> nodePairPre;
	int edgeIDGlobal(-1);
	std::vector<InputedDataType>::const_iterator itrEnd = workVector.end();
	for( std::vector<InputedDataType>::const_iterator itr = workVector.begin(); itr != itrEnd; ++itr ){
		const int elemIDGlobal = (*itr).second.first;
		const int edgeIDLocal = (*itr).second.second;
		const int nodePairCur0 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 0 );
		const int nodePairCur1 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 1 );

		if( ( nodePairCur0 == nodePairPre.first && nodePairCur1 == nodePairPre.second ) ||
			( nodePairCur0 == nodePairPre.second && nodePairCur1 == nodePairPre.first ) ){// Same edge with privious one
		}else{// New edge
			++edgeIDGlobal;
			nodePairPre.first  = nodePairCur0;
			nodePairPre.second = nodePairCur1;
		}
		m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] = edgeIDGlobal;
#ifdef _DEBUG_WRITE
		std::cout << "nodePairCur0 nodePairCur1 m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] : " << nodePairCur0 << " " << nodePairCur1 << " -> " << m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] << std::endl;
#endif
	}

	m_numOfEquation = edgeIDGlobal + 1;
	
	m_hasSetIDsLocal2Global = true;

}

// Calculate array converting global IDs to the ones after degeneration
void Forward3DNonConformingHexaElement0thOrder::calcArrayConvertIDsGlobal2AfterDegenerated(){

	if( !m_hasSetIDsLocal2Global ){
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
		const int nElemOnPlane = m_MeshDataNonConformingHexaElement.getNumElemOnBoundaryPlanes( planeID );
//#ifdef _DEBUG_WRITE
//		std::cout << "planeID " << planeID << std::endl;
//#endif
		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){
			const int elemIDGlobal = m_MeshDataNonConformingHexaElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataNonConformingHexaElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
//#ifdef _DEBUG_WRITE
//				std::cout << "elemIDGlobal edgeIDLocal m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] : " << elemIDGlobal << " " << edgeIDLocal << " " << m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] << std::endl;
//#endif
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] ] = Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE;
			}
		}
	}

	// Planes on which dirichlet boundary condition with zero values is specified
	for( int iPlane = 0; iPlane < 3; ++iPlane ){
		const int planeID = planesWithZeroValues[iPlane];
		const int nElemOnPlane = m_MeshDataNonConformingHexaElement.getNumElemOnBoundaryPlanes( planeID );
		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){
			const int elemIDGlobal = m_MeshDataNonConformingHexaElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataNonConformingHexaElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDGlobal][edgeIDLocal] ] = Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE;
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

	m_hasIDsGlobal2AfterDegenerated[iPol] = true;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquation; ++i ){
		std::cout << "i m_IDsGlobal2AfterDegenerated[iPol][i] : " << i << " " << m_IDsGlobal2AfterDegenerated[iPol][i] << std::endl;
	}
	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 12; ++iEdge ){
			std::cout << "iElem iEdge m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ] : "  << iElem << " "  << iEdge << " " << m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[iElem][iEdge] ] << std::endl;
		}
	}
#endif

}

// Calculate array converting global edge IDs non-zero electric field values specified to the edges
void Forward3DNonConformingHexaElement0thOrder::calcArrayConvertIDGlobal2NonZeroValues(){

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
	OutputFiles::m_logFile << "# Calculating EM field of the boundary planes with 2D forward analysis." << std::endl;
	if( iPol == CommonParameters::EX_POLARIZATION ){// Ex polarization
		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DQuadrilateralElement[MeshData::ZXMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataNonConformingHexaElement );
		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::ZXPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DQuadrilateralElement[MeshData::ZXPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataNonConformingHexaElement );
	}else{// Ey polarization
		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZMinus << " ]-----------------------------" << std::endl;
		m_Fwd2DQuadrilateralElement[MeshData::YZMinus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataNonConformingHexaElement );
		OutputFiles::m_logFile << "#-----------------------------[ Boundary Plane " << MeshData::YZPlus << " ]-----------------------------" << std::endl;
		m_Fwd2DQuadrilateralElement[MeshData::YZPlus]->calcEMFieldsOfBoundaryPlanes( freq, &m_MeshDataNonConformingHexaElement );
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
		const int nElemOnPlane = m_MeshDataNonConformingHexaElement.getNumElemOnBoundaryPlanes( planeID );
		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){
			const int elemIDGlobal = m_MeshDataNonConformingHexaElement.getElemBoundaryPlanes( planeID, iElem );
			const int faceIDLocal = m_MeshDataNonConformingHexaElement.getFaceIDLocalFromElementBoundaryPlanes( planeID, iElem );
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const std::complex<double> val = m_Fwd2DQuadrilateralElement[planeID]->getSolutionFromLocalID( iElem, iEdge );
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[elemIDGlobal][edgeIDLocal], val ) );
			}
		}
	}

	// Top planes on which dirichlet boundary condition with non-zero values ( source field ) is specified
	{
		const double sourceValueElectric = CommonParameters::sourceValueElectric;
		const int nElemOnPlane = m_MeshDataNonConformingHexaElement.getNumElemOnBoundaryPlanes( MeshData::XYMinus );
		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){
			const int elemIDGlobal = m_MeshDataNonConformingHexaElement.getElemBoundaryPlanes( MeshData::XYMinus, iElem );
			const int faceIDLocal = m_MeshDataNonConformingHexaElement.getFaceIDLocalFromElementBoundaryPlanes( MeshData::XYMinus, iElem );
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( faceIDLocal, iEdge );
				const int nodeID0 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 0 );
				const int nodeID1 = m_MeshDataNonConformingHexaElement.getNodeIDGlobalFromElementAndEdge( elemIDGlobal, edgeIDLocal, 1 );
				double val = ( iPol == CommonParameters::EX_POLARIZATION ) ? m_MeshDataNonConformingHexaElement.caldDiffXOfTwoNodes( nodeID0, nodeID1 ) : m_MeshDataNonConformingHexaElement.caldDiffYOfTwoNodes( nodeID0, nodeID1 );
				val /= m_MeshDataNonConformingHexaElement.calcDistanceOfTwoNodes( nodeID0, nodeID1 );
				m_globalID2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[elemIDGlobal][edgeIDLocal], std::complex<double>(sourceValueElectric*val, 0.0) ) );
			}
		}
	}

#ifdef _DEBUG_WRITE
	for ( std::map<int, std::complex<double> >::iterator itr = m_globalID2NonZeroValues.begin(); itr != m_globalID2NonZeroValues.end(); ++itr ){
		std::cout << "key value " << itr->first << " " << itr->second << std::endl;
	}
#endif

}


// Make map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
void Forward3DNonConformingHexaElement0thOrder::makeMapSlaveDofToMasterDofAndFactors(){

	const int iPol = getPolarizationCurrent();// Constraint matrix does not depend on the type of polarization

	if( m_IDsAfterDegenerated2AfterConstrained != NULL ){
		delete [] m_IDsAfterDegenerated2AfterConstrained;
		m_IDsAfterDegenerated2AfterConstrained = NULL;
	}
	m_IDsAfterDegenerated2AfterConstrained = new int[m_numOfEquationDegenerated];
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		m_IDsAfterDegenerated2AfterConstrained[i] = 0;
	}

	if( m_slaveDofToMasterDofAndFactors != NULL ){
		delete [] m_slaveDofToMasterDofAndFactors;
		m_slaveDofToMasterDofAndFactors = NULL;
	}
	m_slaveDofToMasterDofAndFactors = new std::vector< std::pair<int,double> >[m_numOfEquationDegenerated];

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = iElem;
		for( int iFace = 0; iFace < 6; ++iFace ){
			if( !m_MeshDataNonConformingHexaElement.faceSlaveElements(elemID, iFace) ){
				// This face does not have slave edges
				continue;
			}
			// Dofs of master edges
			int dofMasterBeforeDegenerated[4] = { -1, -1, -1, -1 };
			int dofMaster[4] = { -1, -1, -1, -1 };
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFace, iEdge );
				const int dofBeforeDegenerated = m_IDsLocal2Global[elemID][edgeIDLocal];
				dofMasterBeforeDegenerated[iEdge] = dofBeforeDegenerated;
				dofMaster[iEdge] = m_IDsGlobal2AfterDegenerated[iPol][ dofBeforeDegenerated ];
#ifdef _DEBUG_WRITE
				if( dofMaster[iEdge] < 0 ){
					std::cout << "iElem, iEdge, dofMaster :" << iElem << " " << iEdge <<  " " << dofMaster[iEdge] << std::endl;
				}
#endif
			}
			// Face index of neighbor element
			const int iFaceNeib = m_MeshDataNonConformingHexaElement.getFaceIndexOfNeighborElement(iFace);
			const int numNeibElements = m_MeshDataNonConformingHexaElement.getNumNeighborElement(elemID, iFace);
			if( numNeibElements == 2 ){
				// Dofs of slave edges
				int dofSlaves[2][4] = { { -1, -1, -1, -1 }, { -1, -1, -1, -1 } };
				for( int iNeib = 0; iNeib < 2; ++iNeib ){
					const int elemIDNeib = m_MeshDataNonConformingHexaElement.getIDOfNeighborElement(elemID, iFace, iNeib );
					for( int iEdgeNeib = 0; iEdgeNeib < 4; ++iEdgeNeib ){
						const int edgeIDLocalNeib = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFaceNeib, iEdgeNeib );
						dofSlaves[iNeib][iEdgeNeib] = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDNeib][edgeIDLocalNeib] ];
					}
				}
				assert( dofSlaves[0][3] == dofSlaves[1][2] );
				if( dofMaster[0] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][0], dofMaster[0], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[1][0], dofMaster[0], 1.0 );
				}
				if( dofMaster[1] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][1], dofMaster[1], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[1][1], dofMaster[1], 1.0 );
				}
				if( dofMaster[2] >= 0 ){
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][3], dofMaster[2], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				if( dofMaster[3] >= 0 ){
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][3], dofMaster[3], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				// Slave dofs on outer edges of the face
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][0] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[1][0] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[1][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				// Slave dofs on interior edges of the face
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][3] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_INTERIOR_EDGES;
			}else if( numNeibElements == 4 ){
				// Dofs of slave edges
				int dofSlaves[4][4] = { { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 } };
				for( int iNeib = 0; iNeib < 4; ++iNeib ){
					const int elemIDNeib = m_MeshDataNonConformingHexaElement.getIDOfNeighborElement(elemID, iFace, iNeib );
					for( int iEdgeNeib = 0; iEdgeNeib < 4; ++iEdgeNeib ){
						const int edgeIDLocalNeib = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFaceNeib, iEdgeNeib );
						dofSlaves[iNeib][iEdgeNeib] = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDNeib][edgeIDLocalNeib] ];
					}
				}
				assert( dofSlaves[0][1] == dofSlaves[2][0] );
				assert( dofSlaves[1][1] == dofSlaves[3][0] );
				assert( dofSlaves[0][3] == dofSlaves[1][2] );
				assert( dofSlaves[2][3] == dofSlaves[3][2] );
				if( dofMaster[0] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][0], dofMaster[0], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[1][0], dofMaster[0], 1.0 );
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][1], dofMaster[0], 0.5 );
					addMasterDofAndFactorPair( dofSlaves[1][1], dofMaster[0], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				if( dofMaster[1] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[2][1], dofMaster[1], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[3][1], dofMaster[1], 1.0 );
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][1], dofMaster[1], 0.5 );
					addMasterDofAndFactorPair( dofSlaves[1][1], dofMaster[1], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				if( dofMaster[2] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][2], dofMaster[2], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[2][2], dofMaster[2], 1.0 );
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][3], dofMaster[2], 0.5 );
					addMasterDofAndFactorPair( dofSlaves[2][3], dofMaster[2], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				if( dofMaster[3] >= 0 ){
					// Slave dofs on outer edges of the face
					addMasterDofAndFactorPair( dofSlaves[1][3], dofMaster[3], 1.0 );
					addMasterDofAndFactorPair( dofSlaves[3][3], dofMaster[3], 1.0 );
					// Slave dofs on interior edges of the face
					addMasterDofAndFactorPair( dofSlaves[0][3], dofMaster[3], 0.5 );
					addMasterDofAndFactorPair( dofSlaves[2][3], dofMaster[3], 0.5 );
				}else{
					// In this case, MPC constants are calcunated in calcMPCConstants
				}
				// Slave dofs on outer edges of the face
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][0] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[1][0] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[2][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[3][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][2] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[2][2] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[1][3] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[3][3] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_OUTER_EDGES;
				// Slave dofs on interior edges of the face
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_INTERIOR_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[1][1] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_INTERIOR_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[0][3] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_INTERIOR_EDGES;
				m_IDsAfterDegenerated2AfterConstrained[ dofSlaves[2][3] ] = Forward3DNonConformingHexaElement0thOrder::SLAVE_ON_INTERIOR_EDGES;
			}else{
				OutputFiles::m_logFile << "Error : Number of neighbor elements is wrong :  " << numNeibElements << std::endl;
				exit(1);
			}
		}
	}

//#ifdef _DEBUG_WRITE
//	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
//		std::vector< std::pair<int,double> >& vec = m_slaveDofToMasterDofAndFactors[i];
//		const int numMasters = vec.size();
//		assert( numMasters == 1 || numMasters == 2 );
//		double sumFactor(0.0);
//		for( std::vector< std::pair<int,double> >::const_iterator itrVec = vec.begin(); itrVec != vec.end(); ++itrVec ){
//			sumFactor += itrVec->second;
//		}
//		assert( std::abs(sumFactor - 0.5) < m_eps );
//	}
//#endif

	int icount(-1);
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		if( m_IDsAfterDegenerated2AfterConstrained[i] < 0 ){
			// Slave dof
			continue;
		}
		// Master dof
		addMasterDofAndFactorPair( i, i, 1.0 );// Even for master dof, dof and factor are inserted
		m_IDsAfterDegenerated2AfterConstrained[i] = ++icount;
	}
	m_numOfEquationDegeneratedAndConstrained = icount + 1;

	m_hasMadeMapSlaveDofToMasterDofAndFactors = true;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		std::vector< std::pair<int,double> >& vec = m_slaveDofToMasterDofAndFactors[i];
		for( std::vector< std::pair<int,double> >::const_iterator itrVec = vec.begin(); itrVec != vec.end(); ++itrVec ){
			std::cout << "slave master factor : " << i << " " << itrVec->first << " " << itrVec->second << std::endl;
		}
	}
#endif

}

// Calculate MPC constants
void Forward3DNonConformingHexaElement0thOrder::calcMPCConstants(){

	const int iPol = getPolarizationCurrent();// Constraint matrix does not depend on the type of polarization

	if( m_vectorMPCConstants == NULL ){
		m_vectorMPCConstants = new std::complex<double>[m_numOfEquationDegenerated];
	}
	// Zero clear because values are added in the following procedure
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		m_vectorMPCConstants[i] = std::complex<double>(0.0, 0.0);
	}

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = iElem;
		for( int iFace = 0; iFace < 6; ++iFace ){
			if( !m_MeshDataNonConformingHexaElement.faceSlaveElements(elemID, iFace) ){
				// This face does not have slave edges
				continue;
			}
			// Dofs of master edges
			int dofMasterBeforeDegenerated[4] = { -1, -1, -1, -1 };
			int dofMaster[4] = { -1, -1, -1, -1 };
			for( int iEdge = 0; iEdge < 4; ++iEdge ){
				const int edgeIDLocal = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFace, iEdge );
				const int dofBeforeDegenerated = m_IDsLocal2Global[elemID][edgeIDLocal];
				dofMasterBeforeDegenerated[iEdge] = dofBeforeDegenerated;
				dofMaster[iEdge] = m_IDsGlobal2AfterDegenerated[iPol][ dofBeforeDegenerated ];
			}
			// Face index of neighbor element
			const int iFaceNeib = m_MeshDataNonConformingHexaElement.getFaceIndexOfNeighborElement(iFace);
			const int numNeibElements = m_MeshDataNonConformingHexaElement.getNumNeighborElement(elemID, iFace);
			if( numNeibElements == 2 ){
				// Dofs of slave edges
				int dofSlaves[2][4] = { { -1, -1, -1, -1 }, { -1, -1, -1, -1 } };
				for( int iNeib = 0; iNeib < 2; ++iNeib ){
					const int elemIDNeib = m_MeshDataNonConformingHexaElement.getIDOfNeighborElement(elemID, iFace, iNeib );
					for( int iEdgeNeib = 0; iEdgeNeib < 4; ++iEdgeNeib ){
						const int edgeIDLocalNeib = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFaceNeib, iEdgeNeib );
						dofSlaves[iNeib][iEdgeNeib] = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDNeib][edgeIDLocalNeib] ];
					}
				}
				assert( dofSlaves[0][3] == dofSlaves[1][2] );
				if( dofMaster[2] < 0 ){
					if( dofMaster[2] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[2] ];
						m_vectorMPCConstants[ dofSlaves[0][3] ] += value;
					}
				}
				if( dofMaster[3] < 0 ){
					if( dofMaster[3] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[3] ];
						m_vectorMPCConstants[ dofSlaves[0][3] ] += value;
					}
				}
			}else if( numNeibElements == 4 ){
				// Dofs of slave edges
				int dofSlaves[4][4] = { { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 }, { -1, -1, -1, -1 } };
				for( int iNeib = 0; iNeib < 4; ++iNeib ){
					const int elemIDNeib = m_MeshDataNonConformingHexaElement.getIDOfNeighborElement(elemID, iFace, iNeib );
					for( int iEdgeNeib = 0; iEdgeNeib < 4; ++iEdgeNeib ){
						const int edgeIDLocalNeib = m_MeshDataNonConformingHexaElement.getEdgeIDLocalFromFaceIDLocal( iFaceNeib, iEdgeNeib );
						dofSlaves[iNeib][iEdgeNeib] = m_IDsGlobal2AfterDegenerated[iPol][ m_IDsLocal2Global[elemIDNeib][edgeIDLocalNeib] ];
					}
				}
				assert( dofSlaves[0][1] == dofSlaves[2][0] );
				assert( dofSlaves[1][1] == dofSlaves[3][0] );
				assert( dofSlaves[0][3] == dofSlaves[1][2] );
				assert( dofSlaves[2][3] == dofSlaves[3][2] );
				if( dofMaster[0] < 0 ){
					if( dofMaster[0] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[0] ];
						m_vectorMPCConstants[ dofSlaves[0][1] ] += value;
						m_vectorMPCConstants[ dofSlaves[1][1] ] += value;
					}
				}
				if( dofMaster[1] < 0 ){
					if( dofMaster[1] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[1] ];
						m_vectorMPCConstants[ dofSlaves[0][1] ] += value;
						m_vectorMPCConstants[ dofSlaves[1][1] ] += value;
					}
				}
				if( dofMaster[2] < 0 ){
					if( dofMaster[2] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[2] ];
						m_vectorMPCConstants[ dofSlaves[0][3] ] += value;
						m_vectorMPCConstants[ dofSlaves[2][3] ] += value;
					}
				}
				if( dofMaster[3] < 0 ){
					if( dofMaster[3] == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> value = std::complex<double>(0.5, 0.0) * m_globalID2NonZeroValues[ dofMasterBeforeDegenerated[3] ];
						m_vectorMPCConstants[ dofSlaves[0][3] ] += value;
						m_vectorMPCConstants[ dofSlaves[2][3] ] += value;
					}
				}
			}else{
				OutputFiles::m_logFile << "Error : Number of neighbor elements is wrong :  " << numNeibElements << std::endl;
				exit(1);
			}
		}
	}

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numOfEquationDegenerated; ++i ){
		std::cout << "m_vectorMPCConstants[" << i << "] = " << m_vectorMPCConstants[i] << std::endl;
	}
#endif

}

// Add master dof and factor pair to m_slaveDofToMasterDofAndFactors
void Forward3DNonConformingHexaElement0thOrder::addMasterDofAndFactorPair( const int slaveDof, const int masterDof, const double factor  ){

	//std::map< int, std::vector< std::pair<int,double> > >::iterator itr = m_slaveDofToMasterDofAndFactors.find(slaveDof);
	//if( itr == m_slaveDofToMasterDofAndFactors.end() ){
	//	// Pair has not been inserted
	//	std::pair<int,double> pair = std::make_pair( masterDof, factor );
	//	std::vector< std::pair<int,double> >  vec;
	//	vec.push_back(pair);
	//	m_slaveDofToMasterDofAndFactors.insert( std::make_pair( slaveDof, vec ) );
	//}else{
	//	// Pair has already been inserted
	//	std::vector< std::pair<int,double> >& vec = itr->second;
	//	bool found(false);
	//	for( std::vector< std::pair<int,double> >::iterator itrVec = vec.begin(); itrVec != vec.end(); ++itrVec ){
	//		if(itrVec->first == masterDof){
	//			found = true;
	//			break;
	//		}
	//	}
	//	if( !found ){
	//		// Insert only if the master has not been found
	//		itr->second.push_back( std::make_pair( masterDof, factor ) );
	//	}
	//}
	std::vector< std::pair<int,double> >& vec = m_slaveDofToMasterDofAndFactors[slaveDof];;
	bool found(false);
	for( std::vector< std::pair<int,double> >::iterator itrVec = vec.begin(); itrVec != vec.end(); ++itrVec ){
		if(itrVec->first == masterDof){
			found = true;
			break;
		}
	}
	if( !found ){
		// Insert only if the master has not been found
		vec.push_back( std::make_pair( masterDof, factor ) );
	}

}

// Get shape functions of the x direction with respect to the 1st reference element coordinate system
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncX( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			return 0.125 * ( 1.0 + m_etaAtEdge[num] * eta ) * ( 1.0 + m_zetaAtEdge[num] * zeta ) * invJacobMat.mat11;
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			return 0.125 * ( 1.0 + m_zetaAtEdge[num] * zeta ) * ( 1.0 + m_xiAtEdge[num] * xi ) * invJacobMat.mat12;
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			return 0.125 * ( 1.0 + m_xiAtEdge[num] * xi ) * ( 1.0 + m_etaAtEdge[num] * eta ) * invJacobMat.mat13;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of the y direction with respect to the reference element coordinate system
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncY( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			return 0.125 * ( 1.0 + m_etaAtEdge[num] * eta ) * ( 1.0 + m_zetaAtEdge[num] * zeta ) * invJacobMat.mat21;
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			return 0.125 * ( 1.0 + m_zetaAtEdge[num] * zeta ) * ( 1.0 + m_xiAtEdge[num] * xi ) * invJacobMat.mat22;
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			return 0.125 * ( 1.0 + m_xiAtEdge[num] * xi ) * ( 1.0 + m_etaAtEdge[num] * eta ) * invJacobMat.mat23;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of the z direction with respect to the reference element coordinate system
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncZ( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			return 0.125 * ( 1.0 + m_etaAtEdge[num] * eta ) * ( 1.0 + m_zetaAtEdge[num] * zeta ) * invJacobMat.mat31;
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			return 0.125 * ( 1.0 + m_zetaAtEdge[num] * zeta ) * ( 1.0 + m_xiAtEdge[num] * xi ) * invJacobMat.mat32;
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			return 0.125 * ( 1.0 + m_xiAtEdge[num] * xi ) * ( 1.0 + m_etaAtEdge[num] * eta ) * invJacobMat.mat33;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get x component of shape function rotated for 0th order edge-based elements
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncRotatedX( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	double tmp1(0.0);
	double tmp2(0.0);

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			tmp1 = 0.125 * m_etaAtEdge[num]  * ( 1.0 + m_zetaAtEdge[num] * zeta )*( invJacobMat.mat31*invJacobMat.mat22 - invJacobMat.mat21*invJacobMat.mat32 );
			tmp2 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_etaAtEdge[num]  *  eta )*( invJacobMat.mat31*invJacobMat.mat23 - invJacobMat.mat21*invJacobMat.mat33 );
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			tmp1 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]    * xi   )*( invJacobMat.mat32*invJacobMat.mat23 - invJacobMat.mat22*invJacobMat.mat33 );
			tmp2 = 0.125 * m_xiAtEdge[num]   * ( 1.0 + m_zetaAtEdge[num]  * zeta )*( invJacobMat.mat32*invJacobMat.mat21 - invJacobMat.mat22*invJacobMat.mat31 );
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			tmp1 = 0.125 * m_xiAtEdge[num]  * ( 1.0 + m_etaAtEdge[num] * eta )*( invJacobMat.mat33*invJacobMat.mat21 - invJacobMat.mat23*invJacobMat.mat31 );
			tmp2 = 0.125 * m_etaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]  * xi  )*( invJacobMat.mat33*invJacobMat.mat22 - invJacobMat.mat23*invJacobMat.mat32 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}

	return tmp1 + tmp2;

}

// Get y component of shape function rotated for 0th order edge-based elements
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncRotatedY( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	double tmp1(0.0);
	double tmp2(0.0);

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			tmp1 = 0.125 * m_etaAtEdge[num]  * ( 1.0 + m_zetaAtEdge[num] * zeta )*( invJacobMat.mat11*invJacobMat.mat32 - invJacobMat.mat31*invJacobMat.mat12 );
			tmp2 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_etaAtEdge[num]  *  eta )*( invJacobMat.mat11*invJacobMat.mat33 - invJacobMat.mat31*invJacobMat.mat13 );
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			tmp1 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]    * xi   )*( invJacobMat.mat12*invJacobMat.mat33 - invJacobMat.mat32*invJacobMat.mat13 );
			tmp2 = 0.125 * m_xiAtEdge[num]   * ( 1.0 + m_zetaAtEdge[num]  * zeta )*( invJacobMat.mat12*invJacobMat.mat31 - invJacobMat.mat32*invJacobMat.mat11 );
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			tmp1 = 0.125 * m_xiAtEdge[num]  * ( 1.0 + m_etaAtEdge[num] * eta )*( invJacobMat.mat13*invJacobMat.mat31 - invJacobMat.mat33*invJacobMat.mat11 );
			tmp2 = 0.125 * m_etaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]  * xi  )*( invJacobMat.mat13*invJacobMat.mat32 - invJacobMat.mat33*invJacobMat.mat12 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}

	return tmp1 + tmp2;

}

// Get z component of shape function rotated for 0th order edge-based elements
double Forward3DNonConformingHexaElement0thOrder::getShapeFuncRotatedZ( const double xi, const double eta, const double zeta, const int num, const Forward3D::Matrix3x3& invJacobMat ) const{

	double tmp1(0.0);
	double tmp2(0.0);

	switch( num ){
		case 0:// Go through
		case 1:// Go through
		case 2:// Go through
		case 3:
			tmp1 = 0.125 * m_etaAtEdge[num]  * ( 1.0 + m_zetaAtEdge[num] * zeta )*( invJacobMat.mat21*invJacobMat.mat12 - invJacobMat.mat11*invJacobMat.mat22 );
			tmp2 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_etaAtEdge[num]  *  eta )*( invJacobMat.mat21*invJacobMat.mat13 - invJacobMat.mat11*invJacobMat.mat23 );
			break;
		case 4:// Go through
		case 5:// Go through
		case 6:// Go through
		case 7:
			tmp1 = 0.125 * m_zetaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]    * xi   )*( invJacobMat.mat22*invJacobMat.mat13 - invJacobMat.mat12*invJacobMat.mat23 );
			tmp2 = 0.125 * m_xiAtEdge[num]   * ( 1.0 + m_zetaAtEdge[num]  * zeta )*( invJacobMat.mat22*invJacobMat.mat11 - invJacobMat.mat12*invJacobMat.mat21 );
			break;
		case 8:// Go through
		case 9:// Go through
		case 10:// Go through
		case 11:
			tmp1 = 0.125 * m_xiAtEdge[num]  * ( 1.0 + m_etaAtEdge[num] * eta )*( invJacobMat.mat23*invJacobMat.mat11 - invJacobMat.mat13*invJacobMat.mat21 );
			tmp2 = 0.125 * m_etaAtEdge[num] * ( 1.0 + m_xiAtEdge[num]  * xi  )*( invJacobMat.mat23*invJacobMat.mat12 - invJacobMat.mat13*invJacobMat.mat22 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in " << __FUNCTION__ << " : num = " << num << std::endl;
			exit(1);
			break;
	}
	
	return tmp1 + tmp2;
}

// Calculate jacobian matrix of the elements
double Forward3DNonConformingHexaElement0thOrder::calcJacobianMatrix( const int elemID,  const double xi, const double eta, const double zeta, Forward3D::Matrix3x3& JacobMat ) const{

	double xCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double yCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double zCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for( int i = 0; i < 8; ++i ){
		const int nodeID = m_MeshDataNonConformingHexaElement.getNodesOfElements(elemID, i);
		xCoord[i] = m_MeshDataNonConformingHexaElement.getXCoordinatesOfNodes(nodeID);
		yCoord[i] = m_MeshDataNonConformingHexaElement.getYCoordinatesOfNodes(nodeID);
		zCoord[i] = m_MeshDataNonConformingHexaElement.getZCoordinatesOfNodes(nodeID);
	}

	// Zero clear
	JacobMat.mat11 = 0.0;
	JacobMat.mat12 = 0.0;
	JacobMat.mat13 = 0.0;
	JacobMat.mat21 = 0.0;
	JacobMat.mat22 = 0.0;
	JacobMat.mat23 = 0.0;
	JacobMat.mat31 = 0.0;
	JacobMat.mat32 = 0.0;
	JacobMat.mat33 = 0.0;
	for( int i = 0; i < 8; ++i ){
		const double xiNode   = m_xiAtNode[i];
		const double etaNode  = m_etaAtNode[i];
		const double zetaNode = m_zetaAtNode[i];
		const double tmp1 = 0.125 * xiNode   * (1.0 + etaNode  * eta)  * (1.0 + zetaNode * zeta);
		const double tmp2 = 0.125 * etaNode  * (1.0 + zetaNode * zeta) * (1.0 + xiNode * xi);
		const double tmp3 = 0.125 * zetaNode * (1.0 + xiNode * xi)     * (1.0 + etaNode * eta);
		JacobMat.mat11 += tmp1 * xCoord[i];
		JacobMat.mat12 += tmp1 * yCoord[i];
		JacobMat.mat13 += tmp1 * zCoord[i];
		JacobMat.mat21 += tmp2 * xCoord[i];
		JacobMat.mat22 += tmp2 * yCoord[i];
		JacobMat.mat23 += tmp2 * zCoord[i];
		JacobMat.mat31 += tmp3 * xCoord[i];
		JacobMat.mat32 += tmp3 * yCoord[i];
		JacobMat.mat33 += tmp3 * zCoord[i]; 
	}

	const double determinant = JacobMat.mat11 * JacobMat.mat22 * JacobMat.mat33
							 + JacobMat.mat12 * JacobMat.mat23 * JacobMat.mat31
							 + JacobMat.mat13 * JacobMat.mat21 * JacobMat.mat32
							 - JacobMat.mat13 * JacobMat.mat22 * JacobMat.mat31
							 - JacobMat.mat12 * JacobMat.mat21 * JacobMat.mat33
							 - JacobMat.mat11 * JacobMat.mat23 * JacobMat.mat32;

	return determinant;

}

// Calculate inverse of jacobian matrix  multiplied by determinant
void Forward3DNonConformingHexaElement0thOrder::calcInverseOfJacobianMatrix( const Forward3D::Matrix3x3& jacobMat, const double determinant, Forward3D::Matrix3x3& invJacobMat ) const{

	const double invDet = 1.0 / determinant;

	invJacobMat.mat11 = ( jacobMat.mat22 * jacobMat.mat33 - jacobMat.mat23 * jacobMat.mat32 ) * invDet; 
	invJacobMat.mat12 = ( jacobMat.mat13 * jacobMat.mat32 - jacobMat.mat12 * jacobMat.mat33 ) * invDet; 
	invJacobMat.mat13 = ( jacobMat.mat12 * jacobMat.mat23 - jacobMat.mat13 * jacobMat.mat22 ) * invDet; 

	invJacobMat.mat21 = ( jacobMat.mat23 * jacobMat.mat31 - jacobMat.mat21 * jacobMat.mat33 ) * invDet; 
	invJacobMat.mat22 = ( jacobMat.mat11 * jacobMat.mat33 - jacobMat.mat13 * jacobMat.mat31 ) * invDet; 
	invJacobMat.mat23 = ( jacobMat.mat13 * jacobMat.mat21 - jacobMat.mat11 * jacobMat.mat23 ) * invDet; 

	invJacobMat.mat31 = ( jacobMat.mat21 * jacobMat.mat32 - jacobMat.mat22 * jacobMat.mat31 ) * invDet; 
	invJacobMat.mat32 = ( jacobMat.mat12 * jacobMat.mat31 - jacobMat.mat11 * jacobMat.mat32 ) * invDet; 
	invJacobMat.mat33 = ( jacobMat.mat11 * jacobMat.mat22 - jacobMat.mat12 * jacobMat.mat21 ) * invDet; 

}

// Output results of forward calculation to VTK file
void Forward3DNonConformingHexaElement0thOrder::outputResultToVTK() const{

	if( !OutputFiles::m_vtkFile.is_open() ){
		return;
	}
	
	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();

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
	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();

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
void Forward3DNonConformingHexaElement0thOrder::outputResultToBinary( const int iFreq, const int iPol ) const{

	const std::string stringPolarization = iPol == CommonParameters::EX_POLARIZATION ? "ExPol" : "EyPol";
	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();

	//--- Total element number
	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
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
		const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
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
		const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
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
		const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
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

// Add values to right-hand sides matrix consisting of interpolator vectors by taking into account MPC
void Forward3DNonConformingHexaElement0thOrder::addValuesToRhsVectorsByConsideringMPC( const int irow, const int irhs, const std::complex<double>& val ){

	const std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[irow];
	for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
		const int iRowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
		const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
		addValuesToRhsVectors( iRowMod, irhs, valMod );
	}

}

//// Make pair of master dof and factor for a slave dof
//void Forward3DNonConformingHexaElement0thOrder::makePairOfMasterDofAndFactorForASlaveDof( const int slaveDofAfterDegenerated, std::vector< std::pair<int,double> >& masterDofAndFactorsAfterDegenerated ) const{
//
//	std::map< int, std::vector< std::pair<int,double> > > masterDofBeforeDegeneratedAndFactorsAll;
//
//	const int slaveDofBeforeDegenerated = getDofBeforeDegenerationFromDofAfterDegeneration(slaveDofAfterDegenerated);
//
//	std::vector< std::pair<int,double> >& masterDofBeforeDegeneratedAndFactors = masterDofBeforeDegeneratedAndFactorsAll[slaveDofBeforeDegenerated];
//	if( static_cast<int>( masterDofBeforeDegeneratedAndFactors.size() ) == 0 ){
//		// This dof is master 
//		masterDofAndFactorsAfterDegenerated.push_back( std::make_pair(slaveDofAfterDegenerated, 1.0) );
//		return;
//	}
//
//	// This dof depends on some master dofs
//	std::vector< std::pair<int,double> > masterDofAndFactorsBeforeDegenerated;
//	//for( std::vector< std::pair<int,double> >::const_iterator itr = masterDofBeforeDegeneratedAndFactors.begin(); itr != masterDofBeforeDegeneratedAndFactors.end(); ++itr ){
//	//	const int masterDofBeforeDegenerated = itr->first;
//	//	const double factor = itr->second;
//	//	masterDofAndFactorsBeforeDegenerated.push_back( std::make_pair(masterDofBeforeDegenerated, factor) );
//	//}
//	masterDofAndFactorsBeforeDegenerated.push_back( std::make_pair(slaveDofBeforeDegenerated, 1.0) );
//
//	bool found(true);
//	while (found){
//		found = false;		
//		for( std::vector< std::pair<int,double> >::iterator itr = masterDofAndFactorsBeforeDegenerated.begin(); itr != masterDofAndFactorsBeforeDegenerated.end(); ++itr ){
//			const int masterDofBeforeDegeneratedOrg = itr->first;
//			std::vector< std::pair<int,double> > masterDofAndFactorsBeforeDegeneratedAux;
//			makePairOfMasterDofAndFactorForASlaveDofAux( masterDofBeforeDegeneratedOrg, masterDofAndFactorsBeforeDegeneratedAux );
//			if( static_cast<int>( masterDofAndFactorsBeforeDegeneratedAux.size() ) == 0 ){
//				// The dof does not has masters
//				continue;
//			}
//			// The dofhas masters
//			found = true;
//			const double factorOrg = itr->second;
//			std::vector< std::pair<int,double> >::const_iterator itrSub = masterDofAndFactorsBeforeDegeneratedAux.begin();
//			itr->first = itrSub->first;// Replace dof
//			itr->second = factorOrg * itrSub->second;// Replace factor
//			for( ; itrSub != masterDofAndFactorsBeforeDegeneratedAux.end(); ++itrSub ){
//				const int masterDofBeforeDegenerated = itrSub->first;
//				const double factor = factorOrg * itrSub->second;
//				masterDofAndFactorsBeforeDegenerated.push_back( std::make_pair(masterDofBeforeDegenerated, factor) );
//			}
//		}
//	}
//
//	const int iPol = CommonParameters::EX_POLARIZATION;// Constraint matrix does not depend on the type of polarization
//	for( std::vector< std::pair<int,double> >::const_iterator itr = masterDofAndFactorsBeforeDegenerated.begin(); itr != masterDofAndFactorsBeforeDegenerated.end(); ++itr ){
//		const int dofBeforeDegenerated = itr->first;
//		const double factor = itr->second;
//		//std::map<int,int>::const_iterator itrBeforeToAfter = arrayDofBeforeDegeneratedToAfterDegenerated.find(dofBeforeDegenerated);
//		//if( itrBeforeToAfter == arrayDofBeforeDegeneratedToAfterDegenerated.end() ){
//		//	OutputFiles::m_logFile << "Error : Key " << dofBeforeDegenerated << " is not found in arrayDofBeforeDegeneratedToAfterDegenerated." << std::endl;
//		//	exit(1);
//		//}
//		//masterDofAndFactorsAfterDegenerated.push_back( std::make_pair( itrBeforeToAfter->second, factor ) );
//		const int dofAfterDegenerated = m_IDsGlobal2AfterDegenerated[iPol][dofBeforeDegenerated];
//		masterDofAndFactorsAfterDegenerated.push_back( std::make_pair( dofAfterDegenerated, factor ) );
//	}
//
//}
//
//void Forward3DNonConformingHexaElement0thOrder::makePairOfMasterDofAndFactorForASlaveDofAux( const int slaveDofBeforeDegenerated, std::vector< std::pair<int,double> >& masterDofAndFactorsBeforeDegenerated ) const{
//
//	std::map< int, std::vector< std::pair<int,double> > > masterDofBeforeDegeneratedAndFactorsAll;
//
//	std::vector< std::pair<int,double> >& masterDofBeforeDegeneratedAndFactors = masterDofBeforeDegeneratedAndFactorsAll[slaveDofBeforeDegenerated];
//	if( static_cast<int>( masterDofBeforeDegeneratedAndFactors.size() ) == 0 ){
//		return;
//	}
//
//	// This dof depends on some master dofs
//	for( std::vector< std::pair<int,double> >::const_iterator itr = masterDofBeforeDegeneratedAndFactors.begin(); itr != masterDofBeforeDegeneratedAndFactors.end(); ++itr ){
//		const int masterDofBeforeDegenerated = itr->first;
//		const double factor = itr->second;
//		masterDofAndFactorsBeforeDegenerated.push_back( std::make_pair(masterDofBeforeDegenerated, factor) );
//	}
//
//}
//
//// Get dof before degeneration from the after degeneration
//int Forward3DNonConformingHexaElement0thOrder::getDofBeforeDegenerationFromDofAfterDegeneration( const int dofAfterDegeneration ) const{
//
//
//}