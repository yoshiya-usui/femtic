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
#include "Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "ResistivityBlock.h"
#include <algorithm>
#include <assert.h>

// Constructer
Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased( const int planeID, const int iPol ):
	Forward2DQuadrilateralElementEdgeBased( planeID, iPol ),
	m_slaveDofToMasterDofAndFactors(NULL),
	m_hasMadeMapSlaveDofToMasterDofAndFactors(false),
	m_IDsAfterDegenerated2AfterConstrained(NULL),
	m_numEquationDegeneratedAndConstrained(0)
{
	//---------------------------------------------------------------------------
	//--- Calculate integral points and weights of two point Gauss quadrature ---
	//---------------------------------------------------------------------------
	int ip(0);
	for( int i = 0; i < m_numGauss; ++i ){
		for( int j = 0; j < m_numGauss; ++j ){
			m_integralPointXi[ip]   = CommonParameters::abscissas2Point[i];
			m_integralPointEta[ip]  = CommonParameters::abscissas2Point[j];
			m_weights[ip] = CommonParameters::weights2Point[i] * CommonParameters::weights2Point[j];
			++ip;
		}
	}

	// Array of reference coord xi values for each node
	m_xiAtNode[0] = -1.0;
	m_xiAtNode[1] =  1.0;
	m_xiAtNode[2] =  1.0;
	m_xiAtNode[3] = -1.0;

	// Array of reference coord eta values for each node
	m_etaAtNode[0] = -1.0;
	m_etaAtNode[1] = -1.0;
	m_etaAtNode[2] =  1.0;
	m_etaAtNode[3] =  1.0;

	// Array of reference coord xi values for each edge
	m_xiAtEdge[0]  = -9999.999;
	m_xiAtEdge[1]  = -9999.999;
	m_xiAtEdge[2]  = -1.0;
	m_xiAtEdge[3]  =  1.0;

	// Array of reference coord eta values for each edge
	m_etaAtEdge[0]  = -1.0;
	m_etaAtEdge[1]  =  1.0;
	m_etaAtEdge[2]  = -9999.999;
	m_etaAtEdge[3]  = -9999.999;

}

// Destructer
Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::~Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased(){
}

// Calculate EM fields of boundary planes by 2D forward calculcation with 0th order edge element
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataNonConformingHexaElement* const pMeshData ){

	const int imode = calcMode();// TM or TE mode
	if ( imode != CommonParameters::TM_MODE ){
		OutputFiles::m_logFile << "Error : Only TM mode can be treated in Forward2DTriangleElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes ! imode = " << imode << "." << std::endl;
		exit(1);
	}
	if( !m_sourceFieldElectric ){
		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in Forward2DTriangleElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes !" << std::endl;
		exit(1);
	}

	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );

	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();
	if( pAnalysisControl->getBoundaryConditionBottom() != AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		OutputFiles::m_logFile << "Error : When the 0th order edge-based element is used, electric field of the bottom must be zero !" << std::endl;
		exit(1);
	}
	if( m_specifyTEResultToSidesOfEdgeElement ){
		OutputFiles::m_logFile << "Error : Horizontal electric field cannot be specified at the left and right side the 2D model of the 0th order edge-based element !" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Calculate electric field on a boundary plane with 0th order edge-based element." << std::endl;

	//-----------------------------------------------------------------------------
	//--- Set Number of equations and array converting local IDs to global ones ---
	//-----------------------------------------------------------------------------
	if( !m_hasAlreadySetIDsLocal2Global ){
		calcArrayConvertLocalID2Global( pMeshData );
	}

	if( !m_hasMadeMapSlaveDofToMasterDofAndFactors ){
		makeMapSlaveDofToMasterDofAndFactors( pMeshData );
	}

	OutputFiles::m_logFile << "# Number of equation = " << m_numEquations
		<< ", Number of equation after degeneration = " << m_numEquationsDegenerated
		<< ", Number of equation after constraint = " << m_numEquationDegeneratedAndConstrained << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() <= 0 ){
		m_matrix2DAnalysis.setDegreeOfEquation(m_numEquationDegeneratedAndConstrained);
	}

#ifdef _DEBUG_WRITE
	std::cout << "m_numEquations = " << m_numEquations << std::endl;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem m_IDsLocal2Global : " << iElem;
		for( int i = 0; i < 4; ++i ){
			std::cout << " " << m_IDsLocal2Global[iElem][i];
		}
		std::cout << std::endl;
	}
	std::cout << "m_numEquationsDegenerated = " << m_numEquationsDegenerated << std::endl;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
		for( int i = 0; i < 4; ++i ){
			std::cout << " " << m_IDsLocal2GlobalDegenerated[iElem][i];
		}
		std::cout << std::endl;
	}
#endif

	//-----------------------------------------------------------------------
	//--- Set matrix structure and analyze the structure by matrix solver ---
	//-----------------------------------------------------------------------
	if( m_hasMatrixStructureSetAndAnalyzed == false ){
		//// If matrix structure has not been set yet, set matrix structure.
		//OutputFiles::m_logFile << "# Set matrix structure. " << pAnalysisControl->outputElapsedTime() << std::endl;

		////------------------------------------------
		////--- Components due to stiffness matrix ---
		////------------------------------------------
		//for( int iElem = 0; iElem < nElem; ++iElem ){
		//	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){
		//		const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
		//		if( row < 0 ){
		//			continue;
		//		}
		//		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){
		//			const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
		//			if( col < 0 ){
		//				continue;
		//			}
		//			if( col >= row ){// Store only upper triangle part
		//				m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
		//			}
		//		}// iEdge2
		//	}// iEdge1
		//}// iElem]
		//--- Set matrix structure and analyze the structure by matrix solver
		setNonZeroStrucuture(pMeshData);

		//--- Convert the matrix from the triplet format to the CRS format
		if( m_matrix2DAnalysis.hasConvertedToCRSFormat() == false ){
			m_matrix2DAnalysis.convertToCRSFormat();
		}

		//--- Anaysis phase of matrix solver
		OutputFiles::m_logFile << "# Anaysis phase of matrix solver for boundary plane."
			<<  " Polarization : " << m_polarization << " Plane : " << m_planeID << ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
		m_matrix2DAnalysis.analysisPhaseMatrixSolver();
		m_hasMatrixStructureSetAndAnalyzed = true;
	}

	//-------------------------------------------------------
	//--- Set values of matrix and right hand side vector ---
	//-------------------------------------------------------
	OutputFiles::m_logFile << "# Set values of matrix and right hand side vector. " << pAnalysisControl->outputElapsedTime() << std::endl;
	//const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	m_matrix2DAnalysis.zeroClearNonZeroValues();// Zero clear matrix values
	m_matrix2DAnalysis.zeroClearRightHandSideVector();// Zero clear right hand side vector

	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	//--- Calculate omega * mu * sigma
	//	const int elemID = pMeshData->getElemBoundaryPlanes( m_planeID, iElem );
	//	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
	//	const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
	//	const std::complex<double> factor = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form

	//	const double length[4] = { pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 0 ),
	//							   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 1 ),
	//							   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 2 ),
	//							   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 3 ) };

	//	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){
	//		const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
	//		if( row < 0 ){
	//			continue;
	//		}
	//		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){
	//			const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
	//			if( col <= Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE ){
	//				continue;
	//			}
	//			double integral1(0.0);
	//			double integral2(0.0);
	//			for( int ip = 0; ip < m_numIntegralPoints; ++ip ){
	//				const double xi = m_integralPointXi[ip];
	//				const double eta = m_integralPointEta[ip];
	//				Forward2D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	//				const double detJacob = calcJacobianMatrix( pMeshData, iElem, xi, eta, jacobMat );
	//				Forward2D::Matrix2x2 invJacobMat;
	//				calcInverseOfJacobianMatrix( jacobMat, detJacob, invJacobMat );
	//				integral1 += getShapeFuncRotated(xi, eta, iEdge1, invJacobMat) * getShapeFuncRotated(xi, eta, iEdge2, invJacobMat) * detJacob * m_weights[ip];
	//				integral2 += ( getShapeFuncH(xi, eta, iEdge1, invJacobMat) * getShapeFuncH(xi, eta, iEdge2, invJacobMat)
	//					         + getShapeFuncV(xi, eta, iEdge1, invJacobMat) * getShapeFuncV(xi, eta, iEdge2, invJacobMat) ) * detJacob * m_weights[ip];
	//			}
	//			integral1 *= length[iEdge1] * length[iEdge2];
	//			integral2 *= length[iEdge1] * length[iEdge2];
	//			std::complex<double> val = std::complex<double>( integral1 , 0.0 ) - std::complex<double>( integral2 , 0.0 ) * factor;// exp(-i*omega*t) form
	//			if( col == Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
	//				const std::complex<double> nonZeroValue = m_edgesIDGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge2] ];
	//				m_matrix2DAnalysis.addRightHandSideVector( row, -val * nonZeroValue );// Add to right hand side vector
	//			}else if( col >= row ){// Store only upper triangle part
	//				m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
	//			}
	//		}// iEdge2
	//	}// iEdge1		
	//}// iElem
	// 
	//--- Set non-zero values of matrix and right-hande side vector for forward calculation
	setNonZeroValues(freq, pMeshData);

#ifdef _DEBUG_WRITE
	m_matrix2DAnalysis.debugWriteMatrix();
	m_matrix2DAnalysis.debugWriteRightHandSide();
#endif

	//-----------------------------------------------------
	//--- Numrical factorization phase of matrix solver ---
	//-----------------------------------------------------
	OutputFiles::m_logFile << "# Numerical factorization phase of matrix solver for boundary plane."
			<<  " Polarization : " << m_polarization << " Plane : " << m_planeID  << ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix2DAnalysis.factorizationPhaseMatrixSolver();//Numrical factorization phase of matrix solver

	//------------------------------------
	//--- Solve phase of matrix solver ---
	//------------------------------------

	std::complex<double>* solutionConstrained = new std::complex<double>[m_numEquationDegeneratedAndConstrained];
	OutputFiles::m_logFile << "# Solve phase of matrix solver for boundary plane."
			<<  " Polarization : " << m_polarization << " Plane : " << m_planeID <<  ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix2DAnalysis.solvePhaseMatrixSolver(solutionConstrained);//Solve phase of matrix solver

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numEquationDegeneratedAndConstrained; ++i ){
		std::cout << "i solutionConstrained : " << i << " " << solutionConstrained[i] << std::endl;
	}
#endif

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[m_numEquations];

	bool* alreadyFound = new bool[ m_numEquations ];
	for( int i = 0; i < m_numEquations; ++i ){
		alreadyFound[i] = false;
	}

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 4; ++iEdge ){
			const int iNum = m_IDsLocal2Global[iElem][iEdge];
			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsLocal2GlobalDegenerated[iElem][iEdge];
				if( iNumDegenerated == Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					m_solution[iNum] = m_edgesIDGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge] ];
				}else{
					//m_solution[iNum] = solutionDegenerated[iNumDegenerated]; 
					m_solution[iNum] = std::complex<double>(0.0, 0.0);
					const std::vector< std::pair<int,double> >& masterDofAndFactor = m_slaveDofToMasterDofAndFactors[iNumDegenerated];
					for( std::vector< std::pair<int,double> >::const_iterator itr = masterDofAndFactor.begin(); itr != masterDofAndFactor.end(); ++itr ){
						const int dof = m_IDsAfterDegenerated2AfterConstrained[itr->first];
						const std::complex<double> factor = std::complex<double>(itr->second, 0.0);
						m_solution[iNum] += solutionConstrained[dof] * factor;
					}
				}
				alreadyFound[iNum] = true;
			}
		}// iEdge
	}// iElem
	
	delete[] solutionConstrained;
	delete[] alreadyFound;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numEquations; ++i ){
		std::cout << "i m_solution : " << i << " " << m_solution[i] << std::endl;
	}
#endif
	
	output2DResult( freq, pMeshData );
}

// Calculate array converting global element IDs of 2D mesh to global edge IDs of 2D mesh
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcArrayConvertLocalID2Global( const MeshDataNonConformingHexaElement* const pMeshData ){

	typedef std::pair< int, int > NodeIDPair;
	typedef std::pair< int, int > ElemAndEdgeID;
	typedef std::pair< int, std::pair< NodeIDPair, ElemAndEdgeID > > InputedDataType;

	std::vector< InputedDataType > indicatorAndLocalEdgeID;

	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 4; ++iEdge ){
			const int nodeID0 = pMeshData->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, iEdge, 0 );
			const int nodeID1 = pMeshData->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, iEdge, 1 );
			const NodeIDPair nodePairCur = nodeID1 > nodeID0 ? std::make_pair( nodeID0, nodeID1 ) : std::make_pair( nodeID1, nodeID0 );
			const int indicator = ( nodePairCur.first * nodePairCur.second );
			const ElemAndEdgeID elemAndEdgeIDCur = std::make_pair( iElem, iEdge );
			indicatorAndLocalEdgeID.push_back( InputedDataType( indicator, std::pair< NodeIDPair, ElemAndEdgeID >(nodePairCur, elemAndEdgeIDCur) ) );
		}
	}

	std::sort( indicatorAndLocalEdgeID.begin(), indicatorAndLocalEdgeID.end() );

#ifdef _DEBUG_WRITE
	for( std::vector<InputedDataType>::iterator itr = indicatorAndLocalEdgeID.begin(); itr != indicatorAndLocalEdgeID.end(); ++itr ){
		std::cout << "indicatorAndLocalEdgeID : " << (*itr).first << " " << (*itr).second.first.first << " " << (*itr).second.first.second << " " << (*itr).second.second.first << " " << (*itr).second.second.second << std::endl;
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
		m_IDsLocal2Global[iElem] = new int[4];
	}
	//-----------------------------------------

	// Allocate memory to m_IDsLocal2GlobalDegenerated ---
	if( m_IDsLocal2GlobalDegenerated != NULL ){
		const int nElem = sizeof( m_IDsLocal2GlobalDegenerated ) / sizeof( m_IDsLocal2GlobalDegenerated[0] );
		for( int iElem = 0; iElem < nElem; ++iElem ){
			delete [] m_IDsLocal2GlobalDegenerated[iElem];
		}
		delete [] m_IDsLocal2GlobalDegenerated;
		m_IDsLocal2GlobalDegenerated = NULL;
	}

	m_IDsLocal2GlobalDegenerated = new int*[nElem]; 
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_IDsLocal2GlobalDegenerated[iElem] = new int[4];
	}
	//----------------------------------------------------

	int indicatorPre(-1);
	std::pair<int, int> nodePairPre;
	int edgeIDGlobal2D(-1);
	int edgeIDGlobal2DDegenerated(-1);
	bool flag(false);

	std::vector<InputedDataType>::const_iterator itrEnd = indicatorAndLocalEdgeID.end();
	for( std::vector<InputedDataType>::const_iterator itr = indicatorAndLocalEdgeID.begin(); itr != itrEnd; ++itr ){
		const int elemIDGlobal2D = (*itr).second.second.first;
		const int edgeIDLocal2D = (*itr).second.second.second;

		const int nodePairCur0 = pMeshData->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 0 );
		const int nodePairCur1 = pMeshData->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 1 );

		if( (*itr).first == indicatorPre &&
			( ( nodePairCur0 == nodePairPre.first && nodePairCur1 == nodePairPre.second ) ||
			  ( nodePairCur0 == nodePairPre.second && nodePairCur1 == nodePairPre.first ) ) ){
			// Same edge with privious one
			flag = false;
			m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] = edgeIDGlobal2D;
			m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = edgeIDGlobal2DDegenerated;
		}else{// New edge
			if(flag){// If the previous edge is isolated
				const int elemIDGlobal2DPre = (*(itr-1)).second.second.first;
				const int edgeIDLocal2DPre = (*(itr-1)).second.second.second;
				if( isOuterEdge(elemIDGlobal2DPre, edgeIDLocal2DPre, pMeshData ) ){
					// Outer edges
					if( getTypeOfOuterEdgeOfBoundaryPlanes(edgeIDLocal2DPre) == Forward2DQuadrilateralElementEdgeBased::UPPER_EDGE ){
						m_IDsLocal2GlobalDegenerated[elemIDGlobal2DPre][edgeIDLocal2DPre] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE;
					}else{
						m_IDsLocal2GlobalDegenerated[elemIDGlobal2DPre][edgeIDLocal2DPre] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE;
					}
					--edgeIDGlobal2DDegenerated;
				}
			}
			flag = true;
			m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] = ++edgeIDGlobal2D;
			m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = ++edgeIDGlobal2DDegenerated;
			if( itr+1 == itrEnd ){// Last element 
				if( isOuterEdge(elemIDGlobal2D, edgeIDLocal2D, pMeshData ) ){
					// Outer edges
					if( getTypeOfOuterEdgeOfBoundaryPlanes(edgeIDLocal2D) == Forward2DQuadrilateralElementEdgeBased::UPPER_EDGE ){
						m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE;
					}else{
						m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE;
					}
					--edgeIDGlobal2DDegenerated;
				}
			}
			indicatorPre = (*itr).first;
			nodePairPre.first  = nodePairCur0;
			nodePairPre.second = nodePairCur1;
		}

#ifdef _DEBUG_WRITE
		std::cout << "elemIDGlobal2D edgeIDLocal2D nodePairCur0 nodePairCur1 m_IDsLocal2Global :  " 
			<< elemIDGlobal2D << " " << edgeIDLocal2D << " " << nodePairCur0 << " " << nodePairCur1 << " "
			<< m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] << std::endl;
#endif
	}

	const double sourceValueElectric = CommonParameters::sourceValueElectric;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 4; ++iEdge ){
			if( m_IDsLocal2GlobalDegenerated[iElem][iEdge] == Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
				m_edgesIDGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][iEdge], std::complex<double>(sourceValueElectric, 0.0) ) );
			}
		}
	}

#ifdef _DEBUG_WRITE
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 4; ++iEdge ){
			std::cout << "elemIDGlobal2D edgeIDLocal2D m_IDsLocal2Global m_IDsLocal2GlobalDegenerated :  " 
				<< iElem << " " << iEdge << " " << m_IDsLocal2Global[iElem][iEdge] << " "
				<< m_IDsLocal2GlobalDegenerated[iElem][iEdge] << std::endl;
		}
	}
	for( std::map< int, std::complex<double> >::iterator itr = m_edgesIDGlobal2NonZeroValues.begin(); itr != m_edgesIDGlobal2NonZeroValues.end(); ++itr ){
		std::cout << "m_edgesIDGlobal2NonZeroValues first second :  " << itr->first << " " << itr->second << std::endl;
	}
#endif

	m_numEquations = edgeIDGlobal2D + 1;
	m_numEquationsDegenerated = edgeIDGlobal2DDegenerated + 1;
	
	m_hasAlreadySetIDsLocal2Global = true;

}

// Make map converting master dofs after degeneration and MPC factors from slave dof after degeneration 
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::makeMapSlaveDofToMasterDofAndFactors( const MeshDataNonConformingHexaElement* const pMeshData ){

	if( m_IDsAfterDegenerated2AfterConstrained != NULL ){
		delete [] m_IDsAfterDegenerated2AfterConstrained;
		m_IDsAfterDegenerated2AfterConstrained = NULL;
	}
	m_IDsAfterDegenerated2AfterConstrained = new int[m_numEquationsDegenerated];
	for( int i = 0; i < m_numEquationsDegenerated; ++i ){
		m_IDsAfterDegenerated2AfterConstrained[i] = 0;
	}

	if( m_slaveDofToMasterDofAndFactors != NULL ){
		delete [] m_slaveDofToMasterDofAndFactors;
		m_slaveDofToMasterDofAndFactors = NULL;
	}
	m_slaveDofToMasterDofAndFactors = new std::vector< std::pair<int,double> >[m_numEquationsDegenerated];

	int edgeToNeibIndex[4][2] = { {-1, -1}, {-1, -1}, {-1, -1}, {-1, -1} };
	switch (m_planeID){
		case MeshData::YZMinus:
			edgeToNeibIndex[0][0] = 0;
			edgeToNeibIndex[0][1] = 2;
			edgeToNeibIndex[1][0] = 0;
			edgeToNeibIndex[1][1] = 2;
			edgeToNeibIndex[2][0] = 0;
			edgeToNeibIndex[2][1] = 2;
			edgeToNeibIndex[3][0] = 0;
			edgeToNeibIndex[3][1] = 2;
			break;
		case MeshData::YZPlus:
			edgeToNeibIndex[0][0] = 1;
			edgeToNeibIndex[0][1] = 3;
			edgeToNeibIndex[1][0] = 1;
			edgeToNeibIndex[1][1] = 3;
			edgeToNeibIndex[2][0] = 1;
			edgeToNeibIndex[2][1] = 3;
			edgeToNeibIndex[3][0] = 1;
			edgeToNeibIndex[3][1] = 3;
			break;
		case MeshData::ZXMinus:
			edgeToNeibIndex[0][0] = 0;
			edgeToNeibIndex[0][1] = 1;
			edgeToNeibIndex[1][0] = 0;
			edgeToNeibIndex[1][1] = 1;
			edgeToNeibIndex[2][0] = 0;
			edgeToNeibIndex[2][1] = 2;
			edgeToNeibIndex[3][0] = 0;
			edgeToNeibIndex[3][1] = 2;
			break;
		case MeshData::ZXPlus:
			edgeToNeibIndex[0][0] = 2;
			edgeToNeibIndex[0][1] = 3;
			edgeToNeibIndex[1][0] = 2;
			edgeToNeibIndex[1][1] = 3;
			edgeToNeibIndex[2][0] = 1;
			edgeToNeibIndex[2][1] = 3;
			edgeToNeibIndex[3][0] = 1;
			edgeToNeibIndex[3][1] = 3;
			break;
		default:
			OutputFiles::m_logFile << "Error : Plane ID is wrong : " << m_planeID << "." << std::endl;
			exit(1);		
			break;
	}

	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int elemID = pMeshData->getElemBoundaryPlanes( m_planeID, iElem );
		for( int iEdge = 0; iEdge < 4; ++iEdge ){
			const int faceID = getNeighborFaceIndexFromEdgeIndex(iEdge);
			if( !pMeshData->faceSlaveElements(elemID, faceID) ){
				// This face does not have slave edges
				continue;
			}
			if( pMeshData->getNumNeighborElement(elemID, faceID) != 4 ){
				continue;
			}
			// Dofs of master edges
			const int dofMaster = m_IDsLocal2GlobalDegenerated[iElem][iEdge];
			assert( dofMaster >= 0 );
			// Face index of neighbor element
			const int faceIDNeib = pMeshData->getFaceIndexOfNeighborElement(faceID);
			// Dofs of slave edges
			int dofSlaves[2] = { -1, -1 };
			for( int iNeib = 0; iNeib < 2; ++iNeib ){
				const int neibIndex = edgeToNeibIndex[iEdge][iNeib];
				const int elemIDNeib = pMeshData->getIDOfNeighborElement(elemID, faceID, neibIndex);
				int iEdgeNeib = -1;
				switch (iEdge){
					case 0:
						iEdgeNeib = 1;
						break;
					case 1:
						iEdgeNeib = 0;
						break;
					case 2:
						iEdgeNeib = 3;
						break;
					case 3:
						iEdgeNeib = 2;
						break;
					default:
						OutputFiles::m_logFile << "Error : Edge index is worng !! : iEdge = " << iEdge << std::endl;
						exit(1);
						break;
				}
				bool found(false);
				int iElemNeib = 0;
				for( ; iElemNeib < nElem; ++iElemNeib ){
					if( elemIDNeib == pMeshData->getElemBoundaryPlanes( m_planeID, iElemNeib ) ){
						found = true;
						break;
					}
				}
				if(!found){
					OutputFiles::m_logFile << "Error : Element " << elemIDNeib << " is not found in the boundary plane " << m_planeID << " !!" << std::endl;
					exit(1);
				}
				dofSlaves[iNeib] = m_IDsLocal2GlobalDegenerated[iElemNeib][iEdgeNeib];
			}
			addMasterDofAndFactorPair( dofSlaves[0], dofMaster, 1.0 );
			addMasterDofAndFactorPair( dofSlaves[1], dofMaster, 1.0 );
			m_IDsAfterDegenerated2AfterConstrained[dofSlaves[0]] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::SLAVE_DOFS;
			m_IDsAfterDegenerated2AfterConstrained[dofSlaves[1]] = Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::SLAVE_DOFS;
		}
	}

	int icount(-1);
	for( int i = 0; i < m_numEquationsDegenerated; ++i ){
		if( m_IDsAfterDegenerated2AfterConstrained[i] < 0 ){
			// Slave dof
			continue;
		}
		// Master dof
		addMasterDofAndFactorPair( i, i, 1.0 );// Even for master dof, dof and factor are inserted
		m_IDsAfterDegenerated2AfterConstrained[i] = ++icount;
	}
	m_numEquationDegeneratedAndConstrained = icount + 1;

	m_hasMadeMapSlaveDofToMasterDofAndFactors = true;

#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numEquationsDegenerated; ++i ){
		std::vector< std::pair<int,double> >& vec = m_slaveDofToMasterDofAndFactors[i];
		for( std::vector< std::pair<int,double> >::const_iterator itrVec = vec.begin(); itrVec != vec.end(); ++itrVec ){
			std::cout << "slave master factor : " << i << " " << itrVec->first << " " << itrVec->second << std::endl;
		}
	}
#endif

}

// Get type of outer edge
// [note] : You must confirm inputed edge is the outer edge. 
int Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::getTypeOfOuterEdgeOfBoundaryPlanes( const int edgeIDLocal2D ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//Y-Z Plane
		if( edgeIDLocal2D == 2 ){
			return Forward2DQuadrilateralElementEdgeBased::LEFT_EDGE;
		}else if( edgeIDLocal2D == 3 ){
			return Forward2DQuadrilateralElementEdgeBased::RIGHT_EDGE;
		}
	}else if( m_planeID == MeshData::ZXMinus || m_planeID == MeshData::ZXPlus ){//Z-X Plane
		if( edgeIDLocal2D == 2 ){
			return Forward2DQuadrilateralElementEdgeBased::LEFT_EDGE;
		}else if( edgeIDLocal2D == 3 ){
			return Forward2DQuadrilateralElementEdgeBased::RIGHT_EDGE;
		}
	}else{
		OutputFiles::m_logFile << "Error : Wrong plane ID !! m_planeID = " << m_planeID << std::endl;
		exit(1);
	}
	
	if( edgeIDLocal2D == 0 ){
		return Forward2DQuadrilateralElementEdgeBased::UPPER_EDGE;
	}else if( edgeIDLocal2D == 1 ){
		return Forward2DQuadrilateralElementEdgeBased::LOWER_EDGE;
	}

	OutputFiles::m_logFile << "Error : Edge " << edgeIDLocal2D << " of the 2-D element is not outer edge." << std::endl;
	exit(1);

}

// Get shape functions of the horizontal direction with respect to the reference element coordinate system
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::getShapeFuncH( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat ) const{

	switch( num ){
		case 0:// go through
		case 1:
			return 0.25 * ( 1.0 + m_etaAtEdge[num] * eta ) * invJacobMat.mat11;
			break;
		case 2:// go through
		case 3:
			return 0.25 * ( 1.0 + m_xiAtEdge[num] * xi ) * invJacobMat.mat12;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncH : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of the vertical direction with respect to the reference element coordinate system
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::getShapeFuncV( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat ) const{

	switch( num ){
		case 0:// go through
		case 1:
			return 0.25 * ( 1.0 + m_etaAtEdge[num] * eta ) * invJacobMat.mat21;
			break;
		case 2:// go through
		case 3:
			return 0.25 * ( 1.0 + m_xiAtEdge[num] * xi ) * invJacobMat.mat22;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncV : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions rotated with respect to the reference element coordinate system
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::getShapeFuncRotated( const double xi, const double eta, const int num, const Forward2D::Matrix2x2& invJacobMat ) const{

	switch( num ){
		case 0:// go through
		case 1:
			return 0.25 * m_etaAtEdge[num] * (invJacobMat.mat21*invJacobMat.mat12 - invJacobMat.mat11*invJacobMat.mat22 );
			break;
		case 2:// go through
		case 3:
			return 0.25 * m_xiAtEdge[num] * (invJacobMat.mat22*invJacobMat.mat11 - invJacobMat.mat12*invJacobMat.mat21 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotated : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Calculate horizontal electric field
std::complex<double> Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcValueElectricFieldHorizontal( const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	Forward2D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	const double detJacob = calcJacobianMatrix( pMeshData, iElem, xi, eta, jacobMat );
	Forward2D::Matrix2x2 invJacobMat;
	calcInverseOfJacobianMatrix( jacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 4; ++i ){
		const double length = pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncH( xi, eta, i, invJacobMat ) * length, 0.0 );
	}

	return val;

}

// Calculate vertical electric field
std::complex<double> Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcValueElectricFieldVertical( const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	Forward2D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	const double detJacob = calcJacobianMatrix( pMeshData, iElem, xi, eta, jacobMat );
	Forward2D::Matrix2x2 invJacobMat;
	calcInverseOfJacobianMatrix( jacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 4; ++i ){
		const double length = pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncV( xi, eta, i, invJacobMat ) * length, 0.0 );
	}

	return val;

}

// Calculate magnetic field perpendicular to the boundary plane
std::complex<double> Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double xi, const double eta, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	Forward2D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	const double detJacob = calcJacobianMatrix( pMeshData, iElem, xi, eta, jacobMat );
	Forward2D::Matrix2x2 invJacobMat;
	calcInverseOfJacobianMatrix( jacobMat, detJacob, invJacobMat );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 4; ++i ){
		const double length = pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncRotated( xi, eta, i, invJacobMat ) * length, 0.0 );
	}

	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
	const double factor = omega * CommonParameters::mu;

	val /= std::complex<double>(0.0, factor);

	return val;

}

// Calculate jacobian matrix of the elements on the Z-X plane of boundary
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcJacobianMatrixOnZXPlaneOfBoundary( const MeshDataNonConformingHexaElement* const pMeshData, 
	const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const{

	double xCoord[4] = { 0.0, 0.0, 0.0, 0.0 };
	double zCoord[4] = { 0.0, 0.0, 0.0, 0.0 };
	for( int i = 0; i < 4; ++i ){
		xCoord[i] = pMeshData->getCoordXFromElementBoundaryPlanes( m_planeID, elemID2D, i );
		zCoord[i] = pMeshData->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, i );
	}

	// Zero clear
	jacobMat.mat11 = 0.0;
	jacobMat.mat12 = 0.0;
	jacobMat.mat21 = 0.0;
	jacobMat.mat22 = 0.0;
	for( int i = 0; i < 4; ++i ){
		const double xiNode   = m_xiAtNode[i];
		const double etaNode  = m_etaAtNode[i];
		const double tmp1 = 0.25 * xiNode   * (1.0 + etaNode  * eta);
		const double tmp2 = 0.25 * etaNode  * (1.0 + xiNode * xi);
		jacobMat.mat11 += tmp1 * xCoord[i];
		jacobMat.mat12 += tmp1 * zCoord[i];
		jacobMat.mat21 += tmp2 * xCoord[i];
		jacobMat.mat22 += tmp2 * zCoord[i];
	}

	return jacobMat.mat11 * jacobMat.mat22 - jacobMat.mat12 * jacobMat.mat21;

}

// Calculate jacobian matrix of the elements on the Y-Z plane of boundary
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcJacobianMatrixOnYZPlaneOfBoundary( const MeshDataNonConformingHexaElement* const pMeshData,
	const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const{

	double yCoord[4] = { 0.0, 0.0, 0.0, 0.0 };
	double zCoord[4] = { 0.0, 0.0, 0.0, 0.0 };
	for( int i = 0; i < 4; ++i ){
		yCoord[i] = pMeshData->getCoordYFromElementBoundaryPlanes( m_planeID, elemID2D, i );
		zCoord[i] = pMeshData->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, i );
	}

	// Zero clear
	jacobMat.mat11 = 0.0;
	jacobMat.mat12 = 0.0;
	jacobMat.mat21 = 0.0;
	jacobMat.mat22 = 0.0;
	for( int i = 0; i < 4; ++i ){
		const double xiNode   = m_xiAtNode[i];
		const double etaNode  = m_etaAtNode[i];
		const double tmp1 = 0.25 * xiNode   * (1.0 + etaNode  * eta);
		const double tmp2 = 0.25 * etaNode  * (1.0 + xiNode * xi);
		jacobMat.mat11 += tmp1 * yCoord[i];
		jacobMat.mat12 += tmp1 * zCoord[i];
		jacobMat.mat21 += tmp2 * yCoord[i];
		jacobMat.mat22 += tmp2 * zCoord[i];
	}

	return jacobMat.mat11 * jacobMat.mat22 - jacobMat.mat12 * jacobMat.mat21;

}

// Calculate jacobian matrix
double Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcJacobianMatrix( const MeshDataNonConformingHexaElement* const pMeshData,
	const int elemID2D, const double xi, const double eta, Forward2D::Matrix2x2& jacobMat ) const{

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		return calcJacobianMatrixOnYZPlaneOfBoundary( pMeshData, elemID2D, xi, eta, jacobMat );
	}else{//ZX Plane
		return calcJacobianMatrixOnZXPlaneOfBoundary( pMeshData, elemID2D, xi, eta, jacobMat );
	}

}

// Calculate inverse of jacobian matrix  multiplied by determinant
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::calcInverseOfJacobianMatrix( const Forward2D::Matrix2x2& jacobMat, const double determinant, Forward2D::Matrix2x2& invJacobMat ) const{

	const double invDet = 1.0 / determinant;

	invJacobMat.mat11 =   jacobMat.mat22 * invDet; 
	invJacobMat.mat12 = - jacobMat.mat12 * invDet; 
	invJacobMat.mat21 = - jacobMat.mat21 * invDet; 
	invJacobMat.mat22 =   jacobMat.mat11 * invDet; 

}

// Add master dof and factor pair to m_slaveDofToMasterDofAndFactors
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::addMasterDofAndFactorPair( const int slaveDof, const int masterDof, const double factor  ){

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

// Set non-zero strucuture of matrix for forward calculation
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::setNonZeroStrucuture( const MeshDataNonConformingHexaElement* const pMeshData ){

	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();

	// If matrix structure has not been set yet, set matrix structure.
	OutputFiles::m_logFile << "# Set matrix structure. " << pAnalysisControl->outputElapsedTime() << std::endl;

	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){
			const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){
				const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
				if( col < 0 ){
					continue;
				}
				//if( col >= row ){// Store only upper triangle part
				//	m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
				//}
				std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[row];
				std::vector< std::pair<int,double> >& colMasters= m_slaveDofToMasterDofAndFactors[col];
				for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					for( std::vector< std::pair<int,double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol ){
						const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
						if( colMod >= rowMod ){// Store only upper triangle part
							m_matrix2DAnalysis.setStructureByTripletFormat( rowMod, colMod );
						}
					}
				}
			}// iEdge2
		}// iEdge1
	}// iElem

}

// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::setNonZeroValues( const double freq, const MeshDataNonConformingHexaElement* const pMeshData ){

	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	const int nElem = pMeshData->getNumElemOnBoundaryPlanes( m_planeID );
	for( int iElem = 0; iElem < nElem; ++iElem ){
		//--- Calculate omega * mu * sigma
		const int elemID = pMeshData->getElemBoundaryPlanes( m_planeID, iElem );
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		const std::complex<double> factor = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form
		const double length[4] = { pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 0 ),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 1 ),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 2 ),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, 3 ) };
		for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){
			const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
			if( row < 0 ){
				continue;
			}
			for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){
				const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
				if( col <= Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}
				double integral1(0.0);
				double integral2(0.0);
				for( int ip = 0; ip < m_numIntegralPoints; ++ip ){
					const double xi = m_integralPointXi[ip];
					const double eta = m_integralPointEta[ip];
					Forward2D::Matrix2x2 jacobMat = { 0.0, 0.0, 0.0, 0.0 };
					const double detJacob = calcJacobianMatrix( pMeshData, iElem, xi, eta, jacobMat );
					Forward2D::Matrix2x2 invJacobMat;
					calcInverseOfJacobianMatrix( jacobMat, detJacob, invJacobMat );
					integral1 += getShapeFuncRotated(xi, eta, iEdge1, invJacobMat) * getShapeFuncRotated(xi, eta, iEdge2, invJacobMat) * detJacob * m_weights[ip];
					integral2 += ( getShapeFuncH(xi, eta, iEdge1, invJacobMat) * getShapeFuncH(xi, eta, iEdge2, invJacobMat)
						         + getShapeFuncV(xi, eta, iEdge1, invJacobMat) * getShapeFuncV(xi, eta, iEdge2, invJacobMat) ) * detJacob * m_weights[ip];
				}
				integral1 *= length[iEdge1] * length[iEdge2];
				integral2 *= length[iEdge1] * length[iEdge2];
				std::complex<double> val = std::complex<double>( integral1 , 0.0 ) - std::complex<double>( integral2 , 0.0 ) * factor;// exp(-i*omega*t) form
				//if( col == Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
				//	const std::complex<double> nonZeroValue = m_edgesIDGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge2] ];
				//	m_matrix2DAnalysis.addRightHandSideVector( row, -val * nonZeroValue );// Add to right hand side vector
				//}else if( col >= row ){// Store only upper triangle part
				//	m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				//}
				const std::vector< std::pair<int,double> >& rowMasters= m_slaveDofToMasterDofAndFactors[row];
				for( std::vector< std::pair<int,double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow ){
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
					if( col == Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
						const std::complex<double> nonZeroValue = m_edgesIDGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge2] ];
						m_matrix2DAnalysis.addRightHandSideVector( rowMod, -valMod * nonZeroValue );// Add to right hand side vector
					}else{
						const std::vector< std::pair<int,double> >& colMasters= m_slaveDofToMasterDofAndFactors[col];
						for( std::vector< std::pair<int,double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol ){
							const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
							const std::complex<double> valModMod = valMod * std::complex<double>(itrCol->second, 0.0);
							if( colMod >= rowMod ){// Store only upper triangle part
								m_matrix2DAnalysis.addNonZeroValues( rowMod, colMod, valModMod );// Add to matrix
							}
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}// iElem

}

// Get flag specifing whether an 2-D element faces slave element
bool Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::faceSlaveElements( const int iElem, const int iEdge, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	const int elemID = pMeshData->getElemBoundaryPlanes( m_planeID, iElem );
	const int faceID = getNeighborFaceIndexFromEdgeIndex(iEdge);
	if( pMeshData->faceSlaveElements(elemID, faceID) ){
		return true;
	}

	return false;

}

// Get flag specifing whether the inputted element edge is an outer edge
bool Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::isOuterEdge( const int iElem, const int iEdge, const MeshDataNonConformingHexaElement* const pMeshData ) const{

	const int elemID = pMeshData->getElemBoundaryPlanes( m_planeID, iElem );
	const int faceID = getNeighborFaceIndexFromEdgeIndex(iEdge);
	if( pMeshData->isOuterBoundary(elemID, faceID) ){
		return true;
	}

	return false;

}

// Get neighbor face index from edge index
int Forward2DNonConformingQuadrilateralElement0thOrderEdgeBased::getNeighborFaceIndexFromEdgeIndex( const int iEdge ) const{

	switch (m_planeID){
		case MeshData::YZMinus: // Go through
		case MeshData::YZPlus:
			switch (iEdge){
				case 0:
					return 4;
					break;
				case 1:
					return 5;
					break;
				case 2:
					return 2;
					break;
				case 3:
					return 3;
					break;
				default:
					OutputFiles::m_logFile << "Error : Edge ID is wrong : iEdge = " << iEdge << "." << std::endl;
					exit(1);					
					break;
			}
			break;
		case MeshData::ZXMinus: // Go through
		case MeshData::ZXPlus:
			switch (iEdge){
				case 0:
					return 4;
					break;
				case 1:
					return 5;
					break;
				case 2:
					return 0;
					break;
				case 3:
					return 1;
					break;
				default:
					OutputFiles::m_logFile << "Error : Edge ID is wrong : iEdge = " << iEdge << "." << std::endl;
					exit(1);					
					break;
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Plane ID is wrong : " << m_planeID << "." << std::endl;
			exit(1);		
			break;
	}

	return -1;

}