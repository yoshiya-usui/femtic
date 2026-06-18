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
#include "Forward2DTriangleElement0thOrderEdgeBased.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "ResistivityBlock.h"
#include <algorithm>
#include "MeshDataTetraElement.h"
#include <assert.h>

// Constructer
Forward2DTriangleElement0thOrderEdgeBased::Forward2DTriangleElement0thOrderEdgeBased( const int planeID, const int iPol ):
	Forward2DTriangleElementEdgeBased( planeID, iPol )
{
	// U coordinates of integral points
	m_uCoord[0] = 0.5;
	m_uCoord[1] = 0.5;
	m_uCoord[2] = 0.0;

	// V coordinates of integral points
	m_vCoord[0] = 0.0;
	m_vCoord[1] = 0.5;
	m_vCoord[2] = 0.5;

	// Weights of integral points
	m_weight[0] = 0.3333333333333333333;
	m_weight[1] = 0.3333333333333333333;
	m_weight[2] = 0.3333333333333333333;
}

// Destructer
Forward2DTriangleElement0thOrderEdgeBased::~Forward2DTriangleElement0thOrderEdgeBased(){
}

// Calculate EM fields of boundary planes by 2D forward calculcation with 0th order edge element
void Forward2DTriangleElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataTetraElement* const pMeshDataTetraElement ){

	const int imode = calcMode();// TM or TE mode
	if ( imode != CommonParameters::TM_MODE ){
		OutputFiles::m_logFile << "Error : Only TM mode can be treated in Forward2DTriangleElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes ! imode = " << imode << "." << std::endl;
		exit(1);
	}
#ifdef _DEBUG_WRITE
	std::cout << "imode " << imode << std::endl;// For debug
#endif

	if( m_sourceFieldElectric == false ){
		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in Forward2DTriangleElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes !" << std::endl;
		exit(1);
	}
	
	const int nElem = pMeshDataTetraElement->getNumElemOnBoundaryPlanes( m_planeID );

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
	if( m_hasAlreadySetIDsLocal2Global == false ){
		calcArrayConvertLocalID2Global( pMeshDataTetraElement );
	}

	OutputFiles::m_logFile << "# Number of equation = " << m_numEquations
		<< ", Number of equation after degeneration = " << m_numEquationsDegenerated << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() <= 0 ){
		m_matrix2DAnalysis.setDegreeOfEquation( m_numEquationsDegenerated );
	}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	std::cout << "m_numEquations = " << m_numEquations << std::endl;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem m_IDsLocal2Global : " << iElem;
		for( int i = 0; i < 3; ++i ){
			std::cout << " " << m_IDsLocal2Global[iElem][i];
		}
		std::cout << std::endl;
	}

	std::cout << "m_numEquationsDegenerated = " << m_numEquationsDegenerated << std::endl;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
		for( int i = 0; i < 3; ++i ){
			std::cout << " " << m_IDsLocal2GlobalDegenerated[iElem][i];
		}
		std::cout << std::endl;
	}

	for( int iElem = 0; iElem < nElem; ++iElem ){
		std::cout << "iElem m_signInversion : " << iElem;
		for( int i = 0; i < 3; ++i ){
			std::cout << " " << m_signInversion[iElem][i];
		}
		std::cout << std::endl;
	}
#endif
	//----- debug <<<<<

	//-----------------------------------------------------------------------
	//--- Set matrix structure and analyze the structure by matrix solver ---
	//-----------------------------------------------------------------------
	if( m_hasMatrixStructureSetAndAnalyzed == false ){
		// If matrix structure has not been set yet, set matrix structure.

		OutputFiles::m_logFile << "# Set matrix structure. " << pAnalysisControl->outputElapsedTime() << std::endl;

		//------------------------------------------
		//--- Components due to stiffness matrix ---
		//------------------------------------------
		for( int iElem = 0; iElem < nElem; ++iElem ){

			for( int iEdge1 = 0; iEdge1 < 3; ++iEdge1 ){

				const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
				if( row < 0 ){
					continue;
				}

				for( int iEdge2 = 0; iEdge2 < 3; ++iEdge2 ){
					const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
					if( col < 0 ){
						continue;
					}

					if( col >= row ){// Store only upper triangle part
						m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
					}

				}// iEdge2

			}// iEdge1

		}// iElem

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

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	m_matrix2DAnalysis.zeroClearNonZeroValues();// Zero clear matrix values
	m_matrix2DAnalysis.zeroClearRightHandSideVector();// Zero clear right hand side vector

	const double sourceValueElectric = CommonParameters::sourceValueElectric;
	setNonZeroValues(freq, pMeshDataTetraElement);

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	m_matrix2DAnalysis.debugWriteMatrix();
	m_matrix2DAnalysis.debugWriteRightHandSide();
#endif
	//----- debug <<<<<

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
	std::complex<double>* solutionDegenerated = new std::complex<double>[ m_numEquationsDegenerated ];
	OutputFiles::m_logFile << "# Solve phase of matrix solver for boundary plane."
			<<  " Polarization : " << m_polarization << " Plane : " << m_planeID <<  ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix2DAnalysis.solvePhaseMatrixSolver( solutionDegenerated );//Solve phase of matrix solver

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0; i < m_numEquationsDegenerated; ++i ){
		std::cout << "i solutionDegenerated : " << i << " " << solutionDegenerated[i] << std::endl;
	}
#endif
	//----- debug <<<<<

	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[ m_numEquations ];

	bool* alreadyFound = new bool[ m_numEquations ];
	for( int i = 0; i < m_numEquations; ++i ){
		alreadyFound[i] = false;
	}

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 3; ++iEdge ){

			const int iNum = m_IDsLocal2Global[iElem][iEdge];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsLocal2GlobalDegenerated[iElem][iEdge];

				if( iNumDegenerated == DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					//m_solution[iNum] = std::complex<double>( sourceValueElectric, 0.0 );
					m_solution[iNum] = m_edgesIDGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge] ];
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
	for( int i = 0; i < m_numEquations; ++i ){
		std::cout << "i m_solution : " << i << " " << m_solution[i] << std::endl;
	}
#endif
	//----- debug <<<<<

	output2DResult( freq, pMeshDataTetraElement );

}

// Calculate array converting global element IDs of 2D mesh to global edge IDs of 2D mesh
void Forward2DTriangleElement0thOrderEdgeBased::calcArrayConvertLocalID2Global( const MeshDataTetraElement* const pMeshDataTetraElement ){

	//typedef std::pair< int, std::pair< int, int> > inputedDataType;
	typedef std::pair< int, int > NodeIDPair;
	typedef std::pair< int, int > ElemAndEdgeID;
	typedef std::pair< int, std::pair< NodeIDPair, ElemAndEdgeID > > InputedDataType;

	std::vector< InputedDataType > indicatorAndLocalEdgeID;

	const int nElem = pMeshDataTetraElement->getNumElemOnBoundaryPlanes( m_planeID );

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 3; ++iEdge ){
			const int nodeID0 = pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, iEdge, 0 );
			const int nodeID1 = pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, iEdge, 1 );
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
		for( int i = 0; i < m_sizeOfIDsLocal2Global; ++i ){
			delete [] m_IDsLocal2Global[i];
		}
		delete [] m_IDsLocal2Global;
	}
	m_IDsLocal2Global = new int*[nElem]; 
	m_sizeOfIDsLocal2Global = nElem;
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_IDsLocal2Global[iElem] = new int[3];
	}
	//-----------------------------------------

	// Allocate memory to m_IDsLocal2GlobalDegenerated ---
	if( m_IDsLocal2GlobalDegenerated != NULL ){
		for (int i = 0; i < m_sizeOfIDsLocal2Global; ++i) {
			delete [] m_IDsLocal2GlobalDegenerated[i];
		}
		delete [] m_IDsLocal2GlobalDegenerated;
	}
	m_IDsLocal2GlobalDegenerated = new int*[nElem]; 
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_IDsLocal2GlobalDegenerated[iElem] = new int[3];
	}
	//----------------------------------------------------

	// Allocate memory to m_signInversion ---
	if( m_signInversion != NULL ){
		for( int iElem = 0; iElem < nElem; ++iElem ){
			delete [] m_signInversion[iElem];
		}
		delete [] m_signInversion;
	}
	m_signInversion = new bool*[nElem]; 
	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_signInversion[iElem] = new bool[3];
	}
	//----------------------------------------

	int indicatorPre(-1);
	std::pair<int, int> nodePairPre;
	int edgeIDGlobal2D(-1);
	int edgeIDGlobal2DDegenerated(-1);
	bool outerEdge(false);

	std::vector<InputedDataType>::iterator itrEnd = indicatorAndLocalEdgeID.end();
	for( std::vector<InputedDataType>::iterator itr = indicatorAndLocalEdgeID.begin(); itr != itrEnd; ++itr ){

		const int elemIDGlobal2D = (*itr).second.second.first;
		const int edgeIDLocal2D = (*itr).second.second.second;

		const int nodePairCur0 = pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 0 );
		const int nodePairCur1 = pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 1 );

		//const int nodePairCur0 = (*itr).second.first.first;
		//const int nodePairCur1 = (*itr).second.first.second;

//#ifdef _DEBUG_WRITE
//		std::cout << "indicator nodePairCur0 nodePairCur1 nodePairPre.first nodePairPre.second :  " << (*itr).first << " "  << nodePairCur0 << " " << nodePairCur1 << " " << nodePairPre.first << " " << nodePairPre.second << std::endl;
//#endif

		if( (*itr).first == indicatorPre &&
			( ( nodePairCur0 == nodePairPre.first && nodePairCur1 == nodePairPre.second ) ||
			  ( nodePairCur0 == nodePairPre.second && nodePairCur1 == nodePairPre.first ) ) ){// Same edge with privious one

//#ifdef _DEBUG_WRITE
//			std::cout << "OK" << std::endl;
//#endif
			outerEdge = false;

		}else{// New edge

			if( outerEdge ){// If the previous edge is outer edge
				
				const int elemIDGlobal2DPre = (*(itr-1)).second.second.first;
				const int edgeIDLocal2DPre = (*(itr-1)).second.second.second;

				if( getTypeOfOuterEdgeOfBoundaryPlanes( elemIDGlobal2DPre, edgeIDLocal2DPre, pMeshDataTetraElement ) == Forward2DTriangleElementEdgeBased::UPPER_EDGE ){
					m_IDsLocal2GlobalDegenerated[elemIDGlobal2DPre][edgeIDLocal2DPre] = Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE;
				}else{
					m_IDsLocal2GlobalDegenerated[elemIDGlobal2DPre][edgeIDLocal2DPre] = Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE;
				}

				--edgeIDGlobal2DDegenerated;
			}
			if( itr+1 == itrEnd ){// Last element 

				if( getTypeOfOuterEdgeOfBoundaryPlanes( elemIDGlobal2D, edgeIDLocal2D, pMeshDataTetraElement ) == Forward2DTriangleElementEdgeBased::UPPER_EDGE ){
					m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE;
				}else{
					m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_ZERO_VALUE;
				}

				m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] = ++edgeIDGlobal2D;
				m_signInversion[elemIDGlobal2D][edgeIDLocal2D] = nodePairCur1 > nodePairCur0 ? false : true; 

#ifdef _DEBUG_WRITE
				std::cout << "elemIDGlobal2D edgeIDLocal2D nodePairCur0 nodePairCur1 m_IDsLocal2Global m_signInversion :  " << elemIDGlobal2D << " " << edgeIDLocal2D << " " << nodePairCur0 << " " << nodePairCur1 << " " << m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] << " " << m_signInversion[elemIDGlobal2D][edgeIDLocal2D] << std::endl;
#endif

				break;
			}

			++edgeIDGlobal2D;
			++edgeIDGlobal2DDegenerated;

			indicatorPre = (*itr).first;
			nodePairPre.first  = nodePairCur0;
			nodePairPre.second = nodePairCur1;

			outerEdge = true;
		}

		m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] = edgeIDGlobal2D;
		m_IDsLocal2GlobalDegenerated[elemIDGlobal2D][edgeIDLocal2D] = edgeIDGlobal2DDegenerated;


		m_signInversion[elemIDGlobal2D][edgeIDLocal2D] = nodePairCur1 > nodePairCur0 ? false : true; 

#ifdef _DEBUG_WRITE
		std::cout << "elemIDGlobal2D edgeIDLocal2D nodePairCur0 nodePairCur1 m_IDsLocal2Global m_signInversion :  " << elemIDGlobal2D << " " << edgeIDLocal2D << " " << nodePairCur0 << " " << nodePairCur1 << " " << m_IDsLocal2Global[elemIDGlobal2D][edgeIDLocal2D] << " " << m_signInversion[elemIDGlobal2D][edgeIDLocal2D] << std::endl;
#endif

	}

	const double sourceValueElectric = CommonParameters::sourceValueElectric;

	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 3; ++iEdge ){

			if( m_IDsLocal2GlobalDegenerated[iElem][iEdge] == Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
				double val = pMeshDataTetraElement->calcHorizontalCoordDifferenceBoundaryPlanes( m_planeID, iElem, iEdge ) > 0 ? sourceValueElectric : -sourceValueElectric;
				if( m_signInversion[iElem][iEdge] ){
					val *= -1.0;
				}
				m_edgesIDGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][iEdge], std::complex<double>(val, 0.0) ) );
			}

			//if( m_IDsLocal2GlobalDegenerated[iElem][iEdge] == Forward2DTriangleElement0thOrderEdgeBased::DIRICHLET_BOUNDARY_NONZERO_VALUE ){
			//	double val = pMeshDataTetraElement->calcHorizontalCoordDifferenceBoundaryPlanes( m_planeID, iElem, iEdge ) * sourceValueElectric;
			//	if( m_signInversion[iElem][iEdge] ){
			//		val *= -1.0;
			//	}
			//	m_edgesIDGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][iEdge], std::complex<double>(val, 0.0) ) );
			//}

		}
	}

#ifdef _DEBUG_WRITE
	for( std::map< int, std::complex<double> >::iterator itr = m_edgesIDGlobal2NonZeroValues.begin(); itr != m_edgesIDGlobal2NonZeroValues.end(); ++itr ){
		std::cout << "m_edgesIDGlobal2NonZeroValues first second :  " << itr->first << " " << itr->second << std::endl;
	}
#endif

	m_numEquations = edgeIDGlobal2D + 1;
	m_numEquationsDegenerated = edgeIDGlobal2DDegenerated + 1;
	
	m_hasAlreadySetIDsLocal2Global = true;

}


// Get type of outer edge
// [note] : You must confirm inputed edge is the outer edge. 
int Forward2DTriangleElement0thOrderEdgeBased::getTypeOfOuterEdgeOfBoundaryPlanes( const int elemIDGlobal2D,
	const int edgeIDLocal2D, const MeshDataTetraElement* const pMeshDataTetraElement ) const{

	const int elemIDGlobal3D = pMeshDataTetraElement->getElemBoundaryPlanes( m_planeID, elemIDGlobal2D );

	const int nodeIDGlobal3D[2] = { pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 0 ),
									pMeshDataTetraElement->getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( m_planeID, elemIDGlobal2D, edgeIDLocal2D, 1 ) }; 

	const double eps = 1.0e-6;

	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//Y-Z Plane

		if( fabs( pMeshDataTetraElement->getYCoordinatesOfNodes( nodeIDGlobal3D[0] ) - pMeshDataTetraElement->getYCoordinatesOfNodes( nodeIDGlobal3D[1] ) ) < eps ){

			const double yCoordOfCenter = ( pMeshDataTetraElement->getCenterCoord( elemIDGlobal3D ) ).Y;

			if( pMeshDataTetraElement->getYCoordinatesOfNodes( nodeIDGlobal3D[0] ) < yCoordOfCenter ){

				return Forward2DTriangleElementEdgeBased::LEFT_EDGE;

			}else if( pMeshDataTetraElement->getYCoordinatesOfNodes( nodeIDGlobal3D[0] ) > yCoordOfCenter ){ 

				return Forward2DTriangleElementEdgeBased::RIGHT_EDGE;

			}else{
				OutputFiles::m_logFile << "Error : Center locate on the inputed edge !!" << std::endl;
				exit(1);
			}

		}

	}else if( m_planeID == MeshData::ZXMinus || m_planeID == MeshData::ZXPlus ){//Z-X Plane

		if( fabs( pMeshDataTetraElement->getXCoordinatesOfNodes( nodeIDGlobal3D[0] ) - pMeshDataTetraElement->getXCoordinatesOfNodes( nodeIDGlobal3D[1] ) ) < eps ){

			const double xCoordOfCenter = ( pMeshDataTetraElement->getCenterCoord( elemIDGlobal3D ) ).X;

			if( pMeshDataTetraElement->getXCoordinatesOfNodes( nodeIDGlobal3D[0] ) < xCoordOfCenter ){

				return Forward2DTriangleElementEdgeBased::LEFT_EDGE;

			}else if( pMeshDataTetraElement->getXCoordinatesOfNodes( nodeIDGlobal3D[0] ) > xCoordOfCenter ){ 

				return Forward2DTriangleElementEdgeBased::RIGHT_EDGE;

			}else{
				OutputFiles::m_logFile << "Error : Center locate on the inputed edge !!" << std::endl;
				exit(1);
			}

		}

	}else{
		OutputFiles::m_logFile << "Error : Wrong plane ID !! m_planeID = " << m_planeID << std::endl;
		exit(1);
	}
	
	if ( fabs( pMeshDataTetraElement->getZCoordinatesOfNodes( nodeIDGlobal3D[0] ) - pMeshDataTetraElement->getZCoordinatesOfNodes( nodeIDGlobal3D[1] ) ) < eps ){

		const double zCoordOfCenter = ( pMeshDataTetraElement->getCenterCoord( elemIDGlobal3D ) ).Z;

		if( pMeshDataTetraElement->getZCoordinatesOfNodes( nodeIDGlobal3D[0] ) < zCoordOfCenter ){

			return Forward2DTriangleElementEdgeBased::UPPER_EDGE;

		}else if( pMeshDataTetraElement->getZCoordinatesOfNodes( nodeIDGlobal3D[0] ) > zCoordOfCenter ){ 

			return Forward2DTriangleElementEdgeBased::LOWER_EDGE;

		}else{
			OutputFiles::m_logFile << "Error : Center locate on the inputed edge !!" << std::endl;
			exit(1);
		}

	}else{
		OutputFiles::m_logFile << "Error : Inputed edge is parallel to neither horizontal direction nor vertical direction !! elemIDGlobal2D = " << elemIDGlobal2D << ", edgeIDLocal2D= " << edgeIDLocal2D << std::endl;
		exit(1);
	}

	return Forward2DTriangleElementEdgeBased::INNER_EDGE;

}

// Get shape functions of the 1st direction with respect to the reference element coordinate system
inline double Forward2DTriangleElement0thOrderEdgeBased::getShapeFuncReferenceCoordU( const double uLocal, const double vLocal, const int num ) const{

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
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncReferenceCoordU : num = " << num << std::endl;
			exit(1);
			break;
	}


}

// Get shape functions of the 2nd direction with respect to the reference element coordinate system
inline double Forward2DTriangleElement0thOrderEdgeBased::getShapeFuncReferenceCoordV( const double uLocal, const double vLocal, const int num ) const{

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
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncReferenceCoordV : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions rotated with respect to the reference element coordinate system
inline double Forward2DTriangleElement0thOrderEdgeBased::getShapeFuncRotated() const{

	return 2.0;

}

// Calculate horizontal electric field
std::complex<double> Forward2DTriangleElement0thOrderEdgeBased::calcValueElectricFieldHorizontal( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
	//	exit(1);
	//}else if( m_IDsLocal2Global[iElem] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);

	for( int i = 0; i < 3; ++i ){
		const double length = pMeshDataTetraElement->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		if( m_signInversion[iElem][i] ){
			valU -= m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV -= m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}else{
			valU += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}
	}

	Forward2DTriangleElementEdgeBased::JacobianMatrix jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		calcJacobianMatrixOnYZPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}else{//ZX Plane
		calcJacobianMatrixOnZXPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}

	return valU * std::complex<double>(jacobMat.jacob22/detJacob, 0.0) - valV * std::complex<double>(jacobMat.jacob12/detJacob, 0.0);

}

// Calculate vertical electric field
std::complex<double> Forward2DTriangleElement0thOrderEdgeBased::calcValueElectricFieldVertical( const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
	//	exit(1);
	//}else if( m_IDsLocal2Global[iElem] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	std::complex<double> valU(0.0, 0.0);
	std::complex<double> valV(0.0, 0.0);

	for( int i = 0; i < 3; ++i ){
		const double length = pMeshDataTetraElement->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		if( m_signInversion[iElem][i] ){
			valU -= m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV -= m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}else{
			valU += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordU( uCoord, vCoord, i ) * length, 0.0 );
			valV += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( getShapeFuncReferenceCoordV( uCoord, vCoord, i ) * length, 0.0 );
		}
	}
	
	Forward2DTriangleElementEdgeBased::JacobianMatrix jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		calcJacobianMatrixOnYZPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}else{//ZX Plane
		calcJacobianMatrixOnZXPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}

	return valV * std::complex<double>(jacobMat.jacob11/detJacob, 0.0) - valU * std::complex<double>(jacobMat.jacob21/detJacob, 0.0);

}

// Calculate magnetic field perpendicular to the boundary plane
std::complex<double> Forward2DTriangleElement0thOrderEdgeBased::calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double uCoord, const double vCoord, const MeshDataTetraElement* const pMeshDataTetraElement ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
	//	exit(1);
	//}else if( m_IDsLocal2Global[iElem] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 3; ++i ){
		const double length = pMeshDataTetraElement->calcEdgeLengthFromElementAndEdgeBoundaryPlanes( m_planeID, iElem, i );
		if( m_signInversion[iElem][i] ){
			val -= m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( length * getShapeFuncRotated(), 0.0 );
		}else{
			val += m_solution[ m_IDsLocal2Global[iElem][i] ] * std::complex<double>( length * getShapeFuncRotated(), 0.0 );
		}
	}
	
	Forward2DTriangleElementEdgeBased::JacobianMatrix jacobMat = { 0.0, 0.0, 0.0, 0.0 };
	double detJacob(0.0);
	if( m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		calcJacobianMatrixOnYZPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}else{//ZX Plane
		calcJacobianMatrixOnZXPlaneOfBoundary( pMeshDataTetraElement, iElem, jacobMat, detJacob );
	}

	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
	const double factor = omega * CommonParameters::mu;

	val /= std::complex<double>( 0.0, detJacob * factor );

	return val;

}


// Calculate jacobian matrix of the elements on the Z-X plane of boundary
void Forward2DTriangleElement0thOrderEdgeBased::calcJacobianMatrixOnZXPlaneOfBoundary( const MeshDataTetraElement* const pMeshDataTetraElement, const int elemID2D,
	Forward2DTriangleElementEdgeBased::JacobianMatrix& JacobMat, double& determinant ) const{

	JacobMat.jacob11 = pMeshDataTetraElement->getCoordXFromElementBoundaryPlanes( m_planeID, elemID2D, 1 ) - pMeshDataTetraElement->getCoordXFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob12 = pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 1 ) - pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob21 = pMeshDataTetraElement->getCoordXFromElementBoundaryPlanes( m_planeID, elemID2D, 2 ) - pMeshDataTetraElement->getCoordXFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob22 = pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 2 ) - pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );

	determinant = JacobMat.jacob11 * JacobMat.jacob22 - JacobMat.jacob12 * JacobMat.jacob21;
	
}

// Calculate jacobian matrix of the elements on the Y-Z plane of boundary
void Forward2DTriangleElement0thOrderEdgeBased::calcJacobianMatrixOnYZPlaneOfBoundary( const MeshDataTetraElement* const pMeshDataTetraElement, const int elemID2D,
	Forward2DTriangleElementEdgeBased::JacobianMatrix& JacobMat, double& determinant ) const{

	JacobMat.jacob11 = pMeshDataTetraElement->getCoordYFromElementBoundaryPlanes( m_planeID, elemID2D, 1 ) - pMeshDataTetraElement->getCoordYFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob12 = pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 1 ) - pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob21 = pMeshDataTetraElement->getCoordYFromElementBoundaryPlanes( m_planeID, elemID2D, 2 ) - pMeshDataTetraElement->getCoordYFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );
	JacobMat.jacob22 = pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 2 ) - pMeshDataTetraElement->getCoordZFromElementBoundaryPlanes( m_planeID, elemID2D, 0 );

	determinant = JacobMat.jacob11 * JacobMat.jacob22 - JacobMat.jacob12 * JacobMat.jacob21;

}

// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward2DTriangleElement0thOrderEdgeBased::setNonZeroValues(const double freq, const MeshDataTetraElement* const pMeshData){

	if ((AnalysisControl::getInstance())->isAnisotropyConsidered()) {
		setNonZeroValuesForAnisotropicConductivity(freq, pMeshData);
	}
	else {
		setNonZeroValuesForIsotropicConductivity(freq, pMeshData);
	}

}

// Set non-zero values of matrix and right-hande side vector for isotropic conductity
void Forward2DTriangleElement0thOrderEdgeBased::setNonZeroValuesForIsotropicConductivity(const double freq, const MeshDataTetraElement* const pMeshData) {

	const ResistivityBlockIsotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockIsotropic();
	const int nElem = pMeshData->getNumElemOnBoundaryPlanes(m_planeID);
	for (int iElem = 0; iElem < nElem; ++iElem) {
		//--- Calculate Jacobian
		Forward2DTriangleElementEdgeBased::JacobianMatrix jacobMat = { 0.0, 0.0, 0.0, 0.0 };
		double detJacob(0.0);
		if (m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus) {//YZ Plane
			calcJacobianMatrixOnYZPlaneOfBoundary(pMeshData, iElem, jacobMat, detJacob);
		}
		else {//ZX Plane
			calcJacobianMatrixOnZXPlaneOfBoundary(pMeshData, iElem, jacobMat, detJacob);
		}
		const double divDetJacob = 1.0 / detJacob;

		//--- Calculate omega * mu * sigma
		const int elemID = pMeshData->getElemBoundaryPlanes(m_planeID, iElem);
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		const double length[3] = { pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 0),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 1),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 2) };
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		const std::complex<double> factor = std::complex<double>(0.0, omega * CommonParameters::mu * sigma);// exp(-i*omega*t) form

		for (int iEdge1 = 0; iEdge1 < 3; ++iEdge1) {
			const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 3; ++iEdge2) {
				const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
				if (col <= DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1(0.0);
				double integral2(0.0);
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					integral1 += getShapeFuncRotated() * getShapeFuncRotated() * m_weight[ip];
					integral2 += ((getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob22
								- getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob12)
								* (getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob22
								- getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob12)
								+ (getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob11
								- getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob21)
								* (getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob11
								- getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob21)) * m_weight[ip];
				}
				if (m_signInversion[iElem][iEdge1] != m_signInversion[iElem][iEdge2]) {
					integral1 *= -divDetJacob * length[iEdge1] * length[iEdge2];
					integral2 *= -divDetJacob * length[iEdge1] * length[iEdge2];
				}
				else {
					integral1 *= divDetJacob * length[iEdge1] * length[iEdge2];
					integral2 *= divDetJacob * length[iEdge1] * length[iEdge2];
				}
				std::complex<double> val = std::complex<double>(integral1, 0.0) - std::complex<double>(integral2, 0.0) * factor;// exp(-i*omega*t) form
				if (col == DIRICHLET_BOUNDARY_NONZERO_VALUE) {
					const std::complex<double> nonZeroValue = m_edgesIDGlobal2NonZeroValues[m_IDsLocal2Global[iElem][iEdge2]];
					m_matrix2DAnalysis.addRightHandSideVector(row, -val * nonZeroValue);// Add to right hand side vector
				}
				else if (col >= row) {// Store only upper triangle part
					m_matrix2DAnalysis.addNonZeroValues(row, col, val);// Add to matrix
				}
			}// iEdge2
		}// iEdge1		
	}// iElem
}

// Set non-zero values of matrix and right-hande side vector for anisotropic conductity
void Forward2DTriangleElement0thOrderEdgeBased::setNonZeroValuesForAnisotropicConductivity(const double freq, const MeshDataTetraElement* const pMeshData) {

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	const int nElem = pMeshData->getNumElemOnBoundaryPlanes(m_planeID);
	for (int iElem = 0; iElem < nElem; ++iElem) {
		Forward2DTriangleElementEdgeBased::JacobianMatrix jacobMat = { 0.0, 0.0, 0.0, 0.0 };
		double detJacob(0.0);
		if (m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus) {//YZ Plane
			calcJacobianMatrixOnYZPlaneOfBoundary(pMeshData, iElem, jacobMat, detJacob);
		}
		else {//ZX Plane
			calcJacobianMatrixOnZXPlaneOfBoundary(pMeshData, iElem, jacobMat, detJacob);
		}
		const double divDetJacob = 1.0 / detJacob;
		//--- Calculate omega * mu * sigma
		const int elemID = pMeshData->getElemBoundaryPlanes(m_planeID, iElem);
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		const double length[3] = { pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 0),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 1),
								   pMeshData->calcEdgeLengthFromElementAndEdgeBoundaryPlanes(m_planeID, iElem, 2) };
		const std::complex<double> factor = std::complex<double>(0.0, omega * CommonParameters::mu);// exp(-i*omega*t) form
		double conductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcAnisotropicConductivityTensor(pResistivityBlock->getBlockIDFromElemID(elemID), conductivityTensor);
		for (int iEdge1 = 0; iEdge1 < 3; ++iEdge1) {
			const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 3; ++iEdge2) {
				const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
				if (col <= DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1(0.0);
				double integral2(0.0);
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					integral1 += getShapeFuncRotated() * getShapeFuncRotated() * m_weight[ip];
					const double Nh = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob22
									- getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob12;
					const double Nv = getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob11
									- getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge2) * jacobMat.jacob21;
					if (m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus) {//YZ Plane
						integral2 += ((getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob22
									- getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob12)
									* (Nh * conductivityTensor[1][1] + Nv * conductivityTensor[1][2])
									+ (getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob11
									- getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob21)
									* (Nh * conductivityTensor[2][1] + Nv * conductivityTensor[2][2])) * m_weight[ip];
					}
					else {//ZX Plane
						integral2 += ((getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob22
									- getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob12)
									* (Nh * conductivityTensor[0][0] + Nv * conductivityTensor[0][2])
									+ (getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob11
									- getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], iEdge1) * jacobMat.jacob21)
									* (Nh * conductivityTensor[2][0] + Nv * conductivityTensor[2][2])) * m_weight[ip];
					}
				}
				if (m_signInversion[iElem][iEdge1] != m_signInversion[iElem][iEdge2]) {
					integral1 *= -divDetJacob * length[iEdge1] * length[iEdge2];
					integral2 *= -divDetJacob * length[iEdge1] * length[iEdge2];
				}
				else {
					integral1 *= divDetJacob * length[iEdge1] * length[iEdge2];
					integral2 *= divDetJacob * length[iEdge1] * length[iEdge2];
				}
				std::complex<double> val = std::complex<double>(integral1, 0.0) - std::complex<double>(integral2, 0.0) * factor;// exp(-i*omega*t) form
				if (col == DIRICHLET_BOUNDARY_NONZERO_VALUE) {
					const std::complex<double> nonZeroValue = m_edgesIDGlobal2NonZeroValues[m_IDsLocal2Global[iElem][iEdge2]];
					m_matrix2DAnalysis.addRightHandSideVector(row, -val * nonZeroValue);// Add to right hand side vector
				}
				else if (col >= row) {// Store only upper triangle part
					m_matrix2DAnalysis.addNonZeroValues(row, col, val);// Add to matrix
				}
			}// iEdge2
		}// iEdge1		
	}// iElem

}


