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
#include "Forward2DSquareElement0thOrderEdgeBased.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "ResistivityBlock.h"
#include <iostream>
#include <assert.h>

//// Defailt constructer
//Forward2DSquareElement0thOrderEdgeBased::Forward2DSquareElement0thOrderEdgeBased()
//{}

// Constructer
Forward2DSquareElement0thOrderEdgeBased::Forward2DSquareElement0thOrderEdgeBased( const int planeID, const int iPol ):
	Forward2DSquareElementEdgeBased( planeID, iPol )
{}

// Destructer
Forward2DSquareElement0thOrderEdgeBased::~Forward2DSquareElement0thOrderEdgeBased()
{}

// Calculate EM fields of boundary planes by 2D forward calculcation with 0th order edge-based element
void Forward2DSquareElement0thOrderEdgeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ){

	const int imode = calcMode();// TM or TE mode
	if ( imode != CommonParameters::TM_MODE ){
		OutputFiles::m_logFile << "Error : Only TM mode can be treated in calcEMFieldsOfBoundaryPlanes0thOrderEdgeBased ! imode = " << imode << "." << std::endl;
		exit(1);
	}
#ifdef _DEBUG_WRITE
	std::cout << "imode " << imode << std::endl;// For debug
#endif

	if( m_sourceFieldElectric == false ){
		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in calcEMFieldsOfBoundaryPlanes0thOrderEdgeBased !" << std::endl;
		exit(1);
	}
	
	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	if( pAnalysisControl->getBoundaryConditionBottom() != AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		OutputFiles::m_logFile << "Error : When 0th order edge-based element is used, electric field of the bottom must be zero !" << std::endl;
		exit(1);
	}

	if( m_specifyTEResultToSidesOfEdgeElement ){
		OutputFiles::m_logFile << "Error : Horizontal electric field cannot be specified at the left and right side the 2D model of 0th order edge-based element !" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Calculate electric field on a boundary plane with 0th order edge-based element." << std::endl;

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();

	//--- Calculate element division of 2D
	const int numElemW = calcNumElemHorizontal( pMeshDataBrickElement );
	const int numElemH = calcNumElemVertical( pMeshDataBrickElement );
#ifdef _DEBUG_WRITE
	std::cout << "numElemH numElemW " << numElemH << " " << numElemW << std::endl; // For debug
#endif

	//--- Total element number on the boudary plane
	//const int nElem = pMeshDataBrickElement->m_numElemOnBoundaryPlanes[ m_planeID ]; // Total number of elements on the plane
	const int nElem = pMeshDataBrickElement->getNumElemOnBoundaryPlanes( m_planeID ); // Total number of elements on the plane
#ifdef _DEBUG_WRITE
	std::cout << "nElem " << nElem << std::endl;// For debug
#endif

	//----------------------------------
	//--- Set Number of the equation ---
	//----------------------------------
	const int nEq = numElemW * ( numElemH + 1 ) + numElemH * ( numElemW + 1 );// Number of equations
	//const int nAirLayer = pMeshDataBrickElement->m_numAirLayer; // Number of the air layer
	const int nAirLayer = pMeshDataBrickElement->getNumAirLayer(); // Number of the air layer

	int nTmp = nEq - 2 * numElemW;// Exclude horizontal edges on the top and the bottom
	nTmp -= 2 * numElemH;// Exclude vertical edges of the left and right sides
	if( m_specifyTEResultToSidesOfEdgeElement ){// Horizontal electric field is specified at the left and right side
		nTmp -= 2 * ( numElemH - 1 );// Exclude horizontal edges of the left and right sides
	}
	const int nEqDegenerated = nTmp;// Number of equations after degeneration

	OutputFiles::m_logFile << "# Number of equation = " << nEq
		<< ", Number of equation after degeneration = " << nEqDegenerated << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() <= 0 ){
		m_matrix2DAnalysis.setDegreeOfEquation( nEqDegenerated );
	}

	//----------------------------------------------------------------
	//--- Calculate array converting local edge IDs to global ones ---
	//----------------------------------------------------------------
	const int DIRICHLET_BOUNDARY_NONZERO_VALUE = -1;
	const int DIRICHLET_BOUNDARY_ZERO_VALUE    = -2;

	if( m_hasAlreadySetIDsLocal2Global == false ){

		//---
		//--- Calculate array converting local edge IDs to global ones 
		//---
		if( m_IDsLocal2Global != NULL ){// Release memory
			const int num = sizeof( m_IDsLocal2Global ) / sizeof( m_IDsLocal2Global[0] );
			for( int i = 0; i < num; ++i ){
				delete [] m_IDsLocal2Global[i];
			}
			delete [] m_IDsLocal2Global;
			m_IDsLocal2Global = NULL;
		}
		m_IDsLocal2Global = new int*[nElem];

		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int offset = 2 * numElemW + 1;
			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2Global[iElem] = new int[4];

			m_IDsLocal2Global[iElem][0] = ( ih + 1 ) * offset + iw;
			m_IDsLocal2Global[iElem][1] =   ih       * offset + numElemW + iw;
			m_IDsLocal2Global[iElem][2] =   ih       * offset + iw;
			m_IDsLocal2Global[iElem][3] =   ih       * offset + numElemW + iw + 1;
		}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2Global : " << iElem;
			for( int i = 0; i < 4; ++i ){
				std::cout << " " << m_IDsLocal2Global[iElem][i];
			}
			std::cout << std::endl;
		}
#endif
	//----- debug <<<<<

		//---
		//--- Calculate array converting local edge IDs to global ones after degeneration
		//---
		if( m_IDsLocal2GlobalDegenerated != NULL ){// Release memory
			const int num = sizeof( m_IDsLocal2GlobalDegenerated ) / sizeof( m_IDsLocal2GlobalDegenerated[0] );
			for( int i = 0; i < num; ++i ){
				delete [] m_IDsLocal2GlobalDegenerated[i];
			}
			delete [] m_IDsLocal2GlobalDegenerated;
			m_IDsLocal2GlobalDegenerated = NULL;
		}		
		m_IDsLocal2GlobalDegenerated = new int*[nElem];

		for( int iElem = 0; iElem < nElem; ++iElem ){

			int offset = 2 * numElemW - 1;
			if( m_specifyTEResultToSidesOfEdgeElement ){// Horizontal electric field is specified at the left and right side
				offset -= 2;
			}

			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2GlobalDegenerated[iElem] = new int[4];

			if( iw == 0 ){// Elements at the left side
				if( m_specifyTEResultToSidesOfEdgeElement ){// Horizontal electric field is specified at the left and right side
					m_IDsLocal2GlobalDegenerated[iElem][0] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
					m_IDsLocal2GlobalDegenerated[iElem][2] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				}else{
					m_IDsLocal2GlobalDegenerated[iElem][0] = iw +   ih       * offset + numElemW - 1;
					m_IDsLocal2GlobalDegenerated[iElem][2] = iw + ( ih - 1 ) * offset + numElemW - 1;
				}
				m_IDsLocal2GlobalDegenerated[iElem][1] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsLocal2GlobalDegenerated[iElem][3] =   iw       + ih * offset;
			}else if( iw == numElemW - 1 ){// Elements at the right side
				if( m_specifyTEResultToSidesOfEdgeElement ){// Horizontal electric field is specified at the left and right side
					m_IDsLocal2GlobalDegenerated[iElem][0] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
					m_IDsLocal2GlobalDegenerated[iElem][2] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
				}else{
					m_IDsLocal2GlobalDegenerated[iElem][0] = iw +   ih       * offset + numElemW - 1;
					m_IDsLocal2GlobalDegenerated[iElem][2] = iw + ( ih - 1 ) * offset + numElemW - 1;
				}
				m_IDsLocal2GlobalDegenerated[iElem][1] = ( iw - 1 ) + ih * offset;
				m_IDsLocal2GlobalDegenerated[iElem][3] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			}else{
				m_IDsLocal2GlobalDegenerated[iElem][0] = iw +   ih       * offset + numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][1] = ( iw - 1 ) + ih * offset;
				m_IDsLocal2GlobalDegenerated[iElem][2] = iw + ( ih - 1 ) * offset + numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][3] =   iw       + ih * offset;
				if( m_specifyTEResultToSidesOfEdgeElement ){// Horizontal electric field is specified at the left and right side
					m_IDsLocal2GlobalDegenerated[iElem][0] -= 1;
					m_IDsLocal2GlobalDegenerated[iElem][2] += 1;
				}
			}
				
		}
		for( int iElem = 0; iElem < numElemW; ++iElem ){// Elements at the top
			m_IDsLocal2GlobalDegenerated[iElem][2] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
		}
		for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){// Element at the bottom
			m_IDsLocal2GlobalDegenerated[iElem][0] = DIRICHLET_BOUNDARY_ZERO_VALUE;
		}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
			for( int i = 0; i < 4; ++i ){
				std::cout << " " << m_IDsLocal2GlobalDegenerated[iElem][i];
			}
			std::cout << std::endl;
		}
#endif
	//----- debug <<<<<

		m_hasAlreadySetIDsLocal2Global = true;
	}
	
	//-------------------------------------------------------------------
	//--- Calculate non-zero electric field values specified to edges ---
	//-------------------------------------------------------------------
	std::map<int, std::complex<double> > edgesGlobal2NonZeroValues;

	//--- Top of the boundary ( source field )
	const double sourceValueElectric = CommonParameters::sourceValueElectric;
	//const double sourceValueMagnetic = 1;
	for( int iElem = 0; iElem < numElemW; ++iElem ){// Elements at the top
		edgesGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][2], std::complex<double>(sourceValueElectric, 0.0) ) );
	}
	
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

		////---- Output memory required to caluculate the 2D analysis
		//m_matrix2DAnalysis.writeMemoryRequiredByMatrixSolver(); 
		
		m_hasMatrixStructureSetAndAnalyzed = true;

	}

	//---------------------------------------------------------------------------
	//--- Calculate integral points and weights of two point Gauss quadrature ---
	//---------------------------------------------------------------------------
	const int numGauss = 2;
	const int numIntegralPoints = 4;
	double wLocal[4]    = { 0, 0, 0, 0 };
	double hLocal[4]    = { 0, 0, 0, 0 };
	double weights2D[4] = { 0, 0, 0, 0 };
	int ip(0);
	for( int i = 0; i < numGauss; ++i ){
		for( int j = 0; j < numGauss; ++j ){
			wLocal[ip] = CommonParameters::abscissas2Point[i];
			hLocal[ip] = CommonParameters::abscissas2Point[j];
			weights2D[ip] = CommonParameters::weights2Point[i] * CommonParameters::weights2Point[j];
			++ip;
		}
	}

	//-------------------------------------------------------
	//--- Set values of matrix and right hand side vector ---
	//-------------------------------------------------------
	OutputFiles::m_logFile << "# Set values of matrix and right hand side vector. " << pAnalysisControl->outputElapsedTime() << std::endl;
	ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	m_matrix2DAnalysis.zeroClearNonZeroValues();// Zero clear matrix values
	m_matrix2DAnalysis.zeroClearRightHandSideVector();// Zero clear right hand side vector

	for( int iElem = 0; iElem < nElem; ++iElem ){

		//--- Calculate Jacobian
		const double width  = calcWidth( iElem, pMeshDataBrickElement );
		const double height = calcHeight( iElem, pMeshDataBrickElement );
		const double jacobian = 0.25 * width * height;
		const double EpsilonDiffByX = 2.0 / width;
		const double EtaDiffByY     = 2.0 / height;
		const double translationFactor[4] = { EtaDiffByY, EpsilonDiffByX, EtaDiffByY, EpsilonDiffByX };

		//--- Calculate omega * mu * sigma
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		//const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
		const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		//const double factor = omega * CommonParameters::mu * sigma;
		//const std::complex<double> factor = std::complex<double>( omega * omega * CommonParameters::mu * CommonParameters::epsilon, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form
		const std::complex<double> factor = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form

		for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

			const int row = m_IDsLocal2GlobalDegenerated[iElem][iEdge1];
			if( row < 0 ){
				continue;
			}

			for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

				const int col = m_IDsLocal2GlobalDegenerated[iElem][iEdge2];
				if( col <= DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}

				double integral1(0);
				double integral2(0);
				for( ip = 0; ip < numIntegralPoints; ++ip ){
					integral1 += ( getShapeFuncRotated0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge1 ) * translationFactor[iEdge1]
					             * getShapeFuncRotated0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge2 ) * translationFactor[iEdge2]
							     * weights2D[ip] );
					integral2 += (   getShapeFuncHorizontal0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge1 ) * getShapeFuncHorizontal0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge2 )
						           + getShapeFuncVertical0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge1 ) * getShapeFuncVertical0thOrderEdgeBased( wLocal[ip], hLocal[ip], iEdge2 ) ) * weights2D[ip];
				}
				integral1 *= jacobian;
				integral2 *= jacobian;

				//std::complex<double> val = std::complex<double>( integral1 , 0.0 )
				//                         - std::complex<double>( 0.0 , integral2 * factor );// exp(-i*omega*t) form
				std::complex<double> val = std::complex<double>( integral1 , 0.0 )
				                         - std::complex<double>( integral2 , 0.0 ) * factor;// exp(-i*omega*t) form
				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					const std::complex<double> nonZeroValue = edgesGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iEdge2] ];
					m_matrix2DAnalysis.addRightHandSideVector( row, -val * nonZeroValue );// Add to right hand side vector
				}else if( col >= row ){// Store only upper triangle part
					m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				}

			}// iEdge2
		}// iEdge1		

	}// iElem

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
	std::complex<double>* solutionDegenerated = new std::complex<double>[ nEqDegenerated ];
	OutputFiles::m_logFile << "# Solve phase of matrix solver for boundary plane."
			<<  " Polarization : " << m_polarization << " Plane : " << m_planeID <<  ". " 
			<< pAnalysisControl->outputElapsedTime() << std::endl;
	m_matrix2DAnalysis.solvePhaseMatrixSolver( solutionDegenerated );//Solve phase of matrix solver

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
	for( int i = 0; i < nEqDegenerated; ++i ){
		std::cout << "i solutionDegenerated : " << i << " " << solutionDegenerated[i] << std::endl;
	}
#endif
	//----- debug <<<<<

	//std::complex<double>* solution = new std::complex<double>[ nEq ];
	if( m_solution != NULL ){
		delete [] m_solution;
		m_solution = NULL;
	}
	m_solution = new std::complex<double>[ nEq ];
				
	bool* alreadyFound = new bool[ nEq ];
	for( int i = 0; i < nEq; ++i ){
		alreadyFound[i] = false;
	}
	for( int iElem = 0; iElem < nElem; ++iElem ){
		for( int iEdge = 0; iEdge < 4; ++iEdge ){

			const int iNum = m_IDsLocal2Global[iElem][iEdge];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsLocal2GlobalDegenerated[iElem][iEdge];

				if( iNumDegenerated == DIRICHLET_BOUNDARY_ZERO_VALUE ){
					//solution[ iNum ] = std::complex<double>( 0.0, 0.0 );
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					//solution[ iNum ] = edgesGlobal2NonZeroValues[ iNum ];
					m_solution[iNum] = edgesGlobal2NonZeroValues[ iNum ];
				}else{
					//solution[ iNum ] = solutionDegenerated[ iNumDegenerated ];
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
	for( int i = 0; i < nEq; ++i ){
		std::cout << "i m_solution : " << i << " " << m_solution[i] << std::endl;
	}
#endif
	//----- debug <<<<<

//	//----- debug >>>>>
//#ifdef _OUTPUT_2D_RESULT
//	output2DResult( Forward2DSquareElement::EDGE_BASED_ZEROTH_ORDER, iPol, freq, nElem, numElemW );
//#endif
//	//----- debug <<<<<
	output2DResult( Forward2DSquareElement::EDGE_BASED_ZEROTH_ORDER, freq, nElem, numElemW, pMeshDataBrickElement );

	//--- Release memory
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] edgesLocal2Global[iElem];
	//}
	//delete [] edgesLocal2Global;
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] edgesLocal2GlobalDegenerated[iElem];
	//}
	//delete [] edgesLocal2GlobalDegenerated;
	edgesGlobal2NonZeroValues.clear();

}

// Get shape functions of horizontal direction for 0th order edge-based elements
inline double Forward2DSquareElement0thOrderEdgeBased::getShapeFuncHorizontal0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return 0.5 * ( 1 + hLocal );
			break;
		case 1:
			return 0.0;
			break;
		case 2:
			return 0.5 * ( 1 - hLocal );
			break;
		case 3:
			return 0.0;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncHorizontal0thOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of vertical direction for 0th order edge-based elements
inline double Forward2DSquareElement0thOrderEdgeBased::getShapeFuncVertical0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return 0.0;
			break;
		case 1:
			return 0.5 * ( 1 - wLocal );
			break;
		case 2:
			return 0.0;
			break;
		case 3:
			return 0.5 * ( 1 + wLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncVertical0thOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions rotated for 0th order edge-based elements
inline double Forward2DSquareElement0thOrderEdgeBased::getShapeFuncRotated0thOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return - 0.5;
			break;
		case 1:
			return - 0.5;
			break;
		case 2:
			return   0.5;
			break;
		case 3:
			return   0.5;
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncRotated0thOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Calculate horizontal electric field values for 0th order edge-based elements
// [Input] : 1) iPol   : ID of polarization
//           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
//	 	 	 3) wLocal : Local coordinate of horizontal direction of the element
//			 4) hLocal : Local coordinate of vertical direction of the element
// [Output] : Horizontal electric field calculated
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueElectricFieldHorizontal( const int iElem, const double wLocal, const double hLocal ) const{

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
	for( int i = 0; i < 4; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncHorizontal0thOrderEdgeBased( wLocal, hLocal, i );
	}

	return val;
}

// Calculate vertical electric field values for 0th order edge-based elements
// [Input] : 1) iPol   : ID of polarization
//           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
//	 	 	 3) wLocal : Local coordinate of horizontal direction of the element
//			 4) hLocal : Local coordinate of vertical direction of the element
// [Output] : Vertical electric field calculated
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueElectricFieldVertical( const int iElem, const double wLocal, const double hLocal ) const{

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
	for( int i = 0; i < 4; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncVertical0thOrderEdgeBased( wLocal, hLocal, i );
	}

	return val;
}

//// Calculate Magnetic field values for 0th order edge-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) freq   : Frequency
////           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 4) wLocal : Local coordinate of horizontal direction of the element
////			 5) hLocal : Local coordinate of vertical direction of the element
//// [Output] : Magnetic field calculated
//std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueMagneticField( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
//
//	if( m_solution == NULL ){
//		OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
//		exit(1);
//	}else if( m_IDsLocal2Global[iElem] == NULL ){
//		OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
//		exit(1);
//	}
//
//	const double width  = calcWidth( iElem, pMeshDataBrickElement );
//	const double height = calcHeight( iElem, pMeshDataBrickElement );
//	const double EpsilonDiffByX = 2.0 / width;
//	const double EtaDiffByY     = 2.0 / height;
//	const double translationFactor[4] = { EtaDiffByY, EpsilonDiffByX, EtaDiffByY, EpsilonDiffByX };
//
//	std::complex<double> val(0.0, 0.0);
//	for( int i = 0; i < 4; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotated0thOrderEdgeBased( wLocal, hLocal, i ) * translationFactor[i];
//	}
//
//	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
//	const double factor = omega * CommonParameters::mu;
//
//	val /= std::complex<double>(0.0, factor);
//
//	return val;
//}

// Calculate magnetic field perpendicular to the boundary plane
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	//if( m_solution == NULL ){
	//	OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
	//	exit(1);
	//}else if( m_IDsLocal2Global[iElem] == NULL ){
	//	OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
	//	exit(1);
	//}
	assert( m_solution != NULL );
	assert( m_IDsLocal2Global[iElem] != NULL );

	const double width  = calcWidth( iElem, pMeshDataBrickElement );
	const double height = calcHeight( iElem, pMeshDataBrickElement );
	const double EpsilonDiffByX = 2.0 / width;
	const double EtaDiffByY     = 2.0 / height;
	const double translationFactor[4] = { EtaDiffByY, EpsilonDiffByX, EtaDiffByY, EpsilonDiffByX };

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 4; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotated0thOrderEdgeBased( wCoord, hCoord, i ) * translationFactor[i];
	}

	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
	const double factor = omega * CommonParameters::mu;

	val /= std::complex<double>(0.0, factor);

	return val;

}

// Calculate parameter V of Rodi(1976)
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldPerpendicular( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate parameter J of Rodi(1976)
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldHorizontal( iElem, wLocal, hLocal );
}

// Calculate parameter I of Rodi(1976)
std::complex<double> Forward2DSquareElement0thOrderEdgeBased::calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldVertical( iElem, wLocal, hLocal );
}

