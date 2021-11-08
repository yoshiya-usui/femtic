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
#include "Forward2DSquareElement1stOrderEdgeBased.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "ResistivityBlock.h"
#include <iostream>
#include <assert.h>

//// Defailt constructer
//Forward2DSquareElement1stOrderEdgeBased::Forward2DSquareElement1stOrderEdgeBased()
//{}

// Constructer
Forward2DSquareElement1stOrderEdgeBased::Forward2DSquareElement1stOrderEdgeBased( const int planeID, const int iPol ):
	Forward2DSquareElementEdgeBased( planeID, iPol )
{}

// Destructer
Forward2DSquareElement1stOrderEdgeBased::~Forward2DSquareElement1stOrderEdgeBased()
{}

// Calculate EM fields of boundary planes by 2D forward calculcation with 1st order edge-based element
void Forward2DSquareElement1stOrderEdgeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ){

	const int imode = calcMode();// TM or TE mode
	if ( imode != CommonParameters::TM_MODE ){
		OutputFiles::m_logFile << "Error : Only TM mode can be treated in calcEMFieldsOfBoundaryPlanes1stOrderEdgeBased ! imode = " << imode << "." << std::endl;
		exit(1);
	}
#ifdef _DEBUG_WRITE
	std::cout << "imode " << imode << std::endl;// For debug
#endif

	if( m_sourceFieldElectric == false ){
		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in calcEMFieldsOfBoundaryPlanes1stOrderEdgeBased !" << std::endl;
		exit(1);
	}
	
	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	if( pAnalysisControl->getBoundaryConditionBottom() != AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		OutputFiles::m_logFile << "Error : When 1st order edge-based element is used, electric field of the bottom must be zero !" << std::endl;
		exit(1);
	}

	if( m_specifyTEResultToSidesOfEdgeElement ){
		OutputFiles::m_logFile << "Error : Horizontal electric field cannot be specified at the left and right side the 2D model of 1st order edge-based element !" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Calculate electric field on a boundary plane with 1st order edge-based element." << std::endl;

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
	const int nEq = ( 2 * numElemW + 1 ) * ( 2 * numElemH ) + ( 2 * numElemH + 1 ) * ( 2 * numElemW );// Number of equations

	int nTmp = nEq - 4 * numElemW;// Exclude horizontal edges on the top and the bottom
	nTmp -= 4 * numElemH;// Exclude vertical edges of the left and right sides
	const int nEqDegenerated = nTmp;// Number of equations after degeneration

	OutputFiles::m_logFile << "# Number of equation = " << nEq
		<< ", Number of equation after degeneration = " << nEqDegenerated << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() <= 0 ){
		m_matrix2DAnalysis.setDegreeOfEquation( nEqDegenerated );
	}

	//----------------------------------------------------------------
	//--- Calculate array converting local node IDs to global ones ---
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
			const int offset = 8 * numElemW + 2;

			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2Global[iElem] = new int[12];

			m_IDsLocal2Global[iElem][0]  = ih * offset + 2 * iw;
			m_IDsLocal2Global[iElem][1]  = ih * offset + 2 * iw + 1;
			m_IDsLocal2Global[iElem][2]  = ih * offset + 2 * iw + 4 * numElemW + 1;
			m_IDsLocal2Global[iElem][3]  = ih * offset + 2 * iw + 4 * numElemW + 2;
			m_IDsLocal2Global[iElem][4]  = ih * offset + 2 * iw + 8 * numElemW + 2;  
			m_IDsLocal2Global[iElem][5]  = ih * offset + 2 * iw + 8 * numElemW + 3;  
			m_IDsLocal2Global[iElem][6]  = ih * offset + 2 * iw + 2 * numElemW;  
			m_IDsLocal2Global[iElem][7]  = ih * offset + 2 * iw + 6 * numElemW + 1;  
			m_IDsLocal2Global[iElem][8]  = ih * offset + 2 * iw + 2 * numElemW + 1;
			m_IDsLocal2Global[iElem][9]  = ih * offset + 2 * iw + 6 * numElemW + 2;
			m_IDsLocal2Global[iElem][10] = ih * offset + 2 * iw + 2 * numElemW + 2;
			m_IDsLocal2Global[iElem][11] = ih * offset + 2 * iw + 6 * numElemW + 3;  
		}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2Global : " << iElem;
			for( int i = 0; i < 12; ++i ){
				std::cout << " " << m_IDsLocal2Global[iElem][i];
			}
			std::cout << std::endl;
		}
#endif
	//----- debug <<<<<

		//---
		//--- Calculate array converting local node IDs to global ones after degeneration
		//---
		m_IDsLocal2GlobalDegenerated = new int*[nElem];

		for( int iElem = 0; iElem < nElem; ++iElem ){

			const int offset = 8 * numElemW - 2;

			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2GlobalDegenerated[iElem] = new int[12];

			if( iw == 0 ){// Elements at the left side
				m_IDsLocal2GlobalDegenerated[iElem][0]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][1]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][2]  =   ih       * offset + 2 * iw + 2 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][3]  =   ih       * offset + 2 * iw + 2 * numElemW;
				m_IDsLocal2GlobalDegenerated[iElem][4]  =   ih       * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][5]  =   ih       * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][6]  = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsLocal2GlobalDegenerated[iElem][7]  = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsLocal2GlobalDegenerated[iElem][8]  =   ih       * offset + 2 * iw;
				m_IDsLocal2GlobalDegenerated[iElem][9]  =   ih       * offset + 2 * iw + 4 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][10] =   ih       * offset + 2 * iw + 1;
				m_IDsLocal2GlobalDegenerated[iElem][11] =   ih       * offset + 2 * iw + 4 * numElemW;
			}else if( iw == numElemW - 1 ){// Elements at the right side
				m_IDsLocal2GlobalDegenerated[iElem][0]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][1]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][2]  =   ih       * offset + 2 * iw + 2 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][3]  =   ih       * offset + 2 * iw + 2 * numElemW;
				m_IDsLocal2GlobalDegenerated[iElem][4]  =   ih       * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][5]  =   ih       * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][6]  =   ih       * offset + 2 * iw - 1;
				m_IDsLocal2GlobalDegenerated[iElem][7]  =   ih       * offset + 2 * iw + 4 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][8]  =   ih       * offset + 2 * iw;
				m_IDsLocal2GlobalDegenerated[iElem][9]  =   ih       * offset + 2 * iw + 4 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][10] = DIRICHLET_BOUNDARY_ZERO_VALUE;
				m_IDsLocal2GlobalDegenerated[iElem][11] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			}else{
				m_IDsLocal2GlobalDegenerated[iElem][0]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][1]  = ( ih - 1 ) * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][2]  =   ih       * offset + 2 * iw + 2 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][3]  =   ih       * offset + 2 * iw + 2 * numElemW;
				m_IDsLocal2GlobalDegenerated[iElem][4]  =   ih       * offset + 2 * iw + 6 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][5]  =   ih       * offset + 2 * iw + 6 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][6]  =   ih       * offset + 2 * iw - 1;
				m_IDsLocal2GlobalDegenerated[iElem][7]  =   ih       * offset + 2 * iw + 4 * numElemW - 2;
				m_IDsLocal2GlobalDegenerated[iElem][8]  =   ih       * offset + 2 * iw;
				m_IDsLocal2GlobalDegenerated[iElem][9]  =   ih       * offset + 2 * iw + 4 * numElemW - 1;
				m_IDsLocal2GlobalDegenerated[iElem][10] =   ih       * offset + 2 * iw + 1;
				m_IDsLocal2GlobalDegenerated[iElem][11] =   ih       * offset + 2 * iw + 4 * numElemW;
			}
				
		}
		for( int iElem = 0; iElem < numElemW; ++iElem ){// Elements at the top
			m_IDsLocal2GlobalDegenerated[iElem][0] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
			m_IDsLocal2GlobalDegenerated[iElem][1] = DIRICHLET_BOUNDARY_NONZERO_VALUE;
		}
		for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){// Element at the bottom
			m_IDsLocal2GlobalDegenerated[iElem][4] = DIRICHLET_BOUNDARY_ZERO_VALUE;
			m_IDsLocal2GlobalDegenerated[iElem][5] = DIRICHLET_BOUNDARY_ZERO_VALUE;
		}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
			for( int i = 0; i < 12; ++i ){
				std::cout << " " << m_IDsLocal2GlobalDegenerated[iElem][i];
			}
			std::cout << std::endl;
		}
#endif
	//----- debug <<<<<

		m_hasAlreadySetIDsLocal2Global = true;

	}
	
	//-------------------------------------------------------------------
	//--- Calculate non-zero electric field values specified to nodes ---
	//-------------------------------------------------------------------
	std::map<int, std::complex<double> > nodesGlobal2NonZeroValues;
	//--- Top of the boundary ( source field )
	const double sourceValueElectric = CommonParameters::sourceValueElectric;
	//const double sourceValueMagnetic = 1;
	for( int iElem = 0; iElem < numElemW; ++iElem ){// Elements at the top
		nodesGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][0], std::complex<double>(sourceValueElectric, 0.0) ) );
		nodesGlobal2NonZeroValues.insert( std::map<int, std::complex<double> >::value_type( m_IDsLocal2Global[iElem][1], std::complex<double>(sourceValueElectric, 0.0) ) );
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
			for( int iNode1 = 0; iNode1 < 12; ++iNode1 ){
				const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 12; ++iNode2 ){
					const int col = m_IDsLocal2GlobalDegenerated[iElem][iNode2];
					if( col < 0 ){
						continue;
					}

					if( col >= row ){// Store only upper triangle part
						m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
					}

				}// iNode2
			}// iNode1
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
	const int numGauss = 3;
	const int numIntegralPoints = 9;
	double wLocal[9]    = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double hLocal[9]    = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double weights2D[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int ip(0);
	for( int i = 0; i < numGauss; ++i ){
		for( int j = 0; j < numGauss; ++j ){
			wLocal[ip] = CommonParameters::abscissas3Point[i];
			hLocal[ip] = CommonParameters::abscissas3Point[j];
			weights2D[ip] = CommonParameters::weights3Point[i] * CommonParameters::weights3Point[j];
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

		//--- Calculate omega * mu * sigma
		const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
		//const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
		const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
		const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);
		//const std::complex<double> factor = std::complex<double>( omega * omega * CommonParameters::mu * CommonParameters::epsilon, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form
		const std::complex<double> factor = std::complex<double>( 0.0, omega * CommonParameters::mu * sigma );// exp(-i*omega*t) form

		for( int iNode1 = 0; iNode1 < 12; ++iNode1 ){

			const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
			if( row < 0 ){
				continue;
			}

			for( int iNode2 = 0; iNode2 < 12; ++iNode2 ){

				const int col = m_IDsLocal2GlobalDegenerated[iElem][iNode2];
				if( col <= DIRICHLET_BOUNDARY_ZERO_VALUE ){
					continue;
				}

				double integral1(0);
				double integral2(0);
				for( ip = 0; ip < numIntegralPoints; ++ip ){
					integral1 +=   getShapeFuncRotated1stOrderEdgeBased( wLocal[ip], hLocal[ip], width, height, iNode1 )
					             * getShapeFuncRotated1stOrderEdgeBased( wLocal[ip], hLocal[ip], width, height, iNode2 )
							     * weights2D[ip];
					integral2 += ( getShapeFuncHorizontal1stOrderEdgeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncHorizontal1stOrderEdgeBased( wLocal[ip], hLocal[ip], iNode2 )
						         + getShapeFuncVertical1stOrderEdgeBased(   wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncVertical1stOrderEdgeBased(   wLocal[ip], hLocal[ip], iNode2 ) )
								 * weights2D[ip];
				}
				integral1 *= jacobian;
				integral2 *= jacobian;

				std::complex<double> val = std::complex<double>( integral1 , 0.0 )
				                         - std::complex<double>( integral2 , 0.0 ) * factor;// exp(-i*omega*t) form
									
				if( col == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					const std::complex<double> nonZeroValue = nodesGlobal2NonZeroValues[ m_IDsLocal2Global[iElem][iNode2] ];
					m_matrix2DAnalysis.addRightHandSideVector( row, -val * nonZeroValue );// Add to right hand side vector
				}else if( col >= row ){// Store only upper triangle part
					m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				}

			}// iNode2
		}// iNode1		

	}// iElem

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
		for( int iNode = 0; iNode < 12; ++iNode ){

			const int iNum = m_IDsLocal2Global[iElem][iNode];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsLocal2GlobalDegenerated[iElem][iNode];

				if( iNumDegenerated == DIRICHLET_BOUNDARY_ZERO_VALUE ){
					m_solution[iNum] = std::complex<double>( 0.0, 0.0 );
				}else if( iNumDegenerated == DIRICHLET_BOUNDARY_NONZERO_VALUE ){
					m_solution[iNum] = nodesGlobal2NonZeroValues[ iNum ];
				}else{
					m_solution[iNum] = solutionDegenerated[ iNumDegenerated ]; 
				}

				alreadyFound[ iNum ] = true;
			}
			
		}// iNode
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
//	output2DResult( Forward2DSquareElement::EDGE_BASED_FIRST_ORDER, iPol, freq, nElem, numElemW );
//#endif
//	//----- debug <<<<<
	output2DResult( Forward2DSquareElement::EDGE_BASED_FIRST_ORDER, freq, nElem, numElemW, pMeshDataBrickElement );

	//--- Release memory
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] nodesLocal2Global[iElem];
	//}
	//delete [] nodesLocal2Global;
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete [] nodesLocal2GlobalDegenerated[iElem];
	//}
	//delete [] nodesLocal2GlobalDegenerated;
	nodesGlobal2NonZeroValues.clear();

}

// Get shape functions of horizontal direction for 1st order edge-based elements
inline double Forward2DSquareElement1stOrderEdgeBased::getShapeFuncHorizontal1stOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const{

	const double div3 = 1.0 / 3.0;

	switch( num ){
		case 0:
			return - 0.75 *  hLocal        * ( hLocal - 1 ) * ( wLocal - div3 );
			break;
		case 1:
			return   0.75 *  hLocal        * ( hLocal - 1 ) * ( wLocal + div3 );
			break;
		case 2:
			return   1.50 * ( hLocal - 1 ) * ( hLocal + 1 ) * ( wLocal - div3 );
			break;
		case 3:
			return - 1.50 * ( hLocal - 1 ) * ( hLocal + 1 ) * ( wLocal + div3 );
			break;
		case 4:
			return - 0.75 *   hLocal       * ( hLocal + 1 ) * ( wLocal - div3 );
			break;
		case 5:
			return   0.75 *   hLocal       * ( hLocal + 1 ) * ( wLocal + div3 );
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
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncHorizontal1stOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}
}

// Get shape functions of vertical direction for 1st order edge-based elements
inline double Forward2DSquareElement1stOrderEdgeBased::getShapeFuncVertical1stOrderEdgeBased( const double wLocal, const double hLocal, const int num ) const{

	const double div3 = 1.0 / 3.0;

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
			return - 0.75 *  wLocal        * ( wLocal - 1 ) * ( hLocal - div3 );
			break;
		case 7:
			return   0.75 *  wLocal        * ( wLocal - 1 ) * ( hLocal + div3 );
			break;
		case 8:
			return   1.50 * ( wLocal - 1 ) * ( wLocal + 1 ) * ( hLocal - div3 );
			break;
		case 9:
			return - 1.50 * ( wLocal - 1 ) * ( wLocal + 1 ) * ( hLocal + div3 );
			break;
		case 10:
			return - 0.75 *   wLocal       * ( wLocal + 1 ) * ( hLocal - div3 );
			break;
		case 11:
			return   0.75 *   wLocal       * ( wLocal + 1 ) * ( hLocal + div3 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncVertical1stOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}


// Get shape functions rotated for 1st order edge-based elements
inline double Forward2DSquareElement1stOrderEdgeBased::getShapeFuncRotated1stOrderEdgeBased( const double wLocal, const double hLocal, const double width, const double height, const int num ) const{

	const double div3 = 1.0 / 3.0;

	switch( num ){
		case 0:
			return   3.0 / height * ( hLocal - 0.5 ) * ( wLocal - div3 );
			break;
		case 1:
			return - 3.0 / height * ( hLocal - 0.5 ) * ( wLocal + div3 );
			break;
		case 2:
			return - 6.0 / height *  hLocal          * ( wLocal - div3 );
			break;
		case 3:
			return   6.0 / height *  hLocal          * ( wLocal + div3 );
			break;
		case 4:
			return   3.0 / height * ( hLocal + 0.5 ) * ( wLocal - div3 );
			break;
		case 5:
			return - 3.0 / height * ( hLocal + 0.5 ) * ( wLocal + div3 );
			break;
		case 6:
			return - 3.0 / width  * ( wLocal - 0.5 ) * ( hLocal - div3 );
			break;
		case 7:
			return   3.0 / width  * ( wLocal - 0.5 ) * ( hLocal + div3 );
			break;
		case 8:
			return   6.0 / width  *   wLocal         * ( hLocal - div3 );
			break;
		case 9:
			return - 6.0 / width  *   wLocal         * ( hLocal + div3 );
			break;
		case 10:
			return - 3.0 / width  * ( wLocal + 0.5 ) * ( hLocal - div3 );
			break;
		case 11:
			return   3.0 / width  * ( wLocal + 0.5 ) * ( hLocal + div3 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncVertical1stOrderEdgeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Calculate horizontal electric field values for 1st order edge-based elements
// [Input] : 1) iPol   : ID of polarization
//           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
//	 	 	 3) wLocal : Local coordinate of horizontal direction of the element
//			 4) hLocal : Local coordinate of vertical direction of the element
// [Output] : Horizontal electric field calculated
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueElectricFieldHorizontal( const int iElem, const double wLocal, const double hLocal ) const{

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
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncHorizontal1stOrderEdgeBased( wLocal, hLocal, i );
	}

	return val;
}

// Calculate vertical electric field values for 1st order edge-based elements
// [Input] : 1) iPol   : ID of polarization
//           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
//	 	 	 3) wLocal : Local coordinate of horizontal direction of the element
//			 4) hLocal : Local coordinate of vertical direction of the element
// [Output] : Vertical electric field calculated
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueElectricFieldVertical( const int iElem, const double wLocal, const double hLocal ) const{

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
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncVertical1stOrderEdgeBased( wLocal, hLocal, i );
	}

	return val;
}

// Calculate magnetic field values for 1st order edge-based elements
// [Input] : 1) iPol   : ID of polarization
//           2) freq   : Frequency
//           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
//	 	 	 4) wLocal : Local coordinate of horizontal direction of the element
//			 5) hLocal : Local coordinate of vertical direction of the element
// [Output] : Magnetic field calculated
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueMagneticFieldPerpendicular( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

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

	std::complex<double> val(0.0, 0.0);
	for( int i = 0; i < 12; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncRotated1stOrderEdgeBased( wLocal, hLocal, width, height, i );
	}

	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
	const double factor = omega * CommonParameters::mu;

	val /= std::complex<double>(0.0, factor);

	return val;

}

// Calculate parameter V of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldPerpendicular( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate parameter J of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldHorizontal( iElem, wLocal, hLocal );
}

// Calculate parameter I of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderEdgeBased::calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldVertical( iElem, wLocal, hLocal );
}
