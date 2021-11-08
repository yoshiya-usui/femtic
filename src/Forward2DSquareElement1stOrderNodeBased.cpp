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
#include "Forward2DSquareElement1stOrderNodeBased.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "ResistivityBlock.h"
#include <iostream>
#include <assert.h>

//// Defailt constructer
//Forward2DSquareElement1stOrderNodeBased::Forward2DSquareElement1stOrderNodeBased()
//{}

// Constructer
Forward2DSquareElement1stOrderNodeBased::Forward2DSquareElement1stOrderNodeBased( const int planeID, const int iPol ):
	Forward2DSquareElementNodeBased( planeID, iPol )
{}

// Destructer
Forward2DSquareElement1stOrderNodeBased::~Forward2DSquareElement1stOrderNodeBased()
{}

// Calculate EM fields of boundary planes by 2D forward calculcation with 1st order node-based element
void Forward2DSquareElement1stOrderNodeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ){
	
	const int imode = calcMode();// TM or TE mode
	if ( imode != CommonParameters::TE_MODE ){
		OutputFiles::m_logFile << "Error : Only TE mode can be treated in calcEMFieldsOfBoundaryPlanes1stOrderNodeBased ! imode = " << imode << "." << std::endl;
		exit(1);
	}
#ifdef _DEBUG_WRITE
	std::cout << "imode " << imode << std::endl;// For debug
#endif

	if( m_sourceFieldElectric == false ){
		OutputFiles::m_logFile << "Error : Electric field must be specified at the top of the model as source in calcEMFieldsOfBoundaryPlanes1stOrderNodeBased !" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Calculating electric field on a boundary plane with 1st order node-based element." << std::endl;

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	
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
	//--- Set degree of the equation ---
	//----------------------------------
	const int nEq = ( numElemW + 1 ) * ( numElemH + 1 );// Number of equations

	int nTmp = ( numElemW + 1 ) * numElemH;
	if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		// If perfect conductor is specified at the bottom, exclude nodes of bottom side on which Dirichlet condition is given
		nTmp -= ( numElemW + 1 );
	}
	const int nEqDegenerated = nTmp;// Number of equations after excludomg nodes of top side

	OutputFiles::m_logFile << "# Number of equation = " << nEq
		<< ", Number of equation after degeneration = " << nEqDegenerated << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() <= 0 ){
		m_matrix2DAnalysis.setDegreeOfEquation( nEqDegenerated );
	}

	//----------------------------------------------------------------
	//--- Calculate array converting local node IDs to global ones ---
	//----------------------------------------------------------------
	const int DEGENERATED_NODES_ON_TOP_SIDE = -1;
	const int DEGENERATED_NODES_ON_BOTTOM_SIDE = -2;

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

		const int offset = numElemW + 1;
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2Global[iElem] = new int[4];

			m_IDsLocal2Global[iElem][0] = ( ih + 1 ) * offset + iw + 1;
			m_IDsLocal2Global[iElem][1] = ( ih + 1 ) * offset + iw;
			m_IDsLocal2Global[iElem][2] =   ih       * offset + iw;
			m_IDsLocal2Global[iElem][3] =   ih       * offset + iw + 1;
		}

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

		for( int iElem = 0; iElem < nElem; ++iElem ){// Exclude nodes of top side because Dirichlet condition is set to the nodes
			m_IDsLocal2GlobalDegenerated[iElem] = new int[4];

			m_IDsLocal2GlobalDegenerated[iElem][0] = m_IDsLocal2Global[iElem][0] - ( numElemW + 1 );
			m_IDsLocal2GlobalDegenerated[iElem][1] = m_IDsLocal2Global[iElem][1] - ( numElemW + 1 );
			m_IDsLocal2GlobalDegenerated[iElem][2] = m_IDsLocal2Global[iElem][2] - ( numElemW + 1 );
			m_IDsLocal2GlobalDegenerated[iElem][3] = m_IDsLocal2Global[iElem][3] - ( numElemW + 1 );
		}

		for( int iElem = 0; iElem < numElemW; ++iElem ){
			// Exclude nodes of the top side on which Dirichlet condition is given
			m_IDsLocal2GlobalDegenerated[iElem][2] = DEGENERATED_NODES_ON_TOP_SIDE;
			m_IDsLocal2GlobalDegenerated[iElem][3] = DEGENERATED_NODES_ON_TOP_SIDE;
		}
		if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){// Perfect conductor is specified at the bottom,
			// Exclude nodes of the bottom side on which Dirichlet condition is given
			for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
				m_IDsLocal2GlobalDegenerated[iElem][0] = DEGENERATED_NODES_ON_BOTTOM_SIDE;
				m_IDsLocal2GlobalDegenerated[iElem][1] = DEGENERATED_NODES_ON_BOTTOM_SIDE;
			}
		}
	
	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
			for( int i = 0; i < 4; ++i ){
				std::cout << " " <<m_IDsLocal2GlobalDegenerated[iElem][i];
			}
			std::cout << std::endl;
		}
#endif
	//----- debug <<<<<

		m_hasAlreadySetIDsLocal2Global = true;

	}

	////-----------------------------------------------------------------------------------------------------
	////--- Calculate array converting the global node IDs to the ones of full (3D) model and its reverse ---
	////-----------------------------------------------------------------------------------------------------
	//int* nodesGlobal2FullModel = new int[nEq];
	//std::map<int, int> nodesFullModel2Global;
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	for( int i = 0; i < 4; ++i ){
	//		nodesGlobal2FullModel[ nodesLocal2Global[iElem][i] ] = pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4 + i];
	//		nodesFullModel2Global.insert( std::map<int, int>::value_type( pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4 + i], nodesLocal2Global[iElem][i] ) );
	//	}
	//}

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
			for( int iNode1 = 0; iNode1 < 4; ++iNode1 ){
				const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 4; ++iNode2 ){
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

		//-----------------------------------------------
		//--- Components due to the boudary condition ---
		//-----------------------------------------------

		//----- Don't remove because this might need to be restored in the future >>>>>
		////--- Components due to the boundary condition of the top side
		//const int nodesUpperSide[2] = { 2, 3 };
		//for( int iElem = 0; iElem < numElemW; ++iElem ){
		//	for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){
		//		for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){

		//			const int row = nodesLocal2Global[iElem][ nodesUpperSide[iNode1] ];
		//			const int col = nodesLocal2Global[iElem][ nodesUpperSide[iNode2] ];
		//			if( row < 0 || col < 0 ){
		//				continue;
		//			}
		//			if( col >= row ){// Store only upper triangle part
		//				m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
		//			}

		//		}// iNode2
		//	}// iNode1
		//}// iElem
		//----- Don't remove because this might need to be restored in the future <<<<<<

		if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){// One dimensional boundary condition

			//--- Components due to the boudary condition of the bottom side
			const int nodesBottomSide[2] = { 0, 1 };
			for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
				for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){
					const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];
					if( row < 0 ){
						continue;
					}

					for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){										
						const int col = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode2] ];					 
						if( col < 0 ){
							continue;
						}

						if( col >= row ){// Store only upper triangle part
							m_matrix2DAnalysis.setStructureByTripletFormat( row, col );
						}

					}// iNode2
				}// iNode1
			}// iElem

		}

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
	const double sourceValueElectric = CommonParameters::sourceValueElectric;

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
		const double constWLocal = height / width;
		const double constHLocal = width / height;

		//--- Calculate eta and gamma
		std::complex<double> eta = calcEta( imode, freq, iElem );
		std::complex<double> gamma = calcGamma( imode, freq, iElem );

		//const std::complex<double> val1 = gamma * std::complex<double>( height * width / 36 , 0.0 );
		//const std::complex<double> val2 = std::complex<double>( height / 6 / width , 0.0 ) / eta;
		//const std::complex<double> val3 = std::complex<double>( width / 6 / height , 0.0 ) / eta;
		//const std::complex<double> K1 = val1 * std::complex<double>( 4.0 , 0.0 ) + val2 * std::complex<double>( 2.0 , 0.0 ) + val3 * std::complex<double>( 2.0 , 0.0 );
		//const std::complex<double> K2 = val1 * std::complex<double>( 2.0 , 0.0 ) - val2 * std::complex<double>( 2.0 , 0.0 ) + val3;
		//const std::complex<double> K3 = val1 * std::complex<double>( 2.0 , 0.0 ) + val2                                     - val3 * std::complex<double>( 2.0 , 0.0 );
		//const std::complex<double> K4 = val1                                     - val2                                     - val3;
		//for( int i = 0; i < 4; ++i ){
		//	m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][i], nodesLocal2Global[iElem][i], K1 );// Add K1 to matrix
		//}
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][0], nodesLocal2Global[iElem][1], K2 );// Add K2 to matrix
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][2], nodesLocal2Global[iElem][3], K2 );// Add K2 to matrix
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][0], nodesLocal2Global[iElem][2], K3 );// Add K3 to matrix
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][1], nodesLocal2Global[iElem][3], K3 );// Add K3 to matrix
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][0], nodesLocal2Global[iElem][3], K4 );// Add K4 to matrix
		//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][1], nodesLocal2Global[iElem][2], K4 );// Add K4 to matrix

		//----- Don't remove because this might need to be restored in the future >>>>>
		//for( int iNode1 = 0; iNode1 < 4; ++iNode1 ){
		//	for( int iNode2 = 0; iNode2 < 4; ++iNode2 ){

		//		const int row = nodesLocal2Global[iElem][iNode1];
		//		const int col = nodesLocal2Global[iElem][iNode2];
		//			
		//		if( row < 0 || col < 0 ){
		//			continue;
		//		}

		//		if( col >= row ){// Store only upper triangle part
		//			double integral1(0);
		//			double integral2(0);
		//			double integral3(0);
		//			for( ip = 0; ip < numIntegralPoints; ++ip ){
		//				integral1 += getShapeFunc1stOrderNodeBased(             wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc1stOrderNodeBased(             wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
		//				integral2 += getShapeFuncDiffByWLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
		//				integral3 += getShapeFuncDiffByHLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
		//			}
		//			integral1 *= jacobian;
		//			integral2 *= constWLocal;
		//			integral3 *= constHLocal;

		//			std::complex<double> val = std::complex<double>( integral1             , 0.0 ) * gamma  
		//									 + std::complex<double>( integral2 + integral3 , 0.0 ) / eta;

		//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
		//		}
		//			
		//	}// iNode2
		//}// iNode1
		//----- Don't remove because this might need to be restored in the future <<<<<

		for( int iNode1 = 0; iNode1 < 4; ++iNode1 ){

			const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
			if( row < 0 ){
				continue;
			}

			for( int iNode2 = 0; iNode2 < 4; ++iNode2 ){

				const int col = m_IDsLocal2GlobalDegenerated[iElem][iNode2];
				if( col <= DEGENERATED_NODES_ON_BOTTOM_SIDE ){
					continue;
				}

				double integral1(0);
				double integral2(0);
				double integral3(0);
				for( ip = 0; ip < numIntegralPoints; ++ip ){
					integral1 += getShapeFunc1stOrderNodeBased(             wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc1stOrderNodeBased(             wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
					integral2 += getShapeFuncDiffByWLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
					integral3 += getShapeFuncDiffByHLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal1stOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
				}
				integral1 *= jacobian;
				integral2 *= constWLocal;
				integral3 *= constHLocal;

				std::complex<double> val = std::complex<double>( integral1             , 0.0 ) * gamma  
				                         + std::complex<double>( integral2 + integral3 , 0.0 ) / eta;
									
				if( col == DEGENERATED_NODES_ON_TOP_SIDE ){
					m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueElectric, 0.0) );// Add to right hand side vector
				}else if( col >= row ){// Store only upper triangle part
					m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				}
					
			}// iNode2
		}// iNode1
		
	}// iElem

	//------------------------------------------------
	//--- Components due to the boudary conditions ---
	//------------------------------------------------

	//----- Don't remove because this might need to be restored in the future >>>>>
	////--- Components due to the boudary condition of the top side
	//for( int iElem = 0; iElem < numElemW; ++iElem ){
	//		
	//	const int nodesUpperSide[2] = { 2, 3 };

	//	double alpha2(0);
	//	double beta2(0);
	//	if(m_sourceFieldElectric){// Electric field specified at the top of the model as source
	//		alpha2 = 1.0e+20;
	//		beta2 = 1.0e+20 * sourceValueElectric;				
	//	}else{// Magnetic field specified at the top of the model as source
	//		beta2 = -1.0 * sourceValueMagnetic;
	//	}

	//	//--- Calculate width and jacobian
	//	const double width = calcWidth( iElem );
	//	const double jacobian = 0.50 * width;

	//	//const std::complex<double> val1 = std::complex<double>( width * alpha2 / 3 , 0.0 );
	//	//const std::complex<double> val2 = std::complex<double>( width * alpha2 / 6 , 0.0 );
	//	//const std::complex<double> val3 = std::complex<double>( 0.5 * width * beta2 );

	//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][0], nodesLocal2Global[iElem][0], val1 );// Add to matrix
	//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][1], nodesLocal2Global[iElem][1], val1 );// Add to matrix
	//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][0], nodesLocal2Global[iElem][1], val2 );// Add to matrix

	//	//m_matrix2DAnalysis.addRightHandSideVector( nodesLocal2Global[iElem][0], val3 );// Add to right hand side vector
	//	//m_matrix2DAnalysis.addRightHandSideVector( nodesLocal2Global[iElem][1], val3 );// Add to right hand side vector

	//	//--- Calculate contribution to the cofficent matrix
	//	for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){

	//		const int row = nodesLocal2Global[iElem][ nodesUpperSide[iNode1] ];				
	//		if( row < 0 ){
	//			continue;
	//		}

	//		double integral1(0);
	//		for( ip = 0; ip < numGauss; ++ip ){
	//			integral1 += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * CommonParameters::weights2Point[ip];
	//		}
	//		integral1 *= beta2 * jacobian;
	//		std::complex<double> val1( integral1, 0.0 );
	//		m_matrix2DAnalysis.addRightHandSideVector( row, -val1 );// Add to right hand side vector
	//			
	//		if(m_sourceFieldElectric){// Electric field specified at the top of the model as source

	//			for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){										
	//				const int col = nodesLocal2Global[iElem][ nodesUpperSide[iNode2] ];					 
	//				if( col < 0 ){
	//					continue;
	//				}

	//				if( col >= row ){// Store only upper triangle part
	//					double integral2(0);
	//					for( ip = 0; ip < numGauss; ++ip ){
	//						integral2 += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalTopSide, nodesUpperSide[iNode2] ) * CommonParameters::weights2Point[ip];
	//					}
	//					integral2 *= alpha2 * jacobian;

	//					std::complex<double> val2 = std::complex<double>( integral2, 0.0 );

	//					m_matrix2DAnalysis.addNonZeroValues( row, col, -val2 );// Add to matrix
	//				}
	//			}// iNode2

	//		}

	//	}// iNode1

	//}// iElem
	//----- Don't remove because this might need to be restored in the future <<<<<


	//--- Components due to the boudary condition of the bottom side
	if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){// One dimensional boundary condition
		// [Note] : No components due to the boudary condition of the bottom side for right hand side vector because beta1 is equal to zero.

		//----- Don't remove because this might need to be restored in the future >>>>>
//		for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
//			
//			const int nodesBottomSide[2] = { 0, 1 };
//
//			//--- Calculate sigma
//			const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
//			const int blkID = pMeshDataBrickElement->m_elementID2blockID[elemID];		
//			const double sigma = pResistivityBlock->getConductivityValues(blkID);
//
//			//--- Calculate alpha1 ( Inverse number of impedance of half space )
//			const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
//			const double tmp = sigma / ( CommonParameters::mu * omega );
//			const std::complex<double> alpha1 = std::sqrt( std::complex<double>( 0.0, tmp ) ); // alpha1 = 1/Z = H/E = k/omega/mu
//
//			//----- debug >>>>>
//#ifdef _DEBUG_WRITE
//			const std::complex<double> Z = static_cast< std::complex<double> >(1.0)/alpha1;
//			double rho = static_cast<double>( std::pow(std::abs(Z),2.0) ) / ( CommonParameters::mu * omega );
//			double phase = atan2( imag(Z), real(Z) ) * CommonParameters::rad2deg;
//			std::cout << "Z rho phase " << Z << " " << rho << " " << phase <<  std::endl;
//#endif
//			//----- debug <<<<<
//
//			//--- Calculate width and jacobian
//			const double width = calcWidth( iElem );
//			const double jacobian = 0.50 * width;
//
//			//const std::complex<double> val1 = std::complex<double>( width / 3 , 0.0 ) * alpha1;
//			//const std::complex<double> val2 = std::complex<double>( width / 6 , 0.0 ) * alpha1;
//
//			//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][2], nodesLocal2Global[iElem][2], val1 );// Add to matrix
//			//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][3], nodesLocal2Global[iElem][3], val1 );// Add to matrix
//			//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][2], nodesLocal2Global[iElem][3], val2 );// Add to matrix
//
//			//--- Calculate contribution to the cofficent matrix
//			for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){
//				for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){
//
//					const int row = nodesLocal2Global[iElem][ nodesBottomSide[iNode1] ];
//					const int col = nodesLocal2Global[iElem][ nodesBottomSide[iNode2] ];
//
//					if( row < 0 || col < 0 ){
//						continue;
//					}
//
//					if( col >= row ){// Store only upper triangle part
//						double integral(0);
//						for( ip = 0; ip < numGauss; ++ip ){
//							integral += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights2Point[ip];
//						}
//						integral *= jacobian;
//
//						std::complex<double> val = std::complex<double>( integral, 0.0 ) * alpha1;
//
//						m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
//					}
//
//				}// iNode2
//			}// iNode1
//
//		}// iElem
		//----- Don't remove because this might need to be restored in the future <<<<<

		for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
			
			const int nodesBottomSide[2] = { 0, 1 };

			//--- Calculate sigma
			//const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
			const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
			const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

			//--- Calculate alpha1 ( Inverse number of impedance of half space )
			const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
			const double tmp = sigma / ( CommonParameters::mu * omega );
			const std::complex<double> alpha1 = std::sqrt( std::complex<double>( 0.0, tmp ) ); // alpha1 = 1/Z = H/E = k/omega/mu

			//--- Calculate width and jacobian
			const double width = calcWidth( iElem, pMeshDataBrickElement );
			const double jacobian = 0.50 * width;

			//--- Calculate contribution to the cofficent matrix and right hand side vector
			for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){

				const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];				
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){

					const int col = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode2] ];					 
					if( col <= DEGENERATED_NODES_ON_BOTTOM_SIDE ){
						continue;
					}

					double integral(0);
					for( ip = 0; ip < numGauss; ++ip ){
						integral += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights2Point[ip];
					}
					integral *= jacobian;

					std::complex<double> val = std::complex<double>( integral, 0.0 ) * alpha1;

					if( col == DEGENERATED_NODES_ON_TOP_SIDE ){
						m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueElectric, 0.0) );// Add to right hand side vector
					}else if( col >= row ){// Store only upper triangle part
						m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
					}

				}// iNode2
			}// iNode1

		}// iElem

	}else if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){// Perfect conductor at the bottom
		// [Note] : There is no non-zero components because electric field is forced to zero on the bottom by Dirichlet boundary condition

		//----- Don't remove because this might need to be restored in the future >>>>>
		//for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){

		//	const int nodesBottomSide[2] = { 0, 1 };
		//	const double alpha1 = 1.0e+20;

		//	//--- Calculate width and jacobian
		//	const double width = calcWidth( iElem );
		//	const double jacobian = 0.50 * width;

		//	//const std::complex<double> val1 = std::complex<double>( width * alpha1 / 3 , 0.0 );
		//	//const std::complex<double> val2 = std::complex<double>( width * alpha1 / 6 , 0.0 );

		//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][2], nodesLocal2Global[iElem][2], val1 );// Add to matrix
		//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][3], nodesLocal2Global[iElem][3], val1 );// Add to matrix
		//	//m_matrix2DAnalysis.addNonZeroValues( nodesLocal2Global[iElem][2], nodesLocal2Global[iElem][3], val2 );// Add to matrix

		//	//--- Calculate contribution to the cofficent matrix and right hand side vector
		//	for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){
		//		for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){

		//			const int row = nodesLocal2Global[iElem][ nodesBottomSide[iNode1] ];				
		//			const int col = nodesLocal2Global[iElem][ nodesBottomSide[iNode2] ];					 
		//			if( row < 0 || col < 0 ){
		//				continue;
		//			}

		//			if( col >= row ){// Store only upper triangle part
		//				double integral(0);
		//				for( ip = 0; ip < numGauss; ++ip ){
		//					integral += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights2Point[ip];
		//				}
		//				integral *= alpha1 * jacobian;

		//				std::complex<double> val = std::complex<double>( integral, 0.0 );

		//				m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
		//			}

		//		}// iNode2
		//	}// iNode1
		//				
		//}// iElem
		//----- Don't remove because this might need to be restored in the future <<<<<

		//----- Don't remove because this might need to be restored in the future >>>>>
		//for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){

		//	const int nodesBottomSide[2] = { 0, 1 };
		//	const double alpha1 = 1.0e+20;

		//	//--- Calculate width and jacobian
		//	const double width = calcWidth( iElem );
		//	const double jacobian = 0.50 * width;

		//	//--- Calculate contribution to the cofficent matrix and right hand side vector
		//	for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){

		//		const int row = nodesLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];				
		//		if( row < 0 ){
		//			continue;
		//		}

		//		for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){

		//			double integral(0);
		//			for( ip = 0; ip < numGauss; ++ip ){
		//				integral += getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc1stOrderNodeBased( CommonParameters::abscissas2Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights2Point[ip];
		//			}
		//			integral *= alpha1 * jacobian;
		//			std::complex<double> val = std::complex<double>( integral, 0.0 );

		//			const int col = nodesLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode2] ];					 
		//			if( col < 0 ){
		//				m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueElectric, 0.0) );// Add to right hand side vector
		//			}else if( col >= row ){// Store only upper triangle part
		//				m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
		//			}

		//		}// iNode2
		//	}// iNode1
		//				
		//}// iElem
		//----- Don't remove because this might need to be restored in the future <<<<<

	}else{
		OutputFiles::m_logFile << "Error : m_boundaryConditionBottom is wrong !! m_boundaryConditionBottom = " << pAnalysisControl->getBoundaryConditionBottom() << std::endl;
		exit(1);
	}
	
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

	//--- Solve phase of matrix solver
	std::complex<double>* solutionDegenerated = new std::complex<double>[nEqDegenerated];
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
		for( int iNode = 0; iNode < 4; ++iNode ){

			const int iNum = m_IDsLocal2Global[iElem][iNode];

			if( alreadyFound[ iNum ] == false ){
				const int iNumDegenerated = m_IDsLocal2GlobalDegenerated[iElem][iNode];

				if( iNumDegenerated >= 0 ){
					m_solution[iNum] = solutionDegenerated[ iNumDegenerated ]; 
				}

				alreadyFound[ iNum ] = true;
			}
			
		}// iNode
	}// iElem
	for( int i = 0; i < numElemW + 1; ++i ){
		m_solution[i] = std::complex<double>(sourceValueElectric, 0.0);
	}
	if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){// Perfect conductor at the bottom
		for( int i = nEq - 1; i > nEq - 2 - numElemW; --i ){
			m_solution[i] = std::complex<double>(0.0, 0.0);
		}
	}
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
//	output2DResult( Forward2DSquareElement::NODE_BASED_FIRST_ORDER, iPol, freq, nElem, numElemW );	
//#endif
//	//----- debug <<<<<
	output2DResult( Forward2DSquareElement::NODE_BASED_FIRST_ORDER, freq, nElem, numElemW, pMeshDataBrickElement );	

	//--- Release memory
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete[] nodesLocal2Global[iElem];
	//}
	//delete[] nodesLocal2Global;

	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete[] nodesLocal2GlobalDegenerated[iElem];
	//}
	//delete[] nodesLocal2GlobalDegenerated;

	//delete[] nodesGlobal2FullModel;

}

//// Get electric field perpendicular to bondary plane from element ID and coordinate values
//std::complex<double> Forward2DSquareElement1stOrderNodeBased::getElectricFieldPerpendicularToPlane( const int iElem, const double wLocal, const double hLocal ) const{
//
//	return calcValueV1stOrder( iElem, wLocal, hLocal );
//
//}

// Get shape functions of 1st order node-based element
inline double Forward2DSquareElement1stOrderNodeBased::getShapeFunc1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( 1 + wLocal ) * ( 1 + hLocal );
			break;
		case 1:
			return   0.25 * ( 1 - wLocal ) * ( 1 + hLocal );
			break;
		case 2:
			return   0.25 * ( 1 - wLocal ) * ( 1 - hLocal );
			break;
		case 3:
			return   0.25 * ( 1 + wLocal ) * ( 1 - hLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFunc1stOrderNodeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of 1st order node-based element defferentiated by local coordinate of w
inline double Forward2DSquareElement1stOrderNodeBased::getShapeFuncDiffByWLocal1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( 1 + hLocal );
			break;
		case 1:
			return - 0.25 * ( 1 + hLocal );
			break;
		case 2:
			return - 0.25 * ( 1 - hLocal );
			break;
		case 3:
			return   0.25 * ( 1 - hLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncDiffByWLocal1stOrderNodeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions of 1st order node-based element defferentiated by local coordinate of h
inline double Forward2DSquareElement1stOrderNodeBased::getShapeFuncDiffByHLocal1stOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( 1 + wLocal );
			break;
		case 1:
			return   0.25 * ( 1 - wLocal );
			break;
		case 2:
			return - 0.25 * ( 1 - wLocal );
			break;
		case 3:
			return - 0.25 * ( 1 + wLocal );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncDiffByHLocal1stOrderNodeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

//
//// Calculate parameter V of Rodi(1976) for 1st order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 3) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 4) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter V
//std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueV1stOrder( const int iElem, const double wLocal, const double hLocal ) const{
//
//	if( m_solution == NULL ){
//		OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
//		exit(1);
//	}else if( m_IDsLocal2Global[iElem] == NULL ){
//		OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
//		exit(1);
//	}
//
//	std::complex<double> val(0.0, 0.0);
//	for( int i = 0; i < 4; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFunc1stOrderNodeBased( wLocal, hLocal, i );
//	}
//	return val;
//
//}
//
//// Calculate parameter J of Rodi(1976) for 1st order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) freq   : Frequency
////           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 4) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 5) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter J
//std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueJ1stOrder( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
//
//	if( m_solution == NULL ){
//		OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
//		exit(1);
//	}else if( m_IDsLocal2Global[iElem] == NULL ){
//		OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
//		exit(1);
//	}
//
//	std::complex<double> val(0.0, 0.0);
//	for( int i = 0; i < 4; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByHLocal1stOrderNodeBased( wLocal, hLocal, i );
//	}
//
//	const int imode = calcMode();// TM or TE mode
//	std::complex<double> eta = calcEta( imode, freq, iElem );
//	double height = calcHeight( iElem, pMeshDataBrickElement );
//
//	val *= ( static_cast< std::complex<double> >(-2.0/height) / eta );
//
//	return val;
//}
//
//// Calculate parameter I of Rodi(1976) for 1st order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) freq   : Frequency
////           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 4) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 5) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter I
//std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueI1stOrder( const double freq,  const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
//
//	if( m_solution == NULL ){
//		OutputFiles::m_logFile << "Error : m_solution is NULL !! " << std::endl;
//		exit(1);
//	}else if( m_IDsLocal2Global[iElem] == NULL ){
//		OutputFiles::m_logFile << "Error : m_IDsLocal2Global[iElem] is NULL !! " << std::endl;
//		exit(1);
//	}
//
//	std::complex<double> val(0.0, 0.0);
//	for( int i = 0; i < 4; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByWLocal1stOrderNodeBased( wLocal, hLocal, i );
//	}
//
//	const int imode = calcMode();// TM or TE mode
//	std::complex<double> eta = calcEta( imode, freq, iElem );
//	double width = calcWidth( iElem, pMeshDataBrickElement );
//
//	val *= ( static_cast< std::complex<double> >(-2.0/width) / eta );
//
//	return val;
//}

// Calculate parameter V of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldPerpendicular( iElem, wLocal, hLocal );
}

// Calculate parameter J of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldHorizontal( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate parameter I of Rodi(1976)
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldVertical( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate electric field perpendicular to the boundary plane
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueElectricFieldPerpendicular( const int iElem, const double wCoord, const double hCoord ) const{

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
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFunc1stOrderNodeBased( wCoord, hCoord, i );
	}
	return val;

}

// Calculate horizontal magnetic field
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueMagneticFieldHorizontal( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

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
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByHLocal1stOrderNodeBased( wCoord, hCoord, i );
	}

	const int imode = calcMode();// TM or TE mode
	std::complex<double> eta = calcEta( imode, freq, iElem );
	double height = calcHeight( iElem, pMeshDataBrickElement );

	val *= ( static_cast< std::complex<double> >(-2.0/height) / eta );

	return val;

}

// Calculate vertical magnetic field
std::complex<double> Forward2DSquareElement1stOrderNodeBased::calcValueMagneticFieldVertical( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

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
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByWLocal1stOrderNodeBased( wCoord, hCoord, i );
	}

	const int imode = calcMode();// TM or TE mode
	std::complex<double> eta = calcEta( imode, freq, iElem );
	double width = calcWidth( iElem, pMeshDataBrickElement );

	val *= ( static_cast< std::complex<double> >(-2.0/width) / eta );

	return val;

}