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
#include "Forward2DSquareElement2ndOrderNodeBased.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "ResistivityBlock.h"
#include <iostream>
#include <assert.h>

//// Defailt constructer
//Forward2DSquareElement2ndOrderNodeBased::Forward2DSquareElement2ndOrderNodeBased()
//{}

// Constructer
Forward2DSquareElement2ndOrderNodeBased::Forward2DSquareElement2ndOrderNodeBased( const int planeID, const int iPol ):
	Forward2DSquareElementNodeBased( planeID, iPol )
{}

// Destructer
Forward2DSquareElement2ndOrderNodeBased::~Forward2DSquareElement2ndOrderNodeBased()
{}

// Calculate EM fields of boundary planes by 2D forward calculcation
void Forward2DSquareElement2ndOrderNodeBased::calcEMFieldsOfBoundaryPlanes( const double freq, const MeshDataBrickElement* const pMeshDataBrickElement ){

	const int imode = calcMode();// TM or TE mode
	if( imode != CommonParameters::TM_MODE && imode != CommonParameters::TE_MODE ){
		OutputFiles::m_logFile << "Error : Wrong mode in calcEMFieldsOfBoundaryPlanes2ndOrderNodeBased !! imode = " << imode << "." << std::endl;
		exit(1);
	}
#ifdef _DEBUG_WRITE
	std::cout << "imode " << imode << std::endl;// For debug
#endif

	OutputFiles::m_logFile << "# Calculate electric field on a boundary plane with 2nd order node-based element." << std::endl;

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
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
	const int nEq = ( 2 * numElemW + 1 ) * ( 2 * numElemH + 1 ) - numElemW * numElemH;
	int tmpNEq = nEq;
	if( ( imode == CommonParameters::TM_MODE && m_sourceFieldElectric == false ) ||
		( imode == CommonParameters::TE_MODE && m_sourceFieldElectric == true  ) ) {
		// TM mode and magnetic field specified at the top of the model as source
		// TE mode and electric field specified at the top of the model as source
		tmpNEq -= ( 2 * numElemW + 1 );
	}
	if( imode == CommonParameters::TE_MODE && pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		// TE mode and perfect conductor is specified at the bottom
		// => Exclude nodes of bottom side on which Dirichlet condition is given
		tmpNEq -= ( 2 * numElemW + 1 );
	}
	const int nEqDegenerated = tmpNEq;

	OutputFiles::m_logFile << "# Number of equation = " << nEq
		<< ", Number of equation after degeneration = " << nEqDegenerated << std::endl;

	if( m_matrix2DAnalysis.getDegreeOfEquation() == 0 ){
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

		const int offset = 3 * numElemW + 2;// <= ( 2 * numElemW + 1 ) + ( numElemW + 1 )

		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int ih = iElem / numElemW;
			const int iw = iElem - ih * numElemW;

			m_IDsLocal2Global[iElem] = new int[8];

			m_IDsLocal2Global[iElem][0] = ( ih + 1 ) * offset + 2 * iw + 2;
			m_IDsLocal2Global[iElem][1] = ( ih + 1 ) * offset + 2 * iw;
			m_IDsLocal2Global[iElem][2] =   ih       * offset + 2 * iw;
			m_IDsLocal2Global[iElem][3] =   ih       * offset + 2 * iw + 2;
			m_IDsLocal2Global[iElem][4] = ( ih + 1 ) * offset + 2 * iw + 1;
			m_IDsLocal2Global[iElem][5] = m_IDsLocal2Global[iElem][2] + 2 * ( numElemW - iw ) + ( iw + 1 );
			m_IDsLocal2Global[iElem][6] =   ih       * offset + 2 * iw + 1;
			m_IDsLocal2Global[iElem][7] = m_IDsLocal2Global[iElem][5] + 1;
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

		if( ( imode == CommonParameters::TM_MODE && m_sourceFieldElectric == false ) ||
			( imode == CommonParameters::TE_MODE && m_sourceFieldElectric == true  ) ) {
			// TM mode and magnetic field specified at the top of the model as source
			// TE mode and electric field specified at the top of the model as source
	
			for( int iElem = 0; iElem < nElem; ++iElem ){
				m_IDsLocal2GlobalDegenerated[iElem] = new int[8];
				for( int i = 0; i < 8; ++i ){
					m_IDsLocal2GlobalDegenerated[iElem][i] = m_IDsLocal2Global[iElem][i] - ( 2 * numElemW + 1 );
				}

				if( iElem < numElemW ){
					m_IDsLocal2GlobalDegenerated[iElem][2] = DEGENERATED_NODES_ON_TOP_SIDE;
					m_IDsLocal2GlobalDegenerated[iElem][3] = DEGENERATED_NODES_ON_TOP_SIDE;
					m_IDsLocal2GlobalDegenerated[iElem][6] = DEGENERATED_NODES_ON_TOP_SIDE;
				}		
			}

		}else{

			for( int iElem = 0; iElem < nElem; ++iElem ){
				m_IDsLocal2GlobalDegenerated[iElem] = new int[8];
				for( int i = 0; i < 8; ++i ){
					m_IDsLocal2GlobalDegenerated[iElem][i] = m_IDsLocal2Global[iElem][i];
				}
			}

		}
		if( imode == CommonParameters::TE_MODE && pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
			// TE mode and perfect conductor is specified at the bottom
			// => Exclude nodes of bottom side on which Dirichlet condition is given
			for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
				m_IDsLocal2GlobalDegenerated[iElem][0] = DEGENERATED_NODES_ON_BOTTOM_SIDE;
				m_IDsLocal2GlobalDegenerated[iElem][1] = DEGENERATED_NODES_ON_BOTTOM_SIDE;
				m_IDsLocal2GlobalDegenerated[iElem][4] = DEGENERATED_NODES_ON_BOTTOM_SIDE;
			}
		}

	//----- debug >>>>>
#ifdef _DEBUG_WRITE
		for( int iElem = 0; iElem < nElem; ++iElem ){
			std::cout << "iElem m_IDsLocal2GlobalDegenerated : " << iElem;
			for( int i = 0; i < 8; ++i ){
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
	//	for( int i = 4; i < 8; ++i ){
	//		nodesGlobal2FullModel[ nodesLocal2Global[iElem][i] ] = -1;
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

			for( int iNode1 = 0; iNode1 < 8; ++iNode1 ){

				const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 8; ++iNode2 ){

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
		if( imode == CommonParameters::TM_MODE ){ // TM mode

			//--- Components due to the boudary condition of the top side
			if(m_sourceFieldElectric){// Electric field specified at the top of the model as source
				// [Note] : No matrix components due to the boudary condition of the top side because alpha2 is equal to zero.
			}else{// Magnetic field specified at the top of the model as source
				// [Note] : No matrix components due to the Dirichlet boudary condition is specified.

				//----- Don't remove because this might need to be restored in the future >>>>>
				////--- Components due to the boundary condition of the top side
				//int nodesUpperSide[3] = { 2, 3, 6 };
				//for( int iElem = 0; iElem < numElemW; ++iElem ){
				//	for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){
				//		for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

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
				//----- Don't remove because this might need to be restored in the future <<<<<
			}

			//--- Components due to the boudary condition of the bottom side
			if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){
				// One dimensional boundary condition

				int nodesBottomSide[3] = { 0, 1, 4 };
				for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){

					for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

						const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];
						if( row < 0 ){
							continue;
						}

						for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

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
			else if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){// Perfect conductor at the bottom
				// [Note] : No matrix components due to the boudary condition of the bottom side because alpha1 is equal to zero.
			}
			else{
				OutputFiles::m_logFile << "Error : m_boundaryConditionBottom is wrong !! m_boundaryConditionBottom = " << pAnalysisControl->getBoundaryConditionBottom() << std::endl;
				exit(1);
			}				

		}else if ( imode == CommonParameters::TE_MODE ){// TE mode

			//--- Components due to the boudary condition of the top side
			if(m_sourceFieldElectric){// Electric field specified at the top of the model as source
				// [Note] : No matrix components due to the Dirichlet boudary condition is specified.

				////----- Don't remove because this might need to be restored in the future >>>>>
				////--- Components due to the boundary condition of the top side
				//int nodesUpperSide[3] = { 2, 3, 6 };
				//for( int iElem = 0; iElem < numElemW; ++iElem ){
				//	for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){
				//		for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

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
				////----- Don't remove because this might need to be restored in the future <<<<<
			}else{// Magnetic field specified at the top of the model as source
				// [Note] : No matrix components due to the boudary condition of the top side because alpha2 is equal to zero.
			}

			//--- Components due to the boudary condition of the bottom side
			if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){// One dimensional boundary condition
				int nodesBottomSide[3] = { 0, 1, 4 };
				for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){

					for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

						const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];
						if( row < 0 ){
							continue;
						}

						for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

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

		}
		else{
			OutputFiles::m_logFile << "Error : Wrong mode !! imode = " << imode << "." << std::endl;
			exit(1);
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

	//-----------------------------------------------------------------
	//--- Calculate integral points and weights of Gauss quadrature ---
	//-----------------------------------------------------------------
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
	m_matrix2DAnalysis.zeroClearNonZeroValues();// Zero clear matrix values
	m_matrix2DAnalysis.zeroClearRightHandSideVector();// Zero clear right hand side vector

	OutputFiles::m_logFile << "# Start setting values of matrix and right hand side vector." << pAnalysisControl->outputElapsedTime() << std::endl;
	ResistivityBlock* pResistivityBlock = ResistivityBlock::getInstance();

	const double sourceValueMagnetic = 1;
	const double sourceValueElectric = CommonParameters::sourceValueElectric;

	//--- Components due to the boudary conditions
	if( imode == CommonParameters::TM_MODE ){// TM mode

		//------------------------------------------
		//--- Components due to stiffness matrix ---
		//------------------------------------------
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

			//----- Don't remove because this might need to be restored in the future >>>>>
			//for( int iNode1 = 0; iNode1 < 8; ++iNode1 ){
			//	for( int iNode2 = 0; iNode2 < 8; ++iNode2 ){

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
			//				integral1 += getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//				integral2 += getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//				integral3 += getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//			}
			//			integral1 *= jacobian;
			//			integral2 *= constWLocal;
			//			integral3 *= constHLocal;

			//			std::complex<double> val = std::complex<double>( integral1             , 0.0 ) * gamma  
			//										+ std::complex<double>( integral2 + integral3 , 0.0 ) / eta;

			//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
			//		}
			//			
			//	}// iNode2
			//}// iNode1
			//----- Don't remove because this might need to be restored in the future <<<<<

			for( int iNode1 = 0; iNode1 < 8; ++iNode1 ){

				const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 8; ++iNode2 ){

					const int col = m_IDsLocal2GlobalDegenerated[iElem][iNode2];
					if(  m_sourceFieldElectric == true ){
						// Electric field specified at the top of the model as source
						// => Dirichlet boudary condition is NOT specified.
						if( col < 0 ){
							continue;
						}
					}

					double integral1(0);
					double integral2(0);
					double integral3(0);
					for( ip = 0; ip < numIntegralPoints; ++ip ){
						integral1 += getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
						integral2 += getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
						integral3 += getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
					}
					integral1 *= jacobian;
					integral2 *= constWLocal;
					integral3 *= constHLocal;

					std::complex<double> val = std::complex<double>( integral1             , 0.0 ) * gamma  
											 + std::complex<double>( integral2 + integral3 , 0.0 ) / eta;

					if( col < 0 ){
						m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueMagnetic, 0.0) );// Add to right hand side vector
					}else if( col >= row ){// Store only upper triangle part
						m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
					}
					
				}// iNode2
			}// iNode1
		
		}// iElem

		//----- Don't remove because this might need to be restored in the future >>>>>
		////--- Components due to the boudary condition of the top side
		//for( int iElem = 0; iElem < numElemW; ++iElem ){
		//	
		//	const int nodesUpperSide[3] = { 2, 3, 6 };

		//	double alpha2(0);
		//	double beta2(0);
		//	if(m_sourceFieldElectric){// Electric field specified at the top of the model as source
		//		beta2 = -1.0 * sourceValueElectric;
		//	}else{// Magnetic field specified at the top of the model as source
		//		alpha2 = 1.0e+20;
		//		beta2 = 1.0e+20 * sourceValueMagnetic;
		//	}

		//	//--- Calculate width and jacobian
		//	const double width = calcWidth( iElem );
		//	const double jacobian = 0.50 * width;

		//	//--- Calculate contribution to the cofficent matrix
		//	for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

		//		const int row = nodesLocal2Global[iElem][ nodesUpperSide[iNode1] ];				
		//		if( row < 0 ){
		//			continue;
		//		}

		//		double integral1(0);
		//		for( ip = 0; ip < numGauss; ++ip ){
		//			integral1 += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * CommonParameters::weights3Point[ip];
		//		}
		//		integral1 *= beta2 * jacobian;
		//		std::complex<double> val1( integral1, 0.0 );
		//		m_matrix2DAnalysis.addRightHandSideVector( row, -val1 );// Add to right hand side vector
		//		
		//		if(m_sourceFieldElectric == false){// Magnetic field specified at the top of the model as source

		//			for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){										
		//				const int col = nodesLocal2Global[iElem][ nodesUpperSide[iNode2] ];					 
		//				if( col < 0 ){
		//					continue;
		//				}

		//				if( col >= row ){// Store only upper triangle part
		//					double integral2(0);
		//					for( ip = 0; ip < numGauss; ++ip ){
		//						integral2 += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalTopSide, nodesUpperSide[iNode2] ) * CommonParameters::weights3Point[ip];
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

		//---------------------------------------------------------------
		//--- Components due to the boudary condition of the top side ---
		//---------------------------------------------------------------
		if( m_sourceFieldElectric == false ){
			// Magnetic field specified at the top of the model as source
			// => Dirichlet boudary condition is specified
		}else{
			// Electric field specified at the top of the model as source
			for( int iElem = 0; iElem < numElemW; ++iElem ){
			
				const int nodesUpperSide[3] = { 2, 3, 6 };

				double beta2 = -1.0 * sourceValueElectric;

				//--- Calculate width and jacobian
				const double width = calcWidth( iElem, pMeshDataBrickElement );
				const double jacobian = 0.50 * width;

				//--- Calculate contribution to the cofficent matrix
				for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

					//const int row = nodesLocal2Global[iElem][ nodesUpperSide[iNode1] ];				
					const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesUpperSide[iNode1] ];				
					if( row < 0 ){
						continue;
					}

					double integral1(0);
					for( ip = 0; ip < numGauss; ++ip ){
						integral1 += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * CommonParameters::weights3Point[ip];
					}
					integral1 *= beta2 * jacobian;
					std::complex<double> val1( integral1, 0.0 );

					m_matrix2DAnalysis.addRightHandSideVector( row, -val1 );// Add to right hand side vector

				}// iNode1

			}// iElem
		}

		//------------------------------------------------------------------
		//--- Components due to the boudary condition of the bottom side ---
		//------------------------------------------------------------------
		// [Note] : No components due to the boudary condition of the bottom side for right hand side vector because beta1 is equal to zero.

		if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){// One dimensional boundary condition

			for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){
				const int nodesBottomSide[3] = { 0, 1, 4 };

				//--- Calculate sigma
				//const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
				const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
				const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

				//--- Calculate alpha1 ( Inverse number of impedance of half space )
				const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
				const double tmp = ( CommonParameters::mu * omega ) / sigma;
				const std::complex<double> alpha1 = std::sqrt( std::complex<double>( 0.0, -tmp ) ); // alpha1 = Z = E/H = omega*mu/k

				//----- debug >>>>>
#ifdef _DEBUG_WRITE
				const std::complex<double> Z = alpha1;
				double rho = static_cast<double>( std::pow(std::abs(Z),2.0) ) / ( CommonParameters::mu * omega );
				double phase = atan2( imag(Z), real(Z) ) * CommonParameters::rad2deg;
				std::cout << "Z rho phase " << Z << " " << rho << " " << phase <<  std::endl;
#endif
				//----- debug <<<<<

				//--- Calculate Jacobian
				const double width = calcWidth( iElem, pMeshDataBrickElement );
				const double jacobian = 0.50 * width;

				//----- Don't remove because this might need to be restored in the future >>>>>
				//for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){
				//	for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){
				//		const int row = nodesLocal2Global[iElem][ nodesBottomSide[iNode1] ];
				//		const int col = nodesLocal2Global[iElem][ nodesBottomSide[iNode2] ];

				//		if( row < 0 || col < 0 ){
				//			continue;
				//		}

				//		if( col >= row ){// Store only upper triangle part
				//			double integral(0);
				//			for( ip = 0; ip < numGauss; ++ip ){
				//				integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
				//			}
				//			integral *= jacobian;

				//			std::complex<double> val = std::complex<double>( integral, 0.0 ) * alpha1;

				//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				//		}

				//	}// iNode2
				//}// iNode1
				//----- Don't remove because this might need to be restored in the future <<<<<

				for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

					const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];
					if( row < 0 ){
						continue;
					}

					for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

						const int col = m_IDsLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode2] ];
						//if( m_sourceFieldElectric == true ){
						//	// Electric field specified at the top of the model as source
						//	// => Dirichlet boudary condition is NOT specified.
						//	if( col < 0 ){
						//		continue;
						//	}
						//}
						if( col <= DEGENERATED_NODES_ON_BOTTOM_SIDE ){
							continue;
						}
						
						double integral(0);
						for( ip = 0; ip < numGauss; ++ip ){
							integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
						}
						integral *= jacobian;

						std::complex<double> val = std::complex<double>( integral, 0.0 ) * alpha1;

						if( col == DEGENERATED_NODES_ON_TOP_SIDE ){
							m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueMagnetic, 0.0) );// Add to right hand side vector
						}else if( col >= row ){// Store only upper triangle part
							m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
						}

					}// iNode2
				}// iNode1

			}// iElem

		}else if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){// Perfect conductor at the bottom
			// [Note] : No matrix components due to the boudary condition of the bottom side because alpha1 is equal to zero.

		}else{
			OutputFiles::m_logFile << "Error : m_boundaryConditionBottom is wrong !! m_boundaryConditionBottom = " << pAnalysisControl->getBoundaryConditionBottom() << std::endl;
			exit(1);
		}

	}else if ( imode == CommonParameters::TE_MODE ){// TE mode
		
		//------------------------------------------
		//--- Components due to stiffness matrix ---
		//------------------------------------------
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

			//----- Don't remove because this might need to be restored in the future >>>>>
			//for( int iNode1 = 0; iNode1 < 8; ++iNode1 ){
			//	for( int iNode2 = 0; iNode2 < 8; ++iNode2 ){

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
			//				integral1 += getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//				integral2 += getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//				integral3 += getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
			//			}
			//			integral1 *= jacobian;
			//			integral2 *= constWLocal;
			//			integral3 *= constHLocal;

			//			std::complex<double> val = std::complex<double>( integral1             , 0.0 ) * gamma  
			//										+ std::complex<double>( integral2 + integral3 , 0.0 ) / eta;

			//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
			//		}
			//			
			//	}// iNode2
			//}// iNode1
			//----- Don't remove because this might need to be restored in the future <<<<<

			for( int iNode1 = 0; iNode1 < 8; ++iNode1 ){

				const int row = m_IDsLocal2GlobalDegenerated[iElem][iNode1];
				if( row < 0 ){
					continue;
				}

				for( int iNode2 = 0; iNode2 < 8; ++iNode2 ){

					const int col = m_IDsLocal2GlobalDegenerated[iElem][iNode2];
					if( col <= DEGENERATED_NODES_ON_BOTTOM_SIDE ){
						continue;
					}
					//if( m_sourceFieldElectric == false && col == DEGENERATED_NODES_ON_TOP_SIDE ){
					//	// Magnetic field specified at the top of the model as source => Dirichlet boudary condition is NOT specified.
					//	continue;
					//}

					double integral1(0);
					double integral2(0);
					double integral3(0);
					for( ip = 0; ip < numIntegralPoints; ++ip ){
						integral1 += getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFunc2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
						integral2 += getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
						integral3 += getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode1 ) * getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal[ip], hLocal[ip], iNode2 ) * weights2D[ip];
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

		//---------------------------------------------------------------
		//--- Components due to the boudary condition of the top side ---
		//---------------------------------------------------------------
		// [Note] : No matrix components due to the boudary condition of the top side because alpha2 is equal to zero.

		if( m_sourceFieldElectric == true ){
			// Electric field specified at the top of the model as source
			// => Dirichlet boudary condition is specified
		}else{
			// Magnetic field specified at the top of the model as source

			for( int iElem = 0; iElem < numElemW; ++iElem ){
			
				const int nodesUpperSide[3] = { 2, 3, 6 };

				double beta2 = -1.0 * sourceValueMagnetic;

				//--- Calculate width and jacobian
				const double width = calcWidth( iElem, pMeshDataBrickElement );
				const double jacobian = 0.50 * width;

				for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){

					//const int row = nodesLocal2Global[iElem][ nodesUpperSide[iNode1] ];
					const int row = m_IDsLocal2GlobalDegenerated[iElem][ nodesUpperSide[iNode1] ];
					if( row < 0 ){
						continue;
					}

					double integral(0);
					for( ip = 0; ip < numGauss; ++ip ){
						integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalTopSide, nodesUpperSide[iNode1] ) * CommonParameters::weights3Point[ip];
					}
					integral *= beta2 * jacobian;
					std::complex<double> val( integral, 0.0 );

					m_matrix2DAnalysis.addRightHandSideVector( row, -val );// Add to matrix

				}// iNode1

			}// iElem

		}

		//------------------------------------------------------------------
		//--- Components due to the boudary condition of the bottom side ---
		//------------------------------------------------------------------
		if( pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_ONE_DIMENSIONAL ){// One dimensional boundary condition
			// [Note] : No components due to the boudary condition of the bottom side for right hand side vector because beta1 is equal to zero.

			for( int iElem = nElem - 1; iElem > nElem - 1 - numElemW; --iElem ){

				const int nodesBottomSide[3] = { 0, 1, 4 };

				//--- Calculate sigma
				//const int elemID = pMeshDataBrickElement->m_elemBoundaryPlanes[m_planeID][iElem];
				const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
				const double sigma = pResistivityBlock->getConductivityValuesFromElemID(elemID);

				//--- Calculate alpha1 ( Inverse number of impedance of half space )
				const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
				const double tmp = sigma / ( CommonParameters::mu * omega );
				const std::complex<double> alpha1 = std::sqrt( std::complex<double>( 0.0, tmp ) ); // alpha1 = 1/Z = H/E = k/omega/mu

				//----- debug >>>>>
#ifdef _DEBUG_WRITE
				const std::complex<double> Z = static_cast< std::complex<double> >(1.0)/alpha1;
				double rho = static_cast<double>( std::pow(std::abs(Z),2.0) ) / ( CommonParameters::mu * omega );
				double phase = atan2( imag(Z), real(Z) ) * CommonParameters::rad2deg;
				std::cout << "Z rho phase " << Z << " " << rho << " " << phase <<  std::endl;
#endif
				//----- debug <<<<<

				//--- Calculate Jacobian
				const double width = calcWidth( iElem, pMeshDataBrickElement );
				const double jacobian = 0.50 * width;

				//----- Don't remove because this might need to be restored in the future >>>>>
				//for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){
				//	for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

				//		const int row = nodesLocal2Global[iElem][ nodesBottomSide[iNode1] ];
				//		const int col = nodesLocal2Global[iElem][ nodesBottomSide[iNode2] ];					 
				//		if( row < 0 || col < 0 ){
				//			continue;
				//		}

				//		if( col >= row ){// Store only upper triangle part
				//			double integral(0);
				//			for( ip = 0; ip < numGauss; ++ip ){
				//				integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
				//			}
				//			integral *= jacobian;

				//			std::complex<double> val = std::complex<double>( integral, 0.0 ) * alpha1;

				//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
				//		}

				//	}// iNode2
				//}// iNode1
				//----- Don't remove because this might need to be restored in the future <<<<<

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
						//if( m_sourceFieldElectric == false && col == DEGENERATED_NODES_ON_TOP_SIDE ){
						//	// Magnetic field specified at the top of the model as source => Dirichlet boudary condition is NOT specified.
						//	continue;
						//}

						double integral(0);
						for( ip = 0; ip < numGauss; ++ip ){
							integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
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

			//	const int nodesBottomSide[3] = { 0, 1, 4 };
			//	const double alpha1 = 1.0e+20;

			//	//--- Calculate Jacobian
			//	const double width = calcWidth( iElem );
			//	const double jacobian = 0.50 * width;
			//	
			//	//----- Don't remove because this might need to be restored in the future >>>>>
			//	////--- Calculate contribution to the cofficent matrix and right hand side vector
			//	//for( int iNode1 = 0; iNode1 < 3; ++iNode1 ){
			//	//	for( int iNode2 = 0; iNode2 < 3; ++iNode2 ){

			//	//		const int row = nodesLocal2Global[iElem][ nodesBottomSide[iNode1] ];				
			//	//		const int col = nodesLocal2Global[iElem][ nodesBottomSide[iNode2] ];					 
			//	//		if( row < 0 || col < 0 ){
			//	//			continue;
			//	//		}

			//	//		if( col >= row ){// Store only upper triangle part
			//	//			double integral(0);
			//	//			for( ip = 0; ip < numGauss; ++ip ){
			//	//				integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
			//	//			}
			//	//			integral *= alpha1 * jacobian;

			//	//			std::complex<double> val = std::complex<double>( integral, 0.0 );

			//	//			m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
			//	//		}

			//	//	}// iNode2
			//	//}// iNode1
			//	//----- Don't remove because this might need to be restored in the future <<<<<

			//	//--- Calculate contribution to the cofficent matrix and right hand side vector
			//	for( int iNode1 = 0; iNode1 < 2; ++iNode1 ){

			//		const int row = nodesLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode1] ];				
			//		if( row < 0 ){
			//			continue;
			//		}

			//		for( int iNode2 = 0; iNode2 < 2; ++iNode2 ){

			//			const int col = nodesLocal2GlobalDegenerated[iElem][ nodesBottomSide[iNode2] ];
			//			if( m_sourceFieldElectric == false ){
			//				// Magnetic field specified at the top of the model as source
			//				// => Dirichlet boudary condition is NOT specified.
			//				if( col < 0 ){
			//					continue;
			//				}
			//			}

			//			double integral(0);
			//			for( ip = 0; ip < numGauss; ++ip ){
			//				integral += getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode1] ) * getShapeFunc2ndOrderNodeBased( CommonParameters::abscissas3Point[ip], m_hLocalBottomSide, nodesBottomSide[iNode2] ) * CommonParameters::weights3Point[ip];
			//			}
			//			integral *= alpha1 * jacobian;

			//			std::complex<double> val = std::complex<double>( integral, 0.0 );

			//			if( col < 0 ){
			//				m_matrix2DAnalysis.addRightHandSideVector( row, -val * std::complex<double>(sourceValueElectric, 0.0) );// Add to right hand side vector
			//			}else if( col >= row ){// Store only upper triangle part
			//				m_matrix2DAnalysis.addNonZeroValues( row, col, val );// Add to matrix
			//			}

			//		}// iNode2
			//	}// iNode1

			//}// iElem
			//----- Don't remove because this might need to be restored in the future <<<<<

		}else{
			OutputFiles::m_logFile << "Error : m_boundaryConditionBottom is wrong !! m_boundaryConditionBottom = " << pAnalysisControl->getBoundaryConditionBottom() << std::endl;
			exit(1);
		}
		
	}else{
		OutputFiles::m_logFile << "Error : Wrong mode !! imode = " << imode << "." << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	m_matrix2DAnalysis.debugWriteMatrix();
	m_matrix2DAnalysis.debugWriteRightHandSide();
#endif

	//--- Numrical factorization phase of matrix solver
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
		for( int iNode = 0; iNode < 8; ++iNode ){

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

	if( imode == CommonParameters::TM_MODE && m_sourceFieldElectric == false ){
		// TM mode and magnetic field specified at the top of the model as source
		// => Dirichlet boudary condition is specified.
		for( int i = 0; i < 2 * numElemW + 1; ++i ){
			m_solution[i] = std::complex<double>(sourceValueMagnetic, 0.0);
		}
	}else if( imode == CommonParameters::TE_MODE && m_sourceFieldElectric == true  ){
		// TE mode and electric field specified at the top of the model as source
		// => Dirichlet boudary condition is specified.
		for( int i = 0; i < 2 * numElemW + 1; ++i ){
			m_solution[i] = std::complex<double>(sourceValueElectric, 0.0);
		}
	}
	if( imode == CommonParameters::TE_MODE && pAnalysisControl->getBoundaryConditionBottom() == AnalysisControl::BOUNDARY_BOTTOM_PERFECT_CONDUCTOR ){
		// TE mode and perfect conductor specified at the bottom
		for( int i = nEq - 1; i > nEq - 1 - ( 2 * numElemW + 1 ); --i ){
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
//	output2DResult( Forward2DSquareElement::NODE_BASED_SECOND_ORDER, iPol, freq, nElem, numElemW );	
//#endif
//	//----- debug <<<<<
	output2DResult( Forward2DSquareElement::NODE_BASED_SECOND_ORDER, freq, nElem, numElemW, pMeshDataBrickElement );	

	//--- Release memory
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete[] m_IDsLocal2Global[iElem];
	//}
	//delete[] m_IDsLocal2Global;
	//for( int iElem = 0; iElem < nElem; ++iElem ){
	//	delete[] nodesLocal2GlobalDegenerated[iElem];
	//}
	//delete[] nodesLocal2GlobalDegenerated;
	//delete[] nodesGlobal2FullModel;

}

//// Get electric field perpendicular to bondary plane from element ID and coordinate values
//std::complex<double> Forward2DSquareElement2ndOrderNodeBased::getElectricFieldPerpendicularToPlane( const int iElem, const double wLocal, const double hLocal ) const{
//
//	return calcValueV2ndOrder( iElem, wLocal, hLocal );
//
//}

// Get shape functions
inline double Forward2DSquareElement2ndOrderNodeBased::getShapeFunc2ndOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( wLocal + 1 ) * ( hLocal + 1 ) * ( wLocal + hLocal - 1 );
			break;
		case 1:
			return   0.25 * ( wLocal - 1 ) * ( hLocal + 1 ) * ( wLocal - hLocal + 1 );
			break;
		case 2:
			return - 0.25 * ( wLocal - 1 ) * ( hLocal - 1 ) * ( wLocal + hLocal + 1 );
			break;
		case 3:
			return - 0.25 * ( wLocal + 1 ) * ( hLocal - 1 ) * ( wLocal - hLocal - 1 );
			break;
		case 4:
			return - 0.50 * ( wLocal + 1 ) * ( wLocal  - 1 ) * ( hLocal + 1 );
			break;
		case 5:
			return   0.50 * ( wLocal - 1 ) * ( hLocal + 1 ) * ( hLocal - 1 );
			break;
		case 6:
			return   0.50 * ( wLocal + 1 ) * ( wLocal  - 1 ) * ( hLocal - 1 );
			break;
		case 7:
			return - 0.50 * ( wLocal + 1 ) * ( hLocal + 1 ) * ( hLocal - 1 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFunc2ndOrderNodeBased : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions defferentiated by local coordinate of w
inline double Forward2DSquareElement2ndOrderNodeBased::getShapeFuncDiffByWLocal2ndOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( 2 * wLocal + hLocal ) * ( hLocal + 1 );
			break;
		case 1:
			return   0.25 * ( 2 * wLocal - hLocal ) * ( hLocal + 1 );
			break;
		case 2:
			return - 0.25 * ( 2 * wLocal + hLocal ) * ( hLocal - 1 );
			break;
		case 3:
			return - 0.25 * ( 2 * wLocal - hLocal ) * ( hLocal - 1 );
			break;
		case 4:
			return - wLocal * ( hLocal + 1 );
			break;
		case 5:
			return   0.50 * ( hLocal + 1 ) * ( hLocal - 1 );
			break;
		case 6:
			return   wLocal * ( hLocal - 1 );
			break;
		case 7:
			return - 0.50 * ( hLocal + 1 ) * ( hLocal - 1 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncDifferentiatedByWLocal : num = " << num << std::endl;
			exit(1);
			break;
	}

}

// Get shape functions defferentiated by local coordinate of h
inline double Forward2DSquareElement2ndOrderNodeBased::getShapeFuncDiffByHLocal2ndOrderNodeBased( const double wLocal, const double hLocal, const int num ) const{

	switch( num ){
		case 0:
			return   0.25 * ( wLocal + 1 ) * ( wLocal + 2 * hLocal );
			break;
		case 1:
			return   0.25 * ( wLocal - 1 ) * ( wLocal - 2 * hLocal );
			break;
		case 2:
			return - 0.25 * ( wLocal - 1 ) * ( wLocal + 2 * hLocal );
			break;
		case 3:
			return - 0.25 * ( wLocal + 1 ) * ( wLocal - 2 * hLocal );
			break;
		case 4:
			return - 0.50 * ( wLocal + 1 ) * ( wLocal - 1 );
			break;
		case 5:
			return    hLocal * ( wLocal - 1 );
			break;
		case 6:
			return   0.50 * ( wLocal + 1 ) * ( wLocal - 1 );
			break;
		case 7:
			return  - hLocal * ( wLocal + 1 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong number in getShapeFuncDifferentiatedByHLocal : num = " << num << std::endl;
			exit(1);
			break;
	}

}

//
//// Calculate parameter V of Rodi(1976) for 2nd order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 2) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 3) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter V
//std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueV2ndOrder( const int iElem, const double wLocal, const double hLocal ) const{
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
//	for( int i = 0; i < 8; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFunc2ndOrderNodeBased( wLocal, hLocal, i );
//	}
//	return val;
//
//}
//
//// Calculate parameter J of Rodi(1976) for 2nd order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) freq   : Frequency
////           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 4) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 5) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter J
//std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueJ2ndOrder( const double freq,	const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
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
//	for( int i = 0; i < 8; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByHLocal2ndOrderNodeBased( wLocal, hLocal, i );
//	}
//
//	const int imode = calcMode();// TM or TE mode
//	std::complex<double> eta = calcEta( imode, freq, iElem );
//	double height = calcHeight( iElem, pMeshDataBrickElement );
//
//	val *= ( static_cast< std::complex<double> >(-2.0/height) / eta );
//
//	return val;
//
//}
//
//// Calculate parameter I of Rodi(1976) for 2nd order node-based elements
//// [Input] : 1) iPol   : ID of polarization
////           2) freq   : Frequency
////           3) iElem  : Element ID of the boundary plane, that is a numerical sequence beginning with zero
////	 	 	 4) wLocal : Local coordinate of horizontal direction of the element, which corresponds to xi of Rodi(1976)
////			 5) hLocal : Local coordinate of vertical direction of the element, which corresponds to eta of Rodi(1976)
//// [Output] : Vector of parameter I
//std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueI2ndOrder( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
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
//	for( int i = 0; i < 8; ++i ){
//		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByWLocal2ndOrderNodeBased( wLocal, hLocal, i );
//	}
//
//	const int imode = calcMode();// TM or TE mode
//	std::complex<double> eta = calcEta( imode, freq, iElem );
//	double width = calcWidth( iElem, pMeshDataBrickElement );
//
//	val *= ( static_cast< std::complex<double> >(-2.0/width) / eta );
//
//	return val;
//
//}

// Calculate parameter V of Rodi(1976)
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueV( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueElectricFieldPerpendicular( iElem, wLocal, hLocal );
}

// Calculate parameter J of Rodi(1976)
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueJ( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldHorizontal( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate parameter I of Rodi(1976)
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueI( const double freq, const int iElem, const double wLocal, const double hLocal, const MeshDataBrickElement* const pMeshDataBrickElement ) const{
	return calcValueMagneticFieldVertical( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
}

// Calculate electric field perpendicular to the boundary plane
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueElectricFieldPerpendicular( const int iElem, const double wCoord, const double hCoord ) const{

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
	for( int i = 0; i < 8; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFunc2ndOrderNodeBased( wCoord, hCoord, i );
	}
	return val;

}

// Calculate horizontal magnetic field
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueMagneticFieldHorizontal( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

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
	for( int i = 0; i < 8; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByHLocal2ndOrderNodeBased( wCoord, hCoord, i );
	}

	const int imode = calcMode();// TM or TE mode
	std::complex<double> eta = calcEta( imode, freq, iElem );
	double height = calcHeight( iElem, pMeshDataBrickElement );

	val *= ( static_cast< std::complex<double> >(-2.0/height) / eta );

	return val;

}

// Calculate vertical magnetic field
std::complex<double> Forward2DSquareElement2ndOrderNodeBased::calcValueMagneticFieldVertical( const double freq, const int iElem, const double wCoord, const double hCoord, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

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
	for( int i = 0; i < 8; ++i ){
		val += m_solution[ m_IDsLocal2Global[iElem][i] ] * getShapeFuncDiffByWLocal2ndOrderNodeBased( wCoord, hCoord, i );
	}

	const int imode = calcMode();// TM or TE mode
	std::complex<double> eta = calcEta( imode, freq, iElem );
	double width = calcWidth( iElem, pMeshDataBrickElement );

	val *= ( static_cast< std::complex<double> >(-2.0/width) / eta );

	return val;


}