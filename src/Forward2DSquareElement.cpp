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
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <math.h>

#include "CommonParameters.h"
#include "AnalysisControl.h"
#include "MeshDataBrickElement.h"
#include "ResistivityBlock.h"
#include "PARDISOSolver.h"
#include "OutputFiles.h"
#include "Forward2DSquareElement.h"

//// Defailt constructer
//Forward2DSquareElement::Forward2DSquareElement()
//{}

// Constructer
Forward2DSquareElement::Forward2DSquareElement( const int planeID, const int iPol ):
	Forward2D( planeID, iPol )
{}

//Destructer
Forward2DSquareElement::~Forward2DSquareElement()
{}

// Calculate width of element
// [Input] : iElem : Element ID of the boundary plane, that is a numerical sequence beginning with zero
// [Output] : Width of element
double Forward2DSquareElement::calcWidth( const int iElem, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();

	const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
	if(  m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
		return pMeshDataBrickElement->getEdgeLengthY( elemID );
	}else{//XZ Plane
		return pMeshDataBrickElement->getEdgeLengthX( elemID );
	}

	OutputFiles::m_logFile << "Error : Wrong calculation in calcWidth." << std::endl;
	exit(1);

	return 0.0;
}

// Calculate height of element
// [Input] : iElem : Element ID of the boundary plane, that is a numerical sequence beginning with zero
// [Output] : Height of element
double Forward2DSquareElement::calcHeight( const int iElem, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	//	const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();

	//////return pMeshDataBrickElement->getEdgeLengthZ( iElem );

	//const int node1 = pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4    ];
	//const int node3 = pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4 + 2];

	//return std::fabs( pMeshDataBrickElement->m_zCoordinatesOfNodes[node1] - pMeshDataBrickElement->m_zCoordinatesOfNodes[node3] );

	const int elemID = pMeshDataBrickElement->getElemBoundaryPlanes( m_planeID, iElem );
	return pMeshDataBrickElement->getEdgeLengthZ( elemID );

}

// Calculate element division number of horizontal direction
int Forward2DSquareElement::calcNumElemHorizontal( const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	int numElemW = NULL;
	switch( m_planeID ){
		case MeshData::YZMinus:
			// Go through
		case MeshData::YZPlus:
			//numElemW = pMeshDataBrickElement->m_numElemY;
			numElemW = pMeshDataBrickElement->getNumElemY();
			break;
		case MeshData::ZXMinus:
			// Go through
		case MeshData::ZXPlus:
			//numElemW = pMeshDataBrickElement->m_numElemX;
			numElemW = pMeshDataBrickElement->getNumElemX();
			break;
		default:
			OutputFiles::m_logFile << "Error : Unsupported plane in MeshDataBrickElement::calcGlobalDOF " <<  m_planeID << "." << std::endl;
			exit(1);		
			break;
	}

	return numElemW;
}

// Calculate element division number of vertical direction
int Forward2DSquareElement::calcNumElemVertical( const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
	//return pMeshDataBrickElement->m_numElemZ;
	return pMeshDataBrickElement->getNumElemZ();
}

// Output EM field and responses of boundary planes obtained by 2D analysis
void Forward2DSquareElement::output2DResult( const int type, const double freq, const int nElem, const int numElemW, const MeshDataBrickElement* const pMeshDataBrickElement ) const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if( ptrAnalysisControl->getIsOutput2DResult() == false ){
		return;
	}

	//if( !OutputFiles::m_csvFileFor2DFwd.is_open() ){ 
	//	return;
	//}

	const int imode = calcMode();// TM or TE mode

	//----- Input output location
	std::ifstream inFile( "output_2D_result.dat", std::ios::in );
	if( inFile.fail() ){
		OutputFiles::m_logFile << "File open error : output_2D_result.dat !!" << std::endl;
		exit(1);
	}else{

		//----- Read data from output_2D_result.dat
		int ibuf(0);
		//inFile >> ibuf;
		//const int upperOrLower = ibuf;

		//inFile >> ibuf;
		//const int leftOrRight = ibuf;

		const int upperOrLower = 1;

		int leftOrRight(0);
		const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
		const AnalysisControl::UseBackwardOrForwardElement useBackwardOrForwardElement = pAnalysisControl->getUseBackwardOrForwardElement();
		if(  m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
			leftOrRight = useBackwardOrForwardElement.directionY;
		}else{//XZ Plane
			leftOrRight = useBackwardOrForwardElement.directionX;
		}

		inFile >> ibuf;
		const int nCalcPointH = ibuf;

		if( nCalcPointH < 0 ){
			OutputFiles::m_logFile << "Error : Number vertical coordinate for output must be greater than or equals to zero!! nCalcPointH = " << nCalcPointH << std::endl;
			exit(1);
		}

		double* CalcPointH(NULL);
		if( nCalcPointH > 0 ){
			CalcPointH = new double[nCalcPointH];
			for( int i = 0; i < nCalcPointH; ++i ){
				//inFile >> CalcPointH[i];
				double buf(0.0);
				inFile >> buf;
				CalcPointH[i] = buf * 1000.0;
			}
		}

		inFile >> ibuf;
		const int nCalcPointW = ibuf;
		double* CalcPointW(NULL);
		if( nCalcPointW > 0 ){
			CalcPointW = new double[nCalcPointW];
			for( int j = 0; j < nCalcPointW; ++j ){
				//inFile >> CalcPointW[j];
				double buf(0.0);
				inFile >> buf;
				CalcPointW[j] = buf * 1000.0;
			}
		}

		inFile.close();

		//----- Output EM field and response
		if( upperOrLower == 0 ){
			//OutputFiles::m_csvFileFor2DFwd << "Upper,";
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Upper,\n" );
		}else{
			//OutputFiles::m_csvFileFor2DFwd << "Lower,";
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Lower,\n" );
		}

		if( leftOrRight == 0 ){
			//OutputFiles::m_csvFileFor2DFwd << "Left," << std::endl;
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Left,\n" );
		}else{
			//OutputFiles::m_csvFileFor2DFwd << "Right," << std::endl;
			fprintf( OutputFiles::m_csvFileFor2DFwd, "Right,\n" );
		}

		//OutputFiles::m_csvFileFor2DFwd << "Frequency : " << freq << " Polarization : " << iPol << " Plane : "  << m_planeID << std::endl;
		//OutputFiles::m_csvFileFor2DFwd << "Frequency," << freq << "," << std::endl;
		fprintf( OutputFiles::m_csvFileFor2DFwd, "Frequency, %lf\n", freq );
		//OutputFiles::m_csvFileFor2DFwd << "Polarization," << m_polarization << "," << std::endl;
		switch (m_polarization){
			case CommonParameters::EX_POLARIZATION:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Ex-Polarization,\n" );
				break;
			case CommonParameters::EY_POLARIZATION:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Ey-Polarization,\n" );
				break;
			default:
				OutputFiles::m_logFile << "Error : Wrong polarization!! m_polarization = " << m_polarization << std::endl;
				exit(1);
				break;
		}

		//OutputFiles::m_csvFileFor2DFwd << "Plane," << m_planeID << "," << std::endl;
		//fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, %d\n", m_planeID );
		switch (m_planeID){
			case MeshData::YZMinus:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Y-Z(Minus),\n" );
				break;
			case MeshData::YZPlus:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Y-Z(Plus),\n" );
				break;
			case MeshData::ZXMinus:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Z-X(Minus),\n" );
				break;
			case MeshData::ZXPlus:
				fprintf( OutputFiles::m_csvFileFor2DFwd, "Plane, Z-X(Plus),\n" );
				break;
			default:
				OutputFiles::m_logFile << "Error : Wrong plane ID !! m_planeID = " << m_planeID << std::endl;
				exit(1);
				break;
		}

		if(  m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
			//OutputFiles::m_csvFileFor2DFwd << "Y[km],Z[km],real(V),imag(V),real(J),imag(J),real(I),imag(I),real(Z),imag(Z),rho,phase," << std::endl;
			//fprintf( OutputFiles::m_csvFileFor2DFwd, "     Y[km],      Z[km],         real(V),         imag(V),         real(J),         imag(J),         real(I),         imag(I),         real(Z),         imag(Z),             Rho,           Phase,\n" );
			fprintf( OutputFiles::m_csvFileFor2DFwd, "     Y[km],      Z[km],        real(Ex),        imag(Ex),        real(Ey),        imag(Ey),        real(Ez),        imag(Ez),        real(Hx),        imag(Hx),        real(Hy),        imag(Hy),        real(Hz),        imag(Hz),         real(Z),         imag(Z),             Rho,           Phase,\n" );
		}else{//XZ Plane
			//OutputFiles::m_csvFileFor2DFwd << "X[km],Z[km],real(V),imag(V),real(J),imag(J),real(I),imag(I),real(Z),imag(Z),rho,phase," << std::endl;
			//fprintf( OutputFiles::m_csvFileFor2DFwd, "     X[km],      Z[km],         real(V),         imag(V),         real(J),         imag(J),         real(I),         imag(I),         real(Z),         imag(Z),             Rho,           Phase,\n" );
			fprintf( OutputFiles::m_csvFileFor2DFwd, "     X[km],      Z[km],        real(Ex),        imag(Ex),        real(Ey),        imag(Ey),        real(Ez),        imag(Ez),        real(Hx),        imag(Hx),        real(Hy),        imag(Hy),        real(Hz),        imag(Hz),         real(Z),         imag(Z),             Rho,           Phase,\n" );
		}

		int numH = nCalcPointH;

		int numW(0);
		if( nCalcPointW == -1 ){
			numW = numElemW;
		}else{
			numW = nCalcPointW;
		}
	
		//const MeshDataBrickElement* const pMeshDataBrickElement = MeshDataBrickElement::getInstance();
		const MeshData* const pMeshData = ( ( AnalysisControl::getInstance() )->getPointerOfMeshData() );

		for( int i = 0; i < numH; ++i ){

			for( int j = 0; j < numW; ++j ){

				int firstElem = 0;
				int increment = 1;
				if( nCalcPointW == -1 ){
					firstElem = j;
					increment = numElemW;
				}

				for( int iElem = firstElem; iElem < nElem; iElem += increment ){

					//const int node1 = pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4    ];
					//const int node3 = pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem * 4 + 2];
					const int node1 = pMeshData->getNodesOfElementsBoundaryPlanes( m_planeID, iElem, 0 );
					const int node3 = pMeshData->getNodesOfElementsBoundaryPlanes( m_planeID, iElem, 2 ) ;

					double wmin(0);
					double wmax(0);				
					if(  m_planeID == MeshData::YZMinus || m_planeID == MeshData::YZPlus ){//YZ Plane
						//wmin = pMeshDataBrickElement->m_yCoordinatesOfNodes[node3];
						//wmax = pMeshDataBrickElement->m_yCoordinatesOfNodes[node1];
						wmin = pMeshData->getYCoordinatesOfNodes(node3);
						wmax = pMeshData->getYCoordinatesOfNodes(node1);
					}else{//XZ Plane
						//wmin = pMeshDataBrickElement->m_xCoordinatesOfNodes[node3];
						//wmax = pMeshDataBrickElement->m_xCoordinatesOfNodes[node1];
						wmin = pMeshData->getXCoordinatesOfNodes(node3);
						wmax = pMeshData->getXCoordinatesOfNodes(node1);
					}
					//const double hmin = pMeshDataBrickElement->m_zCoordinatesOfNodes[node3];
					//const double hmax = pMeshDataBrickElement->m_zCoordinatesOfNodes[node1];
					const double hmin = pMeshData->getZCoordinatesOfNodes(node3);
					const double hmax = pMeshData->getZCoordinatesOfNodes(node1);

					bool includedH(false);
					if( upperOrLower == 0 ){// Calculate EM field by upper side element if the point is on the elemnet boundary					
						if( CalcPointH[i] > hmin && CalcPointH[i] <= hmax + CommonParameters::EPS ){
							includedH = true;
						}
					}else{// Calculate EM field by lower side element if the point is on the elemnet boundary
						if( CalcPointH[i] >= hmin - CommonParameters::EPS && CalcPointH[i] < hmax ){
							includedH = true;
						}
					}

					bool includedW(false);
					if( nCalcPointW == -1 ){
						includedW = true;
					}else{
						if( leftOrRight == 0 ){// Calculate EM field by left side element if the point is on the elemnet boundary					
							if( CalcPointW[j] > wmin && CalcPointW[j] <= wmax + CommonParameters::EPS ){
								includedW = true;
							}
						}else{// Calculate EM field by right side element if the point is on the right boundary
							if( CalcPointW[j] >= wmin - CommonParameters::EPS && CalcPointW[j] < wmax ){					
								includedW = true;
							}
						}
					}

					if( includedH && includedW ){

						double hcenter = ( hmin + hmax ) * 0.5;
						double calcPointHFinal = CalcPointH[i];
						double hLocal = 2 * ( calcPointHFinal - hcenter ) / calcHeight(iElem, pMeshDataBrickElement);
						if( hLocal > 1.0 ){
							hLocal = 1.0;
						}
						if( hLocal < -1.0 ){
							hLocal = -1.0;
						}

						double wcenter = ( wmin + wmax ) * 0.5;
						double calcPointWFinal(0);
						if( nCalcPointW == -1 ){
							calcPointWFinal = wcenter;
						}else{
							calcPointWFinal = CalcPointW[j];
						}
						double wLocal = 2 * ( calcPointWFinal - wcenter ) / calcWidth(iElem, pMeshDataBrickElement);
						if( wLocal > 1.0 ){
							wLocal = 1.0;
						}
						if( wLocal < -1.0 ){
							wLocal = -1.0;
						}

						//std::complex<double> V(NULL);
						//std::complex<double> J(NULL);
						//std::complex<double> I(NULL);
						//std::complex<double> Z(NULL);

						if( type == Forward2D::NODE_BASED_SECOND_ORDER || type == EDGE_BASED_FIRST_ORDER ){
							if ( imode == CommonParameters::TE_MODE ){
								OutputFiles::m_logFile << "Error : Only TM mode can be treated for edge-based element." << std::endl;
								exit(1);
							}
						}
						//switch( type ){

						//	case NODE_BASED_FIRST_ORDER:						
						//		V = calcValueV1stOrder( iPol,       iElem, wLocal, hLocal );
						//		J = calcValueJ1stOrder( iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		I = calcValueI1stOrder( iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		break;
						//	case NODE_BASED_SECOND_ORDER:
						//		V = calcValueV2ndOrder( iPol,       iElem, wLocal, hLocal );
						//		J = calcValueJ2ndOrder( iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		I = calcValueI2ndOrder( iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		break;
						//	case EDGE_BASED_ZEROTH_ORDER:

						//		V = calcValueMagneticField0thOrderEdgeBased(           iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		J = calcValueElectricFieldHorizontal0thOrderEdgeBased( iPol,       iElem, wLocal, hLocal );
						//		I = calcValueElectricFieldVertical0thOrderEdgeBased(   iPol,       iElem, wLocal, hLocal );
						//		break;
						//	case EDGE_BASED_FIRST_ORDER:
						//		if ( imode == CommonParameters::TE_MODE ){
						//			OutputFiles::m_logFile << "Error : Only TM mode can be treated for 1st order edge-based element." << std::endl;
						//			exit(1);
						//		}
						//		V = calcValueMagneticField1stOrderEdgeBased(           iPol, freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//		J = calcValueElectricFieldHorizontal1stOrderEdgeBased( iPol,       iElem, wLocal, hLocal );
						//		I = calcValueElectricFieldVertical1stOrderEdgeBased(   iPol,       iElem, wLocal, hLocal );
						//		break;
						//	default:
						//		OutputFiles::m_csvFileFor2DFwd << "Warning : Type is unsupported value in output2DResult !! type = " << type << std::endl;
						//		break;
						//}

						//V = calcValueV( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//J = calcValueJ( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );
						//I = calcValueI( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );

						std::complex<double> Ex = calcEx( iElem, wLocal, hLocal, pMeshDataBrickElement ); 
						std::complex<double> Ey = calcEy( iElem, wLocal, hLocal, pMeshDataBrickElement ); 
						std::complex<double> Ez = calcEz( iElem, wLocal, hLocal, pMeshDataBrickElement ); 
						std::complex<double> Hx = calcHx( freq, iElem, wLocal, hLocal, pMeshDataBrickElement ); 
						std::complex<double> Hy = calcHy( freq, iElem, wLocal, hLocal, pMeshDataBrickElement ); 
						std::complex<double> Hz = calcHz( freq, iElem, wLocal, hLocal, pMeshDataBrickElement );


						//if(  m_planeID == MeshData::ZXMinus || m_planeID == MeshData::ZXPlus ){//ZX Plane
						//	if ( imode == CommonParameters::TM_MODE ){
						//		V *= std::complex<double>(-1.0,0.0);
						//	}else{
						//		J *= std::complex<double>(-1.0,0.0);
						//	}
						//}
						//if ( imode == CommonParameters::TM_MODE ){
						//	Z = -J/V;
						//}else if ( imode == CommonParameters::TE_MODE ){
						//	Z = V/J;
						//}

						std::complex<double> Z(0.0,0.0);
						if(  m_planeID == MeshData::ZXMinus || m_planeID == MeshData::ZXPlus ){//ZX Plane
							if( imode == CommonParameters::TM_MODE ){
								Ey *= std::complex<double>(-1.0, 0.0);
								Z = Ex / Hy;
							}else{
								Hx *= std::complex<double>(-1.0, 0.0);
								Z = Ey / Hx;
							}
						}else{//YZ Plane
							if( imode == CommonParameters::TM_MODE ){
								Z = Ey / Hx;
							}else{
								Z = Ex / Hy;
							}
						}

						const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency
						double rho = static_cast<double>( std::pow(std::abs(Z), 2.0) )/ ( CommonParameters::mu * omega );
						double phase = atan2( imag(Z), real(Z) ) * CommonParameters::rad2deg;
						//OutputFiles::m_csvFileFor2DFwd << calcPointWFinal / 1000.0 << "," << calcPointHFinal / 1000.0 << ","
						//   					   << real(V) << "," << imag(V) << ","
						//						   << real(J) << "," << imag(J) << ","
						//						   << real(I) << "," << imag(I) << ","
						//						   << real(Z) << "," << imag(Z) << ","
						//						   << rho << "," << phase << "," << std::endl;

						//fprintf( OutputFiles::m_csvFileFor2DFwd, "%10.3lf, %10.3lf, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
						//	calcPointWFinal/1000.0, calcPointHFinal/1000.0, real(V), imag(V), real(J), imag(J), real(I), imag(I), real(Z), imag(Z), rho, phase );

						fprintf( OutputFiles::m_csvFileFor2DFwd, "%10.3lf, %10.3lf, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
							CalcPointW[j]/1000.0, CalcPointH[i]/1000.0, real(Ex), imag(Ex), real(Ey), imag(Ey), real(Ez), imag(Ez), real(Hx), imag(Hx), real(Hy), imag(Hy), real(Hz), imag(Hz), real(Z), imag(Z), rho, phase );

					}

				}

			}
		}

		if( nCalcPointH != NULL ){
			delete[] CalcPointH;
		}
		if( nCalcPointW != NULL ){
			delete[] CalcPointW;
		}
	}

}

//// Make vtk file of boudary planes with point data
//// [Input] : 1) iPol : Polarization of analysis
////           2) nodesGlobal2FullModel : Array converting global node IDs of boudary plane to the ones of Full (3D) model
////           3) result : solution vector of 2D analysis
//// [Output] : parameter Gamma of Rodi(1976)
//void Forward2DSquareElement::makeVTKBoundaryPlane( const int iPol, const int* nodesGlobal2FullModel, const std::complex<double>* result ){
//
//	const MeshDataBrickElement* pMeshDataBrickElement = MeshDataBrickElement::getInstance();
//
//	 //--- Calculate element division of 2D
//	const int numElemW = calcNumElemHorizontal();
//	const int numElemH = calcNumElemVertical();
//
//	const int nElem = pMeshDataBrickElement->m_numElemOnBoundaryPlanes[m_planeID]; // Total number of elements on the plane
//	const int nEq = ( 2 * numElemW + 1 ) * ( 2 * numElemH + 1 ) - numElemW * numElemH;// Number of equations	
//	const int nNode = ( numElemH + 1 )*( numElemW + 1 );
//
//	// VTK file
//	std::ostringstream fileName;
//	fileName << "BoundaryPlane_Pol" << iPol << "_Plane" << m_planeID  << ".vtk";
//	std::ofstream vtkFile( fileName.str().c_str(), std::ios::out );
//
//	vtkFile << "# vtk DataFile Version 2.0" << std::endl;
//	vtkFile << "BoundaryPlane_Pol" << iPol << "_Plane" << m_planeID  <<  std::endl;
//	vtkFile << "ASCII" << std::endl;
//	vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
//	vtkFile << "POINTS " << nNode << " float" << std::endl;
//
//	//int* nodeFullModel2BoundaryPlane = new int[ nNode ];
//	std::map<int, int> nodeFullModel2BoundaryPlane;
//	std::map<int, int> nodeBoundaryPlane2FullModel;
//	int icount(0);
//	for( int i = 0; i < nEq; ++i ){
//		if( nodesGlobal2FullModel[i] < 0 ){
//			continue;
//		}
//
//		const int iNode = nodesGlobal2FullModel[i];
//
//		vtkFile << pMeshDataBrickElement->m_xCoordinatesOfNodes[iNode] << " " << pMeshDataBrickElement->m_yCoordinatesOfNodes[iNode] << " " << pMeshDataBrickElement->m_zCoordinatesOfNodes[iNode] << std::endl;
//		//nodeFullModel2BoundaryPlane[iNode] = icount;
//		nodeFullModel2BoundaryPlane.insert( std::map<int, int>::value_type( iNode, icount ) );
//		nodeBoundaryPlane2FullModel.insert( std::map<int, int>::value_type( icount, iNode ) );
//		++icount;
//	}
//
//	vtkFile << "CELLS " << nElem << " " << nElem*5 << std::endl;
//
//	for( int iElem = 0; iElem < nElem; ++iElem ){
//		vtkFile << 4 << " "
//				<< nodeFullModel2BoundaryPlane[ pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem*4  ] ] << " " 
//				<< nodeFullModel2BoundaryPlane[ pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem*4+1] ] << " " 
//				<< nodeFullModel2BoundaryPlane[ pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem*4+2] ] << " " 
//				<< nodeFullModel2BoundaryPlane[ pMeshDataBrickElement->m_nodesOfElementsBoundaryPlanes[m_planeID][iElem*4+3] ] << std::endl;
//	}
//
//	vtkFile << "CELL_TYPES " << nElem << std::endl;
//	for( int iElem = 0 ; iElem < nElem; ++iElem ){
//		vtkFile << "9" << std::endl;
//	}
//
//	//vtkFile << "POINT_DATA " << nNode << std::endl;
//	//vtkFile << "SCALARS NodeID(FullModel) float" <<  std::endl;
//	//vtkFile << "LOOKUP_TABLE default" <<  std::endl;
//	//for( int i = 0 ; i < nNode; ++i ){
//	//	vtkFile << nodeBoundaryPlane2FullModel[i] << std::endl;
//	//}
//
//	vtkFile << "POINT_DATA " << nNode << std::endl;
//	vtkFile << "SCALARS data float 4" <<  std::endl;
//	vtkFile << "LOOKUP_TABLE default" <<  std::endl;
//	for( int i = 0 ; i < nNode; ++i ){
//		vtkFile << real( result[i] ) << " "
//			    << imag( result[i] ) << " "
//				<<  std::abs( result[i] ) << " "
//				<<  std::arg( result[i] ) * 180.0 / CommonParameters::PI << std::endl;
//	}
//
//	vtkFile.close();
//}

//// Copy constructer
//Forward2DSquareElement::Forward2DSquareElement(const Forward2DSquareElement& rhs){
//	std::cerr << "Error : Copy constructer of the class Forward2DSquareElement is not implemented." << std::endl;
//	exit(1);
//}
//
//// Copy assignment operator
//Forward2DSquareElement& Forward2DSquareElement::operator=(const Forward2DSquareElement& rhs){
//	std::cerr << "Error : Copy constructer of the class Forward2DSquareElement is not implemented." << std::endl;
//	exit(1);
//}
