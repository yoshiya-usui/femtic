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
#include <stddef.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <map>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "AnalysisControl.h"
#include "MeshDataBrickElement.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"
#include "OutputFiles.h"

// Constructer
MeshDataBrickElement::MeshDataBrickElement():
	m_numElemX(NULL),
	m_numElemY(NULL),
	m_numElemZ(NULL),
	m_numAirLayer(NULL),
	m_edgeLength(NULL)
{
	m_numNodeOneElement = 8;
	m_numEdgeOneElement = 12;
	m_numNodeOnFaceOneElement = 4;
	m_numNeighborElement = 6;

	for ( int i = 0; i < 6; ++i ){
		m_nodesOfElementsBoundaryPlanes[i] = NULL;
	}

}

// Destructer
MeshDataBrickElement::~MeshDataBrickElement(){

	if( m_edgeLength != NULL){
		delete[] m_edgeLength;
		m_edgeLength = NULL;
	}

	for ( int i = 0; i < 6; ++i ){
		m_nodesOfElementsBoundaryPlanes[i] = NULL;
	}

	for ( int i = 0; i < 6; ++i ){
		if( m_nodesOfElementsBoundaryPlanes[i] != NULL){
			delete[] m_nodesOfElementsBoundaryPlanes[i];
			m_nodesOfElementsBoundaryPlanes[i] = NULL;
		}
	}

}

// Copy constructer
MeshDataBrickElement::MeshDataBrickElement(const MeshDataBrickElement& rhs){
	std::cerr << "Error : Copy constructer of the class MeshDataBrickElement is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
MeshDataBrickElement& MeshDataBrickElement::operator=(const MeshDataBrickElement& rhs){
	std::cerr << "Error : Assignment operator of the class MeshDataBrickElement is not implemented." << std::endl;
	exit(1);
}

// Input mesh data from "mesh.dat"
void MeshDataBrickElement::inputMeshData(){

	std::ifstream inFile( "mesh.dat", std::ios::in );
	if( inFile.fail() )
	{
		std::cerr << "File open error : mesh.dat !!" << std::endl;
		exit(1);
	}

	std::string sbuf;
	inFile >> sbuf;

	if( sbuf.substr(0,4).compare("HEXA") != 0 ){
		OutputFiles::m_logFile << "Mesh data written in mesh.dat is different from the ones of hexahedral element !!" << std::endl;
		exit(1);
	}

	int nX(0);
	int nY(0);	
	int nZ(0);
	int nAirLayer(0);

	inFile >> nX;
	inFile >> nY;
	inFile >> nZ;
	inFile >> nAirLayer;
	if( nX <= 0 ){
		std::cerr << "Division number of X direction is less than or equal to zero ! : " << nX << std::endl;
		exit(1);
	}
	if( nY <= 0 ){
		std::cerr << "Division number of Y direction is less than or equal to zero ! : " << nY << std::endl;
		exit(1);
	}
	if( nZ <= 0 ){
		std::cerr << "Division number of Z direction is less than or equal to zero ! : " << nZ << std::endl;
		exit(1);
	}
	if( nAirLayer <= 0 ){
		std::cerr << "Number of the air layer is less than or equal to zero ! : " << nAirLayer << std::endl;
		exit(1);
	}
	
	//setMeshDivisons( nX, nY, nZ );// Set mesh division numbers
	m_numElemX = nX;
	m_numElemY = nY;
	m_numElemZ = nZ;
	m_numElemTotal = nX * nY * nZ;
	m_numNodeTotal = (nX+1) * (nY+1) * (nZ+1);
	m_numElemOnBoundaryPlanes[0] = nY * nZ;
	m_numElemOnBoundaryPlanes[1] = nY * nZ;
	m_numElemOnBoundaryPlanes[2] = nZ * nX;
	m_numElemOnBoundaryPlanes[3] = nZ * nX;
	m_numElemOnBoundaryPlanes[4] = nX * nY;
	m_numElemOnBoundaryPlanes[5] = nX * nY;

	//setNumAirLayer( nAirLayer );// Number of air layer
	m_numAirLayer = nAirLayer;

	if( m_xCoordinatesOfNodes != NULL ){
		delete[] m_xCoordinatesOfNodes;	
	}
	m_xCoordinatesOfNodes = new double[m_numNodeTotal];

	if( m_yCoordinatesOfNodes != NULL ){
		delete[] m_yCoordinatesOfNodes;	
	}
	m_yCoordinatesOfNodes = new double[m_numNodeTotal];

	if( m_zCoordinatesOfNodes != NULL ){
		delete[] m_zCoordinatesOfNodes;	
	}
	m_zCoordinatesOfNodes = new double[m_numNodeTotal];

	const int nNode = ( nX + 1 ) * ( nY + 1 ) * ( nZ + 1 );
	for( int iNode = 0; iNode < nNode; ++iNode ){
		int idum(0);
		inFile >> idum;

		// X coordinates of nodes
		double x(0);
		inFile >> x;
		//setXCoordinatesOfNodes( iNode, x );
		m_xCoordinatesOfNodes[iNode] = x;
		
		// Y coordinates of nodes
		double y(0);
		inFile >> y;
		//setYCoordinatesOfNodes( iNode, y );
		m_yCoordinatesOfNodes[iNode] = y;

		// Z coordinates of nodes
		double z(0);
		inFile >> z;
		//setZCoordinatesOfNodes( iNode, z );
		m_zCoordinatesOfNodes[iNode] = z;

#ifdef _DEBUG_WRITE
		std::cout << idum << " " << x << " " << y << " " << z << std::endl; // For debug
#endif

	}

	const int nElem = nX  * nY * nZ;
	if( m_neighborElements != NULL ){
		delete[] m_neighborElements;
	}
	m_neighborElements = new int[ m_numElemTotal * 6 ];

	if( m_nodesOfElements == NULL ){
		delete[] m_nodesOfElements;
	}
	m_nodesOfElements = new int[ m_numElemTotal * 8 ];

	for( int iElem = 0; iElem < nElem; ++iElem ){
		int idum(0);
		inFile >> idum;

		// IDs of neighbor Elements
		for( int i = 0; i < 6; ++i ){
			int neib(0);
			inFile >> neib;
			//setNeighborElements( iElem, i, neib );
			m_neighborElements[ iElem * 6 + i ] = neib;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << neib << std::endl; // For debug
#endif

		}

		// Nodes of the element
		for( int i = 0; i < 8; ++i ){
			int node(0);
			inFile >> node;
			//setNodesOfElements( iElem, i, node );
			m_nodesOfElements[ iElem * 8 + i ] = node;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << node << std::endl; // For debug
#endif

		}

	}
	

	for( int iPlane = 0; iPlane < 6; ++iPlane ){// Loop of boundary planes

		int nElemOnPlane;
		inFile >> nElemOnPlane;

#ifdef _DEBUG_WRITE
		std::cout << nElemOnPlane << std::endl; // For debug
#endif

		if( nElemOnPlane != getNumElemOnBoundaryPlanes( iPlane ) ){
			OutputFiles::m_logFile << " Error : Total number of elements belonging to the plane "
				<< iPlane << " is inconsistent with division numbers." << std::endl;
			exit(1);
		}

		if( m_elemBoundaryPlanes[iPlane] != NULL ){
			delete [] m_elemBoundaryPlanes[iPlane];
		}
		m_elemBoundaryPlanes[iPlane] = new int[ nElemOnPlane ];

		if( m_nodesOfElementsBoundaryPlanes[iPlane] != NULL ){
			delete [] m_nodesOfElementsBoundaryPlanes[iPlane];	
		}
		m_nodesOfElementsBoundaryPlanes[iPlane] = new int[ nElemOnPlane * 4 ];

		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){		
			// Set elements belonging to the boundary planes
			int elemID(0);
			inFile >> elemID;

#ifdef _DEBUG_WRITE
			std::cout << iElem << " " << elemID << std::endl; // For debug
#endif

			//setElemBoundaryPlanes( iPlane, iElem, elemID );
			m_elemBoundaryPlanes[iPlane][iElem] = elemID;

			// Set nodes of elements belonging to the boundary planes
			for( int iNode = 0; iNode < 4; ++iNode ){		
				int nodeID(0);
				inFile >> nodeID;

#ifdef _DEBUG_WRITE
				std::cout << iNode << " " << nodeID << std::endl; // For debug
#endif

				//setNodesOfElementsBoundaryPlanes( iPlane, iElem, iNode, nodeID );
				m_nodesOfElementsBoundaryPlanes[iPlane][ iElem * 4 + iNode ] = nodeID;
				
			}
		}
		
	}

	inFile.close();

	//// Output mesh data to VTK file
	//if( OutputFiles::m_vtkFile.is_open() ){
	//	outputMeshDataToVTK();
	//}

}

// Find element including a point
// <Input> : locX, locY, locZ
// <Output>: localCoordX, localCoordY, localCoordZ
int MeshDataBrickElement::findElementIncludingPoint( const double locX, const double locY, const double locZ, double& localCoordX, double& localCoordY, double& localCoordZ, 
	const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod) const{

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	const AnalysisControl::UseBackwardOrForwardElement useBackwardOrForwardElement = pAnalysisControl->getUseBackwardOrForwardElement();
	const int useBackwardOrForwardElementDirectionZ = useUpperElem ? AnalysisControl::BACKWARD_ELEMENT : AnalysisControl::FORWARD_ELEMENT;
	const int useforwardOrBackwardElement[3] = { useBackwardOrForwardElement.directionX, useBackwardOrForwardElement.directionY, useBackwardOrForwardElementDirectionZ };

	const double coordinatesOfPoints[3] = { locX, locY, locZ };
	
	// Very small value
	const double eps = 1.0E-8;

	int elementIncludingPoint = 0;

	for( int i = 0; i < 3; ++i ){// Loop for each direction
		int icount = 0;
		while(1){

			const int node0 = m_nodesOfElements[ elementIncludingPoint * 8     ];
			const int node6 = m_nodesOfElements[ elementIncludingPoint * 8 + 6 ];

			const double xyzmin[3] = { m_xCoordinatesOfNodes[node0], m_yCoordinatesOfNodes[node0], m_zCoordinatesOfNodes[node0] };
			const double xyzmax[3] = { m_xCoordinatesOfNodes[node6], m_yCoordinatesOfNodes[node6], m_zCoordinatesOfNodes[node6] };

			if( modLoc ){
				locXMod = 0.5 * ( xyzmin[0] + xyzmax[0] );
				locYMod = 0.5 * ( xyzmin[1] + xyzmax[1] );
			}

			const int elementIncludingPointOld = elementIncludingPoint;
			if( useforwardOrBackwardElement[i] == AnalysisControl::BACKWARD_ELEMENT ){// Calculate EM field by backward side element if the point is on the elemnet boundary
				if( coordinatesOfPoints[i] > xyzmin[i] && coordinatesOfPoints[i] <= xyzmax[i] + eps ){
					break;
				}else if( coordinatesOfPoints[i] > xyzmax[i] + eps ){
					elementIncludingPoint = m_neighborElements[ 6*elementIncludingPointOld + 2*i + 1 ];					
				}else{
					elementIncludingPoint = m_neighborElements[ 6*elementIncludingPointOld + 2*i ];
				}
			}else{// Calculate EM field by forward side element if the point is on the elemnet boundary					
				if( coordinatesOfPoints[i] >= xyzmin[i] - eps && coordinatesOfPoints[i] < xyzmax[i] ){
					break;
				}else if( coordinatesOfPoints[i] >= xyzmax[i] ){
					elementIncludingPoint = m_neighborElements[ 6*elementIncludingPointOld + 2*i + 1 ];
				}else{
					elementIncludingPoint = m_neighborElements[ 6*elementIncludingPointOld + 2*i ];
				}
			}

			if( elementIncludingPoint < 0 || icount >= m_numElemTotal ){
				OutputFiles::m_logFile << " Error : Could not find element including point ( X = " << locX << ", Y= " << locY << ", Z = " << locZ << " )." << std::endl;
				exit(1);
				return -1;
			}

			++icount;
		}
	}

	if( modLoc ){
		getLocalCoordinateValues( elementIncludingPoint, locXMod, locYMod, locZ, localCoordX, localCoordY, localCoordZ );
	}else{
		getLocalCoordinateValues( elementIncludingPoint, locX, locY, locZ, localCoordX, localCoordY, localCoordZ );
	}

	return elementIncludingPoint;

}

// Find element including a point on the surface of the earth
// <Input> : locX, locY
// <Output>: localCoordX, localCoordY, localCoordZ
int MeshDataBrickElement::findElementIncludingPointOnSurface( const double locX, const double locY, double& localCoordX, double& localCoordY, double& localCoordZ,
	const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod ) const{

	return findElementIncludingPoint(locX, locY, 0.0, localCoordX, localCoordY, localCoordZ, useUpperElem, modLoc, locXMod, locYMod );

}

// Find elements including dipole on the surface of the earth
// <Input> : locXStart, locYStart, locXEnd, locYEnd
// <Output>: elements, localCoordXStartPoint, localCoordYStartPoint, localCoordXEndPoint, localCoordYEndPoint
void MeshDataBrickElement::findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd,
		std::vector<int>& elements, std::vector<double>& localCoordXStartPoint, std::vector<double>& localCoordYStartPoint,	std::vector<double>& localCoordXEndPoint, std::vector<double>& localCoordYEndPoint ) const{

	const double start2endX = locXEnd - locXStart;
	const double start2endY = locYEnd - locYStart;
	const double slopeYX = start2endY / start2endX;
	const double slopeXY = start2endX / start2endY;

	double localCoord[3] = { 0.0, 0.0, 0.0 };
	const double thresholdVal = 1.0E-6;
	double dummy(0.0);

	findElementIncludingPoint( locXStart, locYStart, 0.0, localCoord[0], localCoord[1], localCoord[2], false, false, dummy, dummy );
	double locXStartEdited(locXStart);
	double locYStartEdited(locYStart);
	if( fabs( localCoord[0] - 1.0 ) < thresholdVal || fabs( localCoord[0] + 1.0 ) < thresholdVal ){
		const double val = 1.0e-3;
		locXStartEdited = start2endX > 0 ? locXStart + val : locXStart - val;
	}
	if( fabs( localCoord[1] - 1.0 ) < thresholdVal || fabs( localCoord[1] + 1.0 ) < thresholdVal ){
		const double val = 1.0e-3;
		locYStartEdited = start2endY > 0 ? locYStart + val : locYStart - val;
	}
	const int elementIncludingStartPoint = findElementIncludingPoint( locXStartEdited, locYStartEdited, 0.0, localCoord[0], localCoord[1], localCoord[2], false, false, dummy, dummy );
	findElementIncludingPoint( locXEnd, locYEnd, 0.0, localCoord[0], localCoord[1], localCoord[2], false, false, dummy, dummy );

	double locXEndEdited(locXEnd);
	double locYEndEdited(locYEnd);
	if( fabs( start2endX ) < thresholdVal ){// If the vector from start to end point is parallel to a boundary of elements.
		locXEndEdited = locXStartEdited;
	}else{
		if( fabs( localCoord[0] - 1.0 ) < thresholdVal || fabs( localCoord[0] + 1.0 ) < thresholdVal ){
			const double val = 1.0e-3;
			locXEndEdited = start2endX > 0 ? locXEnd - val : locXEnd + val;
		}
	}
	if( fabs( start2endY ) < thresholdVal ){// If the vector from start to end point is parallel to a boundary of elements.
		locYEndEdited = locYStartEdited;
	}else{
		if( fabs( localCoord[1] - 1.0 ) < thresholdVal || fabs( localCoord[1] + 1.0 ) < thresholdVal ){
			const double val = 1.0e-3;
			locYEndEdited = start2endY > 0 ? locYEnd - val : locYEnd + val;
		}
	}
	const int elementIncludingEndPoint   = findElementIncludingPoint( locXEndEdited, locYEndEdited, 0.0, localCoord[0], localCoord[1], localCoord[2], false, false, dummy, dummy );

	// Threshold value
	const double eps = 1.0e-3;

	double XCoordinateOfIntersectionPointPre = locXStart;
	double YCoordinateOfIntersectionPointPre = locYStart;

	for( int iElem = elementIncludingStartPoint; iElem != elementIncludingEndPoint; ){

		const int node0 = m_nodesOfElements[ iElem * 8     ];
		const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];

		coordinateValue minimumCoordinateValue = { m_xCoordinatesOfNodes[node0], m_yCoordinatesOfNodes[node0], m_zCoordinatesOfNodes[node0] };
		coordinateValue maximumCoordinateValue = { m_xCoordinatesOfNodes[node6], m_yCoordinatesOfNodes[node6], m_zCoordinatesOfNodes[node6] };

		// Neighbor element of -X side
		if( m_neighborElements[ 6*iElem ] >= 0 && start2endX < 0.0 ){

			const double XCoordinateOfIntersectionPoint = minimumCoordinateValue.X;
			const double YCoordinateOfIntersectionPoint = slopeYX * ( minimumCoordinateValue.X - locXStart ) + locYStart;

			if( YCoordinateOfIntersectionPoint >= minimumCoordinateValue.Y && YCoordinateOfIntersectionPoint <= maximumCoordinateValue.Y ){
				// Find element including dipole
				if( hypot( ( XCoordinateOfIntersectionPoint - XCoordinateOfIntersectionPointPre ), ( YCoordinateOfIntersectionPoint - YCoordinateOfIntersectionPointPre ) ) > eps ){
					elements.push_back( iElem );
					double localCoordX(0.0);
					double localCoordY(0.0);
					double localCoordZ(0.0);
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPointPre, YCoordinateOfIntersectionPointPre, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXStartPoint.push_back( localCoordX );
					localCoordYStartPoint.push_back( localCoordY );
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPoint, YCoordinateOfIntersectionPoint, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXEndPoint.push_back( localCoordX );
					localCoordYEndPoint.push_back( localCoordY );
				}
				iElem = m_neighborElements[ 6*iElem ];
				XCoordinateOfIntersectionPointPre = XCoordinateOfIntersectionPoint;
				YCoordinateOfIntersectionPointPre = YCoordinateOfIntersectionPoint;
				continue;
			}

		}

		// Neighbor element of +X side
		if( m_neighborElements[ 6*iElem + 1 ] >= 0 &&  start2endX > 0.0 ){

			const double XCoordinateOfIntersectionPoint = maximumCoordinateValue.X;
			const double YCoordinateOfIntersectionPoint = slopeYX * ( maximumCoordinateValue.X - locXStart ) + locYStart;

			if( YCoordinateOfIntersectionPoint >= minimumCoordinateValue.Y && YCoordinateOfIntersectionPoint <= maximumCoordinateValue.Y ){
				// Find element including dipole
				if( hypot( ( XCoordinateOfIntersectionPoint - XCoordinateOfIntersectionPointPre ), ( YCoordinateOfIntersectionPoint - YCoordinateOfIntersectionPointPre ) ) > eps ){
					elements.push_back( iElem );
					double localCoordX(0.0);
					double localCoordY(0.0);
					double localCoordZ(0.0);
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPointPre, YCoordinateOfIntersectionPointPre, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXStartPoint.push_back( localCoordX );
					localCoordYStartPoint.push_back( localCoordY );
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPoint, YCoordinateOfIntersectionPoint, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXEndPoint.push_back( localCoordX );
					localCoordYEndPoint.push_back( localCoordY );
				}
				iElem = m_neighborElements[ 6*iElem + 1 ];
				XCoordinateOfIntersectionPointPre = XCoordinateOfIntersectionPoint;
				YCoordinateOfIntersectionPointPre = YCoordinateOfIntersectionPoint;
				continue;
			}

		}

		// Neighbor element of - Y side
		if( m_neighborElements[ 6*iElem + 2 ] >= 0 && start2endY < 0.0 ){

			const double XCoordinateOfIntersectionPoint = slopeXY * ( minimumCoordinateValue.Y - locYStart ) + locXStart;
			const double YCoordinateOfIntersectionPoint = minimumCoordinateValue.Y;

			if( XCoordinateOfIntersectionPoint >= minimumCoordinateValue.X && XCoordinateOfIntersectionPoint <= maximumCoordinateValue.X ){
				// Find element including dipole
				if( hypot( ( XCoordinateOfIntersectionPoint - XCoordinateOfIntersectionPointPre ), ( YCoordinateOfIntersectionPoint - YCoordinateOfIntersectionPointPre ) ) > eps ){
					elements.push_back( iElem );
					double localCoordX(0.0);
					double localCoordY(0.0);
					double localCoordZ(0.0);
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPointPre, YCoordinateOfIntersectionPointPre, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXStartPoint.push_back( localCoordX );
					localCoordYStartPoint.push_back( localCoordY );
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPoint, YCoordinateOfIntersectionPoint, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXEndPoint.push_back( localCoordX );
					localCoordYEndPoint.push_back( localCoordY );
				}
				iElem = m_neighborElements[ 6*iElem + 2 ];
				XCoordinateOfIntersectionPointPre = XCoordinateOfIntersectionPoint;
				YCoordinateOfIntersectionPointPre = YCoordinateOfIntersectionPoint;
				continue;
			}
		}

		// Neighbor element of + Y side
		if( m_neighborElements[ 6*iElem + 3 ] >= 0 && start2endY > 0.0  ){
			
			const double XCoordinateOfIntersectionPoint = slopeXY * ( maximumCoordinateValue.Y - locYStart ) + locXStart;
			const double YCoordinateOfIntersectionPoint = maximumCoordinateValue.Y;

			if( XCoordinateOfIntersectionPoint >= minimumCoordinateValue.X && XCoordinateOfIntersectionPoint <= maximumCoordinateValue.X ){
				// Find element including dipole
				if( hypot( ( XCoordinateOfIntersectionPoint - XCoordinateOfIntersectionPointPre ), ( YCoordinateOfIntersectionPoint - YCoordinateOfIntersectionPointPre ) ) > eps ){
					elements.push_back( iElem );
					double localCoordX(0.0);
					double localCoordY(0.0);
					double localCoordZ(0.0);
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPointPre, YCoordinateOfIntersectionPointPre, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXStartPoint.push_back( localCoordX );
					localCoordYStartPoint.push_back( localCoordY );
					getLocalCoordinateValues( iElem, XCoordinateOfIntersectionPoint, YCoordinateOfIntersectionPoint, -1.0, localCoordX, localCoordY, localCoordZ );
					localCoordXEndPoint.push_back( localCoordX );
					localCoordYEndPoint.push_back( localCoordY );
				}
				iElem = m_neighborElements[ 6*iElem + 3 ];
				XCoordinateOfIntersectionPointPre = XCoordinateOfIntersectionPoint;
				YCoordinateOfIntersectionPointPre = YCoordinateOfIntersectionPoint;
				continue;
			}
		}

		OutputFiles::m_logFile << " Error : Could not find element including dipole ( " << locXStart << ", " << locYStart << " ) => ( " << locXEnd << ", " << locYEnd << " )." << std::endl;
		exit(1);
					
	}

	elements.push_back( elementIncludingEndPoint );
	double localCoordX(0.0);
	double localCoordY(0.0);
	double localCoordZ(0.0);
	getLocalCoordinateValues( elementIncludingEndPoint, XCoordinateOfIntersectionPointPre, YCoordinateOfIntersectionPointPre, -1.0, localCoordX, localCoordY, localCoordZ );
	localCoordXStartPoint.push_back( localCoordX );
	localCoordYStartPoint.push_back( localCoordY );
	getLocalCoordinateValues( elementIncludingEndPoint, locXEnd, locYEnd, -1.0, localCoordX, localCoordY, localCoordZ );
	localCoordXEndPoint.push_back( localCoordX );
	localCoordYEndPoint.push_back( localCoordY );
	
}

// Get length of the edges parallel to X coordinate
double MeshDataBrickElement::getEdgeLengthX( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node1 = m_nodesOfElements[ iElem * 8 + 1 ];

	return std::fabs( m_xCoordinatesOfNodes[node1] - m_xCoordinatesOfNodes[node0] );

}

// Get length of the edges parallel to Y coordinate
double MeshDataBrickElement::getEdgeLengthY( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node3 = m_nodesOfElements[ iElem * 8 + 3 ];

	return std::fabs( m_yCoordinatesOfNodes[node3] - m_yCoordinatesOfNodes[node0] );

}

// Get length of the edges parallel to Z coordinate
double MeshDataBrickElement::getEdgeLengthZ( const int iElem ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node4 = m_nodesOfElements[ iElem * 8 + 4 ];

	return std::fabs( m_zCoordinatesOfNodes[node4] - m_zCoordinatesOfNodes[node0] );

}

// Get global X coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordX( const int iElem, double localCoordX ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_xCoordinatesOfNodes[node0];
	const double coordMax = m_xCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordX + ( coordMax + coordMin ) * 0.5 ;

}

// Get global Y coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordY( const int iElem, double localCoordY ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_yCoordinatesOfNodes[node0];
	const double coordMax = m_yCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordY + ( coordMax + coordMin ) * 0.5 ;

}

// Get global Z coordinate from local coordinate
double MeshDataBrickElement::calcGlobalCoordZ( const int iElem, double localCoordZ ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	const double coordMin = m_zCoordinatesOfNodes[node0];
	const double coordMax = m_zCoordinatesOfNodes[node6];

	return ( coordMax - coordMin ) * 0.5 * localCoordZ + ( coordMax + coordMin ) * 0.5 ;

}

// Get number of Elements parallel to X direction
int MeshDataBrickElement::getNumElemX() const{
	return m_numElemX;
}

// Get number of Elements parallel to Y direction
int MeshDataBrickElement::getNumElemY() const{
	return m_numElemY;
}

// Get number of Elements parallel to Z direction
int MeshDataBrickElement::getNumElemZ() const{
	return m_numElemZ;
}

// Get number of the air layer
int MeshDataBrickElement::getNumAirLayer() const{
	return m_numAirLayer;
}

// Calculate number of edges of X-Y plane
int MeshDataBrickElement::calcNumEdgesOnXYPlane() const{
	return m_numElemX * ( m_numElemY + 1 ) + ( m_numElemX + 1 ) * m_numElemY;
}

// Calculate number of edges of Y-Z plane
int MeshDataBrickElement::calcNumEdgesOnYZPlane() const{
	return m_numElemY * ( m_numElemZ + 1 ) + ( m_numElemY + 1 ) * m_numElemZ;
}

// Calculate number of edges of Z-X plane
int MeshDataBrickElement::calcNumEdgesOnZXPlane() const{
	return m_numElemZ * ( m_numElemX + 1 ) + ( m_numElemZ + 1 ) * m_numElemX;
}

// Get array of nodes of elements belonging to the boundary planes
int MeshDataBrickElement::getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const{

	//if( iElem < 0 || iElem >= getNumElemOnBoundaryPlanes( iPlane ) ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getNodesOfElementsBoundaryPlanes. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iNode < 0 || iNode >= m_numNodeOnFaceOneElement ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getNodesOfElementsBoundaryPlanes. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < getNumElemOnBoundaryPlanes( iPlane ) );
	assert( iNode >= 0 );
	assert( iNode < m_numNodeOnFaceOneElement );

	return m_nodesOfElementsBoundaryPlanes[iPlane][ m_numNodeOnFaceOneElement*iElem + iNode ];
}

// Get mesh type
int MeshDataBrickElement::getMeshType() const{

	return HEXA;

}

// Decide whether specified elements share same edges
bool MeshDataBrickElement::shareSameEdges( const int elemID1, const int elemID2 ) const{

	//if( elemID1 < 0 || elemID1 >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : elemID1 is out of range in shareSameNodes. elemID1 = " << elemID1 << std::endl;
	//	exit(1);
	//}

	//if( elemID2 < 0 || elemID2 >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : elemID2 is out of range in shareSameNodes. elemID2 = " << elemID2 << std::endl;
	//	exit(1);
	//}
	assert( elemID1 >= 0 );
	assert( elemID1 < m_numElemTotal );
	assert( elemID2 >= 0 );
	assert( elemID2 < m_numElemTotal );

	// Edges parallel to X direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		const int nodeID1[2] = { getNodesOfElements( elemID1, 2*iEdge1 ), getNodesOfElements( elemID1, 2*iEdge1+1 ) };

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			const int nodeID2[2] = { getNodesOfElements( elemID2, 2*iEdge2 ), getNodesOfElements( elemID2, 2*iEdge2+1 ) };

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	// Edges parallel to Y direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		int nodeID1[2];
		if( iEdge1 == 0 || iEdge1 == 2 ){
			nodeID1[0] = getNodesOfElements( elemID1, 2*iEdge1     );
			nodeID1[1] = getNodesOfElements( elemID1, 2*iEdge1 + 3 );
		}else{
			nodeID1[0] = getNodesOfElements( elemID1, 2*(iEdge1-1) + 1 );
			nodeID1[1] = getNodesOfElements( elemID1, 2*(iEdge1-1) + 2 );
		}

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			int nodeID2[2];
			if( iEdge2 == 0 || iEdge2 == 2 ){
				nodeID2[0] = getNodesOfElements( elemID2, 2*iEdge2     );
				nodeID2[1] = getNodesOfElements( elemID2, 2*iEdge2 + 3 );
			}else{
				nodeID2[0] = getNodesOfElements( elemID2, 2*(iEdge2-1) + 1 );
				nodeID2[1] = getNodesOfElements( elemID2, 2*(iEdge2-1) + 2 );
			}

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	// Edges parallel to Z direction
	for( int iEdge1 = 0; iEdge1 < 4; ++iEdge1 ){

		const int nodeID1[2] = { getNodesOfElements( elemID1, iEdge1 ), getNodesOfElements( elemID1, iEdge1+4 ) };

		for( int iEdge2 = 0; iEdge2 < 4; ++iEdge2 ){

			const int nodeID2[2] = { getNodesOfElements( elemID2, iEdge2 ), getNodesOfElements( elemID2, iEdge2+4 ) };

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	return false;

}

// Calculate volume of a specified element
double MeshDataBrickElement::calcVolume( const int elemID ) const{
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	return getEdgeLengthX(elemID) * getEdgeLengthY(elemID) * getEdgeLengthZ(elemID);
}

// Calculate area of face
double MeshDataBrickElement::calcAreaOfFace( const int iElem, const int iFace ) const{

	switch (iFace)
	{
		case 0:// Go through
		case 1:
			// Face perpendicular to x-direction
			return getEdgeLengthY(iElem) * getEdgeLengthZ(iElem);
			break;
		case 2:// Go through
		case 3:
			// Face perpendicular to y-direction
			return getEdgeLengthZ(iElem) * getEdgeLengthX(iElem);
			break;
		case 4:// Go through
		case 5:
			// Face perpendicular to z-direction
			return getEdgeLengthX(iElem) * getEdgeLengthY(iElem);
			break;
		default:
			OutputFiles::m_logFile << " Error : Wrong face ID :" << iFace << std::endl;
			exit(1);
			break;
	}

	return -1.0;

}

// Calculate area of face at bottom of mesh
double MeshDataBrickElement::calcAreaOfFaceAtBottomOfMesh( const int iElem ) const{
	const int elemID = getElemBoundaryPlanes(MeshData::XYPlus, iElem);
	const int iFace = 5;
	return calcAreaOfFace(elemID, iFace);
}

// Output vtk file
void MeshDataBrickElement::outputMeshDataToVTK() const{

	if( !OutputFiles::m_vtkFile.is_open() ){
		return;
	}

	OutputFiles::m_vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	OutputFiles::m_vtkFile << "POINTS " << m_numNodeTotal << " float" << std::endl;

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		OutputFiles::m_vtkFile << m_xCoordinatesOfNodes[iNode] << " "
		   	                   << m_yCoordinatesOfNodes[iNode] << " "
			                   << m_zCoordinatesOfNodes[iNode] << std::endl;
	}

	OutputFiles::m_vtkFile << "CELLS " << m_numElemTotal << " " << m_numElemTotal*9 << std::endl;
	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
		OutputFiles::m_vtkFile << 8 << " "
				  				   << m_nodesOfElements[iElem*8    ] << " " 
				  				   << m_nodesOfElements[iElem*8 + 1] << " " 
				  				   << m_nodesOfElements[iElem*8 + 2] << " " 
				  				   << m_nodesOfElements[iElem*8 + 3] << " " 
				  				   << m_nodesOfElements[iElem*8 + 4] << " " 
				  				   << m_nodesOfElements[iElem*8 + 5] << " " 
				  				   << m_nodesOfElements[iElem*8 + 6] << " " 
				  				   << m_nodesOfElements[iElem*8 + 7] << std::endl;
	}

	OutputFiles::m_vtkFile << "CELL_TYPES " << m_numElemTotal << std::endl;
	for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
		OutputFiles::m_vtkFile << "12" << std::endl;
	}

	OutputFiles::m_vtkFile << "CELL_DATA " << m_numElemTotal << std::endl;

}

// Output mesh data to binary file
void MeshDataBrickElement::outputMeshDataToBinary() const{
	
	std::ofstream fout;
	fout.open( "Mesh.geo", std::ios::out | std::ios::binary | std::ios::trunc );

	float* coordX = new float[m_numNodeTotal];
	float* coordY = new float[m_numNodeTotal];
	float* coordZ = new float[m_numNodeTotal];
	float minX(1.0+20);
	float maxX(-1.0+20);
	float minY(1.0+20);
	float maxY(-1.0+20);
	float minZ(1.0+20);
	float maxZ(-1.0+20);
	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		coordX[iNode] = static_cast<float>( m_xCoordinatesOfNodes[iNode] );
		if( coordX[iNode] < minX ){
			minX = coordX[iNode];
		}
		if( coordX[iNode] > maxX ){
			maxX = coordX[iNode];
		}
		coordY[iNode] = static_cast<float>( m_yCoordinatesOfNodes[iNode] );
		if( coordY[iNode] < minY ){
			minY = coordY[iNode];
		}
		if( coordY[iNode] > maxY ){
			maxY = coordY[iNode];
		}
		coordZ[iNode] = static_cast<float>( m_zCoordinatesOfNodes[iNode] );
		if( coordZ[iNode] < minZ ){
			minZ = coordZ[iNode];
		}
		if( coordZ[iNode] > maxZ ){
			maxZ = coordZ[iNode];
		}
	}

	char line[80];
	strcpy( line, "C Binary" );
	fout.write( line, 80 );

	strcpy( line, "Mesh inputted to FEMTIC" );
	fout.write( line, 80 );

	strcpy( line, "Mesh inputted to FEMTIC" );
	fout.write( line, 80 );

	strcpy( line, "node id off" );
	fout.write( line, 80 );

	strcpy( line, "element id off" );
	fout.write( line, 80 );

	strcpy( line, "extents" );
	fout.write( line, 80 );

	fout.write( (char*) &minX, sizeof( float ) );
	fout.write( (char*) &maxX, sizeof( float ) );
	fout.write( (char*) &minY, sizeof( float ) );
	fout.write( (char*) &maxY, sizeof( float ) );
	fout.write( (char*) &minZ, sizeof( float ) );
	fout.write( (char*) &maxZ, sizeof( float ) );

	strcpy( line, "part" );
	fout.write( line, 80 );

	int ibuf(1);
	fout.write( (char*) &ibuf, sizeof( int ) );

	strcpy( line, "Element data" );
	fout.write( line, 80 );

	strcpy( line, "coordinates" );
	fout.write( line, 80 );

	fout.write( (char*) &m_numNodeTotal, sizeof( int ) );

	fout.write( (char*) coordX, sizeof( float )*m_numNodeTotal );
	fout.write( (char*) coordY, sizeof( float )*m_numNodeTotal );
	fout.write( (char*) coordZ, sizeof( float )*m_numNodeTotal );

	strcpy( line, "hexa8" );
	fout.write( line, 80 );

	fout.write( (char*) &m_numElemTotal, sizeof( int ) );

	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
		for( int i = 0; i < 8; ++i ){
			int ibuf = m_nodesOfElements[iElem*8 + i] + 1;
			fout.write( (char*) &ibuf, sizeof( int ) );
		}
	}

	delete [] coordX;
	delete [] coordY;
	delete [] coordZ;

	fout.close();

}


// Get local coordinate values from coordinate values
// <Input> : iElem, coordX, coordY, coordZ
// <Output>: localCoordX, localCoordY, localCoordZ
void MeshDataBrickElement::getLocalCoordinateValues( const int iElem, const double coordX, const double coordY, const double coordZ,
	double& localCoordX, double& localCoordY, double& localCoordZ ) const{

	const int node0 = m_nodesOfElements[ iElem * 8     ];
	const int node6 = m_nodesOfElements[ iElem * 8 + 6 ];
	coordinateValue coordMin = { m_xCoordinatesOfNodes[node0], m_yCoordinatesOfNodes[node0], m_zCoordinatesOfNodes[node0] };
	coordinateValue coordMax = { m_xCoordinatesOfNodes[node6], m_yCoordinatesOfNodes[node6], m_zCoordinatesOfNodes[node6] };

	localCoordX = 2 * ( coordX - ( coordMin.X + coordMax.X ) * 0.5 ) / std::fabs( coordMax.X - coordMin.X );
	if( localCoordX > 1.0 ){
		localCoordX = 1.0;
	}
	if( localCoordX < -1.0 ){
		localCoordX = -1.0;
	}

	localCoordY = 2 * ( coordY - ( coordMin.Y + coordMax.Y ) * 0.5 ) / std::fabs( coordMax.Y - coordMin.Y );
	if( localCoordY > 1.0 ){
		localCoordY = 1.0;
	}
	if( localCoordY < -1.0 ){
		localCoordY = -1.0;
	}

	localCoordZ = 2 * ( coordZ - ( coordMin.Z + coordMax.Z ) * 0.5 ) / std::fabs( coordMax.Z - coordMin.Z );
	if( localCoordZ > 1.0 ){
		localCoordZ = 1.0;
	}
	if( localCoordZ < -1.0 ){
		localCoordZ = -1.0;
	}

}
