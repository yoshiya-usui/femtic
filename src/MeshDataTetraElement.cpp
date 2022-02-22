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
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "AnalysisControl.h"
#include "MeshDataTetraElement.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "Util.h"

const double MeshDataTetraElement::m_eps = 1.0e-12;

// Constructer
MeshDataTetraElement::MeshDataTetraElement():
	m_numElemOnLandSurface(0),
	m_elemOnLandSurface(NULL),
	m_faceLandSurface(NULL)
{
	m_numNodeOneElement = 4;
	m_numEdgeOneElement = 6;
	m_numNodeOnFaceOneElement = 3;
	m_numNeighborElement = 4;

	for ( int i = 0; i < 6; ++i ){
		m_facesOfElementsBoundaryPlanes[i] = NULL;
	}

	// m_faceID2NodeID
	m_faceID2NodeID[0][0] = 1;
	m_faceID2NodeID[0][1] = 2;
	m_faceID2NodeID[0][2] = 3;

	m_faceID2NodeID[1][0] = 0;
	m_faceID2NodeID[1][1] = 3;
	m_faceID2NodeID[1][2] = 2;

	m_faceID2NodeID[2][0] = 0;
	m_faceID2NodeID[2][1] = 1;
	m_faceID2NodeID[2][2] = 3;

	m_faceID2NodeID[3][0] = 0;
	m_faceID2NodeID[3][1] = 2;
	m_faceID2NodeID[3][2] = 1;

	// m_faceID2EdgeID
	m_faceID2EdgeID[0][0] = 3;
	m_faceID2EdgeID[0][1] = 5;
	m_faceID2EdgeID[0][2] = 4;

	m_faceID2EdgeID[1][0] = 2;
	m_faceID2EdgeID[1][1] = 5;
	m_faceID2EdgeID[1][2] = 1;

	m_faceID2EdgeID[2][0] = 0;
	m_faceID2EdgeID[2][1] = 4;
	m_faceID2EdgeID[2][2] = 2;

	m_faceID2EdgeID[3][0] = 1;
	m_faceID2EdgeID[3][1] = 3;
	m_faceID2EdgeID[3][2] = 0;

	// m_edgeID2NodeID
	m_edgeID2NodeID[0][0] = 0;
	m_edgeID2NodeID[0][1] = 1;

	m_edgeID2NodeID[1][0] = 0;
	m_edgeID2NodeID[1][1] = 2;

	m_edgeID2NodeID[2][0] = 0;
	m_edgeID2NodeID[2][1] = 3;

	m_edgeID2NodeID[3][0] = 1;
	m_edgeID2NodeID[3][1] = 2;

	m_edgeID2NodeID[4][0] = 3;
	m_edgeID2NodeID[4][1] = 1;

	m_edgeID2NodeID[5][0] = 2;
	m_edgeID2NodeID[5][1] = 3;

}

// Destructer
MeshDataTetraElement::~MeshDataTetraElement(){

	for ( int i = 0; i < 6; ++i ){
		if( m_elemBoundaryPlanes[i] != NULL ){
			delete[] m_elemBoundaryPlanes[i];
			m_elemBoundaryPlanes[i] = NULL;
		}
	}

	for ( int i = 0; i < 6; ++i ){
		if( m_facesOfElementsBoundaryPlanes[i] != NULL ){
			delete[] m_facesOfElementsBoundaryPlanes[i];
			m_facesOfElementsBoundaryPlanes[i] = NULL;
		}
	}

	if( m_elemOnLandSurface != NULL ){
		delete[] m_elemOnLandSurface;
		m_elemOnLandSurface = NULL;
	}

	if( m_faceLandSurface != NULL ){
		delete[] m_faceLandSurface;
		m_faceLandSurface = NULL;
	}

}

// Input mesh data from "mesh.dat"
void MeshDataTetraElement::inputMeshData(){

	std::ifstream inFile( "mesh.dat", std::ios::in );
	if( inFile.fail() )
	{
		OutputFiles::m_logFile << "File open error : mesh.dat !!" << std::endl;
		exit(1);
	}

	std::string sbuf;
	inFile >> sbuf;

	if( sbuf.substr(0,5).compare("TETRA") != 0 ){
		OutputFiles::m_logFile << "Mesh data written in mesh.dat is different from the ones of tetrahedral element !!" << std::endl;
		exit(1);
	}

	int ibuf(0);

	inFile >> ibuf;
	if( ibuf > 0 ){
		m_numNodeTotal = ibuf;
	}else{
		OutputFiles::m_logFile << "Total number of nodes is less than or equal to zero ! : " << ibuf << std::endl;
		exit(1);
	}

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

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){

		int idum(0);
		inFile >> idum >> m_xCoordinatesOfNodes[iNode] >> m_yCoordinatesOfNodes[iNode] >> m_zCoordinatesOfNodes[iNode];

#ifdef _DEBUG_WRITE
		std::cout << idum << " " << m_xCoordinatesOfNodes[iNode] << " " << m_yCoordinatesOfNodes[iNode] << " " << m_zCoordinatesOfNodes[iNode] << std::endl; // For debug
#endif

	}

	inFile >> ibuf;
	if( ibuf > 0 ){
		m_numElemTotal = ibuf;
	}else{
		OutputFiles::m_logFile << "Total number of elements is less than or equal to zero ! : " << ibuf << std::endl;
		exit(1);
	}

	if( m_neighborElements != NULL ){
		delete[] m_neighborElements;
	}
	m_neighborElements = new int[ m_numElemTotal * 4 ];

	if( m_nodesOfElements == NULL ){
		delete[] m_nodesOfElements;
	}
	m_nodesOfElements = new int[ m_numElemTotal * m_numNodeOneElement ];

	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){

		int idum(0);
		inFile >> idum;

		// IDs of neighbor Elements
		for( int i = 0; i < 4; ++i ){
			int neib(0);
			inFile >> neib;

			m_neighborElements[ iElem * 4 + i ] = neib;

		}

		// Nodes of the element
		for( int i = 0; i < m_numNodeOneElement; ++i ){
			int node(0);
			inFile >> node;

			m_nodesOfElements[ iElem * m_numNodeOneElement + i ] = node;

		}

	}

#ifdef _DEBUG_WRITE
	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){

		std::cout << iElem << " ";

		// IDs of neighbor Elements
		for( int i = 0; i < 4; ++i ){
			std::cout << m_neighborElements[ iElem * 4 + i ] << " ";
		}

		// Nodes of the element
		for( int i = 0; i < m_numNodeOneElement; ++i ){
			std::cout << m_nodesOfElements[ iElem * m_numNodeOneElement + i ] << " ";
		}

		std::cout << std::endl;
	}
#endif


	for( int iPlane = 0; iPlane < 6; ++iPlane ){// Loop of boundary planes

		int nElemOnPlane;
		inFile >> nElemOnPlane;
		if( nElemOnPlane > 0 ){
			m_numElemOnBoundaryPlanes[iPlane] = nElemOnPlane;
		}else{
			OutputFiles::m_logFile << "Number of faces belonging plane " << iPlane << " is less than or equal to zero ! : " << nElemOnPlane << std::endl;
			exit(1);
		}

#ifdef _DEBUG_WRITE
		std::cout << nElemOnPlane << std::endl; // For debug
#endif

		if( m_elemBoundaryPlanes[iPlane] != NULL ){
			delete [] m_elemBoundaryPlanes[iPlane];
		}
		m_elemBoundaryPlanes[iPlane] = new int[ nElemOnPlane ];

		if( m_facesOfElementsBoundaryPlanes[iPlane] != NULL ){
			delete [] m_facesOfElementsBoundaryPlanes[iPlane];	
		}
		m_facesOfElementsBoundaryPlanes[iPlane] = new int[ nElemOnPlane ];

		// Set elements belonging to the boundary planes
		for( int iElem = 0; iElem < nElemOnPlane; ++iElem ){		

			inFile >> m_elemBoundaryPlanes[iPlane][iElem] >> m_facesOfElementsBoundaryPlanes[iPlane][iElem];

			if( m_elemBoundaryPlanes[iPlane][iElem] < 0 || m_elemBoundaryPlanes[iPlane][iElem] >= m_numElemTotal ){
				OutputFiles::m_logFile << "Element ID of plane " << iPlane << " is out of range !! : " << m_elemBoundaryPlanes[iPlane][iElem] << std::endl;
				exit(1);
			}
			if( m_facesOfElementsBoundaryPlanes[iPlane][iElem] < 0 || m_facesOfElementsBoundaryPlanes[iPlane][iElem] >= 4 ){
				OutputFiles::m_logFile << "Face ID of plane " << iPlane << " is out of range !! : " << m_facesOfElementsBoundaryPlanes[iPlane][iElem] << std::endl;
				exit(1);
			}
			
#ifdef _DEBUG_WRITE
			std::cout << m_elemBoundaryPlanes[iPlane][iElem] << " "  << m_facesOfElementsBoundaryPlanes[iPlane][iElem] << std::endl; // For debug
#endif

		}
		
	}

	inFile >> ibuf;
	if( ibuf > 0 ){
		m_numElemOnLandSurface = ibuf;
	}else{
		OutputFiles::m_logFile << "Total number of faces on the land surface is less than or equal to zero ! : " << ibuf << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
		std::cout << "m_numElemOnLandSurface = " << m_numElemOnLandSurface << std::endl; // For debug
#endif

	if( m_elemOnLandSurface != NULL ){
		delete [] m_elemOnLandSurface;	
	}
	m_elemOnLandSurface = new int[ m_numElemOnLandSurface ];

	if( m_faceLandSurface != NULL ){
		delete [] m_faceLandSurface;	
	}
	m_faceLandSurface = new int[ m_numElemOnLandSurface ];

	// Set faces belonging to the boundary planes
	for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){		

		inFile >> m_elemOnLandSurface[iElem] >> m_faceLandSurface[iElem];

		if( m_elemOnLandSurface[iElem] < 0 || m_elemOnLandSurface[iElem] >= m_numElemTotal ){
			OutputFiles::m_logFile << "Element ID of land surface is out of range !! : " << m_elemOnLandSurface[iElem] << std::endl;
			exit(1);
		}
		if( m_faceLandSurface[iElem] < 0 || m_faceLandSurface[iElem] >= 4 ){
			OutputFiles::m_logFile << "Face ID of land surface is out of range !! : " << m_faceLandSurface[iElem] << std::endl;
			exit(1);
		}
			
#ifdef _DEBUG_WRITE
		std::cout << m_elemOnLandSurface[iElem] << " "  << m_faceLandSurface[iElem] << std::endl; // For debug
#endif

	}

	inFile.close();

}

// Find element including a point
int MeshDataTetraElement::findElementIncludingPoint( const double locX, const double locY, const double locZ, CommonParameters::VolumeCoords& localCoord ) const{

	const CommonParameters::locationXYZ loc = { locX, locY, locZ };

	int elemID = m_elemOnLandSurface[0];

	for( int i = 0; i < m_numElemTotal; ++i ){

		bool found(true);
		for( int iFace = 0; iFace < 4; ++iFace ){

			if( !locateInsideOfFace( elemID, iFace, loc ) ){
				const int neibElemID = m_neighborElements[ 4 * elemID + iFace ];
				found = false;
				if( neibElemID < 0 ){
					continue;
				}
				elemID = neibElemID;
				break;
			}
		}

		if( found ){
			calcVolumeCoordsOfPoint( elemID, loc, localCoord );
			return elemID;
		}

	}

	OutputFiles::m_logFile << " Error : Could not find element including point ( X = " << locX << "[m], Y= " << locY << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	exit(1);

	return -1;

}

int MeshDataTetraElement::findElementIncludingPointOnSurface( const double locX, const double locY, int& faceID, CommonParameters::AreaCoords& localCoord,
	const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod ) const{

	const CommonParameters::locationXY pointCoord = { locX, locY };

	for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){

		int elemID = m_elemOnLandSurface[iElem];
		faceID = m_faceLandSurface[iElem];

		if( useUpperElem ){
			// Find the upper element of the site from the lower element
			const int neibElemID = m_neighborElements[ 4 * elemID + faceID ];
			int neibFaceID = -1;
			for( int iFace = 0; iFace < 4; ++iFace ){
				if( m_neighborElements[ 4 * neibElemID + iFace ] == elemID ){
					neibFaceID = iFace;
					break;
				}
			}
			if( neibFaceID < 0 ){
				OutputFiles::m_logFile << " Error : There is no adjacent element for the Face " << faceID << " of the Element " << elemID << " ." << std::endl;
				exit(1);
			}
			else{
				elemID = neibElemID;
				faceID = neibFaceID;
			}
		}

		const int nodeID0 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][0] );
		const int nodeID1 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][1] );
		const int nodeID2 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][2] );

		const CommonParameters::locationXY nodeCoord0 = { getXCoordinatesOfNodes( nodeID0 ), getYCoordinatesOfNodes( nodeID0 ) };
		const CommonParameters::locationXY nodeCoord1 = { getXCoordinatesOfNodes( nodeID1 ), getYCoordinatesOfNodes( nodeID1 ) };
		const CommonParameters::locationXY nodeCoord2 = { getXCoordinatesOfNodes( nodeID2 ), getYCoordinatesOfNodes( nodeID2 ) };

#ifdef _DEBUG_WRITE
		std::cout << "elemID : " << elemID << std::endl;
		std::cout << "faceID : " << faceID << std::endl;
		std::cout << "locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord0, nodeCoord1 ) : " << locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord0, nodeCoord1 ) << std::endl;
		std::cout << "locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord1, nodeCoord2 ) : " << locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord1, nodeCoord2 ) << std::endl;
		std::cout << "locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord2, nodeCoord0 ) : " << locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord2, nodeCoord0 ) << std::endl;
#endif

		bool locateInElement(false);
		if( useUpperElem ){
			// For upper element
			if( locateLeftOfSegmentOnSeaSurface( pointCoord, nodeCoord0, nodeCoord1 ) &&
				locateLeftOfSegmentOnSeaSurface( pointCoord, nodeCoord1, nodeCoord2 ) &&
				locateLeftOfSegmentOnSeaSurface( pointCoord, nodeCoord2, nodeCoord0 ) ){
				locateInElement = true;
			}
		}
		else{
			// For lower element
			if( locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord0, nodeCoord1 ) &&
				locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord1, nodeCoord2 ) &&
				locateLeftOfSegmentOnLandSurface( pointCoord, nodeCoord2, nodeCoord0 ) ){
				locateInElement = true;
			}
		}
		if( locateInElement ){
			if( modLoc ){
				// Modify the location to the center of the face
				locXMod = ( nodeCoord0.X + nodeCoord1.X + nodeCoord2.X ) / 3.0;
				locYMod = ( nodeCoord0.Y + nodeCoord1.Y + nodeCoord2.Y ) / 3.0;
				const CommonParameters::locationXY pointCoordMod = { locXMod, locYMod };
				calcAreaCoordsOfPointOnLandSurface( elemID, faceID, pointCoordMod, localCoord );
			}
			else{
				calcAreaCoordsOfPointOnLandSurface( elemID, faceID, pointCoord, localCoord );
			}
			return elemID;
		}
	}

	OutputFiles::m_logFile << " Error : Could not find element including point ( X = " << locX << "[m], Y= " << locY << "[m], Z= " << 0.0 <<  "[m] )." << std::endl;
	exit(1);
	return -1;

}

// Find element including a point on the Y-Z plane and return element ID of 2D mesh
int MeshDataTetraElement::findElementIncludingPointOnYZPlaneAndReturnElemID2D( const int iPlane, const double locY, const double locZ, CommonParameters::AreaCoords& localCoord ) const{

	if( iPlane != MeshData::YZMinus && iPlane != MeshData::YZPlus ){
		OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
		exit(1);
	}

	const CommonParameters::locationYZ pointCoord = { locY, locZ };

	const int nElem = m_numElemOnBoundaryPlanes[iPlane];

	for( int iElem = 0; iElem < nElem; ++iElem ){

		//const CommonParameters::locationYZ nodeCoord0 = { getCoordYFromElementBoundaryPlanes( iPlane, iElem, 0 ), getCoordZFromElementBoundaryPlanes( iPlane, iElem, 0 ) };
		//const CommonParameters::locationYZ nodeCoord1 = { getCoordYFromElementBoundaryPlanes( iPlane, iElem, 1 ), getCoordZFromElementBoundaryPlanes( iPlane, iElem, 1 ) };
		//const CommonParameters::locationYZ nodeCoord2 = { getCoordYFromElementBoundaryPlanes( iPlane, iElem, 2 ), getCoordZFromElementBoundaryPlanes( iPlane, iElem, 2 ) };

		if( locateLeftOfSegmentOnYZPlaneOfBoundary( iPlane, iElem, 0, pointCoord ) &&
			locateLeftOfSegmentOnYZPlaneOfBoundary( iPlane, iElem, 1, pointCoord ) &&
			locateLeftOfSegmentOnYZPlaneOfBoundary( iPlane, iElem, 2, pointCoord ) ){

			const CommonParameters::CoordPair point = { pointCoord.Y, pointCoord.Z };

			calcAreaCoordsOfPointOnYZPlaneOfBoundary( iPlane, iElem, point, localCoord );

			return iElem;
		}

	}

	OutputFiles::m_logFile << " Warning : Could not find element including point ( X = " << 0.0 << "[m], Y= " << locY << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	return -1;

}

// Find element including a point on the Z-X plane and return element ID of 2D mesh
int MeshDataTetraElement::findElementIncludingPointOnZXPlaneAndReturnElemID2D( const int iPlane, const double locZ, const double locX, CommonParameters::AreaCoords& localCoord ) const{

	if( iPlane != MeshData::ZXMinus && iPlane != MeshData::ZXPlus ){
		OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
		exit(1);
	}

	const CommonParameters::locationZX pointCoord = { locZ, locX };

	const int nElem = m_numElemOnBoundaryPlanes[iPlane];

	for( int iElem = 0; iElem < nElem; ++iElem ){

		//const CommonParameters::locationZX nodeCoord0 = { getCoordZFromElementBoundaryPlanes( iPlane, iElem, 0 ), getCoordXFromElementBoundaryPlanes( iPlane, iElem, 0 ) };
		//const CommonParameters::locationZX nodeCoord1 = { getCoordZFromElementBoundaryPlanes( iPlane, iElem, 1 ), getCoordXFromElementBoundaryPlanes( iPlane, iElem, 1 ) };
		//const CommonParameters::locationZX nodeCoord2 = { getCoordZFromElementBoundaryPlanes( iPlane, iElem, 2 ), getCoordXFromElementBoundaryPlanes( iPlane, iElem, 2 ) };

		if( locateLeftOfSegmentOnZXPlaneOfBoundary( iPlane, iElem, 0, pointCoord ) &&
			locateLeftOfSegmentOnZXPlaneOfBoundary( iPlane, iElem, 1, pointCoord ) &&
			locateLeftOfSegmentOnZXPlaneOfBoundary( iPlane, iElem, 2, pointCoord ) ){

			const CommonParameters::CoordPair point = { pointCoord.Z, pointCoord.X };

			calcAreaCoordsOfPointOnZXPlaneOfBoundary( iPlane, iElem, point, localCoord );

			return iElem;
		}

	}

	OutputFiles::m_logFile << " Warning : Could not find element including point ( X = " << locX << "[m], Y= " << 0.0 << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	return -1;

}

// Find elements including dipole on the surface of the earth
void MeshDataTetraElement::findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd, std::vector<int>& elements, std::vector<int>& faces,
	std::vector<CommonParameters::AreaCoords>& areaCoordsdStartPoint, std::vector<CommonParameters::AreaCoords>& areaCoordsEndPoint ) const{

	const double thresholdVal = 1.0E-6;

	const CommonParameters::locationXY nodeCoordDipoleStart = { locXStart, locYStart };
	const CommonParameters::locationXY nodeCoordDipoleEnd = { locXEnd, locYEnd };

	const double dipoleLength = hypot( ( locXEnd - locXStart ), ( locYEnd - locYStart ) );
	if( dipoleLength < thresholdVal ){
		OutputFiles::m_logFile << " Error : Length of dipole ( " << dipoleLength << " [m] ) is too small !! " << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	std::cout << "dipoleLength =  " << dipoleLength << std::endl; // For debug
#endif
	
	double dummy(0.0);

	int faceIDStartPoint = 0;
	CommonParameters::AreaCoords areaCoordStart = { 0.0, 0.0, 0.0 }; 
	const int iElemStart = findElementIncludingPointOnSurface( locXStart, locYStart, faceIDStartPoint, areaCoordStart, false, false, dummy, dummy );

	int faceIDEndPoint = 0;
	CommonParameters::AreaCoords areaCoordEnd = { 0.0, 0.0, 0.0 };
	const int iElemEnd = findElementIncludingPointOnSurface( locXEnd, locYEnd, faceIDEndPoint, areaCoordEnd, false, false, dummy, dummy );

#ifdef _DEBUG_WRITE
	std::cout << "iElemStart iElemEnd : " << iElemStart << " " << iElemEnd << std::endl; // For debug
#endif
	
	std::vector<CommonParameters::locationXY> intersectPointsAreadyFound[2];

	for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){

		const int elemID = m_elemOnLandSurface[iElem];
		const int faceID = m_faceLandSurface[iElem];
		
#ifdef _DEBUG_WRITE
		std::cout << "elemID : " << elemID << " " << std::endl;
		std::cout << "faceID : " << faceID << " " << std::endl;
#endif

		const int nodeID0 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][0] );
		const int nodeID1 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][1] );
		const int nodeID2 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][2] );

		const CommonParameters::locationXY nodeCoord0 = { getXCoordinatesOfNodes( nodeID0 ), getYCoordinatesOfNodes( nodeID0 ) };
		const CommonParameters::locationXY nodeCoord1 = { getXCoordinatesOfNodes( nodeID1 ), getYCoordinatesOfNodes( nodeID1 ) };
		const CommonParameters::locationXY nodeCoord2 = { getXCoordinatesOfNodes( nodeID2 ), getYCoordinatesOfNodes( nodeID2 ) };

		//CommonParameters::locationXY intersectPoints[2];
		CommonParameters::locationXY intersectPoints[3];

		if( overlapTwoSegments( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){

#ifdef _DEBUG_WRITE	
			std::cout << "A" << std::endl;
#endif

			const double innerProduct1 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord0 ) / dipoleLength;
			const double innerProduct2 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord1 ) / dipoleLength;

			if( innerProduct1 < innerProduct2 ){
				if( innerProduct1 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord0;
				}
				if( innerProduct2 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord1;
				}
			}else{
				if( innerProduct2 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord1;
				}
				if( innerProduct1 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord0;
				}
			}

		}else if( overlapTwoSegments( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){

#ifdef _DEBUG_WRITE	
			std::cout << "B" << std::endl;
#endif

			const double innerProduct1 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord1 ) / dipoleLength;
			const double innerProduct2 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord2 ) / dipoleLength;

			if( innerProduct1 < innerProduct2 ){
				if( innerProduct1 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord1;
				}
				if( innerProduct2 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord2;
				}
			}else{
				if( innerProduct2 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord2;
				}
				if( innerProduct1 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord1;
				}
			}

		}else if( overlapTwoSegments( nodeCoord2, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){

#ifdef _DEBUG_WRITE	
			std::cout << "C" << std::endl;
#endif

			const double innerProduct1 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord2 ) / dipoleLength;
			const double innerProduct2 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord0 ) / dipoleLength;

			if( innerProduct1 < innerProduct2 ){
				if( innerProduct1 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord2;
				}
				if( innerProduct2 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord0;
				}
			}else{
				if( innerProduct2 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord0;
				}
				if( innerProduct1 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord2;
				}
			}

		}else{

#ifdef _DEBUG_WRITE	
			std::cout << "D" << std::endl;
#endif

			int icount(0);

			if( intersectTwoSegments( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				++icount;
			}

			if( intersectTwoSegments( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				++icount;
			}

			//if( icount < 2 && intersectTwoSegments( nodeCoord2, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
			//	calcCoordOfIntersectionPointOfTwoLines( nodeCoord2, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
			//	++icount;
			//}
			if( intersectTwoSegments( nodeCoord2, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord2, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				++icount;
			}

#ifdef _DEBUG_WRITE	
			std::cout << "icount : " << icount << std::endl;
#endif

			if( elemID == iElemStart && elemID == iElemEnd ){
				intersectPoints[0] = nodeCoordDipoleStart;
				intersectPoints[1] = nodeCoordDipoleEnd;
			}
			else if( icount == 0 ){
				continue;
			}else if( icount == 1 ){

				if( elemID == iElemStart ){
					intersectPoints[1] = nodeCoordDipoleStart;
				}else if( elemID == iElemEnd ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					//OutputFiles::m_logFile << " Error : Number of intersection points of dipole and edges of element " << elemID << " is only one !! " << std::endl;
					//exit(1);
					continue;
				}

			}else if( icount == 3 ){

				double max = calcDistance( intersectPoints[0], intersectPoints[1] );
				int iStart(0);
				int iEnd(1);

				if( calcDistance( intersectPoints[1], intersectPoints[2] ) > max ){
					iStart = 1;
					iEnd = 2;
					max = calcDistance( intersectPoints[iStart], intersectPoints[iEnd] );
				}

				if( calcDistance( intersectPoints[2], intersectPoints[0] ) > max ){
					iStart = 2;
					iEnd = 0;
				}

				if( iStart != 0 || iEnd != 1  ){
					CommonParameters::locationXY tmpStart = { intersectPoints[iStart].X, intersectPoints[iStart].Y };
					CommonParameters::locationXY tmpEnd   = { intersectPoints[iEnd].X,   intersectPoints[iEnd].Y   };
					intersectPoints[0].X = tmpStart.X;
					intersectPoints[0].Y = tmpStart.Y;
					intersectPoints[1].X = tmpEnd.X;
					intersectPoints[1].Y = tmpEnd.Y;
				}

#ifdef _DEBUG_WRITE	
				std::cout << "iStart iEnd : " << iStart << " " << iEnd << std::endl;
				std::cout << "max : " << max << std::endl;
#endif

			}
			
		}

#ifdef _DEBUG_WRITE	
		std::cout << "intersectPoints : " << intersectPoints[0].X << " " << intersectPoints[0].Y << std::endl;
		std::cout << "intersectPoints : " << intersectPoints[1].X << " " << intersectPoints[1].Y << std::endl;
#endif

		const double innerProduct = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[0], intersectPoints[1] );

		if( fabs(innerProduct) < thresholdVal ){// Segment is too short
#ifdef _DEBUG_WRITE	
			std::cout << "Too short" << std::endl;
#endif
			continue;
		}

		if( innerProduct < 0 ){
			const CommonParameters::locationXY temp = intersectPoints[0];
			intersectPoints[0] = intersectPoints[1];
			intersectPoints[1] = temp;
		}

		bool alreadyFound[2] = { false, false };
		for( int i = 0; i < 2; ++i ){
			std::vector<CommonParameters::locationXY>::iterator itrEnd = intersectPointsAreadyFound[i].end();
			for( std::vector<CommonParameters::locationXY>::iterator itr = intersectPointsAreadyFound[i].begin(); itr != itrEnd; ++itr ){
				//const double eps = 1.0e-8;
				if( fabs( (itr->X) - intersectPoints[i].X ) < m_eps && fabs( (itr->Y) - intersectPoints[i].Y ) < m_eps ){
					alreadyFound[i] = true;
				}
			}
		}
		
		if( alreadyFound[0] && alreadyFound[1] ){// Segment has already been found
#ifdef _DEBUG_WRITE	
			std::cout << "Already Found" << std::endl;
#endif
			continue;
		}

		CommonParameters::AreaCoords coords;

		calcAreaCoordsOfPointOnLandSurface( elemID, faceID, intersectPoints[0], coords );
		areaCoordsdStartPoint.push_back( coords );

#ifdef _DEBUG_WRITE	
		std::cout << "start coords : " << coords.coord0 << " " << coords.coord1 << " " << coords.coord2 << std::endl;
#endif

		calcAreaCoordsOfPointOnLandSurface( elemID, faceID, intersectPoints[1], coords );
		areaCoordsEndPoint.push_back( coords );

#ifdef _DEBUG_WRITE	
		std::cout << "end coords : " << coords.coord0 << " " << coords.coord1 << " " << coords.coord2 << std::endl;
#endif

		elements.push_back( elemID );
		faces.push_back( faceID );

#ifdef _DEBUG_WRITE	
		std::cout << "intersectPoints[0] : " << intersectPoints[0].X << " " << intersectPoints[0].Y << " " << std::endl;
		std::cout << "intersectPoints[1] : " << intersectPoints[1].X << " " << intersectPoints[1].Y << " " << std::endl;
#endif

		intersectPointsAreadyFound[0].push_back( intersectPoints[0] );
		intersectPointsAreadyFound[1].push_back( intersectPoints[1] );

	}

	if( elements.empty() ){
		OutputFiles::m_logFile << " Error : Could not find element including dipole ( " << locXStart << ", " << locYStart << " ) => ( " << locXEnd << ", " << locYEnd << " )." << std::endl;
		exit(1);
	}

	double dipoleLengthAccumulated(0.0);
	std::vector<CommonParameters::locationXY>::iterator itr0End = intersectPointsAreadyFound[0].end(); 
	std::vector<CommonParameters::locationXY>::iterator itr0  = intersectPointsAreadyFound[0].begin(); 
	std::vector<CommonParameters::locationXY>::iterator itr1  = intersectPointsAreadyFound[1].begin(); 
	while( itr0 != itr0End ){
		//dipoleLengthAccumulated += hypot( itr1->X - itr0->X, itr1->Y - itr0->Y );
		dipoleLengthAccumulated += calcDistance( *itr0, *itr1 );
		++itr0;
		++itr1;
	}

	if( fabs(dipoleLength - dipoleLengthAccumulated) / fabs(dipoleLength) > 0.01 ){
		OutputFiles::m_logFile << " Warning : Accumulated dipole length (" << dipoleLengthAccumulated
			<< ") is significantly different from the dipole length of the horizontal plane (" << dipoleLength
			<< ")" << std::endl;
	}

}

// Decide whether specified elements share same edges
bool MeshDataTetraElement::shareSameEdges( const int elemID1, const int elemID2 ) const{

	assert( elemID1 >= 0 );
	assert( elemID1 < m_numElemTotal );
	assert( elemID2 >= 0 );
	assert( elemID2 < m_numElemTotal );

	const int iEdgeToNode[6][2] = { {0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3} }; 

	for( int iEdge1 = 0; iEdge1 < 6; ++iEdge1 ){

		const int nodeID1[2] = { getNodesOfElements( elemID1, iEdgeToNode[iEdge1][0] ), getNodesOfElements( elemID1, iEdgeToNode[iEdge1][1] ) };

		for( int iEdge2 = 0; iEdge2 < 6; ++iEdge2 ){

			const int nodeID2[2] = { getNodesOfElements( elemID2, iEdgeToNode[iEdge2][0] ), getNodesOfElements( elemID2, iEdgeToNode[iEdge2][1] ) };

			if( ( nodeID1[0] == nodeID2[0] && nodeID1[1] == nodeID2[1] ) ||
				( nodeID1[0] == nodeID2[1] && nodeID1[1] == nodeID2[0] ) ){
				return true;
			}

		}

	}

	return false;

}

// Calculate volume of a specified element
double MeshDataTetraElement::calcVolume( const int elemID ) const{
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	CommonParameters::locationXYZ nodeCoord[4];
	for( int i = 0; i < 4; ++i ){
		const int nodeID = getNodesOfElements( elemID, i ); 
		nodeCoord[i].X = getXCoordinatesOfNodes( nodeID );
		nodeCoord[i].Y = getYCoordinatesOfNodes( nodeID );
		nodeCoord[i].Z = getZCoordinatesOfNodes( nodeID );
	}

	return calcVolume( nodeCoord[0], nodeCoord[1], nodeCoord[2], nodeCoord[3] );
}

// Output vtk file
void MeshDataTetraElement::outputMeshDataToVTK() const{

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

	OutputFiles::m_vtkFile << "CELLS " << m_numElemTotal << " " << m_numElemTotal*5 << std::endl;
	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
		OutputFiles::m_vtkFile << m_numNodeOneElement << " "
			<< m_nodesOfElements[iElem*m_numNodeOneElement    ] << " " 
			<< m_nodesOfElements[iElem*m_numNodeOneElement + 1] << " " 
			<< m_nodesOfElements[iElem*m_numNodeOneElement + 2] << " " 
			<< m_nodesOfElements[iElem*m_numNodeOneElement + 3] << std::endl;
	}

	OutputFiles::m_vtkFile << "CELL_TYPES " << m_numElemTotal << std::endl;
	for( int iElem = 0 ; iElem < m_numElemTotal; ++iElem ){
		OutputFiles::m_vtkFile << "10" << std::endl;
	}

	OutputFiles::m_vtkFile << "CELL_DATA " << m_numElemTotal << std::endl;

}

// Output mesh data to binary file
void MeshDataTetraElement::outputMeshDataToBinary() const{
	
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

	strcpy( line, "Tetrahedral Mesh" );
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

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		float dbuf = static_cast<float>( m_xCoordinatesOfNodes[iNode] );
		fout.write( (char*) &dbuf, sizeof( float ) );
	}

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		float dbuf = static_cast<float>( m_yCoordinatesOfNodes[iNode] );
		fout.write( (char*) &dbuf, sizeof( float ) );
	}

	for( int iNode = 0; iNode < m_numNodeTotal; ++iNode ){
		float dbuf = static_cast<float>( m_zCoordinatesOfNodes[iNode] );
		fout.write( (char*) &dbuf, sizeof( float ) );
	}

	strcpy( line, "tetra4" );
	fout.write( line, 80 );

	fout.write( (char*) &m_numElemTotal, sizeof( int ) );

	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
		for( int i = 0; i < 4; ++i ){
			int ibuf = m_nodesOfElements[iElem*4 + i] + 1;
			fout.write( (char*) &ibuf, sizeof( int ) );
		}
	}

	delete [] coordX;
	delete [] coordY;
	delete [] coordZ;

	fout.close();

}

// Get array of nodes of elements belonging to the boundary planes
int MeshDataTetraElement::getNodesOfElementsBoundaryPlanes( const int iPlane, const int iElem, const int iNode ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getNodesOfElements. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iPlane < 0 || iPlane >= 6 ){
	//	OutputFiles::m_logFile << " Error : iPlane is out of range in getNodesOfElementsBoundaryPlanes !! : iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}

	//if( iNode < 0 || iNode >= m_numNodeOneElement ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getNodesOfElements. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iPlane >= 0 );
	assert( iPlane < 6 );
	assert( iNode >= 0 );
	assert( iNode < m_numNodeOneElement );

	const int elemID = m_elemBoundaryPlanes[iPlane][iElem];

	const int faceID = m_facesOfElementsBoundaryPlanes[iPlane][iElem];

	return getNodesOfElements( elemID, m_faceID2NodeID[faceID][iNode] );

}

// Get mesh type
int MeshDataTetraElement::getMeshType() const{

	return MeshData::TETRA;

}

// Get local face ID of elements belonging to the boundary planes
int MeshDataTetraElement::getFaceIDLocalFromElementBoundaryPlanes( const int iPlane, const int iElem ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getFaceIDLocalFromElementBoundaryPlanes. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iPlane < 0 || iPlane >= 6 ){
	//	OutputFiles::m_logFile << " Error : iPlane is out of range in getFaceIDLocalFromElementBoundaryPlanes !! : iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iPlane >= 0 );
	assert( iPlane < 6 );

	return m_facesOfElementsBoundaryPlanes[iPlane][iElem];

}

// Get local node ID  from local face ID
int MeshDataTetraElement::getNodeIDLocalFromFaceIDLocal( const int iFace, const int num ) const{

	//if( iFace < 0 || iFace >= 4 ){
	//	OutputFiles::m_logFile << " Error : iFace is out of range in getNodeIDLocalFromFaceIDLocal. iFace = " << iFace << std::endl;
	//	exit(1);
	//}

	//if( num < 0 || num >= 3 ){
	//	OutputFiles::m_logFile << " Error : num is out of range in getNodeIDLocalFromFaceIDLocal. num = " << num << std::endl;
	//	exit(1);
	//}
	assert( iFace >= 0 );
	assert( iFace < 4 );
	assert( num >= 0 );
	assert( num < 3 );

	return m_faceID2NodeID[iFace][num];

}

// Get local node ID from local edge ID
int MeshDataTetraElement::getNodeIDLocalFromEdgeIDLocal( const int iEdge, const int num ) const{

	//if( iEdge < 0 || iEdge >= 6 ){
	//	OutputFiles::m_logFile << " Error : iEdge is out of range in getNodeIDLocalFromElementAndEdge. iEdge = " << iEdge << std::endl;
	//	exit(1);
	//}

	//if( num != 0 && num != 1 ){
	//	OutputFiles::m_logFile << " Error : num must be 0 or 1 in getNodeIDLocalFromElementAndEdge. num = " << num << std::endl;
	//	exit(1);
	//}
	assert( iEdge >= 0 );
	assert( iEdge < 6 );
	assert( num == 0 || num == 1 );

	return m_edgeID2NodeID[iEdge][num];

}

// Get ID of the nodes of specified element and edge
int MeshDataTetraElement::getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getNodeIDFromElementAndEdge. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iEdge < 0 || iEdge >= 6 ){
	//	OutputFiles::m_logFile << " Error : iEdge is out of range in getNodeIDFromElementAndEdge. iEdge = " << iEdge << std::endl;
	//	exit(1);
	//}

	//if( num != 0 && num != 1 ){
	//	OutputFiles::m_logFile << " Error : num must be 0 or 1 in getNodeIDFromElementAndEdge. num = " << num << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iEdge >= 0 );
	assert( iEdge < 6 );
	assert( num == 0 || num == 1 );

	return getNodesOfElements( iElem, m_edgeID2NodeID[iEdge][num] );

}

// Get global node ID of specified element and face
int MeshDataTetraElement::getNodeIDGlobalFromElementAndFace( const int iElem, const int iFace, const int num ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 4 );
	assert( num >= 0 );
	assert( num < 3 );

	return getNodesOfElements( iElem, m_faceID2NodeID[iFace][num] );

}

// Get global node ID of specified element belonging to the boundary planes  
int MeshDataTetraElement::getNodeIDGlobalFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	const int elemID3D = getElemBoundaryPlanes( iPlane, iElem );
	const int faceID3D = getFaceIDLocalFromElementBoundaryPlanes( iPlane, iElem );

	//if( num < 0 || num >= 3 ){
	//	OutputFiles::m_logFile << " Error : num is out of range !! num = " << num << std::endl;
	//	exit(1);
	//}
	assert( num >= 0 );
	assert( num < 3 );

	int comp(num);

	if( iPlane == MeshData::YZMinus || iPlane == MeshData::ZXMinus ){// 2D plane is minus side

		if( num == 1 ){
			comp = 2;
		}else if( num == 2 ){
			comp = 1;
		}

	}

	return getNodesOfElements( elemID3D, m_faceID2NodeID[faceID3D][comp] );

}

// Get X coordinate of node of specified element belonging to the boundary planes  
double MeshDataTetraElement::getCoordXFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getXCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get Y coordinate of node of specified element belonging to the boundary planes  
double MeshDataTetraElement::getCoordYFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getYCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get Z coordinate of node of specified element belonging to the boundary planes  
double MeshDataTetraElement::getCoordZFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getZCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get global node ID from ID of element belonging to the boundary planes and its edge ID
// [note] : node ID is outputed as they make a clockwise turn around +X or +Y direction
int MeshDataTetraElement::getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge, const int num ) const{

	//const int elemID3D = getElemBoundaryPlanes( iPlane, iElem );
	//const int edgeID3D = getEdgeIDLocalFromElementBoundaryPlanes( iPlane, iElem, iEdge );
	//return getNodeIDGlobalFromElementAndEdge( elemID3D, edgeID3D, num );

	//const int elemID3D = getElemBoundaryPlanes( iPlane, iElem );
	//const int faceID3D = getFaceIDLocalFromElementBoundaryPlanes( iPlane, iElem );

	//int comp(0);

	//if( iPlane == MeshData::YZPlus || iPlane == MeshData::ZXPlus ){// node ID is 

	//	if( ( iEdge == 0 && num == 0 ) || ( iEdge == 2 && num == 1 ) ){
	//		comp = 0;
	//	}else if( ( iEdge == 1 && num == 0 ) || ( iEdge == 0 && num == 1 ) ){
	//		comp = 1;
	//	}else if( ( iEdge == 2 && num == 0 ) || ( iEdge == 1 && num == 1 ) ){
	//		comp = 2;
	//	}else{
	//		OutputFiles::m_logFile << " Error : Either iEdge or num is out of range !! iEdge = " << iEdge << ", num = " << num << std::endl;
	//		exit(1);
	//	}

	//}else{

	//	if( ( iEdge == 0 && num == 0 ) || ( iEdge == 2 && num == 1 ) ){
	//		comp = 0;
	//	}else if( ( iEdge == 1 && num == 0 ) || ( iEdge == 0 && num == 1 ) ){
	//		comp = 2;
	//	}else if( ( iEdge == 2 && num == 0 ) || ( iEdge == 1 && num == 1 ) ){
	//		comp = 1;
	//	}else{
	//		OutputFiles::m_logFile << " Error : Either iEdge or num is out of range !! iEdge = " << iEdge << ", num = " << num << std::endl;
	//		exit(1);
	//	}

	//}

	 //	return getNodesOfElements( elemID3D, m_faceID2NodeID[faceID3D][comp] );

	if( ( iEdge == 0 && num == 0 ) || ( iEdge == 2 && num == 1 ) ){
		return getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 0 );
	}else if( ( iEdge == 1 && num == 0 ) || ( iEdge == 0 && num == 1 ) ){
		return getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 1 );
	}else if( ( iEdge == 2 && num == 0 ) || ( iEdge == 1 && num == 1 ) ){
		return getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 2 );
	}else{
		OutputFiles::m_logFile << " Error : Either iEdge or num is out of range !! iEdge = " << iEdge << ", num = " << num << std::endl;
		exit(1);
	}


}

// Get local edge ID from local face ID
int MeshDataTetraElement::getEdgeIDLocalFromFaceIDLocal( const int iFace, const int num ) const{

	//if( iFace < 0 || iFace >= 4 ){
	//	OutputFiles::m_logFile << " Error : iFace is out of range. iFace = " << iFace << std::endl;
	//	exit(1);
	//}

	//if( num < 0 || num >= 3 ){
	//	OutputFiles::m_logFile << " Error : num is out of range. num = " << num << std::endl;
	//	exit(1);
	//}

	assert( iFace >= 0 );
	assert( iFace < 4 );
	assert( num >= 0 );
	assert( num < 3 );

	return m_faceID2EdgeID[iFace][num];

}

#ifdef _DEBUG_WRITE

// Function for debug
void MeshDataTetraElement::testFuction() const{

	//const CommonParameters::locationXY startPointOf1stSegment = { 1.0, 1.0 };
	//const CommonParameters::locationXY endPointOf1stSegment = { 3.0, 3.0 };
	//const CommonParameters::locationXY startPointOf2ndSegment = {  0.0, 0.0 };
	//const CommonParameters::locationXY endPointOf2ndSegment = { 2.0, 2.0 };
	//const bool yesNo = intersectTwoSegments( startPointOf1stSegment, endPointOf1stSegment, startPointOf2ndSegment, endPointOf2ndSegment );
	//std::cout << "intersectTwoSegments : " << yesNo << std::endl;

	//const CommonParameters::locationXY startPointOf1stSegment = { 1.0, 1.0 };
	//const CommonParameters::locationXY endPointOf1stSegment = { 3.0, 3.0 };
	//const CommonParameters::locationXY startPointOf2ndSegment = {  0.0, 0.0 };
	//const CommonParameters::locationXY endPointOf2ndSegment = { 10.0, 10.0 };
	//const bool yesNo = overlapTwoLines( startPointOf1stSegment, endPointOf1stSegment, startPointOf2ndSegment, endPointOf2ndSegment );
	//std::cout << "overlapTwoLines : " << yesNo << std::endl;

	//const CommonParameters::locationXY startPointOf1stSegment = { 1.0, 1.0 };
	//const CommonParameters::locationXY endPointOf1stSegment = { 3.0, 3.0 };
	//const CommonParameters::locationXY startPointOf2ndSegment = {  3.0, 1.0 };
	//const CommonParameters::locationXY endPointOf2ndSegment = { 3.0, 5.0 };
	//CommonParameters::locationXY result = { 0.0, 0.0 };
	//calcCoordOfIntersectionPointOfTwoLines( startPointOf1stSegment, endPointOf1stSegment, startPointOf2ndSegment, endPointOf2ndSegment, result );
	//std::cout << "calcCoordOfIntersectionPointOfTwoLines : " << result.X << " " << result.Y << std::endl;

	//const CommonParameters::locationXY point = { 2.0, 2.0 };
	//const CommonParameters::locationXY startPointOf1stSegment = { 1.0, 1.0 };
	//const CommonParameters::locationXY endPointOf1stSegment = { 3.0, 3.0 };
	//const bool yesNo = locateLeftOfSegmentOnXYPlane( point, startPointOf1stSegment, endPointOf1stSegment );
	//std::cout << "locateLeftOfSegmentOnXYPlane : " << yesNo << std::endl;

}

#endif

// Copy constructer
MeshDataTetraElement::MeshDataTetraElement(const MeshDataTetraElement& rhs){
	std::cerr << "Error : Copy constructer of the class MeshDataTetraElement is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
MeshDataTetraElement& MeshDataTetraElement::operator=(const MeshDataTetraElement& rhs){
	std::cerr << "Error : Assignment operator of the class MeshDataTetraElement is not implemented." << std::endl;
	exit(1);
}

// Function determine if two segments intersect or not
bool MeshDataTetraElement::locateLeftOfSegmentOnLandSurface( const CommonParameters::locationXY& point, 
	const CommonParameters::locationXY& startPointOfSegment, const CommonParameters::locationXY& endPointOfSegment ) const{
		
	if( ( endPointOfSegment.Y - startPointOfSegment.Y )*( point.X - startPointOfSegment.X ) >= ( endPointOfSegment.X - startPointOfSegment.X )*( point.Y - startPointOfSegment.Y ) ){
		return true;
	}

	return false;

}

// Function determine if then inputed point locate at the left of the surface of the upper element
bool MeshDataTetraElement::locateLeftOfSegmentOnSeaSurface( const CommonParameters::locationXY& point, 
	const CommonParameters::locationXY& startPointOfSegment, const CommonParameters::locationXY& endPointOfSegment ) const{
		
	if( ( endPointOfSegment.Y - startPointOfSegment.Y )*( point.X - startPointOfSegment.X ) <= ( endPointOfSegment.X - startPointOfSegment.X )*( point.Y - startPointOfSegment.Y ) ){
		return true;
	}

	return false;

}

// Function determine if then inputed point locate at the left of the segment on the Y-Z plane of boundary
bool MeshDataTetraElement::locateLeftOfSegmentOnYZPlaneOfBoundary( const int iPlane, const int iElem, const int iEdge, const CommonParameters::locationYZ& point ) const{

	//if( iPlane != MeshData::YZMinus && iPlane != MeshData::YZPlus ){
	//	OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iPlane == MeshData::YZMinus || iPlane == MeshData::YZPlus );

	const int nodeID0 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 1 );

	const CommonParameters::locationYZ startPointOfSegment = { getYCoordinatesOfNodes( nodeID0 ), getZCoordinatesOfNodes( nodeID0 ) };
	const CommonParameters::locationYZ endPointOfSegment   = { getYCoordinatesOfNodes( nodeID1 ), getZCoordinatesOfNodes( nodeID1 ) };

	if( ( endPointOfSegment.Z - startPointOfSegment.Z )*( point.Y - startPointOfSegment.Y ) <= ( endPointOfSegment.Y - startPointOfSegment.Y )*( point.Z - startPointOfSegment.Z ) ){
		return true;
	}

	return false;

}

// Function determine if then inputed point locate at the left of the segment on the Z-Y plane of boundary
bool MeshDataTetraElement::locateLeftOfSegmentOnZXPlaneOfBoundary( const int iPlane, const int iElem, const int iEdge, const CommonParameters::locationZX& point ) const{

	//if( iPlane != MeshData::ZXMinus && iPlane != MeshData::ZXPlus ){
	//	OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iPlane == MeshData::ZXMinus || iPlane == MeshData::ZXPlus );

	const int nodeID0 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 1 );

	const CommonParameters::locationZX startPointOfSegment = { getZCoordinatesOfNodes( nodeID0 ), getXCoordinatesOfNodes( nodeID0 ) };
	const CommonParameters::locationZX endPointOfSegment   = { getZCoordinatesOfNodes( nodeID1 ), getXCoordinatesOfNodes( nodeID1 ) };

	if( ( endPointOfSegment.X - startPointOfSegment.X )*( point.Z - startPointOfSegment.Z ) <= ( endPointOfSegment.Z - startPointOfSegment.Z )*( point.X - startPointOfSegment.X ) ){
		return true;
	}

	return false;

}

// Calculate volume of tetrahedron
double MeshDataTetraElement::calcVolume( const CommonParameters::locationXYZ& point1, const CommonParameters::locationXYZ& point2,
	const CommonParameters::locationXYZ& point3, const CommonParameters::locationXYZ& point4 ) const{

	const double val = ( point2.X*point3.Y* point4.Z + point2.Y*point3.Z*point4.X + point2.Z*point3.X*point4.Y - point2.Z*point3.Y*point4.X - point2.X*point3.Z*point4.Y - point2.Y*point3.X*point4.Z )
		             - ( point1.X*point3.Y* point4.Z + point1.Y*point3.Z*point4.X + point1.Z*point3.X*point4.Y - point1.Z*point3.Y*point4.X - point1.X*point3.Z*point4.Y - point1.Y*point3.X*point4.Z )
		             + ( point1.X*point2.Y* point4.Z + point1.Y*point2.Z*point4.X + point1.Z*point2.X*point4.Y - point1.Z*point2.Y*point4.X - point1.X*point2.Z*point4.Y - point1.Y*point2.X*point4.Z )
		             - ( point1.X*point2.Y* point3.Z + point1.Y*point2.Z*point3.X + point1.Z*point2.X*point3.Y - point1.Z*point2.Y*point3.X - point1.X*point2.Z*point3.Y - point1.Y*point2.X*point3.Z );

	return val / 6.0;
}

// Calculate volume coordinates of point on the land surface
void MeshDataTetraElement::calcVolumeCoordsOfPointOnLandSurface( const int elemID, const int faceID, const CommonParameters::locationXY& pointCoord, CommonParameters::VolumeCoords& coords ) const{

	//if( elemID >= m_numElemTotal || elemID < 0 ){
	//	OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
	//	exit(1);
	//}
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	//const CommonParameters::locationXYZ pointCoord3D = { pointCoord.X, pointCoord.Y, 0.0 };

	//calcVolumeCoordsOfPoint( elemID, pointCoord3D, coords );

	const CommonParameters::locationXY pointCoord2D = { pointCoord.X, pointCoord.Y };

	CommonParameters::AreaCoords areaCoord = { 0.0, 0.0, 0.0 };
	calcAreaCoordsOfPointOnLandSurface( elemID, faceID, pointCoord2D, areaCoord );

	calcVolumeCoordFromAreaCoord( faceID, areaCoord, coords );
	
}

// Calculate are of triangle
double MeshDataTetraElement::calcArea(  const CommonParameters::CoordPair& point1, const CommonParameters::CoordPair& point2, const CommonParameters::CoordPair& point3  ) const{

	return 0.5 * fabs( ( point2.first - point1.first ) * ( point3.second - point1.second ) - ( point2.second - point1.second ) * ( point3.first - point1.first ) ); 

}

// Calculate length of edges of elements
double MeshDataTetraElement::calcEdgeLengthFromElementAndEdge( const int iElem, const int iEdge ) const{

	return calcDistanceOfTwoNodes( getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 0 ), getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 1 ) );

}

// Calculate length projected on the horizontal plane of edges of elements
double MeshDataTetraElement::calcEdgeLengthProjectedOnHorizontalPlaneFromElementAndEdge( const int iElem, const int iEdge ) const{

	return calcHorizontalDistanceOfTwoNodes( getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 0 ), getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 1 ) );

}

// Calculate length of edges of elements on boundary planes
double MeshDataTetraElement::calcEdgeLengthFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const{

	const int nodeID0 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 1 );

	const double coordZ0 = getZCoordinatesOfNodes( nodeID0 );
	const double coordZ1 = getZCoordinatesOfNodes( nodeID1 );

	if( iPlane == MeshData::ZXMinus || iPlane == MeshData::ZXPlus ){// Z-X plane
		const double coordX0 = getXCoordinatesOfNodes( nodeID0 );
		const double coordX1 = getXCoordinatesOfNodes( nodeID1 );

#ifdef _DEBUG_WRITE
		std::cout << "coordX0 coordZ0 : " << coordX0 << " " << coordZ0 << std::endl;
		std::cout << "coordX1 coordZ1 : " << coordX1 << " " << coordZ1 << std::endl;
#endif

		return hypot( coordX1 - coordX0, coordZ1 - coordZ0 ); 

	}else if( iPlane == MeshData::YZMinus || iPlane == MeshData::YZPlus ){// Y-Z plane
		const double coordY0 = getYCoordinatesOfNodes( nodeID0 );
		const double coordY1 = getYCoordinatesOfNodes( nodeID1 );

#ifdef _DEBUG_WRITE
		std::cout << "coordY0 coordZ0 : " << coordY0 << " " << coordZ0 << std::endl;
		std::cout << "coordY1 coordZ1 : " << coordY1 << " " << coordZ1 << std::endl;
#endif

		return hypot( coordY1 - coordY0, coordZ1 - coordZ0 );
	}

	OutputFiles::m_logFile << "Error : Wrong plane ID !! : iPlane = " << iPlane << std::endl;
	exit(1);

	return -1;

}

//// Calculate length of edges of element face
//double MeshDataTetraElement::calcEdgeLengthOnElementFace( const int iElem, const int iFace, const int iEdge ) const{
//
//	return calcEdgeLengthFromElementAndEdge( iElem, getEdgeIDLocalFromFaceIDLocal( iFace, iEdge ) );
//
//}

// Calculate horizontal coordinate differences of edges of the elements on boundary planes
double MeshDataTetraElement::calcHorizontalCoordDifferenceBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const{

	const int nodeID0 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 1 );

	if( iPlane == MeshData::ZXMinus || iPlane == MeshData::ZXPlus ){// Z-X plane

		return getXCoordinatesOfNodes( nodeID1 ) - getXCoordinatesOfNodes( nodeID0 ); 

	}else if( iPlane == MeshData::YZMinus || iPlane == MeshData::YZPlus ){// Y-Z plane

		return getYCoordinatesOfNodes( nodeID1 ) - getYCoordinatesOfNodes( nodeID0 ); 

	}

	OutputFiles::m_logFile << "Error : Wrong plane ID !! : iPlane = " << iPlane << std::endl;
	exit(1);

	return -1;

}

// Calculate X coordinate of points on element face
double MeshDataTetraElement::calcXCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iFace < 0 || iFace >= 4 ){
	//	OutputFiles::m_logFile << " Error : iFace is out of range. iFace = " << iFace << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 4 );
	double val(0.0);

	val += getXCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][0] ) ) * areaCoord.coord0; 
	val += getXCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][1] ) ) * areaCoord.coord1; 
	val += getXCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][2] ) ) * areaCoord.coord2; 

	return val;

}

// Calculate Y coordinate of points on element face
double MeshDataTetraElement::calcYCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iFace < 0 || iFace >= 4 ){
	//	OutputFiles::m_logFile << " Error : iFace is out of range. iFace = " << iFace << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 4 );
	double val(0.0);

	val += getYCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][0] ) ) * areaCoord.coord0; 
	val += getYCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][1] ) ) * areaCoord.coord1; 
	val += getYCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][2] ) ) * areaCoord.coord2; 

	return val;

}

// Calculate Z coordinate of points on element face
double MeshDataTetraElement::calcZCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iFace < 0 || iFace >= 4 ){
	//	OutputFiles::m_logFile << " Error : iFace is out of range. iFace = " << iFace << std::endl;
	//	exit(1);
	//}
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 4 );
	double val(0.0);

	val += getZCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][0] ) ) * areaCoord.coord0; 
	val += getZCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][1] ) ) * areaCoord.coord1; 
	val += getZCoordinatesOfNodes( getNodesOfElements( iElem, m_faceID2NodeID[iFace][2] ) ) * areaCoord.coord2; 

	return val;

}

// Calculate volume coordinates from area coordinates
void MeshDataTetraElement::calcVolumeCoordFromAreaCoord( const int iFace, const CommonParameters::AreaCoords& areaCoord, CommonParameters::VolumeCoords& volumeCoord ) const{

	switch (iFace){
		case 0:
			volumeCoord.coord0 = 0.0;
			volumeCoord.coord1 = areaCoord.coord0;
			volumeCoord.coord2 = areaCoord.coord1;
			volumeCoord.coord3 = areaCoord.coord2;
			break;
		case 1:
			volumeCoord.coord1 = 0.0;
			volumeCoord.coord0 = areaCoord.coord0;
			volumeCoord.coord3 = areaCoord.coord1;
			volumeCoord.coord2 = areaCoord.coord2;
			break;
		case 2:
			volumeCoord.coord2 = 0.0;
			volumeCoord.coord0 = areaCoord.coord0;
			volumeCoord.coord1 = areaCoord.coord1;
			volumeCoord.coord3 = areaCoord.coord2;
			break;
		case 3:
			volumeCoord.coord3 = 0.0;
			volumeCoord.coord0 = areaCoord.coord0;
			volumeCoord.coord2 = areaCoord.coord1;
			volumeCoord.coord1 = areaCoord.coord2;
			break;
		default:
			OutputFiles::m_logFile << " Error : iFace is out of range. iFace = " << iFace << std::endl;
			exit(1);
			break;
	}

}

// Calculate area projected on X-Y plane with sign of triangle from area coordinate.
// If rotation direction is plus Z, area value is plus. Otherwise, area value is minus.
double MeshDataTetraElement::calcAreaOnXYPlaneWithSignFromAreaCoords( const int elemID, const int faceID,
	const CommonParameters::AreaCoords& coord0, const CommonParameters::AreaCoords& coord1, const CommonParameters::AreaCoords& coord2 ) const{

	const CommonParameters::CoordPair coordPair0 = { calcXCoordOfPointOnFace( elemID, faceID, coord0 ), calcYCoordOfPointOnFace( elemID, faceID, coord0 ) };
	const CommonParameters::CoordPair coordPair1 = { calcXCoordOfPointOnFace( elemID, faceID, coord1 ), calcYCoordOfPointOnFace( elemID, faceID, coord1 ) };
	const CommonParameters::CoordPair coordPair2 = { calcXCoordOfPointOnFace( elemID, faceID, coord2 ), calcYCoordOfPointOnFace( elemID, faceID, coord2 ) };

	const double val = calcArea( coordPair0, coordPair1, coordPair2 );
	const double outerProduct = ( coordPair1.first - coordPair0.first ) * ( coordPair2.second - coordPair0.second ) - ( coordPair2.first - coordPair0.first ) * ( coordPair1.second - coordPair0.second );

	if( outerProduct > 0.0 ){
		return val;
	}else{
		return -val;
	}	

}

// Calculate area with sign of triangle from area coordinate.
// If rotation direction is plus Z, area value is plus. Otherwise, area value is minus.
double MeshDataTetraElement::calcAreaWithSignFromAreaCoords( const int elemID, const int faceID,
	const CommonParameters::AreaCoords& coord0, const CommonParameters::AreaCoords& coord1, const CommonParameters::AreaCoords& coord2 ) const{

	const CommonParameters::locationXYZ p0 = { calcXCoordOfPointOnFace( elemID, faceID, coord0 ), calcYCoordOfPointOnFace( elemID, faceID, coord0 ), calcZCoordOfPointOnFace( elemID, faceID, coord0 ) };
	const CommonParameters::locationXYZ p1 = { calcXCoordOfPointOnFace( elemID, faceID, coord1 ), calcYCoordOfPointOnFace( elemID, faceID, coord1 ), calcZCoordOfPointOnFace( elemID, faceID, coord1 ) };
	const CommonParameters::locationXYZ p2 = { calcXCoordOfPointOnFace( elemID, faceID, coord2 ), calcYCoordOfPointOnFace( elemID, faceID, coord2 ), calcZCoordOfPointOnFace( elemID, faceID, coord2 ) };

	const CommonParameters::Vector3D vec0 = { p1.X - p0.X, p1.Y - p0.Y, p1.Z - p0.Z };
	const CommonParameters::Vector3D vec1 = { p2.X - p0.X, p2.Y - p0.Y, p2.Z - p0.Z };

	const CommonParameters::Vector3D vecOut = calcOuterProduct( vec0, vec1 );

	if( vecOut.Z > 0.0 ){
		return  0.5 * sqrt( vecOut.X * vecOut.X + vecOut.Y * vecOut.Y + vecOut.Z * vecOut.Z );
	}else{
		return -0.5 * sqrt( vecOut.X * vecOut.X + vecOut.Y * vecOut.Y + vecOut.Z * vecOut.Z );
	}	

}

// Calculate area of face
double MeshDataTetraElement::calcAreaOfFace( const int iElem, const int iFace ) const{

	const int nodeID[3] = {
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 0 ),
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 1 ),
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 2 )
	};

	const CommonParameters::Vector3D vec0 = {
		getXCoordinatesOfNodes(nodeID[1]) - getXCoordinatesOfNodes(nodeID[0]),  
		getYCoordinatesOfNodes(nodeID[1]) - getYCoordinatesOfNodes(nodeID[0]),  
		getZCoordinatesOfNodes(nodeID[1]) - getZCoordinatesOfNodes(nodeID[0])
	};

	const CommonParameters::Vector3D vec1 = {
		getXCoordinatesOfNodes(nodeID[2]) - getXCoordinatesOfNodes(nodeID[0]),  
		getYCoordinatesOfNodes(nodeID[2]) - getYCoordinatesOfNodes(nodeID[0]),  
		getZCoordinatesOfNodes(nodeID[2]) - getZCoordinatesOfNodes(nodeID[0])
	};
	
	const CommonParameters::Vector3D vecOut = calcOuterProduct(vec0, vec1);

	return  0.5 * fabs( sqrt( vecOut.X * vecOut.X + vecOut.Y * vecOut.Y + vecOut.Z * vecOut.Z ) );

}

// Calculate area of face at bottom of mesh
double MeshDataTetraElement::calcAreaOfFaceAtBottomOfMesh( const int iElem ) const{
	const int elemID = getElemBoundaryPlanes(MeshData::XYPlus, iElem);
	const int iFace = getFaceIDLocalFromElementBoundaryPlanes(MeshData::XYPlus, iElem);
	return calcAreaOfFace(elemID, iFace);
}

//// Calculate volume coordinates of point on YZ plane
//void MeshDataTetraElement::calcVolumeCoordsOfPointOnYZPlane( const int elemID, const int faceID, const CommonParameters::locationYZ& pointCoord, CommonParameters::VolumeCoords& coords ) const{
//
//	if( elemID >= m_numElemTotal || elemID < 0 ){
//		OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
//		exit(1);
//	}
//
//	const CommonParameters::locationXYZ pointCoord3D = { 0.0, pointCoord.Y, pointCoord.Z };
//
//	calcVolumeCoordsOfPoint( elemID, pointCoord3D, coords );
//
//}
//
//// Calculate volume coordinates of point on ZX plane
//void MeshDataTetraElement::calcVolumeCoordsOfPointOnZXPlane( const int elemID, const int faceID, const CommonParameters::locationZX& pointCoord, CommonParameters::VolumeCoords& coords ) const{
//
//	if( elemID >= m_numElemTotal || elemID < 0 ){
//		OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
//		exit(1);
//	}
//
//	const CommonParameters::locationXYZ pointCoord3D = { pointCoord.X, 0.0, pointCoord.Z };
//
//	calcVolumeCoordsOfPoint( elemID, pointCoord3D, coords );
//
//}

// Calculate volume coordinates of point
void MeshDataTetraElement::calcVolumeCoordsOfPoint( const int elemID, const CommonParameters::locationXYZ& pointCoord, CommonParameters::VolumeCoords& coords ) const{

	//if( elemID >= m_numElemTotal || elemID < 0 ){
	//	OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
	//	exit(1);
	//}
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	CommonParameters::locationXYZ nodeCoord[4];
	for( int i = 0; i < 4; ++i ){
		const int nodeID = getNodesOfElements( elemID, i ); 
		nodeCoord[i].X = getXCoordinatesOfNodes( nodeID );
		nodeCoord[i].Y = getYCoordinatesOfNodes( nodeID );
		nodeCoord[i].Z = getZCoordinatesOfNodes( nodeID );
	}

	const double volTotal = calcVolume( nodeCoord[0], nodeCoord[1], nodeCoord[2], nodeCoord[3] );

	coords.coord0 = calcVolume(   pointCoord, nodeCoord[1], nodeCoord[2], nodeCoord[3] ) / volTotal;
	coords.coord1 = calcVolume( nodeCoord[0],   pointCoord, nodeCoord[2], nodeCoord[3] ) / volTotal;
	coords.coord2 = calcVolume( nodeCoord[0], nodeCoord[1],   pointCoord, nodeCoord[3] ) / volTotal;
	coords.coord3 = calcVolume( nodeCoord[0], nodeCoord[1], nodeCoord[2],   pointCoord ) / volTotal;

	if( coords.coord0 < - m_eps || coords.coord1 < - m_eps || coords.coord2 < - m_eps || coords.coord3 < - m_eps ){
		OutputFiles::m_logFile << " Error : Volume coordinate of element " << elemID << " is negative for the point (X,Y,Z) = ( " << pointCoord.X << ", " << pointCoord.Y << ", " << pointCoord.Z << " ) ." << std::endl;
		exit(1);
	}

}

// Calculate area coordinates of point on the land surface
void MeshDataTetraElement::calcAreaCoordsOfPointOnLandSurface( const int elemID, const int faceID, const CommonParameters::locationXY& pointCoord, CommonParameters::AreaCoords& coords ) const{

	//if( elemID >= m_numElemTotal || elemID < 0 ){
	//	OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
	//	exit(1);
	//}

	//if( faceID < 0 || faceID >= 4 ){
	//	OutputFiles::m_logFile << "Error : ID of face is out of range !! : faceID = " << faceID << std::endl;
	//	exit(1);
	//}
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );
	assert( faceID >= 0 );
	assert( faceID < 4 );

	CommonParameters::CoordPair nodeCoord[3];
	for( int i = 0; i < 3; ++i ){
		const int nodeID = getNodesOfElements( elemID, m_faceID2NodeID[faceID][i] );
		nodeCoord[i].first  = getXCoordinatesOfNodes( nodeID );
		nodeCoord[i].second = getYCoordinatesOfNodes( nodeID );
	}

	const CommonParameters::CoordPair point = { pointCoord.X, pointCoord.Y };

	const double areaTotal = calcArea( nodeCoord[0], nodeCoord[1], nodeCoord[2] );

	coords.coord0 = calcArea( point, nodeCoord[1], nodeCoord[2] ) / areaTotal;
	coords.coord1 = calcArea( nodeCoord[0], point, nodeCoord[2] ) / areaTotal;
	coords.coord2 = calcArea( nodeCoord[0], nodeCoord[1], point ) / areaTotal;

	if( coords.coord0 < - m_eps || coords.coord1 < - m_eps || coords.coord2 < - m_eps ){
		OutputFiles::m_logFile << " Error : Area coordinate of face " << faceID << " of element " << elemID << " is negative for the point (X,Y) = ( " << pointCoord.X << ", " << pointCoord.Y << " ) ." << std::endl;
		exit(1);
	}

}

// Calculate area coordinates of the specified point on the Y-Z plane of boundary
void MeshDataTetraElement::calcAreaCoordsOfPointOnYZPlaneOfBoundary( const int iPlane, const int iElem, const CommonParameters::CoordPair& point, CommonParameters::AreaCoords& coords ) const{

	//if( iPlane != MeshData::YZMinus && iPlane != MeshData::YZPlus ){
	//	OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iPlane == MeshData::YZMinus || iPlane == MeshData::YZPlus );

	const int nodeID0 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 1 );
	const int nodeID2 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 2 );

	const CommonParameters::CoordPair nodeCoord0 = { getYCoordinatesOfNodes(nodeID0), getZCoordinatesOfNodes(nodeID0) };
	const CommonParameters::CoordPair nodeCoord1 = { getYCoordinatesOfNodes(nodeID1), getZCoordinatesOfNodes(nodeID1) };
	const CommonParameters::CoordPair nodeCoord2 = { getYCoordinatesOfNodes(nodeID2), getZCoordinatesOfNodes(nodeID2) };

	const double areaTotal = calcArea( nodeCoord0, nodeCoord1, nodeCoord2 );

	coords.coord0 = calcArea( point, nodeCoord1, nodeCoord2 ) / areaTotal;
	coords.coord1 = calcArea( nodeCoord0, point, nodeCoord2 ) / areaTotal;
	coords.coord2 = calcArea( nodeCoord0, nodeCoord1, point ) / areaTotal;

	if( coords.coord0 < - m_eps || coords.coord1 < - m_eps || coords.coord2 < - m_eps ){
		OutputFiles::m_logFile << " Error : Area coordinate is negative for the point (Y,Z) = ( " << point.first << ", " << point.second << " ) ." << std::endl;
		exit(1);
	}

}

// Calculate area coordinates of the specified point on the Z-X plane of boundary
void MeshDataTetraElement::calcAreaCoordsOfPointOnZXPlaneOfBoundary( const int iPlane, const int iElem, const CommonParameters::CoordPair& point, CommonParameters::AreaCoords& coords ) const{

	//if( iPlane != MeshData::ZXMinus && iPlane != MeshData::ZXPlus ){
	//	OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iPlane == MeshData::ZXMinus || iPlane == MeshData::ZXPlus );

	const int nodeID0 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 1 );
	const int nodeID2 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 2 );

	const CommonParameters::CoordPair nodeCoord0 = { getZCoordinatesOfNodes(nodeID0), getXCoordinatesOfNodes(nodeID0) };
	const CommonParameters::CoordPair nodeCoord1 = { getZCoordinatesOfNodes(nodeID1), getXCoordinatesOfNodes(nodeID1) };
	const CommonParameters::CoordPair nodeCoord2 = { getZCoordinatesOfNodes(nodeID2), getXCoordinatesOfNodes(nodeID2) };

	const double areaTotal = calcArea( nodeCoord0, nodeCoord1, nodeCoord2 );

	coords.coord0 = calcArea( point, nodeCoord1, nodeCoord2 ) / areaTotal;
	coords.coord1 = calcArea( nodeCoord0, point, nodeCoord2 ) / areaTotal;
	coords.coord2 = calcArea( nodeCoord0, nodeCoord1, point ) / areaTotal;

	if( coords.coord0 < - m_eps || coords.coord1 < - m_eps || coords.coord2 < - m_eps ){
		OutputFiles::m_logFile << " Error : Area coordinate is negative for the point (Z,X) = ( " << point.first << ", " << point.second << " ) ." << std::endl;
		exit(1);
	}

}

//// Calculate are of triangle projected on the X-Y plane
//double MeshDataTetraElement::calcAreaOfTriangleOnXYPlane( const CommonParameters::locationXY& vertex1, const CommonParameters::locationXY& vertex2,
//	const CommonParameters::locationXY& vertex3 ) const{
//
//	return 0.5 * fabs( ( vertex2.X - vertex1.X )*( vertex3.Y - vertex1.Y ) - ( vertex2.Y - vertex1.Y )*( vertex3.X - vertex1.X ) );
//
//}

//// Calulate normal vector of element face
//CommonParameters::Vector3D MeshDataTetraElement::calulateNormalVectorOfElementFace( const int elemID, const int faceID ) const{
//
//	if( elemID >= m_numElemTotal || elemID < 0 ){
//		OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
//		exit(1);
//	}
//
//	const int nodeID0 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][0] );
//	const int nodeID1 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][1] );
//	const int nodeID2 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][2] );
//
//	const CommonParameters::locationXYZ nodeCoord0 = { getXCoordinatesOfNodes( nodeID0 ), getYCoordinatesOfNodes( nodeID0 ), getZCoordinatesOfNodes( nodeID0 ) };
//	const CommonParameters::locationXYZ nodeCoord1 = { getXCoordinatesOfNodes( nodeID1 ), getYCoordinatesOfNodes( nodeID1 ), getZCoordinatesOfNodes( nodeID1 ) };
//	const CommonParameters::locationXYZ nodeCoord2 = { getXCoordinatesOfNodes( nodeID2 ), getYCoordinatesOfNodes( nodeID2 ), getZCoordinatesOfNodes( nodeID2 ) };
//
//	const CommonParameters::Vector3D vec1 = { nodeCoord1.X - nodeCoord0.X, nodeCoord1.Y - nodeCoord0.Y, nodeCoord1.Z - nodeCoord0.Z };
//	const CommonParameters::Vector3D vec2 = { nodeCoord2.X - nodeCoord1.X, nodeCoord2.Y - nodeCoord1.Y, nodeCoord2.Z - nodeCoord1.Z };
//
//	return calcOuterProduct( vec1, vec2 );
//
//}

// Decide whether specified point locate inside of face
//bool MeshDataTetraElement::locateInsideOfFace( const int elemID, const int faceID, const double locX, const double locY, const double locZ ) const{
bool MeshDataTetraElement::locateInsideOfFace( const int elemID, const int faceID, const CommonParameters::locationXYZ& loc ) const{

	//if( elemID >= m_numElemTotal || elemID < 0 ){
	//	OutputFiles::m_logFile << "Error : ID of element is out of range !! : elemID = " << elemID << std::endl;
	//	exit(1);
	//}
	assert( elemID >= 0 );
	assert( elemID < m_numElemTotal );

	CommonParameters::locationXYZ nodeCoord[3];

	for( int i = 0; i < 3; ++i ){
		const int nodeID = getNodesOfElements( elemID, m_faceID2NodeID[faceID][i] );
		nodeCoord[i].X =getXCoordinatesOfNodes( nodeID );
		nodeCoord[i].Y =getYCoordinatesOfNodes( nodeID );
		nodeCoord[i].Z =getZCoordinatesOfNodes( nodeID );
	}

	//const CommonParameters::Vector3D vec1 = { nodeCoord1.X - nodeCoord0.X, nodeCoord1.Y - nodeCoord0.Y, nodeCoord1.Z - nodeCoord0.Z };
	//const CommonParameters::Vector3D vec2 = { nodeCoord2.X - nodeCoord1.X, nodeCoord2.Y - nodeCoord1.Y, nodeCoord2.Z - nodeCoord1.Z };
	//const CommonParameters::Vector3D vec3 = { loc.X - nodeCoord0.X, loc.Y - nodeCoord0.Y, locZ - nodeCoord0.Z };

	const CommonParameters::Vector3D vec1 = calcVector3D( nodeCoord[0], nodeCoord[1] );
	const CommonParameters::Vector3D vec2 = calcVector3D( nodeCoord[1], nodeCoord[2] );
	const CommonParameters::Vector3D vec3 = calcVector3D( nodeCoord[0], loc );

	if( calcInnerProduct( calcOuterProduct( vec1, vec2 ), vec3 ) <= 0.0 ){
		return true;
	}

	return false;

}

//// Decide whether specified point is included in element
//bool MeshDataTetraElement::includePointInElement( const int elemID, const double locX, const double locY, const double locZ ) const{
//
//	for( int i = 0; i < 4; ++i ){
//		if( !locateInsideOfFace( elemID, i, locX, locY, locZ ) ){
//			return false;
//		}
//	}
//
//	return true;
//
//}
