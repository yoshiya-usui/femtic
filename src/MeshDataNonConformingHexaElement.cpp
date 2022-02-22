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
#include <assert.h>

#include "AnalysisControl.h"
#include "MeshDataNonConformingHexaElement.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"
#include "OutputFiles.h"
#include "Util.h"

// Constructer
MeshDataNonConformingHexaElement::MeshDataNonConformingHexaElement():
	m_neighborElementsForNonConformingHexa(NULL),
	m_numElemOnLandSurface(0),
	m_elemOnLandSurface(NULL),
	m_faceLandSurface(NULL)
{
	m_numNodeOneElement = 8;
	m_numEdgeOneElement = 12;
	m_numNodeOnFaceOneElement = 4;
	m_numNeighborElement = -1;// This variable is not used

	for ( int i = 0; i < 6; ++i ){
		m_facesOfElementsBoundaryPlanes[i] = NULL;
	}

	// Array converting from face ID to node ID
	m_faceID2NodeID[0][0] = 0;
	m_faceID2NodeID[0][1] = 3;
	m_faceID2NodeID[0][2] = 7;
	m_faceID2NodeID[0][3] = 4;

	m_faceID2NodeID[1][0] = 1;
	m_faceID2NodeID[1][1] = 2;
	m_faceID2NodeID[1][2] = 6;
	m_faceID2NodeID[1][3] = 5;

	m_faceID2NodeID[2][0] = 0;
	m_faceID2NodeID[2][1] = 1;
	m_faceID2NodeID[2][2] = 5;
	m_faceID2NodeID[2][3] = 4;

	m_faceID2NodeID[3][0] = 3;
	m_faceID2NodeID[3][1] = 2;
	m_faceID2NodeID[3][2] = 6;
	m_faceID2NodeID[3][3] = 7;

	m_faceID2NodeID[4][0] = 0;
	m_faceID2NodeID[4][1] = 1;
	m_faceID2NodeID[4][2] = 2;
	m_faceID2NodeID[4][3] = 3;

	m_faceID2NodeID[5][0] = 4;
	m_faceID2NodeID[5][1] = 5;
	m_faceID2NodeID[5][2] = 6;
	m_faceID2NodeID[5][3] = 7;

	// Array converting from face ID to edge ID
	m_faceID2EdgeID[0][0] = 4;
	m_faceID2EdgeID[0][1] = 5;
	m_faceID2EdgeID[0][2] = 8;
	m_faceID2EdgeID[0][3] = 10;

	m_faceID2EdgeID[1][0] = 6;
	m_faceID2EdgeID[1][1] = 7;
	m_faceID2EdgeID[1][2] = 9;
	m_faceID2EdgeID[1][3] = 11;

	m_faceID2EdgeID[2][0] = 0;
	m_faceID2EdgeID[2][1] = 2;
	m_faceID2EdgeID[2][2] = 8;
	m_faceID2EdgeID[2][3] = 9;

	m_faceID2EdgeID[3][0] = 1;
	m_faceID2EdgeID[3][1] = 3;
	m_faceID2EdgeID[3][2] = 10;
	m_faceID2EdgeID[3][3] = 11;

	m_faceID2EdgeID[4][0] = 0;
	m_faceID2EdgeID[4][1] = 1;
	m_faceID2EdgeID[4][2] = 4;
	m_faceID2EdgeID[4][3] = 6;

	m_faceID2EdgeID[5][0] = 2;
	m_faceID2EdgeID[5][1] = 3;
	m_faceID2EdgeID[5][2] = 5;
	m_faceID2EdgeID[5][3] = 7;

	// Array converting from edge ID to node ID
	m_edgeID2NodeID[0][0] = 0;
	m_edgeID2NodeID[0][1] = 1;

	m_edgeID2NodeID[1][0] = 3;
	m_edgeID2NodeID[1][1] = 2;

	m_edgeID2NodeID[2][0] = 4;
	m_edgeID2NodeID[2][1] = 5;

	m_edgeID2NodeID[3][0] = 7;
	m_edgeID2NodeID[3][1] = 6;

	m_edgeID2NodeID[4][0] = 0;
	m_edgeID2NodeID[4][1] = 3;

	m_edgeID2NodeID[5][0] = 4;
	m_edgeID2NodeID[5][1] = 7;

	m_edgeID2NodeID[6][0] = 1;
	m_edgeID2NodeID[6][1] = 2;

	m_edgeID2NodeID[7][0] = 5;
	m_edgeID2NodeID[7][1] = 6;

	m_edgeID2NodeID[8][0] = 0;
	m_edgeID2NodeID[8][1] = 4;

	m_edgeID2NodeID[9][0] = 1;
	m_edgeID2NodeID[9][1] = 5;

	m_edgeID2NodeID[10][0] = 3;
	m_edgeID2NodeID[10][1] = 7;

	m_edgeID2NodeID[11][0] = 2;
	m_edgeID2NodeID[11][1] = 6;

	// Calculate integral points and weights of two point Gauss quadrature
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

}

// Destructer
MeshDataNonConformingHexaElement::~MeshDataNonConformingHexaElement(){

	//if( m_subElements != NULL ){
	//	delete[] m_subElements;
	//	m_subElements = NULL;
	//}

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

	if( m_neighborElementsForNonConformingHexa != NULL ){
		delete[] m_neighborElementsForNonConformingHexa;
		m_neighborElementsForNonConformingHexa = NULL;
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

// Copy constructer
MeshDataNonConformingHexaElement::MeshDataNonConformingHexaElement(const MeshDataNonConformingHexaElement& rhs){
	std::cerr << "Error : Copy constructer of the class MeshDataNonConformingHexaElement is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
MeshDataNonConformingHexaElement& MeshDataNonConformingHexaElement::operator=(const MeshDataNonConformingHexaElement& rhs){
	std::cerr << "Error : Assignment operator of the class MeshDataNonConformingHexaElement is not implemented." << std::endl;
	exit(1);
}

// Input mesh data from "mesh.dat"
void MeshDataNonConformingHexaElement::inputMeshData(){

	std::ifstream inFile( "mesh.dat", std::ios::in );
	if( inFile.fail() )
	{
		std::cerr << "File open error : mesh.dat !!" << std::endl;
		exit(1);
	}

	std::string sbuf;
	inFile >> sbuf;
	if( sbuf.substr(0,5).compare("DHEXA") != 0 ){
		OutputFiles::m_logFile << "Mesh data written in mesh.dat is different from those for nonconforming hexahedral element !!" << std::endl;
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
		assert( idum == iNode ); 
	}

	inFile >> ibuf;
	if( ibuf > 0 ){
		m_numElemTotal = ibuf;
	}else{
		OutputFiles::m_logFile << "Total number of elements is less than or equal to zero ! : " << ibuf << std::endl;
		exit(1);
	}

	if( m_nodesOfElements == NULL ){
		delete[] m_nodesOfElements;
	}
	m_nodesOfElements = new int[ m_numElemTotal * m_numNodeOneElement ];

	if( m_neighborElementsForNonConformingHexa != NULL ){
		delete[] m_neighborElementsForNonConformingHexa;
	}
	m_neighborElementsForNonConformingHexa = new std::vector<int>[ m_numElemTotal * 6 ];

	for( int iElem = 0; iElem < m_numElemTotal; ++iElem ){
		int idum(0);
		inFile >> idum;
		assert( idum == iElem ); 
		// Nodes of the element
		for( int i = 0; i < m_numNodeOneElement; ++i ){
			int node(0);
			inFile >> node;
			m_nodesOfElements[ iElem * m_numNodeOneElement + i ] = node;
		}
		// IDs of neighbor Elements
		for( int i = 0; i < 6; ++i ){
			int nFace(-1);
			inFile >> nFace;
			for( int iFace = 0; iFace < nFace; ++iFace ){
				int ibuf(-1);
				inFile >> ibuf;
				m_neighborElementsForNonConformingHexa[ iElem * 6 + i ].push_back(ibuf);
			}
		}
	}

	// Check whether side element-faces are parallel to Z-X or Y-Z plane
	checkWhetherSideFaceIsParallelToZXOrYZPlane();

	for( int iPlane = 0; iPlane < 6; ++iPlane ){// Loop of boundary planes
		int nElemOnPlane;
		inFile >> nElemOnPlane;
		if( nElemOnPlane > 0 ){
			m_numElemOnBoundaryPlanes[iPlane] = nElemOnPlane;
		}else{
			OutputFiles::m_logFile << "Number of faces belonging plane " << iPlane << " is less than or equal to zero ! : " << nElemOnPlane << std::endl;
			exit(1);
		}

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
			if( m_facesOfElementsBoundaryPlanes[iPlane][iElem] < 0 || m_facesOfElementsBoundaryPlanes[iPlane][iElem] >= 6 ){
				OutputFiles::m_logFile << "Face ID of plane " << iPlane << " is out of range !! : " << m_facesOfElementsBoundaryPlanes[iPlane][iElem] << std::endl;
				exit(1);
			}
		}
	}

	inFile >> ibuf;
	if( ibuf > 0 ){
		m_numElemOnLandSurface = ibuf;
	}else{
		OutputFiles::m_logFile << "Total number of faces on the land surface is less than or equal to zero ! : " << ibuf << std::endl;
		exit(1);
	}

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
		if( m_faceLandSurface[iElem] < 0 || m_faceLandSurface[iElem] >= 6 ){
			OutputFiles::m_logFile << "Face ID of land surface is out of range !! : " << m_faceLandSurface[iElem] << std::endl;
			exit(1);
		}
	}

	inFile.close();

}

// Find element including a point
int MeshDataNonConformingHexaElement::findElementIncludingPoint( const double locX, const double locY, const double locZ, double& xi, double& eta, double& zeta ) const{

	const int numElemTotal = getNumElemTotal();
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){
		if( isLocatedInTheElement( locX, locY, locZ, iElem ) ){
			calcLocalCoordinates( iElem, locX, locY, locZ, xi, eta, zeta );
			return iElem;
		}
	}

	OutputFiles::m_logFile << " Error : Could not find element including point ( X = " << locX << "[m], Y= " << locY << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	exit(1);
	return -1;

}

// Find elements including point on the surface of the earth
int MeshDataNonConformingHexaElement::findElementIncludingPointOnSurface( const double locX, const double locY, int& faceID, double& xi, double& eta, double& zeta,
	const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod ) const{

	const double eps = 1.0e-9;
	if( useUpperElem ){
		for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){
			const int elemID = m_elemOnLandSurface[iElem];
			faceID = m_faceLandSurface[iElem];
			assert( faceID == 4 );
			// Find the upper element of the site from the lower element
			std::vector<int>& neibElem = m_neighborElementsForNonConformingHexa[6 * elemID + faceID];
			for( std::vector<int>::const_iterator itrNeib = neibElem.begin(); itrNeib != neibElem.end(); ++itrNeib ){
				const int neibElemID = *itrNeib;
				for( int iFace = 0; iFace < 6; ++iFace ){
					std::vector<int>& neibNeibElem = m_neighborElementsForNonConformingHexa[6 * neibElemID + iFace];
					for( std::vector<int>::const_iterator itrNeibNeib = neibNeibElem.begin(); itrNeibNeib != neibNeibElem.end(); ++itrNeibNeib ){
						if( *itrNeibNeib == elemID ){
							assert( iFace == 5 );
							calcHorizontalLocalCoordinates( neibElemID, locX, locY, xi, eta );
							if( xi <= 1.0 + eps && xi >= -1.0 - eps && eta <= 1.0 + eps && eta >= -1.0 - eps ){
								// Point locates in the element 
								if( xi > 1.0 ){
									xi = 1.0;
								}
								if( eta > 1.0 ){
									eta = 1.0;
								}
								if( xi < -1.0 ){
									xi = -1.0;
								}
								if( eta < -1.0 ){
									eta = -1.0;
								}
								zeta = 1.0;
								if( modLoc ){
									xi = 0.0;
									eta = 0.0;
									locXMod = calcXCoordOfPointOnFace(neibElemID, iFace, xi, eta);
									locYMod = calcYCoordOfPointOnFace(neibElemID, iFace, xi, eta);
								}
								faceID = iFace;
								return neibElemID;
							}
						}
					}
				}
			}
		}
	}else{
		for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){
			const int elemID = m_elemOnLandSurface[iElem];
			faceID = m_faceLandSurface[iElem];
			assert( faceID == 4 );		
			calcHorizontalLocalCoordinates( elemID, locX, locY, xi, eta );
			if( xi <= 1.0 + eps && xi >= -1.0 - eps && eta <= 1.0 + eps && eta >= -1.0 - eps ){
				// Point locates in the element 
				if( xi > 1.0 ){
					xi = 1.0;
				}
				if( eta > 1.0 ){
					eta = 1.0;
				}
				if( xi < -1.0 ){
					xi = -1.0;
				}
				if( eta < -1.0 ){
					eta = -1.0;
				}
				zeta = -1.0;
				if( modLoc ){
					xi = 0.0;
					eta = 0.0;
					locXMod = calcXCoordOfPointOnFace(elemID, faceID, xi, eta);
					locYMod = calcYCoordOfPointOnFace(elemID, faceID, xi, eta);
				}
				return elemID;
			}
		}
	}
	
	OutputFiles::m_logFile << " Error : Could not find element including point ( X = " << locX << "[m], Y= " << locY << "[m] )." << std::endl;
	exit(1);
	return -1;

}

// Find element including a point on the Y-Z plane and return element ID of 2D mesh
void MeshDataNonConformingHexaElement::findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd,
	std::vector<int>& elements, std::vector<double>& localCoordXStartPoint, std::vector<double>& localCoordYStartPoint,
	std::vector<double>& localCoordXEndPoint, std::vector<double>& localCoordYEndPoint ) const{

	const double thresholdVal = 1.0E-6;

	const CommonParameters::locationXY nodeCoordDipoleStart = { locXStart, locYStart };
	const CommonParameters::locationXY nodeCoordDipoleEnd = { locXEnd, locYEnd };

	const double dipoleLength = hypot( ( locXEnd - locXStart ), ( locYEnd - locYStart ) );
	if( dipoleLength < thresholdVal ){
		OutputFiles::m_logFile << " Error : Length of dipole ( " << dipoleLength << " [m] ) is too small !! " << std::endl;
		exit(1);
	}

	double dummy(0.0);
	int faceID(0);

	double localCoordStart[3] = { 0.0, 0.0, 0.0 };
	const int iElemStart = findElementIncludingPointOnSurface( locXStart, locYStart, faceID,
		localCoordStart[0], localCoordStart[1], localCoordStart[2], false, false, dummy, dummy );
	assert(faceID == 4);

	double localCoordEnd[3] = { 0.0, 0.0, 0.0 };
	const int iElemEnd = findElementIncludingPointOnSurface( locXEnd, locYEnd, faceID,
		localCoordEnd[0], localCoordEnd[1], localCoordEnd[2], false, false, dummy, dummy );
	assert(faceID == 4);

	std::vector<CommonParameters::locationXY> intersectPointsAreadyFound[2];
	for( int iElem = 0; iElem < m_numElemOnLandSurface; ++iElem ){
		const int elemID = m_elemOnLandSurface[iElem];
		const int faceID = m_faceLandSurface[iElem];

		const int nodeID0 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][0] );
		const int nodeID1 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][1] );
		const int nodeID2 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][2] );
		const int nodeID3 = getNodesOfElements( elemID, m_faceID2NodeID[faceID][3] );

		const CommonParameters::locationXY nodeCoord0 = { getXCoordinatesOfNodes( nodeID0 ), getYCoordinatesOfNodes( nodeID0 ) };
		const CommonParameters::locationXY nodeCoord1 = { getXCoordinatesOfNodes( nodeID1 ), getYCoordinatesOfNodes( nodeID1 ) };
		const CommonParameters::locationXY nodeCoord2 = { getXCoordinatesOfNodes( nodeID2 ), getYCoordinatesOfNodes( nodeID2 ) };
		const CommonParameters::locationXY nodeCoord3 = { getXCoordinatesOfNodes( nodeID3 ), getYCoordinatesOfNodes( nodeID3 ) };

		CommonParameters::locationXY intersectPoints[4];
		if( overlapTwoSegments( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){// node0 - node1
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
		}else if( overlapTwoSegments( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){// node1 - node2
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
		}else if( overlapTwoSegments( nodeCoord2, nodeCoord3, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){// node2 - node3
			const double innerProduct1 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord2 ) / dipoleLength;
			const double innerProduct2 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord3 ) / dipoleLength;
			if( innerProduct1 < innerProduct2 ){
				if( innerProduct1 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord2;
				}
				if( innerProduct2 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord3;
				}
			}else{
				if( innerProduct2 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord3;
				}
				if( innerProduct1 > dipoleLength ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					intersectPoints[1] = nodeCoord2;
				}
			}
		}else if( overlapTwoSegments( nodeCoord3, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){// node3 - node0
			const double innerProduct1 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord3 ) / dipoleLength;
			const double innerProduct2 = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, nodeCoordDipoleStart, nodeCoord0 ) / dipoleLength;
			if( innerProduct1 < innerProduct2 ){
				if( innerProduct1 < 0.0 ){
					intersectPoints[0] = nodeCoordDipoleStart;
				}else{
					intersectPoints[0] = nodeCoord3;
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
					intersectPoints[1] = nodeCoord3;
				}
			}	
		}else{
			int icount(0);
			if( intersectTwoSegments( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord0, nodeCoord1, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				++icount;
			}
			if( intersectTwoSegments( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord1, nodeCoord2, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				bool duplicated(false);
				for( int i = 0; i < icount; ++i ){
					if( calcDistance(intersectPoints[icount], intersectPoints[i]) < 1.0e-9 ){
						duplicated = true;
					}
				}
				if(!duplicated){
					++icount;
				}
			}
			if( intersectTwoSegments( nodeCoord2, nodeCoord3, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord2, nodeCoord3, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				bool duplicated(false);
				for( int i = 0; i < icount; ++i ){
					if( calcDistance(intersectPoints[icount], intersectPoints[i]) < 1.0e-9 ){
						duplicated = true;
					}
				}
				if(!duplicated){
					++icount;
				}
			}
			if( intersectTwoSegments( nodeCoord3, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd ) ){
				calcCoordOfIntersectionPointOfTwoLines( nodeCoord3, nodeCoord0, nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[icount] );
				bool duplicated(false);
				for( int i = 0; i < icount; ++i ){
					if( calcDistance(intersectPoints[icount], intersectPoints[i]) < 1.0e-9 ){
						duplicated = true;
					}
				}
				if(!duplicated){
					++icount;
				}
			}
			if( elemID == iElemStart && elemID == iElemEnd ){
				intersectPoints[0] = nodeCoordDipoleStart;
				intersectPoints[1] = nodeCoordDipoleEnd;
			}else if( icount == 0 ){
				continue;
			}else if( icount == 1 ){
				if( elemID == iElemStart ){
					intersectPoints[1] = nodeCoordDipoleStart;
				}else if( elemID == iElemEnd ){
					intersectPoints[1] = nodeCoordDipoleEnd;
				}else{
					continue;
				}
			}else if( icount > 2 ){
				OutputFiles::m_logFile << " Error : Number of intersection points is larger than two: " << icount << std::endl;
				exit(1);
			}

		}
		const double innerProduct = calcInnerProduct2D( nodeCoordDipoleStart, nodeCoordDipoleEnd, intersectPoints[0], intersectPoints[1] );
		if( fabs(innerProduct) < thresholdVal ){
			// Segment is too short
			continue;
		}
		if( innerProduct < 0 ){
			const CommonParameters::locationXY temp = intersectPoints[0];
			intersectPoints[0] = intersectPoints[1];
			intersectPoints[1] = temp;
		}
		bool alreadyFound = false;
		const int numSegments = static_cast<int>( intersectPointsAreadyFound[0].size() );
		for( int iSeg = 0; iSeg < numSegments; ++iSeg ){
			const CommonParameters::locationXY coordStart = intersectPointsAreadyFound[0][iSeg];
			const CommonParameters::locationXY coordEnd = intersectPointsAreadyFound[1][iSeg];
			if( does1stSegmentContain2ndSegment(coordStart, coordEnd, intersectPoints[0], intersectPoints[1]) ){
				alreadyFound = true;
			}
			if( does1stSegmentContain2ndSegment(intersectPoints[0], intersectPoints[1], coordStart, coordEnd) ){
				alreadyFound = true;
			}
		}
		if( alreadyFound ){
			// Segment has already been found
			continue;
		}
		CommonParameters::locationXY localCoord = {0.0, 0.0};
		calcHorizontalLocalCoordinates( elemID, intersectPoints[0].X, intersectPoints[0].Y, localCoord.X, localCoord.Y );
		localCoordXStartPoint.push_back(localCoord.X);
		localCoordYStartPoint.push_back(localCoord.Y);
		calcHorizontalLocalCoordinates( elemID, intersectPoints[1].X, intersectPoints[1].Y, localCoord.X, localCoord.Y );
		localCoordXEndPoint.push_back(localCoord.X);
		localCoordYEndPoint.push_back(localCoord.Y);
		elements.push_back( elemID );
		for( int i = 0; i < 2; ++i ){
			intersectPointsAreadyFound[i].push_back( intersectPoints[i] );
		}
	}

	if( elements.empty() ){
		OutputFiles::m_logFile << " Error : Could not find element including dipole ( " << locXStart << ", " << locYStart << " ) => ( " << locXEnd << ", " << locYEnd << " )." << std::endl;
		exit(1);
	}

	// For check
	double dipoleLengthAccumulated(0.0);
	std::vector<CommonParameters::locationXY>::iterator itr0End = intersectPointsAreadyFound[0].end(); 
	std::vector<CommonParameters::locationXY>::iterator itr0  = intersectPointsAreadyFound[0].begin(); 
	std::vector<CommonParameters::locationXY>::iterator itr1  = intersectPointsAreadyFound[1].begin(); 
	while( itr0 != itr0End ){
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

// Find element including a point on the Y-Z plane and return element ID of 2D mesh
int MeshDataNonConformingHexaElement::findElementIncludingPointOnYZPlaneAndReturnElemID2D( const int iPlane, const double locY, const double locZ, double& xi, double& eta ) const{

	if( iPlane != MeshData::YZMinus && iPlane != MeshData::YZPlus ){
		OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
		exit(1);
	}

	const double eps = 1.0e-9;
	const int nElem = getNumElemOnBoundaryPlanes(iPlane);
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int nodeID0 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 0 );
		const int nodeID2 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 2 );

		const double h = locY;
		const double v = locZ;
		const double hMin = getYCoordinatesOfNodes(nodeID0);
		const double hMax = getYCoordinatesOfNodes(nodeID2);
		const double vMin = getZCoordinatesOfNodes(nodeID0);
		const double vMax = getZCoordinatesOfNodes(nodeID2);

		const double hMid = 0.5 * (hMin + hMax);
		const double vMid = 0.5 * (vMin + vMax);
		xi   = 2.0 * (h - hMid) / (hMax - hMin);
		eta  = 2.0 * (v - vMid) / (vMax - vMin);

		if( xi <= 1.0 + eps && xi >= -1.0 - eps && eta <= 1.0 + eps && eta >= -1.0 - eps ){
			if( xi > 1.0 ){
				xi = 1.0;
			}
			if( eta > 1.0 ){
				eta = 1.0;
			}
			if( xi < -1.0 ){
				xi = -1.0;
			}
			if( eta < -1.0 ){
				eta = -1.0;
			}
			return iElem;
		}
	}

	OutputFiles::m_logFile << " Warning : Could not find element including point ( Y= " << locY << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	return -1;

}

// Find element including a point on the Z-X plane and return element ID of 2D mesh
int MeshDataNonConformingHexaElement::findElementIncludingPointOnZXPlaneAndReturnElemID2D( const int iPlane, const double locX, const double locZ, double& xi, double& eta ) const{

	if( iPlane != MeshData::ZXMinus && iPlane != MeshData::ZXPlus ){
		OutputFiles::m_logFile << " Error : iPlane is wrong !! iPlane = " << iPlane << std::endl;
		exit(1);
	}

	const double eps = 1.0e-9;
	const int nElem = getNumElemOnBoundaryPlanes(iPlane);
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int nodeID0 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 0 );
		const int nodeID2 = getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, 2 );

		const double h = locX;
		const double v = locZ;
		const double hMin = getXCoordinatesOfNodes(nodeID0);
		const double hMax = getXCoordinatesOfNodes(nodeID2);
		const double vMin = getZCoordinatesOfNodes(nodeID0);
		const double vMax = getZCoordinatesOfNodes(nodeID2);

		const double hMid = 0.5 * (hMin + hMax);
		const double vMid = 0.5 * (vMin + vMax);
		xi   = 2.0 * (h - hMid) / (hMax - hMin);
		eta  = 2.0 * (v - vMid) / (vMax - vMin);

		if( xi <= 1.0 + eps && xi >= -1.0 - eps && eta <= 1.0 + eps && eta >= -1.0 - eps ){
			if( xi > 1.0 ){
				xi = 1.0;
			}
			if( eta > 1.0 ){
				eta = 1.0;
			}
			if( xi < -1.0 ){
				xi = -1.0;
			}
			if( eta < -1.0 ){
				eta = -1.0;
			}
			return iElem;
		}
	}

	OutputFiles::m_logFile << " Warning : Could not find element including point ( X= " << locX << "[m], Z= " << locZ <<  "[m] )." << std::endl;
	return -1;

}

// Get mesh type
int MeshDataNonConformingHexaElement::getMeshType() const{

	return MeshData::NONCONFORMING_HEXA;

}

// Get ID of a neighbor element
int MeshDataNonConformingHexaElement::getIDOfNeighborElement( const int iElem, const int iFace, const int num ) const{

	assert( num >= 0 );
	assert( num < getNumNeighborElement(iElem, iFace) );

	return m_neighborElementsForNonConformingHexa[iElem * 6 + iFace][num];

}

// Get flag specifing whether an element face has slave faces
bool MeshDataNonConformingHexaElement::faceSlaveElements( const int iElem, const int iFace ) const{

	return getNumNeighborElement(iElem, iFace) > 1;

}

// Get flag specifing whether an element face is outer boundary
bool MeshDataNonConformingHexaElement::isOuterBoundary( const int iElem, const int iFace ) const{

	return getNumNeighborElement(iElem, iFace) == 0;

}

// Get number of neighbor elements for an element-face
int MeshDataNonConformingHexaElement::getNumNeighborElement( const int iElem, const int iFace ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 6 );

	return static_cast<int>( m_neighborElementsForNonConformingHexa[iElem * 6 + iFace].size() );

}

// Get local face ID of elements belonging to the boundary planes
int MeshDataNonConformingHexaElement::getFaceIDLocalFromElementBoundaryPlanes( const int iPlane, const int iElem ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iPlane >= 0 );
	assert( iPlane < 6 );

	return m_facesOfElementsBoundaryPlanes[iPlane][iElem];

}

// Get global node ID of specified element and edge
int MeshDataNonConformingHexaElement::getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iEdge >= 0 );
	assert( iEdge < 12 );
	assert( num == 0 || num == 1 );

	return getNodesOfElements( iElem, m_edgeID2NodeID[iEdge][num] );

}

// Get global node ID of specified element and face
int MeshDataNonConformingHexaElement::getNodeIDGlobalFromElementAndFace( const int iElem, const int iFace, const int num ) const{

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iFace >= 0 );
	assert( iFace < 6 );
	assert( num >= 0 );
	assert( num < 4 );

	return getNodesOfElements( iElem, m_faceID2NodeID[iFace][num] );

}

// Get global node ID of specified element belonging to the boundary planes  
int MeshDataNonConformingHexaElement::getNodeIDGlobalFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	const int elemID3D = getElemBoundaryPlanes( iPlane, iElem );
	const int faceID3D = getFaceIDLocalFromElementBoundaryPlanes( iPlane, iElem );

	assert( num >= 0 );
	assert( num < 4 );

	//return getNodesOfElements( elemID3D, m_faceID2NodeID[faceID3D][num] );
	return getNodeIDGlobalFromElementAndFace( elemID3D, faceID3D, num );

}

// Get global node ID from ID of element belonging to the boundary planes and its edge index
int MeshDataNonConformingHexaElement::getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge, const int num ) const{

	assert( iEdge >= 0 || iEdge < 4 );
	assert( num == 0 || num == 1 );

	const int elemID3D = getElemBoundaryPlanes( iPlane, iElem );
	const int faceID = getFaceIDLocalFromElementBoundaryPlanes( iPlane, iElem );
	const int edgeID = getEdgeIDLocalFromFaceIDLocal( faceID, iEdge );

	//return getNodesOfElements( elemID3D, m_edgeID2NodeID[edgeID][num] );
	return getNodeIDGlobalFromElementAndEdge( elemID3D, edgeID, num );

}

// Get X coordinate of node of specified element belonging to the boundary planes  
double MeshDataNonConformingHexaElement::getCoordXFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getXCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get Y coordinate of node of specified element belonging to the boundary planes  
double MeshDataNonConformingHexaElement::getCoordYFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getYCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get Z coordinate of node of specified element belonging to the boundary planes  
double MeshDataNonConformingHexaElement::getCoordZFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const{

	return getZCoordinatesOfNodes( getNodeIDGlobalFromElementBoundaryPlanes( iPlane, iElem, num ) );

}

// Get local edge ID from local face ID
int MeshDataNonConformingHexaElement::getEdgeIDLocalFromFaceIDLocal( const int iFace, const int num ) const{

	assert( iFace >= 0 );
	assert( iFace < 6 );
	assert( num >= 0 );
	assert( num < 4 );

	return m_faceID2EdgeID[iFace][num];

}

// Decide whether specified elements share same edges
bool MeshDataNonConformingHexaElement::shareSameEdges( const int elemID1, const int elemID2 ) const{

	OutputFiles::m_logFile << "Error : MeshDataTetraElement::shareSameEdges is not implemented" << std::endl;
	exit(1);

}

// Calculate volume of a specified element
double MeshDataNonConformingHexaElement::calcVolume( const int elemID ) const{

	double volume(0.0);
	for( int ip = 0; ip < 8; ++ip ){
		const double xi = m_integralPointXi[ip];
		const double eta = m_integralPointEta[ip];
		const double zeta = m_integralPointZeta[ip];
		volume += calcDeterminantOfJacobianMatrix( elemID, xi, eta, zeta ) * m_weights[ip];
	}

	return volume;

}

// Output vtk file
void MeshDataNonConformingHexaElement::outputMeshDataToVTK() const{

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
void MeshDataNonConformingHexaElement::outputMeshDataToBinary() const{

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

// Get array of nodes of elements belonging to the boundary planes
int MeshDataNonConformingHexaElement::getNodesOfElementsBoundaryPlanes( const int iPlane, const int iElem, const int iNode ) const{

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

// Calculate horizontal coordinate differences of edges of the elements on boundary planes
double MeshDataNonConformingHexaElement::calcHorizontalCoordDifferenceBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const{

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

// Interpolate x coordinate on top or bottom face from local coordinate of horizontal plane
double MeshDataNonConformingHexaElement::calcXCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const{

	assert( iFace == 4 || iFace == 5 );

	// Array of reference coord xi values for each node
	const double xiAtNode[4] = { -1.0, 1.0, 1.0, -1.0 };

	// Array of reference coord eta values for each node
	const double etaAtNode[4] = { -1.0, -1.0, 1.0, 1.0 };

	double coordOut(0.0);
	for( int i = 0; i < 4; ++i ){
		const int nodeID = getNodeIDGlobalFromElementAndFace(iElem, iFace, i );
		const double coordAtNode = getXCoordinatesOfNodes(nodeID);
		coordOut += 0.25 * (1.0 + xi * xiAtNode[i]) * (1.0 + eta * etaAtNode[i]) * coordAtNode;
	}
	return coordOut;

}

// Interpolate y coordinate on top or bottom face from local coordinate of horizontal plane
double MeshDataNonConformingHexaElement::calcYCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const{

	assert( iFace == 4 || iFace == 5 );

	// Array of reference coord xi values for each node
	const double xiAtNode[4] = { -1.0, 1.0, 1.0, -1.0 };

	// Array of reference coord eta values for each node
	const double etaAtNode[4] = { -1.0, -1.0, 1.0, 1.0 };

	double coordOut(0.0);
	for( int i = 0; i < 4; ++i ){
		const int nodeID = getNodeIDGlobalFromElementAndFace(iElem, iFace, i );
		const double coordAtNode = getYCoordinatesOfNodes(nodeID);
		coordOut += 0.25 * (1.0 + xi * xiAtNode[i]) * (1.0 + eta * etaAtNode[i]) * coordAtNode;
	}
	return coordOut;

}

// Interpolate z coordinate on top or bottom face from local coordinate of horizontal plane
double MeshDataNonConformingHexaElement::calcZCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const{

	assert( iFace == 4 || iFace == 5 );

	// Array of reference coord xi values for each node
	const double xiAtNode[4] = { -1.0, 1.0, 1.0, -1.0 };

	// Array of reference coord eta values for each node
	const double etaAtNode[4] = { -1.0, -1.0, 1.0, 1.0 };

	double coordOut(0.0);
	for( int i = 0; i < 4; ++i ){
		const int nodeID = getNodeIDGlobalFromElementAndFace(iElem, iFace, i );
		const double coordAtNode = getZCoordinatesOfNodes(nodeID);
		coordOut += 0.25 * (1.0 + xi * xiAtNode[i]) * (1.0 + eta * etaAtNode[i]) * coordAtNode;
	}
	return coordOut;

}

// Calculate length of edges of elements
double MeshDataNonConformingHexaElement::calcEdgeLengthFromElementAndEdge( const int iElem, const int iEdge ) const{

	return calcDistanceOfTwoNodes( getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 0 ), getNodeIDGlobalFromElementAndEdge( iElem, iEdge, 1 ) );

}

// Calculate length of edges of elements on boundary planes
double MeshDataNonConformingHexaElement::calcEdgeLengthFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const{

	const int nodeID0 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 0 );
	const int nodeID1 = getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( iPlane, iElem, iEdge, 1 );

	const double coordZ0 = getZCoordinatesOfNodes( nodeID0 );
	const double coordZ1 = getZCoordinatesOfNodes( nodeID1 );

	if( iPlane == MeshData::ZXMinus || iPlane == MeshData::ZXPlus ){// Z-X plane
		const double coordX0 = getXCoordinatesOfNodes( nodeID0 );
		const double coordX1 = getXCoordinatesOfNodes( nodeID1 );
		return hypot( coordX1 - coordX0, coordZ1 - coordZ0 ); 
	}else if( iPlane == MeshData::YZMinus || iPlane == MeshData::YZPlus ){// Y-Z plane
		const double coordY0 = getYCoordinatesOfNodes( nodeID0 );
		const double coordY1 = getYCoordinatesOfNodes( nodeID1 );
		return hypot( coordY1 - coordY0, coordZ1 - coordZ0 );
	}

	OutputFiles::m_logFile << "Error : Wrong plane ID !! : iPlane = " << iPlane << std::endl;
	exit(1);

	return -1;

}

// Get face index of neighbor element
double MeshDataNonConformingHexaElement::getEdgeLengthX( const int iElem ) const{

	const int node0 = getNodesOfElements( iElem, 0 );
	const int node2 = getNodesOfElements( iElem, 2 );

	return caldDiffXOfTwoNodes( node0, node2 );

}

// Get length of the edges parallel to Y coordinate
double MeshDataNonConformingHexaElement::getEdgeLengthY( const int iElem ) const{

	const int node0 = getNodesOfElements( iElem, 0 );
	const int node2 = getNodesOfElements( iElem, 2 );

	return caldDiffYOfTwoNodes( node0, node2 );

}

// Get face index of neighbor element
int MeshDataNonConformingHexaElement::getFaceIndexOfNeighborElement( const int iFace ) const{

	switch (iFace){
		case 0:
			return 1;
			break;
		case 1:
			return 0;
			break;
		case 2:
			return 3;
			break;
		case 3:
			return 2;
			break;
		case 4:
			return 5;
			break;
		case 5:
			return 4;
			break;
		default:
			OutputFiles::m_logFile << "Error : Face ID is worng !! : iFace = " << iFace << std::endl;
			exit(1);
			break;
	}

	return -1;
}

// Calculate area of face
double MeshDataNonConformingHexaElement::calcAreaOfFace( const int iElem, const int iFace ) const{

	const int nodeID[4] = {
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 0 ),
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 1 ),
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 2 ),
		getNodeIDGlobalFromElementAndFace( iElem, iFace, 3 )
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
	
	const CommonParameters::Vector3D vec2 = {
		getXCoordinatesOfNodes(nodeID[1]) - getXCoordinatesOfNodes(nodeID[3]),
		getYCoordinatesOfNodes(nodeID[1]) - getYCoordinatesOfNodes(nodeID[3]),
		getZCoordinatesOfNodes(nodeID[1]) - getZCoordinatesOfNodes(nodeID[3])
	};

	const CommonParameters::Vector3D vec3 = {
		getXCoordinatesOfNodes(nodeID[2]) - getXCoordinatesOfNodes(nodeID[3]),  
		getYCoordinatesOfNodes(nodeID[2]) - getYCoordinatesOfNodes(nodeID[3]),  
		getZCoordinatesOfNodes(nodeID[2]) - getZCoordinatesOfNodes(nodeID[3])
	};

	CommonParameters::Vector3D vecOut = calcOuterProduct(vec0, vec1);
	double area = 0.5 * fabs( sqrt( vecOut.X * vecOut.X + vecOut.Y * vecOut.Y + vecOut.Z * vecOut.Z ) );

	vecOut = calcOuterProduct(vec2, vec3);
	area += 0.5 * fabs( sqrt( vecOut.X * vecOut.X + vecOut.Y * vecOut.Y + vecOut.Z * vecOut.Z ) );

	return area;

}

// Calculate area of face at bottom of mesh
double MeshDataNonConformingHexaElement::calcAreaOfFaceAtBottomOfMesh( const int iElem ) const{
	const int elemID = getElemBoundaryPlanes(MeshData::XYPlus, iElem);
	const int iFace = getFaceIDLocalFromElementBoundaryPlanes(MeshData::XYPlus, iElem);
	return calcAreaOfFace(elemID, iFace);
}

// Check whether side element-faces are parallel to Z-X or Y-Z plane
void MeshDataNonConformingHexaElement::checkWhetherSideFaceIsParallelToZXOrYZPlane() const{

	const double eps = 1.0e-6;
	const int numElemTotal = getNumElemTotal();
	for( int iElem = 0; iElem < numElemTotal; ++iElem ){
		for( int iFace = 0; iFace < 4; ++iFace ){
			const int nodeIDs[4] = {
				getNodeIDGlobalFromElementAndFace(iElem, iFace, 0 ), 
				getNodeIDGlobalFromElementAndFace(iElem, iFace, 1 ), 
				getNodeIDGlobalFromElementAndFace(iElem, iFace, 2 ), 
				getNodeIDGlobalFromElementAndFace(iElem, iFace, 3 ) };
			if( iFace == 0 || iFace == 1 ){
				// Check whether the side element-face is parallel to the Y-Z plane
				const double x[4] = {
					getXCoordinatesOfNodes(nodeIDs[0]),
					getXCoordinatesOfNodes(nodeIDs[1]),
					getXCoordinatesOfNodes(nodeIDs[2]),
					getXCoordinatesOfNodes(nodeIDs[3]) };
				if( fabs(x[1] - x[0]) > eps || fabs(x[2] - x[0]) > eps || fabs(x[3] - x[0]) > eps ){
					OutputFiles::m_logFile << "Error : Face " << iFace << " of element " << iElem << " is not parallel to the Y-Z plane." << std::endl;
					exit(1);
				}
			}else{
				// Check whether the side element-face is parallel to the Z-X plane
				const double y[4] = {
					getYCoordinatesOfNodes(nodeIDs[0]),
					getYCoordinatesOfNodes(nodeIDs[1]),
					getYCoordinatesOfNodes(nodeIDs[2]),
					getYCoordinatesOfNodes(nodeIDs[3]) };
				if( fabs(y[1] - y[0]) > eps || fabs(y[2] - y[0]) > eps || fabs(y[3] - y[0]) > eps ){
					OutputFiles::m_logFile << "Error : Face " << iFace << " of element " << iElem << " is not parallel to the Z-X plane." << std::endl;
					exit(1);
				}
			}
		}
	}

}

// Check whether the specified point is located in the specified element
bool MeshDataNonConformingHexaElement::isLocatedInTheElement( const double x, const double y, const double z, const int iElem ) const{

	const double eps = 1.0e-9;
	double xi(0.0);
	double eta(0.0);
	double zeta(0.0);
	calcLocalCoordinates(iElem, x, y, z, xi, eta, zeta);
	if( xi > 1.0 + eps || xi < -1.0 - eps || eta > 1.0 + eps || eta < -1.0 - eps || zeta > 1.0 + eps || zeta < -1.0 - eps ){
		return false;
	}else{
		if( xi > 1.0 ){
			xi = 1.0;
		}
		if( eta > 1.0 ){
			eta = 1.0;
		}
		if( zeta > 1.0 ){
			zeta = 1.0;
		}
		if( xi < -1.0 ){
			xi = -1.0;
		}
		if( eta < -1.0 ){
			eta = -1.0;
		}
		if( zeta < -1.0 ){
			zeta = -1.0;
		}
		return true;
	}

}

// Calculate local coordinates
void MeshDataNonConformingHexaElement::calcLocalCoordinates( const int iElem, const double x, const double y, const double z, double& xi, double& eta, double& zeta ) const{

	calcHorizontalLocalCoordinates(iElem, x, y, xi, eta);

	const double zMin = calcZCoordOfPointOnFace(iElem, 4, xi, eta);
	const double zMax = calcZCoordOfPointOnFace(iElem, 5, xi, eta);

	const double zMid = 0.5 * (zMin + zMax);
	zeta = 2.0 * (z - zMid) / (zMax - zMin);

}

// Calculate horizontal local coordinates
void MeshDataNonConformingHexaElement::calcHorizontalLocalCoordinates( const int iElem, const double x, const double y, double& xi, double& eta ) const{

	const int nodeID0 = getNodesOfElements(iElem, 0);
	const int nodeID2 = getNodesOfElements(iElem, 2);

	const double xMin = getXCoordinatesOfNodes(nodeID0);
	const double yMin = getYCoordinatesOfNodes(nodeID0);
	const double xMax = getXCoordinatesOfNodes(nodeID2);
	const double yMax = getYCoordinatesOfNodes(nodeID2);

	const double xMid = 0.5 * (xMin + xMax);
	const double yMid = 0.5 * (yMin + yMax);
	xi   = 2.0 * (x - xMid) / (xMax - xMin);
	eta  = 2.0 * (y - yMid) / (yMax - yMin); 

}

// Calculate determinant of jacobian matrix of the elements
double MeshDataNonConformingHexaElement::calcDeterminantOfJacobianMatrix( const int iElem,  const double xi, const double eta, const double zeta ) const{

	double xCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double yCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	double zCoord[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for( int i = 0; i < 8; ++i ){
		const int nodeID = getNodesOfElements(iElem, i);
		xCoord[i] = getXCoordinatesOfNodes(nodeID);
		yCoord[i] = getYCoordinatesOfNodes(nodeID);
		zCoord[i] = getZCoordinatesOfNodes(nodeID);
	}

	// Zero clear
	CommonParameters::DoubleMatrix3x3 JacobMat;
	JacobMat.comp11 = 0.0;
	JacobMat.comp12 = 0.0;
	JacobMat.comp13 = 0.0;
	JacobMat.comp21 = 0.0;
	JacobMat.comp22 = 0.0;
	JacobMat.comp23 = 0.0;
	JacobMat.comp31 = 0.0;
	JacobMat.comp32 = 0.0;
	JacobMat.comp33 = 0.0;
	for( int i = 0; i < 8; ++i ){
		const double xiNode   = m_xiAtNode[i];
		const double etaNode  = m_etaAtNode[i];
		const double zetaNode = m_zetaAtNode[i];
		const double tmp1 = 0.125 * xiNode   * (1.0 + etaNode  * eta)  * (1.0 + zetaNode * zeta);
		const double tmp2 = 0.125 * etaNode  * (1.0 + zetaNode * zeta) * (1.0 + xiNode * xi);
		const double tmp3 = 0.125 * zetaNode * (1.0 + xiNode * xi)     * (1.0 + etaNode * eta);
		JacobMat.comp11 += tmp1 * xCoord[i];
		JacobMat.comp12 += tmp1 * yCoord[i];
		JacobMat.comp13 += tmp1 * zCoord[i];
		JacobMat.comp21 += tmp2 * xCoord[i];
		JacobMat.comp22 += tmp2 * yCoord[i];
		JacobMat.comp23 += tmp2 * zCoord[i];
		JacobMat.comp31 += tmp3 * xCoord[i];
		JacobMat.comp32 += tmp3 * yCoord[i];
		JacobMat.comp33 += tmp3 * zCoord[i]; 
	}

	const double determinant = JacobMat.comp11 * JacobMat.comp22 * JacobMat.comp33
							 + JacobMat.comp12 * JacobMat.comp23 * JacobMat.comp31
							 + JacobMat.comp13 * JacobMat.comp21 * JacobMat.comp32
							 - JacobMat.comp13 * JacobMat.comp22 * JacobMat.comp31
							 - JacobMat.comp12 * JacobMat.comp21 * JacobMat.comp33
							 - JacobMat.comp11 * JacobMat.comp23 * JacobMat.comp32;

	return determinant;

}
