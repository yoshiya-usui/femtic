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

#include "MeshData.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"
#include "OutputFiles.h"

// Constructer
MeshData::MeshData():
	m_numElemTotal(NULL),
	m_numNodeTotal(NULL),
	m_numNodeOneElement(8),
	m_numEdgeOneElement(12),
	m_numNodeOnFaceOneElement(4),
	m_numNeighborElement(6),
	m_xCoordinatesOfNodes(NULL),
	m_yCoordinatesOfNodes(NULL),
	m_zCoordinatesOfNodes(NULL),
	m_neighborElements(NULL),
	m_nodesOfElements(NULL)
{

	for ( int i = 0; i < 6; ++i ){
		m_numElemOnBoundaryPlanes[i] = NULL;
	}
	
	for ( int i = 0; i < 6; ++i ){
		m_elemBoundaryPlanes[i] = NULL;
	}

}

// Destructer
MeshData::~MeshData(){

	if( m_xCoordinatesOfNodes != NULL){
		delete[] m_xCoordinatesOfNodes;
		m_xCoordinatesOfNodes = NULL;
	}

	if( m_yCoordinatesOfNodes != NULL){
		delete[] m_yCoordinatesOfNodes;
		m_yCoordinatesOfNodes = NULL;
	}

	if( m_zCoordinatesOfNodes != NULL){
		delete[] m_zCoordinatesOfNodes;
		m_zCoordinatesOfNodes = NULL;
	}

	if( m_neighborElements != NULL){
		delete[] m_neighborElements;
		m_neighborElements = NULL;
	}

	if( m_nodesOfElements != NULL){
		delete[] m_nodesOfElements;
		m_nodesOfElements = NULL;
	}

	for ( int i = 0; i < 6; ++i ){
		if( m_elemBoundaryPlanes[i] != NULL){
			delete[] m_elemBoundaryPlanes[i];
			m_elemBoundaryPlanes[i] = NULL;
		}
	}

}

// Copy constructer
MeshData::MeshData(const MeshData& rhs){
	std::cerr << "Error : Copy constructer of the class MeshData is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
MeshData& MeshData::operator=(const MeshData& rhs){
	std::cerr << "Error : Assignment operator of the class MeshData is not implemented." << std::endl;
	exit(1);
}

// Get tolal number of elements
int MeshData::getNumElemTotal() const{
	return m_numElemTotal;
}

// Get tolal number of nodes
int MeshData::getNumNodeTotal() const{
	return m_numNodeTotal;
}

// Get total number of elements belonging to the boundary planes
int MeshData::getNumElemOnBoundaryPlanes( const int iPlane ) const{

	//if( iPlane < 0 || iPlane >= 6 ){
	//	OutputFiles::m_logFile << "iPlane is out of range in getNumElemOnBoundaryPlanes. iPlane = " << iPlane << std::endl;
	//	exit(1);
	//}
	assert( iPlane >= 0 );
	assert( iPlane < 6);

	return m_numElemOnBoundaryPlanes[iPlane];
}

// Get X coordinates of node
double MeshData::getXCoordinatesOfNodes( const int iNode ) const{

	//if( iNode < 0 || iNode >= m_numNodeTotal ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getXCoordinatesOfNodes. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iNode >= 0 );
	assert( iNode < m_numNodeTotal );

	return m_xCoordinatesOfNodes[iNode];

}

// Get Y coordinates of node
double MeshData::getYCoordinatesOfNodes( const int iNode ) const{

	//if( iNode < 0 || iNode >= m_numNodeTotal ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getYCoordinatesOfNodes. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iNode >= 0 );
	assert( iNode < m_numNodeTotal );

	return m_yCoordinatesOfNodes[iNode];
}

// Get Z coordinates of node
double MeshData::getZCoordinatesOfNodes( const int iNode ) const{

	//if( iNode < 0 || iNode >= m_numNodeTotal ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getZCoordinatesOfNodes. iNode = " << iNode << std::endl;
	//	exit(1);
	//}
	assert( iNode >= 0 );
	assert( iNode < m_numNodeTotal );

	return m_zCoordinatesOfNodes[iNode];
}

// Get ID of the Node composing specified element
int MeshData::getNodesOfElements( const int iElem, const int iNode ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getNodesOfElements. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	//if( iNode < 0 || iNode >= m_numNodeOneElement ){
	//	OutputFiles::m_logFile << " Error : iNode is out of range in getNodesOfElements. iNode = " << iNode << std::endl;
	//	exit(1);
	//}

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( iNode >= 0 );
	assert( iNode < m_numNodeOneElement );

	return m_nodesOfElements[ m_numNodeOneElement*iElem + iNode ];
}

// Get ID of the element belonging to the boundary planes
int MeshData::getElemBoundaryPlanes( const int iPlane, const int iElem ) const{

	//if( iElem < 0 || iElem >= getNumElemOnBoundaryPlanes( iPlane ) ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getElemBoundaryPlanes. iElem = " << iElem << std::endl;
	//	exit(1);
	//}

	assert( iElem >= 0 );
	assert( iElem < getNumElemOnBoundaryPlanes( iPlane ) );

	return m_elemBoundaryPlanes[iPlane][iElem];
}

// Get ID of neighbor Elements
int MeshData::getIDOfNeighborElement( const int iElem, const int num ) const{

	//if( iElem < 0 || iElem >= m_numElemTotal ){
	//	OutputFiles::m_logFile << " Error : iElem is out of range in getIDOfNeighborElement. iElem = " << iElem << std::endl;
	//	exit(1);
	//}
	//if( num < 0 || num >= m_numNeighborElement ){
	//	OutputFiles::m_logFile << " Error : num is out of range in getIDOfNeighborElement. num = " << num << std::endl;
	//	exit(1);
	//}

	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal );
	assert( num >= 0 );
	assert( num < m_numNeighborElement );

	return m_neighborElements[ iElem * m_numNeighborElement + num ];

}

// Get number of neighbor elements of one element
int MeshData::getNumNeighborElement() const{

	return m_numNeighborElement;

}

// Calculate distance of two nodes
double MeshData::calcDistanceOfTwoNodes( const int nodeID0,  const int nodeID1 ) const{

	double val = pow( getXCoordinatesOfNodes(nodeID1) - getXCoordinatesOfNodes(nodeID0) , 2.0 ); 
	val += pow( getYCoordinatesOfNodes(nodeID1) - getYCoordinatesOfNodes(nodeID0) , 2.0 ); 
	val += pow( getZCoordinatesOfNodes(nodeID1) - getZCoordinatesOfNodes(nodeID0) , 2.0 ); 

	return sqrt( val );
	
}

// Calculate horizotal distance of two nodes
double MeshData::calcHorizontalDistanceOfTwoNodes( const int nodeID0,  const int nodeID1 ) const{

	double val = pow( getXCoordinatesOfNodes(nodeID1) - getXCoordinatesOfNodes(nodeID0) , 2.0 ); 
	val += pow( getYCoordinatesOfNodes(nodeID1) - getYCoordinatesOfNodes(nodeID0) , 2.0 ); 

	return sqrt( val );
	
}

// Calculate distance of two nodes along X direction
double MeshData::caldDiffXOfTwoNodes( const int nodeID0,  const int nodeID1 ) const{

	return getXCoordinatesOfNodes(nodeID1) - getXCoordinatesOfNodes(nodeID0);

}

// Calculate distance of two nodes along Y direction
double MeshData::caldDiffYOfTwoNodes( const int nodeID0,  const int nodeID1 ) const{

	return getYCoordinatesOfNodes(nodeID1) - getYCoordinatesOfNodes(nodeID0);

}

// Calculate distance of two nodes along Z direction
double MeshData::caldDiffZOfTwoNodes( const int nodeID0,  const int nodeID1 ) const{

	return getZCoordinatesOfNodes(nodeID1) - getZCoordinatesOfNodes(nodeID0);

}

// Decide whether specified elements share same nodes
bool MeshData::shareSameNodes( const int elemID1, const int elemID2 ) const{

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

	for( int iNod1 = 0; iNod1 < m_numNodeOneElement; ++iNod1 ){
		
		const int nodeID = getNodesOfElements( elemID1, iNod1 );

		for( int iNod2 = 0; iNod2 < m_numNodeOneElement; ++iNod2 ){

			if( nodeID == getNodesOfElements( elemID2, iNod2 ) ){
				return true;
			}

		}

	}

	return false;

}

// Calculate coordinate of the center of a specified element
CommonParameters::locationXYZ MeshData::getCenterCoord( const int iElem ) const{
	assert( iElem >= 0 );
	assert( iElem < m_numElemTotal);

	CommonParameters::locationXYZ val = { 0.0, 0.0, 0.0 };

	for( int i = 0; i < m_numNodeOneElement; ++i ){
		val.X += getXCoordinatesOfNodes( getNodesOfElements( iElem, i ) );
		val.Y += getYCoordinatesOfNodes( getNodesOfElements( iElem, i ) );
		val.Z += getZCoordinatesOfNodes( getNodesOfElements( iElem, i ) );
	}

	val.X /= static_cast<double>(m_numNodeOneElement);
	val.Y /= static_cast<double>(m_numNodeOneElement);
	val.Z /= static_cast<double>(m_numNodeOneElement);

	return val;
}

// Calculate difference of the centers of the specified two element
CommonParameters::locationXYZ MeshData::calDiffOfCenters( const int iElem1, const int iElem2 ) const{
	const CommonParameters::locationXYZ coord1 = getCenterCoord(iElem1);
	const CommonParameters::locationXYZ coord2 = getCenterCoord(iElem2);
	const CommonParameters::locationXYZ coordOut = { coord1.X - coord2.X, coord1.Y - coord2.Y, coord1.Z - coord2.Z };
	return coordOut;
}

// Calculate distanceof two points
double MeshData::calcDistance( const CommonParameters::locationXY& point0,  const CommonParameters::locationXY& point1 ) const{

	return hypot( point1.X - point0.X, point1.Y - point0.Y );

}

//// Decide whether specified elements share same edges
//bool MeshData::shareSameEdges( const int elemID1, const int elemID2 ) const{
//
//	if( elemID1 < 0 || elemID1 >= m_numElemTotal ){
//		OutputFiles::m_logFile << " Error : elemID1 is out of range in shareSameNodes. elemID1 = " << elemID1 << std::endl;
//		exit(1);
//	}
//
//	if( elemID2 < 0 || elemID2 >= m_numElemTotal ){
//		OutputFiles::m_logFile << " Error : elemID2 is out of range in shareSameNodes. elemID2 = " << elemID2 << std::endl;
//		exit(1);
//	}
//
//
//
//}
