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
#include "Util.h"

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
bool MeshData::does1stSegmentContain2ndSegment( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
	const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const{

	const double dx1st = endPointOf1stSegment.X - startPointOf1stSegment.X;
	const double dy1st = endPointOf1stSegment.Y - startPointOf1stSegment.Y;
	const double outerProduct1 = dx1st * ( startPointOf2ndSegment.Y - startPointOf1stSegment.Y ) - dy1st * ( startPointOf2ndSegment.X - startPointOf1stSegment.X );
	const double outerProduct2 = dx1st * ( endPointOf2ndSegment.Y - startPointOf1stSegment.Y ) - dy1st * ( endPointOf2ndSegment.X - startPointOf1stSegment.X );

	const double eps = 1.0e-9;
	if( fabs(outerProduct1) < eps && fabs(outerProduct2) < eps ){
		// Two segmenta are located along a same line
		const double innerProduct1 = ( startPointOf2ndSegment.X - startPointOf1stSegment.X ) * ( endPointOf1stSegment.X - startPointOf1stSegment.X )
								+ ( startPointOf2ndSegment.Y - startPointOf1stSegment.Y ) * ( endPointOf1stSegment.Y - startPointOf1stSegment.Y );
		const double innerProduct2 = ( endPointOf2ndSegment.X - startPointOf1stSegment.X ) * ( endPointOf1stSegment.X - startPointOf1stSegment.X )
								+ ( endPointOf2ndSegment.Y - startPointOf1stSegment.Y ) * ( endPointOf1stSegment.Y - startPointOf1stSegment.Y );
		const double denominator = pow(endPointOf1stSegment.X - startPointOf1stSegment.X, 2) + pow(endPointOf1stSegment.Y - startPointOf1stSegment.Y, 2);
		if( innerProduct1 / denominator > - eps && innerProduct1 / denominator < 1.0 + eps &&
			innerProduct2 / denominator > - eps && innerProduct2 / denominator < 1.0 + eps ){
			return true;
		}
	}

	return false;
}

// Function determine if two segments intersect or not
bool MeshData::intersectTwoSegments( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
	const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const{

	const double val1 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * ( startPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y - startPointOf2ndSegment.Y );
	const double val2 = ( endPointOf1stSegment.Y - startPointOf1stSegment.Y ) * (   endPointOf2ndSegment.X - startPointOf1stSegment.X ) + ( endPointOf1stSegment.X - startPointOf1stSegment.X ) * ( startPointOf1stSegment.Y -   endPointOf2ndSegment.Y );

	const double eps = 1.0e-9;

	if( calcDistance(startPointOf1stSegment, startPointOf2ndSegment) < eps ||
		calcDistance(endPointOf1stSegment, startPointOf2ndSegment) < eps ||
		calcDistance(startPointOf1stSegment, endPointOf2ndSegment) < eps ||
		calcDistance(endPointOf1stSegment, endPointOf2ndSegment) < eps ){
		// Two segment share a node
		return true;
	}

	if( val1*val2 < eps ){
		const double val3 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * ( startPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y - startPointOf1stSegment.Y );
		const double val4 = ( endPointOf2ndSegment.Y - startPointOf2ndSegment.Y ) * (   endPointOf1stSegment.X - startPointOf2ndSegment.X ) + ( endPointOf2ndSegment.X - startPointOf2ndSegment.X ) * ( startPointOf2ndSegment.Y -   endPointOf1stSegment.Y );
		if( fabs(val1*val2) < eps && fabs(val3*val4) < eps ){
			return false;
		}else if( val3*val4 < eps ){
			return true;
		}
		//if( val3*val4 < eps ){
		//	return true;
		//}
	}

	return false;

}

// Function determine if two lines overlap or not
bool MeshData::overlapTwoLines( const CommonParameters::locationXY& coord1stLine1, const CommonParameters::locationXY& coord1stLine2,
	const CommonParameters::locationXY& coord2ndLine1, const CommonParameters::locationXY& coord2ndLine2 ) const{

	const double eps = 1.0e-12;
	if( fabs( coord1stLine2.X - coord1stLine1.X ) < eps ){// If the first line is parallel to Y direction
		if( fabs( coord2ndLine2.X - coord2ndLine1.X ) < eps ){// If the second line is also parallel to Y direction
			if( fabs( coord1stLine1.X - coord2ndLine1.X ) < eps && fabs( coord1stLine2.X - coord2ndLine2.X ) < eps ){
				return true;
			}else{
				return false;
			}
		}else{// If the second line is not parallel to Y direction
			return false;
		}
	}

	if( fabs( coord1stLine2.Y - coord1stLine1.Y ) < eps ){// If the first line is parallel to X direction
		if( fabs( coord2ndLine2.Y - coord2ndLine1.Y ) < eps ){// If the second line is also parallel to X direction
			if( fabs( coord1stLine1.Y - coord2ndLine1.Y ) < eps && fabs( coord1stLine2.Y - coord2ndLine2.Y ) < eps ){
				return true;
			}else{
				return false;
			}
		}else{// If the second line is not parallel to X direction
			return false;
		}
	}

	const double val1 = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X ) - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y );

	const double val2 = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X )*coord1stLine1.X
					  - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y )*coord2ndLine1.X  
		  			  + ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.X - coord2ndLine1.X )*( coord2ndLine1.Y - coord1stLine1.Y ); 

	if( fabs(val1) < eps && fabs(val2) < eps ){
		return true;
	}

	return false;

}

// Function determine if two segments overlap or not
bool MeshData::overlapTwoSegments( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
	const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const{

	if( !overlapTwoLines( startPointOf1stSegment, endPointOf1stSegment, startPointOf2ndSegment, endPointOf2ndSegment ) ){
		// Two lines don't overlap
		return false;
	}
	
	const double innerProduct1 = calcInnerProduct2D( startPointOf1stSegment, endPointOf1stSegment, startPointOf1stSegment, endPointOf1stSegment );
	const double innerProduct2 = calcInnerProduct2D( startPointOf1stSegment, endPointOf1stSegment, startPointOf1stSegment, startPointOf2ndSegment );
	const double innerProduct3 = calcInnerProduct2D( startPointOf1stSegment, endPointOf1stSegment, startPointOf1stSegment, endPointOf2ndSegment );

	if( ( innerProduct2 < 0.0 && innerProduct3 < 0.0 ) || ( innerProduct2 > innerProduct1 && innerProduct3 > innerProduct1 ) ){
		return false;
	}

	return true;

}

// Calculate inner product of two vectors
double MeshData::calcInnerProduct2D( const CommonParameters::locationXY& startCoordOf1stVec, const CommonParameters::locationXY& endCoordOf1stVec,
	const CommonParameters::locationXY& startCoordOf2ndVec, const CommonParameters::locationXY& endCoordOf2ndVec) const{

	CommonParameters::Vector3D vec1 = { endCoordOf1stVec.X - startCoordOf1stVec.X, endCoordOf1stVec.Y - startCoordOf1stVec.Y, 0.0 };
	CommonParameters::Vector3D vec2 = { endCoordOf2ndVec.X - startCoordOf2ndVec.X, endCoordOf2ndVec.Y - startCoordOf2ndVec.Y, 0.0 };

	return calcInnerProduct( vec1, vec2 );

}

// Calculate coordinates of intersection point of two lines
void MeshData::calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::locationXY& coord1stLine1, const CommonParameters::locationXY& coord1stLine2,
	const CommonParameters::locationXY& coord2ndLine1, const CommonParameters::locationXY& coord2ndLine2, CommonParameters::locationXY& coordIntersectionPoint) const{

	const double eps = 1.0e-9;
	if( calcDistance(coord1stLine1, coord2ndLine1) < eps || calcDistance(coord1stLine1, coord2ndLine2) < eps ){
		coordIntersectionPoint = coord1stLine1;
		return;
	}
	if( calcDistance(coord1stLine2, coord2ndLine1) < eps || calcDistance(coord1stLine2, coord2ndLine2) < eps ){
		coordIntersectionPoint = coord1stLine2;
		return;
	}

	const double temp1 = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X ) - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y );
	if( fabs( temp1 ) < eps ){
		OutputFiles::m_logFile << " Error : Divide by zero in calculating X coordinate of intersection point of two lines !!" << std::endl;
		exit(1);
	}

	coordIntersectionPoint.X = ( coord1stLine2.Y - coord1stLine1.Y )*( coord2ndLine2.X - coord2ndLine1.X )*coord1stLine1.X
							 - ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.Y - coord2ndLine1.Y )*coord2ndLine1.X  
							 + ( coord1stLine2.X - coord1stLine1.X )*( coord2ndLine2.X - coord2ndLine1.X )*( coord2ndLine1.Y - coord1stLine1.Y );
	coordIntersectionPoint.X /= temp1;

	const double temp2 = coord1stLine2.X - coord1stLine1.X;
	const double temp3 = coord2ndLine2.X - coord2ndLine1.X;

	if( fabs( temp2 ) < 1.0e-8 && fabs( temp3 ) < 1.0e-8 ){
		OutputFiles::m_logFile << " Error : Divide by zero in calculating Y coordinate of intersection point of two lines !!" << std::endl;
		exit(1);
	}

	if( fabs( temp2 ) > fabs( temp3 ) ){
		coordIntersectionPoint.Y = ( coord1stLine2.Y - coord1stLine1.Y )/temp2*( coordIntersectionPoint.X - coord1stLine1.X ) + coord1stLine1.Y;
	}else{
		coordIntersectionPoint.Y = ( coord2ndLine2.Y - coord2ndLine1.Y )/temp3*( coordIntersectionPoint.X - coord2ndLine1.X ) + coord2ndLine1.Y;
	}

	return;

}
