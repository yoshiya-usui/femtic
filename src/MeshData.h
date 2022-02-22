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
#ifndef DBLDEF_MESHDATA
#define DBLDEF_MESHDATA

#include <vector>
#include "CommonParameters.h"

// Class of FEM mesh for brick element
class MeshData{

public:

	enum BoundaryPlanes{
		YZMinus = 0,
		YZPlus,
		ZXMinus,
		ZXPlus,
		XYMinus,
		XYPlus,
	};

	enum MeshType{
		HEXA = 0,
		TETRA,
		NONCONFORMING_HEXA,
	};

	struct coordinateValue{
		double X;
		double Y;
		double Z;
	};

	struct coordinateValueXY{
		double X;
		double Y;
	};

	// Constructer
	MeshData();

	// Destructer
	virtual ~MeshData();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData() = 0;

	// Get tolal number of elements
	int getNumElemTotal() const;

	// Get tolal number of nodes
	int getNumNodeTotal() const;

	// Get total number of elements belonging to the boundary planes
	int getNumElemOnBoundaryPlanes( const int iPlane ) const;

	// Get X coordinates of node
	double getXCoordinatesOfNodes( const int iNode ) const;

	// Get Y coordinates of node
	double getYCoordinatesOfNodes( const int iNode ) const;

	// Get Z coordinates of node
	double getZCoordinatesOfNodes( const int iNode ) const;

	// Get ID of the Node composing specified element
	int getNodesOfElements( const int iElem, const int iNode ) const;

	// Get ID of the element belonging to the boundary planes
	int getElemBoundaryPlanes( const int iPlane, const int iElem ) const;

	// Get ID of neighbor Elements
	int getIDOfNeighborElement( const int iElem, const int num ) const;

	// Get number of neighbor elements of one element
	int getNumNeighborElement() const;

	// Get mesh type
	virtual int getMeshType() const = 0;

	// Calculate distance of two nodes
	double calcDistanceOfTwoNodes( const int nodeID0,  const int nodeID1 ) const;

	// Calculate horizontal distance of two nodes
	double calcHorizontalDistanceOfTwoNodes( const int nodeID0,  const int nodeID1 ) const;

	// Calculate distance of two nodes along X direction
	double caldDiffXOfTwoNodes( const int nodeID0,  const int nodeID1 ) const;

	// Calculate distance of two nodes along Y direction
	double caldDiffYOfTwoNodes( const int nodeID0,  const int nodeID1 ) const;

	// Calculate distance of two nodes along Z direction
	double caldDiffZOfTwoNodes( const int nodeID0,  const int nodeID1 ) const;

	// Decide whether specified elements share same nodes
	virtual bool shareSameNodes( const int elemID1, const int elemID2 ) const;

	// Calculate coordinate of the center of a specified element
	virtual CommonParameters::locationXYZ getCenterCoord( const int iElem ) const;

	// Calculate difference of the centers of the specified two element
	CommonParameters::locationXYZ calDiffOfCenters( const int iElem1, const int iElem2 ) const;

	// Decide whether specified elements share same edges
	virtual bool shareSameEdges( const int elemID1, const int elemID2 ) const = 0;

	// Get ID of the nodes of elements belonging to the boundary planes
	virtual int getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const = 0;

	// Calculate volume of a specified element
	virtual double calcVolume( const int elemID ) const = 0;

	// Output mesh data to VTK file
	virtual void outputMeshDataToVTK() const = 0;

	// Output mesh data to binary file
	virtual void outputMeshDataToBinary() const = 0;

	// Calculate area of face
	virtual double calcAreaOfFace( const int iElem, const int iFace ) const = 0;

	// Calculate area of face at bottom of mesh
	virtual double calcAreaOfFaceAtBottomOfMesh( const int iElem ) const = 0;

protected:

	// Copy constructer
	MeshData(const MeshData& rhs);

	// Copy assignment operator
	MeshData& operator=(const MeshData& rhs);

	// Total number of elements
	int m_numElemTotal;

	// Total number of nodes
	int m_numNodeTotal;

	// Number of nodes belonging to one element
	int m_numNodeOneElement;

	// Number of edges belonging to one element
	int m_numEdgeOneElement;

	// Number of nodes on a face of one element
	int m_numNodeOnFaceOneElement;

	// Number of neighbor elements of one element
	int m_numNeighborElement;

	// Total number of elements belonging to the boundary planes
	int m_numElemOnBoundaryPlanes[6];

	// Array of the X coordinates of nodes
	double* m_xCoordinatesOfNodes;

	// Array of the Y coordinates of nodes
	double* m_yCoordinatesOfNodes;

	// Array of the Z coordinates of nodes
	double* m_zCoordinatesOfNodes;

	// Array of IDs of neighbor Elements
	int* m_neighborElements;

	// Array of nodes composing each element
	int* m_nodesOfElements;

	// Array of elements belonging to the boundary planes
	//   m_elemBoundaryPlane[0] : Y-Z Plane ( Minus Side )
	//   m_elemBoundaryPlane[1] : Y-Z Plane ( Plus Side  )
	//   m_elemBoundaryPlane[2] : Z-X Plane ( Minus Side )
	//   m_elemBoundaryPlane[3] : Z-X Plane ( Plus Side  )
	//   m_elemBoundaryPlane[4] : X-Y Plane ( Minus Side )
	//   m_elemBoundaryPlane[5] : X-Y Plane ( Plus Side  )
	int* m_elemBoundaryPlanes[6];

	// Calculate distanceof two points
	double calcDistance( const CommonParameters::locationXY& point0,  const CommonParameters::locationXY& point1 ) const;

	// Function determine if the 1st segment contains the 2nd segment
	bool does1stSegmentContain2ndSegment( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
		const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const;

	// Function determine if two segments intersect or not
	bool intersectTwoSegments( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
		const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const;

	// Function determine if two lines overlap or not
	bool overlapTwoLines( const CommonParameters::locationXY& coord1stLine1, const CommonParameters::locationXY& coord1stLine2,
		const CommonParameters::locationXY& coord2ndLine1, const CommonParameters::locationXY& coord2ndLine2 ) const;

	// Function determine if two segments overlap or not
	bool overlapTwoSegments( const CommonParameters::locationXY& startPointOf1stSegment, const CommonParameters::locationXY& endPointOf1stSegment,
		const CommonParameters::locationXY& startPointOf2ndSegment, const CommonParameters::locationXY& endPointOf2ndSegment ) const;

	// Calculate inner product of two vectors
	double calcInnerProduct2D( const CommonParameters::locationXY& startCoordOf1stVec, const CommonParameters::locationXY& endCoordOf1stVec,
		const CommonParameters::locationXY& startCoordOf2ndVec, const CommonParameters::locationXY& endCoordOf2ndVec) const;

	// Calculate coordinates of intersection point of two lines
	void calcCoordOfIntersectionPointOfTwoLines( const CommonParameters::locationXY& coord1stLine1, const CommonParameters::locationXY& coord1stLine2,
		const CommonParameters::locationXY& coord2ndLine1, const CommonParameters::locationXY& coord2ndLine2, CommonParameters::locationXY& coordIntersectionPoint) const;

};

#endif
