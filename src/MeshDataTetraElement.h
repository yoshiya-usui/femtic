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
#ifndef DBLDEF_MESHDATA_TETRA_ELEMENT
#define DBLDEF_MESHDATA_TETRA_ELEMENT

#include <vector>
#include "MeshData.h"
#include "CommonParameters.h"

// Class of FEM mesh for tetrahedral element
class MeshDataTetraElement : public MeshData {

public:

	// Default constructer
	MeshDataTetraElement();

	// Destructer
	virtual ~MeshDataTetraElement();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData();

	// Find element including a point
	int findElementIncludingPoint( const double locX, const double locY, const double locZ, CommonParameters::VolumeCoords& localCoord ) const;

	// Find element including a point on the surface of the earth
	int findElementIncludingPointOnSurface( const double locX, const double locY, int& faceID, CommonParameters::AreaCoords& localCoord, const bool useUpperElem,
		const bool modLoc, double& locXMod, double& locYMod ) const;

	// Find element including a point on the Y-Z plane and return element ID of 2D mesh
	int findElementIncludingPointOnYZPlaneAndReturnElemID2D( const int iPlane, const double locY, const double locZ, CommonParameters::AreaCoords& localCoord ) const;

	// Find element including a point on the Z-X plane and return element ID of 2D mesh
	int findElementIncludingPointOnZXPlaneAndReturnElemID2D( const int iPlane, const double locZ, const double locX, CommonParameters::AreaCoords& localCoord ) const;

	// Find elements including dipole on the surface of the earth
	void findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd, std::vector<int>& elements, std::vector<int>& faces,
		std::vector<CommonParameters::AreaCoords>& areaCoordsdStartPoint, std::vector<CommonParameters::AreaCoords>& areaCoordsdEndPoint ) const;

	// Decide whether specified elements share same edges
	virtual bool shareSameEdges( const int elemID1, const int elemID2 ) const;

	// Calculate volume of a specified element
	virtual double calcVolume( const int elemID ) const;

	// Output mesh data to VTK file
	virtual void outputMeshDataToVTK() const;
	
	// Output mesh data to binary file
	virtual void outputMeshDataToBinary() const;

	// Get ID of the nodes of elements belonging to the boundary planes
	virtual int getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const;

	// Get mesh type
	int getMeshType() const;

	// Get local face ID of elements belonging to the boundary planes
	int getFaceIDLocalFromElementBoundaryPlanes( const int iPlane, const int iElem ) const;

	// Get local node ID  from local face ID
	int getNodeIDLocalFromFaceIDLocal( const int iFace, const int num ) const;

	// Get local node ID from local edge ID
	int getNodeIDLocalFromEdgeIDLocal( const int iEdge, const int num ) const;

	// Get global node ID of specified element and edge
	int getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const;

	// Get global node ID of specified element and face
	int getNodeIDGlobalFromElementAndFace( const int iElem, const int iFace, const int num ) const;

	// Get global node ID of specified element belonging to the boundary planes  
	int getNodeIDGlobalFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get X coordinate of node of specified element belonging to the boundary planes  
	double getCoordXFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Y coordinate of node of specified element belonging to the boundary planes  
	double getCoordYFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Z coordinate of node of specified element belonging to the boundary planes  
	double getCoordZFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get global node ID from ID of element belonging to the boundary planes and its edge ID
	// [note] : node ID is outputed as they make a clockwise turn around +X or +Y direction
	int getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge, const int num ) const;

	// Get local edge ID from local face ID
	int getEdgeIDLocalFromFaceIDLocal( const int iFace, const int num ) const;

	// Calculate length of edges of elements
	double calcEdgeLengthFromElementAndEdge( const int iElem, const int iEdge ) const;

	// Calculate length projected on the horizontal plane of edges of elements
	double calcEdgeLengthProjectedOnHorizontalPlaneFromElementAndEdge( const int iElem, const int iEdge ) const;

	// Calculate length of edges of elements on boundary planes
	double calcEdgeLengthFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Calculate horizontal coordinate differences of edges of the elements on boundary planes
	double calcHorizontalCoordDifferenceBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Calculate X coordinate of points on element face
	double calcXCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const;

	// Calculate Y coordinate of points on element face
	double calcYCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const;

	// Calculate Z coordinate of points on element face
	double calcZCoordOfPointOnFace( const int iElem, const int iFace, const CommonParameters::AreaCoords& areaCoord ) const;

	// Calculate volume coordinates from area coordinates
	void calcVolumeCoordFromAreaCoord( const int iFace, const CommonParameters::AreaCoords& areaCoord, CommonParameters::VolumeCoords& volumeCoord ) const;

	// Calculate area projected on X-Y plane with sign of triangle from area coordinate.
	// If rotation direction is plus Z, area value is plus. Otherwise, area value is minus.
	double calcAreaOnXYPlaneWithSignFromAreaCoords( const int elemID, const int faceID, const CommonParameters::AreaCoords& coord0, const CommonParameters::AreaCoords& coord1, const CommonParameters::AreaCoords& coord2 ) const;

	// Calculate area with sign of triangle from area coordinate.
	// If rotation direction is plus Z, area value is plus. Otherwise, area value is minus.
	double calcAreaWithSignFromAreaCoords( const int elemID, const int faceID, const CommonParameters::AreaCoords& coord0, const CommonParameters::AreaCoords& coord1, const CommonParameters::AreaCoords& coord2 ) const;

	// Calculate area of face
	double calcAreaOfFace( const int iElem, const int iFace ) const;

	// Calculate area of face at bottom of mesh
	double calcAreaOfFaceAtBottomOfMesh( const int iElem ) const;

#ifdef _DEBUG_WRITE
	// Function for debug
	void testFuction() const;
#endif

private:

	// Copy constructer
	MeshDataTetraElement(const MeshDataTetraElement& rhs);

	// Copy assignment operator
	MeshDataTetraElement& operator=(const MeshDataTetraElement& rhs);

	static const double m_eps;

	// Array of faces of elements belonging to the boundary planes
	//   m_facesOfElementsBoundaryPlanes[0] : Y-Z Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[1] : Y-Z Plane ( Plus Side  )
	//   m_facesOfElementsBoundaryPlanes[2] : Z-X Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[3] : Z-X Plane ( Plus Side  )
	//   m_facesOfElementsBoundaryPlanes[4] : X-Y Plane ( Minus Side )
	//   m_facesOfElementsBoundaryPlanes[5] : X-Y Plane ( Plus Side  )
	int* m_facesOfElementsBoundaryPlanes[6];

	// Number of elements belonging to the land surface
	int m_numElemOnLandSurface;

	// Array of elements belonging to the land surface
	int* m_elemOnLandSurface;

	// Array of faces belonging to the land surface
	int* m_faceLandSurface;

	// Array converting from face ID to node ID
	int m_faceID2NodeID[4][3];

	// Array converting from face ID to edge ID
	int m_faceID2EdgeID[4][3];

	// Array converting from edge ID to node ID
	int m_edgeID2NodeID[6][2];

	// Function determine if two segments intersect or not
	bool locateLeftOfSegmentOnLandSurface( const CommonParameters::locationXY& point, 
		const CommonParameters::locationXY& startPointOfSegment, const CommonParameters::locationXY& endPointOfSegment ) const;

	// Function determine if then inputed point locate at the left of the segment on the surface of the upper element
	bool locateLeftOfSegmentOnSeaSurface( const CommonParameters::locationXY& point, 
		const CommonParameters::locationXY& startPointOfSegment, const CommonParameters::locationXY& endPointOfSegment ) const;

	// Function determine if then inputed point locate at the left of the segment on the Y-Z plane of boundary
	bool locateLeftOfSegmentOnYZPlaneOfBoundary( const int iPlane, const int iElem, const int iEdge, const CommonParameters::locationYZ& point ) const;

	// Function determine if then inputed point locate at the left of the segment on the Z-X plane of boundary
	bool locateLeftOfSegmentOnZXPlaneOfBoundary( const int iPlane, const int iElem, const int iEdge, const CommonParameters::locationZX& point ) const;

	// Calculate volume of tetrahedron
	double calcVolume( const CommonParameters::locationXYZ& point1, const CommonParameters::locationXYZ& point2,
		const CommonParameters::locationXYZ& point3, const CommonParameters::locationXYZ& point4 ) const;

	// Calculate volume coordinates of point on the land surface
	void calcVolumeCoordsOfPointOnLandSurface( const int elemID, const int faceID, const CommonParameters::locationXY& pointCoord, CommonParameters::VolumeCoords& coords ) const;

	// Calculate volume coordinates of the nputed point
	void calcVolumeCoordsOfPoint( const int elemID, const CommonParameters::locationXYZ& pointCoord, CommonParameters::VolumeCoords& coords ) const;

	// Calculate area coordinates of point on the land surface
	void calcAreaCoordsOfPointOnLandSurface( const int elemID, const int faceID, const CommonParameters::locationXY& pointCoord, CommonParameters::AreaCoords& coords ) const;

	// Calculate area of triangle from two dimensinal coordinates
	double calcArea( const CommonParameters::CoordPair& point1, const CommonParameters::CoordPair& point2, const CommonParameters::CoordPair& point3 ) const;

	// Calculate area coordinates of the specified point on the Y-Z plane of boundary
	void calcAreaCoordsOfPointOnYZPlaneOfBoundary( const int iPlane, const int iElem, const CommonParameters::CoordPair& point, CommonParameters::AreaCoords& coords ) const;

	// Calculate area coordinates of the specified point on the Z-X plane of boundary
	void calcAreaCoordsOfPointOnZXPlaneOfBoundary( const int iPlane, const int iElem, const CommonParameters::CoordPair& point, CommonParameters::AreaCoords& coords ) const;

	// Decide whether specified point locate inside of face
	bool locateInsideOfFace( const int elemID, const int faceID, const CommonParameters::locationXYZ& loc ) const;

};
#endif
