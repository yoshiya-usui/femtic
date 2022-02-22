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
#ifndef DBLDEF_MESHDATA_NONCONFORMING_HEXA_ELEMENT
#define DBLDEF_MESHDATA_NONCONFORMING_HEXA_ELEMENT

#include <vector>
#include "MeshData.h"

// Class of FEM mesh for brick element
class MeshDataNonConformingHexaElement : public MeshData {

public:

	// Constructer
	MeshDataNonConformingHexaElement();

	// Destructer
	virtual ~MeshDataNonConformingHexaElement();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData();

	// Find element including a point
	int findElementIncludingPoint( const double locX, const double locY, const double locZ, double& xi, double& eta, double& zeta ) const;

	// Find elements including point on the surface of the earth
	int findElementIncludingPointOnSurface( const double locX, const double locY, int& faceID, double& xi, double& eta, double& zeta,
		const bool useUpperElem, const bool modLoc, double& locXMod, double& locYMod ) const;

	// Find element including a point on the Y-Z plane and return element ID of 2D mesh
	void findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd,
		std::vector<int>& elements, std::vector<double>& localCoordXStartPoint, std::vector<double>& localCoordYStartPoint,	std::vector<double>& localCoordXEndPoint, std::vector<double>& localCoordYEndPoint ) const;

	// Find element including a point on the Y-Z plane and return element ID of 2D mesh
	int findElementIncludingPointOnYZPlaneAndReturnElemID2D( const int iPlane, const double locY, const double locZ, double& xi, double& eta ) const;

	// Find element including a point on the Z-X plane and return element ID of 2D mesh
	int findElementIncludingPointOnZXPlaneAndReturnElemID2D( const int iPlane, const double locX, const double locZ, double& xi, double& eta ) const;

	// Get mesh type
	int getMeshType() const;

	// Get ID of a neighbor element
	int getIDOfNeighborElement( const int iElem, const int iFace, const int num ) const;

	// Get number of neighbor elements for an element-face
	int getNumNeighborElement( const int iElem, const int iFace ) const;

	// Get flag specifing whether an element face has slave faces
	bool faceSlaveElements( const int iElem, const int iFace ) const;

	// Get flag specifing whether an element face is outer boundary
	bool isOuterBoundary( const int iElem, const int iFace ) const;

	// Get local face ID of elements belonging to the boundary planes
	int getFaceIDLocalFromElementBoundaryPlanes( const int iPlane, const int iElem ) const;

	// Get global node ID of specified element and edge
	int getNodeIDGlobalFromElementAndEdge( const int iElem, const int iEdge, const int num ) const;

	// Get global node ID of specified element and face
	int getNodeIDGlobalFromElementAndFace( const int iElem, const int iFace, const int num ) const;

	// Get global node ID of specified element belonging to the boundary planes  
	int getNodeIDGlobalFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get global node ID from ID of element belonging to the boundary planes and its edge index
	int getNodeIDGlobalFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge, const int num ) const;

	// Get X coordinate of node of specified element belonging to the boundary planes  
	double getCoordXFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Y coordinate of node of specified element belonging to the boundary planes  
	double getCoordYFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get Z coordinate of node of specified element belonging to the boundary planes  
	double getCoordZFromElementBoundaryPlanes( const int iPlane, const int iElem, const int num ) const;

	// Get local edge ID from local face ID
	int getEdgeIDLocalFromFaceIDLocal( const int iFace, const int num ) const;

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

	// Calculate horizontal coordinate differences of edges of the elements on boundary planes
	double calcHorizontalCoordDifferenceBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Interpolate x coordinate on top or bottom face from local coordinate of horizontal plane
	double calcXCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Interpolate y coordinate on top or bottom face from local coordinate of horizontal plane
	double calcYCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Interpolate z coordinate on top or bottom face from local coordinate of horizontal plane
	double calcZCoordOfPointOnFace( const int iElem, const int iFace, const double xi, const double eta ) const;

	// Calculate length of edges of elements
	double calcEdgeLengthFromElementAndEdge( const int iElem, const int iEdge ) const;

	// Calculate length of edges of elements on boundary planes
	double calcEdgeLengthFromElementAndEdgeBoundaryPlanes( const int iPlane, const int iElem, const int iEdge ) const;

	// Get face index of neighbor element
	double getEdgeLengthX( const int iElem ) const;

	// Get length of the edges parallel to Y coordinate
	double getEdgeLengthY( const int iElem ) const;

	// Get face index of neighbor element
	int getFaceIndexOfNeighborElement( const int iFace ) const;

	// Calculate area of face
	double calcAreaOfFace( const int iElem, const int iFace ) const;

	// Calculate area of face at bottom of mesh
	double calcAreaOfFaceAtBottomOfMesh( const int iElem ) const;

private:

	// Copy constructer
	MeshDataNonConformingHexaElement(const MeshDataNonConformingHexaElement& rhs);

	// Copy assignment operator
	MeshDataNonConformingHexaElement& operator=(const MeshDataNonConformingHexaElement& rhs);

	// Array of IDs of neighbor Elements
	std::vector<int>* m_neighborElementsForNonConformingHexa;

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
	int m_faceID2NodeID[6][4];

	// Array converting from face ID to edge ID
	int m_faceID2EdgeID[6][4];

	// Array converting from edge ID to node ID
	int m_edgeID2NodeID[12][2];

	const static int m_numGauss = 2;

	const static int m_numIntegralPoints = m_numGauss * m_numGauss * m_numGauss;

	double m_integralPointXi[m_numIntegralPoints];

	double m_integralPointEta[m_numIntegralPoints];

	double m_integralPointZeta[m_numIntegralPoints];

	double m_weights[m_numIntegralPoints];

	// Array of reference coord xi values for each node
	double m_xiAtNode[8];

	// Array of reference coord eta values for each node
	double m_etaAtNode[8];

	// Array of reference coord zeta values for each node
	double m_zetaAtNode[8];

	// Check whether side element-faces are parallel to Z-X or Y-Z plane
	void checkWhetherSideFaceIsParallelToZXOrYZPlane() const;

	// Check whether the specified point is located in the specified element
	bool isLocatedInTheElement( const double x, const double y, const double z, const int iElem ) const; 
	
	// Calculate local coordinates
	void calcLocalCoordinates( const int iElem, const double x, const double y, const double z, double& xi, double& eta, double& zeta ) const;

	// Calculate horizontal local coordinates
	void calcHorizontalLocalCoordinates( const int iElem, const double x, const double y, double& xi, double& eta ) const;

	// Calculate determinant of jacobian matrix of the elements
	double calcDeterminantOfJacobianMatrix( const int iElem,  const double xi, const double eta, const double zeta ) const;

};

#endif
