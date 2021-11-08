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
#ifndef DBLDEF_MESHDATA_BRICKELEMENT
#define DBLDEF_MESHDATA_BRICKELEMENT

#include <vector>
#include "MeshData.h"

// Class of FEM mesh for brick element
class MeshDataBrickElement : public MeshData {

public:

	// Constructer
	MeshDataBrickElement();

	// Destructer
	virtual ~MeshDataBrickElement();

	// Input mesh data from "mesh.dat"
	virtual void inputMeshData();

	// Find element including a point
	int findElementIncludingPoint( const double locX, const double locY, const double locZ,
		double& localCoordX, double& localCoordY, double& localCoordZ, const bool useUpperElem,
		const bool modLoc, double& locXMod, double& locYMod ) const;

	// Find elements including point on the surface of the earth
	int findElementIncludingPointOnSurface( const double locX, const double locY,
		double& localCoordX, double& localCoordY, double& localCoordZ, const bool useUpperElem,
		const bool modLoc, double& locXMod, double& locYMod ) const;

	// Find elements including dipole on the surface of the earth
	void findElementsIncludingDipoleOnSurface( const double locXStart, const double locYStart, const double locXEnd, const double locYEnd,
		std::vector<int>& elements, std::vector<double>& localCoordXStartPoint, std::vector<double>& localCoordYStartPoint,	std::vector<double>& localCoordXEndPoint, std::vector<double>& localCoordYEndPoint ) const;

	// Get length of the edges parallel to X coordinate
	double getEdgeLengthX( const int iElem ) const;

	// Get length of the edges parallel to Y coordinate
	double getEdgeLengthY( const int iElem ) const;

	// Get length of the edges parallel to Z coordinate
	double getEdgeLengthZ( const int iElem ) const;

	// Get global X coordinate from local coordinate
	double calcGlobalCoordX( const int iElem, double localCoordX ) const;

	// Get global Y coordinate from local coordinate
	double calcGlobalCoordY( const int iElem, double localCoordY ) const;

	// Get global Z coordinate from local coordinate
	double calcGlobalCoordZ( const int iElem, double localCoordZ ) const;

	// Get number of Elements parallel to X direction
	int getNumElemX() const;

	// Get number of Elements parallel to Y direction
	int getNumElemY() const;

	// Get number of Elements parallel to Z direction
	int getNumElemZ() const;

	// Get number of the air layer
	int getNumAirLayer() const;	

	// Calculate number of edges of X-Y plane
	int calcNumEdgesOnXYPlane() const;

	// Calculate number of edges of Y-Z plane
	int calcNumEdgesOnYZPlane() const;

	// Calculate number of edges of Z-X plane
	int calcNumEdgesOnZXPlane() const;

	// Get ID of the nodes of elements belonging to the boundary planes
	virtual int getNodesOfElementsBoundaryPlanes(  const int iPlane, const int iElem, const int iNode ) const;

	// Get mesh type
	virtual int getMeshType() const;

	// Decide whether specified elements share same edges
	virtual bool shareSameEdges( const int elemID1, const int elemID2 ) const;

	// Calculate volume of a specified element
	virtual double calcVolume( const int elemID ) const;

	// Calculate area of face
	virtual double calcAreaOfFace( const int iElem, const int iFace ) const;

	// Calculate area of face at bottom of mesh
	virtual double calcAreaOfFaceAtBottomOfMesh( const int iElem ) const;

	// Output mesh data to VTK file
	virtual void outputMeshDataToVTK() const;

	// Output mesh data to binary file
	virtual void outputMeshDataToBinary() const;

private:

	// Copy constructer
	MeshDataBrickElement(const MeshDataBrickElement& rhs);

	// Copy assignment operator
	MeshDataBrickElement& operator=(const MeshDataBrickElement& rhs);

	// Number of Elements parallel to X direction
	int m_numElemX;

	// Number of Elements parallel to Y direction
	int m_numElemY;

	// Number of Elements parallel to Z direction
	int m_numElemZ;

	// Number of the air layer
	int m_numAirLayer;

	// Array of edge lenghts
	double* m_edgeLength;

	// Array of nodes of elements belonging to the boundary planes
	//   m_nodesOfElementsBoundaryPlanes[0] : Y-Z Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[1] : Y-Z Plane ( Plus Side  )
	//   m_nodesOfElementsBoundaryPlanes[2] : Z-X Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[3] : Z-X Plane ( Plus Side  )
	//   m_nodesOfElementsBoundaryPlanes[4] : X-Y Plane ( Minus Side )
	//   m_nodesOfElementsBoundaryPlanes[5] : X-Y Plane ( Plus Side  )
	int* m_nodesOfElementsBoundaryPlanes[6];

	// Get local coordinate values from coordinate values
	virtual void getLocalCoordinateValues( const int iElem, const double coordX, const double coordY, const double coordZ,
		double& localCoordX, double& localCoordY, double& localCoordZ ) const;

};

#endif
