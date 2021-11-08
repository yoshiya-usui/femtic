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
#include <string>
#include <assert.h>

#include "ObservedDataStationPoint.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"
#include "OutputFiles.h"

// Constructer
ObservedDataStationPoint::ObservedDataStationPoint():
	ObservedDataStation(),
	m_elementIncludingStation(0),
	m_faceIncludingStation(0),
	m_useUpperElementForInterpolationOfEMField(false),
	m_rhsVectorIDOfHx(0),
	m_rhsVectorIDOfHy(0)
{
	m_location.X = 0.0;
	m_location.Y = 0.0;

	m_localCoordinateValues.X = 0.0;
	m_localCoordinateValues.Y = 0.0;
	m_localCoordinateValues.Z = 0.0;

	for( int i = 0; i < 2; ++i ){
		m_HxCalculated[i] = std::complex<double>(0.0,0.0);
		m_HyCalculated[i] = std::complex<double>(0.0,0.0);
	}

	m_volumeCoordinateValues.coord0 = 0.0;
	m_volumeCoordinateValues.coord1 = 0.0;
	m_volumeCoordinateValues.coord2 = 0.0;
	m_volumeCoordinateValues.coord3 = 0.0;

	m_areaCoordinateValues.coord0 = 0.0;
	m_areaCoordinateValues.coord1 = 0.0;
	m_areaCoordinateValues.coord2 = 0.0;

}

// Destructer
ObservedDataStationPoint::~ObservedDataStationPoint(){
}

// Find element including a point
void ObservedDataStationPoint::findElementIncludingStation(){

	const bool modLoc = ( AnalysisControl::getInstance() )->getIsObsLocMovedToCenter();

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh

		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();

		double locXMod( m_location.X );
		double locYMod( m_location.Y );
		m_elementIncludingStation = ptrMeshDataTetraElement->findElementIncludingPointOnSurface( m_location.X, m_location.Y, m_faceIncludingStation, m_areaCoordinateValues,
			useUpperElementForInterpolationOfEMField(),	modLoc, locXMod, locYMod );
		if(modLoc){
			m_location.X = locXMod;
			m_location.Y = locYMod;
		}

		ptrMeshDataTetraElement->calcVolumeCoordFromAreaCoord( m_faceIncludingStation, m_areaCoordinateValues, m_volumeCoordinateValues );

#ifdef _DEBUG_WRITE
		std::cout << " m_stationID m_elementIncludingStation m_faceIncludingStation m_localCoordinateValues.x m_localCoordinateValues.y m_localCoordinateValues.z " 
			<< m_stationID << " " << m_elementIncludingStation << " " << m_faceIncludingStation << " "
			<< m_volumeCoordinateValues.coord0 << " " << m_volumeCoordinateValues.coord1 << " " << m_volumeCoordinateValues.coord2 << " " << m_volumeCoordinateValues.coord3 << std::endl;
#endif
	}else if( meshType == MeshData::HEXA ){// Hexa mesh

		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();

		double locXMod( m_location.X );
		double locYMod( m_location.Y );
		m_elementIncludingStation = ptrMeshDataBrickElement->findElementIncludingPointOnSurface( m_location.X, m_location.Y,
			m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, useUpperElementForInterpolationOfEMField(), modLoc, locXMod, locYMod );
		if(modLoc){
			m_location.X = locXMod;
			m_location.Y = locYMod;
		}

#ifdef _DEBUG_WRITE
		std::cout << " m_stationID m_elementIncludingStation m_localCoordinateValues.x m_localCoordinateValues.y m_localCoordinateValues.z " 
			<< m_stationID << " " << m_elementIncludingStation << " "
			<< m_localCoordinateValues.X << " " << m_localCoordinateValues.Y << " " << m_localCoordinateValues.Z << std::endl;
#endif

	}else if( meshType == MeshData::NONCONFORMING_HEXA ){

		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();

		double locXMod( m_location.X );
		double locYMod( m_location.Y );
		m_elementIncludingStation = ptrMeshDataNonConformingHexaElement->findElementIncludingPointOnSurface( m_location.X, m_location.Y, m_faceIncludingStation,
			m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, useUpperElementForInterpolationOfEMField(), modLoc, locXMod, locYMod );
		if(modLoc){
			m_location.X = locXMod;
			m_location.Y = locYMod;
		}

	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Get caluculated value of Hx
std::complex<double> ObservedDataStationPoint::getHxCalculated( const int iPol ) const{
	return m_HxCalculated[iPol];
}

// Get caluculated value of Hy
std::complex<double> ObservedDataStationPoint::getHyCalculated( const int iPol ) const{
	return m_HyCalculated[iPol];
}

// Get right-hand side vector ID of Hx
int ObservedDataStationPoint::getRhsVectorIDOfHx() const{
	return m_rhsVectorIDOfHx;
}

// Get right-hand side vector ID of Hy
int ObservedDataStationPoint::getRhsVectorIDOfHy() const{
	return m_rhsVectorIDOfHy;
}

// Calulate horizontal magnetic field
void ObservedDataStationPoint::calculateHorizontalMagneticField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfHx, const int rhsVectorIDOfHy ){

	const int iPol = ptrForward3D->getPolarizationCurrent();

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		m_HxCalculated[iPol] = ptrForward3D->calcValueMagneticFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_HyCalculated[iPol] = ptrForward3D->calcValueMagneticFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
	}else if( meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA ){
		m_HxCalculated[iPol] = ptrForward3D->calcValueMagneticFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_HyCalculated[iPol] = ptrForward3D->calcValueMagneticFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

	// For inversion
	m_rhsVectorIDOfHx = rhsVectorIDOfHx;
	m_rhsVectorIDOfHy = rhsVectorIDOfHy;
}

// Initialize horizontal magnetic field
void ObservedDataStationPoint::initializeHorizontalMagneticField( const int iPol ){

	m_HxCalculated[iPol] = std::complex<double>(0.0,0.0);
	m_HyCalculated[iPol] = std::complex<double>(0.0,0.0);

}

// Calulate interpolator vector of horizontal magnetic field
void ObservedDataStationPoint::calcInterpolatorVectorOfHorizontalMagneticField( Forward3D* const ptrForward3D ){

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHx );
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHy );
	}else if( meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA ){
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHx );
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHy );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Get location of point
const CommonParameters::locationXY& ObservedDataStationPoint::getLocationOfPoint() const{

	return m_location;

}

// Get Z coordinate of the point
double ObservedDataStationPoint::getZCoordOfPoint() const{

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();
		return ptrMeshDataTetraElement->calcZCoordOfPointOnFace( m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues );
	}else if( meshType == MeshData::HEXA ){
		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();
		return ptrMeshDataBrickElement->calcGlobalCoordZ( m_elementIncludingStation, m_localCoordinateValues.Z );
	}else if( meshType == MeshData::NONCONFORMING_HEXA ){
		const MeshDataNonConformingHexaElement* const ptrMeshDataHexaElement =  ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();
		return ptrMeshDataHexaElement->calcZCoordOfPointOnFace( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Get flag specifing whether the EM field is interpolated from the values of the upper element
bool ObservedDataStationPoint::useUpperElementForInterpolationOfEMField() const{

	return m_useUpperElementForInterpolationOfEMField;

}

// Set flag specifing whether the EM field is interpolated from the values of the upper element
void ObservedDataStationPoint::setFlagUseUpperElementForInterpolationOfEMField( const bool useUpperElem ){

	m_useUpperElementForInterpolationOfEMField = useUpperElem;

}


