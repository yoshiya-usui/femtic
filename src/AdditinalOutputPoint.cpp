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
#include <fstream>
#include <iostream>

#include "AnalysisControl.h"
#include "AdditinalOutputPoint.h"
#include "OutputFiles.h"
#include "CommonParameters.h"

// Constructer
AdditinalOutputPoint::AdditinalOutputPoint():
	m_elementIncludingStation(0)
{
	m_location.X = 0.0;
	m_location.Y = 0.0;
	m_location.Z = 0.0;

	m_localCoordinateValues.X = 0.0;
	m_localCoordinateValues.Y = 0.0;
	m_localCoordinateValues.Z = 0.0;

	m_volumeCoordinateValues.coord0 = 0.0;
	m_volumeCoordinateValues.coord1 = 0.0;
	m_volumeCoordinateValues.coord2 = 0.0;
	m_volumeCoordinateValues.coord3 = 0.0;

	for( int iPol = 0; iPol < 2; ++iPol ){
		m_ExCalculated[iPol] =NULL;
		m_EyCalculated[iPol] =NULL;
		m_EzCalculated[iPol] =NULL;
		m_HxCalculated[iPol] =NULL;
		m_HyCalculated[iPol] =NULL;
		m_HzCalculated[iPol] =NULL;
	}

}

// Destructer
AdditinalOutputPoint::~AdditinalOutputPoint(){

	for( int iPol = 0; iPol < 2; ++iPol ){

		if( m_ExCalculated[iPol] != NULL ){
			delete [] m_ExCalculated[iPol];
			m_ExCalculated[iPol] =NULL;
		}

		if( m_EyCalculated[iPol] != NULL ){
			delete [] m_EyCalculated[iPol];
			m_EyCalculated[iPol] =NULL;
		}

		if( m_EzCalculated[iPol] != NULL ){
			delete [] m_EzCalculated[iPol];
			m_EzCalculated[iPol] =NULL;
		}

		if( m_HxCalculated[iPol] != NULL ){
			delete [] m_HxCalculated[iPol];
			m_HxCalculated[iPol] =NULL;
		}

		if( m_HyCalculated[iPol] != NULL ){
			delete [] m_HyCalculated[iPol];
			m_HyCalculated[iPol] =NULL;
		}

		if( m_HzCalculated[iPol] != NULL ){
			delete [] m_HzCalculated[iPol];
			m_HzCalculated[iPol] =NULL;
		}

	}

}

// Read data from input file
void AdditinalOutputPoint::inputObservedData( std::ifstream& inFile ){

	double dbuf(0.0);
	inFile >> dbuf;
	m_location.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.Z = dbuf * CommonParameters::convKilometerToMeter;

#ifdef _DEBUG_WRITE
	std::cout << m_location.X << " "
		      << m_location.Y << " "
			  << m_location.Z << std::endl;
#endif

}

// Find element including station
void AdditinalOutputPoint::findElementIncludingStation(){

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();
		m_elementIncludingStation = ptrMeshDataTetraElement->findElementIncludingPoint( m_location.X, m_location.Y, m_location.Z, m_volumeCoordinateValues );
#ifdef _DEBUG_WRITE
		std::cout << "m_elementIncludingStation m_localCoordinateValues.coord0 m_localCoordinateValues.coord1 m_localCoordinateValues.coord2 m_localCoordinateValues.coord3 "
			<< m_elementIncludingStation << " "
			<< m_volumeCoordinateValues.coord0 << " " << m_volumeCoordinateValues.coord1 << " " << m_volumeCoordinateValues.coord2 << " " << m_volumeCoordinateValues.coord3 << std::endl;
#endif
	}else if( meshType == MeshData::HEXA ){// Hexa mesh
		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();
		double dummy(0.0);
		m_elementIncludingStation = ptrMeshDataBrickElement->findElementIncludingPoint( m_location.X, m_location.Y, m_location.Z,
			m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, false, false, dummy, dummy );
	}else if( meshType == MeshData::NONCONFORMING_HEXA ){
		const MeshDataNonConformingHexaElement* const ptrMeshDataHexaElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();
		m_elementIncludingStation = ptrMeshDataHexaElement->findElementIncludingPoint( m_location.X, m_location.Y, m_location.Z,
			m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Calulate EM field
void AdditinalOutputPoint::calculateEMField( const int ifreq, const Forward3D* const ptrForward3D ){

	if( m_elementIncludingStation < 0 ){
		return;
	}

	const int iPol = ptrForward3D->getPolarizationCurrent();

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		m_ExCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_EyCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_EzCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldZDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_HxCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_HyCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
		m_HzCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldZDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
	}else if( meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA ){
		m_ExCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_EyCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_EzCalculated[iPol][ifreq] = ptrForward3D->calcValueElectricFieldZDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_HxCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_HyCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_HzCalculated[iPol][ifreq] = ptrForward3D->calcValueMagneticFieldZDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Initialize EM field
void AdditinalOutputPoint::initializeEMfield( const int nfreq ){

	for( int iPol = 0; iPol < 2; ++iPol ){
		for( int ifreq = 0; ifreq < nfreq; ++ifreq ){
			m_ExCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
			m_EyCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
			m_EzCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
			m_HxCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
			m_HyCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
			m_HzCalculated[iPol][ifreq] = std::complex<double>(0.0,0.0);
		}
	}

}

// Allcate memory for EM field
void AdditinalOutputPoint::allocateMemoryForCalculatedValues( const int nfreq ){

	for( int iPol = 0; iPol < 2; ++iPol ){
		m_ExCalculated[iPol] = new std::complex<double>[nfreq];
		m_EyCalculated[iPol] = new std::complex<double>[nfreq];
		m_EzCalculated[iPol] = new std::complex<double>[nfreq];
		m_HxCalculated[iPol] = new std::complex<double>[nfreq];
		m_HyCalculated[iPol] = new std::complex<double>[nfreq];
		m_HzCalculated[iPol] = new std::complex<double>[nfreq];
	}

}
			
// Output calculated values of EM field
void AdditinalOutputPoint::outputCalculatedValues( const int nfreq, const double* const freq ) const{

	if( m_elementIncludingStation < 0 ){
		return;
	}

	for( int ifreq = 0; ifreq < nfreq; ++ifreq ){

		fprintf( OutputFiles::m_csvFile, " %15.6lf, %15.6lf, %15.6lf, %15.6e,",
			m_location.X*0.001, m_location.Y*0.001, m_location.Z*0.001, freq[ifreq] );
		for( int iPol = 0; iPol < 2; ++iPol ){
 			fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
				m_ExCalculated[iPol][ifreq].real(), m_ExCalculated[iPol][ifreq].imag(),
				m_EyCalculated[iPol][ifreq].real(), m_EyCalculated[iPol][ifreq].imag(),
				m_EzCalculated[iPol][ifreq].real(), m_EzCalculated[iPol][ifreq].imag(),
				m_HxCalculated[iPol][ifreq].real(), m_HxCalculated[iPol][ifreq].imag(),
				m_HyCalculated[iPol][ifreq].real(), m_HyCalculated[iPol][ifreq].imag(),
				m_HzCalculated[iPol][ifreq].real(), m_HzCalculated[iPol][ifreq].imag() );
		}
 		fprintf( OutputFiles::m_csvFile, "\n");

	}

}
