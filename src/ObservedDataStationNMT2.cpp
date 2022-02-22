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
#include <string>
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "ObservedDataStationNMT2.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationNMT2::ObservedDataStationNMT2():
	ObservedDataStation(),
	m_ZxxObserved(NULL),
	m_ZxyObserved(NULL),
	m_ZyxObserved(NULL),
	m_ZyyObserved(NULL),
	m_ZxxSD(NULL),
	m_ZxySD(NULL),
	m_ZyxSD(NULL),
	m_ZyySD(NULL),
	m_ZxxCalculated(NULL),
	m_ZxyCalculated(NULL),
	m_ZyxCalculated(NULL),
	m_ZyyCalculated(NULL),
	m_ZxxResidual(NULL),
	m_ZxyResidual(NULL),
	m_ZyxResidual(NULL),
	m_ZyyResidual(NULL),
	m_dataIDOfZxx(NULL),
	m_dataIDOfZxy(NULL),
	m_dataIDOfZyx(NULL),
	m_dataIDOfZyy(NULL)
{
	for( int i = 0; i < 2 ; ++i ){
		m_location[i].startPoint.X = 0.0;
		m_location[i].startPoint.Y = 0.0;
		m_location[i].endPoint.X = 0.0;
		m_location[i].endPoint.Y = 0.0;
		m_voltageCalculated[i][0] = std::complex<double>(0.0,0.0);
		m_voltageCalculated[i][1] = std::complex<double>(0.0,0.0);
		m_numElementsIncludingDipole[i] = 0;
		m_elementsIncludingDipole[i] =NULL;
		m_facesIncludingDipole[i] =NULL;
		m_localCoordinateValuesStartPoint[i] = NULL;
		m_localCoordinateValuesEndPoint[i] = NULL;
		m_areaCoordinateValuesStartPoint[i] = NULL;
		m_areaCoordinateValuesEndPoint[i] = NULL;
		m_rhsVectorIDOfVoltageDifference[i] = 0;
	}

}

// Destructer
ObservedDataStationNMT2::~ObservedDataStationNMT2(){

	if( m_ZxxObserved != NULL){
		delete[] m_ZxxObserved;
		m_ZxxObserved = NULL;
	}

	if( m_ZxyObserved != NULL){
		delete[] m_ZxyObserved;
		m_ZxyObserved = NULL;
	}

	if( m_ZyxObserved != NULL){
		delete[] m_ZyxObserved;
		m_ZyxObserved = NULL;
	}

	if( m_ZyyObserved != NULL){
		delete[] m_ZyyObserved;
		m_ZyyObserved = NULL;
	}

	if( m_ZxxSD != NULL){
		delete[] m_ZxxSD;
		m_ZxxSD = NULL;
	}

	if( m_ZxySD != NULL){
		delete[] m_ZxySD;
		m_ZxySD = NULL;
	}

	if( m_ZyxSD != NULL){
		delete[] m_ZyxSD;
		m_ZyxSD = NULL;
	}

	if( m_ZyySD != NULL){
		delete[] m_ZyySD;
		m_ZyySD = NULL;
	}

	if( m_ZxxCalculated != NULL){
		delete[] m_ZxxCalculated;
		m_ZxxCalculated = NULL;
	}

	if( m_ZxyCalculated != NULL){
		delete[] m_ZxyCalculated;
		m_ZxyCalculated = NULL;
	}

	if( m_ZyxCalculated != NULL){
		delete[] m_ZyxCalculated;
		m_ZyxCalculated = NULL;
	}

	if( m_ZyyCalculated != NULL){
		delete[] m_ZyyCalculated;
		m_ZyyCalculated = NULL;
	}

	if( m_ZxxResidual != NULL){
		delete[] m_ZxxResidual;
		m_ZxxResidual = NULL;
	}

	if( m_ZxyResidual != NULL){
		delete[] m_ZxyResidual;
		m_ZxyResidual = NULL;
	}

	if( m_ZyxResidual != NULL){
		delete[] m_ZyxResidual;
		m_ZyxResidual = NULL;
	}

	if( m_ZyyResidual != NULL){
		delete[] m_ZyyResidual;
		m_ZyyResidual = NULL;
	}

	for( int i = 0; i < 2; ++i ){

		if( m_elementsIncludingDipole[i] != NULL ){
			delete[] m_elementsIncludingDipole[i];
			m_elementsIncludingDipole[i] = NULL;
		}

		if( m_facesIncludingDipole[i] != NULL ){
			delete[] m_facesIncludingDipole[i];
			m_facesIncludingDipole[i] = NULL;
		}

		if( m_localCoordinateValuesStartPoint[i] != NULL ){
			delete[] m_localCoordinateValuesStartPoint[i];
			m_localCoordinateValuesStartPoint[i] = NULL;
		}

		if( m_localCoordinateValuesEndPoint[i] != NULL ){
			delete[] m_localCoordinateValuesEndPoint[i];
			m_localCoordinateValuesEndPoint[i] = NULL;
		}

		if( m_areaCoordinateValuesStartPoint[i] != NULL ){
			delete[] m_areaCoordinateValuesStartPoint[i];
			m_areaCoordinateValuesStartPoint[i] = NULL;
		}

		if( m_areaCoordinateValuesEndPoint[i] != NULL ){
			delete[] m_areaCoordinateValuesEndPoint[i];
			m_areaCoordinateValuesEndPoint[i] = NULL;
		}

	}

	if( m_dataIDOfZxx != NULL){
		delete[] m_dataIDOfZxx;
		m_dataIDOfZxx = NULL;
	}

	if( m_dataIDOfZxy != NULL){
		delete[] m_dataIDOfZxy;
		m_dataIDOfZxy = NULL;
	}

	if( m_dataIDOfZyx != NULL){
		delete[] m_dataIDOfZyx;
		m_dataIDOfZyx = NULL;
	}

	if( m_dataIDOfZyy != NULL){
		delete[] m_dataIDOfZyy;
		m_dataIDOfZyy = NULL;
	}

}

// Read data from input file
void ObservedDataStationNMT2::inputObservedData( std::ifstream& inFile ){

	inFile >> m_stationID;

	inFile >> m_IDOfMagneticFieldStation;

	OutputFiles::m_logFile << "# " << std::setw(15) << std::left << m_stationID << std::setw(15) << std::left << m_IDOfMagneticFieldStation << std::endl;

	double dbuf(0.0);
	inFile >> dbuf;
	m_location[0].startPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[0].startPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[0].endPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[0].endPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[1].startPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[1].startPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[1].endPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location[1].endPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> m_numOfFrequency;
	const int nFreq = m_numOfFrequency;
	if( nFreq > 0 ){
		m_freq = new double[nFreq];
		m_ZxxObserved = new std::complex<double>[nFreq];
		m_ZxyObserved = new std::complex<double>[nFreq];
		m_ZyxObserved = new std::complex<double>[nFreq];
		m_ZyyObserved = new std::complex<double>[nFreq];
		m_ZxxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZxySD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyySD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZxxCalculated = new std::complex<double>[nFreq];
		m_ZxyCalculated = new std::complex<double>[nFreq];
		m_ZyxCalculated = new std::complex<double>[nFreq];
		m_ZyyCalculated = new std::complex<double>[nFreq];
		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1 >> dbuf2;
			m_ZxxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_ZxyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_ZyxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_ZyyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> m_ZxxSD[i].realPart;
			inFile >> m_ZxxSD[i].imagPart;
			inFile >> m_ZxySD[i].realPart;
			inFile >> m_ZxySD[i].imagPart;
			inFile >> m_ZyxSD[i].realPart;
			inFile >> m_ZyxSD[i].imagPart;
			inFile >> m_ZyySD[i].realPart;
			inFile >> m_ZyySD[i].imagPart;
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << " NMT2 " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location[0].startPoint.X << " "
				<< m_location[0].startPoint.Y << " "
				<< m_location[0].endPoint.X << " "
				<< m_location[0].endPoint.Y << " "
				<< m_location[1].startPoint.X << " "
				<< m_location[1].startPoint.Y << " "
				<< m_location[1].endPoint.X << " "
				<< m_location[1].endPoint.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for( int i = 0; i < m_numOfFrequency; ++i ){
		std::cout << m_freq[i] << " "
					<< m_ZxxObserved[i] << " "
					<< m_ZxyObserved[i] << " "
					<< m_ZyxObserved[i] << " "
					<< m_ZyyObserved[i] << " "
					<< m_ZxxSD[i].realPart << " "
					<< m_ZxxSD[i].imagPart << " "
					<< m_ZxySD[i].realPart << " "
					<< m_ZxySD[i].imagPart << " "
					<< m_ZyxSD[i].realPart << " "
					<< m_ZyxSD[i].imagPart << " "
					<< m_ZyySD[i].realPart << " "
					<< m_ZyySD[i].imagPart << std::endl;
	}
#endif

}

// Find elements including dipole
void ObservedDataStationNMT2::findElementsIncludingDipoles(){

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh

		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();

		for( int idipole = 0; idipole < 2; ++idipole ){

			std::vector<CommonParameters::AreaCoords> localCoordStartPoint;
			std::vector<CommonParameters::AreaCoords> localCoordEndPoint;

			std::vector<int> elementsIncludingDipole;
			std::vector<int> facesIncludingDipole;

			ptrMeshDataTetraElement->findElementsIncludingDipoleOnSurface( m_location[idipole].startPoint.X, m_location[idipole].startPoint.Y,
				m_location[idipole].endPoint.X, m_location[idipole].endPoint.Y,
				elementsIncludingDipole, facesIncludingDipole, localCoordStartPoint, localCoordEndPoint );

			m_numElementsIncludingDipole[idipole] = static_cast<int>( elementsIncludingDipole.size() );
			if( m_numElementsIncludingDipole[idipole] > 0 ){
				m_elementsIncludingDipole[idipole] = new int[ m_numElementsIncludingDipole[idipole] ];
				m_facesIncludingDipole[idipole] = new int[ m_numElementsIncludingDipole[idipole] ];
				m_areaCoordinateValuesStartPoint[idipole] = new CommonParameters::AreaCoords[ m_numElementsIncludingDipole[idipole] ];
				m_areaCoordinateValuesEndPoint[idipole] = new CommonParameters::AreaCoords[ m_numElementsIncludingDipole[idipole] ];
			}

			for( int i = 0; i < m_numElementsIncludingDipole[idipole]; ++i ){
				m_elementsIncludingDipole[idipole][i] = elementsIncludingDipole[i];
				m_facesIncludingDipole[idipole][i] = facesIncludingDipole[i];
				m_areaCoordinateValuesStartPoint[idipole][i] = localCoordStartPoint[i];
				m_areaCoordinateValuesEndPoint[idipole][i] = localCoordEndPoint[i];
			}

		}

#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		for( int i = 0; i < m_numElementsIncludingDipole[0]; ++i ){
			std::cout << m_elementsIncludingDipole[0][i] << " " << m_facesIncludingDipole[0][i] << " " << m_areaCoordinateValuesStartPoint[0][i].coord0  << " " << m_areaCoordinateValuesStartPoint[0][i].coord1 << " " << m_areaCoordinateValuesStartPoint[0][i].coord2 << " " << std::endl;
			std::cout << m_elementsIncludingDipole[0][i] << " " << m_facesIncludingDipole[0][i] << " " << m_areaCoordinateValuesEndPoint[0][i].coord0  << " " << m_areaCoordinateValuesEndPoint[0][i].coord1 << " " << m_areaCoordinateValuesStartPoint[0][i].coord2 << " " << std::endl;
		}
		for( int i = 0; i < m_numElementsIncludingDipole[1]; ++i ){
			std::cout << m_elementsIncludingDipole[1][i] << " " << m_facesIncludingDipole[1][i] << " " << m_areaCoordinateValuesStartPoint[1][i].coord0  << " " << m_areaCoordinateValuesStartPoint[1][i].coord1 << " " << m_areaCoordinateValuesStartPoint[1][i].coord2 << std::endl;
			std::cout << m_elementsIncludingDipole[1][i] << " " << m_facesIncludingDipole[1][i] << " " << m_areaCoordinateValuesEndPoint[1][i].coord0  << " " << m_areaCoordinateValuesEndPoint[1][i].coord1 << " " << m_areaCoordinateValuesStartPoint[1][i].coord2 << std::endl;
		}
#endif
	}else if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::NONCONFORMING_HEXA ){// Non-conforming deformed hexahedral mesh
		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();
		for( int idipole = 0; idipole < 2; ++idipole ){
			std::vector<double> localCoordXStartPoint;
			std::vector<double> localCoordYStartPoint;
			std::vector<double> localCoordXEndPoint;
			std::vector<double> localCoordYEndPoint;
			std::vector<int> elementsIncludingDipole;
			ptrMeshDataNonConformingHexaElement->findElementsIncludingDipoleOnSurface(
				m_location[idipole].startPoint.X, m_location[idipole].startPoint.Y, m_location[idipole].endPoint.X, m_location[idipole].endPoint.Y,
				elementsIncludingDipole, localCoordXStartPoint, localCoordYStartPoint, localCoordXEndPoint, localCoordYEndPoint );
			m_numElementsIncludingDipole[idipole] = static_cast<int>( elementsIncludingDipole.size() );
			if( m_numElementsIncludingDipole[idipole] > 0 ){
				m_elementsIncludingDipole[idipole]        = new int[ m_numElementsIncludingDipole[idipole] ];
				m_localCoordinateValuesStartPoint[idipole] = new CommonParameters::locationXY[ m_numElementsIncludingDipole[idipole] ];
				m_localCoordinateValuesEndPoint[idipole]   = new CommonParameters::locationXY[ m_numElementsIncludingDipole[idipole] ];
			}
			for( int ielem = 0; ielem < m_numElementsIncludingDipole[idipole]; ++ielem ){
				m_elementsIncludingDipole[idipole][ielem] = elementsIncludingDipole[ielem];
				m_localCoordinateValuesStartPoint[idipole][ielem].X = localCoordXStartPoint[ielem];
				m_localCoordinateValuesStartPoint[idipole][ielem].Y = localCoordYStartPoint[ielem];
				m_localCoordinateValuesEndPoint[idipole][ielem].X   = localCoordXEndPoint[ielem];
				m_localCoordinateValuesEndPoint[idipole][ielem].Y   = localCoordYEndPoint[ielem];
			}
		}
#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		std::cout << "m_elementsIncludingDipole[0] m_localCoordinateValuesStartPoint[0].x m_localCoordinateValuesStartPoint[0].y m_localCoordinateValuesEndPoint[0].x m_localCoordinateValuesEndPoint[0].y" << std::endl;
		for( int ielem = 0; ielem < m_numElementsIncludingDipole[0]; ++ielem ){
			std::cout << m_elementsIncludingDipole[0][ielem] << " " << m_localCoordinateValuesStartPoint[0][ielem].X << " " << m_localCoordinateValuesStartPoint[0][ielem].Y << " " << m_localCoordinateValuesEndPoint[0][ielem].X << " " << m_localCoordinateValuesEndPoint[0][ielem].Y << std::endl;
		}
		std::cout << "m_elementsIncludingDipole[1] m_localCoordinateValuesStartPoint[1].x m_localCoordinateValuesStartPoint[1].y m_localCoordinateValuesEndPoint[1].x m_localCoordinateValuesEndPoint[1].y" << std::endl;
		for( int ielem = 0; ielem < m_numElementsIncludingDipole[1]; ++ielem ){
			std::cout << m_elementsIncludingDipole[1][ielem] << " " << m_localCoordinateValuesStartPoint[1][ielem].X << " " << m_localCoordinateValuesStartPoint[1][ielem].Y << " " << m_localCoordinateValuesEndPoint[1][ielem].X << " " << m_localCoordinateValuesEndPoint[1][ielem].Y << std::endl;
		}
#endif
	}else{// Hexa mesh

		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();

		for( int idipole = 0; idipole < 2; ++idipole ){
			std::vector<double> localCoordXStartPoint;
			std::vector<double> localCoordYStartPoint;
			std::vector<double> localCoordXEndPoint;
			std::vector<double> localCoordYEndPoint;
			std::vector<int> elementsIncludingDipole;
			ptrMeshDataBrickElement->findElementsIncludingDipoleOnSurface(
				m_location[idipole].startPoint.X, m_location[idipole].startPoint.Y, m_location[idipole].endPoint.X, m_location[idipole].endPoint.Y,
				elementsIncludingDipole, localCoordXStartPoint, localCoordYStartPoint, localCoordXEndPoint, localCoordYEndPoint );

			m_numElementsIncludingDipole[idipole] = static_cast<int>( elementsIncludingDipole.size() );
			if( m_numElementsIncludingDipole[idipole] > 0 ){
				m_elementsIncludingDipole[idipole]         = new int[ m_numElementsIncludingDipole[idipole] ];
				m_localCoordinateValuesStartPoint[idipole] = new CommonParameters::locationXY[ m_numElementsIncludingDipole[idipole] ];
				m_localCoordinateValuesEndPoint[idipole]   = new CommonParameters::locationXY[ m_numElementsIncludingDipole[idipole] ];
			}

			for( int ielem = 0; ielem < m_numElementsIncludingDipole[idipole]; ++ielem ){
				m_elementsIncludingDipole[idipole][ielem] = elementsIncludingDipole[ielem];
				m_localCoordinateValuesStartPoint[idipole][ielem].X = localCoordXStartPoint[ielem];
				m_localCoordinateValuesStartPoint[idipole][ielem].Y = localCoordYStartPoint[ielem];
				m_localCoordinateValuesEndPoint[idipole][ielem].X   = localCoordXEndPoint[ielem];
				m_localCoordinateValuesEndPoint[idipole][ielem].Y   = localCoordYEndPoint[ielem];
			}
		}
	
#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		std::cout << "m_elementsIncludingDipole[0] m_localCoordinateValuesStartPoint[0].x m_localCoordinateValuesStartPoint[0].y m_localCoordinateValuesEndPoint[0].x m_localCoordinateValuesEndPoint[0].y" << std::endl;
		for( int ielem = 0; ielem < m_numElementsIncludingDipole[0]; ++ielem ){
			std::cout << m_elementsIncludingDipole[0][ielem] << " " << m_localCoordinateValuesStartPoint[0][ielem].X << " " << m_localCoordinateValuesStartPoint[0][ielem].Y << " " << m_localCoordinateValuesEndPoint[0][ielem].X << " " << m_localCoordinateValuesEndPoint[0][ielem].Y << std::endl;
		}
		std::cout << "m_elementsIncludingDipole[1] m_localCoordinateValuesStartPoint[1].x m_localCoordinateValuesStartPoint[1].y m_localCoordinateValuesEndPoint[1].x m_localCoordinateValuesEndPoint[1].y" << std::endl;
		for( int ielem = 0; ielem < m_numElementsIncludingDipole[1]; ++ielem ){
			std::cout << m_elementsIncludingDipole[1][ielem] << " " << m_localCoordinateValuesStartPoint[1][ielem].X << " " << m_localCoordinateValuesStartPoint[1][ielem].Y << " " << m_localCoordinateValuesEndPoint[1][ielem].X << " " << m_localCoordinateValuesEndPoint[1][ielem].Y << std::endl;
		}
#endif

	}

}

// Calulate difference of voltages
void ObservedDataStationNMT2::calculateVoltageDifferences( const Forward3D* const ptrForward3D, const int rhsVectorIDOfVoltageDifference1st, const int rhsVectorIDOfVoltageDifference2nd ){

	const int iPol = ptrForward3D->getPolarizationCurrent();

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh
		for( int iDipole = 0; iDipole < 2; ++iDipole ){
			m_voltageCalculated[iDipole][iPol] = ptrForward3D->calcVoltageDifference( m_numElementsIncludingDipole[iDipole], m_elementsIncludingDipole[iDipole], m_facesIncludingDipole[iDipole],
				m_areaCoordinateValuesStartPoint[iDipole], m_areaCoordinateValuesEndPoint[iDipole] );
		}
	}else{// Hexa mesh
		for( int iDipole = 0; iDipole < 2; ++iDipole ){
			m_voltageCalculated[iDipole][iPol] = ptrForward3D->calcVoltageDifference( m_numElementsIncludingDipole[iDipole], m_elementsIncludingDipole[iDipole],
				m_localCoordinateValuesStartPoint[iDipole], m_localCoordinateValuesEndPoint[iDipole] );
		}
	}


#ifdef _DEBUG_WRITE
	std::cout << "iPol V1 V2 : " << iPol << " " << m_voltageCalculated[0][iPol] << " " << m_voltageCalculated[1][iPol] << std::endl;
#endif

	// For inversion
	m_rhsVectorIDOfVoltageDifference[0] = rhsVectorIDOfVoltageDifference1st;
	m_rhsVectorIDOfVoltageDifference[1] = rhsVectorIDOfVoltageDifference2nd;

}

// Calulate Impedance tensor
//void ObservedDataStationNMT2::calculateImpedanceTensor( const int freqIDAmongThisPE, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){
void ObservedDataStationNMT2::calculateImpedanceTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };
	const std::complex<double> det = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];

	const double lengthOfDipole1 = hypot( m_location[0].endPoint.X - m_location[0].startPoint.X, m_location[0].endPoint.Y - m_location[0].startPoint.Y );
	const double lengthOfDipole2 = hypot( m_location[1].endPoint.X - m_location[1].startPoint.X, m_location[1].endPoint.Y - m_location[1].startPoint.Y );

	const double angleOfDipole1 = atan2( m_location[0].endPoint.Y - m_location[0].startPoint.Y , m_location[0].endPoint.X - m_location[0].startPoint.X );
	const double angleOfDipole2 = atan2( m_location[1].endPoint.Y - m_location[1].startPoint.Y , m_location[1].endPoint.X - m_location[1].startPoint.X );


	std::complex<double> ElectricFieldDipole1[2] = { m_voltageCalculated[0][0] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 ) , m_voltageCalculated[0][1] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 ) };
	std::complex<double> ElectricFieldDipole2[2] = { m_voltageCalculated[1][0] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 ) , m_voltageCalculated[1][1] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 ) };

	std::complex<double> ExCalculated[2];
	std::complex<double> EyCalculated[2];
	for( int iPol = 0; iPol < 2 ; ++iPol ){
		ExCalculated[iPol] =   ElectricFieldDipole1[iPol] * static_cast< std::complex<double> >( sin( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) - ElectricFieldDipole2[iPol] * static_cast< std::complex<double> >( sin( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) );
		EyCalculated[iPol] = - ElectricFieldDipole1[iPol] * static_cast< std::complex<double> >( cos( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) + ElectricFieldDipole2[iPol] * static_cast< std::complex<double> >( cos( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) );
	}

	m_ZxxCalculated[freqIDThisPEInSta] = ( ExCalculated[0]*HyCalculated[1] - ExCalculated[1]*HyCalculated[0] ) / det;
	m_ZxyCalculated[freqIDThisPEInSta] = ( ExCalculated[1]*HxCalculated[0] - ExCalculated[0]*HxCalculated[1] ) / det;
	m_ZyxCalculated[freqIDThisPEInSta] = ( EyCalculated[0]*HyCalculated[1] - EyCalculated[1]*HyCalculated[0] ) / det;
	m_ZyyCalculated[freqIDThisPEInSta] = ( EyCalculated[1]*HxCalculated[0] - EyCalculated[0]*HxCalculated[1] ) / det;

	//m_ZxxResidual[freqIDThisPEInSta].realPart = std::fabs( ( m_ZxxCalculated[freqIDThisPEInSta].real() - m_ZxxObserved[freqIDGlobalInSta].real() ) / m_ZxxObserved[freqIDGlobalInSta].real() );
	//m_ZxxResidual[freqIDThisPEInSta].imagPart = std::fabs( ( m_ZxxCalculated[freqIDThisPEInSta].imag() - m_ZxxObserved[freqIDGlobalInSta].imag() ) / m_ZxxObserved[freqIDGlobalInSta].imag() );
	//m_ZxyResidual[freqIDThisPEInSta].realPart = std::fabs( ( m_ZxyCalculated[freqIDThisPEInSta].real() - m_ZxyObserved[freqIDGlobalInSta].real() ) / m_ZxyObserved[freqIDGlobalInSta].real() );
	//m_ZxyResidual[freqIDThisPEInSta].imagPart = std::fabs( ( m_ZxyCalculated[freqIDThisPEInSta].imag() - m_ZxyObserved[freqIDGlobalInSta].imag() ) / m_ZxyObserved[freqIDGlobalInSta].imag() );
	//m_ZyxResidual[freqIDThisPEInSta].realPart = std::fabs( ( m_ZyxCalculated[freqIDThisPEInSta].real() - m_ZyxObserved[freqIDGlobalInSta].real() ) / m_ZyxObserved[freqIDGlobalInSta].real() );
	//m_ZyxResidual[freqIDThisPEInSta].imagPart = std::fabs( ( m_ZyxCalculated[freqIDThisPEInSta].imag() - m_ZyxObserved[freqIDGlobalInSta].imag() ) / m_ZyxObserved[freqIDGlobalInSta].imag() );
	//m_ZyyResidual[freqIDThisPEInSta].realPart = std::fabs( ( m_ZyyCalculated[freqIDThisPEInSta].real() - m_ZyyObserved[freqIDGlobalInSta].real() ) / m_ZyyObserved[freqIDGlobalInSta].real() );
	//m_ZyyResidual[freqIDThisPEInSta].imagPart = std::fabs( ( m_ZyyCalculated[freqIDThisPEInSta].imag() - m_ZyyObserved[freqIDGlobalInSta].imag() ) / m_ZyyObserved[freqIDGlobalInSta].imag() );
	m_ZxxResidual[freqIDThisPEInSta].realPart = ( m_ZxxObserved[freqIDGlobalInSta].real() - m_ZxxCalculated[freqIDThisPEInSta].real() ) / m_ZxxSD[freqIDGlobalInSta].realPart;
	m_ZxxResidual[freqIDThisPEInSta].imagPart = ( m_ZxxObserved[freqIDGlobalInSta].imag() - m_ZxxCalculated[freqIDThisPEInSta].imag() ) / m_ZxxSD[freqIDGlobalInSta].imagPart;
	m_ZxyResidual[freqIDThisPEInSta].realPart = ( m_ZxyObserved[freqIDGlobalInSta].real() - m_ZxyCalculated[freqIDThisPEInSta].real() ) / m_ZxySD[freqIDGlobalInSta].realPart;
	m_ZxyResidual[freqIDThisPEInSta].imagPart = ( m_ZxyObserved[freqIDGlobalInSta].imag() - m_ZxyCalculated[freqIDThisPEInSta].imag() ) / m_ZxySD[freqIDGlobalInSta].imagPart;
	m_ZyxResidual[freqIDThisPEInSta].realPart = ( m_ZyxObserved[freqIDGlobalInSta].real() - m_ZyxCalculated[freqIDThisPEInSta].real() ) / m_ZyxSD[freqIDGlobalInSta].realPart;
	m_ZyxResidual[freqIDThisPEInSta].imagPart = ( m_ZyxObserved[freqIDGlobalInSta].imag() - m_ZyxCalculated[freqIDThisPEInSta].imag() ) / m_ZyxSD[freqIDGlobalInSta].imagPart;
	m_ZyyResidual[freqIDThisPEInSta].realPart = ( m_ZyyObserved[freqIDGlobalInSta].real() - m_ZyyCalculated[freqIDThisPEInSta].real() ) / m_ZyySD[freqIDGlobalInSta].realPart;
	m_ZyyResidual[freqIDThisPEInSta].imagPart = ( m_ZyyObserved[freqIDGlobalInSta].imag() - m_ZyyCalculated[freqIDThisPEInSta].imag() ) / m_ZyySD[freqIDGlobalInSta].imagPart;

#ifdef _DEBUG_WRITE
	std::cout << "freqIDThisPEInSta Zxx(NMT) Zxy(NMT) Zyx(NMT) Zyy(NMT) : " << freqIDThisPEInSta << " " << m_ZxxCalculated[freqIDThisPEInSta] << " " << m_ZxyCalculated[freqIDThisPEInSta] << " " << m_ZyxCalculated[freqIDThisPEInSta] << " " << m_ZyyCalculated[freqIDThisPEInSta] << std::endl;
#endif

	// For inversion
	//ObservedData* const ptrObservedData = ObservedData::getInstance();
	//this->m_dataIDOfZxx[freqIDThisPEInSta].realPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZxx[freqIDThisPEInSta].imagPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZxy[freqIDThisPEInSta].realPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZxy[freqIDThisPEInSta].imagPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZyx[freqIDThisPEInSta].realPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZyx[freqIDThisPEInSta].imagPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//this->m_dataIDOfZyy[freqIDThisPEInSta].realPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );	
	//this->m_dataIDOfZyy[freqIDThisPEInSta].imagPart = ptrObservedData->incrementNumObservedDataThisPE( ifreq );

	this->m_dataIDOfZxx[freqIDThisPEInSta].realPart = icount++;
	this->m_dataIDOfZxx[freqIDThisPEInSta].imagPart = icount++;
	this->m_dataIDOfZxy[freqIDThisPEInSta].realPart = icount++;
	this->m_dataIDOfZxy[freqIDThisPEInSta].imagPart = icount++;
	this->m_dataIDOfZyx[freqIDThisPEInSta].realPart = icount++;
	this->m_dataIDOfZyx[freqIDThisPEInSta].imagPart = icount++;
	this->m_dataIDOfZyy[freqIDThisPEInSta].realPart = icount++;	
	this->m_dataIDOfZyy[freqIDThisPEInSta].imagPart = icount++;

}

// Initialize difference of voltages
void ObservedDataStationNMT2::initializeVoltageDifferences( const int iPol ){

	m_voltageCalculated[0][iPol] = std::complex<double>(0.0,0.0);
	m_voltageCalculated[1][iPol] = std::complex<double>(0.0,0.0);

}

// Initialize Impedance tensor and errors
void ObservedDataStationNMT2::initializeImpedanceTensorsAndErrors(){

	//for( int i = 0; i < m_numOfFrequency; ++i ){
	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_ZxxResidual[i].realPart = 0.0;
		m_ZxxResidual[i].imagPart = 0.0;
		m_ZxyResidual[i].realPart = 0.0;
		m_ZxyResidual[i].imagPart = 0.0;
		m_ZyxResidual[i].realPart = 0.0;
		m_ZyxResidual[i].imagPart = 0.0;
		m_ZyyResidual[i].realPart = 0.0;
		m_ZyyResidual[i].imagPart = 0.0;

		m_ZxxCalculated[i] = std::complex<double>(0.0,0.0);
		m_ZxyCalculated[i] = std::complex<double>(0.0,0.0);
		m_ZyxCalculated[i] = std::complex<double>(0.0,0.0);
		m_ZyyCalculated[i] = std::complex<double>(0.0,0.0);
	}

}


// Allocate memory for the calculated Impedance and errors
void ObservedDataStationNMT2::allocateMemoryForCalculatedValues(){

	if( m_ZxxCalculated != NULL){
		delete[] m_ZxxCalculated;
		m_ZxxCalculated = NULL;
	}

	if( m_ZxyCalculated != NULL){
		delete[] m_ZxyCalculated;
		m_ZxyCalculated = NULL;
	}

	if( m_ZyxCalculated != NULL){
		delete[] m_ZyxCalculated;
		m_ZyxCalculated = NULL;
	}

	if( m_ZyyCalculated != NULL){
		delete[] m_ZyyCalculated;
		m_ZyyCalculated = NULL;
	}

	if( m_ZxxResidual != NULL){
		delete[] m_ZxxResidual;
		m_ZxxResidual = NULL;
	}

	if( m_ZxyResidual != NULL){
		delete[] m_ZxyResidual;
		m_ZxyResidual = NULL;
	}

	if( m_ZyxResidual != NULL){
		delete[] m_ZyxResidual;
		m_ZyxResidual = NULL;
	}

	if( m_ZyyResidual != NULL){
		delete[] m_ZyyResidual;
		m_ZyyResidual = NULL;
	}

	if( m_dataIDOfZxx != NULL){
		delete[] m_dataIDOfZxx;
		m_dataIDOfZxx = NULL;
	}

	if( m_dataIDOfZxy != NULL){
		delete[] m_dataIDOfZxy;
		m_dataIDOfZxy = NULL;
	}

	if( m_dataIDOfZyx != NULL){
		delete[] m_dataIDOfZyx;
		m_dataIDOfZyx = NULL;
	}

	if( m_dataIDOfZyy != NULL){
		delete[] m_dataIDOfZyy;
		m_dataIDOfZyy = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){
		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_ZxxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZxyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZyxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZyyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZxxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZxyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZyxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_ZyyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfZxx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfZxy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfZyx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfZyy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_ZxxCalculated[i] = std::complex<double>(0.0,0.0);
			m_ZxyCalculated[i] = std::complex<double>(0.0,0.0);
			m_ZyxCalculated[i] = std::complex<double>(0.0,0.0);
			m_ZyyCalculated[i] = std::complex<double>(0.0,0.0);
			m_ZxxResidual[i].realPart = 0.0;
			m_ZxxResidual[i].imagPart = 0.0;
			m_ZxyResidual[i].realPart = 0.0;
			m_ZxyResidual[i].imagPart = 0.0;
			m_ZyxResidual[i].realPart = 0.0;
			m_ZyxResidual[i].imagPart = 0.0;
			m_ZyyResidual[i].realPart = 0.0;
			m_ZyyResidual[i].imagPart = 0.0;
			m_dataIDOfZxx[i].realPart = -1;
			m_dataIDOfZxx[i].imagPart = -1;
			m_dataIDOfZxy[i].realPart = -1;
			m_dataIDOfZxy[i].imagPart = -1;
			m_dataIDOfZyx[i].realPart = -1;
			m_dataIDOfZyx[i].imagPart = -1;
			m_dataIDOfZyy[i].realPart = -1;
			m_dataIDOfZyy[i].imagPart = -1;
		}

	}
}

// Output calculated values of impedance tensors
void ObservedDataStationNMT2::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		//OutputFiles::m_csvFile << m_stationID << "," << m_freq[*itr] << "," << 
		//	m_ZxxCalculated[icount].real() << "," << m_ZxxCalculated[icount].imag() << "," << 
		//	m_ZxyCalculated[icount].real() << "," << m_ZxyCalculated[icount].imag() << "," << 
		//	m_ZyxCalculated[icount].real() << "," << m_ZyxCalculated[icount].imag() << "," << 
		//	m_ZyyCalculated[icount].real() << "," << m_ZyyCalculated[icount].imag() << "," << 
		//	m_ZxxResidual[icount].realPart << "," << m_ZxxResidual[icount].imagPart << "," <<
		//	m_ZxyResidual[icount].realPart << "," << m_ZxyResidual[icount].imagPart << "," <<
		//	m_ZyxResidual[icount].realPart << "," << m_ZyxResidual[icount].imagPart << "," <<
		//	m_ZyyResidual[icount].realPart << "," << m_ZyyResidual[icount].imagPart << "," <<
		//	m_ZxxObserved[*itr].real() << "," << m_ZxxObserved[*itr].imag() << "," << 
		//	m_ZxyObserved[*itr].real() << "," << m_ZxyObserved[*itr].imag() << "," << 
		//	m_ZyxObserved[*itr].real() << "," << m_ZyxObserved[*itr].imag() << "," << 
		//	m_ZyyObserved[*itr].real() << "," << m_ZyyObserved[*itr].imag() << "," << 
		//	m_ZxxSD[*itr].realPart << "," << m_ZxxSD[*itr].imagPart << "," << 
		//	m_ZxySD[*itr].realPart << "," << m_ZxySD[*itr].imagPart << "," << 
		//	m_ZyxSD[*itr].realPart << "," << m_ZyxSD[*itr].imagPart << "," << 
		//	m_ZyySD[*itr].realPart << "," << m_ZyySD[*itr].imagPart << ","<< std::endl;	

		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_stationID, m_freq[*itr], m_ZxxCalculated[icount].real(), m_ZxxCalculated[icount].imag(), m_ZxyCalculated[icount].real(), m_ZxyCalculated[icount].imag(),
			m_ZyxCalculated[icount].real(), m_ZyxCalculated[icount].imag(), m_ZyyCalculated[icount].real(), m_ZyyCalculated[icount].imag());
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_ZxxResidual[icount].realPart, m_ZxxResidual[icount].imagPart, m_ZxyResidual[icount].realPart, m_ZxyResidual[icount].imagPart,
			m_ZyxResidual[icount].realPart, m_ZyxResidual[icount].imagPart, m_ZyyResidual[icount].realPart, m_ZyyResidual[icount].imagPart );
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_ZxxObserved[*itr].real(),m_ZxxObserved[*itr].imag(),m_ZxyObserved[*itr].real(),m_ZxyObserved[*itr].imag(),
			m_ZyxObserved[*itr].real(),m_ZyxObserved[*itr].imag(),m_ZyyObserved[*itr].real(),m_ZyyObserved[*itr].imag() );
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n", 
			m_ZxxSD[*itr].realPart,m_ZxxSD[*itr].imagPart,m_ZxySD[*itr].realPart,m_ZxySD[*itr].imagPart,
			m_ZyxSD[*itr].realPart,m_ZyxSD[*itr].imagPart,m_ZyySD[*itr].realPart,m_ZyySD[*itr].imagPart );

		++icount;
	}

}

// Calulate interpolator vector of voltage difference
void ObservedDataStationNMT2::calcInterpolatorVectorOfVoltageDifference( Forward3D* const ptrForward3D ){

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh
		for( int iDipole = 0; iDipole < 2; ++iDipole ){
			ptrForward3D->calcInterpolatorVectorOfVoltageDifference( m_numElementsIncludingDipole[iDipole], m_elementsIncludingDipole[iDipole], m_facesIncludingDipole[iDipole],
				m_areaCoordinateValuesStartPoint[iDipole], m_areaCoordinateValuesEndPoint[iDipole], m_rhsVectorIDOfVoltageDifference[iDipole] );
		}
	}else{// Hexa mesh
		for( int iDipole = 0; iDipole < 2; ++iDipole ){
			ptrForward3D->calcInterpolatorVectorOfVoltageDifference( m_numElementsIncludingDipole[iDipole], m_elementsIncludingDipole[iDipole],
				m_localCoordinateValuesStartPoint[iDipole], m_localCoordinateValuesEndPoint[iDipole], m_rhsVectorIDOfVoltageDifference[iDipole] );
		}
	}

}


// Calulate sensitivity matrix of Impedance tensors
void ObservedDataStationNMT2::calculateSensitivityMatrix( const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMFieldExPol,
	const std::complex<double>* const derivativesOfEMFieldEyPol,
	double* const sensitivityMatrix, const bool forceSDToOne ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];
	const std::complex<double> divDet = std::complex<double>(1.0,0.0) / det;
	const std::complex<double> divDet2 = divDet * divDet;

	const long long rhsVectorIDOfHx = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	const long long rhsVectorIDOfHy = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());

	const double lengthOfDipole1 = hypot( m_location[0].endPoint.X - m_location[0].startPoint.X, m_location[0].endPoint.Y - m_location[0].startPoint.Y );
	const double lengthOfDipole2 = hypot( m_location[1].endPoint.X - m_location[1].startPoint.X, m_location[1].endPoint.Y - m_location[1].startPoint.Y );

	const double angleOfDipole1 = atan2( m_location[0].endPoint.Y - m_location[0].startPoint.Y , m_location[0].endPoint.X - m_location[0].startPoint.X );
	const double angleOfDipole2 = atan2( m_location[1].endPoint.Y - m_location[1].startPoint.Y , m_location[1].endPoint.X - m_location[1].startPoint.X );

	const std::complex<double> ElectricFieldDipole1[2] = { m_voltageCalculated[0][0] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 ) , m_voltageCalculated[0][1] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 ) };
	const std::complex<double> ElectricFieldDipole2[2] = { m_voltageCalculated[1][0] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 ) , m_voltageCalculated[1][1] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 ) };

	std::complex<double> ExCalculated[2];
	std::complex<double> EyCalculated[2];
	for( int iPol = 0; iPol < 2 ; ++iPol ){
		ExCalculated[iPol] =   ElectricFieldDipole1[iPol] * static_cast< std::complex<double> >( sin( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) - ElectricFieldDipole2[iPol] * static_cast< std::complex<double> >( sin( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) );
		EyCalculated[iPol] = - ElectricFieldDipole1[iPol] * static_cast< std::complex<double> >( cos( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) + ElectricFieldDipole2[iPol] * static_cast< std::complex<double> >( cos( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) );
	}

	CommonParameters::DoubleComplexValues ZxxSD = m_ZxxSD[freqIDGlobalInSta];
	CommonParameters::DoubleComplexValues ZxySD = m_ZxySD[freqIDGlobalInSta];
	CommonParameters::DoubleComplexValues ZyxSD = m_ZyxSD[freqIDGlobalInSta];
	CommonParameters::DoubleComplexValues ZyySD = m_ZyySD[freqIDGlobalInSta];
	if( forceSDToOne ){
		// force erro to one
		ZxxSD.realPart = 1.0;
		ZxySD.realPart = 1.0;
		ZyxSD.realPart = 1.0;
		ZyySD.realPart = 1.0;
		ZxxSD.imagPart = 1.0;
		ZxySD.imagPart = 1.0;
		ZyxSD.imagPart = 1.0;
		ZyySD.imagPart = 1.0;
	}

	const long long nBlkNotFixed = static_cast<long long>( ( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed() );
	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const std::complex<double> dEdmDipole1[2] = {
			derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference[0]) + imdl ] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 ) ,
			derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference[0]) + imdl ] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole1 )
		};

		const std::complex<double> dEdmDipole2[2] = {
			derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference[1]) + imdl ] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 ) ,
			derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference[1]) + imdl ] * static_cast< std::complex<double> >( -1.0 / lengthOfDipole2 )
		};
		
		const std::complex<double> dExdm[2] = {
			dEdmDipole1[0] * static_cast< std::complex<double> >( sin( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) - dEdmDipole2[0] * static_cast< std::complex<double> >( sin( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) ),
			dEdmDipole1[1] * static_cast< std::complex<double> >( sin( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) - dEdmDipole2[1] * static_cast< std::complex<double> >( sin( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) )
		};

		const std::complex<double> dEydm[2] = {
			- dEdmDipole1[0] * static_cast< std::complex<double> >( cos( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) + dEdmDipole2[0] * static_cast< std::complex<double> >( cos( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) ),
			- dEdmDipole1[1] * static_cast< std::complex<double> >( cos( angleOfDipole2 ) / sin( angleOfDipole2 - angleOfDipole1 ) ) + dEdmDipole2[1] * static_cast< std::complex<double> >( cos( angleOfDipole1 ) / sin( angleOfDipole2 - angleOfDipole1 ) )
		};
		
		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] *   HyCalculated[1]
										 + derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] *   HxCalculated[0]
										 - derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] *   HyCalculated[0]
										 - derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] *   HxCalculated[1];

		// dZxx/dm
		const std::complex<double> workXX1	= dExdm[0] * HyCalculated[1]
											+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * ExCalculated[0]
											- dExdm[1] * HyCalculated[0]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * ExCalculated[1];
		
		const std::complex<double> workXX2	= ExCalculated[0]*HyCalculated[1] - ExCalculated[1]*HyCalculated[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).real() / ZxxSD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).imag() / ZxxSD.imagPart;

		// dZxy/dm
		const std::complex<double> workXY1	= dExdm[1] * HxCalculated[0]
											+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * ExCalculated[1]
											- dExdm[0] * HxCalculated[1]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * ExCalculated[0];

		const std::complex<double> workXY2	= ExCalculated[1]*HxCalculated[0] - ExCalculated[0]*HxCalculated[1];
		 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).real() / ZxySD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).imag() / ZxySD.imagPart;

		// dZyx/dm
		const std::complex<double> workYX1	= dEydm[0] * HyCalculated[1]
											+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * EyCalculated[0]
											- dEydm[1] * HyCalculated[0]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * EyCalculated[1];

		const std::complex<double> workYX2	= EyCalculated[0]*HyCalculated[1] - EyCalculated[1]*HyCalculated[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + imdl ] = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).real() / ZyxSD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + imdl ] = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).imag() / ZyxSD.imagPart;

		// dZyy/dm
		const std::complex<double> workYY1	= dEydm[1] * HxCalculated[0]
											+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * EyCalculated[1]
											- dEydm[0] * HxCalculated[1]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * EyCalculated[0];

		const std::complex<double> workYY2	= EyCalculated[1]*HxCalculated[0] - EyCalculated[0]*HxCalculated[1];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + imdl ] = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).real() / ZyySD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + imdl ] = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).imag() / ZyySD.imagPart;
	}
	
}

// Calculate data vector of this PE
void ObservedDataStationNMT2::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	vector[ offset + m_dataIDOfZxx[freqIDThisPEInSta].realPart ] = m_ZxxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfZxx[freqIDThisPEInSta].imagPart ] = m_ZxxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfZxy[freqIDThisPEInSta].realPart ] = m_ZxyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfZxy[freqIDThisPEInSta].imagPart ] = m_ZxyResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfZyx[freqIDThisPEInSta].realPart ] = m_ZyxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfZyx[freqIDThisPEInSta].imagPart ] = m_ZyxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfZyy[freqIDThisPEInSta].realPart ] = m_ZyyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfZyy[freqIDThisPEInSta].imagPart ] = m_ZyyResidual[freqIDThisPEInSta].imagPart;
}

// Calulate L2 norm of misfit
double ObservedDataStationNMT2::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_ZxxResidual[ifreq].realPart , 2 );
		misfit += pow( m_ZxxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_ZxyResidual[ifreq].realPart , 2 );
		misfit += pow( m_ZxyResidual[ifreq].imagPart , 2 );
		misfit += pow( m_ZyxResidual[ifreq].realPart , 2 );
		misfit += pow( m_ZyxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_ZyyResidual[ifreq].realPart , 2 );
		misfit += pow( m_ZyyResidual[ifreq].imagPart , 2 );
	}

	return misfit;

}

// Get location of the station
const CommonParameters::locationDipole& ObservedDataStationNMT2::getLocationOfStation( const int iDipole ) const{

	assert( iDipole == 0 || iDipole == 1 );

	return m_location[ iDipole ];

}

// Get Z coordinate of the point
double ObservedDataStationNMT2::getZCoordOfPoint( const int iDipole , const int num ) const{

	assert( iDipole == 0 || iDipole == 1 );
	assert( num == 0 || num == 1 );

	int iElem(-1);
	int iFace(-1);
	CommonParameters::AreaCoords areaCoord = { -1.0, -1.0, -1.0 };

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh

		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();

		double dummy(0.0);
		if( num == 0 ){
			iElem = ptrMeshDataTetraElement->findElementIncludingPointOnSurface( m_location[iDipole].startPoint.X, m_location[iDipole].startPoint.Y, iFace, areaCoord, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataTetraElement->findElementIncludingPointOnSurface( m_location[iDipole].endPoint.X, m_location[iDipole].endPoint.Y, iFace, areaCoord, false, false, dummy, dummy );
		}
		return ptrMeshDataTetraElement->calcZCoordOfPointOnFace( iElem, iFace, areaCoord );

	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::NONCONFORMING_HEXA ){

		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();

		int faceID(0);
		double dummy(0.0);
		double localCoordX(0.0);
		double localCoordY(0.0);
		double localCoordZ(0.0);
		if( num == 0 ){
			iElem = ptrMeshDataNonConformingHexaElement->findElementIncludingPointOnSurface( 
				m_location[iDipole].startPoint.X, m_location[iDipole].startPoint.Y, faceID, localCoordX, localCoordY, localCoordZ, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataNonConformingHexaElement->findElementIncludingPointOnSurface( 
				m_location[iDipole].endPoint.X, m_location[iDipole].endPoint.Y, faceID, localCoordX, localCoordY, localCoordZ, false, false, dummy, dummy );
		}
		return ptrMeshDataNonConformingHexaElement->calcZCoordOfPointOnFace(iElem, 4, localCoordX, localCoordY);

	}
	else{

		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();

		double dummy(0.0);
		double localCoordZ(0.0);
		if( num == 0 ){
			iElem = ptrMeshDataBrickElement->findElementIncludingPointOnSurface( m_location[iDipole].startPoint.X, m_location[iDipole].startPoint.Y, dummy, dummy, localCoordZ, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataBrickElement->findElementIncludingPointOnSurface( m_location[iDipole].endPoint.X,   m_location[iDipole].endPoint.Y, dummy, dummy, localCoordZ, false, false, dummy, dummy );
		}
		return ptrMeshDataBrickElement->calcGlobalCoordZ( iElem, localCoordZ );

	}

}

