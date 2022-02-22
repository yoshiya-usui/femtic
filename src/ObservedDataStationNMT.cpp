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

#include "ObservedDataStationNMT.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationNMT::ObservedDataStationNMT():
	ObservedDataStation(),
	m_YxObserved(NULL),
	m_YyObserved(NULL),
	m_YxSD(NULL),
	m_YySD(NULL),
	m_YxCalculated(NULL),
	m_YyCalculated(NULL),
	m_YxResidual(NULL),
	m_YyResidual(NULL),
	m_numElementsIncludingDipole(0),
	m_elementsIncludingDipole(NULL),
	m_facesIncludingDipole(NULL),
	m_localCoordinateValuesStartPoint(NULL),
	m_localCoordinateValuesEndPoint(NULL),
	m_areaCoordinateValuesStartPoint(NULL),
	m_areaCoordinateValuesEndPoint(NULL),
	m_rhsVectorIDOfVoltageDifference(0),
	m_dataIDOfYx(NULL),
	m_dataIDOfYy(NULL)
{

	m_location.startPoint.X = 0.0;
	m_location.startPoint.Y = 0.0;
	m_location.endPoint.X = 0.0;
	m_location.endPoint.Y = 0.0;
	
	for( int i = 0; i < 2; ++i ){
		m_voltageCalculated[i] = std::complex<double>(0.0,0.0);
	}

}

// Destructer
ObservedDataStationNMT::~ObservedDataStationNMT(){

	if( m_YxObserved != NULL){
		delete[] m_YxObserved;
		m_YxObserved = NULL;
	}

	if( m_YyObserved != NULL){
		delete[] m_YyObserved;
		m_YyObserved = NULL;
	}

	if( m_YxSD != NULL){
		delete[] m_YxSD;
		m_YxSD = NULL;
	}

	if( m_YySD != NULL){
		delete[] m_YySD;
		m_YySD = NULL;
	}

	if( m_YxCalculated != NULL){
		delete[] m_YxCalculated;
		m_YxCalculated = NULL;
	}

	if( m_YyCalculated != NULL){
		delete[] m_YyCalculated;
		m_YyCalculated = NULL;
	}

	if( m_YxResidual != NULL){
		delete[] m_YxResidual;
		m_YxResidual = NULL;
	}

	if( m_YyResidual != NULL){
		delete[] m_YyResidual;
		m_YyResidual = NULL;
	}

	if( m_elementsIncludingDipole != NULL){
		delete[] m_elementsIncludingDipole;
		m_elementsIncludingDipole = NULL;
	}

	if( m_facesIncludingDipole != NULL){
		delete[] m_facesIncludingDipole;
		m_facesIncludingDipole = NULL;
	}

	if( m_localCoordinateValuesStartPoint != NULL){
		delete[] m_localCoordinateValuesStartPoint;
		m_localCoordinateValuesStartPoint = NULL;
	}

	if( m_localCoordinateValuesEndPoint != NULL){
		delete[] m_localCoordinateValuesEndPoint;
		m_localCoordinateValuesEndPoint = NULL;
	}

	if( m_areaCoordinateValuesStartPoint != NULL){
		delete[] m_areaCoordinateValuesStartPoint;
		m_areaCoordinateValuesStartPoint = NULL;
	}

	if( m_areaCoordinateValuesEndPoint != NULL){
		delete[] m_areaCoordinateValuesEndPoint;
		m_areaCoordinateValuesEndPoint = NULL;
	}

	if( m_dataIDOfYx != NULL){
		delete[] m_dataIDOfYx;
		m_dataIDOfYx = NULL;
	}

	if( m_dataIDOfYy != NULL){
		delete[] m_dataIDOfYy;
		m_dataIDOfYy = NULL;
	}

}

// Read data from input file
void ObservedDataStationNMT::inputObservedData( std::ifstream& inFile ){

	inFile >> m_stationID;

	inFile >> m_IDOfMagneticFieldStation;

	OutputFiles::m_logFile << "# " << std::setw(15) << std::left << m_stationID << std::setw(15) << std::left << m_IDOfMagneticFieldStation << std::endl;

	double dbuf(0.0);
	inFile >> dbuf;
	m_location.startPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.startPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.endPoint.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.endPoint.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> m_numOfFrequency;
	const int nFreq = m_numOfFrequency;
	if( nFreq > 0 ){
		m_freq = new double[nFreq];
		m_YxObserved = new std::complex<double>[nFreq];
		m_YyObserved = new std::complex<double>[nFreq];
		m_YxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_YySD = new CommonParameters::DoubleComplexValues[nFreq];
		//m_YxCalculated = new std::complex<double>[nFreq];
		//m_YyCalculated = new std::complex<double>[nFreq];
		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1 >> dbuf2;
			m_YxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_YyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> m_YxSD[i].realPart;
			inFile >> m_YxSD[i].imagPart;
			inFile >> m_YySD[i].realPart;
			inFile >> m_YySD[i].imagPart;
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << " NMT " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.startPoint.X << " "
			    << m_location.startPoint.Y << " "
				<< m_location.endPoint.X << " "
				<< m_location.endPoint.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for( int i = 0; i < m_numOfFrequency; ++i ){
		std::cout << m_freq[i] << " "
				    << m_YxObserved[i] << " "
				    << m_YyObserved[i] << " "
					<< m_YxSD[i].realPart << " "
					<< m_YxSD[i].imagPart << " "
					<< m_YySD[i].realPart << " "
					<< m_YySD[i].imagPart << std::endl;
	}
#endif

}

// Find elements including dipole
void ObservedDataStationNMT::findElementsIncludingDipole(){

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra/Prism mesh

		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();

		std::vector<CommonParameters::AreaCoords> localCoordStartPoint;
		std::vector<CommonParameters::AreaCoords> localCoordEndPoint;

		std::vector<int> elementsIncludingDipole;
		std::vector<int> facesIncludingDipole;

		ptrMeshDataTetraElement->findElementsIncludingDipoleOnSurface( m_location.startPoint.X, m_location.startPoint.Y, m_location.endPoint.X, m_location.endPoint.Y,
			elementsIncludingDipole, facesIncludingDipole, localCoordStartPoint, localCoordEndPoint );

		m_numElementsIncludingDipole = static_cast<int>( elementsIncludingDipole.size() );
		if( m_numElementsIncludingDipole > 0 ){
			m_elementsIncludingDipole = new int[m_numElementsIncludingDipole];
			m_facesIncludingDipole = new int[m_numElementsIncludingDipole];
			m_areaCoordinateValuesStartPoint = new CommonParameters::AreaCoords[m_numElementsIncludingDipole];
			m_areaCoordinateValuesEndPoint = new CommonParameters::AreaCoords[m_numElementsIncludingDipole];
		}

		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			m_elementsIncludingDipole[i] = elementsIncludingDipole[i];
			m_facesIncludingDipole[i] = facesIncludingDipole[i];
			m_areaCoordinateValuesStartPoint[i] = localCoordStartPoint[i];
			m_areaCoordinateValuesEndPoint[i] = localCoordEndPoint[i];
		}

#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			std::cout << m_elementsIncludingDipole[i] << " " << m_facesIncludingDipole[i] << " " << m_areaCoordinateValuesStartPoint[i].coord0  << " " << m_areaCoordinateValuesStartPoint[i].coord1 << " " << m_areaCoordinateValuesStartPoint[i].coord2 << std::endl;
			std::cout << m_elementsIncludingDipole[i] << " " << m_facesIncludingDipole[i] << " " << m_areaCoordinateValuesStartPoint[i].coord0  << " " << m_areaCoordinateValuesStartPoint[i].coord1 << " " << m_areaCoordinateValuesStartPoint[i].coord2 << std::endl;
		}
#endif
	}else if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::NONCONFORMING_HEXA ){// Non-conforming deformed hexahedral mesh

		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataNonConformingHexaElement();

		std::vector<double> localCoordXStartPoint;
		std::vector<double> localCoordYStartPoint;
		std::vector<double> localCoordXEndPoint;
		std::vector<double> localCoordYEndPoint;
		std::vector<int> elementsIncludingDipole;

		ptrMeshDataNonConformingHexaElement->findElementsIncludingDipoleOnSurface( m_location.startPoint.X, m_location.startPoint.Y, m_location.endPoint.X, m_location.endPoint.Y,
			elementsIncludingDipole, localCoordXStartPoint, localCoordYStartPoint, localCoordXEndPoint, localCoordYEndPoint );

		m_numElementsIncludingDipole = static_cast<int>( elementsIncludingDipole.size() );
		if( m_numElementsIncludingDipole > 0 ){
			m_elementsIncludingDipole = new int[m_numElementsIncludingDipole];
			m_localCoordinateValuesStartPoint = new CommonParameters::locationXY[m_numElementsIncludingDipole];
			m_localCoordinateValuesEndPoint   = new CommonParameters::locationXY[m_numElementsIncludingDipole];
		}
	
		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			m_elementsIncludingDipole[i] = elementsIncludingDipole[i];
			m_localCoordinateValuesStartPoint[i].X = localCoordXStartPoint[i];
			m_localCoordinateValuesStartPoint[i].Y = localCoordYStartPoint[i];
			m_localCoordinateValuesEndPoint[i].X = localCoordXEndPoint[i];
			m_localCoordinateValuesEndPoint[i].Y = localCoordYEndPoint[i];
		}

#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		std::cout << "m_elementsIncludingDipole m_localCoordinateValuesStartPoint[i].x m_localCoordinateValuesStartPoint[i].y m_localCoordinateValuesEndPoint[i].x m_localCoordinateValuesEndPoint[i].y" << std::endl;
		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			std::cout << m_elementsIncludingDipole[i] << " " << m_localCoordinateValuesStartPoint[i].X  << " " << m_localCoordinateValuesStartPoint[i].Y << " " << m_localCoordinateValuesEndPoint[i].X << " " << m_localCoordinateValuesEndPoint[i].Y << std::endl;
		}
#endif
	}else{// Hexa mesh

		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();

		std::vector<double> localCoordXStartPoint;
		std::vector<double> localCoordYStartPoint;
		std::vector<double> localCoordXEndPoint;
		std::vector<double> localCoordYEndPoint;

		std::vector<int> elementsIncludingDipole;
		ptrMeshDataBrickElement->findElementsIncludingDipoleOnSurface( m_location.startPoint.X, m_location.startPoint.Y, m_location.endPoint.X, m_location.endPoint.Y,
			elementsIncludingDipole, localCoordXStartPoint, localCoordYStartPoint, localCoordXEndPoint, localCoordYEndPoint );

		m_numElementsIncludingDipole = static_cast<int>( elementsIncludingDipole.size() );
		if( m_numElementsIncludingDipole > 0 ){
			m_elementsIncludingDipole = new int[m_numElementsIncludingDipole];
			m_localCoordinateValuesStartPoint = new CommonParameters::locationXY[m_numElementsIncludingDipole];
			m_localCoordinateValuesEndPoint   = new CommonParameters::locationXY[m_numElementsIncludingDipole];
		}
	
		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			m_elementsIncludingDipole[i] = elementsIncludingDipole[i];
			m_localCoordinateValuesStartPoint[i].X = localCoordXStartPoint[i];
			m_localCoordinateValuesStartPoint[i].Y = localCoordYStartPoint[i];
			m_localCoordinateValuesEndPoint[i].X = localCoordXEndPoint[i];
			m_localCoordinateValuesEndPoint[i].Y = localCoordYEndPoint[i];
		}

#ifdef _DEBUG_WRITE
		std::cout << "m_stationID " << m_stationID << std::endl;
		std::cout << "m_elementsIncludingDipole m_localCoordinateValuesStartPoint[i].x m_localCoordinateValuesStartPoint[i].y m_localCoordinateValuesEndPoint[i].x m_localCoordinateValuesEndPoint[i].y" << std::endl;
		for( int i = 0; i < m_numElementsIncludingDipole; ++i ){
			std::cout << m_elementsIncludingDipole[i] << " " << m_localCoordinateValuesStartPoint[i].X  << " " << m_localCoordinateValuesStartPoint[i].Y << " " << m_localCoordinateValuesEndPoint[i].X << " " << m_localCoordinateValuesEndPoint[i].Y << std::endl;
		}
#endif

	}

}

// Calulate difference of voltages
void ObservedDataStationNMT::calculateVoltageDifferences( const Forward3D* const ptrForward3D, const int rhsVectorIDOfVoltageDifference ){

	const int iPol = ptrForward3D->getPolarizationCurrent();

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh
		m_voltageCalculated[iPol] = ptrForward3D->calcVoltageDifference( m_numElementsIncludingDipole, m_elementsIncludingDipole, m_facesIncludingDipole, m_areaCoordinateValuesStartPoint, m_areaCoordinateValuesEndPoint );
	}else{// Hexa mesh
		m_voltageCalculated[iPol] = ptrForward3D->calcVoltageDifference( m_numElementsIncludingDipole, m_elementsIncludingDipole, m_localCoordinateValuesStartPoint, m_localCoordinateValuesEndPoint );
	}

#ifdef _DEBUG_WRITE
	std::cout << "iPol V : " << iPol << " " << m_voltageCalculated[iPol] << std::endl;
#endif

	// For inversion
	m_rhsVectorIDOfVoltageDifference = rhsVectorIDOfVoltageDifference;
}

// Calulate Network-MT response
//void ObservedDataStationNMT::calculateNetworkMTResponse( const int freqIDAmongThisPE, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){
void ObservedDataStationNMT::calculateNetworkMTResponse( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];

	m_YxCalculated[freqIDThisPEInSta] = ( m_voltageCalculated[0]*HyCalculated[1] - m_voltageCalculated[1]*HyCalculated[0] ) / det;
	m_YyCalculated[freqIDThisPEInSta] = ( m_voltageCalculated[1]*HxCalculated[0] - m_voltageCalculated[0]*HxCalculated[1] ) / det;

	m_YxResidual[freqIDThisPEInSta].realPart = ( m_YxObserved[freqIDGlobalInSta].real() - m_YxCalculated[freqIDThisPEInSta].real() ) / m_YxSD[freqIDGlobalInSta].realPart;
	m_YxResidual[freqIDThisPEInSta].imagPart = ( m_YxObserved[freqIDGlobalInSta].imag() - m_YxCalculated[freqIDThisPEInSta].imag() ) / m_YxSD[freqIDGlobalInSta].imagPart;
	m_YyResidual[freqIDThisPEInSta].realPart = ( m_YyObserved[freqIDGlobalInSta].real() - m_YyCalculated[freqIDThisPEInSta].real() ) / m_YySD[freqIDGlobalInSta].realPart;
	m_YyResidual[freqIDThisPEInSta].imagPart = ( m_YyObserved[freqIDGlobalInSta].imag() - m_YyCalculated[freqIDThisPEInSta].imag() ) / m_YySD[freqIDGlobalInSta].imagPart;

#ifdef _DEBUG_WRITE
	std::cout << "ifreq Yx Yy: " << freqIDThisPEInSta << " " << m_YxCalculated[freqIDThisPEInSta] << " " << m_YyCalculated[freqIDThisPEInSta]<< std::endl;
#endif

	// For inversion
	//ObservedData* const ptrObservedData = ObservedData::getInstance();
	//m_dataIDOfYx[freqIDThisPEInSta] = ptrObservedData->incrementNumObservedDataThisPE( freqIDThisPEInSta );
	//m_dataIDOfYy[freqIDThisPEInSta] = ptrObservedData->incrementNumObservedDataThisPE( freqIDThisPEInSta );

	m_dataIDOfYx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfYx[freqIDThisPEInSta].imagPart = icount++;
	m_dataIDOfYy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfYy[freqIDThisPEInSta].imagPart = icount++;
}

// Initialize difference of voltages
void ObservedDataStationNMT::initializeVoltageDifferences( const int iPol ){

	m_voltageCalculated[iPol] = std::complex<double>(0.0,0.0);

}

// Initialize Network-MT responses and errors
void ObservedDataStationNMT::initializeNetworkMTResponsesAndErrors(){
	
	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_YxResidual[i].realPart = 0.0;
		m_YxResidual[i].imagPart = 0.0;
		m_YyResidual[i].realPart = 0.0;
		m_YyResidual[i].imagPart = 0.0;

		m_YxCalculated[i] = std::complex<double>(0.0,0.0);
		m_YyCalculated[i] = std::complex<double>(0.0,0.0);
	}

}

// Allocate memory for the calculated values of Network-MT responses and errors
void ObservedDataStationNMT::allocateMemoryForCalculatedValues(){
	
	if( m_YxCalculated != NULL){
		delete[] m_YxCalculated;
		m_YxCalculated = NULL;
	}

	if( m_YyCalculated != NULL){
		delete[] m_YyCalculated;
		m_YyCalculated = NULL;
	}

	if( m_YxResidual != NULL){
		delete[] m_YxResidual;
		m_YxResidual = NULL;
	}

	if( m_YyResidual != NULL){
		delete[] m_YyResidual;
		m_YyResidual = NULL;
	}

	if( m_dataIDOfYx != NULL){
		delete[] m_dataIDOfYx;
		m_dataIDOfYx = NULL;
	}

	if( m_dataIDOfYy != NULL){
		delete[] m_dataIDOfYy;
		m_dataIDOfYy = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){
		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_YxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_YyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_YxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_YyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfYx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfYy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_YxCalculated[i] = std::complex<double>(0.0,0.0);
			m_YyCalculated[i] = std::complex<double>(0.0,0.0);;
			m_YxResidual[i].realPart = 0.0;
			m_YxResidual[i].imagPart = 0.0;
			m_YyResidual[i].realPart = 0.0;
			m_YyResidual[i].imagPart = 0.0;
			m_dataIDOfYx[i].realPart = -1;
			m_dataIDOfYx[i].imagPart = -1;
			m_dataIDOfYy[i].realPart = -1;
			m_dataIDOfYy[i].imagPart = -1;
		}

	}

}

// Output calculated values of Network-MT responses
void ObservedDataStationNMT::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		//OutputFiles::m_csvFile << m_stationID << "," << m_freq[*itr] << "," << 
		//	m_YxCalculated[icount].real() << "," << m_YxCalculated[icount].imag() << "," << 
		//	m_YyCalculated[icount].real() << "," << m_YyCalculated[icount].imag() << "," << 
		//	m_YxResidual[icount].realPart << "," << m_YxResidual[icount].imagPart << "," << 
		//	m_YyResidual[icount].realPart << "," << m_YyResidual[icount].imagPart << "," << 
		//	m_YxObserved[*itr].real() << "," << m_YxObserved[*itr].imag() << "," << 
		//	m_YyObserved[*itr].real() << "," << m_YyObserved[*itr].imag() << "," << 
		//	m_YxSD[*itr].realPart << "," << m_YxSD[*itr].imagPart << "," << 
		//	m_YySD[*itr].realPart << "," << m_YySD[*itr].imagPart << std::endl;	

		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_stationID, m_freq[*itr], m_YxCalculated[icount].real(), m_YxCalculated[icount].imag(), m_YyCalculated[icount].real(), m_YyCalculated[icount].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_YxResidual[icount].realPart, m_YxResidual[icount].imagPart, m_YyResidual[icount].realPart, m_YyResidual[icount].imagPart ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_YxObserved[*itr].real() , m_YxObserved[*itr].imag(), m_YyObserved[*itr].real(), m_YyObserved[*itr].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_YxSD[*itr].realPart, m_YxSD[*itr].imagPart, m_YySD[*itr].realPart, m_YySD[*itr].imagPart ); 

		++icount;
	}
}

// Calulate interpolator vector of voltage difference
void ObservedDataStationNMT::calcInterpolatorVectorOfVoltageDifference( Forward3D* const ptrForward3D ){

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh
		ptrForward3D->calcInterpolatorVectorOfVoltageDifference( m_numElementsIncludingDipole, m_elementsIncludingDipole, m_facesIncludingDipole, m_areaCoordinateValuesStartPoint, m_areaCoordinateValuesEndPoint, m_rhsVectorIDOfVoltageDifference );
	}else{// Hexa mesh
		ptrForward3D->calcInterpolatorVectorOfVoltageDifference( m_numElementsIncludingDipole, m_elementsIncludingDipole, m_localCoordinateValuesStartPoint, m_localCoordinateValuesEndPoint, m_rhsVectorIDOfVoltageDifference );
	}

}

// Calulate sensitivity matrix of Impedance tensors
void ObservedDataStationNMT::calculateSensitivityMatrix( const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMFieldExPol,
	const std::complex<double>* const derivativesOfEMFieldEyPol,
	double* const sensitivityMatrix ) const{

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
	const long long nBlkNotFixed = static_cast<long long>(( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed());

	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * HyCalculated[1]
									+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * HxCalculated[0]
									- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * HyCalculated[0]
									- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * HxCalculated[1];

		// dYx/dm
		const std::complex<double> workX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference) + imdl ] * HyCalculated[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * m_voltageCalculated[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference) + imdl ] * HyCalculated[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * m_voltageCalculated[1];
		
		const std::complex<double> workX2	= m_voltageCalculated[0]*HyCalculated[1] - m_voltageCalculated[1]*HyCalculated[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfYx[freqIDThisPEInSta].realPart) + imdl ] = ( workX1 * divDet - work1 * workX2 * divDet2 ).real() / m_YxSD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfYx[freqIDThisPEInSta].imagPart) + imdl ] = ( workX1 * divDet - work1 * workX2 * divDet2 ).imag() / m_YxSD[freqIDGlobalInSta].imagPart;

		// dYy/dm
		const std::complex<double> workY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference) + imdl ] * HxCalculated[0]
										+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * m_voltageCalculated[1]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfVoltageDifference) + imdl ] * HxCalculated[1]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * m_voltageCalculated[0];

		const std::complex<double> workY2	= m_voltageCalculated[1]*HxCalculated[0] - m_voltageCalculated[0]*HxCalculated[1];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfYy[freqIDThisPEInSta].realPart) + imdl ] = ( workY1 * divDet - work1 * workY2 * divDet2 ).real() / m_YySD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfYy[freqIDThisPEInSta].imagPart) + imdl ] = ( workY1 * divDet - work1 * workY2 * divDet2 ).imag() / m_YySD[freqIDGlobalInSta].imagPart;

	}

}

// Calculate data vector of this PE
void ObservedDataStationNMT::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	vector[ offset + m_dataIDOfYx[freqIDThisPEInSta].realPart ] = m_YxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfYx[freqIDThisPEInSta].imagPart ] = m_YxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfYy[freqIDThisPEInSta].realPart ] = m_YyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfYy[freqIDThisPEInSta].imagPart ] = m_YyResidual[freqIDThisPEInSta].imagPart;

}

// Calulate sum of square of misfit
double ObservedDataStationNMT::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_YxResidual[ifreq].realPart , 2 );
		misfit += pow( m_YxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_YyResidual[ifreq].realPart , 2 );
		misfit += pow( m_YyResidual[ifreq].imagPart , 2 );
	}

	return misfit;

}

// Get location of the station
const CommonParameters::locationDipole& ObservedDataStationNMT::getLocationOfStation() const{

	return m_location;

}

// Get Z coordinate of the point
double ObservedDataStationNMT::getZCoordOfPoint( const int num ) const{

	assert( num == 0 || num == 1 );

	int iElem(-1);
	int iFace(-1);
	CommonParameters::AreaCoords areaCoord = { -1.0, -1.0, -1.0 };

	if( ( AnalysisControl::getInstance() )->getTypeOfMesh() == MeshData::TETRA ){// Tetra mesh

		const MeshDataTetraElement* const ptrMeshDataTetraElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataTetraElement();

		double dummy(0.0);
		if( num == 0 ){
			iElem = ptrMeshDataTetraElement->findElementIncludingPointOnSurface( m_location.startPoint.X, m_location.startPoint.Y, iFace, areaCoord, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataTetraElement->findElementIncludingPointOnSurface( m_location.endPoint.X, m_location.endPoint.Y, iFace, areaCoord, false, false, dummy, dummy );
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
				m_location.startPoint.X, m_location.startPoint.Y, faceID, localCoordX, localCoordY, localCoordZ, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataNonConformingHexaElement->findElementIncludingPointOnSurface( 
				m_location.endPoint.X, m_location.endPoint.Y, faceID, localCoordX, localCoordY, localCoordZ, false, false, dummy, dummy );
		}
		return ptrMeshDataNonConformingHexaElement->calcZCoordOfPointOnFace(iElem, 4, localCoordX, localCoordY);

	}
	else{

		const MeshDataBrickElement* const ptrMeshDataBrickElement = ( AnalysisControl::getInstance() )->getPointerOfMeshDataBrickElement();

		double dummy(0.0);
		double localCoordZ(0.0);
		if( num == 0 ){
			iElem = ptrMeshDataBrickElement->findElementIncludingPointOnSurface( m_location.startPoint.X, m_location.startPoint.Y, dummy, dummy, localCoordZ, false, false, dummy, dummy );
		}else{
			iElem = ptrMeshDataBrickElement->findElementIncludingPointOnSurface( m_location.endPoint.X, m_location.endPoint.Y, dummy, dummy, localCoordZ, false, false, dummy, dummy );
		}
		return ptrMeshDataBrickElement->calcGlobalCoordZ( iElem, localCoordZ );

	}

}
