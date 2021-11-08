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

#include "ObservedDataStationVTF.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationVTF::ObservedDataStationVTF():
	ObservedDataStationPoint(),
	m_TzxObserved(NULL),
	m_TzyObserved(NULL),
	m_TzxSD(NULL),
	m_TzySD(NULL),
	m_TzxCalculated(NULL),
	m_TzyCalculated(NULL),
	m_TzxResidual(NULL),
	m_TzyResidual(NULL),
	m_rhsVectorIDOfHz(0),
	m_dataIDOfTzx(NULL),
	m_dataIDOfTzy(NULL)
{
	for( int i = 0; i < 2; ++i ){
		m_HzCalculated[i] = std::complex<double>(0.0,0.0);
	}
}

// Destructer
ObservedDataStationVTF::~ObservedDataStationVTF(){

	if( m_TzxObserved != NULL){
		delete[] m_TzxObserved;
		m_TzxObserved = NULL;
	}

	if( m_TzyObserved != NULL){
		delete[] m_TzyObserved;
		m_TzyObserved = NULL;
	}

	if( m_TzxSD != NULL){
		delete[] m_TzxSD;
		m_TzxSD = NULL;
	}

	if( m_TzySD != NULL){
		delete[] m_TzySD;
		m_TzySD = NULL;
	}

	if( m_TzxCalculated != NULL){
		delete[] m_TzxCalculated;
		m_TzxCalculated = NULL;
	}

	if( m_TzyCalculated != NULL){
		delete[] m_TzyCalculated;
		m_TzyCalculated = NULL;
	}

	if( m_TzxResidual != NULL){
		delete[] m_TzxResidual;
		m_TzxResidual = NULL;
	}

	if( m_dataIDOfTzx != NULL){
		delete[] m_dataIDOfTzx;
		m_dataIDOfTzx = NULL;
	}

	if( m_dataIDOfTzy != NULL){
		delete[] m_dataIDOfTzy;
		m_dataIDOfTzy = NULL;
	}

}

// Read data from input file
void ObservedDataStationVTF::inputObservedData( std::ifstream& inFile ){

	inFile >> m_stationID;

	inFile >> m_IDOfMagneticFieldStation;

	OutputFiles::m_logFile << "# " << std::setw(15) << std::left << m_stationID << std::setw(18) << std::left << m_IDOfMagneticFieldStation;

	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();
	// Type of owner element
	int ownerType(-1);
	std::string ownerElemType;
	if( pAnalysisControl->isTypeOfOwnerElementSetIndivisually() ){
		// Owner element type of each site is specified
		inFile >> ownerType;
	}
	else{
		ownerType = pAnalysisControl->getTypeOfOwnerElement();
	}
	switch(ownerType){
		case AnalysisControl::USE_LOWER_ELEMENT:
			m_useUpperElementForInterpolationOfEMField = false;
			ownerElemType = "Lower";
			break;
		case AnalysisControl::USE_UPPER_ELEMENT:
			m_useUpperElementForInterpolationOfEMField = true;
			ownerElemType = "Upper";
			break;
		default:
			OutputFiles::m_logFile << std::endl << "Error : Unknown type of owner element : " << ownerType << std::endl;
			exit(1);
			break;
	}
	OutputFiles::m_logFile << std::left << ownerElemType << std::endl;

	double dbuf(0.0);
	inFile >> dbuf;
	m_location.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> m_numOfFrequency;
	const int nFreq = m_numOfFrequency;
	if( nFreq > 0 ){
		m_freq = new double[nFreq];
		m_TzxObserved = new std::complex<double>[nFreq];
		m_TzyObserved = new std::complex<double>[nFreq];
		m_TzxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_TzySD = new CommonParameters::DoubleComplexValues[nFreq];
		//m_TzxCalculated = new std::complex<double>[nFreq];
		//m_TzyCalculated = new std::complex<double>[nFreq];
		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1 >> dbuf2;
			m_TzxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_TzyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> m_TzxSD[i].realPart;
			inFile >> m_TzxSD[i].imagPart;
			inFile >> m_TzySD[i].realPart;
			inFile >> m_TzySD[i].imagPart;
		}
	}

#ifdef _DEBUG_WRITE
	std::cout <<  " VTF " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.X << " " << m_location.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for( int i = 0; i < m_numOfFrequency; ++i ){
		std::cout << m_freq[i] << " "
				    << m_TzxObserved[i] << " "
				    << m_TzyObserved[i] << " "
					<< m_TzxSD[i].realPart << " "
					<< m_TzxSD[i].imagPart << " "
					<< m_TzySD[i].realPart << " "
					<< m_TzySD[i].imagPart << std::endl;
	}
#endif

}

// Calulate vertical magnetic field
void ObservedDataStationVTF::calculateVerticalMagneticField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfHz ){

	const int iPol = ptrForward3D->getPolarizationCurrent();

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		m_HzCalculated[iPol] = ptrForward3D->calcValueMagneticFieldZDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
	}else if( meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA ){
		m_HzCalculated[iPol] = ptrForward3D->calcValueMagneticFieldZDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

	// For inversion
	m_rhsVectorIDOfHz = rhsVectorIDOfHz;

}

// Calulate vertical magnetic field transfer function
//void ObservedDataStationVTF::calculateVTF( const int freqIDAmongThisPE ){
void ObservedDataStationVTF::calculateVTF( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];

	m_TzxCalculated[freqIDThisPEInSta] = ( m_HzCalculated[0]*HyCalculated[1] - m_HzCalculated[1]*HyCalculated[0] ) / det;
	m_TzyCalculated[freqIDThisPEInSta] = ( m_HzCalculated[1]*HxCalculated[0] - m_HzCalculated[0]*HxCalculated[1] ) / det;

	m_TzxResidual[freqIDThisPEInSta].realPart = ( m_TzxObserved[freqIDGlobalInSta].real() - m_TzxCalculated[freqIDThisPEInSta].real() ) / m_TzxSD[freqIDGlobalInSta].realPart;
	m_TzxResidual[freqIDThisPEInSta].imagPart = ( m_TzxObserved[freqIDGlobalInSta].imag() - m_TzxCalculated[freqIDThisPEInSta].imag() ) / m_TzxSD[freqIDGlobalInSta].imagPart;
	m_TzyResidual[freqIDThisPEInSta].realPart = ( m_TzyObserved[freqIDGlobalInSta].real() - m_TzyCalculated[freqIDThisPEInSta].real() ) / m_TzySD[freqIDGlobalInSta].realPart;
	m_TzyResidual[freqIDThisPEInSta].imagPart = ( m_TzyObserved[freqIDGlobalInSta].imag() - m_TzyCalculated[freqIDThisPEInSta].imag() ) / m_TzySD[freqIDGlobalInSta].imagPart;

#ifdef _DEBUG_WRITE
	std::cout << "freqIDThisPEInSta Tzx Tzy : " << freqIDThisPEInSta << " " << m_TzxCalculated[freqIDThisPEInSta] << " " << m_TzyCalculated[freqIDThisPEInSta] << std::endl;
#endif

	// For inversion
	//ObservedData* const ptrObservedData = ObservedData::getInstance();
	//m_dataIDOfTzx[freqIDAmongThisPE] = ptrObservedData->incrementNumObservedDataThisPE( ifreq );
	//m_dataIDOfTzy[freqIDAmongThisPE] = ptrObservedData->incrementNumObservedDataThisPE( ifreq );

	m_dataIDOfTzx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTzx[freqIDThisPEInSta].imagPart = icount++;
	m_dataIDOfTzy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTzy[freqIDThisPEInSta].imagPart = icount++;

}

// Initialize vertical magnetic field
void ObservedDataStationVTF::initializeVerticalMagneticField( const int iPol ){

	m_HzCalculated[iPol] = std::complex<double>(0.0,0.0);

}

// Initialize vertical magnetic field transfer functions and errors
void ObservedDataStationVTF::initializeVTFsAndErrors(){

	//for( int i = 0; i < m_numOfFrequency; ++i ){
	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_TzxResidual[i].realPart = 0.0;
		m_TzxResidual[i].imagPart = 0.0;
		m_TzyResidual[i].realPart = 0.0;
		m_TzyResidual[i].imagPart = 0.0;

		m_TzxCalculated[i] = std::complex<double>(0.0,0.0);
		m_TzyCalculated[i] = std::complex<double>(0.0,0.0);
	}

}

// Allocate memory for the calculated values of vertical magnetic field transfer functions and errors
void ObservedDataStationVTF::allocateMemoryForCalculatedValues(){

	if( m_TzxCalculated != NULL){
		delete[] m_TzxCalculated;
		m_TzxCalculated = NULL;
	}

	if( m_TzyCalculated != NULL){
		delete[] m_TzyCalculated;
		m_TzyCalculated = NULL;
	}

	if( m_TzxResidual != NULL){
		delete[] m_TzxResidual;
		m_TzxResidual = NULL;
	}

	if( m_TzyResidual != NULL){
		delete[] m_TzyResidual;
		m_TzyResidual = NULL;
	}

	if( m_dataIDOfTzx != NULL){
		delete[] m_dataIDOfTzx;
		m_dataIDOfTzx = NULL;
	}

	if( m_dataIDOfTzy != NULL){
		delete[] m_dataIDOfTzy;
		m_dataIDOfTzy = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){

		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_TzxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TzyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TzxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_TzyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//m_dataIDOfTzx = new int[m_numOfFreqCalculatedByThisStaAndPE];
		//m_dataIDOfTzy = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTzx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTzy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_TzxCalculated[i] = std::complex<double>(0.0,0.0);
			m_TzyCalculated[i] = std::complex<double>(0.0,0.0);
			m_TzxResidual[i].realPart = 0.0;
			m_TzxResidual[i].imagPart = 0.0;
			m_TzyResidual[i].realPart = 0.0;
			m_TzyResidual[i].imagPart = 0.0;
			m_dataIDOfTzx[i].realPart = -1;
			m_dataIDOfTzx[i].imagPart = -1;
			m_dataIDOfTzy[i].realPart = -1;
			m_dataIDOfTzy[i].imagPart = -1;
		}

	}

}

// Output calculated values of vertical magnetic field transfer functions
void ObservedDataStationVTF::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		//OutputFiles::m_csvFile << m_stationID << "," << m_freq[*itr] << "," << 
		//	m_TzxCalculated[icount].real() << "," << m_TzxCalculated[icount].imag() << "," << 
		//	m_TzyCalculated[icount].real() << "," << m_TzyCalculated[icount].imag() << "," << 
		//	m_TzxResidual[icount].realPart << "," << m_TzxResidual[icount].imagPart << "," << 
		//	m_TzyResidual[icount].realPart << "," << m_TzyResidual[icount].imagPart << "," << 
		//	m_TzxObserved[*itr].real() << "," << m_TzxObserved[*itr].imag() << "," << 
		//	m_TzyObserved[*itr].real() << "," << m_TzyObserved[*itr].imag() << "," << 
		//	m_TzxSD[*itr].realPart << "," << m_TzxSD[*itr].imagPart << "," << 
		//	m_TzySD[*itr].realPart << "," << m_TzySD[*itr].imagPart << ","<< std::endl;	

		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_stationID, m_freq[*itr], m_TzxCalculated[icount].real(), m_TzxCalculated[icount].imag(), m_TzyCalculated[icount].real(), m_TzyCalculated[icount].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_TzxResidual[icount].realPart, m_TzxResidual[icount].imagPart, m_TzyResidual[icount].realPart, m_TzyResidual[icount].imagPart ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_TzxObserved[*itr].real(), m_TzxObserved[*itr].imag(), m_TzyObserved[*itr].real(), m_TzyObserved[*itr].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_TzxSD[*itr].realPart, m_TzxSD[*itr].imagPart, m_TzySD[*itr].realPart, m_TzySD[*itr].imagPart ); 

		++icount;
	}
}

// Calulate interpolator vector of vertical magnetic field
void ObservedDataStationVTF::calcInterpolatorVectorOfVerticalMagneticField( Forward3D* const ptrForward3D ){

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldZDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHz );
	}else if( meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA ){
		ptrForward3D->calcInterpolatorVectorOfMagneticFieldZDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHz );
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Calulate sensitivity matrix of vertical magnetic field transfer functions
void ObservedDataStationVTF::calculateSensitivityMatrix( const double freq, const int nModel,
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

	const long long rhsVectorIDOfHx = static_cast<long long>( ptrStationOfMagneticField->getRhsVectorIDOfHx() );
	const long long rhsVectorIDOfHy = static_cast<long long>( ptrStationOfMagneticField->getRhsVectorIDOfHy() );
	const long long nBlkNotFixed = static_cast<long long>( ( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed() );

	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * HyCalculated[1]
										 + derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * HxCalculated[0]
										 - derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx + imdl ] * HyCalculated[0]
										 - derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy + imdl ] * HxCalculated[1];

		// dTzx/dm
		const std::complex<double> workXX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHz) + imdl ] *   HyCalculated[1]
										    + derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_HzCalculated[0]
										    - derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHz) + imdl ] *   HyCalculated[0]
										    - derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_HzCalculated[1];
		
		const std::complex<double> workXX2	= m_HzCalculated[0]*HyCalculated[1]	- m_HzCalculated[1]*HyCalculated[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTzx[freqIDThisPEInSta].realPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).real() / m_TzxSD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTzx[freqIDThisPEInSta].imagPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).imag() / m_TzxSD[freqIDGlobalInSta].imagPart;

		// dTzy/dm
		const std::complex<double> workXY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHz) + imdl ] *   HxCalculated[0]
											+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_HzCalculated[1]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHz) + imdl ] *   HxCalculated[1]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_HzCalculated[0];

		const std::complex<double> workXY2	= m_HzCalculated[1]*HxCalculated[0]	- m_HzCalculated[0]*HxCalculated[1];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTzy[freqIDThisPEInSta].realPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).real() / m_TzySD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTzy[freqIDThisPEInSta].imagPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).imag() / m_TzySD[freqIDGlobalInSta].imagPart;

	}

}

// Calculate data vector of this PE
void ObservedDataStationVTF::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	vector[ offset + m_dataIDOfTzx[freqIDThisPEInSta].realPart ] = m_TzxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTzx[freqIDThisPEInSta].imagPart ] = m_TzxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfTzy[freqIDThisPEInSta].realPart ] = m_TzyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTzy[freqIDThisPEInSta].imagPart ] = m_TzyResidual[freqIDThisPEInSta].imagPart;
}

// Calulate L2 norm of misfit
double ObservedDataStationVTF::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_TzxResidual[ifreq].realPart , 2 );
		misfit += pow( m_TzxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_TzyResidual[ifreq].realPart , 2 );
		misfit += pow( m_TzyResidual[ifreq].imagPart , 2 );
	}

	return misfit;

}

// Get VTK
bool ObservedDataStationVTF::getVTF( const double freq, std::complex<double>& Tzx, std::complex<double>& Tzy ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return false;
	}

	Tzx = m_TzxCalculated[freqIDThisPEInSta];
	Tzy = m_TzyCalculated[freqIDThisPEInSta];

	return true;

}

