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

#include "ObservedDataStationHTF.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationHTF::ObservedDataStationHTF():
	ObservedDataStationPoint(),
	m_TxxObserved(NULL),
	m_TxyObserved(NULL),
	m_TyxObserved(NULL),
	m_TyyObserved(NULL),
	m_TxxSD(NULL),
	m_TxySD(NULL),
	m_TyxSD(NULL),
	m_TyySD(NULL),
	m_TxxCalculated(NULL),
	m_TxyCalculated(NULL),
	m_TyxCalculated(NULL),
	m_TyyCalculated(NULL),
	m_TxxResidual(NULL),
	m_TxyResidual(NULL),
	m_TyxResidual(NULL),
	m_TyyResidual(NULL),
	m_dataIDOfTxx(NULL),
	m_dataIDOfTxy(NULL),
	m_dataIDOfTyx(NULL),
	m_dataIDOfTyy(NULL)
{
}

// Destructer
ObservedDataStationHTF::~ObservedDataStationHTF(){

	if( m_TxxObserved != NULL){
		delete[] m_TxxObserved;
		m_TxxObserved = NULL;
	}

	if( m_TxyObserved != NULL){
		delete[] m_TxyObserved;
		m_TxyObserved = NULL;
	}

	if( m_TyxObserved != NULL){
		delete[] m_TyxObserved;
		m_TyxObserved = NULL;
	}

	if( m_TyyObserved != NULL){
		delete[] m_TyyObserved;
		m_TyyObserved = NULL;
	}

	if( m_TxxSD != NULL){
		delete[] m_TxxSD;
		m_TxxSD = NULL;
	}

	if( m_TxySD != NULL){
		delete[] m_TxySD;
		m_TxySD = NULL;
	}

	if( m_TyxSD != NULL){
		delete[] m_TyxSD;
		m_TyxSD = NULL;
	}

	if( m_TyySD != NULL){
		delete[] m_TyySD;
		m_TyySD = NULL;
	}

	if( m_TxxCalculated != NULL){
		delete[] m_TxxCalculated;
		m_TxxCalculated = NULL;
	}

	if( m_TxyCalculated != NULL){
		delete[] m_TxyCalculated;
		m_TxyCalculated = NULL;
	}

	if( m_TyxCalculated != NULL){
		delete[] m_TyxCalculated;
		m_TyxCalculated = NULL;
	}

	if( m_TyyCalculated != NULL){
		delete[] m_TyyCalculated;
		m_TyyCalculated = NULL;
	}
	
	if( m_TxxResidual != NULL){
		delete[] m_TxxResidual;
		m_TxxResidual = NULL;
	}
	
	if( m_TxyResidual != NULL){
		delete[] m_TxyResidual;
		m_TxyResidual = NULL;
	}
	
	if( m_TyxResidual != NULL){
		delete[] m_TyxResidual;
		m_TyxResidual = NULL;
	}
	
	if( m_TyyResidual != NULL){
		delete[] m_TyyResidual;
		m_TyyResidual = NULL;
	}

	if( m_dataIDOfTxx != NULL){
		delete[] m_dataIDOfTxx;
		m_dataIDOfTxx = NULL;
	}

	if( m_dataIDOfTxy != NULL){
		delete[] m_dataIDOfTxy;
		m_dataIDOfTxy = NULL;
	}

	if( m_dataIDOfTyx != NULL){
		delete[] m_dataIDOfTyx;
		m_dataIDOfTyx = NULL;
	}

	if( m_dataIDOfTyy != NULL){
		delete[] m_dataIDOfTyy;
		m_dataIDOfTyy = NULL;
	}

}

// Read data from input file
void ObservedDataStationHTF::inputObservedData( std::ifstream& inFile ){

	inFile >> m_stationID;

	inFile >> m_IDOfMagneticFieldStation;

	if( m_stationID == m_IDOfMagneticFieldStation ){
		OutputFiles::m_logFile << "HTF station (" << m_stationID << ") must be different from its reference station (" << m_IDOfMagneticFieldStation << ") !!" << std::endl;
		exit(1);
	}
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
		m_TxxObserved = new std::complex<double>[nFreq];
		m_TxyObserved = new std::complex<double>[nFreq];
		m_TyxObserved = new std::complex<double>[nFreq];
		m_TyyObserved = new std::complex<double>[nFreq];
		m_TxxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_TxySD = new CommonParameters::DoubleComplexValues[nFreq];
		m_TyxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_TyySD = new CommonParameters::DoubleComplexValues[nFreq];
		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1 >> dbuf2;
			m_TxxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_TxyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_TyxObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> dbuf1 >> dbuf2;
			m_TyyObserved[i] = std::complex<double>(dbuf1,dbuf2);
			inFile >> m_TxxSD[i].realPart;
			inFile >> m_TxxSD[i].imagPart;
			inFile >> m_TxySD[i].realPart;
			inFile >> m_TxySD[i].imagPart;
			inFile >> m_TyxSD[i].realPart;
			inFile >> m_TyxSD[i].imagPart;
			inFile >> m_TyySD[i].realPart;
			inFile >> m_TyySD[i].imagPart;
		}
	}

#ifdef _DEBUG_WRITE
	std::cout <<  " HTF " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.X << " " << m_location.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for( int i = 0; i < m_numOfFrequency; ++i ){
		std::cout << m_freq[i] << " "
				    << m_TxxObserved[i] << " "
				    << m_TxyObserved[i] << " "
				    << m_TyxObserved[i] << " "
				    << m_TyyObserved[i] << " "
					<< m_TxxSD[i].realPart << " "
					<< m_TxxSD[i].imagPart << " "
					<< m_TxySD[i].realPart << " "
					<< m_TxySD[i].imagPart << " "
					<< m_TyxSD[i].realPart << " "
					<< m_TyxSD[i].imagPart << " "
					<< m_TyySD[i].realPart << " "
					<< m_TyySD[i].imagPart << std::endl;
	}
#endif

}

// Calulate horizontal magnetic field transfer function
void ObservedDataStationHTF::calculateHTF( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalRef[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalRef[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalRef[0]*HyCalRef[1] - HxCalRef[1]*HyCalRef[0];

	m_TxxCalculated[freqIDThisPEInSta] = ( m_HxCalculated[0]*HyCalRef[1] - m_HxCalculated[1]*HyCalRef[0] ) / det;
	m_TxyCalculated[freqIDThisPEInSta] = ( m_HxCalculated[1]*HxCalRef[0] - m_HxCalculated[0]*HxCalRef[1] ) / det;
	m_TyxCalculated[freqIDThisPEInSta] = ( m_HyCalculated[0]*HyCalRef[1] - m_HyCalculated[1]*HyCalRef[0] ) / det;
	m_TyyCalculated[freqIDThisPEInSta] = ( m_HyCalculated[1]*HxCalRef[0] - m_HyCalculated[0]*HxCalRef[1] ) / det;

	m_TxxResidual[freqIDThisPEInSta].realPart = ( m_TxxObserved[freqIDGlobalInSta].real() - m_TxxCalculated[freqIDThisPEInSta].real() ) / m_TxxSD[freqIDGlobalInSta].realPart;
	m_TxxResidual[freqIDThisPEInSta].imagPart = ( m_TxxObserved[freqIDGlobalInSta].imag() - m_TxxCalculated[freqIDThisPEInSta].imag() ) / m_TxxSD[freqIDGlobalInSta].imagPart;
	m_TxyResidual[freqIDThisPEInSta].realPart = ( m_TxyObserved[freqIDGlobalInSta].real() - m_TxyCalculated[freqIDThisPEInSta].real() ) / m_TxySD[freqIDGlobalInSta].realPart;
	m_TxyResidual[freqIDThisPEInSta].imagPart = ( m_TxyObserved[freqIDGlobalInSta].imag() - m_TxyCalculated[freqIDThisPEInSta].imag() ) / m_TxySD[freqIDGlobalInSta].imagPart;
	m_TyxResidual[freqIDThisPEInSta].realPart = ( m_TyxObserved[freqIDGlobalInSta].real() - m_TyxCalculated[freqIDThisPEInSta].real() ) / m_TyxSD[freqIDGlobalInSta].realPart;
	m_TyxResidual[freqIDThisPEInSta].imagPart = ( m_TyxObserved[freqIDGlobalInSta].imag() - m_TyxCalculated[freqIDThisPEInSta].imag() ) / m_TyxSD[freqIDGlobalInSta].imagPart;
	m_TyyResidual[freqIDThisPEInSta].realPart = ( m_TyyObserved[freqIDGlobalInSta].real() - m_TyyCalculated[freqIDThisPEInSta].real() ) / m_TyySD[freqIDGlobalInSta].realPart;
	m_TyyResidual[freqIDThisPEInSta].imagPart = ( m_TyyObserved[freqIDGlobalInSta].imag() - m_TyyCalculated[freqIDThisPEInSta].imag() ) / m_TyySD[freqIDGlobalInSta].imagPart;

#ifdef _DEBUG_WRITE
	std::cout << "freqIDThisPEInSta Txx Txy Tyx Tyy : " << freqIDThisPEInSta << " " << m_TxxCalculated[freqIDThisPEInSta] << " " << m_TxyCalculated[freqIDThisPEInSta] << " " << m_TyxCalculated[freqIDThisPEInSta] << " " << m_TyyCalculated[freqIDThisPEInSta] << std::endl;
#endif

	// For inversion
	m_dataIDOfTxx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTxx[freqIDThisPEInSta].imagPart = icount++;
	m_dataIDOfTxy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTxy[freqIDThisPEInSta].imagPart = icount++;
	m_dataIDOfTyx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTyx[freqIDThisPEInSta].imagPart = icount++;
	m_dataIDOfTyy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfTyy[freqIDThisPEInSta].imagPart = icount++;

}

// Initialize horizontal magnetic field transfer functions and errors
void ObservedDataStationHTF::initializeHTFsAndErrors(){

	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_TxxResidual[i].realPart = 0.0;
		m_TxxResidual[i].imagPart = 0.0;
		m_TxyResidual[i].realPart = 0.0;
		m_TxyResidual[i].imagPart = 0.0;
		m_TyxResidual[i].realPart = 0.0;
		m_TyxResidual[i].imagPart = 0.0;
		m_TyyResidual[i].realPart = 0.0;
		m_TyyResidual[i].imagPart = 0.0;

		m_TxxCalculated[i] = std::complex<double>(0.0,0.0);
		m_TxyCalculated[i] = std::complex<double>(0.0,0.0);
		m_TyxCalculated[i] = std::complex<double>(0.0,0.0);
		m_TyyCalculated[i] = std::complex<double>(0.0,0.0);
	}

}

// Allocate memory for the calculated values of horizontal magnetic field transfer functions and errors
void ObservedDataStationHTF::allocateMemoryForCalculatedValues(){

	if( m_TxxCalculated != NULL){
		delete[] m_TxxCalculated;
		m_TxxCalculated = NULL;
	}

	if( m_TxyCalculated != NULL){
		delete[] m_TxyCalculated;
		m_TxyCalculated = NULL;
	}

	if( m_TyxCalculated != NULL){
		delete[] m_TyxCalculated;
		m_TyxCalculated = NULL;
	}

	if( m_TyyCalculated != NULL){
		delete[] m_TyyCalculated;
		m_TyyCalculated = NULL;
	}
	
	if( m_TxxResidual != NULL){
		delete[] m_TxxResidual;
		m_TxxResidual = NULL;
	}
	
	if( m_TxyResidual != NULL){
		delete[] m_TxyResidual;
		m_TxyResidual = NULL;
	}
	
	if( m_TyxResidual != NULL){
		delete[] m_TyxResidual;
		m_TyxResidual = NULL;
	}
	
	if( m_TyyResidual != NULL){
		delete[] m_TyyResidual;
		m_TyyResidual = NULL;
	}

	if( m_dataIDOfTxx != NULL){
		delete[] m_dataIDOfTxx;
		m_dataIDOfTxx = NULL;
	}

	if( m_dataIDOfTxy != NULL){
		delete[] m_dataIDOfTxy;
		m_dataIDOfTxy = NULL;
	}

	if( m_dataIDOfTyx != NULL){
		delete[] m_dataIDOfTyx;
		m_dataIDOfTyx = NULL;
	}

	if( m_dataIDOfTyy != NULL){
		delete[] m_dataIDOfTyy;
		m_dataIDOfTyy = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){

		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_TxxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TxyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TyxCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TyyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_TxxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_TxyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_TyxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_TyyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTxx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTxy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTyx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfTyy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_TxxCalculated[i] = std::complex<double>(0.0,0.0);
			m_TxyCalculated[i] = std::complex<double>(0.0,0.0);
			m_TyxCalculated[i] = std::complex<double>(0.0,0.0);
			m_TyyCalculated[i] = std::complex<double>(0.0,0.0);
			m_TxxResidual[i].realPart = 0.0;
			m_TxxResidual[i].imagPart = 0.0;
			m_TxyResidual[i].realPart = 0.0;
			m_TxyResidual[i].imagPart = 0.0;
			m_TyxResidual[i].realPart = 0.0;
			m_TyxResidual[i].imagPart = 0.0;
			m_TyyResidual[i].realPart = 0.0;
			m_TyyResidual[i].imagPart = 0.0;
			m_dataIDOfTxx[i].realPart = -1;
			m_dataIDOfTxx[i].imagPart = -1;
			m_dataIDOfTxy[i].realPart = -1;
			m_dataIDOfTxy[i].imagPart = -1;
			m_dataIDOfTyx[i].realPart = -1;
			m_dataIDOfTyx[i].imagPart = -1;
			m_dataIDOfTyy[i].realPart = -1;
			m_dataIDOfTyy[i].imagPart = -1;
		}

	}

}

// Output calculated values of horizontal magnetic field transfer functions
void ObservedDataStationHTF::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){

		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_stationID, m_freq[*itr], 
			m_TxxCalculated[icount].real(), m_TxxCalculated[icount].imag(), m_TxyCalculated[icount].real(), m_TxyCalculated[icount].imag(),
			m_TyxCalculated[icount].real(), m_TyxCalculated[icount].imag(), m_TyyCalculated[icount].real(), m_TyyCalculated[icount].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_TxxResidual[icount].realPart, m_TxxResidual[icount].imagPart, m_TxyResidual[icount].realPart, m_TxyResidual[icount].imagPart,
			m_TyxResidual[icount].realPart, m_TyxResidual[icount].imagPart, m_TyyResidual[icount].realPart, m_TyyResidual[icount].imagPart ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_TxxObserved[*itr].real(), m_TxxObserved[*itr].imag(), m_TxyObserved[*itr].real(), m_TxyObserved[*itr].imag(),
			m_TyxObserved[*itr].real(), m_TyxObserved[*itr].imag(), m_TyyObserved[*itr].real(), m_TyyObserved[*itr].imag() ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_TxxSD[*itr].realPart, m_TxxSD[*itr].imagPart, m_TxySD[*itr].realPart, m_TxySD[*itr].imagPart,
			m_TyxSD[*itr].realPart, m_TyxSD[*itr].imagPart, m_TyySD[*itr].realPart, m_TyySD[*itr].imagPart ); 

		++icount;
	}
}

// Calulate sensitivity matrix of horizontal magnetic field transfer functions
void ObservedDataStationHTF::calculateSensitivityMatrix( const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMFieldExPol,
	const std::complex<double>* const derivativesOfEMFieldEyPol,
	double* const sensitivityMatrix ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalRef[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalRef[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalRef[0]*HyCalRef[1] - HxCalRef[1]*HyCalRef[0];
	const std::complex<double> divDet = std::complex<double>(1.0,0.0) / det;
	const std::complex<double> divDet2 = divDet * divDet;

	const long long rhsVectorIDOfHxRef = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	const long long rhsVectorIDOfHyRef = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());
	const long long nBlkNotFixed = static_cast<long long>( ( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed() );

	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){
		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] *   HyCalRef[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] *   HxCalRef[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] *   HyCalRef[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] *   HxCalRef[1];

		// dTxx/dm
		const std::complex<double> workXX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHx) + imdl ] *   HyCalRef[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] * m_HxCalculated[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHx) + imdl ] *   HyCalRef[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] * m_HxCalculated[1];
		
		const std::complex<double> workXX2	= m_HxCalculated[0]*HyCalRef[1]	- m_HxCalculated[1]*HyCalRef[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTxx[freqIDThisPEInSta].realPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).real() / m_TxxSD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTxx[freqIDThisPEInSta].imagPart) + imdl ] = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).imag() / m_TxxSD[freqIDGlobalInSta].imagPart;

		// dTxy/dm
		const std::complex<double> workXY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHx) + imdl ] *   HxCalRef[0]
										+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] * m_HxCalculated[1]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHx) + imdl ] *   HxCalRef[1]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] * m_HxCalculated[0];

		const std::complex<double> workXY2	= m_HxCalculated[1]*HxCalRef[0]	- m_HxCalculated[0]*HxCalRef[1];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTxy[freqIDThisPEInSta].realPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).real() / m_TxySD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTxy[freqIDThisPEInSta].imagPart) + imdl ] = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).imag() / m_TxySD[freqIDGlobalInSta].imagPart;

		// dTyx/dm
		const std::complex<double> workYX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHy) + imdl ] *   HyCalRef[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] * m_HyCalculated[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHy) + imdl ] *   HyCalRef[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHyRef + imdl ] * m_HyCalculated[1];

		const std::complex<double> workYX2	= m_HyCalculated[0]*HyCalRef[1]	- m_HyCalculated[1]*HyCalRef[0];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTyx[freqIDThisPEInSta].realPart) + imdl ] = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).real() / m_TyxSD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTyx[freqIDThisPEInSta].imagPart) + imdl ] = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).imag() / m_TyxSD[freqIDGlobalInSta].imagPart;

		// dTyy/dm
		const std::complex<double> workYY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHy) + imdl ] *   HxCalRef[0]
										+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] * m_HyCalculated[1]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfHy) + imdl ] *   HxCalRef[1]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHxRef + imdl ] * m_HyCalculated[0];

		const std::complex<double> workYY2	= m_HyCalculated[1]*HxCalRef[0]	- m_HyCalculated[0]*HxCalRef[1];

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTyy[freqIDThisPEInSta].realPart) + imdl ] = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).real() / m_TyySD[freqIDGlobalInSta].realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfTyy[freqIDThisPEInSta].imagPart) + imdl ] = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).imag() / m_TyySD[freqIDGlobalInSta].imagPart;
	}

}

// Calculate data vector of this PE
void ObservedDataStationHTF::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	vector[ offset + m_dataIDOfTxx[freqIDThisPEInSta].realPart ] = m_TxxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTxx[freqIDThisPEInSta].imagPart ] = m_TxxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfTxy[freqIDThisPEInSta].realPart ] = m_TxyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTxy[freqIDThisPEInSta].imagPart ] = m_TxyResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfTyx[freqIDThisPEInSta].realPart ] = m_TyxResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTyx[freqIDThisPEInSta].imagPart ] = m_TyxResidual[freqIDThisPEInSta].imagPart;
	vector[ offset + m_dataIDOfTyy[freqIDThisPEInSta].realPart ] = m_TyyResidual[freqIDThisPEInSta].realPart;
	vector[ offset + m_dataIDOfTyy[freqIDThisPEInSta].imagPart ] = m_TyyResidual[freqIDThisPEInSta].imagPart;
}

// Calulate L2 norm of misfit
double ObservedDataStationHTF::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_TxxResidual[ifreq].realPart , 2 );
		misfit += pow( m_TxxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_TxyResidual[ifreq].realPart , 2 );
		misfit += pow( m_TxyResidual[ifreq].imagPart , 2 );
		misfit += pow( m_TyxResidual[ifreq].realPart , 2 );
		misfit += pow( m_TyxResidual[ifreq].imagPart , 2 );
		misfit += pow( m_TyyResidual[ifreq].realPart , 2 );
		misfit += pow( m_TyyResidual[ifreq].imagPart , 2 );
	}

	return misfit;

}
