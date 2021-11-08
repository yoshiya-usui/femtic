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

#include "ObservedDataStationPT.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationPT::ObservedDataStationPT():
	ObservedDataStationPoint(),
	m_PTxxObserved(NULL),
	m_PTxyObserved(NULL),
	m_PTyxObserved(NULL),
	m_PTyyObserved(NULL),
	m_PTxxSD(NULL),
	m_PTxySD(NULL),
	m_PTyxSD(NULL),
	m_PTyySD(NULL),
	m_PTxxCalculated(NULL),
	m_PTxyCalculated(NULL),
	m_PTyxCalculated(NULL),
	m_PTyyCalculated(NULL),
	m_PTxxResidual(NULL),
	m_PTxyResidual(NULL),
	m_PTyxResidual(NULL),
	m_PTyyResidual(NULL),
	m_rhsVectorIDOfEx(0),
	m_rhsVectorIDOfEy(0),
	m_dataIDOfPTxx(NULL),
	m_dataIDOfPTxy(NULL),
	m_dataIDOfPTyx(NULL),
	m_dataIDOfPTyy(NULL),
	m_typeOfElectricField(AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD)
{
	for( int i = 0; i < 2; ++i ){
		m_ExCalculated[i] = std::complex<double>(0.0,0.0);
		m_EyCalculated[i] = std::complex<double>(0.0,0.0);
	}
}

// Destructer
ObservedDataStationPT::~ObservedDataStationPT(){

	if( m_PTxxObserved != NULL){
		delete[] m_PTxxObserved;
		m_PTxxObserved = NULL;
	}

	if( m_PTxyObserved != NULL){
		delete[] m_PTxyObserved;
		m_PTxyObserved = NULL;
	}

	if( m_PTyxObserved != NULL){
		delete[] m_PTyxObserved;
		m_PTyxObserved = NULL;
	}

	if( m_PTyyObserved != NULL){
		delete[] m_PTyyObserved;
		m_PTyyObserved = NULL;
	}

	if( m_PTxxSD != NULL){
		delete[] m_PTxxSD;
		m_PTxxSD = NULL;
	}

	if( m_PTxySD != NULL){
		delete[] m_PTxySD;
		m_PTxySD = NULL;
	}

	if( m_PTyxSD != NULL){
		delete[] m_PTyxSD;
		m_PTyxSD = NULL;
	}

	if( m_PTyySD != NULL){
		delete[] m_PTyySD;
		m_PTyySD = NULL;
	}

	if( m_PTxxCalculated != NULL){
		delete[] m_PTxxCalculated;
		m_PTxxCalculated = NULL;
	}

	if( m_PTxyCalculated != NULL){
		delete[] m_PTxyCalculated;
		m_PTxyCalculated = NULL;
	}

	if( m_PTyxCalculated != NULL){
		delete[] m_PTyxCalculated;
		m_PTyxCalculated = NULL;
	}

	if( m_PTyyCalculated != NULL){
		delete[] m_PTyyCalculated;
		m_PTyyCalculated = NULL;
	}

	if( m_PTxxResidual != NULL){
		delete[] m_PTxxResidual;
		m_PTxxResidual = NULL;
	}

	if( m_PTxyResidual != NULL){
		delete[] m_PTxyResidual;
		m_PTxyResidual = NULL;
	}

	if( m_PTyxResidual != NULL){
		delete[] m_PTyxResidual;
		m_PTyxResidual = NULL;
	}

	if( m_PTyyResidual != NULL){
		delete[] m_PTyyResidual;
		m_PTyyResidual = NULL;
	}

	if( m_dataIDOfPTxx != NULL){
		delete[] m_dataIDOfPTxx;
		m_dataIDOfPTxx = NULL;
	}

	if( m_dataIDOfPTxy != NULL){
		delete[] m_dataIDOfPTxy;
		m_dataIDOfPTxy = NULL;
	}

	if( m_dataIDOfPTyx != NULL){
		delete[] m_dataIDOfPTyx;
		m_dataIDOfPTyx = NULL;
	}

	if( m_dataIDOfPTyy != NULL){
		delete[] m_dataIDOfPTyy;
		m_dataIDOfPTyy = NULL;
	}

}

// Read data from input file
void ObservedDataStationPT::inputObservedData( std::ifstream& inFile ){

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
	OutputFiles::m_logFile << std::setw(15) << std::left << ownerElemType;

	// Type of electric field
	if( pAnalysisControl->isTypeOfElectricFieldSetIndivisually() ){
		// Electric field type of each site is specified
		inFile >> m_typeOfElectricField;
	}else{
		m_typeOfElectricField = pAnalysisControl->getTypeOfElectricField();
	}

	if( pAnalysisControl->getTypeOfMesh() == MeshData::HEXA ){
		if(	m_typeOfElectricField != AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD ){
			OutputFiles::m_logFile << std::endl << "Warning : Horizontal electric field must be used for hexahedral mesh." << std::endl;
		}
		m_typeOfElectricField = AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD;
	}
	std::string elecType;
	switch( m_typeOfElectricField ){
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			elecType = "Tangential";
			break;
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			elecType = "Horizontal";
			break;
		default:
			OutputFiles::m_logFile << std::endl << "Error : Unknown type of the electric field : " << m_typeOfElectricField << std::endl;
			exit(1);
			break;
	}
	OutputFiles::m_logFile << std::left << elecType << std::endl;

	double dbuf(0.0);
	inFile >> dbuf;
	m_location.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.Y = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> m_numOfFrequency;
	const int nFreq = m_numOfFrequency;
	if( nFreq > 0 ){
		m_freq = new double[nFreq];
		m_PTxxObserved = new double[nFreq];
		m_PTxyObserved = new double[nFreq];
		m_PTyxObserved = new double[nFreq];
		m_PTyyObserved = new double[nFreq];
		m_PTxxSD = new double[nFreq];
		m_PTxySD = new double[nFreq];
		m_PTyxSD = new double[nFreq];
		m_PTyySD = new double[nFreq];
		//m_PTxxCalculated = new double[nFreq];
		//m_PTxyCalculated = new double[nFreq];
		//m_PTyxCalculated = new double[nFreq];
		//m_PTyyCalculated = new double[nFreq];
		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			inFile >> m_PTxxObserved[i];
			inFile >> m_PTxyObserved[i];
			inFile >> m_PTyxObserved[i];
			inFile >> m_PTyyObserved[i];
			inFile >> m_PTxxSD[i];
			inFile >> m_PTxySD[i];
			inFile >> m_PTyxSD[i];
			inFile >> m_PTyySD[i];
		}
	}

#ifdef _DEBUG_WRITE
	std::cout <<  " PT " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.X << " " << m_location.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for( int i = 0; i < m_numOfFrequency; ++i ){
		std::cout << m_freq[i] << " "
				    << m_PTxxObserved[i] << " "
				    << m_PTxyObserved[i] << " "
				    << m_PTyxObserved[i] << " "
				    << m_PTyyObserved[i] << " "
					<< m_PTxxSD[i] << " "
					<< m_PTxySD[i] << " "
					<< m_PTyxSD[i] << " "
					<< m_PTyySD[i] << std::endl;
	}
#endif

}

// Calulate electric field
void ObservedDataStationPT::calculateElectricField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy ){

	const int iPol = ptrForward3D->getPolarizationCurrent();

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		switch( getTypeOfElectricField() ){
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
				m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3 );
				break;
			case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
				m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldTangentialX( m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2 );
				m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldTangentialY( m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2 );
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
				exit(1);
				break;
		}
	}else if( meshType == MeshData::HEXA ){
		m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else if( meshType == MeshData::NONCONFORMING_HEXA ){
		switch( getTypeOfElectricField() ){
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
				m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
				break;
			case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
				m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldTangentialXFromAllEdges( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
				m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldTangentialYFromAllEdges( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
				exit(1);
				break;
		}
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

	// For inversion
	m_rhsVectorIDOfEx = rhsVectorIDOfEx;
	m_rhsVectorIDOfEy = rhsVectorIDOfEy;

}

// Calulate Impedance tensor
void ObservedDataStationPT::calculatePhaseTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> d = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];

	std::complex<double> Zxx = ( m_ExCalculated[0]*HyCalculated[1] - m_ExCalculated[1]*HyCalculated[0] ) / d;
	std::complex<double> Zxy = ( m_ExCalculated[1]*HxCalculated[0] - m_ExCalculated[0]*HxCalculated[1] ) / d;
	std::complex<double> Zyx = ( m_EyCalculated[0]*HyCalculated[1] - m_EyCalculated[1]*HyCalculated[0] ) / d;
	std::complex<double> Zyy = ( m_EyCalculated[1]*HxCalculated[0] - m_EyCalculated[0]*HxCalculated[1] ) / d;

	const double det = Zxx.real()*Zyy.real() - Zxy.real()*Zyx.real();

	m_PTxxCalculated[freqIDThisPEInSta] = ( Zyy.real()*Zxx.imag() - Zxy.real()*Zyx.imag() )/ det;
	m_PTxyCalculated[freqIDThisPEInSta] = ( Zyy.real()*Zxy.imag() - Zxy.real()*Zyy.imag() )/ det;
	m_PTyxCalculated[freqIDThisPEInSta] = ( Zxx.real()*Zyx.imag() - Zyx.real()*Zxx.imag() )/ det;
	m_PTyyCalculated[freqIDThisPEInSta] = ( Zxx.real()*Zyy.imag() - Zyx.real()*Zxy.imag() )/ det;

	m_PTxxResidual[freqIDThisPEInSta] = ( m_PTxxObserved[freqIDGlobalInSta] - m_PTxxCalculated[freqIDThisPEInSta] ) / m_PTxxSD[freqIDGlobalInSta];
	m_PTxyResidual[freqIDThisPEInSta] = ( m_PTxyObserved[freqIDGlobalInSta] - m_PTxyCalculated[freqIDThisPEInSta] ) / m_PTxySD[freqIDGlobalInSta];
	m_PTyxResidual[freqIDThisPEInSta] = ( m_PTyxObserved[freqIDGlobalInSta] - m_PTyxCalculated[freqIDThisPEInSta] ) / m_PTyxSD[freqIDGlobalInSta];
	m_PTyyResidual[freqIDThisPEInSta] = ( m_PTyyObserved[freqIDGlobalInSta] - m_PTyyCalculated[freqIDThisPEInSta] ) / m_PTyySD[freqIDGlobalInSta];

#ifdef _DEBUG_WRITE
	std::cout << "ifreq Pxx Pxy Pyx Pyy: " << freqIDThisPEInSta << " " << m_PTxxCalculated[freqIDThisPEInSta] << " " << m_PTxyCalculated[freqIDThisPEInSta]<< " " << m_PTyxCalculated[freqIDThisPEInSta] << " " << m_PTyyCalculated[freqIDThisPEInSta] << std::endl;
	std::cout << "obs cal sd res: " << m_PTxxObserved[freqIDGlobalInSta] << " " << m_PTxxCalculated[freqIDThisPEInSta] << " " << m_PTxxSD[freqIDGlobalInSta] << " " << m_PTxxResidual[freqIDThisPEInSta] << std::endl;
#endif

	m_dataIDOfPTxx[freqIDThisPEInSta] = icount++;
	m_dataIDOfPTxy[freqIDThisPEInSta] = icount++;
	m_dataIDOfPTyx[freqIDThisPEInSta] = icount++;
	m_dataIDOfPTyy[freqIDThisPEInSta] = icount++;
}

// Initialize electric field
void ObservedDataStationPT::initializeElectricField( const int iPol ){

	m_ExCalculated[iPol] = std::complex<double>(0.0,0.0);
	m_EyCalculated[iPol] = std::complex<double>(0.0,0.0);

}

// Initialize Impedance tensor and errors
void ObservedDataStationPT::initializePhaseTensorsAndErrors(){

	//for( int i = 0; i < m_numOfFrequency; ++i ){
	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_PTxxResidual[i] = 0.0;
		m_PTxyResidual[i] = 0.0;
		m_PTyxResidual[i] = 0.0;
		m_PTyyResidual[i] = 0.0;

		m_PTxxCalculated[i] = 0.0;
		m_PTxyCalculated[i] = 0.0;
		m_PTyxCalculated[i] = 0.0;
		m_PTyyCalculated[i] = 0.0;
	}

}

// Allocate memory for the calculated values of phase tensors and errors
void ObservedDataStationPT::allocateMemoryForCalculatedValues(){

	if( m_PTxxCalculated != NULL){
		delete[] m_PTxxCalculated;
		m_PTxxCalculated = NULL;
	}

	if( m_PTxyCalculated != NULL){
		delete[] m_PTxyCalculated;
		m_PTxyCalculated = NULL;
	}

	if( m_PTyxCalculated != NULL){
		delete[] m_PTyxCalculated;
		m_PTyxCalculated = NULL;
	}

	if( m_PTyyCalculated != NULL){
		delete[] m_PTyyCalculated;
		m_PTyyCalculated = NULL;
	}

	if( m_PTxxResidual != NULL){
		delete[] m_PTxxResidual;
		m_PTxxResidual = NULL;
	}

	if( m_PTxyResidual != NULL){
		delete[] m_PTxyResidual;
		m_PTxyResidual = NULL;
	}

	if( m_PTyxResidual != NULL){
		delete[] m_PTyxResidual;
		m_PTyxResidual = NULL;
	}

	if( m_PTyyResidual != NULL){
		delete[] m_PTyyResidual;
		m_PTyyResidual = NULL;
	}
	
	if( m_dataIDOfPTxx != NULL){
		delete[] m_dataIDOfPTxx;
		m_dataIDOfPTxx = NULL;
	}

	if( m_dataIDOfPTxy != NULL){
		delete[] m_dataIDOfPTxy;
		m_dataIDOfPTxy = NULL;
	}

	if( m_dataIDOfPTyx != NULL){
		delete[] m_dataIDOfPTyx;
		m_dataIDOfPTyx = NULL;
	}

	if( m_dataIDOfPTyy != NULL){
		delete[] m_dataIDOfPTyy;
		m_dataIDOfPTyy = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){
		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_PTxxCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTxyCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTyxCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTyyCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTxxResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTxyResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTyxResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PTyyResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPTxx = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPTxy = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPTyx = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPTyy = new int[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_PTxxCalculated[i] = 0.0;
			m_PTxyCalculated[i] = 0.0;
			m_PTyxCalculated[i] = 0.0;
			m_PTyyCalculated[i] = 0.0;
			m_PTxxResidual[i] = 0.0;
			m_PTxyResidual[i] = 0.0;
			m_PTyxResidual[i] = 0.0;
			m_PTyyResidual[i] = 0.0;
			m_dataIDOfPTxx[i] = -1;
			m_dataIDOfPTxy[i] = -1;
			m_dataIDOfPTyx[i] = -1;
			m_dataIDOfPTyy[i] = -1;
		}
	}

}

// Output calculated values of phase tensors
void ObservedDataStationPT::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		//OutputFiles::m_csvFile << m_stationID << "," << m_freq[*itr] << "," << 
		//	m_PTxxCalculated[icount] << "," << m_PTxyCalculated[icount] << "," << 
		//	m_PTyxCalculated[icount] << "," << m_PTyyCalculated[icount] << "," << 
		//	m_PTxxResidual[icount] << "," << m_PTxyResidual[icount] << "," << 
		//	m_PTyxResidual[icount] << "," << m_PTyyResidual[icount] << "," << 
		//	m_PTxxObserved[*itr] << "," << m_PTxyObserved[*itr] << "," << 
		//	m_PTyxObserved[*itr] << "," << m_PTyyObserved[*itr] << "," << 
		//	m_PTxxSD[*itr] << "," << m_PTxySD[*itr] << "," << 
		//	m_PTyxSD[*itr] << "," << m_PTyySD[*itr] << std::endl;	

		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_stationID, m_freq[*itr], m_PTxxCalculated[icount], m_PTxyCalculated[icount], m_PTyxCalculated[icount], m_PTyyCalculated[icount] ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_PTxxResidual[icount], m_PTxyResidual[icount], m_PTyxResidual[icount], m_PTyyResidual[icount] ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,",
			m_PTxxObserved[*itr], m_PTxyObserved[*itr], m_PTyxObserved[*itr], m_PTyyObserved[*itr] ); 
		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_PTxxSD[*itr], m_PTxySD[*itr], m_PTyxSD[*itr], m_PTyySD[*itr] ); 

		++icount;
	}
}

// Calulate interpolator vector of electric field
void ObservedDataStationPT::calcInterpolatorVectorOfElectricField( Forward3D* const ptrForward3D ){

	const int meshType = ( AnalysisControl::getInstance() )->getTypeOfMesh();
	if( meshType == MeshData::TETRA ){// Tetra mesh
		switch( getTypeOfElectricField() ){
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfEx );
				ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection( m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfEy );
				break;
			case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialX( m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2, m_rhsVectorIDOfEx );
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialY( m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2, m_rhsVectorIDOfEy );
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
				exit(1);
				break;
		}
	}else if( meshType == MeshData::HEXA ){
		ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx );
		ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy );
	}else if( meshType == MeshData::NONCONFORMING_HEXA ){
		switch( getTypeOfElectricField() ){
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx );
				ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy );
				break;
			case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx );
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx );
				break;
			default:
				OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
				exit(1);
				break;
		}
	}else{
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Calulate sensitivity matrix of phase tensors
void ObservedDataStationPT::calculateSensitivityMatrix( const double freq, const int nModel,
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

	const std::complex<double> d = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];
	const std::complex<double> Zxx = ( m_ExCalculated[0]*HyCalculated[1] - m_ExCalculated[1]*HyCalculated[0] ) / d;
	const std::complex<double> Zxy = ( m_ExCalculated[1]*HxCalculated[0] - m_ExCalculated[0]*HxCalculated[1] ) / d;
	const std::complex<double> Zyx = ( m_EyCalculated[0]*HyCalculated[1] - m_EyCalculated[1]*HyCalculated[0] ) / d;
	const std::complex<double> Zyy = ( m_EyCalculated[1]*HxCalculated[0] - m_EyCalculated[0]*HxCalculated[1] ) / d;

#ifdef _DEBUG_WRITE
		std::cout << "Zxx, Zxy, Zyx, Zyy, factor : " << Zxx << " " << Zxy << " " << Zyx << " " << Zyy << " " << ( Zxx.real()*Zyy.real() - Zxy.real()*Zyx.real() ) << std::endl;
#endif

	const long long rhsVectorIDOfHx = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	const long long rhsVectorIDOfHy = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());
	const long long nBlkNotFixed = static_cast<long long>(( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed());

	const double factor = 1.0 / ( Zxx.real()*Zyy.real() - Zxy.real()*Zyx.real() );

	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] *   HyCalculated[1]
										 + derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] *   HxCalculated[0]
										 - derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] *   HyCalculated[0]
										 - derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] *   HxCalculated[1];

		// dZxx/dm
		const std::complex<double> workXX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] *   HyCalculated[1]
											+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_ExCalculated[0]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] *   HyCalculated[0]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_ExCalculated[1];
		
		const std::complex<double> workXX2	= m_ExCalculated[0]*HyCalculated[1]	- m_ExCalculated[1]*HyCalculated[0];

		const CommonParameters::DoubleComplexValues dZxxdm = {
			( workXX1 * divDet - work1 * workXX2 * divDet2 ).real(),
			( workXX1 * divDet - work1 * workXX2 * divDet2 ).imag()
		};

		// dZxy/dm
		const std::complex<double> workXY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] *   HxCalculated[0]
											+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_ExCalculated[1]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] *   HxCalculated[1]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_ExCalculated[0];

		const std::complex<double> workXY2	= m_ExCalculated[1]*HxCalculated[0]	- m_ExCalculated[0]*HxCalculated[1];

		const CommonParameters::DoubleComplexValues dZxydm = {
			( workXY1 * divDet - work1 * workXY2 * divDet2 ).real(),
			( workXY1 * divDet - work1 * workXY2 * divDet2 ).imag()
		};

		// dZyx/dm
		const std::complex<double> workYX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HyCalculated[1]
											+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_EyCalculated[0]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HyCalculated[0]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_EyCalculated[1];

		const std::complex<double> workYX2	= m_EyCalculated[0]*HyCalculated[1]	- m_EyCalculated[1]*HyCalculated[0];

		const CommonParameters::DoubleComplexValues dZyxdm = {
			( workYX1 * divDet - work1 * workYX2 * divDet2 ).real(),
			( workYX1 * divDet - work1 * workYX2 * divDet2 ).imag()
		};

		// dZyy/dm
		const std::complex<double> workYY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HxCalculated[0]
											+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_EyCalculated[1]
											- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HxCalculated[1]
											- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_EyCalculated[0];

		const std::complex<double> workYY2	= m_EyCalculated[1]*HxCalculated[0]	- m_EyCalculated[0]*HxCalculated[1];

		const CommonParameters::DoubleComplexValues dZyydm = {
			( workYY1 * divDet - work1 * workYY2 * divDet2 ).real(),
			( workYY1 * divDet - work1 * workYY2 * divDet2 ).imag()
		};

		const double factor2 = factor * factor * ( dZxxdm.realPart * Zyy.real() + Zxx.real() * dZyydm.realPart - dZxydm.realPart * Zyx.real() - Zxy.real() * dZyxdm.realPart );

		// dPTxx/dm
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTxx[freqIDThisPEInSta]) + imdl ] =
			  ( dZyydm.realPart * Zxx.imag() + Zyy.real() * dZxxdm.imagPart - dZxydm.realPart * Zyx.imag() - Zxy.real() * dZyxdm.imagPart ) * factor
			- ( Zyy.real() * Zxx.imag() - Zxy.real() * Zyx.imag() ) * factor2; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTxx[freqIDThisPEInSta]) + imdl ] /= m_PTxxSD[freqIDGlobalInSta];

		// dPTxy/dm
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTxy[freqIDThisPEInSta]) + imdl ] =
			  ( dZyydm.realPart * Zxy.imag() + Zyy.real() * dZxydm.imagPart - dZxydm.realPart * Zyy.imag() - Zxy.real() * dZyydm.imagPart ) * factor
			- ( Zyy.real() * Zxy.imag() - Zxy.real() * Zyy.imag() ) * factor2; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTxy[freqIDThisPEInSta]) + imdl ] /= m_PTxySD[freqIDGlobalInSta];

		// dPTyx/dm
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTyx[freqIDThisPEInSta]) + imdl ] =
			  ( dZxxdm.realPart * Zyx.imag() + Zxx.real() * dZyxdm.imagPart - dZyxdm.realPart * Zxx.imag() - Zyx.real() * dZxxdm.imagPart ) * factor
			- ( Zxx.real() * Zyx.imag() - Zyx.real() * Zxx.imag() ) * factor2; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTyx[freqIDThisPEInSta]) + imdl ] /= m_PTyxSD[freqIDGlobalInSta];

		// dPTyy/dm
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTyy[freqIDThisPEInSta]) + imdl ] =
			  ( dZxxdm.realPart * Zyy.imag() + Zxx.real() * dZyydm.imagPart - dZyxdm.realPart * Zxy.imag() - Zyx.real() * dZxydm.imagPart ) * factor
			- ( Zxx.real() * Zyy.imag() - Zyx.real() * Zxy.imag() ) * factor2; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPTyy[freqIDThisPEInSta]) + imdl ] /= m_PTyySD[freqIDGlobalInSta];

#ifdef _DEBUG_WRITE
		std::cout << "senseMat: "
			<< sensitivityMatrix[ nModel*m_dataIDOfPTxx[freqIDThisPEInSta] + imdl ] << " " 
			<< sensitivityMatrix[ nModel*m_dataIDOfPTxy[freqIDThisPEInSta] + imdl ] << " " 
			<< sensitivityMatrix[ nModel*m_dataIDOfPTyx[freqIDThisPEInSta] + imdl ] << " " 
			<< sensitivityMatrix[ nModel*m_dataIDOfPTyy[freqIDThisPEInSta] + imdl ] << std::endl;
#endif

	}

}

// Calculate data vector of this PE
void ObservedDataStationPT::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	vector[ offset + m_dataIDOfPTxx[freqIDThisPEInSta] ] = m_PTxxResidual[freqIDThisPEInSta];
	vector[ offset + m_dataIDOfPTxy[freqIDThisPEInSta] ] = m_PTxyResidual[freqIDThisPEInSta];
	vector[ offset + m_dataIDOfPTyx[freqIDThisPEInSta] ] = m_PTyxResidual[freqIDThisPEInSta];
	vector[ offset + m_dataIDOfPTyy[freqIDThisPEInSta] ] = m_PTyyResidual[freqIDThisPEInSta];
}

// Calulate error sum of squares
double ObservedDataStationPT::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_PTxxResidual[ifreq], 2 );
		misfit += pow( m_PTxyResidual[ifreq], 2 );
		misfit += pow( m_PTyxResidual[ifreq], 2 );
		misfit += pow( m_PTyyResidual[ifreq], 2 );
	}

	return misfit;

}

// Get type of the electric field used to calculate response functions
int ObservedDataStationPT::getTypeOfElectricField() const{

	return m_typeOfElectricField;

}

// Set type of the electric field used to calculate response functions
void ObservedDataStationPT::setTypeOfElectricField( const int type ){

	m_typeOfElectricField = type;

}