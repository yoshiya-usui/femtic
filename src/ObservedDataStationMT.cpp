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

#include "ObservedDataStationMT.h"
#include "ObservedDataStationPoint.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationMT::ObservedDataStationMT():
	ObservedDataStationPoint(),
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
	m_rhsVectorIDOfEx(0),
	m_rhsVectorIDOfEy(0),
	m_dataIDOfZxx(NULL),
	m_dataIDOfZxy(NULL),
	m_dataIDOfZyx(NULL),
	m_dataIDOfZyy(NULL),
	m_fixDistortionMatrix(true),
	m_typeOfElectricField(AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD),
	m_arrayDistortionMatrixDifferences(NULL),
	m_arrayGainsAndRotations(NULL)
{

	for( int i = 0; i < 2; ++i ){
		m_ExCalculated[i] = std::complex<double>(0.0,0.0);
		m_EyCalculated[i] = std::complex<double>(0.0,0.0);
	}

	if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::DISTORTION_TYPE_UNDEFINED ){
		OutputFiles::m_logFile << "Error : Type of distortion must be defined before instantiation of ObservedDataStationMT !!" << std::endl;
		exit(1);
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		m_arrayDistortionMatrixDifferences = new DistortionMatrixDifferences;
		for( int i = 0; i < 4; ++i ){
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[i] = 0.0;
			m_arrayDistortionMatrixDifferences->distortionMatrixDifference[i] = 0.0;
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[i] = 0.0;
			m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[i] = -1;
		}
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
			 ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		m_arrayGainsAndRotations = new GainsAndRotations;
		for( int i = 0; i < 4; ++i ){
			m_arrayGainsAndRotations->gainsAndRotationsPre[i] = 0.0;
			m_arrayGainsAndRotations->gainsAndRotations[i] = 0.0;
			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[i] = 0.0;
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[i] = -1;
		}
	}

}

// Destructer
ObservedDataStationMT::~ObservedDataStationMT(){

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

	if( m_arrayDistortionMatrixDifferences != NULL){
		delete m_arrayDistortionMatrixDifferences;
		m_arrayDistortionMatrixDifferences = NULL;
	}

	if( m_arrayGainsAndRotations != NULL){
		delete m_arrayGainsAndRotations;
		m_arrayGainsAndRotations = NULL;
	}

}

// Read data from input file
void ObservedDataStationMT::inputObservedData( std::ifstream& inFile ){

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
		m_ZxxObserved = new std::complex<double>[nFreq];
		m_ZxyObserved = new std::complex<double>[nFreq];
		m_ZyxObserved = new std::complex<double>[nFreq];
		m_ZyyObserved = new std::complex<double>[nFreq];
		m_ZxxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZxySD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyySD = new CommonParameters::DoubleComplexValues[nFreq];
		//m_ZxxCalculated = new std::complex<double>[nFreq];
		//m_ZxyCalculated = new std::complex<double>[nFreq];
		//m_ZyxCalculated = new std::complex<double>[nFreq];
		//m_ZyyCalculated = new std::complex<double>[nFreq];

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
	std::cout << " MT " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.X << " " << m_location.Y << std::endl;
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

// Calulate electric field
void ObservedDataStationMT::calculateElectricField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy ){

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
	}else if( meshType == MeshData::HEXA ){// Hexa mesh
		m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
		m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z );
	}else if( meshType == MeshData::NONCONFORMING_HEXA ){
		switch( getTypeOfElectricField() ){
			case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
				m_ExCalculated[iPol] = ptrForward3D->calcValueElectricFieldXDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z  );
				m_EyCalculated[iPol] = ptrForward3D->calcValueElectricFieldYDirection( m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z  );
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
void ObservedDataStationMT::calculateImpedanceTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated( 0 ), ptrStationOfMagneticField->getHxCalculated( 1 ) };
	const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated( 0 ), ptrStationOfMagneticField->getHyCalculated( 1 ) };

	const std::complex<double> det = HxCalculated[0]*HyCalculated[1] - HxCalculated[1]*HyCalculated[0];

	const std::complex<double> ZxxUndist = ( m_ExCalculated[0]*HyCalculated[1] - m_ExCalculated[1]*HyCalculated[0] ) / det;
	const std::complex<double> ZxyUndist = ( m_ExCalculated[1]*HxCalculated[0] - m_ExCalculated[0]*HxCalculated[1] ) / det;
	const std::complex<double> ZyxUndist = ( m_EyCalculated[0]*HyCalculated[1] - m_EyCalculated[1]*HyCalculated[0] ) / det;
	const std::complex<double> ZyyUndist = ( m_EyCalculated[1]*HxCalculated[0] - m_EyCalculated[0]*HxCalculated[1] ) / det;

	std::complex<double> cxx( 1.0, 0.0 );
	std::complex<double> cxy( 0.0, 0.0 );
	std::complex<double> cyx( 0.0, 0.0 );
	std::complex<double> cyy( 1.0, 0.0 );
	const int typeOfDistortion = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		cxx = std::complex<double>( m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXX] + 1.0 , 0.0 );
		cxy = std::complex<double>( m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXY]       , 0.0 );
		cyx = std::complex<double>( m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYX]       , 0.0 );
		cyy = std::complex<double>( m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYY] + 1.0 , 0.0 );
	}
	else if( typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
		const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );
		const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
		const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];

		cxx = std::complex<double>(   gX * cos( betaX ), 0.0 );
		cxy = std::complex<double>( - gY * sin( betaY ), 0.0 );
		cyx = std::complex<double>(   gX * sin( betaX ), 0.0 );
		cyy = std::complex<double>(   gY * cos( betaY ), 0.0 );
	}
	else if( typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
		const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );

		cxx = std::complex<double>(   gX, 0.0 );
		cxy = std::complex<double>(  0.0, 0.0 );
		cyx = std::complex<double>(  0.0, 0.0 );
		cyy = std::complex<double>(   gY, 0.0 );
	}

	m_ZxxCalculated[freqIDThisPEInSta] = cxx * ZxxUndist + cxy * ZyxUndist;
	m_ZxyCalculated[freqIDThisPEInSta] = cxx * ZxyUndist + cxy * ZyyUndist;
	m_ZyxCalculated[freqIDThisPEInSta] = cyx * ZxxUndist + cyy * ZyxUndist;
	m_ZyyCalculated[freqIDThisPEInSta] = cyx * ZxyUndist + cyy * ZyyUndist;

	m_ZxxResidual[freqIDThisPEInSta].realPart = ( m_ZxxObserved[freqIDGlobalInSta].real() - m_ZxxCalculated[freqIDThisPEInSta].real() ) / m_ZxxSD[freqIDGlobalInSta].realPart;
	m_ZxxResidual[freqIDThisPEInSta].imagPart = ( m_ZxxObserved[freqIDGlobalInSta].imag() - m_ZxxCalculated[freqIDThisPEInSta].imag() ) / m_ZxxSD[freqIDGlobalInSta].imagPart;
	m_ZxyResidual[freqIDThisPEInSta].realPart = ( m_ZxyObserved[freqIDGlobalInSta].real() - m_ZxyCalculated[freqIDThisPEInSta].real() ) / m_ZxySD[freqIDGlobalInSta].realPart;
	m_ZxyResidual[freqIDThisPEInSta].imagPart = ( m_ZxyObserved[freqIDGlobalInSta].imag() - m_ZxyCalculated[freqIDThisPEInSta].imag() ) / m_ZxySD[freqIDGlobalInSta].imagPart;
	m_ZyxResidual[freqIDThisPEInSta].realPart = ( m_ZyxObserved[freqIDGlobalInSta].real() - m_ZyxCalculated[freqIDThisPEInSta].real() ) / m_ZyxSD[freqIDGlobalInSta].realPart;
	m_ZyxResidual[freqIDThisPEInSta].imagPart = ( m_ZyxObserved[freqIDGlobalInSta].imag() - m_ZyxCalculated[freqIDThisPEInSta].imag() ) / m_ZyxSD[freqIDGlobalInSta].imagPart;
	m_ZyyResidual[freqIDThisPEInSta].realPart = ( m_ZyyObserved[freqIDGlobalInSta].real() - m_ZyyCalculated[freqIDThisPEInSta].real() ) / m_ZyySD[freqIDGlobalInSta].realPart;
	m_ZyyResidual[freqIDThisPEInSta].imagPart = ( m_ZyyObserved[freqIDGlobalInSta].imag() - m_ZyyCalculated[freqIDThisPEInSta].imag() ) / m_ZyySD[freqIDGlobalInSta].imagPart;

#ifdef _DEBUG_WRITE
	std::cout << "freqIDThisPEInSta Zxx Zxy Zyx Zyy : " << freqIDThisPEInSta << " " << m_ZxxCalculated[freqIDThisPEInSta] << " " << m_ZxyCalculated[freqIDThisPEInSta] << " " << m_ZyxCalculated[freqIDThisPEInSta] << " " << m_ZyyCalculated[freqIDThisPEInSta] << std::endl;
#endif

	// For inversion
	m_dataIDOfZxx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfZxx[freqIDThisPEInSta].imagPart = icount++;

	m_dataIDOfZxy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfZxy[freqIDThisPEInSta].imagPart = icount++;

	m_dataIDOfZyx[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfZyx[freqIDThisPEInSta].imagPart = icount++;

	m_dataIDOfZyy[freqIDThisPEInSta].realPart = icount++;
	m_dataIDOfZyy[freqIDThisPEInSta].imagPart = icount++;

#ifdef _DEBUG_WRITE
	std::cout << "m_dataIDOfZxx[freqIDThisPEInSta].realPart : " << m_dataIDOfZxx[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZxx[freqIDThisPEInSta].imagPart : " << m_dataIDOfZxx[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZxy[freqIDThisPEInSta].realPart : " << m_dataIDOfZxy[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZxy[freqIDThisPEInSta].imagPart : " << m_dataIDOfZxy[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZyx[freqIDThisPEInSta].realPart : " << m_dataIDOfZyx[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZyx[freqIDThisPEInSta].imagPart : " << m_dataIDOfZyx[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZyy[freqIDThisPEInSta].realPart : " << m_dataIDOfZyy[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZyy[freqIDThisPEInSta].imagPart : " << m_dataIDOfZyy[freqIDThisPEInSta].imagPart << std::endl;
#endif

}

// Initialize electric field
void ObservedDataStationMT::initializeElectricField( const int iPol ){

	m_ExCalculated[iPol] = std::complex<double>(0.0,0.0);
	m_EyCalculated[iPol] = std::complex<double>(0.0,0.0);

}

// Initialize Impedance tensor and errors
void ObservedDataStationMT::initializeImpedanceTensorsAndErrors(){

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
void ObservedDataStationMT::allocateMemoryForCalculatedValues(){

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

// Output calculated Impedance tensors
void ObservedDataStationMT::outputCalculatedValues() const{

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
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

// Calulate interpolator vector of electric field
void ObservedDataStationMT::calcInterpolatorVectorOfElectricField( Forward3D* const ptrForward3D ){

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
	}else if( meshType == MeshData::HEXA){// Hexa mesh
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
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges( m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy );
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

// Calulate sensitivity matrix of Impedance tensors
void ObservedDataStationMT::calculateSensitivityMatrix( const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMFieldExPol,
	const std::complex<double>* const derivativesOfEMFieldEyPol,
	double* const sensitivityMatrix, const bool forceSDToOne ) const{

	//const int freqID = getFreqIDs( freq );
	//const int freqIDAmongThisStationAndPE = getFreqIDsAmongThisPE( freq );

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

//	//----- debug >>>>>
//#ifdef _DEBUG_WRITE
//	std::cout << "divDet " << divDet << std::endl;
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Ex imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*m_rhsVectorIDOfEx + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Ey imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*m_rhsVectorIDOfEy + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Hx imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*rhsVectorIDOfHx   + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Hy imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*rhsVectorIDOfHy   + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Ex imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*m_rhsVectorIDOfEx + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Ey imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*m_rhsVectorIDOfEy + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Hx imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*rhsVectorIDOfHx   + imdl ] << std::endl;
//	}
//	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
//		std::cout << "Hy imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*rhsVectorIDOfHy   + imdl ] << std::endl;
//	}
//#endif
//	//----- debug <<<<<

	//const double baseOfStaticShift = ( ObservedData::getInstance() )->getBaseOfStaticShift();

	//const double factor( pow( 10.0, m_staticShiftFactor ) );

	double cxx = 1.0;
	double cxy = 0.0;
	double cyx = 0.0;
	double cyy = 1.0;
	if( !m_fixDistortionMatrix ){ // Distortion matrix is not fixed
		const int typeOfDistortion = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
		if( typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			assert( m_arrayDistortionMatrixDifferences != NULL );
			cxx = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXX] + 1.0;
			cxy = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXY];
			cyx = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYX];
			cyy = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYY] + 1.0;
		}
		else if( typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
			assert( m_arrayGainsAndRotations != NULL );
			const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
			const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );
			const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
			const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];
			cxx =   gX * cos( betaX );
			cxy = - gY * sin( betaY );
			cyx =   gX * sin( betaX );
			cyy =   gY * cos( betaY );
		}
		else if( typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY ){
			assert( m_arrayGainsAndRotations != NULL );
			const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
			const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );
			cxx =  gX;
			cxy = 0.0;
			cyx = 0.0;
			cyy =  gY;
		}
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

	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const std::complex<double> work1 = derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * HyCalculated[1]
									+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * HxCalculated[0]
									- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * HyCalculated[0]
									- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * HxCalculated[1];

		// dZxx/dm
		const std::complex<double> workXX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] * HyCalculated[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_ExCalculated[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] * HyCalculated[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_ExCalculated[1];

		const std::complex<double> workXX2	= m_ExCalculated[0]*HyCalculated[1]	- m_ExCalculated[1]*HyCalculated[0];

		const double dZxxRealUndist = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).real();  
		const double dZxxImagUndist = ( workXX1 * divDet - work1 * workXX2 * divDet2 ).imag();

		// dZxy/dm
		const std::complex<double> workXY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] * HxCalculated[0]
										+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_ExCalculated[1]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl ] * HxCalculated[1]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_ExCalculated[0];

		const std::complex<double> workXY2	= m_ExCalculated[1]*HxCalculated[0]	- m_ExCalculated[0]*HxCalculated[1];

		const double dZxyRealUndist = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).real();
		const double dZxyImagUndist = ( workXY1 * divDet - work1 * workXY2 * divDet2 ).imag();

		// dZyx/dm
		const std::complex<double> workYX1	= derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HyCalculated[1]
										+ derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_EyCalculated[0]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HyCalculated[0]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHy   + imdl ] * m_EyCalculated[1];

		const std::complex<double> workYX2	= m_EyCalculated[0]*HyCalculated[1]	- m_EyCalculated[1]*HyCalculated[0];

		const double dZyxRealUndist = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).real();
		const double dZyxImagUndist = ( workYX1 * divDet - work1 * workYX2 * divDet2 ).imag();

		// dZyy/dm
		const std::complex<double> workYY1	= derivativesOfEMFieldEyPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HxCalculated[0]
										+ derivativesOfEMFieldExPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_EyCalculated[1]
										- derivativesOfEMFieldExPol[ nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl ] *   HxCalculated[1]
										- derivativesOfEMFieldEyPol[ nBlkNotFixed * rhsVectorIDOfHx   + imdl ] * m_EyCalculated[0];

		const std::complex<double> workYY2	= m_EyCalculated[1]*HxCalculated[0]	- m_EyCalculated[0]*HxCalculated[1];

		const double dZyyRealUndist = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).real();
		const double dZyyImagUndist = ( workYY1 * divDet - work1 * workYY2 * divDet2 ).imag();

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + imdl ] = ( cxx * dZxxRealUndist + cxy * dZyxRealUndist ) / ZxxSD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + imdl ] = ( cxx * dZxyRealUndist + cxy * dZyyRealUndist ) / ZxySD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + imdl ] = ( cyx * dZxxRealUndist + cyy * dZyxRealUndist ) / ZyxSD.realPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + imdl ] = ( cyx * dZxyRealUndist + cyy * dZyyRealUndist ) / ZyySD.realPart;

		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + imdl ] = ( cxx * dZxxImagUndist + cxy * dZyxImagUndist ) / ZxxSD.imagPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + imdl ] = ( cxx * dZxyImagUndist + cxy * dZyyImagUndist ) / ZxySD.imagPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + imdl ] = ( cyx * dZxxImagUndist + cyy * dZyxImagUndist ) / ZyxSD.imagPart;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + imdl ] = ( cyx * dZxyImagUndist + cyy * dZyyImagUndist ) / ZyySD.imagPart;

		if( !m_fixDistortionMatrix ){ // Distortion matrix is not fixed

			const std::complex<double> ZxxUndist = ( m_ExCalculated[0]*HyCalculated[1] - m_ExCalculated[1]*HyCalculated[0] ) / det;
			const std::complex<double> ZxyUndist = ( m_ExCalculated[1]*HxCalculated[0] - m_ExCalculated[0]*HxCalculated[1] ) / det;
			const std::complex<double> ZyxUndist = ( m_EyCalculated[0]*HyCalculated[1] - m_EyCalculated[1]*HyCalculated[0] ) / det;
			const std::complex<double> ZyyUndist = ( m_EyCalculated[1]*HxCalculated[0] - m_EyCalculated[0]*HxCalculated[1] ) / det;

			if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
				const long long ID_Cxx = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CXX]);
				const long long ID_Cxy = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CXY]);
				const long long ID_Cyx = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CYX]);
				const long long ID_Cyy = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CYY]);
				if( ID_Cxx >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx ] = ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx ] = ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx ] = ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx ] = ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx ] = 0.0;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx ] = 0.0;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx ] = 0.0;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx ] = 0.0;// For Im(Zyy)
				}
				if( ID_Cxy >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy ] = ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy ] = ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy ] = ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy ] = ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy ] = 0.0;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy ] = 0.0;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy ] = 0.0;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy ] = 0.0;// For Im(Zyy)
				}
				if( ID_Cyx >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx ] = 0.0;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx ] = 0.0;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx ] = 0.0;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx ] = 0.0;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx ] = ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx ] = ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx ] = ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx ] = ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
				if( ID_Cyy >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy ] = 0.0;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy ] = 0.0;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy ] = 0.0;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy ] = 0.0;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy ] = ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy ] = ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy ] = ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy ] = ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
			}
			else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){

				const double ln10 = log(10.0);

				const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
				const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );
				const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
				const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];

				const long long ID_GainX = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN]);
				const long long ID_GainY = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN]);
				const long long ID_RotX  = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_ROTATION]);  
				const long long ID_RotY  = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_ROTATION]);

				if( ID_GainX >= 0 ){ 
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * cos( betaX ) * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * cos( betaX ) * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * cos( betaX ) * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * cos( betaX ) * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * sin( betaX ) * ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * sin( betaX ) * ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * sin( betaX ) * ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * sin( betaX ) * ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
				if( ID_GainY >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] = - gY * ln10 * sin( betaY ) * ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] = - gY * ln10 * sin( betaY ) * ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] = - gY * ln10 * sin( betaY ) * ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] = - gY * ln10 * sin( betaY ) * ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * cos( betaY ) * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * cos( betaY ) * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * cos( betaY ) * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * cos( betaY ) * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
				if( ID_RotX >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX  ] = - gX *        sin( betaX ) * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX  ] = - gX *        sin( betaX ) * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX  ] = - gX *        sin( betaX ) * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX  ] = - gX *        sin( betaX ) * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX  ] =   gX *        cos( betaX ) * ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX  ] =   gX *        cos( betaX ) * ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX  ] =   gX *        cos( betaX ) * ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX  ] =   gX *        cos( betaX ) * ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
				if( ID_RotY >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY  ] = - gY *        cos( betaY ) * ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY  ] = - gY *        cos( betaY ) * ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY  ] = - gY *        cos( betaY ) * ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY  ] = - gY *        cos( betaY ) * ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY  ] = - gY *        sin( betaY ) * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY  ] = - gY *        sin( betaY ) * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY  ] = - gY *        sin( betaY ) * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY  ] = - gY *        sin( betaY ) * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
			}
			else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){

				const double ln10 = log(10.0);

				const double gX = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN] );
				const double gY = pow( 10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN] );

				const long long ID_GainX = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN]);
				const long long ID_GainY = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN]);
				if( ID_GainX >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   gX * ln10 * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   0.0;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   0.0;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX ] =   0.0;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX ] =   0.0;// For Im(Zyy)
				}
				if( ID_GainY >= 0 ){
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   0.0;// For Re(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   0.0;// For Im(Zxx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   0.0;// For Re(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   0.0;// For Im(Zxy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
					sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY ] =   gY * ln10 * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
				}
			}
		}
	}

}

// Calculate data vector of this PE
void ObservedDataStationMT::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

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
double ObservedDataStationMT::calculateErrorSumOfSquaresThisPE() const{

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

// Copy current distortion parameters to previous ones
void ObservedDataStationMT::copyDistortionParamsCurToPre( const int iComp ){

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp] = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

/// Get flag specifing whether distortion matrix are fixed or not
bool ObservedDataStationMT::doesFixDistortionMatrix() const{

	return m_fixDistortionMatrix;

}

// Get type of the electric field used to calculate response functions
int ObservedDataStationMT::getTypeOfElectricField() const{

	return m_typeOfElectricField;

}

// Set flag specifing whether distortion matrix are fixed or not
void ObservedDataStationMT::setFixDistortionMatrix( const bool doesFix ){

	m_fixDistortionMatrix = doesFix;

}

// Set type of the electric field used to calculate response functions
void ObservedDataStationMT::setTypeOfElectricField( const int type ){

	m_typeOfElectricField = type;

}

// Set distortion parameters of previous iteration
void ObservedDataStationMT::setDistortionParamsPre( const int iComp, const double val ){

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp] = val;
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = val;
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = val;
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}


}

// Set distortion parameters
void ObservedDataStationMT::setDistortionParams( const int iComp, const double val ){

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp] = val;
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );
		m_arrayGainsAndRotations->gainsAndRotations[iComp] = val;
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );
		m_arrayGainsAndRotations->gainsAndRotations[iComp] = val;
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Set ID of distortion parameters
void ObservedDataStationMT::setIDOfDistortionParams( const int iComp, const int ID ){

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		if( doesFixDistortionMatrix() ){
			m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] = -1;
			OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
		}else{// Distortion matrix is not fixed
			m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] = ID;
		}

	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );
		if( doesFixDistortionMatrix() ){
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = -1;
			OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
		}else{// Distortion matrix is not fixed
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = ID;
		}
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );
		if( doesFixDistortionMatrix() ){
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = -1;
			OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
		}else{// Distortion matrix is not fixed
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = ID;
		}
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Set full updated value of distortion parameters
void ObservedDataStationMT::setDistortionParamsUpdatedFull( const int iComp, const double val ){

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		if( doesFixDistortionMatrix() ){// Not change distortion matrix
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] = m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
		}else{
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] = val;
		}
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );
		if( doesFixDistortionMatrix() ){// Not change distortion matrix
			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
		}else{
			if( iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION ){
				// Rotation angle is bounded in from -pi/2 to pi/2 radians
				if( val > 0.5 * CommonParameters::PI ){
					m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = 0.5 * CommonParameters::PI;
				}
				else if( val < - 0.5 * CommonParameters::PI ){
					m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = - 0.5 * CommonParameters::PI;
				}
				else{
					m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
				}
			}
			else{
				m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
			}
		}

	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );
		if( doesFixDistortionMatrix() ){// Not change distortion matrix
			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
		}else{
			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
		}
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}
	
}

// Update distortion parameters
void ObservedDataStationMT::updateDistortionParams( const double dampingFactor ){

	if( doesFixDistortionMatrix() ){// distortion matrix if fixed
		return;
	}

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();

	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );

		for( int iComp = 0; iComp < 4; ++iComp ){
			if( m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] >= 0 ){
				m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp] = 
					dampingFactor * m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] + ( 1.0 - dampingFactor ) * m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
			}
		}

	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );

		for( int iComp = 0; iComp < 4; ++iComp ){
			if( m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] >= 0 ){
				if( iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION ){
					// Rotation angle is bounded in from -pi/2 to pi/2 radians
					// y = tan(beta)
					// beta = Arctan(y)
					// dy/dbeta = 1 / cos(beta)^2
					const double yUpdatedFull = tan(m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp]);
					const double yPre = tan(m_arrayGainsAndRotations->gainsAndRotationsPre[iComp]);
					const double yNew = dampingFactor * yUpdatedFull + ( 1.0 - dampingFactor ) * yPre;
					const double betaNew = atan(yNew);
					m_arrayGainsAndRotations->gainsAndRotations[iComp] = betaNew;
				}else{
					m_arrayGainsAndRotations->gainsAndRotations[iComp] = 
						dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] + ( 1.0 - dampingFactor ) * m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
				}
			}
		}		

	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		if( m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN] >= 0 ){
			m_arrayGainsAndRotations->gainsAndRotations[EX_GAIN] = 
				dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[EX_GAIN] + ( 1.0 - dampingFactor ) * m_arrayGainsAndRotations->gainsAndRotationsPre[EX_GAIN];
		}
		if( m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN] >= 0 ){
			m_arrayGainsAndRotations->gainsAndRotations[EY_GAIN] = 
				dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[EY_GAIN] + ( 1.0 - dampingFactor ) * m_arrayGainsAndRotations->gainsAndRotationsPre[EY_GAIN];
		}
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Get distortion parameters of previous iteration
double ObservedDataStationMT::getDistortionParamsPre( const int iComp ) const{

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );

		return m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );

		return m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );

		return m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Get distortion parameters
double ObservedDataStationMT::getDistortionParams( const int iComp ) const{

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );

		return m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );

		return m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );

		return m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Get ID of distortion parameters
int ObservedDataStationMT::getIDOfDistortionParams( const int iComp ) const{

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );

		return m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );

		return m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );

		return m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp];
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}

// Get full updated value of distortion parameters
double ObservedDataStationMT::getDistortionParamsUpdatedFull( const int iComp ) const{

	const int type = ( AnalysisControl::getInstance() )->getTypeOfDistortion();
	if( type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		assert( m_arrayDistortionMatrixDifferences != NULL );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );
		assert( iComp >= COMPONENT_ID_CXX );
		assert( iComp <= COMPONENT_ID_CYY );

		return m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp >= EX_GAIN );
		assert( iComp <= EY_ROTATION );

		return m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp];
	}
	else if( type == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		assert( m_arrayGainsAndRotations != NULL );
		assert( iComp == EX_GAIN || iComp == EY_GAIN );

		return m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp];
	}
	else{
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}

}
