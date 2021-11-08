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

#include "ObservedDataStationApparentResistivityAndPhase.h"
#include "ObservedDataStationMT.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"
#include "Util.h"

// Constructer
ObservedDataStationApparentResistivityAndPhase::ObservedDataStationApparentResistivityAndPhase():
	ObservedDataStationMT(),
	m_apparentResistivityXXObserved(NULL),
	m_apparentResistivityXYObserved(NULL),
	m_apparentResistivityYXObserved(NULL),
	m_apparentResistivityYYObserved(NULL),
	m_PhaseXXObserved(NULL),
	m_PhaseXYObserved(NULL),
	m_PhaseYXObserved(NULL),
	m_PhaseYYObserved(NULL),
	m_apparentResistivityXXSD(NULL),
	m_apparentResistivityXYSD(NULL),
	m_apparentResistivityYXSD(NULL),
	m_apparentResistivityYYSD(NULL),
	m_PhaseXXSD(NULL),
	m_PhaseXYSD(NULL),
	m_PhaseYXSD(NULL),
	m_PhaseYYSD(NULL),
	m_apparentResistivityXXCalculated(NULL),
	m_apparentResistivityXYCalculated(NULL),
	m_apparentResistivityYXCalculated(NULL),
	m_apparentResistivityYYCalculated(NULL),
	m_PhaseXXCalculated(NULL),
	m_PhaseXYCalculated(NULL),
	m_PhaseYXCalculated(NULL),
	m_PhaseYYCalculated(NULL),
	m_apparentResistivityXXResidual(NULL),
	m_apparentResistivityXYResidual(NULL),
	m_apparentResistivityYXResidual(NULL),
	m_apparentResistivityYYResidual(NULL),
	m_PhaseXXResidual(NULL),
	m_PhaseXYResidual(NULL),
	m_PhaseYXResidual(NULL),
	m_PhaseYYResidual(NULL),
	m_dataIDOfApparentResistivityXX(NULL),
	m_dataIDOfApparentResistivityXY(NULL),
	m_dataIDOfApparentResistivityYX(NULL),
	m_dataIDOfApparentResistivityYY(NULL),
	m_dataIDOfPhaseXX(NULL),
	m_dataIDOfPhaseXY(NULL),
	m_dataIDOfPhaseYX(NULL),
	m_dataIDOfPhaseYY(NULL)
{
}

// Destructer
ObservedDataStationApparentResistivityAndPhase::~ObservedDataStationApparentResistivityAndPhase(){

	if( m_apparentResistivityXXObserved != NULL){
		delete[] m_apparentResistivityXXObserved;
		m_apparentResistivityXXObserved = NULL;
	}

	if( m_apparentResistivityXYObserved != NULL){
		delete[] m_apparentResistivityXYObserved;
		m_apparentResistivityXYObserved = NULL;
	}

	if( m_apparentResistivityYXObserved != NULL){
		delete[] m_apparentResistivityYXObserved;
		m_apparentResistivityYXObserved = NULL;
	}

	if( m_apparentResistivityYYObserved != NULL){
		delete[] m_apparentResistivityYYObserved;
		m_apparentResistivityYYObserved = NULL;
	}

	if( m_PhaseXXObserved != NULL){
		delete[] m_PhaseXXObserved;
		m_PhaseXXObserved = NULL;
	}

	if( m_PhaseXYObserved != NULL){
		delete[] m_PhaseXYObserved;
		m_PhaseXYObserved = NULL;
	}

	if( m_PhaseYXObserved != NULL){
		delete[] m_PhaseYXObserved;
		m_PhaseYXObserved = NULL;
	}

	if( m_PhaseYYObserved != NULL){
		delete[] m_PhaseYYObserved;
		m_PhaseYYObserved = NULL;
	}

	if( m_apparentResistivityXXSD != NULL){
		delete[] m_apparentResistivityXXSD;
		m_apparentResistivityXXSD = NULL;
	}

	if( m_apparentResistivityXYSD != NULL){
		delete[] m_apparentResistivityXYSD;
		m_apparentResistivityXYSD = NULL;
	}

	if( m_apparentResistivityYXSD != NULL){
		delete[] m_apparentResistivityYXSD;
		m_apparentResistivityYXSD = NULL;
	}

	if( m_apparentResistivityYYSD != NULL){
		delete[] m_apparentResistivityYYSD;
		m_apparentResistivityYYSD = NULL;
	}

	if( m_PhaseXXSD != NULL){
		delete[] m_PhaseXXSD;
		m_PhaseXXSD = NULL;
	}

	if( m_PhaseXYSD != NULL){
		delete[] m_PhaseXYSD;
		m_PhaseXYSD = NULL;
	}

	if( m_PhaseYXSD != NULL){
		delete[] m_PhaseYXSD;
		m_PhaseYXSD = NULL;
	}

	if( m_PhaseYYSD != NULL){
		delete[] m_PhaseYYSD;
		m_PhaseYYSD = NULL;
	}

	if( m_apparentResistivityXXCalculated != NULL){
		delete[] m_apparentResistivityXXCalculated;
		m_apparentResistivityXXCalculated = NULL;
	}

	if( m_apparentResistivityXYCalculated != NULL){
		delete[] m_apparentResistivityXYCalculated;
		m_apparentResistivityXYCalculated = NULL;
	}

	if( m_apparentResistivityYXCalculated != NULL){
		delete[] m_apparentResistivityYXCalculated;
		m_apparentResistivityYXCalculated = NULL;
	}

	if( m_apparentResistivityYYCalculated != NULL){
		delete[] m_apparentResistivityYYCalculated;
		m_apparentResistivityYYCalculated = NULL;
	}

	if( m_PhaseXXCalculated != NULL){
		delete[] m_PhaseXXCalculated;
		m_PhaseXXCalculated = NULL;
	}

	if( m_PhaseXYCalculated != NULL){
		delete[] m_PhaseXYCalculated;
		m_PhaseXYCalculated = NULL;
	}

	if( m_PhaseYXCalculated != NULL){
		delete[] m_PhaseYXCalculated;
		m_PhaseYXCalculated = NULL;
	}

	if( m_PhaseYYCalculated != NULL){
		delete[] m_PhaseYYCalculated;
		m_PhaseYYCalculated = NULL;
	}

	if( m_apparentResistivityXXResidual != NULL){
		delete[] m_apparentResistivityXXResidual;
		m_apparentResistivityXXResidual = NULL;
	}

	if( m_apparentResistivityXYResidual != NULL){
		delete[] m_apparentResistivityXYResidual;
		m_apparentResistivityXYResidual = NULL;
	}

	if( m_apparentResistivityYXResidual != NULL){
		delete[] m_apparentResistivityYXResidual;
		m_apparentResistivityYXResidual = NULL;
	}

	if( m_apparentResistivityYYResidual != NULL){
		delete[] m_apparentResistivityYYResidual;
		m_apparentResistivityYYResidual = NULL;
	}

	if( m_PhaseXXResidual != NULL){
		delete[] m_PhaseXXResidual;
		m_PhaseXXResidual = NULL;
	}

	if( m_PhaseXYResidual != NULL){
		delete[] m_PhaseXYResidual;
		m_PhaseXYResidual = NULL;
	}

	if( m_PhaseYXResidual != NULL){
		delete[] m_PhaseYXResidual;
		m_PhaseYXResidual = NULL;
	}

	if( m_PhaseYYResidual != NULL){
		delete[] m_PhaseYYResidual;
		m_PhaseYYResidual = NULL;
	}

	if( m_dataIDOfApparentResistivityXX != NULL){
		delete[] m_dataIDOfApparentResistivityXX;
		m_dataIDOfApparentResistivityXX = NULL;
	}

	if( m_dataIDOfApparentResistivityXY != NULL){
		delete[] m_dataIDOfApparentResistivityXY;
		m_dataIDOfApparentResistivityXY = NULL;
	}

	if( m_dataIDOfApparentResistivityYX != NULL){
		delete[] m_dataIDOfApparentResistivityYX;
		m_dataIDOfApparentResistivityYX = NULL;
	}

	if( m_dataIDOfApparentResistivityYY != NULL){
		delete[] m_dataIDOfApparentResistivityYY;
		m_dataIDOfApparentResistivityYY = NULL;
	}

	if( m_dataIDOfPhaseXX != NULL){
		delete[] m_dataIDOfPhaseXX;
		m_dataIDOfPhaseXX = NULL;
	}

	if( m_dataIDOfPhaseXY != NULL){
		delete[] m_dataIDOfPhaseXY;
		m_dataIDOfPhaseXY = NULL;
	}

	if( m_dataIDOfPhaseYX != NULL){
		delete[] m_dataIDOfPhaseYX;
		m_dataIDOfPhaseYX = NULL;
	}

	if( m_dataIDOfPhaseYY != NULL){
		delete[] m_dataIDOfPhaseYY;
		m_dataIDOfPhaseYY = NULL;
	}

}

// Read data from input file
void ObservedDataStationApparentResistivityAndPhase::inputObservedData( std::ifstream& inFile ){

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
		int ibuf(-1);
		inFile >> ibuf;
		setTypeOfElectricField(ibuf);
	}else{
		setTypeOfElectricField( pAnalysisControl->getTypeOfElectricField() );
	}

	if( pAnalysisControl->getTypeOfMesh() == MeshData::HEXA ){
		if(	getTypeOfElectricField() != AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD ){
			OutputFiles::m_logFile << std::endl << "Warning : Horizontal electric field must be used for hexahedral mesh." << std::endl;
		}
		setTypeOfElectricField( AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD );
	}
	std::string elecType;
	switch( getTypeOfElectricField() ){
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			elecType = "Tangential";
			break;
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			elecType = "Horizontal";
			break;
		default:
			OutputFiles::m_logFile << std::endl << "Error : Unknown type of the electric field : " << getTypeOfElectricField() << std::endl;
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
		m_apparentResistivityXXObserved = new double[nFreq];
		m_apparentResistivityXYObserved = new double[nFreq];
		m_apparentResistivityYXObserved = new double[nFreq];
		m_apparentResistivityYYObserved = new double[nFreq];
		m_PhaseXXObserved = new double[nFreq];
		m_PhaseXYObserved = new double[nFreq];
		m_PhaseYXObserved = new double[nFreq];
		m_PhaseYYObserved = new double[nFreq];
		m_apparentResistivityXXSD = new double[nFreq];
		m_apparentResistivityXYSD = new double[nFreq];
		m_apparentResistivityYXSD = new double[nFreq];
		m_apparentResistivityYYSD = new double[nFreq];
		m_PhaseXXSD = new double[nFreq];
		m_PhaseXYSD = new double[nFreq];
		m_PhaseYXSD = new double[nFreq];
		m_PhaseYYSD = new double[nFreq];

		// Their arrays are used in some functions of the base class
		m_ZxxObserved = new std::complex<double>[nFreq];
		m_ZxyObserved = new std::complex<double>[nFreq];
		m_ZyxObserved = new std::complex<double>[nFreq];
		m_ZyyObserved = new std::complex<double>[nFreq];
		m_ZxxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZxySD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyxSD = new CommonParameters::DoubleComplexValues[nFreq];
		m_ZyySD = new CommonParameters::DoubleComplexValues[nFreq];

		for( int i = 0; i < nFreq; ++i ){
			inFile >> m_freq[i];
			double dbuf1(0.0);
			double dbuf2(0.0);
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityXXObserved[i] = dbuf1;
			m_PhaseXXObserved[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityXYObserved[i] = dbuf1;
			m_PhaseXYObserved[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityYXObserved[i] = dbuf1;
			m_PhaseYXObserved[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityYYObserved[i] = dbuf1;
			m_PhaseYYObserved[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityXXSD[i] = dbuf1;
			m_PhaseXXSD[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityXYSD[i] = dbuf1;
			m_PhaseXYSD[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityYXSD[i] = dbuf1;
			m_PhaseYXSD[i] = dbuf2;
			inFile >> dbuf1 >> dbuf2;
			m_apparentResistivityYYSD[i] = dbuf1;
			m_PhaseYYSD[i] = dbuf2;
			
			if( m_apparentResistivityXXSD[i] > 0.0 && m_PhaseXXSD[i] > 0.0 ){
				calcImpedanceTensorComponentFromApparentResistivityAndPhase( m_freq[i], m_apparentResistivityXXObserved[i], m_apparentResistivityXXSD[i], m_PhaseXXObserved[i], m_PhaseXXSD[i], m_ZxxObserved[i], m_ZxxSD[i] );
			}else{
				m_ZxxObserved[i] = std::complex<double>(0.0, 0.0);
				m_ZxxSD[i].realPart = 1.0e10;
				m_ZxxSD[i].imagPart = 1.0e10;
			}
			if( m_apparentResistivityXYSD[i] > 0.0 && m_PhaseXYSD[i] > 0.0 ){
				calcImpedanceTensorComponentFromApparentResistivityAndPhase( m_freq[i], m_apparentResistivityXYObserved[i], m_apparentResistivityXYSD[i], m_PhaseXYObserved[i], m_PhaseXYSD[i], m_ZxyObserved[i], m_ZxySD[i] );
			}else{
				m_ZxyObserved[i] = std::complex<double>(0.0, 0.0);
				m_ZxySD[i].realPart = 1.0e10;
				m_ZxySD[i].imagPart = 1.0e10;
			}
			if( m_apparentResistivityYXSD[i] > 0.0 && m_PhaseYXSD[i] > 0.0 ){
				calcImpedanceTensorComponentFromApparentResistivityAndPhase( m_freq[i], m_apparentResistivityYXObserved[i], m_apparentResistivityYXSD[i], m_PhaseYXObserved[i], m_PhaseYXSD[i], m_ZyxObserved[i], m_ZyxSD[i] );
			}else{
				m_ZyxObserved[i] = std::complex<double>(0.0, 0.0);
				m_ZyxSD[i].realPart = 1.0e10;
				m_ZyxSD[i].imagPart = 1.0e10;
			}
			if( m_apparentResistivityYYSD[i] > 0.0 && m_PhaseYYSD[i] > 0.0 ){
				calcImpedanceTensorComponentFromApparentResistivityAndPhase( m_freq[i], m_apparentResistivityYYObserved[i], m_apparentResistivityYYSD[i], m_PhaseYYObserved[i], m_PhaseYYSD[i], m_ZyyObserved[i], m_ZyySD[i] );
			}else{
				m_ZyyObserved[i] = std::complex<double>(0.0, 0.0);
				m_ZyySD[i].realPart = 1.0e10;
				m_ZyySD[i].imagPart = 1.0e10;
			}
		}
	}

}

// Calulate Impedance tensor
void ObservedDataStationApparentResistivityAndPhase::calculateApparentResistivityAndPhase( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount ){

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );
	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}
	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	ObservedDataStationMT::calculateImpedanceTensor( freq, ptrStationOfMagneticField, icount );

	const double omega = 2.0 * CommonParameters::PI * freq;

	const std::complex<double> Zxx = m_ZxxCalculated[freqIDThisPEInSta];
	const std::complex<double> Zxy = m_ZxyCalculated[freqIDThisPEInSta];
	const std::complex<double> Zyx = m_ZyxCalculated[freqIDThisPEInSta];
	const std::complex<double> Zyy = m_ZyyCalculated[freqIDThisPEInSta];
	m_apparentResistivityXXCalculated[freqIDThisPEInSta] = std::norm( Zxx ) / ( omega * CommonParameters::mu );
	m_apparentResistivityXYCalculated[freqIDThisPEInSta] = std::norm( Zxy ) / ( omega * CommonParameters::mu );
	m_apparentResistivityYXCalculated[freqIDThisPEInSta] = std::norm( Zyx ) / ( omega * CommonParameters::mu );
	m_apparentResistivityYYCalculated[freqIDThisPEInSta] = std::norm( Zyy ) / ( omega * CommonParameters::mu );
	m_PhaseXXCalculated[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2( Zxx.imag(), Zxx.real() );
	m_PhaseXYCalculated[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2( Zxy.imag(), Zxy.real() );
	m_PhaseYXCalculated[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2( Zyx.imag(), Zyx.real() );
	m_PhaseYYCalculated[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2( Zyy.imag(), Zyy.real() );

	// Zero clear
	m_apparentResistivityXXResidual[freqIDThisPEInSta] = 0.0;
	m_apparentResistivityXYResidual[freqIDThisPEInSta] = 0.0;
	m_apparentResistivityYXResidual[freqIDThisPEInSta] = 0.0;
	m_apparentResistivityYYResidual[freqIDThisPEInSta] = 0.0;
	m_PhaseXXResidual[freqIDThisPEInSta] = 0.0;
	m_PhaseXYResidual[freqIDThisPEInSta] = 0.0;
	m_PhaseYXResidual[freqIDThisPEInSta] = 0.0;
	m_PhaseYYResidual[freqIDThisPEInSta] = 0.0;

	if( m_apparentResistivityXXSD[freqIDGlobalInSta] > 0.0 ){
		m_apparentResistivityXXResidual[freqIDThisPEInSta] =
			log10( m_apparentResistivityXXObserved[freqIDGlobalInSta] / m_apparentResistivityXXCalculated[freqIDThisPEInSta] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XX );
	}
	if( m_apparentResistivityXYSD[freqIDGlobalInSta] > 0.0 ){
		m_apparentResistivityXYResidual[freqIDThisPEInSta] =
			log10( m_apparentResistivityXYObserved[freqIDGlobalInSta] / m_apparentResistivityXYCalculated[freqIDThisPEInSta] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XY );
	}
	if( m_apparentResistivityYXSD[freqIDGlobalInSta] > 0.0 ){
		m_apparentResistivityYXResidual[freqIDThisPEInSta] =
			log10( m_apparentResistivityYXObserved[freqIDGlobalInSta] / m_apparentResistivityYXCalculated[freqIDThisPEInSta] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YX );
	}
	if( m_apparentResistivityYYSD[freqIDGlobalInSta] > 0.0 ){
		m_apparentResistivityYYResidual[freqIDThisPEInSta] =
			log10( m_apparentResistivityYYObserved[freqIDGlobalInSta] / m_apparentResistivityYYCalculated[freqIDThisPEInSta] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YY );
	}

	if( m_PhaseXXSD[freqIDGlobalInSta] > 0.0 ){	
		double phaseObs = m_PhaseXXObserved[freqIDGlobalInSta];
		if( Zxx.real() < 0.0 && m_ZxxObserved[freqIDGlobalInSta].real() < 0.0 && Zxx.imag() * m_ZxxObserved[freqIDGlobalInSta].imag() < 0.0 ){
			if( Zxx.imag() > 0.0 ){
				phaseObs +=  360.0;
			}else{
				phaseObs -=  360.0;
			}
		}
		m_PhaseXXResidual[freqIDThisPEInSta] = ( phaseObs - m_PhaseXXCalculated[freqIDThisPEInSta] ) / m_PhaseXXSD[freqIDGlobalInSta];
	}
	if( m_PhaseXYSD[freqIDGlobalInSta] > 0.0 ){
		double phaseObs = m_PhaseXYObserved[freqIDGlobalInSta];
		if( Zxy.real() < 0.0 && m_ZxyObserved[freqIDGlobalInSta].real() < 0.0 && Zxy.imag() * m_ZxyObserved[freqIDGlobalInSta].imag() < 0.0 ){
			if( Zxy.imag() > 0.0 ){
				phaseObs += 360.0;
			}else{
				phaseObs -= 360.0;
			}
		}
		m_PhaseXYResidual[freqIDThisPEInSta] = ( phaseObs - m_PhaseXYCalculated[freqIDThisPEInSta] ) / m_PhaseXYSD[freqIDGlobalInSta];
	}
	if( m_PhaseYXSD[freqIDGlobalInSta] > 0.0 ){
		double phaseObs = m_PhaseYXObserved[freqIDGlobalInSta];
		if( Zyx.real() < 0.0 && m_ZyxObserved[freqIDGlobalInSta].real() < 0.0 && Zyx.imag() * m_ZyxObserved[freqIDGlobalInSta].imag() < 0.0 ){
			if( Zyx.imag() > 0.0 ){
				phaseObs += 360.0;
			}else{
				phaseObs -= 360.0;
			}
		}
		m_PhaseYXResidual[freqIDThisPEInSta] = ( phaseObs - m_PhaseYXCalculated[freqIDThisPEInSta] ) / m_PhaseYXSD[freqIDGlobalInSta];
	}
	if( m_PhaseYYSD[freqIDGlobalInSta] > 0.0 ){
		double phaseObs = m_PhaseYYObserved[freqIDGlobalInSta];
		if( Zyy.real() < 0.0 && m_ZyyObserved[freqIDGlobalInSta].real() < 0.0 && Zyy.imag() * m_ZyyObserved[freqIDGlobalInSta].imag() < 0.0 ){
			if( Zyy.imag() > 0.0 ){
				phaseObs += 360.0;
			}else{
				phaseObs -= 360.0;
			}
		}
		m_PhaseYYResidual[freqIDThisPEInSta] = ( phaseObs - m_PhaseYYCalculated[freqIDThisPEInSta] ) / m_PhaseYYSD[freqIDGlobalInSta];
	}

	// Share data ID with the ones for impedance tensors
	m_dataIDOfApparentResistivityXX[freqIDThisPEInSta] = m_dataIDOfZxx[freqIDThisPEInSta].realPart;
	m_dataIDOfApparentResistivityXY[freqIDThisPEInSta] = m_dataIDOfZxy[freqIDThisPEInSta].realPart;
	m_dataIDOfApparentResistivityYX[freqIDThisPEInSta] = m_dataIDOfZyx[freqIDThisPEInSta].realPart;
	m_dataIDOfApparentResistivityYY[freqIDThisPEInSta] = m_dataIDOfZyy[freqIDThisPEInSta].realPart;
	m_dataIDOfPhaseXX[freqIDThisPEInSta] = m_dataIDOfZxx[freqIDThisPEInSta].imagPart;
	m_dataIDOfPhaseXY[freqIDThisPEInSta] = m_dataIDOfZxy[freqIDThisPEInSta].imagPart;
	m_dataIDOfPhaseYX[freqIDThisPEInSta] = m_dataIDOfZyx[freqIDThisPEInSta].imagPart;
	m_dataIDOfPhaseYY[freqIDThisPEInSta] = m_dataIDOfZyy[freqIDThisPEInSta].imagPart;

}

// Initialize apparent resistivitu, phase and their errors
void ObservedDataStationApparentResistivityAndPhase::initializeApparentResistivityPhaseAndErrors(){

	ObservedDataStationMT::initializeImpedanceTensorsAndErrors();

	//for( int i = 0; i < m_numOfFrequency; ++i ){
	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		m_apparentResistivityXXCalculated[i] = 0.0;
		m_apparentResistivityXYCalculated[i] = 0.0;
		m_apparentResistivityYXCalculated[i] = 0.0;
		m_apparentResistivityYYCalculated[i] = 0.0;
		m_PhaseXXCalculated[i] = 0.0;
		m_PhaseXYCalculated[i] = 0.0;
		m_PhaseYXCalculated[i] = 0.0;
		m_PhaseYYCalculated[i] = 0.0;
		m_apparentResistivityXXResidual[i] = 0.0;
		m_apparentResistivityXYResidual[i] = 0.0;
		m_apparentResistivityYXResidual[i] = 0.0;
		m_apparentResistivityYYResidual[i] = 0.0;
		m_PhaseXXResidual[i] = 0.0;
		m_PhaseXYResidual[i] = 0.0;
		m_PhaseYXResidual[i] = 0.0;
		m_PhaseYYResidual[i] = 0.0;
	}

}

// Allocate memory for the calculated Impedance and errors
void ObservedDataStationApparentResistivityAndPhase::allocateMemoryForCalculatedValues(){

	ObservedDataStationMT::allocateMemoryForCalculatedValues();

	if( m_apparentResistivityXXCalculated != NULL){
		delete[] m_apparentResistivityXXCalculated;
		m_apparentResistivityXXCalculated = NULL;
	}

	if( m_apparentResistivityXYCalculated != NULL){
		delete[] m_apparentResistivityXYCalculated;
		m_apparentResistivityXYCalculated = NULL;
	}

	if( m_apparentResistivityYXCalculated != NULL){
		delete[] m_apparentResistivityYXCalculated;
		m_apparentResistivityYXCalculated = NULL;
	}

	if( m_apparentResistivityYYCalculated != NULL){
		delete[] m_apparentResistivityYYCalculated;
		m_apparentResistivityYYCalculated = NULL;
	}

	if( m_PhaseXXCalculated != NULL){
		delete[] m_PhaseXXCalculated;
		m_PhaseXXCalculated = NULL;
	}

	if( m_PhaseXYCalculated != NULL){
		delete[] m_PhaseXYCalculated;
		m_PhaseXYCalculated = NULL;
	}

	if( m_PhaseYXCalculated != NULL){
		delete[] m_PhaseYXCalculated;
		m_PhaseYXCalculated = NULL;
	}

	if( m_PhaseYYCalculated != NULL){
		delete[] m_PhaseYYCalculated;
		m_PhaseYYCalculated = NULL;
	}

	if( m_apparentResistivityXXResidual != NULL){
		delete[] m_apparentResistivityXXResidual;
		m_apparentResistivityXXResidual = NULL;
	}

	if( m_apparentResistivityXYResidual != NULL){
		delete[] m_apparentResistivityXYResidual;
		m_apparentResistivityXYResidual = NULL;
	}

	if( m_apparentResistivityYXResidual != NULL){
		delete[] m_apparentResistivityYXResidual;
		m_apparentResistivityYXResidual = NULL;
	}

	if( m_apparentResistivityYYResidual != NULL){
		delete[] m_apparentResistivityYYResidual;
		m_apparentResistivityYYResidual = NULL;
	}

	if( m_PhaseXXResidual != NULL){
		delete[] m_PhaseXXResidual;
		m_PhaseXXResidual = NULL;
	}

	if( m_PhaseXYResidual != NULL){
		delete[] m_PhaseXYResidual;
		m_PhaseXYResidual = NULL;
	}

	if( m_PhaseYXResidual != NULL){
		delete[] m_PhaseYXResidual;
		m_PhaseYXResidual = NULL;
	}

	if( m_PhaseYYResidual != NULL){
		delete[] m_PhaseYYResidual;
		m_PhaseYYResidual = NULL;
	}

	if( m_dataIDOfApparentResistivityXX != NULL){
		delete[] m_dataIDOfApparentResistivityXX;
		m_dataIDOfApparentResistivityXX = NULL;
	}

	if( m_dataIDOfApparentResistivityXY != NULL){
		delete[] m_dataIDOfApparentResistivityXY;
		m_dataIDOfApparentResistivityXY = NULL;
	}

	if( m_dataIDOfApparentResistivityYX != NULL){
		delete[] m_dataIDOfApparentResistivityYX;
		m_dataIDOfApparentResistivityYX = NULL;
	}

	if( m_dataIDOfApparentResistivityYY != NULL){
		delete[] m_dataIDOfApparentResistivityYY;
		m_dataIDOfApparentResistivityYY = NULL;
	}

	if( m_dataIDOfPhaseXX != NULL){
		delete[] m_dataIDOfPhaseXX;
		m_dataIDOfPhaseXX = NULL;
	}

	if( m_dataIDOfPhaseXY != NULL){
		delete[] m_dataIDOfPhaseXY;
		m_dataIDOfPhaseXY = NULL;
	}

	if( m_dataIDOfPhaseYX != NULL){
		delete[] m_dataIDOfPhaseYX;
		m_dataIDOfPhaseYX = NULL;
	}

	if( m_dataIDOfPhaseYY != NULL){
		delete[] m_dataIDOfPhaseYY;
		m_dataIDOfPhaseYY = NULL;
	}

	if( m_numOfFreqCalculatedByThisStaAndPE > 0 ){

		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_apparentResistivityXXCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityXYCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityYXCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityYYCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseXXCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseXYCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseYXCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseYYCalculated = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityXXResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityXYResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityYXResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_apparentResistivityYYResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseXXResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseXYResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseYXResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_PhaseYYResidual = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfApparentResistivityXX = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfApparentResistivityXY = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfApparentResistivityYX = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfApparentResistivityYY = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPhaseXX = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPhaseXY = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPhaseYX = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfPhaseYY = new int[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
			m_apparentResistivityXXCalculated[i] = 0.0;
			m_apparentResistivityXYCalculated[i] = 0.0;
			m_apparentResistivityYXCalculated[i] = 0.0;
			m_apparentResistivityYYCalculated[i] = 0.0;
			m_PhaseXXCalculated[i] = 0.0;
			m_PhaseXYCalculated[i] = 0.0;
			m_PhaseYXCalculated[i] = 0.0;
			m_PhaseYYCalculated[i] = 0.0;
			m_apparentResistivityXXResidual[i] = 0.0;
			m_apparentResistivityXYResidual[i] = 0.0;
			m_apparentResistivityYXResidual[i] = 0.0;
			m_apparentResistivityYYResidual[i] = 0.0;
			m_PhaseXXResidual[i] = 0.0;
			m_PhaseXYResidual[i] = 0.0;
			m_PhaseYXResidual[i] = 0.0;
			m_PhaseYYResidual[i] = 0.0;
			m_dataIDOfApparentResistivityXX[i] = -1;
			m_dataIDOfApparentResistivityXY[i] = -1;
			m_dataIDOfApparentResistivityYX[i] = -1;
			m_dataIDOfApparentResistivityYY[i] = -1;
			m_dataIDOfPhaseXX[i] = -1;
			m_dataIDOfPhaseXY[i] = -1;
			m_dataIDOfPhaseYX[i] = -1;
			m_dataIDOfPhaseYY[i] = -1;
		}
	}

}

// Output calculated Impedance tensors
void ObservedDataStationApparentResistivityAndPhase::outputCalculatedValues() const{

//#ifdef _DEBUG_WRITE
//	ObservedDataStationMT::outputCalculatedValues();
//#endif
	const bool useImpedanceTensorInsteadOfPhase = ( AnalysisControl::getInstance()->getApparentResistivityAndPhaseTreatmentOption() == AnalysisControl::USE_Z_IF_SIGN_OF_RE_Z_DIFFER );

	int icount(0);
	for( std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		fprintf( OutputFiles::m_csvFile, "%10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_stationID, m_freq[*itr], 
			m_apparentResistivityXXCalculated[icount], m_PhaseXXCalculated[icount],
			m_apparentResistivityXYCalculated[icount], m_PhaseXYCalculated[icount],
			m_apparentResistivityYXCalculated[icount], m_PhaseYXCalculated[icount],
			m_apparentResistivityYYCalculated[icount], m_PhaseYYCalculated[icount]
			);
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_apparentResistivityXXResidual[icount], m_PhaseXXResidual[icount],
			m_apparentResistivityXYResidual[icount], m_PhaseXYResidual[icount],
			m_apparentResistivityYXResidual[icount], m_PhaseYXResidual[icount],
			m_apparentResistivityYYResidual[icount], m_PhaseYYResidual[icount]
			);
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_apparentResistivityXXObserved[*itr], m_PhaseXXObserved[*itr],
			m_apparentResistivityXYObserved[*itr], m_PhaseXYObserved[*itr],
			m_apparentResistivityYXObserved[*itr], m_PhaseYXObserved[*itr],
			m_apparentResistivityYYObserved[*itr], m_PhaseYYObserved[*itr]
			);
 		fprintf( OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,", 
			m_apparentResistivityXXSD[*itr], m_PhaseXXSD[*itr],
			m_apparentResistivityXYSD[*itr], m_PhaseXYSD[*itr],
			m_apparentResistivityYXSD[*itr], m_PhaseYXSD[*itr],
			m_apparentResistivityYYSD[*itr], m_PhaseYYSD[*itr]
			);
		if(useImpedanceTensorInsteadOfPhase){
	 		fprintf( OutputFiles::m_csvFile, "%3s,%3s,%3s,%3s,\n",
			isUsedImpedanceTensorFromFreqIDs( icount, ObservedDataStationMT::XX ) ? "*" : "",
			isUsedImpedanceTensorFromFreqIDs( icount, ObservedDataStationMT::XY ) ? "*" : "",
			isUsedImpedanceTensorFromFreqIDs( icount, ObservedDataStationMT::YX ) ? "*" : "",
			isUsedImpedanceTensorFromFreqIDs( icount, ObservedDataStationMT::YY ) ? "*" : ""
			);
		}else{
	 		fprintf( OutputFiles::m_csvFile, "\n");
		}
		++icount;
	}

}

// Calulate sensitivity matrix of apparent resistivity and phase
void ObservedDataStationApparentResistivityAndPhase::calculateSensitivityMatrix( const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMFieldExPol,
	const std::complex<double>* const derivativesOfEMFieldEyPol,
	double* const sensitivityMatrix ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );
	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}
	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];

	ObservedDataStationMT::calculateSensitivityMatrix( freq, nModel, ptrStationOfMagneticField, derivativesOfEMFieldExPol, derivativesOfEMFieldEyPol, sensitivityMatrix, true );

	const std::complex<double> Zxx = m_ZxxCalculated[freqIDThisPEInSta];
	const std::complex<double> Zxy = m_ZxyCalculated[freqIDThisPEInSta];
	const std::complex<double> Zyx = m_ZyxCalculated[freqIDThisPEInSta];
	const std::complex<double> Zyy = m_ZyyCalculated[freqIDThisPEInSta];
	const double eps = std::max(std::max(std::norm(Zxx), std::norm(Zxy)), std::max(std::norm(Zyx), std::norm(Zyy))) * 1.0e-20;

	const bool useImpedanceTensorInsteadOfPhase = ( AnalysisControl::getInstance()->getApparentResistivityAndPhaseTreatmentOption() == AnalysisControl::USE_Z_IF_SIGN_OF_RE_Z_DIFFER );

	const double ln10 = log(10.0);
	const long long nBlkNotFixed = static_cast<long long>( ( ResistivityBlock::getInstance() )->getNumResistivityBlockNotFixed() );
	for( long long imdl = 0; imdl < nBlkNotFixed; ++imdl ){

		const double dZxxRe = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + imdl ];
		const double dZxyRe = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + imdl ];
		const double dZyxRe = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + imdl ];
		const double dZyyRe = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + imdl ];
		const double dZxxIm = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + imdl ];
		const double dZxyIm = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + imdl ];
		const double dZyxIm = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + imdl ];
		const double dZyyIm = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + imdl ];

		// Zero clear
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + imdl ] = 0.0;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + imdl ] = 0.0; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + imdl ] = 0.0; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + imdl ] = 0.0;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + imdl ] = 0.0;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + imdl ] = 0.0;
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + imdl ] = 0.0; 
		sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + imdl ] = 0.0; 

		// Sensitivity is kept to be zero if corresponding error is negative
		if( m_apparentResistivityXXSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + imdl ] = dZxxRe / m_ZxxSD[freqIDGlobalInSta].realPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + imdl ] =
					2.0 / ln10 / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxRe + Zxx.imag() * dZxxIm ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XX );
			}
		}
		if( m_apparentResistivityXYSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + imdl ] = dZxyRe / m_ZxySD[freqIDGlobalInSta].realPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + imdl ] = 
					2.0 / ln10 / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyRe + Zxy.imag() * dZxyIm ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XY );
			}
		}
		if( m_apparentResistivityYXSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + imdl ] = dZyxRe / m_ZyxSD[freqIDGlobalInSta].realPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + imdl ] = 
					2.0 / ln10 / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxRe + Zyx.imag() * dZyxIm ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YX );
			}
		}
		if( m_apparentResistivityYYSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + imdl ] = dZyyRe / m_ZyySD[freqIDGlobalInSta].realPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + imdl ] =
					2.0 / ln10 / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyRe + Zyy.imag() * dZyyIm ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YY );
			}
		}
		if( m_PhaseXXSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + imdl ] = dZxxIm / m_ZxxSD[freqIDGlobalInSta].imagPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + imdl ] =
					CommonParameters::rad2deg / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxIm - Zxx.imag() * dZxxRe ) / m_PhaseXXSD[freqIDGlobalInSta];
			}
		}
		if( m_PhaseXYSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + imdl ] = dZxyIm / m_ZxySD[freqIDGlobalInSta].imagPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + imdl ] =
					CommonParameters::rad2deg / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyIm - Zxy.imag() * dZxyRe ) / m_PhaseXYSD[freqIDGlobalInSta];
			}
		}
		if( m_PhaseYXSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + imdl ] = dZyxIm / m_ZyxSD[freqIDGlobalInSta].imagPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + imdl ] = 
					CommonParameters::rad2deg / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxIm - Zyx.imag() * dZyxRe ) / m_PhaseYXSD[freqIDGlobalInSta];
			}
		}
		if( m_PhaseYYSD[freqIDGlobalInSta] > 0.0 ){
			if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + imdl ] = dZyyIm / m_ZyySD[freqIDGlobalInSta].imagPart;
			}else{
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + imdl ] = 
					CommonParameters::rad2deg / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyIm - Zyy.imag() * dZyyRe ) / m_PhaseYYSD[freqIDGlobalInSta];
			}
		}
	}

	if( !doesFixDistortionMatrix() ){ // Distortion matrix is not fixed

		if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
			
			assert( COMPONENT_ID_CXX + COMPONENT_ID_CXY + COMPONENT_ID_CYX + COMPONENT_ID_CYY == 6 ); 

			double dZxxRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxxIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxyRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxyIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyxRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyxIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyyRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyyIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			for( int i = 0; i < 4; ++i ){
				const long long ID = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[i]);
				if( ID < 0 ){
					continue;
				}
				dZxxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZxyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
			}

			for( int i = 0; i < 4; ++i ){
				const long long ID = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[i]);
				if( ID < 0 ){
					continue;
				}
				// Zero clear 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				// Sensitivity is kept to be zero if corresponding error is negative
				if( m_apparentResistivityXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxRe[i] / m_ZxxSD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxRe[i] + Zxx.imag() * dZxxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XX );
					}
				}
				if( m_apparentResistivityXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyRe[i] / m_ZxySD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyRe[i] + Zxy.imag() * dZxyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XY );
					}
				}
				if( m_apparentResistivityYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxRe[i] / m_ZyxSD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxRe[i] + Zyx.imag() * dZyxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YX );
					}
				}
				if( m_apparentResistivityYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyRe[i] / m_ZyySD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyRe[i] + Zyy.imag() * dZyyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YY );
					}
				}
				if( m_PhaseXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxIm[i] / m_ZxxSD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxIm[i] - Zxx.imag() * dZxxRe[i] ) / m_PhaseXXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyIm[i] / m_ZxySD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyIm[i] - Zxy.imag() * dZxyRe[i] ) / m_PhaseXYSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxIm[i] / m_ZyxSD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] =
							CommonParameters::rad2deg / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxIm[i] - Zyx.imag() * dZyxRe[i] ) / m_PhaseYXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyIm[i] / m_ZyySD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyIm[i] - Zyy.imag() * dZyyRe[i] ) / m_PhaseYYSD[freqIDGlobalInSta];
					}
				}
			}

		}
		else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
			
			assert( EX_GAIN + EY_GAIN + EX_ROTATION + EY_ROTATION == 6 ); 

			double dZxxRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxxIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxyRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZxyIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyxRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyxIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyyRe[4] = { -1.0, -1.0, -1.0, -1.0 };
			double dZyyIm[4] = { -1.0, -1.0, -1.0, -1.0 };
			for( int i = 0; i < 4; ++i ){
				const long long ID = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[i]);
				if( ID < 0 ){
					continue;
				}
				dZxxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZxyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
			}

			for( int i = 0; i < 4; ++i ){
				const long long ID = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[i]);
				if( ID < 0 ){
					continue;
				}
				// Zero clear 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				// Sensitivity is kept to be zero if corresponding error is negative
				if( m_apparentResistivityXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxRe[i] / m_ZxxSD[freqIDGlobalInSta].realPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxRe[i] + Zxx.imag() * dZxxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XX );
					}
				}
				if( m_apparentResistivityXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyRe[i] / m_ZxySD[freqIDGlobalInSta].realPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyRe[i] + Zxy.imag() * dZxyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XY );
					}
				}
				if( m_apparentResistivityYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxRe[i] / m_ZyxSD[freqIDGlobalInSta].realPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxRe[i] + Zyx.imag() * dZyxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YX );
					}
				}
				if( m_apparentResistivityYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyRe[i] / m_ZyySD[freqIDGlobalInSta].realPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyRe[i] + Zyy.imag() * dZyyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YY );
					}
				}
				if( m_PhaseXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxIm[i] / m_ZxxSD[freqIDGlobalInSta].imagPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxIm[i] - Zxx.imag() * dZxxRe[i] ) / m_PhaseXXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyIm[i] / m_ZxySD[freqIDGlobalInSta].imagPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyIm[i] - Zxy.imag() * dZxyRe[i] ) / m_PhaseXYSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxIm[i] / m_ZyxSD[freqIDGlobalInSta].imagPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] =
							CommonParameters::rad2deg / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxIm[i] - Zyx.imag() * dZyxRe[i] ) / m_PhaseYXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyIm[i] / m_ZyySD[freqIDGlobalInSta].imagPart; 
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyIm[i] - Zyy.imag() * dZyyRe[i] ) / m_PhaseYYSD[freqIDGlobalInSta];
					}
				}
			}
		}
		else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
			
			assert( EX_GAIN + EY_GAIN == 1 ); 

			double dZxxRe[2] = { -1.0, -1.0 };
			double dZxxIm[2] = { -1.0, -1.0 };
			double dZxyRe[2] = { -1.0, -1.0 };
			double dZxyIm[2] = { -1.0, -1.0 };
			double dZyxRe[2] = { -1.0, -1.0 };
			double dZyxIm[2] = { -1.0, -1.0 };
			double dZyyRe[2] = { -1.0, -1.0 };
			double dZyyIm[2] = { -1.0, -1.0 };
			for( int i = 0; i < 2; ++i ){
				const long long ID = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[i]);
				if( ID < 0 ){
					continue;
				}
				dZxxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZxyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZxyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyxRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyxIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
				dZyyRe[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID ];
				dZyyIm[i] = sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID ];
			}

			for( int i = 0; i < 2; ++i ){
				const long long ID = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[i]);
				if( ID < 0 ){
					continue;
				}
				// Zero clear 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0; 
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 0.0;
				// Sensitivity is kept to be zero if corresponding error is negative
				if( m_apparentResistivityXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxRe[i] / m_ZxxSD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxRe[i] + Zxx.imag() * dZxxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XX );
					}
				}
				if( m_apparentResistivityXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyRe[i] / m_ZxySD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyRe[i] + Zxy.imag() * dZxyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::XY );
					}
				}
				if( m_apparentResistivityYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxRe[i] / m_ZyxSD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxRe[i] + Zyx.imag() * dZyxIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YX );
					}
				}
				if( m_apparentResistivityYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyRe[i] / m_ZyySD[freqIDGlobalInSta].realPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfApparentResistivityYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							2.0 / ln10 / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyRe[i] + Zyy.imag() * dZyyIm[i] ) / calcLog10ErrorOfApparentResistivity( freqIDGlobalInSta, ObservedDataStationMT::YY );
					}
				}
				if( m_PhaseXXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxxIm[i] / m_ZxxSD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxx), eps ) * ( Zxx.real() * dZxxIm[i] - Zxx.imag() * dZxxRe[i] ) / m_PhaseXXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseXYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZxyIm[i] / m_ZxySD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseXY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zxy), eps ) * ( Zxy.real() * dZxyIm[i] - Zxy.imag() * dZxyRe[i] ) / m_PhaseXYSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYXSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyxIm[i] / m_ZyxSD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYX[freqIDThisPEInSta]) + nBlkNotFixed + ID ] =
							CommonParameters::rad2deg / std::max( std::norm(Zyx), eps ) * ( Zyx.real() * dZyxIm[i] - Zyx.imag() * dZyxRe[i] ) / m_PhaseYXSD[freqIDGlobalInSta];
					}
				}
				if( m_PhaseYYSD[freqIDGlobalInSta] > 0.0 ){
					if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = dZyyIm[i] / m_ZyySD[freqIDGlobalInSta].imagPart;
					}else{
						sensitivityMatrix[ static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfPhaseYY[freqIDThisPEInSta]) + nBlkNotFixed + ID ] = 
							CommonParameters::rad2deg / std::max( std::norm(Zyy), eps ) * ( Zyy.real() * dZyyIm[i] - Zyy.imag() * dZyyRe[i] ) / m_PhaseYYSD[freqIDGlobalInSta];
					}
				}
			}
		}
	}

}

// Calculate data vector of this PE
void ObservedDataStationApparentResistivityAndPhase::calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const{

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE( freq );

	if( freqIDThisPEInSta < 0 ){// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[ freqIDThisPEInSta ];
	const bool useImpedanceTensorInsteadOfPhase = ( AnalysisControl::getInstance()->getApparentResistivityAndPhaseTreatmentOption() == AnalysisControl::USE_Z_IF_SIGN_OF_RE_Z_DIFFER );

	if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XX ) ){
		OutputFiles::m_logFile << "Notice : Zxx is used instead of apparent resistivity and phase (Station ID : " << getStationID() << ", Frequency : " << freq << " [Hz])." << std::endl;
		vector[ offset + m_dataIDOfApparentResistivityXX[freqIDThisPEInSta] ] = m_ZxxResidual[freqIDThisPEInSta].realPart;
		vector[ offset + m_dataIDOfPhaseXX[freqIDThisPEInSta] ] = m_ZxxResidual[freqIDThisPEInSta].imagPart;
	}else{
		vector[ offset + m_dataIDOfApparentResistivityXX[freqIDThisPEInSta] ] = m_apparentResistivityXXResidual[freqIDThisPEInSta];
		vector[ offset + m_dataIDOfPhaseXX[freqIDThisPEInSta] ] = m_PhaseXXResidual[freqIDThisPEInSta];
	}
	if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::XY ) ){
		OutputFiles::m_logFile << "Notice : Zxy is used instead of apparent resistivity and phase (Station ID : " << getStationID() << ", Frequency : " << freq << " [Hz])." << std::endl;
		vector[ offset + m_dataIDOfApparentResistivityXY[freqIDThisPEInSta] ] = m_ZxyResidual[freqIDThisPEInSta].realPart;
		vector[ offset + m_dataIDOfPhaseXY[freqIDThisPEInSta] ] = m_ZxyResidual[freqIDThisPEInSta].imagPart;
	}else{
		vector[ offset + m_dataIDOfApparentResistivityXY[freqIDThisPEInSta] ] = m_apparentResistivityXYResidual[freqIDThisPEInSta];
		vector[ offset + m_dataIDOfPhaseXY[freqIDThisPEInSta] ] = m_PhaseXYResidual[freqIDThisPEInSta];
	}
	if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YX ) ){
		OutputFiles::m_logFile << "Notice : Zyx is used instead of apparent resistivity and phase (Station ID : " << getStationID() << ", Frequency : " << freq << " [Hz])." << std::endl;
		vector[ offset + m_dataIDOfApparentResistivityYX[freqIDThisPEInSta] ] = m_ZyxResidual[freqIDThisPEInSta].realPart;
		vector[ offset + m_dataIDOfPhaseYX[freqIDThisPEInSta] ] = m_ZyxResidual[freqIDThisPEInSta].imagPart;
	}else{
		vector[ offset + m_dataIDOfApparentResistivityYX[freqIDThisPEInSta] ] = m_apparentResistivityYXResidual[freqIDThisPEInSta];
		vector[ offset + m_dataIDOfPhaseYX[freqIDThisPEInSta] ] = m_PhaseYXResidual[freqIDThisPEInSta];
	}
	if( useImpedanceTensorInsteadOfPhase && isUsedImpedanceTensorFromFreqIDs( freqIDThisPEInSta, ObservedDataStationMT::YY ) ){
		OutputFiles::m_logFile << "Notice : Zyy is used instead of apparent resistivity and phase (Station ID : " << getStationID() << ", Frequency : " << freq << " [Hz])." << std::endl;
		vector[ offset + m_dataIDOfApparentResistivityYY[freqIDThisPEInSta] ] = m_ZyyResidual[freqIDThisPEInSta].realPart;
		vector[ offset + m_dataIDOfPhaseYY[freqIDThisPEInSta] ] = m_ZyyResidual[freqIDThisPEInSta].imagPart;
	}else{
		vector[ offset + m_dataIDOfApparentResistivityYY[freqIDThisPEInSta] ] = m_apparentResistivityYYResidual[freqIDThisPEInSta];
		vector[ offset + m_dataIDOfPhaseYY[freqIDThisPEInSta] ] = m_PhaseYYResidual[freqIDThisPEInSta];
	}
}

// Calulate L2 norm of misfit
double ObservedDataStationApparentResistivityAndPhase::calculateErrorSumOfSquaresThisPE() const{

	double misfit(0.0);

	for( int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq ){
		misfit += pow( m_apparentResistivityXXResidual[ifreq], 2 );
		misfit += pow( m_apparentResistivityXYResidual[ifreq], 2 );
		misfit += pow( m_apparentResistivityYXResidual[ifreq], 2 );
		misfit += pow( m_apparentResistivityYYResidual[ifreq], 2 );
		misfit += pow( m_PhaseXXResidual[ifreq], 2 );
		misfit += pow( m_PhaseXYResidual[ifreq], 2 );
		misfit += pow( m_PhaseYXResidual[ifreq], 2 );
		misfit += pow( m_PhaseYYResidual[ifreq], 2 );
	}

	return misfit;

}

// Apparent resistivity is ignored for all frequencies
bool ObservedDataStationApparentResistivityAndPhase::isApparentResistivityIgnoredForAllFrequencies( const int iComp ) const{

	switch (iComp)
	{
	case ObservedDataStationMT::XX:
		for( int iFreq = 0; iFreq <	m_numOfFrequency; ++iFreq ){
			if( m_apparentResistivityXXSD[iFreq] >= 0.0 ){
				return false;
			}
		}
		break;
	case ObservedDataStationMT::XY:
		for( int iFreq = 0; iFreq <	m_numOfFrequency; ++iFreq ){
			if( m_apparentResistivityXYSD[iFreq] >= 0.0 ){
				return false;
			}
		}
		break;
	case ObservedDataStationMT::YX:
		for( int iFreq = 0; iFreq <	m_numOfFrequency; ++iFreq ){
			if( m_apparentResistivityYXSD[iFreq] >= 0.0 ){
				return false;
			}
		}
		break;
	case ObservedDataStationMT::YY:
		for( int iFreq = 0; iFreq <	m_numOfFrequency; ++iFreq ){
			if( m_apparentResistivityYYSD[iFreq] >= 0.0 ){
				return false;
			}
		}
		break;
	default:
		OutputFiles::m_logFile << "Error : Type of index of apparent resistivity component is wrong. : " << iComp << std::endl;
		exit(1);
		break;
	}

	return true;

}

// Calculate log 10 error of apparent resistivity
double ObservedDataStationApparentResistivityAndPhase::calcLog10ErrorOfApparentResistivity( const int freqIDGlobalInSta, const int iComp ) const{

	switch (iComp)
	{
	case ObservedDataStationMT::XX:
		return m_apparentResistivityXXSD[freqIDGlobalInSta] / m_apparentResistivityXXObserved[freqIDGlobalInSta] / log(10);
	case ObservedDataStationMT::XY:
		return m_apparentResistivityXYSD[freqIDGlobalInSta] / m_apparentResistivityXYObserved[freqIDGlobalInSta] / log(10);
	case ObservedDataStationMT::YX:
		return m_apparentResistivityYXSD[freqIDGlobalInSta] / m_apparentResistivityYXObserved[freqIDGlobalInSta] / log(10);
	case ObservedDataStationMT::YY:
		return m_apparentResistivityYYSD[freqIDGlobalInSta] / m_apparentResistivityYYObserved[freqIDGlobalInSta] / log(10);
	default:
		OutputFiles::m_logFile << "Error : Type of index of apparent resistivity component is wrong. : " << iComp << std::endl;
		exit(1);
		break;
	}

	double dummy(-1);	
	return dummy;

}

// Determine whether to use impedance tensor instead of phase by checking sign differece of real part of impedance tensor 
bool ObservedDataStationApparentResistivityAndPhase::isUsedImpedanceTensor( const double phaseObs, const double phaseError, const double phaseCalc ) const{

	const double phaseObsRad = phaseObs * CommonParameters::deg2rad; 
	const double phaseCalcRad = phaseCalc * CommonParameters::deg2rad;
	if( cos(phaseObsRad) * cos(phaseCalcRad) < 0 ){
		return true;
	}else{
		return false;
	}

}

// Determine whether to use impedance tensor instead of phase by checking sign differece of real part of impedance tensor from frequeny IDs
bool ObservedDataStationApparentResistivityAndPhase::isUsedImpedanceTensorFromFreqIDs( const int freqIDThisPEInSta, const int iComp ) const{

	assert( freqIDThisPEInSta >= 0 );
	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[freqIDThisPEInSta];

	switch (iComp)
	{
	case ObservedDataStationMT::XX:
		return isUsedImpedanceTensor( m_PhaseXXObserved[freqIDGlobalInSta], m_PhaseXXSD[freqIDGlobalInSta], m_PhaseXXCalculated[freqIDThisPEInSta] );
	case ObservedDataStationMT::XY:
		return isUsedImpedanceTensor( m_PhaseXYObserved[freqIDGlobalInSta], m_PhaseXYSD[freqIDGlobalInSta], m_PhaseXYCalculated[freqIDThisPEInSta] );
	case ObservedDataStationMT::YX:
		return isUsedImpedanceTensor( m_PhaseYXObserved[freqIDGlobalInSta], m_PhaseYXSD[freqIDGlobalInSta], m_PhaseYXCalculated[freqIDThisPEInSta] );
	case ObservedDataStationMT::YY:
		return isUsedImpedanceTensor( m_PhaseYYObserved[freqIDGlobalInSta], m_PhaseYYSD[freqIDGlobalInSta], m_PhaseYYCalculated[freqIDThisPEInSta] );
	default:
		OutputFiles::m_logFile << "Error : Type of index of apparent resistivity component is wrong. : " << iComp << std::endl;
		exit(1);
		break;
	}

	return false;

}
