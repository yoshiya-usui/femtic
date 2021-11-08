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
#ifndef DBLDEF_OBSERVED_DATA_HTF
#define DBLDEF_OBSERVED_DATA_HTF

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of HTF station
class ObservedDataStationHTF: public ObservedDataStationPoint{
	public:
		// Constructer
		explicit ObservedDataStationHTF();

		// Destructer
		~ObservedDataStationHTF();
			
		// Read data from input file
		void inputObservedData( std::ifstream& inFile );

		// Calulate horizontal magnetic field transfer function
		void calculateHTF( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount );

		// Initialize horizontal magnetic field transfer functions and errors
		void initializeHTFsAndErrors();

		// Allocate memory for the calculated values of horizontal magnetic field transfer functions and errors
		void allocateMemoryForCalculatedValues();

		// Output calculated values of horizontal magnetic field transfer functions
		void outputCalculatedValues() const;

		// Calulate sensitivity matrix of HTF
		void calculateSensitivityMatrix( const double freq, const int nModel,
			const ObservedDataStationPoint* const ptrStationOfMagneticField,
			const std::complex<double>* const derivativesOfEMFieldExPol,
			const std::complex<double>* const derivativesOfEMFieldEyPol,
			double* const sensitivityMatrix ) const;
	
		// Calculate data vector of this PE
		void calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const;

		// Calulate sum of square of misfit
		double calculateErrorSumOfSquaresThisPE() const;

	private:
		std::complex<double>* m_TxxObserved;
		std::complex<double>* m_TxyObserved;
		std::complex<double>* m_TyxObserved;
		std::complex<double>* m_TyyObserved;

		CommonParameters::DoubleComplexValues* m_TxxSD;
		CommonParameters::DoubleComplexValues* m_TxySD;
		CommonParameters::DoubleComplexValues* m_TyxSD;
		CommonParameters::DoubleComplexValues* m_TyySD;

		std::complex<double>* m_TxxCalculated;
		std::complex<double>* m_TxyCalculated;
		std::complex<double>* m_TyxCalculated;
		std::complex<double>* m_TyyCalculated;

		CommonParameters::DoubleComplexValues* m_TxxResidual;
		CommonParameters::DoubleComplexValues* m_TxyResidual;
		CommonParameters::DoubleComplexValues* m_TyxResidual;
		CommonParameters::DoubleComplexValues* m_TyyResidual;

		CommonParameters::InitComplexValues* m_dataIDOfTxx;
		CommonParameters::InitComplexValues* m_dataIDOfTxy;
		CommonParameters::InitComplexValues* m_dataIDOfTyx;
		CommonParameters::InitComplexValues* m_dataIDOfTyy;

		// Copy constructer
		ObservedDataStationHTF(const ObservedDataStationHTF& rhs);

		// Copy assignment operator
		ObservedDataStationHTF& operator=(const ObservedDataStationHTF& rhs);

};

#endif
