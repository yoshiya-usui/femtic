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
#ifndef DBLDEF_OBSERVED_DATA_VTF
#define DBLDEF_OBSERVED_DATA_VTF

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of VTF station
class ObservedDataStationVTF: public ObservedDataStationPoint{
	public:
		// Constructer
		explicit ObservedDataStationVTF();

		// Destructer
		~ObservedDataStationVTF();
			
		// Read data from input file
		void inputObservedData( std::ifstream& inFile );
			
		// Calulate vertical magnetic field
		void calculateVerticalMagneticField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfHz );

		// Calulate vertical magnetic field transfer function
		void calculateVTF( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount );

		// Initialize vertical magnetic field
		void initializeVerticalMagneticField( const int iPol );

		// Initialize vertical magnetic field transfer functions and errors
		void initializeVTFsAndErrors();

		// Allocate memory for the calculated values of vertical magnetic field transfer functions and errors
		void allocateMemoryForCalculatedValues();

		// Output calculated values of vertical magnetic field transfer functions
		void outputCalculatedValues() const;

		// Calulate interpolator vector of vertical magnetic field
		void calcInterpolatorVectorOfVerticalMagneticField( Forward3D* const ptrForward3D );

		// Calulate sensitivity matrix of VTF
		void calculateSensitivityMatrix( const double freq, const int nModel,
			const ObservedDataStationPoint* const ptrStationOfMagneticField,
			const std::complex<double>* const derivativesOfEMFieldExPol,
			const std::complex<double>* const derivativesOfEMFieldEyPol,
			double* const sensitivityMatrix ) const;
	
		// Calculate data vector of this PE
		void calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const;

		// Calulate sum of square of misfit
		double calculateErrorSumOfSquaresThisPE() const;

		// Get VTK
		bool getVTF( const double freq, std::complex<double>& Tzx, std::complex<double>& Tzy ) const;

	private:
		std::complex<double>* m_TzxObserved;
		std::complex<double>* m_TzyObserved;

		CommonParameters::DoubleComplexValues* m_TzxSD;
		CommonParameters::DoubleComplexValues* m_TzySD;

		std::complex<double>* m_TzxCalculated;
		std::complex<double>* m_TzyCalculated;

		CommonParameters::DoubleComplexValues* m_TzxResidual;
		CommonParameters::DoubleComplexValues* m_TzyResidual;

		std::complex<double> m_HzCalculated[2];

		//int m_columnNumberOfHzInRhsMatrix;
		int m_rhsVectorIDOfHz;

		//int* m_dataIDOfTzx;
		//int* m_dataIDOfTzy;
		CommonParameters::InitComplexValues* m_dataIDOfTzx;
		CommonParameters::InitComplexValues* m_dataIDOfTzy;

		// Copy constructer
		ObservedDataStationVTF(const ObservedDataStationVTF& rhs);

		// Copy assignment operator
		ObservedDataStationVTF& operator=(const ObservedDataStationVTF& rhs);

};

#endif
