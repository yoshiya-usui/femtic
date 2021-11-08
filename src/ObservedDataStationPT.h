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
#ifndef DBLDEF_OBSERVED_DATA_PT
#define DBLDEF_OBSERVED_DATA_PT

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of PT station
class ObservedDataStationPT: public ObservedDataStationPoint{
	public:
		// Constructer
		explicit ObservedDataStationPT();

		// Destructer
		~ObservedDataStationPT();

		// Read data from input file
		void inputObservedData( std::ifstream& inFile );

		// Calulate electric field
		void calculateElectricField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy );

		// Calulate phase tensor
		void calculatePhaseTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount );

		// Initialize electric field
		void initializeElectricField( const int iPol );

		// Initialize phase tensors and errors
		void initializePhaseTensorsAndErrors();

		// Allocate memory for the calculated values of phase tensors and errors
		void allocateMemoryForCalculatedValues();

		// Output calculated values of phase tensors
		void outputCalculatedValues() const;
			
		// Calulate interpolator vector of electric field
		void calcInterpolatorVectorOfElectricField( Forward3D* const ptrForward3D );

		// Calulate sensitivity matrix of phase tensor
		void calculateSensitivityMatrix( const double freq, const int nModel,
			const ObservedDataStationPoint* const ptrStationOfMagneticField,
			const std::complex<double>* const derivativesOfEMFieldExPol,
			const std::complex<double>* const derivativesOfEMFieldEyPol,
			double* const sensitivityMatrix ) const;
	
		// Calculate data vector of this PE
		void calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const;

		// Calulate sum of square of misfit
		double calculateErrorSumOfSquaresThisPE() const;

		// Get type of the electric field used to calculate response functions
		int getTypeOfElectricField() const;

		// Set type of the electric field used to calculate response functions
		void setTypeOfElectricField( const int type );

	private:
		double* m_PTxxObserved;
		double* m_PTxyObserved;
		double* m_PTyxObserved;
		double* m_PTyyObserved;

		double* m_PTxxSD;
		double* m_PTxySD;
		double* m_PTyxSD;
		double* m_PTyySD;

		double* m_PTxxCalculated;
		double* m_PTxyCalculated;
		double* m_PTyxCalculated;
		double* m_PTyyCalculated;

		double* m_PTxxResidual;
		double* m_PTxyResidual;
		double* m_PTyxResidual;
		double* m_PTyyResidual;

		std::complex<double> m_ExCalculated[2];
		std::complex<double> m_EyCalculated[2];

		int m_rhsVectorIDOfEx;
		int m_rhsVectorIDOfEy;

		int* m_dataIDOfPTxx;
		int* m_dataIDOfPTxy;
		int* m_dataIDOfPTyx;
		int* m_dataIDOfPTyy;

		// Type of the electric field used to calculated response functions
		int m_typeOfElectricField;

		// Copy constructer
		ObservedDataStationPT(const ObservedDataStationPT& rhs);

		// Copy assignment operator
		ObservedDataStationPT& operator=(const ObservedDataStationPT& rhs);

};

#endif
