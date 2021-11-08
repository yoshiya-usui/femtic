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
#ifndef DBLDEF_OBSERVED_DATA_MT
#define DBLDEF_OBSERVED_DATA_MT

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataStationMT: public ObservedDataStationPoint{
	public:
		enum ImpedanceTensorComponent{
			XX = 0,
			XY,
			YX,
			YY
		};

		enum ComponentIDOfDistortionMatrix{
			COMPONENT_ID_CXX = 0,
			COMPONENT_ID_CXY,
			COMPONENT_ID_CYX,
			COMPONENT_ID_CYY
		};

		enum ComponentIDOfGainsAndRotations{
			EX_GAIN = 0,
			EY_GAIN,
			EX_ROTATION,
			EY_ROTATION
		};

		// Constructer
		explicit ObservedDataStationMT();

		// Destructer
		virtual ~ObservedDataStationMT();

		// Read data from input file
		virtual void inputObservedData( std::ifstream& inFile );

		// Calulate electric field
		void calculateElectricField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy );

		// Calulate Impedance tensor
		void calculateImpedanceTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount );

		// Initialize electric field
		void initializeElectricField( const int iPol );

		// Initialize Impedance tensor and errors
		void initializeImpedanceTensorsAndErrors();

		// Allocate memory for the calculated Impedance tensors and errors
		virtual void allocateMemoryForCalculatedValues();

		// Output calculated Impedance tensors
		virtual void outputCalculatedValues() const;

		// Calulate interpolator vector of electric field
		void calcInterpolatorVectorOfElectricField( Forward3D* const ptrForward3D );

		// Calulate sensitivity matrix
		virtual void calculateSensitivityMatrix( const double freq, const int nModel,
			const ObservedDataStationPoint* const ptrStationOfMagneticField,
			const std::complex<double>* const derivativesOfEMFieldExPol,
			const std::complex<double>* const derivativesOfEMFieldEyPol,
			double* const sensitivityMatrix, const bool forceSDToOne = false ) const;

		// Calculate data vector of this PE
		virtual void calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const;

		// Calulate sum of square of misfit
		virtual double calculateErrorSumOfSquaresThisPE() const;

		// Copy current distortion parameters to previous ones
		void copyDistortionParamsCurToPre( const int iComp );

		// Get flag specifing whether distortion matrix are fixed or not
		bool doesFixDistortionMatrix() const;

		// Get type of the electric field used to calculate response functions
		int getTypeOfElectricField() const;

		// Set flag specifing whether distortion matrix are fixed or not
		void setFixDistortionMatrix( const bool doesFix );

		// Set type of the electric field used to calculate response functions
		void setTypeOfElectricField( const int type );

		// Set distortion parameters of previous iteration
		void setDistortionParamsPre( const int iComp, const double val );

		// Set distortion parameters 
		void setDistortionParams( const int iComp, const double val );

		// Set ID of distortion parameters
		void setIDOfDistortionParams( const int iComp, const int ID );
			
		// Set full updated value of distortion parameters 
		void setDistortionParamsUpdatedFull( const int iComp, const double val );

		// Update distortion parameters
		void updateDistortionParams( const double dampingFactor );

		// Get distortion parameters of previous iteration
		double getDistortionParamsPre( const int iComp ) const;

		// Get distortion parameters
		double getDistortionParams( const int iComp ) const;

		// Get ID of distortion parameters
		int getIDOfDistortionParams( const int iComp ) const;

		// Get full updated value of distortion parameters
		double getDistortionParamsUpdatedFull( const int iComp ) const;

	protected:
		struct DistortionMatrixDifferences{
			// Previous distortion matrix difference
			double distortionMatrixDifferencePre[4];

			// Distortion matrix difference
			double distortionMatrixDifference[4];

			// Distorsion matrix difference obtained by inversion which is the ones fully updated ( damping factor = 1 )
			double distortionMatrixDifferenceUpdatedFull[4];

			// Component ID of distortion matrix difference. This ID is -1 if distortion matrix is fixed or this station has no data
			int IDsOfDistortionMatrixDifference[4];
		};

		struct GainsAndRotations{
			// Previous gains and rotations of distortion matrix
			double gainsAndRotationsPre[4];

			// Gains and rotations of distortion matrix
			double gainsAndRotations[4];

			// Gains and rotations of distortion matrix obtained by inversion which is the ones fully updated ( damping factor = 1 )
			double gainsAndRotationsUpdatedFull[4];

			// Component ID of gains and rotations of distortion matrix. This ID is -1 if distortion matrix is fixed or this station has no data
			int IDsOfGainsAndRotations[4];
		};

		std::complex<double>* m_ZxxObserved;
		std::complex<double>* m_ZxyObserved;
		std::complex<double>* m_ZyxObserved;
		std::complex<double>* m_ZyyObserved;

		CommonParameters::DoubleComplexValues* m_ZxxSD;
		CommonParameters::DoubleComplexValues* m_ZxySD;
		CommonParameters::DoubleComplexValues* m_ZyxSD;
		CommonParameters::DoubleComplexValues* m_ZyySD;

		std::complex<double>* m_ZxxCalculated;
		std::complex<double>* m_ZxyCalculated;
		std::complex<double>* m_ZyxCalculated;
		std::complex<double>* m_ZyyCalculated;

		CommonParameters::DoubleComplexValues* m_ZxxResidual;
		CommonParameters::DoubleComplexValues* m_ZxyResidual;
		CommonParameters::DoubleComplexValues* m_ZyxResidual;
		CommonParameters::DoubleComplexValues* m_ZyyResidual;

		CommonParameters::InitComplexValues* m_dataIDOfZxx;
		CommonParameters::InitComplexValues* m_dataIDOfZxy;
		CommonParameters::InitComplexValues* m_dataIDOfZyx;
		CommonParameters::InitComplexValues* m_dataIDOfZyy;

		// Arrays of distortion matrix differences
		DistortionMatrixDifferences* m_arrayDistortionMatrixDifferences;

		// Arrays of gains and rotations of distortion matrix
		GainsAndRotations* m_arrayGainsAndRotations;

	private:
		std::complex<double> m_ExCalculated[2];
		std::complex<double> m_EyCalculated[2];

		int m_rhsVectorIDOfEx;
		int m_rhsVectorIDOfEy;

		// Flag specifing fix distortion matrix or not
		bool m_fixDistortionMatrix;

		// Type of the electric field used to calculated response functions
		int m_typeOfElectricField;

		// Copy constructer
		ObservedDataStationMT(const ObservedDataStationMT& rhs);

		// Copy assignment operator
		ObservedDataStationMT& operator=(const ObservedDataStationMT& rhs);

};

#endif
