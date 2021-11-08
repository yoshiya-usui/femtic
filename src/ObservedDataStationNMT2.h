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
#ifndef DBLDEF_OBSERVED_DATA_NMT2
#define DBLDEF_OBSERVED_DATA_NMT2

#include <vector>
#include <complex>

#include "ObservedDataStation.h"
#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of NMT station ( triangle area )
class ObservedDataStationNMT2: public ObservedDataStation{
	public:
		// Constructer
		explicit ObservedDataStationNMT2();

		// Destructer
		virtual ~ObservedDataStationNMT2();

		// Read data from input file
		virtual void inputObservedData( std::ifstream& inFile );

		// Find elements including dipole
		void findElementsIncludingDipoles();

		// Calulate difference of voltages
		void calculateVoltageDifferences( const Forward3D* const ptrForward3D, const int rhsVectorIDOfVoltageDifference1st, const int rhsVectorIDOfVoltageDifference2nd );

		// Calulate Impedance tensor
		void calculateImpedanceTensor( const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount );

		// Initialize difference of voltages
		void initializeVoltageDifferences( const int iPol );

		// Initialize Impedance tensor and errors
		void initializeImpedanceTensorsAndErrors();

		// Allocate memory for the calculated impedance tensors and errors
		virtual void allocateMemoryForCalculatedValues();

		// Output calculated values of impedance tensors
		virtual void outputCalculatedValues() const;
			
		// Calulate interpolator vector of voltage difference
		void calcInterpolatorVectorOfVoltageDifference( Forward3D* const ptrForward3D );

		// Calulate sensitivity matrix of Impedance tensors
		virtual void calculateSensitivityMatrix( const double freq, const int nModel,
			const ObservedDataStationPoint* const ptrStationOfMagneticField,
			const std::complex<double>* const derivativesOfEMFieldExPol,
			const std::complex<double>* const derivativesOfEMFieldEyPol,
			double* const sensitivityMatrix, const bool forceSDToOne = false ) const;

		// Calculate data vector of this PE
		virtual void calculateResidualVectorOfDataThisPE( const double freq, const int offset, double* vector ) const;

		// Calulate sum of square of misfit
		virtual double calculateErrorSumOfSquaresThisPE() const;

		// Get location of the station
		const CommonParameters::locationDipole& getLocationOfStation( const int iDipole ) const;

		// Get Z coordinate of the point
		double getZCoordOfPoint( const int iDipole , const int num ) const;

	protected:
		enum ImpedanceTensorComponentNMT2{
			XX = 0,
			XY,
			YX,
			YY
		};

		// Location of the station
		CommonParameters::locationDipole m_location[2];

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

	private:
		std::complex<double> m_voltageCalculated[2][2];

		// Total number of elemenst including dipole
		int m_numElementsIncludingDipole[2];

		// Elements including dipoles
		int* m_elementsIncludingDipole[2];

		// Faces including dipoles
		int* m_facesIncludingDipole[2];

		// local coordinate values of the location of start point of line segment
		CommonParameters::locationXY* m_localCoordinateValuesStartPoint[2];

		// local coordinate values of the location of end point of line segment
		CommonParameters::locationXY* m_localCoordinateValuesEndPoint[2];

		//// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		//CommonParameters::VolumeCoords* m_volumeCoordinateValuesStartPoint[2];

		//// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		//CommonParameters::VolumeCoords* m_volumeCoordinateValuesEndPoint[2];

		// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		CommonParameters::AreaCoords* m_areaCoordinateValuesStartPoint[2];

		// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		CommonParameters::AreaCoords* m_areaCoordinateValuesEndPoint[2];

		//int m_columnNumberOfVoltageDifferenceInRhsMatrix[2];
		int m_rhsVectorIDOfVoltageDifference[2];

		// Copy constructer
		ObservedDataStationNMT2(const ObservedDataStationNMT2& rhs);

		// Copy assignment operator
		ObservedDataStationNMT2& operator=(const ObservedDataStationNMT2& rhs);
};

#endif
