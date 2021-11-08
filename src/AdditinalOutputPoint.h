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
#ifndef DBLDEF_ADDITIONAL_OUTPUT_POINT
#define DBLDEF_ADDITIONAL_OUTPUT_POINT

#include <vector>
#include <complex>

#include "Forward3D.h"
#include "CommonParameters.h"

// Class of additional output points
class AdditinalOutputPoint{
	public:
		// Constructer
		explicit AdditinalOutputPoint();

		// Destructer
		~AdditinalOutputPoint();

		// Read data from input file
		void inputObservedData( std::ifstream& inFile );

		// Find element including a point
		void findElementIncludingStation();

		// Calulate EM field
		void calculateEMField( const int ifreq, const Forward3D* const ptrForward3D );

		// Initialize EM field
		void initializeEMfield( const int nfreq );

		// Allocate memory for EM field
		void allocateMemoryForCalculatedValues( const int nfreq );

		// Output calculated values of EM field
		void outputCalculatedValues( const int nfreq, const double* const freq ) const;

	private:
		// Location of the station
		CommonParameters::locationXYZ m_location;

		// Element including the station
		int m_elementIncludingStation;

		// local coordinate values of the location of the station
		CommonParameters::locationXYZ m_localCoordinateValues;

		// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		CommonParameters::VolumeCoords m_volumeCoordinateValues;

		std::complex<double>* m_ExCalculated[2];
		std::complex<double>* m_EyCalculated[2];
		std::complex<double>* m_EzCalculated[2];
		std::complex<double>* m_HxCalculated[2];
		std::complex<double>* m_HyCalculated[2];
		std::complex<double>* m_HzCalculated[2];

		// Copy constructer
		AdditinalOutputPoint(const AdditinalOutputPoint& rhs);

		// Copy assignment operator
		AdditinalOutputPoint& operator=(const AdditinalOutputPoint& rhs);

};

#endif
