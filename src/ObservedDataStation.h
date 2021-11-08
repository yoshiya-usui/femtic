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
#ifndef DBLDEF_OBSERVED_DATA_STATION
#define DBLDEF_OBSERVED_DATA_STATION

#include <vector>
#include <complex>

#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of each station
class ObservedDataStation{
	public:
		// Constructer
		explicit ObservedDataStation();

		// Destructer
		~ObservedDataStation();

		// Get ID of station
		int getStationID() const;

		// Get ID of the station where magnetic field is observed
		int getIDOfMagneticFieldStation() const;

		// Get total number of frequencies
		int getTotalNumberOfFrequency() const;

		// Get frequencies at which observed value exists
		double getFrequencyValues( const int num ) const;

		// Find and return frequency IDs ( consecutive number in this station ) from frequency value.
		// If the specified frequency is not included, return -1.
		int getFreqIDs( const double freq ) const;

		// Set up frequencies calculated by this PE, at which observed value exists
		void setupFrequenciesCalculatedByThisPE( const int nFreqCalculatedByThisPE, const double* freqCalculatedByThisPE );

		// Find and return frequency IDs among the ones calculated by this PE ( consecutive number in this PE and station ) from frequency value.
		// If the specified frequency is not included, return -1.
		int getFreqIDsAmongThisPE( const double freq ) const;

	protected:
		// ID of station
		int m_stationID;

		// ID of the station where magnetic field is observed
		int m_IDOfMagneticFieldStation;

		// Total number of frequencies
		int m_numOfFrequency;

		// Frequencies at which observed value exists
		double* m_freq;

		// Total number of frequencies calculated by this PE
		int m_numOfFreqCalculatedByThisStaAndPE;

		// IDs of frequencies calculated by this PE, at which observed value exists
		std::vector<int> m_freqIDsAmongThisStationCalculatedByThisPE;

	private:
		// Copy constructer
		ObservedDataStation(const ObservedDataStation& rhs);

		// Copy assignment operator
		ObservedDataStation& operator=(const ObservedDataStation& rhs);

};

#endif
