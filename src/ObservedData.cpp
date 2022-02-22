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

#include "mpi.h"

#include "AnalysisControl.h"
#include "ObservedData.h"
#include "MeshDataBrickElement.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

#include "Util.h"

// Return the the instance of the class
ObservedData* ObservedData::getInstance(){
   	static ObservedData instance;// The only instance
  	return &instance;
}

// Constructer
ObservedData::ObservedData():
	m_numAllStations(0),
	m_numStationsMT(0),
	m_numStationsApparentResistivityAndPhase(0),
	m_numStationsHTF(0),
	m_numStationsVTF(0),
	m_numStationsPT(0),
	m_numStationsNMT(0),
	m_numStationsNMT2(0),
	m_numStationsNMT2ApparentResistivityAndPhase(0),
	m_numAdditinalOutputPoint(0),
	m_numInterpolatorVectors(0),
	m_numOfFrequenciesCalculatedByThisPE(0),
	m_IDsOfFrequenciesCalculatedByThisPE(NULL),
	m_valuesOfFrequenciesCalculatedByThisPE(NULL),
	m_observedStationMT(NULL),
	m_observedStationApparentResistivityAndPhase(NULL),
	m_observedStationHTF(NULL),
	m_observedStationVTF(NULL),
	m_observedStationPT(NULL),
	m_observedStationNMT(NULL),
	m_observedStationNMT2(NULL),
	m_observedStationNMT2ApparentResistivityAndPhase(NULL),
	m_additinalOutputPoint(NULL),
	m_numObservedDataThisPEAccumulated(NULL)
{
	//for( int i = 0; i < 2; ++i ){
	//	m_derivativesOfEMField[i] = NULL;
	//}
}

//Destructer
ObservedData::~ObservedData(){
	
	if( m_observedStationMT != NULL ){
		delete [] m_observedStationMT;
		m_observedStationMT = NULL;
	}

	if( m_observedStationApparentResistivityAndPhase != NULL ){
		delete [] m_observedStationApparentResistivityAndPhase;
		m_observedStationApparentResistivityAndPhase = NULL;
	}

	if( m_observedStationHTF != NULL ){
		delete [] m_observedStationHTF;
		m_observedStationHTF = NULL;
	}

	if( m_observedStationVTF != NULL ){
		delete [] m_observedStationVTF;
		m_observedStationVTF = NULL;
	}

	if( m_observedStationPT != NULL ){
		delete [] m_observedStationPT;
		m_observedStationPT = NULL;
	}

	if( m_observedStationNMT != NULL ){
		delete [] m_observedStationNMT;
		m_observedStationNMT = NULL;
	}

	if( m_observedStationNMT2 != NULL ){
		delete [] m_observedStationNMT2;
		m_observedStationNMT2 = NULL;
	}

	if( m_observedStationNMT2ApparentResistivityAndPhase != NULL ){
		delete [] m_observedStationNMT2ApparentResistivityAndPhase;
		m_observedStationNMT2ApparentResistivityAndPhase = NULL;
	}

	if( m_additinalOutputPoint != NULL ){
		delete [] m_additinalOutputPoint;
		m_additinalOutputPoint = NULL;
	}

	if( m_numObservedDataThisPEAccumulated != NULL ){
		delete [] m_numObservedDataThisPEAccumulated;
		m_numObservedDataThisPEAccumulated = NULL;
	}

	//for( int i = 0; i < 2; ++i ){
	//	if( m_derivativesOfEMField[i] != NULL ){
	//		delete [] m_derivativesOfEMField[i];
	//		m_derivativesOfEMField[i] = NULL;
	//	}
	//}

	if( m_IDsOfFrequenciesCalculatedByThisPE != NULL ){
		delete [] m_IDsOfFrequenciesCalculatedByThisPE;
		m_IDsOfFrequenciesCalculatedByThisPE = NULL;
	}

	if( m_valuesOfFrequenciesCalculatedByThisPE != NULL ){
		delete [] m_valuesOfFrequenciesCalculatedByThisPE;
		m_valuesOfFrequenciesCalculatedByThisPE = NULL;
	}


}

// Input mesh data from "observe.dat"
void ObservedData::inputObservedData(){

	std::ifstream inFile( "observe.dat", std::ios::in );
	if( inFile.fail() ){
		OutputFiles::m_logFile << "File open error : observe.dat !!" << std::endl;
		exit(1);
	}

	const double convKilometerToMeter = 1000.0;

	while(!inFile.eof())
	{
		std::string strbuf;
		inFile >> strbuf;

		int ibuf(0);
		if( strbuf == "MT" || strbuf == "IMP" ){
			inFile >> m_numStationsMT;
			if( m_numStationsMT > 0 ){
				m_observedStationMT = new ObservedDataStationMT[m_numStationsMT];
			}

			OutputFiles::m_logFile << "# Number of the stations of the impedance tensor : " << m_numStationsMT << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)  Owner Element  Electric field" << std::endl;

			for( int i = 0; i < m_numStationsMT; ++i ){
				
				m_observedStationMT[i].inputObservedData( inFile );

				const int stationID = m_observedStationMT[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::MT ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationMT[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationMT[i].getFrequencyValues(ifrq) );
				}

			}
		}else if( strbuf == "APP_RES_AND_PHS" ){
			inFile >> m_numStationsApparentResistivityAndPhase;
			if( m_numStationsApparentResistivityAndPhase > 0 ){
				m_observedStationApparentResistivityAndPhase = new ObservedDataStationApparentResistivityAndPhase[m_numStationsApparentResistivityAndPhase];
			}

			OutputFiles::m_logFile << "# Number of the stations of the apparent resistivity and phase : " << m_numStationsApparentResistivityAndPhase << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)  Owner Element  Electric field" << std::endl;

			for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
				
				m_observedStationApparentResistivityAndPhase[i].inputObservedData( inFile );

				const int stationID = m_observedStationApparentResistivityAndPhase[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::APP_RES_AND_PHS ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationApparentResistivityAndPhase[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationApparentResistivityAndPhase[i].getFrequencyValues(ifrq) );
				}
			}
		}else if( strbuf == "HTF" ){

			inFile >> m_numStationsHTF;
			if( m_numStationsHTF > 0 ){
				m_observedStationHTF = new ObservedDataStationHTF[m_numStationsHTF];
			}

			OutputFiles::m_logFile << "# Number of the stations of the horizontal magnetic transfer function : " << m_numStationsHTF << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)  Owner Element" << std::endl;

			for( int i = 0; i < m_numStationsHTF; ++i ){

				m_observedStationHTF[i].inputObservedData( inFile );

				const int stationID = m_observedStationHTF[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::HTF ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationHTF[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationHTF[i].getFrequencyValues(ifrq) );
				}
				
			}

		}else if( strbuf == "VTF" ){

			inFile >> m_numStationsVTF;
			if( m_numStationsVTF > 0 ){
				m_observedStationVTF = new ObservedDataStationVTF[m_numStationsVTF];
			}

			OutputFiles::m_logFile << "# Number of the stations of the vertical magnetic transfer function : " << m_numStationsVTF << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)  Owner Element" << std::endl;

			for( int i = 0; i < m_numStationsVTF; ++i ){
				
				m_observedStationVTF[i].inputObservedData( inFile );

				const int stationID = m_observedStationVTF[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::VTF ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationVTF[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationVTF[i].getFrequencyValues(ifrq) );
				}


			}

		}else if( strbuf == "PT" ){

			inFile >> m_numStationsPT;
			if( m_numStationsPT > 0 ){
				m_observedStationPT = new ObservedDataStationPT[m_numStationsPT];
			}

			OutputFiles::m_logFile << "# Number of the stations of the  phase tensor : " << m_numStationsPT << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)  Owner Element  Electric field" << std::endl;

			for( int i = 0; i < m_numStationsPT; ++i ){
				
				m_observedStationPT[i].inputObservedData( inFile );

				const int stationID = m_observedStationPT[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::PT ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationPT[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationPT[i].getFrequencyValues(ifrq) );
				}

			}
		}else if( strbuf == "NMT" ){

			inFile >> m_numStationsNMT;
			if( m_numStationsNMT > 0 ){
				m_observedStationNMT = new ObservedDataStationNMT[m_numStationsNMT];
			}

			OutputFiles::m_logFile << "# Number of the stations of the Network-MT response function : " << m_numStationsNMT << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)" << std::endl;

			for( int i = 0; i < m_numStationsNMT; ++i ){
				
				m_observedStationNMT[i].inputObservedData( inFile );

				const int stationID = m_observedStationNMT[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::NMT ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationNMT[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationNMT[i].getFrequencyValues(ifrq) );
				}

			}

		}else if( strbuf == "NMT2" ){

			inFile >> m_numStationsNMT2;
			if( m_numStationsNMT2 > 0 ){
				m_observedStationNMT2 = new ObservedDataStationNMT2[m_numStationsNMT2];
			}

			OutputFiles::m_logFile << "# Number of the stations of the Network-MT impedance tensor : " << m_numStationsNMT2 << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)" << std::endl;

			for( int i = 0; i < m_numStationsNMT2; ++i ){
				
				m_observedStationNMT2[i].inputObservedData( inFile );

				const int stationID = m_observedStationNMT2[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::NMT2 ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationNMT2[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationNMT2[i].getFrequencyValues(ifrq) );
				}

			}
		}else if( strbuf == "NMT2_APP_RES_AND_PHS" ){
			inFile >> m_numStationsNMT2ApparentResistivityAndPhase;
			if( m_numStationsNMT2ApparentResistivityAndPhase > 0 ){
				m_observedStationNMT2ApparentResistivityAndPhase = new ObservedDataStationNMT2ApparentResistivityAndPhase[m_numStationsNMT2ApparentResistivityAndPhase];
			}

			OutputFiles::m_logFile << "# Number of the stations of the Network-MT apparent resistivity and phase : " << m_numStationsNMT2ApparentResistivityAndPhase << std::endl;
			OutputFiles::m_logFile << "# Station ID     Station ID (Ref)" << std::endl;

			for( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
				
				m_observedStationNMT2ApparentResistivityAndPhase[i].inputObservedData( inFile );

				const int stationID = m_observedStationNMT2ApparentResistivityAndPhase[i].getStationID();
				if( m_stationID2type.find( stationID ) != m_stationID2type.end() ){
					OutputFiles::m_logFile << "Error : Statin ID " <<  stationID << " is already exist !!" << std::endl;
					exit(1);
				}
				m_stationID2type.insert( std::map<int, int>::value_type( stationID, ObservedData::NMT2_APP_RES_AND_PHS ) ); 
				m_stationID2IDAmongEachStationType.insert( std::map<int, int>::value_type( stationID, i ) ); 

				const int nFreq = m_observedStationNMT2ApparentResistivityAndPhase[i].getTotalNumberOfFrequency();
				for( int ifrq = 0; ifrq < nFreq; ++ifrq ){
					checkAndAddNewFrequency( m_observedStationNMT2ApparentResistivityAndPhase[i].getFrequencyValues(ifrq) );
				}
			}
		}else if( strbuf == "ADDITIONAL_OUTPUT_POINTS" ){

			inFile >> m_numAdditinalOutputPoint;

			OutputFiles::m_logFile << "# Number of additional points for output : " << m_numAdditinalOutputPoint << std::endl;

			if( m_numAdditinalOutputPoint > 0 ){
				m_additinalOutputPoint = new AdditinalOutputPoint[m_numAdditinalOutputPoint];
			}

			for( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
				m_additinalOutputPoint[i].inputObservedData( inFile );
			}

		}else if( strbuf == "END" ){
			break;
		}else{
			OutputFiles::m_logFile << "Error : Wrong type of station. " << strbuf << std::endl;
			exit(1);
		}

	}
	
	inFile.close();

	m_numAllStations = m_numStationsMT + m_numStationsApparentResistivityAndPhase + m_numStationsHTF + m_numStationsVTF 
		+ m_numStationsPT + m_numStationsNMT + m_numStationsNMT2 + m_numStationsNMT2ApparentResistivityAndPhase;

	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	OutputFiles::m_logFile << "# Total number of stations : " << m_numAllStations << " ." << std::endl;

	if( m_numStationsMT > 0 ){
		OutputFiles::m_logFile << "#  - Number of impedance tensor stations : " << m_numStationsMT << " ." << std::endl;
	}
	if( m_numStationsApparentResistivityAndPhase > 0 ){
		OutputFiles::m_logFile << "#  - Number of apparent resistivity & phase stations : " << m_numStationsApparentResistivityAndPhase << " ." << std::endl;
	}
	if( m_numStationsHTF > 0 ){
		OutputFiles::m_logFile << "#  - Number of HTF stations : " << m_numStationsHTF << " ." << std::endl;
	}
	if( m_numStationsVTF > 0 ){
		OutputFiles::m_logFile << "#  - Number of VTF stations : " << m_numStationsVTF << " ." << std::endl;
	}
	if( m_numStationsPT > 0 ){
		OutputFiles::m_logFile << "#  - Number of PT stations : " << m_numStationsPT << " ." << std::endl;
	}
	if( m_numStationsNMT > 0 ){
		OutputFiles::m_logFile << "#  - Number of NMT stations : " << m_numStationsNMT << " ." << std::endl;
	}
	if( m_numStationsNMT2 > 0 ){
		OutputFiles::m_logFile << "#  - Number of NMT2 stations (impedance tensor) : " << m_numStationsNMT2 << " ." << std::endl;
	}
	if( m_numStationsNMT2ApparentResistivityAndPhase > 0 ){
		OutputFiles::m_logFile << "#  - Number of NMT2 stations (apparent resistivity & phase) : " << m_numStationsNMT2ApparentResistivityAndPhase << " ." << std::endl;
	}
	if( m_numAdditinalOutputPoint > 0 ){
		OutputFiles::m_logFile << "# Number of additional output points : " << m_numAdditinalOutputPoint << " ." << std::endl;
	}

	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	//m_stationID2TypeAndIDAmongEachStationType = new typeAndIDAmongEachStationType[ m_stationID2type.size() ];

	//int icount(0);
	//for( std::map<int, int>::iterator itr = m_stationID2type.begin(); itr != m_stationID2type.end(); ++itr, ++icount ){
	//	m_stationID2TypeAndIDAmongEachStationType[icount].type = (*itr).second;
	//}
	//icount = 0;
	//for( std::map<int, int>::iterator itr = m_stationID2IDAmongEachStationType.begin(); itr != m_stationID2IDAmongEachStationType.end(); ++itr, ++icount ){
	//	m_stationID2TypeAndIDAmongEachStationType[icount].IDAmongEachStationType = (*itr).second;
	//}

	// Sort all frequencies for which observed data exists
	std::sort( m_frequencyAll.begin(), m_frequencyAll.end());

	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	OutputFiles::m_logFile << "# Forward calculation is performed for the frequencies/periods below." << std::endl;
	OutputFiles::m_logFile << "# " << std::setw(7) << "Freq#" << std::setw(15) << "Frequency[Hz]" << std::setw(15) << "Period[sec]" << std::endl;
	int iFreq(0);
	for( std::vector<double>::iterator it = m_frequencyAll.begin(); it != m_frequencyAll.end(); ++it, ++iFreq ){
		OutputFiles::m_logFile << "# " << std::setw(7) << iFreq << std::setw(15) << *it << std::setw(15) << 1.0 / *it << std::endl; 
	}
	OutputFiles::m_logFile << "#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;

#ifdef _DEBUG_WRITE
	//for( int i = 0; i < m_numAllStations; ++i ){
	//	std::cout << i << " " << m_stationID2TypeAndIDAmongEachStationType[i].type << " " << m_stationID2TypeAndIDAmongEachStationType[i].IDAmongEachStationType  << std::endl;
	//}
	for( std::map<int, int>::iterator itr = m_stationID2type.begin(); itr != m_stationID2type.end(); ++itr ){
		std::cout << (*itr).first << " " << (*itr).second << std::endl;
	}
	for( std::map<int, int>::iterator itr = m_stationID2IDAmongEachStationType.begin(); itr != m_stationID2IDAmongEachStationType.end(); ++itr ){
		std::cout << (*itr).first << " " << (*itr).second << std::endl;
	}
#endif

}

// Input data of distortion matrix from "distortion_iterX.dat"
void ObservedData::inputDistortionMatrixData(){

	const int iterInit = ( AnalysisControl::getInstance() )->getIterationNumInit();

	std::ostringstream inputFile;
	inputFile << "distortion_iter" << iterInit << ".dat";
	std::ifstream inFile( inputFile.str().c_str(), std::ios::in );

	if( inFile.fail() ){
		OutputFiles::m_logFile << "File open error : " << inputFile.str().c_str() << " !!" << std::endl;
		exit(1);
	}

	int numMTSta(0);
	inFile >> numMTSta;
	if( numMTSta != m_numStationsMT + m_numStationsApparentResistivityAndPhase ){
		OutputFiles::m_logFile << "Error : Number of MT station written in " << inputFile.str().c_str() << " ( " << numMTSta << " ) is not equal to the one written in observe.dat ( " << m_numStationsMT + m_numStationsApparentResistivityAndPhase << " ) ." << std::endl;
		exit(1);
	}

	if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){

		for ( int iSta = 0; iSta < numMTSta; ++iSta ){
		
			int StaID(0);
			double distMat[4] = { 0.0, 0.0, 0.0, 0.0 };
			int fixDistortionMatrix;

			inFile >> StaID >> distMat[0] >> distMat[1] >> distMat[2] >> distMat[3] >> fixDistortionMatrix;

			//distMat[ObservedData::COMPONENT_ID_CXX] -= 1.0;
			//distMat[ObservedData::COMPONENT_ID_CYY] -= 1.0;

			bool foundMT(false);
			for ( int indexMT = 0; indexMT < m_numStationsMT; ++indexMT ){
				if( StaID == m_observedStationMT[indexMT].getStationID() ){
					for( int iComp = 0; iComp < 4; ++iComp ){
						m_observedStationMT[indexMT].setDistortionParamsPre( iComp, distMat[iComp] );
						m_observedStationMT[indexMT].setDistortionParams(    iComp, distMat[iComp] );
					}
					m_observedStationMT[indexMT].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundMT = true;
					break;
				}
			}
			bool foundAppResAndPhs(false);
			for ( int indexAppResAndPhs = 0; indexAppResAndPhs < m_numStationsApparentResistivityAndPhase; ++indexAppResAndPhs ){
				if( StaID == m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].getStationID() ){
					for( int iComp = 0; iComp < 4; ++iComp ){
						m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParamsPre( iComp, distMat[iComp] );
						m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParams(    iComp, distMat[iComp] );
					}
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundAppResAndPhs = true;
					break;
				}
			}
			if( !foundMT && !foundAppResAndPhs ){
				OutputFiles::m_logFile << "Error : MT station with ID " << StaID << " is not found in observe.dat." << std::endl;
				exit(1);
			}

		}

	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){

		for ( int iSta = 0; iSta < numMTSta; ++iSta ){
		
			int StaID(0);
			double distParam[4] = { 0.0, 0.0, 0.0, 0.0 };
			int fixDistortionMatrix;

			inFile >> StaID >> distParam[0] >> distParam[1] >> distParam[2] >> distParam[3] >> fixDistortionMatrix;

			if( distParam[ObservedDataStationMT::EX_ROTATION] > 90.0 ){
				OutputFiles::m_logFile << "Warning : Rotation angle of distortion is forced to be less than or equal to 90 (deg.)" << std::endl;
				distParam[ObservedDataStationMT::EX_ROTATION] = 90.0;
			}else if( distParam[ObservedDataStationMT::EX_ROTATION] < -90.0 ){
				OutputFiles::m_logFile << "Warning : Rotation angle of distortion is forced to be larger than or equal to -90 (deg.)" << std::endl;
				distParam[ObservedDataStationMT::EX_ROTATION] = -90.0;
			}
			if( distParam[ObservedDataStationMT::EY_ROTATION] > 90.0 ){
				OutputFiles::m_logFile << "Warning : Rotation angle of distortion is forced to be less than or equal to 90 (deg.)" << std::endl;
				distParam[ObservedDataStationMT::EY_ROTATION] = 90.0;
			}else if( distParam[ObservedDataStationMT::EY_ROTATION] < -90.0 ){
				OutputFiles::m_logFile << "Warning : Rotation angle of distortion is forced to be larger than or equal to -90 (deg.)" << std::endl;
				distParam[ObservedDataStationMT::EY_ROTATION] = -90.0;
			}

			distParam[ObservedDataStationMT::EX_ROTATION] *= CommonParameters::deg2rad;
			distParam[ObservedDataStationMT::EY_ROTATION] *= CommonParameters::deg2rad;

			bool foundMT(false);
			for ( int indexMT = 0; indexMT < m_numStationsMT; ++indexMT ){
				if( StaID == m_observedStationMT[indexMT].getStationID() ){
					for( int iComp = 0; iComp < 4; ++iComp ){
						m_observedStationMT[indexMT].setDistortionParamsPre( iComp, distParam[iComp] );
						m_observedStationMT[indexMT].setDistortionParams(    iComp, distParam[iComp] );
					}
					m_observedStationMT[indexMT].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundMT = true;
					break;
				}
			}
			bool foundAppResAndPhs(false);
			for ( int indexAppResAndPhs = 0; indexAppResAndPhs < m_numStationsApparentResistivityAndPhase; ++indexAppResAndPhs ){
				if( StaID == m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].getStationID() ){
					for( int iComp = 0; iComp < 4; ++iComp ){
						m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParamsPre( iComp, distParam[iComp] );
						m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParams(    iComp, distParam[iComp] );
					}
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundAppResAndPhs = true;
					break;
				}
			}
			if( !foundMT && !foundAppResAndPhs ){
				OutputFiles::m_logFile << "Error : MT station with ID " << StaID << " is not found in observe.dat." << std::endl;
				exit(1);
			}
		}

	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		for ( int iSta = 0; iSta < numMTSta; ++iSta ){		
			int StaID(0);
			int fixDistortionMatrix;

			double Gx(0.0);
			double Gy(0.0);
			inFile >> StaID >> Gx >> Gy >> fixDistortionMatrix;

			bool foundMT(false);
			for ( int indexMT = 0; indexMT < m_numStationsMT; ++indexMT ){
				if( StaID == m_observedStationMT[indexMT].getStationID() ){
					// Gx
					m_observedStationMT[indexMT].setDistortionParamsPre( ObservedDataStationMT::EX_GAIN, Gx );
					m_observedStationMT[indexMT].setDistortionParams(    ObservedDataStationMT::EX_GAIN, Gx );
					// Gy
					m_observedStationMT[indexMT].setDistortionParamsPre( ObservedDataStationMT::EY_GAIN, Gy );
					m_observedStationMT[indexMT].setDistortionParams(    ObservedDataStationMT::EY_GAIN, Gy );
					m_observedStationMT[indexMT].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundMT = true;
					break;
				}
			}
			bool foundAppResAndPhs(false);
			for ( int indexAppResAndPhs = 0; indexAppResAndPhs < m_numStationsApparentResistivityAndPhase; ++indexAppResAndPhs ){
				if( StaID == m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].getStationID() ){
					// Gx
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParamsPre( ObservedDataStationMT::EX_GAIN, Gx );
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParams(    ObservedDataStationMT::EX_GAIN, Gx );
					// Gy
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParamsPre( ObservedDataStationMT::EY_GAIN, Gy );
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setDistortionParams(    ObservedDataStationMT::EY_GAIN, Gy );
					m_observedStationApparentResistivityAndPhase[indexAppResAndPhs].setFixDistortionMatrix( fixDistortionMatrix == 1 ? true : false );
					foundAppResAndPhs = true;
					break;
				}
			}
			if( !foundMT && !foundAppResAndPhs ){
				OutputFiles::m_logFile << "Error : MT station with ID " << StaID << " is not found in observe.dat." << std::endl;
				exit(1);
			}
		}
	}
	else{
		OutputFiles::m_logFile << "Type of distortion is wrong !!" << std::endl;
		exit(1);
	}

	// Calculate total number of components of distortion matrix difference whose value is not fixed
	calcNumDistortionParamsNotFixed();

	inFile.close();
	
}

// Get total number of different frquencies
int ObservedData::getTotalNumberOfDifferenetFrequencies() const{

	return static_cast<int>( m_frequencyAll.size() );

}

// Get frequency value
double ObservedData::getFrequencyValue( const int num ) const{

	//if( num < 0 || num >= getTotalNumberOfDifferenetFrequencies() ){
	//	OutputFiles::m_logFile << "Error : num is out of range at getFrequencyValue !! num = " << num << " ." << std::endl;
	//	exit(1);
	//}
	assert( num >= 0 ); 
	assert( num < getTotalNumberOfDifferenetFrequencies() );

	return m_frequencyAll[num];

}

// Find element including each station
void ObservedData::findElementIncludingEachStation(){

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].findElementIncludingStation();
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].findElementIncludingStation();
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].findElementIncludingStation();
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].findElementIncludingStation();
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].findElementIncludingStation();
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].findElementsIncludingDipole();
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].findElementsIncludingDipoles();
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].findElementsIncludingDipoles();
	}

	for ( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
		m_additinalOutputPoint[i].findElementIncludingStation();
	}

}

// Add new frequency to the array of all frequencies for which observed data exists after checking 
// wheter inputed frequency have not been contained
void ObservedData::checkAndAddNewFrequency( const double freq ){

	const double EPS = 1.0E-8;
	const double thresholdFrequency = 1e-8;

	if( freq < thresholdFrequency ){
		OutputFiles::m_logFile << "Error : Inputed frequency " << freq << " is less than threshold value " << thresholdFrequency << " ." << std::endl;
		exit(1);
	}

	for( std::vector<double>::iterator it = m_frequencyAll.begin(); it != m_frequencyAll.end(); ++it ){
		if( std::fabs( freq - *it ) < EPS ){
			return;
		}
	}

	m_frequencyAll.push_back( freq );

}
	
// Check whether magnetic field value of the specified station is used to calculate response functions of the specified frequency
bool ObservedData::isUseMagneticFieldOfTheStation( const double freq, const int staID ){

	for ( int i = 0; i < m_numStationsMT; ++i ){
		if( staID == m_observedStationMT[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		if( staID == m_observedStationApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		if( staID == m_observedStationHTF[i].getIDOfMagneticFieldStation() ||
			staID == m_observedStationHTF[i].getStationID() ){
			if( m_observedStationHTF[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		if( staID == m_observedStationVTF[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationVTF[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		if( staID == m_observedStationPT[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationPT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		if( staID == m_observedStationNMT[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationNMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		if( staID == m_observedStationNMT2[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationNMT2[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		if( staID == m_observedStationNMT2ApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ){
			if( m_observedStationNMT2ApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
				return true; 
			}
		}
	}

	return false;
}

// Get instance of the station where magnetic field is observed
ObservedDataStationPoint* ObservedData::getInstanceOfMagneticStation( const int IDOfMagSta ) const{

	std::map<int, int>::const_iterator itr = m_stationID2IDAmongEachStationType.find( IDOfMagSta );
	if( itr == m_stationID2IDAmongEachStationType.end() ){
		OutputFiles::m_logFile << "Error : Reference station ( " << IDOfMagSta << " ) cannnot be found !!" << std::endl;
		exit(1);
	}

	std::map<int, int>::const_iterator itr2 = m_stationID2type.find( IDOfMagSta );
	if( itr2 == m_stationID2type.end() ){
		OutputFiles::m_logFile << "Error : Reference station ( " << IDOfMagSta << " ) cannnot be found !!" << std::endl;
		exit(1);
	}

	const int num = (*itr).second;

	switch( (*itr2).second ){
		case MT:
			return &m_observedStationMT[num];
			break;
		case APP_RES_AND_PHS:
			return &m_observedStationApparentResistivityAndPhase[num];
			break;
		case HTF:
			return &m_observedStationHTF[num];
			break;
		case VTF:
			return &m_observedStationVTF[num];
			break;
		case PT:
			return &m_observedStationPT[num];
			break;
		case NMT:
			// go through
		case NMT2:
			// go through
		case NMT2_APP_RES_AND_PHS:
			OutputFiles::m_logFile << "Error : Network-MT station cannnot have magnetic data !! IDOfMagSta = " << IDOfMagSta << std::endl;
			exit(1);
			break;
		default:
			OutputFiles::m_logFile << "Error : Wrong type ID of station !!" << (*itr2).second << " ." << std::endl;
			exit(1);
			break;
	}

}

// Calculate total number of components of distortion parameters whose value is not fixed
void ObservedData::calcNumDistortionParamsNotFixed(){

	int icount(0);
	if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			if( m_observedStationMT[i].getTotalNumberOfFrequency() > 0 && !m_observedStationMT[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				for( int j = 0; j < 4; ++j ){
					m_observedStationMT[i].setIDOfDistortionParams( j, icount );
					if( ObservedDataStationMT::COMPONENT_ID_CXX || ObservedDataStationMT::COMPONENT_ID_CXY ||
						ObservedDataStationMT::COMPONENT_ID_CYX || ObservedDataStationMT::COMPONENT_ID_CYY ){
						m_typesOfDistortionParamsNotFixed.push_back(j);
					}else{
						OutputFiles::m_logFile << "Error : Component number of distortion matrix is wrong !!" << j << std::endl;
						exit(1);
					}
					++icount;
				}
			} 
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			if( m_observedStationApparentResistivityAndPhase[i].getTotalNumberOfFrequency() > 0 && !m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				for( int j = 0; j < 4; ++j ){
					m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( j, icount );
					if( ObservedDataStationMT::COMPONENT_ID_CXX || ObservedDataStationMT::COMPONENT_ID_CXY ||
						ObservedDataStationMT::COMPONENT_ID_CYX || ObservedDataStationMT::COMPONENT_ID_CYY ){
						m_typesOfDistortionParamsNotFixed.push_back(j);
					}else{
						OutputFiles::m_logFile << "Error : Component number of distortion matrix is wrong !!" << j << std::endl;
						exit(1);
					}
					++icount;
				}
			} 
		}
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			if( m_observedStationMT[i].getTotalNumberOfFrequency() > 0 && !m_observedStationMT[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				for( int j = 0; j < 4; ++j ){
					m_observedStationMT[i].setIDOfDistortionParams( j, icount );
					if( ObservedDataStationMT::EX_GAIN || ObservedDataStationMT::EY_GAIN ||
						ObservedDataStationMT::EX_ROTATION || ObservedDataStationMT::EY_ROTATION ){
						m_typesOfDistortionParamsNotFixed.push_back(j);
					}else{
						OutputFiles::m_logFile << "Error : Component number of gains or rotations is wrong !!" << j << std::endl;
						exit(1);
					}
					++icount;
				}
			} 
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			const double eps = 1.0 * CommonParameters::deg2rad;
			if( m_observedStationApparentResistivityAndPhase[i].getTotalNumberOfFrequency() > 0 && !m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				// Gx
				if( m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::XX ) &&
					m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::XY ) &&
					fabs(m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION)) < eps &&
					fabs(m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION)) < eps ){
					OutputFiles::m_logFile << "Notice : Gx of MT station (ID = " << m_observedStationApparentResistivityAndPhase[i].getStationID()
						<< ") is not estimated because XX and XY components of apparent resistivity are not used as data and rotations are nearly zeros." << std::endl;
				}else{
					m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EX_GAIN, icount );
					m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EX_GAIN);
					++icount;
				}
				// Gy
				if( m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::YX ) &&
					m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::YY ) &&
					fabs(m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION)) < eps &&
					fabs(m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION)) < eps ){
					OutputFiles::m_logFile << "Notice : Gy of MT station (ID = " << m_observedStationApparentResistivityAndPhase[i].getStationID()
						<< ") is not estimated because YX and YY components of apparent resistivity are not used as data and rotations are nearly zeros." << std::endl;
				}else{
					m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EY_GAIN, icount );
					m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EY_GAIN);
					++icount;
				}
				// Rotation X
				m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EX_ROTATION, icount );
				m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EX_ROTATION);
				++icount;
				// Rotation Y
				m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EY_ROTATION, icount );
				m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EY_ROTATION);
				++icount;
			} 
		}
	}
	else if( ( AnalysisControl::getInstance() )->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			if( m_observedStationMT[i].getTotalNumberOfFrequency() > 0 && !m_observedStationMT[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				// Gx
				m_observedStationMT[i].setIDOfDistortionParams( ObservedDataStationMT::EX_GAIN, icount );
				m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EX_GAIN);
				++icount;
				// Gy
				m_observedStationMT[i].setIDOfDistortionParams( ObservedDataStationMT::EY_GAIN, icount );
				m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EY_GAIN);
				++icount;
			} 
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			if( m_observedStationApparentResistivityAndPhase[i].getTotalNumberOfFrequency() > 0 && !m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ){
				// Exclude station without data and fixed distortion matrix
				// Gx
				if( m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::XX ) &&
					m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::XY ) ){
					OutputFiles::m_logFile << "Notice : Gx of MT station (ID = " << m_observedStationApparentResistivityAndPhase[i].getStationID()
						<< ") is not estimated because XX and XY components of apparent resistivity are not used as data." << std::endl;
				}else{
					m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EX_GAIN, icount );
					m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EX_GAIN);
					++icount;
				}
				// Gy
				if( m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::YX ) &&
					m_observedStationApparentResistivityAndPhase[i].isApparentResistivityIgnoredForAllFrequencies( ObservedDataStationMT::YY ) ){
					OutputFiles::m_logFile << "Notice : Gy of MT station (ID = " << m_observedStationApparentResistivityAndPhase[i].getStationID()
						<< ") is not estimated because YX and YY components of apparent resistivity are not used as data." << std::endl;
				}else{
					m_observedStationApparentResistivityAndPhase[i].setIDOfDistortionParams( ObservedDataStationMT::EY_GAIN, icount );
					m_typesOfDistortionParamsNotFixed.push_back(ObservedDataStationMT::EY_GAIN);
					++icount;
				}
			}
		}
	}
	else{
		OutputFiles::m_logFile << "Error : Wrong type of distortion." << std::endl;
		exit(1);
	}

}

//// Get station type from statin ID 
//int ObservedData::getTypeFromStationID( const int staID ) const{
//	return m_stationID2TypeAndIDAmongEachStationType[staID].type;
//}
//
//// Get IDs among each station type of all stations from statin ID 
//int ObservedData::getIDAmongEachStationTypeFromStationID( const int staID ) const{
//	return m_stationID2TypeAndIDAmongEachStationType[staID].IDAmongEachStationType;
//}

// Calulate EM field of all stations
void ObservedData::calculateEMFieldOfAllStations( const Forward3D* const ptrForward3D, const double freq, const int iPol, const int ifreq ){

	// Initialization
	m_numInterpolatorVectors = 0;

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].initializeElectricField( iPol );
		m_observedStationMT[i].initializeHorizontalMagneticField( iPol );

		if( m_observedStationMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfEx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfEy = incrementNumInterpolatorVectors();
			m_observedStationMT[i].calculateElectricField( ptrForward3D, rhsVectorIDOfEx, rhsVectorIDOfEy );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationMT[i].getStationID() ) ){
			const int rhsVectorIDOfHx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfHy = incrementNumInterpolatorVectors();
			m_observedStationMT[i].calculateHorizontalMagneticField( ptrForward3D, rhsVectorIDOfHx, rhsVectorIDOfHy );
		}
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].initializeElectricField( iPol );
		m_observedStationApparentResistivityAndPhase[i].initializeHorizontalMagneticField( iPol );

		if( m_observedStationApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfEx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfEy = incrementNumInterpolatorVectors();
			m_observedStationApparentResistivityAndPhase[i].calculateElectricField( ptrForward3D, rhsVectorIDOfEx, rhsVectorIDOfEy );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationApparentResistivityAndPhase[i].getStationID() ) ){
			const int rhsVectorIDOfHx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfHy = incrementNumInterpolatorVectors();
			m_observedStationApparentResistivityAndPhase[i].calculateHorizontalMagneticField( ptrForward3D, rhsVectorIDOfHx, rhsVectorIDOfHy );
		}
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].initializeHorizontalMagneticField( iPol );

		if( m_observedStationHTF[i].getFreqIDsAmongThisPE( freq ) >= 0 || 
			isUseMagneticFieldOfTheStation( freq, m_observedStationHTF[i].getStationID() ) ){
			const int rhsVectorIDOfHx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfHy = incrementNumInterpolatorVectors();
			m_observedStationHTF[i].calculateHorizontalMagneticField( ptrForward3D, rhsVectorIDOfHx, rhsVectorIDOfHy );
		}
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].initializeVerticalMagneticField( iPol );
		m_observedStationVTF[i].initializeHorizontalMagneticField( iPol );

		if( m_observedStationVTF[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfHz = incrementNumInterpolatorVectors();
			m_observedStationVTF[i].calculateVerticalMagneticField( ptrForward3D, rhsVectorIDOfHz );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationVTF[i].getStationID() ) ){
			const int rhsVectorIDOfHx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfHy = incrementNumInterpolatorVectors();
			m_observedStationVTF[i].calculateHorizontalMagneticField( ptrForward3D, rhsVectorIDOfHx, rhsVectorIDOfHy );
		}
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].initializeElectricField( iPol );
		m_observedStationPT[i].initializeHorizontalMagneticField( iPol );

		if( m_observedStationPT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfEx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfEy = incrementNumInterpolatorVectors();
			m_observedStationPT[i].calculateElectricField( ptrForward3D, rhsVectorIDOfEx, rhsVectorIDOfEy );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationPT[i].getStationID() ) ){
			const int rhsVectorIDOfHx = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfHy = incrementNumInterpolatorVectors();
			m_observedStationPT[i].calculateHorizontalMagneticField( ptrForward3D, rhsVectorIDOfHx, rhsVectorIDOfHy );
		}
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].initializeVoltageDifferences( iPol );

		if( m_observedStationNMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfVoltageDifference = incrementNumInterpolatorVectors();
			m_observedStationNMT[i].calculateVoltageDifferences( ptrForward3D, rhsVectorIDOfVoltageDifference );
		}
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].initializeVoltageDifferences( iPol );

		if( m_observedStationNMT2[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfVoltageDifference1st = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfVoltageDifference2nd = incrementNumInterpolatorVectors();
			m_observedStationNMT2[i].calculateVoltageDifferences( ptrForward3D, rhsVectorIDOfVoltageDifference1st, rhsVectorIDOfVoltageDifference2nd );
		}
	}
	
	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].initializeVoltageDifferences( iPol );

		if( m_observedStationNMT2ApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			const int rhsVectorIDOfVoltageDifference1st = incrementNumInterpolatorVectors();
			const int rhsVectorIDOfVoltageDifference2nd = incrementNumInterpolatorVectors();
			m_observedStationNMT2ApparentResistivityAndPhase[i].calculateVoltageDifferences( ptrForward3D, rhsVectorIDOfVoltageDifference1st, rhsVectorIDOfVoltageDifference2nd );
		}
	}

	for ( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
		m_additinalOutputPoint[i].calculateEMField( ifreq, ptrForward3D );
	}

}

// Calulate response function of all stations
void ObservedData::calculateResponseFunctionOfAllStations( const int freqIDAmongThisPE ){
	

	const double freq = getValuesOfFrequenciesCalculatedByThisPE( freqIDAmongThisPE );

	if( m_numObservedDataThisPEAccumulated == NULL ){
		OutputFiles::m_logFile << "Error : m_numObservedDataThisPEAccumulated is NULL" << std::endl;
		exit(1);
	}


	// Initialization
	int icount = 0;

#ifdef _DEBUG_WRITE
	std::cout << "freq = " << freq << std::endl;
#endif

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].calculateImpedanceTensor( freq, getInstanceOfMagneticStation( m_observedStationMT[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].calculateApparentResistivityAndPhase( freq, getInstanceOfMagneticStation( m_observedStationApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].calculateHTF( freq, getInstanceOfMagneticStation( m_observedStationHTF[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].calculateVTF( freq, getInstanceOfMagneticStation( m_observedStationVTF[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].calculatePhaseTensor( freq, getInstanceOfMagneticStation( m_observedStationPT[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].calculateNetworkMTResponse( freq, getInstanceOfMagneticStation( m_observedStationNMT[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].calculateImpedanceTensor( freq, getInstanceOfMagneticStation( m_observedStationNMT2[i].getIDOfMagneticFieldStation() ), icount );
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].calculateApparentResistivityAndPhase( freq, getInstanceOfMagneticStation( m_observedStationNMT2ApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ), icount );
	}

	m_numObservedDataThisPEAccumulated[freqIDAmongThisPE+1] = m_numObservedDataThisPEAccumulated[freqIDAmongThisPE] + icount;

#ifdef _DEBUG_WRITE
	std::cout << "freqIDAmongThisPE m_numObservedDataThisPEAccumulated[freqIDAmongThisPE+1] : " << freqIDAmongThisPE << " " << m_numObservedDataThisPEAccumulated[freqIDAmongThisPE+1] << std::endl;
#endif
	
}

// Initialize response functions and errors of all stations
void ObservedData::initializeResponseFunctionsAndErrorsOfAllStations(){

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].initializeImpedanceTensorsAndErrors();
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].initializeApparentResistivityPhaseAndErrors();
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].initializeHTFsAndErrors();
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].initializeVTFsAndErrors();
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].initializePhaseTensorsAndErrors();
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].initializeNetworkMTResponsesAndErrors();
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].initializeImpedanceTensorsAndErrors();
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].initializeApparentResistivityPhaseAndErrors();
	}

	for ( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
		m_additinalOutputPoint[i].initializeEMfield( getNumOfFrequenciesCalculatedByThisPE() ); 
	}

}

// Allocate memory for the calculated values and errors of all stations 
void ObservedData::allocateMemoryForCalculatedValuesOfAllStations(){

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	const int nFreqCalculatedByThisPE = getNumOfFrequenciesCalculatedByThisPE();

	double* freqCalculatedByThisPE = new double[nFreqCalculatedByThisPE];
	for( int i = 0; i < nFreqCalculatedByThisPE; ++i ){
		freqCalculatedByThisPE[i] = getValuesOfFrequenciesCalculatedByThisPE(i);
	}

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationMT[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationApparentResistivityAndPhase[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationHTF[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationVTF[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationPT[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationNMT[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationNMT2[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].setupFrequenciesCalculatedByThisPE( nFreqCalculatedByThisPE, freqCalculatedByThisPE );
		m_observedStationNMT2ApparentResistivityAndPhase[i].allocateMemoryForCalculatedValues();
	}

	for ( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
		m_additinalOutputPoint[i].allocateMemoryForCalculatedValues( getNumOfFrequenciesCalculatedByThisPE() );
	}

}

// Output calculated values of all stations
void ObservedData::outputCalculatedValuesOfAllStations() const{

	const bool useImpedanceTensorInsteadOfPhase = ( AnalysisControl::getInstance()->getApparentResistivityAndPhaseTreatmentOption() == AnalysisControl::USE_Z_IF_SIGN_OF_RE_Z_DIFFER );

	if( m_numStationsMT > 0 ){
		fprintf( OutputFiles::m_csvFile, "MT,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Cal,     Im(Zxx)_Cal,     Re(Zxy)_Cal,     Im(Zxy)_Cal,     Re(Zyx)_Cal,     Im(Zyx)_Cal,     Re(Zyy)_Cal,     Im(Zyy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Res,     Im(Zxx)_Res,     Re(Zxy)_Res,     Im(Zxy)_Res,     Re(Zyx)_Res,     Im(Zyx)_Res,     Re(Zyy)_Res,     Im(Zyy)_Res,");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Obs,     Im(Zxx)_Obs,     Re(Zxy)_Obs,     Im(Zxy)_Obs,     Re(Zyx)_Obs,     Im(Zyx)_Obs,     Re(Zyy)_Obs,     Im(Zyy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "      Re(Zxx)_SD,      Im(Zxx)_SD,      Re(Zxy)_SD,      Im(Zxy)_SD,      Re(Zyx)_SD,      Im(Zyx)_SD,      Re(Zyy)_SD,      Im(Zyy)_SD,\n");

		for ( int i = 0; i < m_numStationsMT; ++i ){
			m_observedStationMT[i].outputCalculatedValues();
		}

	}

	if( m_numStationsApparentResistivityAndPhase > 0 ){
		fprintf( OutputFiles::m_csvFile, "APP_RES_AND_PHS,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Cal,    Phase_xx_Cal,  App.Res_xy_Cal,    Phase_xy_Cal,  App.Res_yx_Cal,    Phase_yx_Cal,  App.Res_yy_Cal,    Phase_yy_Cal,");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Res,    Phase_xx_Res,  App.Res_xy_Res,    Phase_xy_Res,  App.Res_yx_Res,    Phase_yx_Res,  App.Res_yy_Res,    Phase_yy_Res,");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Obs,    Phase_xx_Obs,  App.Res_xy_Obs,    Phase_xy_Obs,  App.Res_yx_Obs,    Phase_yx_Obs,  App.Res_yy_Obs,    Phase_yy_Obs,");
		fprintf( OutputFiles::m_csvFile, "   App.Res_xx_SD,     Phase_xx_SD,   App.Res_xy_SD,     Phase_xy_SD,   App.Res_yx_SD,     Phase_yx_SD,   App.Res_yy_SD,     Phase_yy_SD,");
		if(useImpedanceTensorInsteadOfPhase){
			fprintf( OutputFiles::m_csvFile, "Zxx,Zxy,Zyx,Zyy,\n");
		}else{
			fprintf( OutputFiles::m_csvFile, "\n");
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			m_observedStationApparentResistivityAndPhase[i].outputCalculatedValues();
		}

	}

	if( m_numStationsHTF > 0 ){
		fprintf( OutputFiles::m_csvFile, "HTF,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "     Re(Txx)_Cal,     Im(Txx)_Cal,     Re(Txy)_Cal,     Im(Txy)_Cal,     Re(Tyx)_Cal,     Im(Tyx)_Cal,     Re(Tyy)_Cal,     Im(Tyy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "     Re(Txx)_Res,     Im(Txx)_Res,     Re(Txy)_Res,     Im(Txy)_Res,     Re(Tyx)_Res,     Im(Tyx)_Res,     Re(Tyy)_Res,     Im(Tyy)_Res,");
		fprintf( OutputFiles::m_csvFile, "     Re(Txx)_Obs,     Im(Txx)_Obs,     Re(Txy)_Obs,     Im(Txy)_Obs,     Re(Tyx)_Obs,     Im(Tyx)_Obs,     Re(Tyy)_Obs,     Im(Tyy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "      Re(Txx)_SD,      Im(Txx)_SD,      Re(Txy)_SD,      Im(Txy)_SD,      Re(Tyx)_SD,      Im(Tyx)_SD,      Re(Tyy)_SD,      Im(Tyy)_SD,\n");

		for ( int i = 0; i < m_numStationsHTF; ++i ){
			m_observedStationHTF[i].outputCalculatedValues();
		}

	}

	if( m_numStationsVTF > 0 ){
		fprintf( OutputFiles::m_csvFile, "VTF,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "     Re(Tzx)_Cal,     Im(Tzx)_Cal,     Re(Tzy)_Cal,     Im(Tzy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "     Re(Tzx)_Res,     Im(Tzx)_Res,     Re(Tzy)_Res,     Im(Tzy)_Res,");
		fprintf( OutputFiles::m_csvFile, "     Re(Tzx)_Obs,     Im(Tzx)_Obs,     Re(Tzy)_Obs,     Im(Tzy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "      Re(Tzx)_SD,      Im(Tzx)_SD,      Re(Tzy)_SD,      Im(Tzy)_SD,\n");

		for ( int i = 0; i < m_numStationsVTF; ++i ){
			m_observedStationVTF[i].outputCalculatedValues();
		}

	}

	if( m_numStationsPT > 0 ){
		fprintf( OutputFiles::m_csvFile, "PT,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "    Re(PTxx)_Cal,    Re(PTxy)_Cal,    Re(PTyx)_Cal,    Re(PTyy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "    Re(PTxx)_Res,    Re(PTxy)_Res,    Re(PTyx)_Res,    Re(PTyy)_Res,");
		fprintf( OutputFiles::m_csvFile, "    Re(PTxx)_Obs,    Re(PTxy)_Obs,    Re(PTyx)_Obs,    Re(PTyy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "     Re(PTxx)_SD,     Re(PTxy)_SD,     Re(PTyx)_SD,     Re(PTyy)_SD,\n");

		for ( int i = 0; i < m_numStationsPT; ++i ){
			m_observedStationPT[i].outputCalculatedValues();
		}
	}

	if( m_numStationsNMT > 0 ){
		fprintf( OutputFiles::m_csvFile, "NMT,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "      Re(Yx)_Cal,      Im(Yx)_Cal,      Re(Yy)_Cal,      Im(Yy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "      Re(Yx)_Res,      Im(Yx)_Res,      Re(Yy)_Res,      Im(Yy)_Res,");
		fprintf( OutputFiles::m_csvFile, "      Re(Yx)_Obs,      Im(Yx)_Obs,      Re(Yy)_Obs,      Im(Yy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "       Re(Yx)_SD,       Im(Yx)_SD,       Re(Yy)_SD,       Im(Yy)_SD,\n");

		for ( int i = 0; i < m_numStationsNMT; ++i ){
			m_observedStationNMT[i].outputCalculatedValues();
		}

	}

	if( m_numStationsNMT2 > 0 ){
		fprintf( OutputFiles::m_csvFile, "NMT2,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Cal,     Im(Zxx)_Cal,     Re(Zxy)_Cal,     Im(Zxy)_Cal,     Re(Zyx)_Cal,     Im(Zyx)_Cal,     Re(Zyy)_Cal,     Im(Zyy)_Cal,");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Res,     Im(Zxx)_Res,     Re(Zxy)_Res,     Im(Zxy)_Res,     Re(Zyx)_Res,     Im(Zyx)_Res,     Re(Zyy)_Res,     Im(Zyy)_Res,");
		fprintf( OutputFiles::m_csvFile, "     Re(Zxx)_Obs,     Im(Zxx)_Obs,     Re(Zxy)_Obs,     Im(Zxy)_Obs,     Re(Zyx)_Obs,     Im(Zyx)_Obs,     Re(Zyy)_Obs,     Im(Zyy)_Obs,");
		fprintf( OutputFiles::m_csvFile, "      Re(Zxx)_SD,      Im(Zxx)_SD,      Re(Zxy)_SD,      Im(Zxy)_SD,      Re(Zyx)_SD,      Im(Zyx)_SD,      Re(Zyy)_SD,      Im(Zyy)_SD,\n");

		for ( int i = 0; i < m_numStationsNMT2; ++i ){
			m_observedStationNMT2[i].outputCalculatedValues();
		}

	}

	if( m_numStationsNMT2ApparentResistivityAndPhase > 0 ){
		fprintf( OutputFiles::m_csvFile, "NMT2_APP_RES_AND_PHS,\n");
		fprintf( OutputFiles::m_csvFile, "     StaID,        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Cal,    Phase_xx_Cal,  App.Res_xy_Cal,    Phase_xy_Cal,  App.Res_yx_Cal,    Phase_yx_Cal,  App.Res_yy_Cal,    Phase_yy_Cal,");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Res,    Phase_xx_Res,  App.Res_xy_Res,    Phase_xy_Res,  App.Res_yx_Res,    Phase_yx_Res,  App.Res_yy_Res,    Phase_yy_Res,");
		fprintf( OutputFiles::m_csvFile, "  App.Res_xx_Obs,    Phase_xx_Obs,  App.Res_xy_Obs,    Phase_xy_Obs,  App.Res_yx_Obs,    Phase_yx_Obs,  App.Res_yy_Obs,    Phase_yy_Obs,");
		fprintf( OutputFiles::m_csvFile, "   App.Res_xx_SD,     Phase_xx_SD,   App.Res_xy_SD,     Phase_xy_SD,   App.Res_yx_SD,     Phase_yx_SD,   App.Res_yy_SD,     Phase_yy_SD,");
		if(useImpedanceTensorInsteadOfPhase){
			fprintf( OutputFiles::m_csvFile, "Zxx,Zxy,Zyx,Zyy,\n");
		}else{
			fprintf( OutputFiles::m_csvFile, "\n");
		}

		for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
			m_observedStationNMT2ApparentResistivityAndPhase[i].outputCalculatedValues();
		}
	}

	if( m_numAdditinalOutputPoint > 0 ){
		fprintf( OutputFiles::m_csvFile, "ADDITINAL_OUTPUT_POINT,\n");
		fprintf( OutputFiles::m_csvFile, "           X[km],           Y[km],           Z[km],        Freq[Hz],");
		fprintf( OutputFiles::m_csvFile, "    Re(Ex)_ExPol,    Im(Ex)_ExPol,    Re(Ey)_ExPol,    Im(Ey)_ExPol,    Re(Ez)_ExPol,    Im(Ez)_ExPol,");
		fprintf( OutputFiles::m_csvFile, "    Re(Hx)_ExPol,    Im(Hx)_ExPol,    Re(Hy)_ExPol,    Im(Hy)_ExPol,    Re(Hz)_ExPol,    Im(Hz)_ExPol,");
		fprintf( OutputFiles::m_csvFile, "    Re(Ex)_EyPol,    Im(Ex)_EyPol,    Re(Ey)_EyPol,    Im(Ey)_EyPol,    Re(Ez)_EyPol,    Im(Ez)_EyPol,");
		fprintf( OutputFiles::m_csvFile, "    Re(Hx)_EyPol,    Im(Hx)_EyPol,    Re(Hy)_EyPol,    Im(Hy)_EyPol,    Re(Hz)_EyPol,    Im(Hz)_EyPol,\n");

		for ( int i = 0; i < m_numAdditinalOutputPoint; ++i ){
			m_additinalOutputPoint[i].outputCalculatedValues( getNumOfFrequenciesCalculatedByThisPE(), m_valuesOfFrequenciesCalculatedByThisPE );
		}

	}

	fflush(OutputFiles::m_csvFile);

}

// Calulate right-hand sides matrix consisting of interpolator vectors
void ObservedData::calculateInterpolatorVectors( Forward3D* const ptrForward3D, const double freq ){

	// Initialize right-hand sides matrix for calculating derivatives of EM field
	ptrForward3D->initializeRhsVectors( m_numInterpolatorVectors );

	for ( int i = 0; i < m_numStationsMT; ++i ){
		if( m_observedStationMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationMT[i].calcInterpolatorVectorOfElectricField( ptrForward3D );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationMT[i].getStationID() ) ){
			m_observedStationMT[i].calcInterpolatorVectorOfHorizontalMagneticField( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		if( m_observedStationApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationApparentResistivityAndPhase[i].calcInterpolatorVectorOfElectricField( ptrForward3D );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationApparentResistivityAndPhase[i].getStationID() ) ){
			m_observedStationApparentResistivityAndPhase[i].calcInterpolatorVectorOfHorizontalMagneticField( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		if( m_observedStationHTF[i].getFreqIDsAmongThisPE( freq ) >= 0 || 
			isUseMagneticFieldOfTheStation( freq, m_observedStationHTF[i].getStationID() ) ){
			m_observedStationHTF[i].calcInterpolatorVectorOfHorizontalMagneticField( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		if( m_observedStationVTF[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationVTF[i].calcInterpolatorVectorOfVerticalMagneticField( ptrForward3D );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationVTF[i].getStationID() ) ){
			m_observedStationVTF[i].calcInterpolatorVectorOfHorizontalMagneticField( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		if( m_observedStationPT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationPT[i].calcInterpolatorVectorOfElectricField( ptrForward3D );
		}
		if( isUseMagneticFieldOfTheStation( freq, m_observedStationPT[i].getStationID() ) ){
			m_observedStationPT[i].calcInterpolatorVectorOfHorizontalMagneticField( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		if( m_observedStationNMT[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationNMT[i].calcInterpolatorVectorOfVoltageDifference( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		if( m_observedStationNMT2[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationNMT2[i].calcInterpolatorVectorOfVoltageDifference( ptrForward3D );
		}
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		if( m_observedStationNMT2ApparentResistivityAndPhase[i].getFreqIDsAmongThisPE( freq ) >= 0 ){
			m_observedStationNMT2ApparentResistivityAndPhase[i].calcInterpolatorVectorOfVoltageDifference( ptrForward3D );
		}
	}

}

// Calulate sensitivity matrix
void ObservedData::calculateSensitivityMatrix( const std::complex<double>* const derivativesOfEMFieldExPol, const std::complex<double>* const derivativesOfEMFieldEyPol,
	const double freq, const int nModel, double* const sensitivityMatrix ) const{

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationMT[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		m_observedStationHTF[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationHTF[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		m_observedStationVTF[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationVTF[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		m_observedStationPT[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationPT[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		m_observedStationNMT[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationNMT[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		m_observedStationNMT2[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationNMT2[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		m_observedStationNMT2ApparentResistivityAndPhase[i].calculateSensitivityMatrix( freq, nModel,
			getInstanceOfMagneticStation( m_observedStationNMT2ApparentResistivityAndPhase[i].getIDOfMagneticFieldStation() ),
			derivativesOfEMFieldExPol,
			derivativesOfEMFieldEyPol,
			sensitivityMatrix );
	}

}

// Calculate data vector of this PE
void ObservedData::calculateResidualVectorOfDataThisPE( double* const vector ) const{

	const int nfreq = getNumOfFrequenciesCalculatedByThisPE();

	for( int ifreq = 0; ifreq < nfreq; ++ifreq ){

		const double freq = getValuesOfFrequenciesCalculatedByThisPE( ifreq );

		if( m_numObservedDataThisPEAccumulated == NULL ){
			OutputFiles::m_logFile << "Error : m_numObservedDataThisPEAccumulated is NULL" << std::endl;
			exit(1);
		}

		const int offset = m_numObservedDataThisPEAccumulated[ifreq];

		for ( int i = 0; i < m_numStationsMT; ++i ){
			m_observedStationMT[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			m_observedStationApparentResistivityAndPhase[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsHTF; ++i ){
			m_observedStationHTF[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsVTF; ++i ){
			m_observedStationVTF[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsPT; ++i ){
			m_observedStationPT[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsNMT; ++i ){
			m_observedStationNMT[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsNMT2; ++i ){
			m_observedStationNMT2[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

		for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
			m_observedStationNMT2ApparentResistivityAndPhase[i].calculateResidualVectorOfDataThisPE( freq, offset, vector );
		}

	}

}

// Get number of interpolator vectors
int ObservedData::getNumInterpolatorVectors() const{
	return m_numInterpolatorVectors;
}

// Get total number of observed data
int ObservedData::getNumObservedDataThisPETotal() const{
	return m_numObservedDataThisPEAccumulated[ m_numOfFrequenciesCalculatedByThisPE ];
}

// Get number of observed data
int ObservedData::getNumObservedDataThisPE( const int ifreq ) const{

	//if( ifreq < 0 || ifreq > m_numOfFrequenciesCalculatedByThisPE ){
	//	OutputFiles::m_logFile << "Error : ifreq is out of range. ifreq = " << ifreq << std::endl;
	//	exit(1);
	//}
	assert( ifreq >= 0 ); 
	assert( ifreq <= m_numOfFrequenciesCalculatedByThisPE );

	//return m_numObservedDataThisPE[ifreq];
	return m_numObservedDataThisPEAccumulated[ifreq+1] - m_numObservedDataThisPEAccumulated[ifreq];

}

// Get number of observed data accumulated
int ObservedData::getNumObservedDataThisPEAccumulated( const int num ) const{

	//if( num < 0 || num > m_numOfFrequenciesCalculatedByThisPE + 1 ){
	//	OutputFiles::m_logFile << "Error : num is out of range. num = " << num << std::endl;
	//	exit(1);
	//}
	assert( num >= 0 ); 
	assert( num <= m_numOfFrequenciesCalculatedByThisPE + 1 );

	return m_numObservedDataThisPEAccumulated[num];

}

// Increment number of interpolator vectors
int ObservedData::incrementNumInterpolatorVectors(){
	return m_numInterpolatorVectors++;
}

// Calculate frequencies treated by thie PE
void ObservedData::calcFrequenciesCalculatedByThisPE(){

	const int nfreq = getTotalNumberOfDifferenetFrequencies();
	//if( nfreq <= 0 ){
	//	OutputFiles::m_logFile << "Error : Total number of frequencies is less than zero !!" << std::endl;
	//	exit(1);
	//}
	assert( nfreq > 0 );

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	const int totalPE = ptrAnalysisControl->getTotalPE();
	const int myPE = ptrAnalysisControl->getMyPE();

	const int quotient = nfreq / totalPE;
	const int remainder = nfreq % totalPE;

//	if( myPE < remainder ){
//		m_numOfFrequenciesCalculatedByThisPE = quotient + 1;
//	}else{
//		m_numOfFrequenciesCalculatedByThisPE = quotient;
//	}
//	m_IDsOfFrequenciesCalculatedByThisPE = new int[m_numOfFrequenciesCalculatedByThisPE];
//	m_valuesOfFrequenciesCalculatedByThisPE = new double[m_numOfFrequenciesCalculatedByThisPE];
//
//#ifdef _DEBUG_WRITE
//	std::cout << "m_numOfFrequenciesCalculatedByThisPE = " << m_numOfFrequenciesCalculatedByThisPE << std::endl; // For debug
//#endif
//
//	// Calculate array containing IDs and values of frequencies caluculated by this PE
//	int icount(0);
//	for( int ifreq = myPE; ifreq < nfreq; ifreq += totalPE ){
//
//		const double frquencyValue = getFrequencyValue(ifreq);
//
//		m_IDsOfFrequenciesCalculatedByThisPE[icount] = ifreq;
//		m_valuesOfFrequenciesCalculatedByThisPE[icount] = frquencyValue;
//		
//		++icount;
//	}
//
//	if( m_numObservedDataThisPEAccumulated != NULL ){
//		delete [] m_numObservedDataThisPEAccumulated;
//		m_numObservedDataThisPEAccumulated = NULL;
//	}
//
//	m_numObservedDataThisPEAccumulated = new int[nfreq+1];
//	for( int i = 0; i < nfreq + 1; ++i ){
//		m_numObservedDataThisPEAccumulated[i] = 0;
//	}


	int* index = new int[ totalPE + 1 ];
	index[0] = 0;
	for( int i = 0; i < totalPE; ++i ){
		if( i < totalPE - remainder ){
			index[i+1] = quotient;
		}else{
			index[i+1] = quotient + 1;
		}
	}
	for( int i = 0; i < totalPE; ++i ){
		index[i+1] += index[i];
	}

	m_numOfFrequenciesCalculatedByThisPE = index[myPE+1] - index[myPE];

#ifdef _DEBUG_WRITE
	for( int i = 0; i < totalPE + 1; ++i ){
		std::cout << "i index = " << i << " " << index[i] << std::endl; // For debug
	}
	std::cout << "m_numOfFrequenciesCalculatedByThisPE = " << m_numOfFrequenciesCalculatedByThisPE << std::endl; // For debug

#endif

	m_IDsOfFrequenciesCalculatedByThisPE = new int[m_numOfFrequenciesCalculatedByThisPE];
	m_valuesOfFrequenciesCalculatedByThisPE = new double[m_numOfFrequenciesCalculatedByThisPE];

	// Calculate array containing IDs and values of frequencies caluculated by this PE
	int icount(0);
	for( int ifreq = index[myPE]; ifreq < index[myPE+1]; ++ifreq ){

		const double frquencyValue = getFrequencyValue(ifreq);

		m_IDsOfFrequenciesCalculatedByThisPE[icount] = ifreq;
		m_valuesOfFrequenciesCalculatedByThisPE[icount] = frquencyValue;
		
		++icount;
	}

	delete [] index;

	if( m_numObservedDataThisPEAccumulated != NULL ){
		delete [] m_numObservedDataThisPEAccumulated;
		m_numObservedDataThisPEAccumulated = NULL;
	}

	m_numObservedDataThisPEAccumulated = new int[nfreq+1];
	for( int i = 0; i < nfreq + 1; ++i ){
		m_numObservedDataThisPEAccumulated[i] = 0;
	}

}

// Get Total number of frequencies calculated by this PE
int ObservedData::getNumOfFrequenciesCalculatedByThisPE() const{
	return m_numOfFrequenciesCalculatedByThisPE;
}
		
// Get IDs of Frequencies calculated by this PE
int ObservedData::getIDsOfFrequenciesCalculatedByThisPE( const int num ) const{
	
	//if( num >= m_numOfFrequenciesCalculatedByThisPE || num < 0 ){
	//	OutputFiles::m_logFile << "Error : num is out range of getIDsOfFrequenciesCalculatedByThisPE !!" << std::endl;
	//	exit(1);
	//}
	assert( num >= 0 );
	assert( num < m_numOfFrequenciesCalculatedByThisPE );

	return m_IDsOfFrequenciesCalculatedByThisPE[num];

}

// Get values of Frequencies calculated by this PE
double ObservedData::getValuesOfFrequenciesCalculatedByThisPE( const int num ) const{

	//if( num >= m_numOfFrequenciesCalculatedByThisPE || num < 0 ){
	//	OutputFiles::m_logFile << "Error : num is out range of m_valuesOfFrequenciesCalculatedByThisPE !!" << std::endl;
	//	exit(1);
	//}
	assert( num >= 0 );
	assert( num < m_numOfFrequenciesCalculatedByThisPE );

	return m_valuesOfFrequenciesCalculatedByThisPE[num];
}

// Calculate error sum of squares
double ObservedData::calculateErrorSumOfSquaresThisPE() const{

	//const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	//if( ptrAnalysisControl->getTotalPE() > 1 ){
	//	OutputFiles::m_logFile << "calculateErrorSumOfSquaresThisPE() can be used only for serial calculation !!" << std::endl;
	//	exit(1);
	//} 

	double misfit(0.0);

	for ( int i = 0; i < m_numStationsMT; ++i ){
		misfit += m_observedStationMT[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		misfit += m_observedStationApparentResistivityAndPhase[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		misfit += m_observedStationHTF[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		misfit += m_observedStationVTF[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		misfit += m_observedStationPT[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		misfit += m_observedStationNMT[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		misfit += m_observedStationNMT2[i].calculateErrorSumOfSquaresThisPE();
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		misfit += m_observedStationNMT2ApparentResistivityAndPhase[i].calculateErrorSumOfSquaresThisPE();
	}

	return misfit;

}

// Calculate sum of squares of distortion matrix 
double ObservedData::calculateSumSquareOfDistortionMatrixComplexity() const{

	assert( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE );

	double norm(0.0);

	for( int i = 0; i < m_numStationsMT; ++i ){
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXX) , 2.0 );
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXY) , 2.0 );
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYX) , 2.0 );
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYY) , 2.0 );
	}

	for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXX) , 2.0 );
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXY) , 2.0 );
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYX) , 2.0 );
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYY) , 2.0 );
	}

#ifdef _DEBUG_WRITE
	std::cout << "Dist_NORM " << norm << std::endl;
#endif

	return norm;

}

// Calculate sum of squares of distortion matrix gains
double ObservedData::calculateSumSquareOfDistortionMatrixGains() const{

	assert( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS || 
			(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY );

	double norm(0.0);

	for( int i = 0; i < m_numStationsMT; ++i ){
		if( m_observedStationMT[i].getIDOfDistortionParams(ObservedDataStationMT::EX_GAIN) >= 0 ){
			norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EX_GAIN) , 2.0 );
		}
		if( m_observedStationMT[i].getIDOfDistortionParams(ObservedDataStationMT::EY_GAIN) >= 0 ){
			norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EY_GAIN) , 2.0 );
		}
	}

	for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		if( m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(ObservedDataStationMT::EX_GAIN) >= 0 ){
			norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_GAIN) , 2.0 );
		}
		if( m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(ObservedDataStationMT::EY_GAIN) >= 0 ){
			norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_GAIN) , 2.0 );
		}
	}

	return norm;

}

// Calculate sum of squares of distortion matrix rotations
double ObservedData::calculateSumSquareOfDistortionMatrixRotations() const{

	assert( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS );

	double norm(0.0);

	for( int i = 0; i < m_numStationsMT; ++i ){
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION) , 2.0 );
		norm += pow( m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION) , 2.0 );
	}

	for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION) , 2.0 );
		norm += pow( m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION) , 2.0 );
	}

	return norm;

}

// Get total number of the distortion parameters whose value is not fixed
int ObservedData::getNumDistortionParamsNotFixed() const{

	return static_cast<int>( m_typesOfDistortionParamsNotFixed.size() );

}

// Get types of distortion parameters whose value is not fixed
int ObservedData::getTypesOfDistortionParamsNotFixed( const int iParamsNotFixed ) const{

	return m_typesOfDistortionParamsNotFixed[iParamsNotFixed];

}

// Copy the distortion parameters not fixed to vector
void ObservedData::copyDistortionParamsNotFixedToVector( double* vector ) const{

	if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){

		for ( int i = 0; i < m_numStationsMT; ++i ){
			for( int j = 0; j < 4; ++j ){
				const int ID = m_observedStationMT[i].getIDOfDistortionParams(j);
				if( ID >= 0 ){
					vector[ ID ] = m_observedStationMT[i].getDistortionParams(j);
				}
			} 
		}

		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			for( int j = 0; j < 4; ++j ){
				const int ID = m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(j);
				if( ID >= 0 ){
					vector[ ID ] = m_observedStationApparentResistivityAndPhase[i].getDistortionParams(j);
				}
			} 
		}
	}

}

// Copy previous distortion parameters not fixed to vector
void ObservedData::copyDistortionParamsNotFixedPreToVector( double* vector ) const{

	if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			for( int j = 0; j < 4; ++j ){
				const int ID = m_observedStationMT[i].getIDOfDistortionParams(j);
				if( ID >= 0 ){
					vector[ ID ] = m_observedStationMT[i].getDistortionParamsPre(j);
				}
			} 
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			for( int j = 0; j < 4; ++j ){
				const int ID = m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(j);
				if( ID >= 0 ){
					vector[ ID ] = m_observedStationApparentResistivityAndPhase[i].getDistortionParamsPre(j);
				}
			} 
		}
	}

}

// Copy current distortion parameters to previous values
void ObservedData::copyDistortionParamsCurToPre() const{

	if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			for( int iComp = 0; iComp < 4; ++iComp ){
				const int index = m_observedStationMT[i].getIDOfDistortionParams(iComp);
				if( index >= 0 ){
					m_observedStationMT[i].copyDistortionParamsCurToPre(iComp); 
				}
			}
		}
		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			for( int iComp = 0; iComp < 4; ++iComp ){
				const int index = m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(iComp);
				if( index >= 0 ){
					m_observedStationApparentResistivityAndPhase[i].copyDistortionParamsCurToPre(iComp); 
				}
			}
		}
	}

}

// Calculate full updated values of distortion parameters
void ObservedData::calcDistortionParamsUpdatedFullFromIncrements( const double* const distortionParamIncre ){

	const double epsForDYdBeta = 0.01;
	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if( ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ||
		ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
		ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){
		for ( int i = 0; i < m_numStationsMT; ++i ){
			for( int iComp = 0; iComp < 4; ++iComp ){
				const int index = m_observedStationMT[i].getIDOfDistortionParams(iComp);
				if( index >= 0 ){
					if( ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS &&
						( iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION ) ){
						// Rotation angle is bounded in from -pi/2 to pi/2 radians
						// y = tan(beta)
						// beta = Arctan(y)
						// dy/dbeta = 1 / cos(beta)^2
						const double betaPre = m_observedStationMT[i].getDistortionParamsPre(iComp);
						const double dYdBeta = 1.0 / std::max( cos(betaPre) * cos(betaPre), epsForDYdBeta );
						const double newY = tan(betaPre) + dYdBeta * distortionParamIncre[index];
						const double newBeta = atan(newY);
						m_observedStationMT[i].setDistortionParamsUpdatedFull( iComp, newBeta );
					}else{
						const double newDistortionParam = m_observedStationMT[i].getDistortionParamsPre(iComp) + distortionParamIncre[index];
						m_observedStationMT[i].setDistortionParamsUpdatedFull( iComp, newDistortionParam ); 
					}
				}
			}
		}

		for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			for( int iComp = 0; iComp < 4; ++iComp ){
				const int index = m_observedStationApparentResistivityAndPhase[i].getIDOfDistortionParams(iComp);
				if( index >= 0 ){
					if( ptrAnalysisControl->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS &&
						( iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION ) ){
						// Rotation angle is bounded in from -pi/2 to pi/2 radians
						// y = tan(beta)
						// beta = Arctan(y)
						// dy/dbeta = 1 / cos(beta)^2
						const double betaPre = m_observedStationApparentResistivityAndPhase[i].getDistortionParamsPre(iComp);
						const double dYdBeta = 1.0 / std::max( cos(betaPre) * cos(betaPre), epsForDYdBeta );
						const double newY = tan(betaPre) + dYdBeta * distortionParamIncre[index];
						const double newBeta = atan(newY);
						m_observedStationApparentResistivityAndPhase[i].setDistortionParamsUpdatedFull( iComp, newBeta ); 
					}else{
						const double newDistortionParam = m_observedStationApparentResistivityAndPhase[i].getDistortionParamsPre(iComp) + distortionParamIncre[index];
						m_observedStationApparentResistivityAndPhase[i].setDistortionParamsUpdatedFull( iComp, newDistortionParam ); 
					}
				}
			}
		}
	}

}

// Update distortion parameters
void ObservedData::updateDistortionParams(){

	const double stepLengthDampingFactor = (  AnalysisControl::getInstance() )->getStepLengthDampingFactorCur();

	for ( int i = 0; i < m_numStationsMT; ++i ){
		m_observedStationMT[i].updateDistortionParams( stepLengthDampingFactor );
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		m_observedStationApparentResistivityAndPhase[i].updateDistortionParams( stepLengthDampingFactor );
	}

}

// Output distortion parametersto file
void ObservedData::outputDistortionParams( const int iterNum ) const{

	std::ostringstream fileName;
	fileName << "distortion_iter" << iterNum << ".dat";

	FILE *fp;
	if( (fp = fopen( fileName.str().c_str(), "w")) == NULL ) {
		OutputFiles::m_logFile  << "File open error !! : " << fileName.str() << std::endl;
		exit(1);
	}

	fprintf(fp, "%10d\n", m_numStationsMT + m_numStationsApparentResistivityAndPhase );

	if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE ){

		for( int i = 0; i < m_numStationsMT; ++i ){
			fprintf(fp, "%10d%15lf%15lf%15lf%15lf%10d\n",  m_observedStationMT[i].getStationID(),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXX),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXY),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYX),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYY),
				m_observedStationMT[i].doesFixDistortionMatrix() ? 1 : 0 );
		}
		for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			fprintf(fp, "%10d%15lf%15lf%15lf%15lf%10d\n",  m_observedStationApparentResistivityAndPhase[i].getStationID(),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXX),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CXY),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYX),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::COMPONENT_ID_CYY),
				m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ? 1 : 0 );
		}

	}
	else if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ){

		for( int i = 0; i < m_numStationsMT; ++i ){
			fprintf(fp, "%10d%15lf%15lf%15lf%15lf%10d\n",  m_observedStationMT[i].getStationID(),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EX_GAIN),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EY_GAIN),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION) * CommonParameters::rad2deg,
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION) * CommonParameters::rad2deg,
				m_observedStationMT[i].doesFixDistortionMatrix() ? 1 : 0 );
		}
		for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			fprintf(fp, "%10d%15lf%15lf%15lf%15lf%10d\n",  m_observedStationApparentResistivityAndPhase[i].getStationID(),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_GAIN),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_GAIN),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_ROTATION) * CommonParameters::rad2deg,
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_ROTATION) * CommonParameters::rad2deg,
				m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ? 1 : 0 );
		}

	}
	else if( (AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY ){

		for( int i = 0; i < m_numStationsMT; ++i ){
			fprintf(fp, "%10d%15lf%15lf%10d\n",  m_observedStationMT[i].getStationID(),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EX_GAIN),
				m_observedStationMT[i].getDistortionParams(ObservedDataStationMT::EY_GAIN),
				m_observedStationMT[i].doesFixDistortionMatrix() ? 1 : 0 );
		}
		for( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
			fprintf(fp, "%10d%15lf%15lf%10d\n",  m_observedStationApparentResistivityAndPhase[i].getStationID(),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EX_GAIN),
				m_observedStationApparentResistivityAndPhase[i].getDistortionParams(ObservedDataStationMT::EY_GAIN),
				m_observedStationApparentResistivityAndPhase[i].doesFixDistortionMatrix() ? 1 : 0 );
		}

	}
	else{
		OutputFiles::m_logFile << "Error : Wrong type of distortion matrix !!" << std::endl;
		exit(1);
	}
	
	fclose(fp);

}

// Output information of locations of observed stations to vtk file
void ObservedData::outputLocationsOfObservedStationsToVtk() const{

	if( !OutputFiles::m_vtkFileForObservedStation.is_open() ){
		return;
	}

	const int numPoints = m_numStationsMT + m_numStationsApparentResistivityAndPhase + m_numStationsHTF + m_numStationsVTF 
		+ m_numStationsPT + m_numStationsNMT * 2 + m_numStationsNMT2 * 4 + m_numStationsNMT2ApparentResistivityAndPhase * 4;

	OutputFiles::m_vtkFileForObservedStation << "DATASET UNSTRUCTURED_GRID" << std::endl;
	OutputFiles::m_vtkFileForObservedStation << "POINTS " << numPoints << " float" << std::endl;

	OutputFiles::m_vtkFileForObservedStation.precision(6);
	for ( int i = 0; i < m_numStationsMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationMT[i].getLocationOfPoint() ).X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationMT[i].getLocationOfPoint() ).Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific <<  m_observedStationMT[i].getZCoordOfPoint() << std::endl;
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationApparentResistivityAndPhase[i].getLocationOfPoint() ).X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationApparentResistivityAndPhase[i].getLocationOfPoint() ).Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific <<  m_observedStationApparentResistivityAndPhase[i].getZCoordOfPoint() << std::endl;
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationHTF[i].getLocationOfPoint() ).X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationHTF[i].getLocationOfPoint() ).Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific <<  m_observedStationHTF[i].getZCoordOfPoint() << std::endl;
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationVTF[i].getLocationOfPoint() ).X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationVTF[i].getLocationOfPoint() ).Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific <<  m_observedStationVTF[i].getZCoordOfPoint() << std::endl;
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationPT[i].getLocationOfPoint() ).X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationPT[i].getLocationOfPoint() ).Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific <<  m_observedStationPT[i].getZCoordOfPoint() << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT[i].getLocationOfStation() ).startPoint.X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT[i].getLocationOfStation() ).startPoint.Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT[i].getZCoordOfPoint(0) << std::endl;
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT[i].getLocationOfStation() ).endPoint.X << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT[i].getLocationOfStation() ).endPoint.Y << " ";
		OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT[i].getZCoordOfPoint(1) << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		for ( int iDip = 0; iDip < 2; ++iDip ){ 
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2[i].getLocationOfStation( iDip ) ).startPoint.X << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2[i].getLocationOfStation( iDip ) ).startPoint.Y << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT2[i].getZCoordOfPoint(iDip, 0) << std::endl;
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2[i].getLocationOfStation( iDip ) ).endPoint.X << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2[i].getLocationOfStation( iDip ) ).endPoint.Y << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT2[i].getZCoordOfPoint(iDip, 1) << std::endl;
		}
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		for ( int iDip = 0; iDip < 2; ++iDip ){ 
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2ApparentResistivityAndPhase[i].getLocationOfStation( iDip ) ).startPoint.X << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2ApparentResistivityAndPhase[i].getLocationOfStation( iDip ) ).startPoint.Y << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT2ApparentResistivityAndPhase[i].getZCoordOfPoint(iDip, 0) << std::endl;
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2ApparentResistivityAndPhase[i].getLocationOfStation( iDip ) ).endPoint.X << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << ( m_observedStationNMT2ApparentResistivityAndPhase[i].getLocationOfStation( iDip ) ).endPoint.Y << " ";
			OutputFiles::m_vtkFileForObservedStation << std::setw(15) << std::scientific << m_observedStationNMT2ApparentResistivityAndPhase[i].getZCoordOfPoint(iDip, 1) << std::endl;
		}
	}

	const int numCells  = m_numStationsMT + m_numStationsApparentResistivityAndPhase + m_numStationsHTF + m_numStationsVTF
		+ m_numStationsPT + m_numStationsNMT + m_numStationsNMT2 * 2 + m_numStationsNMT2ApparentResistivityAndPhase * 2;
	const int sizeOfCellData = numCells + numPoints;

	OutputFiles::m_vtkFileForObservedStation << "CELLS " << numCells << " " << sizeOfCellData << std::endl;

	int icount(0);
	for ( int i = 0; i < m_numStationsMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 2 << " " << icount++ << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 2 << " " << icount++ << " " << icount++ << std::endl;
		OutputFiles::m_vtkFileForObservedStation << 2 << " " << icount++ << " " << icount++ << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 2 << " " << icount++ << " " << icount++ << std::endl;
		OutputFiles::m_vtkFileForObservedStation << 2 << " " << icount++ << " " << icount++ << std::endl;
	}

	OutputFiles::m_vtkFileForObservedStation << "CELL_TYPES " << numCells << std::endl;
	for ( int i = 0; i < m_numStationsMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << std::endl;
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << std::endl;
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << std::endl;
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << std::endl;
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 1 << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 3 << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 3 << std::endl;
		OutputFiles::m_vtkFileForObservedStation << 3 << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << 3 << std::endl;
		OutputFiles::m_vtkFileForObservedStation << 3 << std::endl;
	}

	OutputFiles::m_vtkFileForObservedStation << "CELL_DATA " << numCells << std::endl;

	OutputFiles::m_vtkFileForObservedStation << "SCALARS station_ID int" << std::endl;
	OutputFiles::m_vtkFileForObservedStation << "LOOKUP_TABLE default" <<  std::endl;
	for ( int i = 0; i < m_numStationsMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationMT[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationApparentResistivityAndPhase[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationHTF[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationVTF[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationPT[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationNMT[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationNMT2[i].getStationID() << std::endl;
		OutputFiles::m_vtkFileForObservedStation << m_observedStationNMT2[i].getStationID() << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationNMT2ApparentResistivityAndPhase[i].getStationID() << std::endl;
		OutputFiles::m_vtkFileForObservedStation << m_observedStationNMT2ApparentResistivityAndPhase[i].getStationID() << std::endl;
	}

	OutputFiles::m_vtkFileForObservedStation << "SCALARS station_type int" << std::endl;
	OutputFiles::m_vtkFileForObservedStation << "LOOKUP_TABLE default" <<  std::endl;

	for ( int i = 0; i < m_numStationsMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::MT << std::endl;
	}

	for ( int i = 0; i < m_numStationsApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::APP_RES_AND_PHS << std::endl;
	}

	for ( int i = 0; i < m_numStationsHTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::HTF << std::endl;
	}

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::VTF << std::endl;
	}

	for ( int i = 0; i < m_numStationsPT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::PT << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::NMT << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::NMT2 << std::endl;
		OutputFiles::m_vtkFileForObservedStation << ObservedData::NMT2 << std::endl;
	}

	for ( int i = 0; i < m_numStationsNMT2ApparentResistivityAndPhase; ++i ){
		OutputFiles::m_vtkFileForObservedStation << ObservedData::NMT2_APP_RES_AND_PHS << std::endl;
		OutputFiles::m_vtkFileForObservedStation << ObservedData::NMT2_APP_RES_AND_PHS<< std::endl;
	}

	OutputFiles::m_vtkFileForObservedStation.close();
}

// Output induction arrow to vtk file
void ObservedData::outputInductionArrowToVtk( const int iterNum ) const{

	if( m_numStationsVTF <= 0 ){
		return;
	}

	// Get process ID and total process number
	int myPE(0);
	MPI_Comm_rank ( MPI_COMM_WORLD, &myPE );

	// Open VTK file
	std::ostringstream vtkFileName;
	vtkFileName << "induction_arrow_" << myPE << "_iter" << iterNum << ".vtk";

	std::ofstream vtkFileForInductionArrow( vtkFileName.str().c_str(), std::ios::out );
	if( vtkFileForInductionArrow.fail() )
	{
		OutputFiles::m_logFile << "File open error !! : " << vtkFileName.str() << std::endl;
		exit(1);
	}
	vtkFileForInductionArrow << "# vtk DataFile Version 2.0" << std::endl;
	vtkFileForInductionArrow << "Induction_Arrow_Iter" <<  iterNum << "_PE" << myPE << std::endl;
	vtkFileForInductionArrow << "ASCII" << std::endl;

	//-----------------------
	//--- Output location ---
	//-----------------------
	// Points
	const int numPoints = m_numStationsVTF;
	
	vtkFileForInductionArrow << "DATASET UNSTRUCTURED_GRID" << std::endl;
	vtkFileForInductionArrow << "POINTS " << numPoints << " float" << std::endl;

	for ( int i = 0; i < m_numStationsVTF; ++i ){
		vtkFileForInductionArrow << std::setw(15) << std::scientific << ( m_observedStationVTF[i].getLocationOfPoint() ).X << " ";
		vtkFileForInductionArrow << std::setw(15) << std::scientific << ( m_observedStationVTF[i].getLocationOfPoint() ).Y << " ";
		vtkFileForInductionArrow << std::setw(15) << std::scientific <<  m_observedStationVTF[i].getZCoordOfPoint() << std::endl;
	}

	// Cells
	const int numCells  = m_numStationsVTF;
	const int sizeOfCellData = numCells + numPoints;

	vtkFileForInductionArrow << "CELLS " << numCells << " " << sizeOfCellData << std::endl;

	int icount(0);
	for ( int i = 0; i < m_numStationsVTF; ++i ){
		vtkFileForInductionArrow << 1 << " " << icount++ << std::endl;
	}

	vtkFileForInductionArrow << "CELL_TYPES " << numCells << std::endl;
	for ( int i = 0; i < m_numStationsVTF; ++i ){
		vtkFileForInductionArrow << 1 << std::endl;
	}

	OutputFiles::m_vtkFileForObservedStation << "CELL_DATA " << numCells << std::endl;

	OutputFiles::m_vtkFileForObservedStation << "SCALARS Station_ID int" << std::endl;
	OutputFiles::m_vtkFileForObservedStation << "LOOKUP_TABLE default" <<  std::endl;
	for ( int i = 0; i < m_numStationsVTF; ++i ){
		OutputFiles::m_vtkFileForObservedStation << m_observedStationVTF[i].getStationID() << std::endl;
	}

	//-------------------------------
	//--- Output induction arrows ---
	//-------------------------------
	vtkFileForInductionArrow << "POINT_DATA " << numPoints << std::endl;

	const int nfreq = getNumOfFrequenciesCalculatedByThisPE();
	for( int ifreq = 0; ifreq < nfreq; ++ifreq ){

		const double freq = getValuesOfFrequenciesCalculatedByThisPE( ifreq );

		//--- Real Arrow -----
		vtkFileForInductionArrow << "VECTORS Real_Arrow_" << freq <<"(Hz) float" << std::endl;
		for ( int i = 0; i < m_numStationsVTF; ++i ){
			std::complex<double> Tzx(0.0,0.0);
			std::complex<double> Tzy(0.0,0.0);

			if( m_observedStationVTF[i].getVTF( freq, Tzx, Tzy ) ){
				vtkFileForInductionArrow << std::setw(15) << - Tzx.real() << std::setw(15) << - Tzy.real() << std::setw(15) << 0.0 << std::endl;
			}
			else{
				vtkFileForInductionArrow << std::setw(15) << 0.0 << std::setw(15) << 0.0 << std::setw(15) << 0.0 << std::endl;
			}
		}
		//---------------------

		//--- Imaginary Arrow -----
		vtkFileForInductionArrow << "VECTORS Imag_Arrow_" << freq <<"(Hz) float" << std::endl;
		for ( int i = 0; i < m_numStationsVTF; ++i ){
			std::complex<double> Tzx(0.0,0.0);
			std::complex<double> Tzy(0.0,0.0);
			if( m_observedStationVTF[i].getVTF( freq, Tzx, Tzy ) ){
				vtkFileForInductionArrow << std::setw(15) << - Tzx.imag() << std::setw(15) << - Tzy.imag() << std::setw(15) << 0.0 << std::endl;
			}
			else{
				vtkFileForInductionArrow << std::setw(15) << 0.0 << std::setw(15) << 0.0 << std::setw(15) << 0.0 << std::endl;
			}
		}
		//---------------------

	}

	vtkFileForInductionArrow.close();

}

// Find and return frequency IDs from frequency value
int  ObservedData::getFreqIDs( const double freq ) const{

	const double EPS = 1.0E-8;

	int icount(0);
	for( std::vector<double>::const_iterator itr = m_frequencyAll.begin(); itr != m_frequencyAll.end(); ++itr ){
		if( std::fabs( freq - *itr ) < EPS ){
			return icount;
		}
		++icount;
	}

	return -1;

}
