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
#include "ResistivityBlock.h"
#include "MeshDataBrickElement.h"
#include "MeshDataNonConformingHexaElement.h"
#include "OutputFiles.h"
#include "ObservedData.h"
#ifdef _ANISOTOROPY
#include "Util.h"
#endif
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <algorithm>

// Return the the instance of the class
ResistivityBlock* ResistivityBlock::getInstance(){
   	static ResistivityBlock instance;// The only instance
  	return &instance;
}

// Constructer
ResistivityBlock::ResistivityBlock():
	m_elementID2blockID(NULL),
	m_blockID2modelID(NULL),
	m_modelID2blockID(NULL),
	m_numResistivityBlockTotal(0),
	m_numResistivityBlockNotFixed(0),
	m_resistivityValues(NULL),
	m_resistivityValuesPre(NULL),
	m_resistivityValuesUpdatedFull(NULL),
	m_resistivityValuesMin(NULL),
	m_resistivityValuesMax(NULL),
	m_weightingConstants(NULL),
	m_fixResistivityValues(NULL),
	m_isolated(NULL),
	m_rougheningMatrix(),
	m_includeBottomResistivity(false),
	m_blockID2Elements(NULL),
	m_bottomResistivity(1.0),
	m_roughningFactorAtBottom(1.0),
	m_addSmallValueToDiagonals(false),
	m_smallValueAddedToDiagonals(0.0),
	m_minDistanceToBounds(0.01),
	m_inverseDistanceWeightingFactor(1.0),
	m_typeBoundConstraints(ResistivityBlock::SIMPLE_BOUND_CONSTRAINING)
{}

// Destructer
ResistivityBlock::~ResistivityBlock(){

	if( m_elementID2blockID != NULL){
		delete[] m_elementID2blockID;
		m_elementID2blockID = NULL;
	}

	if( m_blockID2modelID != NULL){
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}

	if( m_modelID2blockID != NULL){
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}

	if( m_resistivityValues != NULL){
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}

	if( m_resistivityValuesPre != NULL){
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}

	if( m_resistivityValuesUpdatedFull != NULL){
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}

	if( m_resistivityValuesMin != NULL){
		delete[] m_resistivityValuesMin;
		m_resistivityValuesMin = NULL;
	}

	if( m_resistivityValuesMax != NULL){
		delete[] m_resistivityValuesMax;
		m_resistivityValuesMax = NULL;
	}

	if( m_weightingConstants != NULL){
		delete[] m_weightingConstants;
		m_weightingConstants = NULL;
	}

	if( m_fixResistivityValues != NULL){
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}

	if( m_isolated != NULL){
		delete[] m_isolated;
		m_isolated = NULL;
	}

	if( m_blockID2Elements != NULL){
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}

	m_rougheningMatrix.releaseMemory();

}

// Read data of resisitivity block model from input file
void ResistivityBlock::inputResisitivityBlock(){

	const int iterInit = ( AnalysisControl::getInstance() )->getIterationNumInit();
	std::ostringstream inputFile;
	inputFile << "resistivity_block_iter" << iterInit << ".dat";
	std::ifstream inFile( inputFile.str().c_str(), std::ios::in );

	if( inFile.fail() )
	{
		OutputFiles::m_logFile << "File open error : " << inputFile.str().c_str() << " !!" << std::endl;
		exit(1);
	}

	int nElem(0);
	inFile >> nElem;
	if( m_elementID2blockID != NULL ){
		delete[] m_elementID2blockID;
		m_elementID2blockID = NULL;
	}
	m_elementID2blockID = new int[ nElem ];

	int nBlk(0);
	inFile >> nBlk;	
	m_numResistivityBlockTotal = nBlk;

	if( m_blockID2modelID != NULL){
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}
	m_blockID2modelID = new int[ m_numResistivityBlockTotal ];

	if( m_resistivityValues != NULL ){
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}
	m_resistivityValues = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesPre != NULL){
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}
	m_resistivityValuesPre = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesUpdatedFull != NULL){
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}
	m_resistivityValuesUpdatedFull = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesMin != NULL){
		delete[] m_resistivityValuesMin;
		m_resistivityValuesMin = NULL;
	}
	m_resistivityValuesMin = new double[ m_numResistivityBlockTotal ];

	if( m_resistivityValuesMax != NULL){
		delete[] m_resistivityValuesMax;
		m_resistivityValuesMax = NULL;
	}
	m_resistivityValuesMax = new double[ m_numResistivityBlockTotal ];

	if( m_weightingConstants != NULL){
		delete[] m_weightingConstants;
		m_weightingConstants = NULL;
	}
	m_weightingConstants = new double[ m_numResistivityBlockTotal ];

	if( m_fixResistivityValues != NULL ){
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}
	m_fixResistivityValues = new bool[ m_numResistivityBlockTotal ];

	if( m_isolated != NULL){
		delete[] m_isolated;
		m_isolated = NULL;
	}
	m_isolated = new bool[ m_numResistivityBlockTotal ];

	for( int i = 0; i < m_numResistivityBlockTotal; ++i ){
		m_resistivityValues[i] = 0.0;
		m_resistivityValuesPre[i] = 0.0;
		m_resistivityValuesUpdatedFull[i] = 0.0;
		m_resistivityValuesMin[i] = 0.0;
		m_resistivityValuesMax[i] = 0.0;
		m_weightingConstants[i] = 1.0;
		m_fixResistivityValues[i] = false;
		m_isolated[i] = false;
	}

#ifdef _DEBUG_WRITE
	std::cout << nElem << " " << m_numResistivityBlockTotal << std::endl; // For debug
#endif
	
	for( int iElem = 0; iElem < nElem; ++iElem ){
		int idum(0);
		int iblk(0);// Resistivity block ID
		inFile >> idum >> iblk;
		m_elementID2blockID[ iElem ] = iblk;

		if( iblk >= m_numResistivityBlockTotal || iblk < 0 )	{
			OutputFiles::m_logFile << "Error : Resistivity block ID " << iblk << " of element " << iElem << " is improper !!" << std::endl;
			exit(1);
		}		

#ifdef _DEBUG_WRITE
		std::cout << iElem << " " << iblk << std::endl; // For debug
#endif
	}
	
	const bool dataSpaceInversionMethod = ( ( AnalysisControl::getInstance() )->getInversionMethod() == Inversion::GAUSS_NEWTON_DATA_SPECE );
	m_numResistivityBlockNotFixed = 0;
#ifdef _ANISOTOROPY
	int counterAnisotropicBlock(0);
#endif
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		int idum(0);
		int itype(0);

#ifdef _ANISOTOROPY
		inFile >> idum;
		if( idum != iBlk ){
			OutputFiles::m_logFile << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		bool anisotropy(false);
		if( ( AnalysisControl::getInstance() )->isAnisotropyConsidered() ){
			inFile >> idum;
			if( idum != 0 ){
				anisotropy = true;
			}
		}
		if( anisotropy ){
			m_mapBlockIDWithAnisotropyToIndex.insert( std::make_pair(iBlk, counterAnisotropicBlock) );
			++counterAnisotropicBlock;
			m_resistivityValues[iBlk] = -9999.999;
			CommonParameters::Vector3D resistivity;
			inFile >> resistivity.X >> resistivity.Y >> resistivity.Z;
			if( resistivity.X <= 0.0 ){
				OutputFiles::m_logFile << "Error : Resistivity component XX of block " << iBlk << " is less than or equal to zero !! : " << resistivity.X << std::endl;
				exit(1);
			}
			if( resistivity.Y <= 0.0 ){
				OutputFiles::m_logFile << "Error : Resistivity component YY of block " << iBlk << " is less than or equal to zero !! : " << resistivity.Y << std::endl;
				exit(1);
			}
			if( resistivity.Z <= 0.0 ){
				OutputFiles::m_logFile << "Error : Resistivity component ZZ of block " << iBlk << " is less than or equal to zero !! : " << resistivity.Z << std::endl;
				exit(1);
			}
			m_resistivityValuesAxialAnisotropy.push_back(resistivity);
			m_resistivityValuesAxialAnisotropyPre.push_back(resistivity);
			m_resistivityValuesAxialAnisotropyFull.push_back(resistivity);
			double strike(0.0);
			double dip(0.0);
			double slant(0.0);
			inFile >> strike >> dip >> slant;
			strike *= CommonParameters::deg2rad;
			dip *= CommonParameters::deg2rad;
			slant *= CommonParameters::deg2rad;
			m_axialAnisotropyStrileAngle.push_back(strike);
			m_axialAnisotropyDipAngle.push_back(dip);
			m_axialAnisotropySlantAngle.push_back(slant);
		}else{
			inFile >> m_resistivityValues[iBlk];
			if( m_resistivityValues[iBlk] <= 0.0 ){
				OutputFiles::m_logFile << "Error : Resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValues[iBlk] << std::endl;
				exit(1);
			}
		}
		inFile >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk] >> m_weightingConstants[iBlk] >> itype;
#else
		inFile >> idum >> m_resistivityValues[iBlk] >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk] >> m_weightingConstants[iBlk] >> itype;
		if( idum != iBlk ){
			OutputFiles::m_logFile << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		if( m_resistivityValues[iBlk] <= 0.0 ){
			OutputFiles::m_logFile << "Error : Resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValues[iBlk] << std::endl;
			exit(1);
		}
#endif
		if( m_resistivityValuesMin[iBlk] <= 0.0 ){
			OutputFiles::m_logFile << "Error : Minimum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMin[iBlk] << std::endl;
			exit(1);
		}
		if( m_resistivityValuesMax[iBlk] <= 0.0 ){
			OutputFiles::m_logFile << "Error : Maximum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMax[iBlk] << std::endl;
			exit(1);
		}
#ifdef _ANISOTOROPY
		if( anisotropy ){
			const CommonParameters::Vector3D resistivity = m_resistivityValuesAxialAnisotropy.back();
			if( m_resistivityValuesMax[iBlk] < resistivity.X ){
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity ( " << resistivity.X << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMax[iBlk] < resistivity.Y ){
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity ( " << resistivity.Y << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMax[iBlk] < resistivity.Z ){
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity ( " << resistivity.Z << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMin[iBlk] > resistivity.X ){
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity ( " << resistivity.X << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMin[iBlk] > resistivity.Y ){
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity ( " << resistivity.Y << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMin[iBlk] > resistivity.Z ){
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity ( " << resistivity.Z << " )." << std::endl;
				exit(1);
			}
		}else{
			if( m_resistivityValuesMax[iBlk] < m_resistivityValues[iBlk] ){
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
				exit(1);
			}
			if( m_resistivityValuesMin[iBlk] > m_resistivityValues[iBlk] ){
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
				exit(1);
			}
		}
#else
		if( m_resistivityValuesMax[iBlk] < m_resistivityValues[iBlk] ){
			OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
		if( m_resistivityValuesMin[iBlk] > m_resistivityValues[iBlk] ){
			OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
#endif
		if( m_weightingConstants[iBlk] <= 0.0 ){
			OutputFiles::m_logFile << "Error : Weighting constant of block " << iBlk << " is less than or equal to zero !! : " << m_weightingConstants[iBlk] << std::endl;
			exit(1);
		}
		if( dataSpaceInversionMethod && itype == FREE_AND_ISOLATED){
			OutputFiles::m_logFile << "Error : Resistivity of isolated block must be fixed when data space inverson method is selected !!" << std::endl;
			exit(1);
		}

		switch(itype){
			case ResistivityBlock::FREE_AND_CONSTRAINED:// Go through
			case ResistivityBlock::FREE_AND_ISOLATED:
				m_blockID2modelID[iBlk] = m_numResistivityBlockNotFixed++;
				break;
			case ResistivityBlock::FIXED_AND_ISOLATED:// Go through
			case ResistivityBlock::FIXED_AND_CONSTRAINED:
				m_fixResistivityValues[iBlk] = true;
				m_blockID2modelID[iBlk] = -1;
				break;
			default:
				OutputFiles::m_logFile << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
				exit(1);
				break;
		}

		switch(itype){
			case ResistivityBlock::FREE_AND_CONSTRAINED:// Go through
			case ResistivityBlock::FIXED_AND_CONSTRAINED:
				m_isolated[iBlk] = false;
				break;
			case ResistivityBlock::FIXED_AND_ISOLATED:// Go through
			case ResistivityBlock::FREE_AND_ISOLATED:
				m_isolated[iBlk] = true;
				break;
			default:
				OutputFiles::m_logFile << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
				exit(1);
				break;
		}

#ifdef _DEBUG_WRITE
		std::cout << std::setw(5) << iBlk << std::setw(15) << m_resistivityValues[iBlk] << std::setw(15)	<< m_resistivityValuesMin[iBlk]	<< std::setw(15) << m_resistivityValuesMax[iBlk]
			<< std::setw(15) << m_weightingConstants[iBlk] <<  std::setw(5) << m_fixResistivityValues[iBlk]  <<  std::setw(5) << m_blockID2modelID[iBlk] << std::endl; // For debug
#endif

	}

	if( !m_fixResistivityValues[0] ){
		OutputFiles::m_logFile << "Error : Resistivity block 0 must be the air. And, its resistivity must be fixed." << std::endl;
		exit(1);
	}
	
	inFile.close();

	memcpy( m_resistivityValuesPre, m_resistivityValues, sizeof(double)*(m_numResistivityBlockTotal) );

	if( m_modelID2blockID != NULL){
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}
	m_modelID2blockID = new int[ m_numResistivityBlockNotFixed ];

	int icount(0);
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){

		if( !m_fixResistivityValues[iBlk] ){
			m_modelID2blockID[icount] = iBlk;
			++icount;
		}

	}

	if( icount != m_numResistivityBlockNotFixed ){
		OutputFiles::m_logFile << "Error : icount is not equal to m_numResistivityBlockNotFixed. icount = " << icount << " m_numResistivityBlockNotFixed = " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	for( int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl ){
		std::cout << " iMdl m_modelID2blockID[iMdl] : " << iMdl << " " << m_modelID2blockID[iMdl] << std::endl;
	}
#endif

	if( m_blockID2Elements != NULL){
		delete[] m_blockID2Elements;
		m_blockID2Elements = NULL;
	}
	m_blockID2Elements= new std::vector< std::pair<int,double> >[m_numResistivityBlockTotal];

	for( int iElem = 0; iElem < nElem; ++iElem ){
		m_blockID2Elements[ m_elementID2blockID[ iElem ] ].push_back( std::make_pair(iElem,1.0) );
	}
#ifndef _LINUX
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		m_blockID2Elements[ iBlk ].shrink_to_fit();
	}
#endif

	if( m_numResistivityBlockNotFixed <= 0 ){
		OutputFiles::m_logFile << "Error : Total number of modifiable resisitivity value is zero or negative !! : " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		const int num = static_cast<int>( m_blockID2Elements[ iBlk ].size() ); 
		for( int i = 0; i < num; ++i ){
			std::cout << " m_blockID2Elements[ " << iBlk << " ][ " << i << "] : " << m_blockID2Elements[iBlk][i].first << std::endl;
		}
	}
#endif

}

// Get resistivity values from resisitivity block ID
double ResistivityBlock::getResistivityValuesFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValues[ iblk ];
}

// Get previous resistivity values from resisitivity block ID
double ResistivityBlock::getResistivityValuesPreFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValuesPre[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValuesPre[ iblk ];

}

//// Get resisitivity block ID from element ID
//inline int ResistivityBlock::getBlockIDFromElemID( const int ielem ) const{
//	return m_elementID2blockID[ ielem ];
//}

// Get conductivity values from resisitivity block ID
double ResistivityBlock::getConductivityValuesFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 0.0;
	}

	return 1.0/m_resistivityValues[ iblk ];
}

// Get resistivity values from element ID
double ResistivityBlock::getResistivityValuesFromElemID( const int ielem ) const{

	//const int iblk = m_elementID2blockID[ ielem ];
	const int iblk = getBlockIDFromElemID(ielem);
	return getResistivityValuesFromBlockID( iblk );
}

// Get conductivity values from element ID
double ResistivityBlock::getConductivityValuesFromElemID( const int ielem ) const{

	//const int iblk = m_elementID2blockID[ ielem ];
	const int iblk = getBlockIDFromElemID(ielem);
	return getConductivityValuesFromBlockID( iblk );
}

// Get model ID from block ID
int ResistivityBlock::getModelIDFromBlockID( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_blockID2modelID[ iblk ];
}

// Get block ID from model ID
int ResistivityBlock::getBlockIDFromModelID( const int imdl ) const{

	//if( imdl < 0 || imdl >= m_numResistivityBlockNotFixed ){
	//	OutputFiles::m_logFile << "Error : Specified model ID is out of range. imdl = " << imdl << std::endl;
	//	exit(1);
	//}
	assert( imdl >= 0 );
	assert( imdl < m_numResistivityBlockNotFixed );

	return m_modelID2blockID[ imdl ];

}

// Get total number of resistivity blocks
int ResistivityBlock::getNumResistivityBlockTotal() const{
	return m_numResistivityBlockTotal;
}

// Get number of resistivity blocks whose resistivity values are fixed
int ResistivityBlock::getNumResistivityBlockNotFixed() const{
	return m_numResistivityBlockNotFixed;
}

// Get flag specifing whether bottom resistivity is included in roughning
bool ResistivityBlock::includeBottomResistivity() const{
	return m_includeBottomResistivity;
}

// Get resistivity of the bottom of the model
double ResistivityBlock::getBottomResistivity() const{
	return m_bottomResistivity;
}

// Get roughning factor at the bottom of the model
double ResistivityBlock::getRoughningFactorAtBottom() const{
	return m_roughningFactorAtBottom;
}

// Get flag specifing whether small value is added to diagonals
bool ResistivityBlock::getFlagAddSmallValueToDiagonals() const{
	return m_addSmallValueToDiagonals;
}

// Set small value added to bottom
double ResistivityBlock::getSmallValueAddedToDiagonals() const{
	return m_smallValueAddedToDiagonals;
}

// Get minimum distance between current resistivity and resistivity bounds in common logarithm scale
double ResistivityBlock::getMinDistanceToBounds() const{
	return m_minDistanceToBounds;
}

// Get type of bound constraints
int ResistivityBlock::getTypeBoundConstraints() const{
	return m_typeBoundConstraints;
}

// Get positive real factor of inverse distance weighting
double ResistivityBlock::getInverseDistanceWeightingFactor() const{
	return m_inverseDistanceWeightingFactor;
}

// Set flag specifing whether bottom resistivity is included in roughning
void ResistivityBlock::setFlagIncludeBottomResistivity( bool include ){
	m_includeBottomResistivity = include;
}

// Set resistivity of the bottom of the model
void ResistivityBlock::setBottomResistivity(  const double resistivity ){
	m_bottomResistivity = resistivity;
}

// Set roughning factor at the bottom of the model
void ResistivityBlock::setRoughningFactorAtBottom( const double factor ){
	m_roughningFactorAtBottom = factor;
}

// Set flag specifing whether small value is added to diagonals
void ResistivityBlock::setFlagAddSmallValueToDiagonals( const bool flag ){
	m_addSmallValueToDiagonals = flag;
}

// Set small value added to bottom
void ResistivityBlock::setSmallValueAddedToDiagonals( const double value ){
	m_smallValueAddedToDiagonals = value;
}

// Set minimum distance to resistivity bounds in common logarithm scale
void ResistivityBlock::setMinDistanceToBounds( const double distance ){
	assert( distance > 0.0 );
	m_minDistanceToBounds = distance;
}

// Set type of bound constraints
void ResistivityBlock::setTypeBoundConstraints( const int type ){
	if( type != SIMPLE_BOUND_CONSTRAINING && type != TRANSFORMING_METHOD ){
		OutputFiles::m_logFile << "Error : Wrong type of bound constraining metohd !! : " << type << " ." << std::endl;
		exit(1);
	}
	m_typeBoundConstraints = type;
}

// Set positive real factor of inverse distance weighting
void ResistivityBlock::setInverseDistanceWeightingFactor( const double factor ){
	m_inverseDistanceWeightingFactor = factor;
}

// Get flag specifing whether resistivity value of each block is fixed or not
bool ResistivityBlock::isFixedResistivityValue( const int iblk ) const{

	//if( iblk < 0 || iblk >= m_numResistivityBlockTotal ){
	//	OutputFiles::m_logFile << "Error : Specified block ID is out of range. iblk = " << iblk << std::endl;
	//	exit(1);
	//}
	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_fixResistivityValues[iblk];
}

// Get flag specifing whether resistivity block is excluded from roughing matrix
bool ResistivityBlock::isolated( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_isolated[iblk];

}

// Calculate volume of the specified resistivity block
double ResistivityBlock::calcVolumeOfBlock( int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	double volume(0.0);

	const std::vector< std::pair<int,double> >& blk2Elem = getBlockID2Elements(iblk);
	for(  std::vector< std::pair<int,double> >::const_iterator itr = blk2Elem.begin(); itr != blk2Elem.end(); ++itr ){
		const int elemID = itr->first;
		volume += ptrMeshData->calcVolume( elemID );
	}

	return volume;
}

// Calculate pre-degenerated roughning matrix
void ResistivityBlock::calcRougheningMatrix(){

	m_rougheningMatrix.setDegreeOfEquation( m_numResistivityBlockTotal );

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	switch ( ptrAnalysisControl->geTypeOfRoughningMatrix() )
	{
		case AnalysisControl::USE_ELEMENTS_SHARE_FACES:
			// Calculate roughening matrix with using elements share faces 
			calcRougheningMatrixUsingElementsShareFaces( 1.0 );
			break;
		case AnalysisControl::USER_DEFINED_ROUGHNING:
			// Calculate roughening matrix from user-defined roughning factor
			calcRougheningMatrixUserDefined( 1.0 );
			break;
		case AnalysisControl::USE_RESISTIVITY_BLOCKS_SHARE_FACES:
			// Calculate roughening matrix with using resistivity blocks share faces 
			calcRougheningMatrixUsingResistivityBlocksShareFaces( 1.0 );
			break;
		case AnalysisControl::USE_ELEMENTS_SHARE_FACES_AREA_VOL_RATIO:
			// Calculate roughening matrix with using elements share faces (weighting by area-volume ratio)
			calcRougheningMatrixUsingElementsShareFacesWeightingByAreaVolumeRatio( 1.0 );
			break;
		default:
			OutputFiles::m_logFile << "Error : Type of the method creating roughning matrix is wrong !!" << std::endl;
			exit(1);
			break;
	}

	// Output roughning matrix before converting to CRS format
	if( ptrAnalysisControl->getIsRougheningMatrixOutputted() && ptrAnalysisControl->getMyPE() == 0 ){
		// Only PE 0 output roughning factor
		m_rougheningMatrix.outputRougheningMatrix();
	}

	if( includeBottomResistivity() ){
		// Calculate roughning matrix from user-defined roughning factor
		addBottomResistivityContribution();
	}

	if( getFlagAddSmallValueToDiagonals() ){
		if( ptrAnalysisControl->getInversionMethod() == Inversion::GAUSS_NEWTON_DATA_SPECE && 
			ptrAnalysisControl->getTypeOfDataSpaceAlgorithm() == AnalysisControl::NEW_DATA_SPACE_ALGORITHM_USING_INV_RTR_MATRIX ){
			// In this case, small values are added to diagonals of [R]T[R] matrix
		}else{
			// Add small value to diagonals of roughning matrix
			addSmallValueToDiagonals();
		}
	}

	m_rougheningMatrix.convertToCRSFormat();

}

// Calculate roughning matrix degenerated for laplacian filter
void ResistivityBlock::calcRougheningMatrixDegeneratedForLaplacianFilter( DoubleSparseMatrix& rougheningMatrixDegenerated, const double factor ) const{

	const int nRow = m_rougheningMatrix.getNumRows();
	for( int iRow = 0; iRow < nRow; ++iRow ){
		if( isFixedResistivityValue(iRow) ){
			continue;
		}
		const int iRowDeg = getModelIDFromBlockID( iRow );
		rougheningMatrixDegenerated.addRightHandSideVector( iRowDeg, m_rougheningMatrix.getRightHandSideVector(iRow) * factor );
		const int nonZeroEnd = m_rougheningMatrix.getRowIndexCRS(iRow+1);
		for( int iNonZero = m_rougheningMatrix.getRowIndexCRS(iRow); iNonZero < nonZeroEnd; ++iNonZero ){
			const int iCol = m_rougheningMatrix.getColumnsCRS(iNonZero);
			const double val = m_rougheningMatrix.getValueCRS(iNonZero);
			if( isFixedResistivityValue(iCol) ){
				rougheningMatrixDegenerated.addRightHandSideVector( iRowDeg, - val * log10(getResistivityValuesPreFromBlockID(iCol)) * factor );
			}
			else{
				const int iColDeg = getModelIDFromBlockID( iCol );
				rougheningMatrixDegenerated.setStructureAndAddValueByTripletFormat( iRowDeg, iColDeg, val*factor ); 
			}
		}
	}

}

// Calculate roughning matrix degenerated for difference filter
void ResistivityBlock::calcRougheningMatrixDegeneratedForDifferenceFilter( const double factor,
	std::vector< std::pair<int,int> >& nonZeroCols, std::vector<double>& matValues, std::vector<double>& rhsValues ) const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const int pValue = ptrAnalysisControl->getDegreeOfLpOptimization();
	const double diffLog10RhoMin = ptrAnalysisControl->getLowerLimitOfDifflog10RhoForLpOptimization();
	const double diffLog10RhoMax = ptrAnalysisControl->getUpperLimitOfDifflog10RhoForLpOptimization();

	const int nRow = m_rougheningMatrix.getNumRows();
	for( int iRow = 0; iRow < nRow; ++iRow ){
		if( isFixedResistivityValue(iRow) ){
			continue;
		}
		const double log10rho = log10(getResistivityValuesFromBlockID(iRow));
		const int iRowDeg = getModelIDFromBlockID( iRow );
		const int nonZeroEnd = m_rougheningMatrix.getRowIndexCRS(iRow+1);
		double diffComps(0.0);
		for( int iNonZero = m_rougheningMatrix.getRowIndexCRS(iRow); iNonZero < nonZeroEnd; ++iNonZero ){
			const int iCol = m_rougheningMatrix.getColumnsCRS(iNonZero);
			const double val = m_rougheningMatrix.getValueCRS(iNonZero);
			diffComps += val;
			if( iRow == iCol ){
				continue;
			}
			const double log10rhoNeib = log10(getResistivityValuesFromBlockID(iCol));
			double diffLog10rho = fabs( log10rho - log10rhoNeib );
			if( diffLog10rho < diffLog10RhoMin ){
				diffLog10rho = diffLog10RhoMin;
			}
			if( diffLog10rho > diffLog10RhoMax ){
				diffLog10rho = diffLog10RhoMax;
			}
			const double weight = sqrt( 0.5 * pValue * pow( diffLog10rho, pValue - 2 ) );
			if( isFixedResistivityValue(iCol) ){
				nonZeroCols.push_back( std::make_pair(iRowDeg, -1) );
				matValues.push_back( -val * factor * weight );
				rhsValues.push_back( -val * log10rhoNeib * factor * weight );
			}else{
				const int iColDeg = getModelIDFromBlockID( iCol );
				nonZeroCols.push_back( std::make_pair(iRowDeg, iColDeg) );
				matValues.push_back( -val * factor * weight );
				rhsValues.push_back( 0.0 );
			}
		}
		// It is assumed that original right hand vector contains only bottom resistivity
		if( fabs(diffComps) > CommonParameters::EPS ){
			const double log10rhoNeib = log10(m_bottomResistivity);
			double diffLog10rho = fabs( log10rho - log10rhoNeib );
			if( diffLog10rho < diffLog10RhoMin ){
				diffLog10rho = diffLog10RhoMin;
			}
			if( diffLog10rho > diffLog10RhoMax ){
				diffLog10rho = diffLog10RhoMax;
			}
			const double weight = sqrt( 0.5 * pValue * pow( diffLog10rho, pValue - 2 ) );
			nonZeroCols.push_back( std::make_pair(iRowDeg, -1) );
			matValues.push_back( diffComps * factor * weight );
			rhsValues.push_back( m_rougheningMatrix.getRightHandSideVector(iRow) * factor * weight );
		}
	}

}

// Calculate array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 )
// from common logarithm
void ResistivityBlock::calctResistivityUpdatedFullFromLog10ResistivityIncres( const double* const log10resistivity ){

	const int typeBoundConstraints = getTypeBoundConstraints();

	int iMdl(0);
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		if( isFixedResistivityValue(iBlk) ){// If resistivity value is fixed
			m_resistivityValuesUpdatedFull[iBlk] = m_resistivityValuesPre[iBlk];
			continue;
		}
		assert( iMdl == getModelIDFromBlockID(iBlk) );
		if( typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD ){
			const double jacobian = calcDerivativeXWithRespectToLog10Resistivity(iBlk);// dx / dm
			const double xUpdatedFull = calcTransformedModelParameterFromResistivity(iBlk, m_resistivityValuesPre[iBlk]) + jacobian * log10resistivity[iMdl];
			m_resistivityValuesUpdatedFull[iBlk] = calcResistivityFromTransformedModelParameter(iBlk, xUpdatedFull);
		}
		else if( typeBoundConstraints == ResistivityBlock::SIMPLE_BOUND_CONSTRAINING ){
			const double log10ResistivityIncreUpdatedFull = log10( m_resistivityValuesPre[iBlk] ) + log10resistivity[iMdl];
			m_resistivityValuesUpdatedFull[iBlk] = pow( 10.0, log10ResistivityIncreUpdatedFull );
		}
		else{
			OutputFiles::m_logFile << "Error : Wrong type of bound constraining metohd !! : " << m_typeBoundConstraints << " ." << std::endl;
			exit(1);
		}
		++iMdl;
	}
	
	assert( iMdl == m_numResistivityBlockNotFixed );
}

// Copy derivatives of logarithm of resistivities with respect to transformed model parameter x
void ResistivityBlock::copyDerivativeLog10ResistivityWithRespectToX( double* derivs ) const{

	const int numBlockNotFixed = getNumResistivityBlockNotFixed();
	for( int iMdl = 0; iMdl < numBlockNotFixed; ++iMdl ){
		derivs[iMdl] = calcDerivativeLog10ResistivityWithRespectToX( getBlockIDFromModelID(iMdl) );
	}

}

// Change resistivity values
void ResistivityBlock::updateResistivityValues(){

	const double stepLengthDampingFactor = ( AnalysisControl::getInstance() )->getStepLengthDampingFactorCur();

	const int typeBoundConstraints = getTypeBoundConstraints();
	if( typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD ){
		for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
			if( isFixedResistivityValue(iBlk) ){// If resistivity value is fixed
				continue;
			}
			const double xUpdatedFull = calcTransformedModelParameterFromResistivity(iBlk, m_resistivityValuesUpdatedFull[iBlk]);
			const double xPre = calcTransformedModelParameterFromResistivity(iBlk, m_resistivityValuesPre[iBlk]);
			const double xNew = stepLengthDampingFactor * xUpdatedFull + ( 1.0 - stepLengthDampingFactor ) * xPre;
			m_resistivityValues[iBlk] = calcResistivityFromTransformedModelParameter(iBlk, xNew);
		}	
	}
	else if( typeBoundConstraints == ResistivityBlock::SIMPLE_BOUND_CONSTRAINING ){
		for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
			if( isFixedResistivityValue(iBlk) ){// If resistivity value is fixed
				continue;
			}
			const double newValueLog10 = stepLengthDampingFactor * log10( m_resistivityValuesUpdatedFull[iBlk] ) + ( 1.0 - stepLengthDampingFactor ) * log10( m_resistivityValuesPre[iBlk] );
			m_resistivityValues[iBlk] = pow( 10.0, newValueLog10 );
			if( m_resistivityValues[iBlk] < m_resistivityValuesMin[iBlk] ){
				OutputFiles::m_logFile << "Warning : Updated resistivity value ( " << m_resistivityValues[iBlk] << " [Ohm-m] ) of block " << iBlk << " is lower than the minimum value ( " << m_resistivityValuesMin[iBlk] <<  " [Ohm-m] ). Its resistivity is set to be the minimum value." << std::endl;
				m_resistivityValues[iBlk] = m_resistivityValuesMin[iBlk];
			}
			else if( m_resistivityValues[iBlk] > m_resistivityValuesMax[iBlk] ){
				OutputFiles::m_logFile << "Warning : Updated resistivity value ( " << m_resistivityValues[iBlk] << " [Ohm-m] ) of block " << iBlk << " is higher the maximum value ( " << m_resistivityValuesMax[iBlk] <<  " [Ohm-m] ). Its resistivity is set to be the maximum value." << std::endl;
				m_resistivityValues[iBlk] = m_resistivityValuesMax[iBlk];
			}
		}
	}
	else{
		OutputFiles::m_logFile << "Error : Wrong type of bound constraining metohd !! : " << m_typeBoundConstraints << " ." << std::endl;
		exit(1);
	}

}


// Output resistivity data to VTK file
void ResistivityBlock::outputResistivityDataToVTK() const{

	if( !OutputFiles::m_vtkFile.is_open() ){
		return;
	}

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	if( pAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK ) ){
		//MeshDataBrickElement* pMeshDataBrickElement = MeshDataBrickElement::getInstance();	
		const int nElem = ( ( AnalysisControl::getInstance() )->getPointerOfMeshData() )->getNumElemTotal();

		OutputFiles::m_vtkFile << "SCALARS blockIDs int" <<  std::endl;
		OutputFiles::m_vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){
			OutputFiles::m_vtkFile << m_elementID2blockID[iElem] << std::endl;
		}

	}

}

// Output resistivity data to binary file
void ResistivityBlock::outputResistivityDataToBinary() const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	
	if( ptrAnalysisControl->getMyPE() != 0 ){// Only PE 0 output sensitivity
		return;
	}

	std::ofstream fout;
	fout.open( "BlockIDs", std::ios::out | std::ios::binary | std::ios::trunc );

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Block IDs";
	strcpy( line, ossTitle.str().c_str() );
	fout.write( line, 80 );

	strcpy( line, "part" );
	fout.write( line, 80 );

	int ibuf(1);
	fout.write( (char*) &ibuf, sizeof( int ) );

	if( ptrAnalysisControl->getTypeOfMesh() ==  MeshData::TETRA ){
		strcpy( line, "tetra4" );
	}
	else{
		strcpy( line, "hexa8" );
	}
	fout.write( line, 80 );

	const MeshData* const ptrMeshData =ptrAnalysisControl->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();

	for( int iElem = 0 ; iElem < nElem; ++iElem ){
		float dbuf = static_cast<float>( m_elementID2blockID[iElem] );
		fout.write( (char*) &dbuf, sizeof( float ) );
	}

	fout.close();

}


// Get arrays of elements belonging to each resistivity block
//const std::vector<int>& ResistivityBlock::getBlockID2Elements( const int iBlk ) const{
const std::vector< std::pair<int,double> >&  ResistivityBlock::getBlockID2Elements( const int iBlk ) const{

	return m_blockID2Elements[iBlk];

}

#ifdef _ANISOTOROPY
//// Calculate anisotropy coefficient rotated
//void ResistivityBlock::calcAisotropyCoefficientRotated( const int iBlk, CommonParameters::Vector3D& coeffX, CommonParameters::Vector3D& coeffY, CommonParameters::Vector3D& coeffZ ) const{
//
//	const double angle = - m_anosotoropyStrileAngle[iBlk] * CommonParameters::deg2rad;
//
//	const CommonParameters::Vector3D tmpX = {
//		cos(angle) * m_anisotoropyCoeff[iBlk].X,
//		sin(angle) * m_anisotoropyCoeff[iBlk].Y,
//		0.0
//	};
//
//	const CommonParameters::Vector3D tmpY = {
//		-sin(angle) * m_anisotoropyCoeff[iBlk].X,
//		 cos(angle) * m_anisotoropyCoeff[iBlk].Y,
//		 0.0
//	};
//
//	const CommonParameters::Vector3D tmpZ = {
//		0.0,
//		0.0,
//		m_anisotoropyCoeff[iBlk].Z
//	};
//
//	coeffX.X =   cos(angle) * tmpX.X + sin(angle) * tmpX.Y;
//	coeffX.Y = - sin(angle) * tmpX.X + cos(angle) * tmpX.Y;
//	coeffX.Z = 0.0;
//
//	coeffY.X =   cos(angle) * tmpY.X + sin(angle) * tmpY.Y;
//	coeffY.Y = - sin(angle) * tmpY.X + cos(angle) * tmpY.Y;
//	coeffY.Z = 0.0;
//
//	coeffZ.X = 0.0;
//	coeffZ.Y = 0.0;
//	coeffZ.Z = tmpZ.Z;
//	
//}

// Calculate anisotropic conductivity tensor
//void ResistivityBlock::calcAisotropicConductivityTensor( const int blockID, CommonParameters::Vector3D& matX, CommonParameters::Vector3D& matY, CommonParameters::Vector3D& matZ ) const{
//
//	assert( ( AnalysisControl::getInstance() )->isAnisotropyConsidered() );
//
//	const std::map<int, int>::const_iterator itr = m_mapBlockIDWithAnisotropyToIndex.find(blockID);
//	if( itr == m_mapBlockIDWithAnisotropyToIndex.end() ){
//		// When the conducitivity of the specified block is isotropic
//		const double sigma = 1.0 / getResistivityValuesFromBlockID(blockID);
//		matX.X = sigma;
//		matX.Y = 0.0;
//		matX.Z = 0.0;
//		matY.X = 0.0;
//		matY.Y = sigma;
//		matY.Z = 0.0;
//		matZ.X = 0.0;
//		matZ.Y = 0.0;
//		matZ.Z = sigma;
//		return;
//	}
//
//	// Index of anisotropic block
//	const int iBlk = itr->second;
//
//	const double strikeAngle = - m_axialAnisotropyStrileAngle[iBlk] * CommonParameters::deg2rad;
//
//	// Rotation around z-axis
//	CommonParameters::Vector3D tmpX = {
//		cos(strikeAngle) / m_resistivityValuesAxialAnisotropy[iBlk].X,
//		sin(strikeAngle) / m_resistivityValuesAxialAnisotropy[iBlk].Y,
//		0.0
//	};
//
//	CommonParameters::Vector3D tmpY = {
//		-sin(strikeAngle) / m_resistivityValuesAxialAnisotropy[iBlk].X,
//		 cos(strikeAngle) / m_resistivityValuesAxialAnisotropy[iBlk].Y,
//		 0.0
//	};
//
//	CommonParameters::Vector3D tmpZ = {
//		0.0,
//		0.0,
//		1.0 / m_resistivityValuesAxialAnisotropy[iBlk].Z
//	};
//
//	matX.X =   cos(strikeAngle) * tmpX.X + sin(strikeAngle) * tmpX.Y;
//	matX.Y = - sin(strikeAngle) * tmpX.X + cos(strikeAngle) * tmpX.Y;
//	matX.Z = 0.0;
//
//	matY.X =   cos(strikeAngle) * tmpY.X + sin(strikeAngle) * tmpY.Y;
//	matY.Y = - sin(strikeAngle) * tmpY.X + cos(strikeAngle) * tmpY.Y;
//	matY.Z = 0.0;
//
//	matZ.X = 0.0;
//	matZ.Y = 0.0;
//	matZ.Z = tmpZ.Z;
//
//	// Rotation around x'-axis
//	const double dipAngle = - m_axialAnisotropyDipAngle[iBlk] * CommonParameters::deg2rad;;
//
//	tmpX.X = matX.X;
//	tmpX.Y = matX.Y;
//	tmpX.Z = matX.Z;
//
//	tmpY.X =   cos(dipAngle) * matY.X + sin(dipAngle) * matZ.X;
//	tmpY.Y =   cos(dipAngle) * matY.Y + sin(dipAngle) * matZ.Y;
//	tmpY.Z =   cos(dipAngle) * matY.Z + sin(dipAngle) * matZ.Z;
//
//	tmpZ.X = - sin(dipAngle) * matY.X + cos(dipAngle) * matZ.X;
//	tmpZ.Y = - sin(dipAngle) * matY.Y + cos(dipAngle) * matZ.Y;
//	tmpZ.Z = - sin(dipAngle) * matY.Z + cos(dipAngle) * matZ.Z;
//
//	matX.X = tmpX.X;
//	matX.Y =   cos(dipAngle) * tmpX.Y + sin(dipAngle) * tmpX.Z;
//	matX.Z = - sin(dipAngle) * tmpX.Y + cos(dipAngle) * tmpX.Z;
//
//	matY.X = tmpY.X;
//	matY.Y =   cos(dipAngle) * tmpY.Y + sin(dipAngle) * tmpY.Z;
//	matY.Z = - sin(dipAngle) * tmpY.Y + cos(dipAngle) * tmpY.Z;
//
//	matZ.X = tmpZ.X;
//	matZ.Y =   cos(dipAngle) * tmpZ.Y + sin(dipAngle) * tmpZ.Z;
//	matZ.Z = - sin(dipAngle) * tmpZ.Y + cos(dipAngle) * tmpZ.Z;
//
//}
void ResistivityBlock::calcAisotropicConductivityTensor( const int blockID, double conductivityTensor[3][3] ) const{

	assert( ( AnalysisControl::getInstance() )->isAnisotropyConsidered() );

	// Zero clear
	for( int row = 0; row < 3; ++row ){
		for( int col = 0; col < 3; ++col ){
			conductivityTensor[row][col] = 0.0;
		}
	}

	const std::map<int, int>::const_iterator itr = m_mapBlockIDWithAnisotropyToIndex.find(blockID);

	if( itr == m_mapBlockIDWithAnisotropyToIndex.end() ){
		// When the conducitivity of the specified block is isotropic
		const double sigma = 1.0 / getResistivityValuesFromBlockID(blockID);
		for( int row = 0; row < 3; ++row ){
			conductivityTensor[row][row] = sigma;
		}
		return;
	}

	// Index of anisotropic block
	const int iBlk = itr->second;

	conductivityTensor[0][0] = 1.0 / m_resistivityValuesAxialAnisotropy[iBlk].X;
	conductivityTensor[1][1] = 1.0 / m_resistivityValuesAxialAnisotropy[iBlk].Y;
	conductivityTensor[2][2] = 1.0 / m_resistivityValuesAxialAnisotropy[iBlk].Z;

	// Rotation around z-axis
	const double strike = m_axialAnisotropyStrileAngle[iBlk];
	double rotationTensor[3][3];
	rotationTensor[0][0] =  cos(strike);
	rotationTensor[0][1] =  sin(strike);
	rotationTensor[0][2] =  0.0;
	rotationTensor[1][0] = -sin(strike);
	rotationTensor[1][1] =  cos(strike);
	rotationTensor[1][2] =  0.0;
	rotationTensor[2][0] =  0.0;
	rotationTensor[2][1] =  0.0;
	rotationTensor[2][2] =  1.0;
	rotateTensor( conductivityTensor, rotationTensor );

	// Rotation around x'-axis
	const double dip = m_axialAnisotropyDipAngle[iBlk];
	rotationTensor[0][0] =  1.0;
	rotationTensor[0][1] =  0.0;
	rotationTensor[0][2] =  0.0;
	rotationTensor[1][0] =  0.0;
	rotationTensor[1][1] =  cos(dip);
	rotationTensor[1][2] =  sin(dip);
	rotationTensor[2][0] =  0.0;
	rotationTensor[2][1] = -sin(dip);
	rotationTensor[2][2] =  cos(dip);
	rotateTensor( conductivityTensor, rotationTensor );

	// Rotation around z'-axis
	const double slant = m_axialAnisotropySlantAngle[iBlk];
	rotationTensor[0][0] =  cos(slant);
	rotationTensor[0][1] =  sin(slant);
	rotationTensor[0][2] =  0.0;
	rotationTensor[1][0] = -sin(slant);
	rotationTensor[1][1] =  cos(slant);
	rotationTensor[1][2] =  0.0;
	rotationTensor[2][0] =  0.0;
	rotationTensor[2][1] =  0.0;
	rotationTensor[2][2] =  1.0;
	rotateTensor( conductivityTensor, rotationTensor );

}
#endif

// Output resistivity values to VTK file
void ResistivityBlock::outputResistivityValuesToVTK() const{

	if( !OutputFiles::m_vtkFile.is_open() ){
		return;
	}

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if( ptrAnalysisControl->doesOutputToVTK( AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK ) ){

		//const MeshDataBrickElement* const ptrMeshDataBrickElement = MeshDataBrickElement::getInstance();	
		//const int nElem = ptrMeshDataBrickElement->getNumElemTotal();
		const int nElem = ( ( AnalysisControl::getInstance() )->getPointerOfMeshData() )->getNumElemTotal();

		OutputFiles::m_vtkFile << "SCALARS Resistivity[Ohm-m] float" << std::endl;
		OutputFiles::m_vtkFile << "LOOKUP_TABLE default" <<  std::endl;
		for( int iElem = 0 ; iElem < nElem; ++iElem ){ 
			OutputFiles::m_vtkFile << m_resistivityValues[ m_elementID2blockID[ iElem ] ] << std::endl;
		}

	}

}

// Output resistivity values to binary file
void ResistivityBlock::outputResistivityValuesToBinary( const int iterNum ) const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	if( ptrAnalysisControl->getMyPE() != 0 ){// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Resistivity.iter" << iterNum;
	std::ofstream fout;
	fout.open( oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc );

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Resistivity[Ohm-m]";
	strcpy( line, ossTitle.str().c_str() );
	fout.write( line, 80 );

	strcpy( line, "part" );
	fout.write( line, 80 );

	int ibuf(1);
	fout.write( (char*) &ibuf, sizeof( int ) );

	if( ptrAnalysisControl->getTypeOfMesh() ==  MeshData::TETRA ){
		strcpy( line, "tetra4" );
	}
	else{
		strcpy( line, "hexa8" );
	}
	fout.write( line, 80 );

	const MeshData* const ptrMeshData =ptrAnalysisControl->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();

	for( int iElem = 0 ; iElem < nElem; ++iElem ){
		float dbuf = static_cast<float>( m_resistivityValues[ m_elementID2blockID[ iElem ] ] );
		fout.write( (char*) &dbuf, sizeof( float ) );
	}

	fout.close();

}

// Output data of resisitivity block model to file
void ResistivityBlock::outputResisitivityBlock( const int iterNum ) const{

	//const MeshDataBrickElement* const ptrMeshDataBrickElement = MeshDataBrickElement::getInstance();

	std::ostringstream fileName;
	//fileName << "resistivity_model" << iterNum << ".dat";
	//fileName << "resistivity_model_iter" << iterNum << ".dat";
	fileName << "resistivity_block_iter" << iterNum << ".dat";

	FILE *fp;
	if( (fp = fopen( fileName.str().c_str(), "w")) == NULL ) {
		OutputFiles::m_logFile  << "File open error !! : " << fileName.str() << std::endl;
		exit(1);
	}

	//const int numElemTotal = ptrMeshDataBrickElement->getNumElemTotal();
	const int numElemTotal = ( ( AnalysisControl::getInstance() )->getPointerOfMeshData() )->getNumElemTotal();
	fprintf(fp, "%10d%10d\n",numElemTotal, m_numResistivityBlockTotal );

	for( int iElem = 0; iElem < numElemTotal; ++iElem ){
		fprintf(fp, "%10d%10d\n", iElem, m_elementID2blockID[iElem] );
	}

	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", iBlk, "     ",
			m_resistivityValues[iBlk], m_resistivityValuesMin[iBlk], m_resistivityValuesMax[iBlk], m_weightingConstants[iBlk],
			getTypeOfResistivityBlock(m_fixResistivityValues[iBlk], m_isolated[iBlk]) );
	}

	fclose(fp);

}

// Copy resistivity values not fixed to vector
void ResistivityBlock::copyResistivityValuesNotFixedToVectorLog10( double* vector ) const{

	for( int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl ){
		vector[iMdl] = log10( m_resistivityValues[ m_modelID2blockID[ iMdl ] ] );
	}

}

// Copy common logarithm of previous free resistivity values to vector
void ResistivityBlock::copyResistivityValuesNotFixedPreToVectorLog10( double* vector ) const{

	for( int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl ){
		vector[iMdl] = log10( m_resistivityValuesPre[ m_modelID2blockID[ iMdl ] ] );
	}

}

// Copy common logarithm of current free resistivity values to previous ones
void ResistivityBlock::copyResistivityValuesNotFixedCurToPre() const{

	// Convert current resistivity vector to previous one
	memcpy( m_resistivityValuesPre, m_resistivityValues, sizeof(double)*(m_numResistivityBlockTotal) );

}

// Calculate model roughness for laplacian filter
double ResistivityBlock::calcModelRoughnessForLaplacianFilter() const{

	double* modelVec = new double[m_numResistivityBlockTotal];

	copyResistivityValuesToVectorLog10( modelVec );

	const double roughness = m_rougheningMatrix.calcModelRoughness( modelVec );

	delete[] modelVec;

	return roughness;

}

// Calculate model roughness for difference filter
double ResistivityBlock::calcModelRoughnessForDifferenceFilter() const{

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const int pValue = ptrAnalysisControl->getDegreeOfLpOptimization();
	const double diffLog10RhoMin = ptrAnalysisControl->getLowerLimitOfDifflog10RhoForLpOptimization();
	const double diffLog10RhoMax = ptrAnalysisControl->getUpperLimitOfDifflog10RhoForLpOptimization();

	double roughness(0.0);
	const int nRow = m_rougheningMatrix.getNumRows();
	for( int iRow = 0; iRow < nRow; ++iRow ){
		const int nonZeroEnd = m_rougheningMatrix.getRowIndexCRS(iRow+1);
		const double log10rho = log10(getResistivityValuesFromBlockID(iRow));
		double diffComps(0.0);
		for( int iNonZero = m_rougheningMatrix.getRowIndexCRS(iRow); iNonZero < nonZeroEnd; ++iNonZero ){
			const int iCol = m_rougheningMatrix.getColumnsCRS(iNonZero);
			const double val = m_rougheningMatrix.getValueCRS(iNonZero);
			diffComps += val;
			if( iRow == iCol ){
				continue;
			}
			const double log10rhoNeib = log10(getResistivityValuesFromBlockID(iCol));
			double diffLog10rho = fabs( log10rho - log10rhoNeib );
			if( diffLog10rho < diffLog10RhoMin ){
				diffLog10rho = diffLog10RhoMin;
			}
			if( diffLog10rho > diffLog10RhoMax ){
				diffLog10rho = diffLog10RhoMax;
			}
			const double weight = 0.5 * pValue * pow( diffLog10rho, pValue - 2 );
			roughness += pow( ( log10rho - log10rhoNeib ) * val, 2 ) * weight;
		}
		if( fabs(diffComps) > CommonParameters::EPS ){
			// It is assumed that original right hand vector contains only bottom resistivity
			const double log10rhoNeib = log10(m_bottomResistivity);
			double diffLog10rho = fabs( log10rho - log10rhoNeib );
			if( diffLog10rho < diffLog10RhoMin ){
				diffLog10rho = diffLog10RhoMin;
			}
			if( diffLog10rho > diffLog10RhoMax ){
				diffLog10rho = diffLog10RhoMax;
			}
			const double weight = 0.5 * pValue * pow( diffLog10rho, pValue - 2 );
			roughness += pow( ( log10rho - log10rhoNeib ) * diffComps, 2 ) * weight;
		}
	}
	return roughness;

}

// Calculate model roughness at the bottom
double ResistivityBlock::calcModelRoughnessAtBottom() const{

	if( !includeBottomResistivity() ){
		return 0.0;
	}

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();
	const int numElemBot = ptrMeshData->getNumElemOnBoundaryPlanes(MeshData::XYPlus);

	std::vector< std::pair<int,double> > diffs;
	diffs.reserve(numElemBot);

	for( int i = 0; i < numElemBot; ++i ){
		const int iElem = ptrMeshData->getElemBoundaryPlanes(MeshData::XYPlus, i);
		const int iBlk = getBlockIDFromElemID( iElem );
		if( isolated(iBlk) ){
			continue;
		}
		const double diff = m_roughningFactorAtBottom * ( log10(m_bottomResistivity) - log10(m_resistivityValues[iBlk]) );
		diffs.push_back( std::make_pair(iBlk, diff) );
	}

	sort( diffs.begin(), diffs.end() );

	double roughness(0.0);
	double diff(0.0);
	int iBlkPre(-1);
	for( std::vector< std::pair<int,double> >::const_iterator itr = diffs.begin(); itr != diffs.end(); ++itr ){
		if( itr->first != iBlkPre ){
			roughness += diff * diff;
			diff = itr->second;
			iBlkPre = itr->first;
		}
		else{
			diff += itr->second;
		}
	}
	roughness += diff * diff;

	return roughness;

}

// Copy common logarithm of resistivity values
void ResistivityBlock::copyResistivityValuesToVectorLog10( double* vector ) const{

	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		vector[iBlk] = log10( m_resistivityValues[iBlk] );
	}

}

// Calculate weighting factor
double ResistivityBlock::calWeightingFactor( const double alphaX, const double alphaY, const double alphaZ, const int iElem1, const int iElem2 ) const{
	const CommonParameters::locationXYZ diff = ( (AnalysisControl::getInstance())->getPointerOfMeshData() )->calDiffOfCenters(iElem1, iElem2);
	const double numer = sqrt( pow(alphaX*diff.X, 2) + pow(alphaY*diff.Y, 2) + pow(alphaZ*diff.Z, 2) ) * 1.0e-3;
	const double denom = sqrt( pow(diff.X, 2) + pow(diff.Y, 2) + pow(diff.Z, 2) ) * 1.0e-3;
	return numer / pow( denom, m_inverseDistanceWeightingFactor );
}

// Get type of resistivity block
int ResistivityBlock::getTypeOfResistivityBlock( const bool fixed, const bool isolated ) const{
	if( fixed ){
		if( isolated ){
			return ResistivityBlock::FIXED_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlock::FIXED_AND_CONSTRAINED;
		}
	}
	else{// Free
		if( isolated ){
			return ResistivityBlock::FREE_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlock::FREE_AND_CONSTRAINED;
		}
	}
}

// Calculate roughning matrix with using elements share faces 
void ResistivityBlock::calcRougheningMatrixUsingElementsShareFaces( const double factor ){

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	const int nElem = ptrMeshData->getNumElemTotal();

	const double alphaX = ( AnalysisControl::getInstance() )->getAlphaWeight(0);
	const double alphaY = ( AnalysisControl::getInstance() )->getAlphaWeight(1);
	const double alphaZ = ( AnalysisControl::getInstance() )->getAlphaWeight(2);

	if( ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA ){
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement(); 
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			for( int iFace = 0; iFace < 6; ++iFace ){
				if( ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace) ){
					continue;
				}
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				for( int iNeigh = 0; iNeigh < numNeibs; ++iNeigh ){
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if( iElemNeib < 0 ){
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
					if( isolated(iBlkNeib) ){
						continue;
					}
					if( iBlk != iBlkNeib ){
						const double weightedFactor = factor * calWeightingFactor(alphaX, alphaY, alphaZ, iElem, iElemNeib) / static_cast<double>(numNeibs);
						m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk,     weightedFactor ); 
						m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlkNeib, -weightedFactor );
					}
				}
			}
		}
	}else{
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			const int numNeighborElement = ptrMeshData->getNumNeighborElement();
			for( int iNeigh = 0; iNeigh < numNeighborElement; ++iNeigh ){
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement( iElem, iNeigh );
				if( iElemNeib < 0 ){
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
				if( isolated(iBlkNeib) ){
					continue;
				}
				if( iBlk != iBlkNeib ){
					const double weightedFactor = factor * calWeightingFactor(alphaX, alphaY, alphaZ, iElem, iElemNeib);
					m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk,     weightedFactor ); 
					m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlkNeib, -weightedFactor );
				}
			}
		}
	}

}

// Calculate roughening matrix with using elements share faces using area-volume ratio as weights
void ResistivityBlock::calcRougheningMatrixUsingElementsShareFacesWeightingByAreaVolumeRatio( const double factor ){

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	const int nElem = ptrMeshData->getNumElemTotal();

	// Calculate volume of each parameter cell (resistivity block)
	const int numResistivityBlocks = getNumResistivityBlockTotal();
	double* volumes = new double[numResistivityBlocks];
	for( int iBlk = 0; iBlk < numResistivityBlocks; ++iBlk ){
		// Zero clear
		volumes[iBlk] = 0.0;
	}
	for( int iElem = 0; iElem < nElem; ++iElem ){
		const int iBlk = getBlockIDFromElemID( iElem );
		volumes[iBlk] += ptrMeshData->calcVolume(iElem); 
	}
#ifdef _DEBUG_WRITE
	for( int iBlk = 0; iBlk < numResistivityBlocks; ++iBlk ){
		std::cout << iBlk << " " << volumes[iBlk] << std::endl;
	}
#endif
	const double inverseDistanceWeightingFactor = getInverseDistanceWeightingFactor();
	if( ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA ){
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement(); 
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			for( int iFace = 0; iFace < 6; ++iFace ){
				if( ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace) ){
					continue;
				}
				const double area = ptrMeshData->calcAreaOfFace(iElem, iFace);
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				const double volume = volumes[iBlk];
				const double areaVolRatio = pow( area / volume / static_cast<double>(numNeibs), inverseDistanceWeightingFactor );
				for( int iNeigh = 0; iNeigh < numNeibs; ++iNeigh ){
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if( iElemNeib < 0 ){
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
					if( isolated(iBlkNeib) ){
						continue;
					}
					if( iBlk != iBlkNeib ){
						const double weightedFactor = factor * areaVolRatio;
						m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk,     weightedFactor ); 
						m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlkNeib, -weightedFactor );
					}
				}
			}
		}
	}else{
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			const int numNeighborElement = ptrMeshData->getNumNeighborElement();
			for( int iFace = 0; iFace < numNeighborElement; ++iFace ){
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement( iElem, iFace );
				if( iElemNeib < 0 ){
					continue;
				}
				const double area = ptrMeshData->calcAreaOfFace(iElem, iFace);
				const double volume = volumes[iBlk];
				const double areaVolRatio = pow( area / volume, inverseDistanceWeightingFactor );
				const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
				if( isolated(iBlkNeib) ){
					continue;
				}
				if( iBlk != iBlkNeib ){
					const double weightedFactor = factor * areaVolRatio;
					m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk,     weightedFactor ); 
					m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlkNeib, -weightedFactor );
				}
			}
		}
	}

	delete [] volumes;

}

// Calculate roughning matrix with using resistivty blocks share faces
void ResistivityBlock::calcRougheningMatrixUsingResistivityBlocksShareFaces( const double factor ){

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	const int nElem = ptrMeshData->getNumElemTotal();

	const double alphaX = ( AnalysisControl::getInstance() )->getAlphaWeight(0);
	const double alphaY = ( AnalysisControl::getInstance() )->getAlphaWeight(1);
	const double alphaZ = ( AnalysisControl::getInstance() )->getAlphaWeight(2);

	std::map< std::pair<int, int>, int> blockPairToCounter;
	std::map< std::pair<int, int>, double> blockPairToWeight;

	if( ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA ){
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement(); 
		for( int iElem = 0; iElem < nElem; ++iElem ){
			int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			for( int iFace = 0; iFace < 6; ++iFace ){
				if( ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace) ){
					continue;
				}
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				for( int iNeigh = 0; iNeigh < numNeibs; ++iNeigh ){
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if( iElemNeib < 0 ){
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
					if( isolated(iBlkNeib) ){
						continue;
					}
					if( iBlk != iBlkNeib ){
						const double weightedFactor = factor * calWeightingFactor(alphaX, alphaY, alphaZ, iElem, iElemNeib);
						const std::pair<int, int> blockPair = std::make_pair(iBlk, iBlkNeib);
						// Counter
						std::map< std::pair<int, int>, int>::iterator itrCounter = blockPairToCounter.find(blockPair);
						if( itrCounter != blockPairToCounter.end() ){
							// Found
							itrCounter->second += 1;
						}else{
							// Not found
							blockPairToCounter.insert( std::make_pair( blockPair, 1 ) );
						}
						// Weight
						std::map< std::pair<int, int>, double>::iterator itrWeight = blockPairToWeight.find(blockPair);
						if( itrWeight != blockPairToWeight.end() ){
							// Found
							itrWeight->second += weightedFactor;
						}else{
							// Not found
							blockPairToWeight.insert( std::make_pair( blockPair, weightedFactor ) );
						}
						assert( (itrCounter != blockPairToCounter.end()) == (itrWeight != blockPairToWeight.end()) );
					}
				}
			}
		}
	}else{
		for( int iElem = 0; iElem < nElem; ++iElem ){
			int iBlk = getBlockIDFromElemID( iElem );
			if( isolated(iBlk) ){
				continue;
			}
			const int numNeighborElement = ptrMeshData->getNumNeighborElement();
			for( int iNeigh = 0; iNeigh < numNeighborElement; ++iNeigh ){
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement( iElem, iNeigh );
				if( iElemNeib < 0 ){
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID( iElemNeib );
				if( isolated(iBlkNeib) ){
					continue;
				}
				if( iBlk != iBlkNeib ){
					const double weightedFactor = factor * calWeightingFactor(alphaX, alphaY, alphaZ, iElem, iElemNeib);
					const std::pair<int, int> blockPair = std::make_pair(iBlk, iBlkNeib);
					// Counter
					std::map< std::pair<int, int>, int>::iterator itrCounter = blockPairToCounter.find(blockPair);
					if( itrCounter != blockPairToCounter.end() ){
						// Found
						itrCounter->second += 1;
					}else{
						// Not found
						blockPairToCounter.insert( std::make_pair( blockPair, 1 ) );
					}
					// Weight
					std::map< std::pair<int, int>, double>::iterator itrWeight = blockPairToWeight.find(blockPair);
					if( itrWeight != blockPairToWeight.end() ){
						// Found
						itrWeight->second += weightedFactor;
					}else{
						// Not found
						blockPairToWeight.insert( std::make_pair( blockPair, weightedFactor ) );
					}
					assert( (itrCounter != blockPairToCounter.end()) == (itrWeight != blockPairToWeight.end()) );
				}
			}
		}
	}

	for( std::map< std::pair<int, int>, int>::const_iterator itrCounter = blockPairToCounter.begin();
		itrCounter != blockPairToCounter.end(); ++itrCounter ){
		const std::pair<int, int> blockPair = itrCounter->first;
		std::map< std::pair<int, int>, double>::const_iterator itrWeight = blockPairToWeight.find(blockPair);
		if( itrWeight == blockPairToWeight.end() ){
			// Not found
			OutputFiles::m_logFile << "Error : Weight correspondig to block pair (" << blockPair.first << ", " << blockPair.second << ") is not found !!" << std::endl;
			exit(1);
		}
		const int iBlk = blockPair.first;
		const int iBlkNeib = blockPair.second;
		const int counter = itrCounter->second;
		const double weight = itrWeight->second;
		const double weightedFactorAveraged = itrWeight->second / static_cast<double>(itrCounter->second);
		m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk,     weightedFactorAveraged );
		m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlkNeib, -weightedFactorAveraged );
	}

}

// Calculate roughning matrix from user-defined roughning factor
void ResistivityBlock::calcRougheningMatrixUserDefined( const double factor ){

	// Read user-defined roughening matrix
	const std::string fileName = "roughening_matrix.dat";
	std::ifstream ifs( fileName.c_str(), std::ios::in );

	if( ifs.fail() ){
		OutputFiles::m_logFile << "File open error : " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	OutputFiles::m_logFile << "# Read user-defined roughening matrix from " << fileName.c_str() << "." << std::endl;

	int ibuf(0);
	ifs >> ibuf;
	const int numBlock(ibuf);
	if( numBlock <= 0 ){
		OutputFiles::m_logFile << "Error : Total number of resistivity blocks must be positive !! : " << numBlock << std::endl;
		exit(1);
	}

	for( int iBlock = 0 ; iBlock < numBlock; ++iBlock ){
		ifs >> ibuf;
		if( iBlock != ibuf ){
			OutputFiles::m_logFile << "Error : Resistivity block numbers must be numbered consecutively from zero !!" << std::endl;
			exit(1);
		}

		ifs >> ibuf;
		const int numNonzeros(ibuf);
		std::vector< std::pair<int, double> > blockIDAndFactor;
		blockIDAndFactor.resize(numNonzeros);
		for( int innz = 0 ; innz < numNonzeros; ++innz ){
			ifs >> ibuf;
			blockIDAndFactor[innz].first = ibuf;
		}
		for( int innz = 0 ; innz < numNonzeros; ++innz ){
			double dbuf(0.0);
			ifs >> dbuf;
			blockIDAndFactor[innz].second = dbuf;
		}
		for( int innz = 0 ; innz < numNonzeros; ++innz ){
			m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlock, blockIDAndFactor[innz].first, blockIDAndFactor[innz].second );
		}
	}

	ifs.close();

}

// Add contribution of bottom resistivity to roughning matrix
void ResistivityBlock::addBottomResistivityContribution(){

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	// Calculate volume of each parameter cell (resistivity block)
	const int nElem = ptrMeshData->getNumElemTotal();
	const int numResistivityBlocks = getNumResistivityBlockTotal();
	double* volumes(NULL);
	const int typeOfRoughningMatrix = ptrAnalysisControl->geTypeOfRoughningMatrix();
	if( typeOfRoughningMatrix == AnalysisControl::USE_ELEMENTS_SHARE_FACES_AREA_VOL_RATIO ){
		volumes = new double[numResistivityBlocks];
		for( int iBlk = 0; iBlk < numResistivityBlocks; ++iBlk ){
			// Zero clear
			volumes[iBlk] = 0.0;
		}
		for( int iElem = 0; iElem < nElem; ++iElem ){
			const int iBlk = getBlockIDFromElemID( iElem );
			volumes[iBlk] += ptrMeshData->calcVolume(iElem); 
		}
	}

	const double inverseDistanceWeightingFactor = getInverseDistanceWeightingFactor();
	const int numElemBot = ptrMeshData->getNumElemOnBoundaryPlanes(MeshData::XYPlus);
	std::map< int, int> blockoCounter;
	std::map< int, double> blockToWeight;
	for( int iElemBot = 0; iElemBot < numElemBot; ++iElemBot ){
		const int iElem = ptrMeshData->getElemBoundaryPlanes(MeshData::XYPlus, iElemBot);
		const int iBlk = getBlockIDFromElemID( iElem );
		if( isolated(iBlk) ){
			continue;
		}
		double factor = m_roughningFactorAtBottom;
		if( typeOfRoughningMatrix == AnalysisControl::USE_ELEMENTS_SHARE_FACES_AREA_VOL_RATIO ){
			const double areaVolRatio = pow(ptrMeshData->calcAreaOfFaceAtBottomOfMesh(iElemBot) / volumes[iBlk], inverseDistanceWeightingFactor);
			factor *= areaVolRatio;
		}
		if( typeOfRoughningMatrix == AnalysisControl::USE_RESISTIVITY_BLOCKS_SHARE_FACES ){
			// Counter
			std::map<int, int>::iterator itrCounter = blockoCounter.find(iBlk);
			if( itrCounter != blockoCounter.end() ){
				// Found
				itrCounter->second += 1;
			}else{
				// Not found
				blockoCounter.insert( std::make_pair( iBlk, 1 ) );
			}
			// Weight
			std::map< int, double>::iterator itrWeight = blockToWeight.find(iBlk);
			if( itrWeight != blockToWeight.end() ){
				// Found
				itrWeight->second += factor;
			}else{
				// Not found
				blockToWeight.insert( std::make_pair( iBlk, factor ) );
			}
			assert( (itrCounter != blockoCounter.end()) == (itrWeight != blockToWeight.end()) );
		}else{
			m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk, factor );
			m_rougheningMatrix.addRightHandSideVector( iBlk, factor * log10(m_bottomResistivity) );
		}
	}
	if( typeOfRoughningMatrix == AnalysisControl::USE_RESISTIVITY_BLOCKS_SHARE_FACES ){
		for( std::map<int, int>::const_iterator itrCounter = blockoCounter.begin();
			itrCounter != blockoCounter.end(); ++itrCounter ){
			const int iBlk = itrCounter->first;
			std::map<int, double>::const_iterator itrWeight = blockToWeight.find(iBlk);
			if( itrWeight == blockToWeight.end() ){
				// Not found
				OutputFiles::m_logFile << "Error : Weight correspondig to block (" << iBlk << ") is not found !!" << std::endl;
				exit(1);
			}
			const int counter = itrCounter->second;
			const double weight = itrWeight->second;
			const double weightedFactorAveraged = itrWeight->second / static_cast<double>(itrCounter->second);
			m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk, weightedFactorAveraged );
			m_rougheningMatrix.addRightHandSideVector( iBlk, weightedFactorAveraged * log10(m_bottomResistivity) );
		}
	}

	if( volumes != NULL ){
		delete [] volumes;
	}

}

// Add small value to diagonals
void ResistivityBlock::addSmallValueToDiagonals(){

	const double value = getSmallValueAddedToDiagonals();
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk, value );
	}

}

// Calculate transformed model parameter x from resistivity
double ResistivityBlock::calcTransformedModelParameterFromResistivity( const int iBlk, const double resistivity ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	return calcTransformedModelParameterFromLog10Resistivity(iBlk, log10(resistivity));
}

// Calculate transformed model parameter x from common logarithm of resistivity
double ResistivityBlock::calcTransformedModelParameterFromLog10Resistivity( const int iBlk, const double log10Resistivity ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	const double m = log10Resistivity;
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];

	const double eps = 1.0e-5;
	const double val1 = b - m < eps ? eps : b - m;
	const double val2 = m - a < eps ? eps : m - a;
	return log( val2 / val1 ) / n;
}

// Calculate resistivity from transformed model parameter x
double ResistivityBlock::calcResistivityFromTransformedModelParameter( const int iBlk, const double x ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	return pow( 10.0, calcLog10ResistivityFromTransformedModelParameter(iBlk, x) );
}

// Calculate common logarithm of resistivity from transformed model parameter x
double ResistivityBlock::calcLog10ResistivityFromTransformedModelParameter( const int iBlk, const double x ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];
	return 0.5 * ( b - a ) * tanh( 0.5 * n * x ) + 0.5 * ( b + a );
}

// Calculate derivative of logarithm of resistivity with respect to transformed model parameter x
double ResistivityBlock::calcDerivativeLog10ResistivityWithRespectToX( const int iBlk ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	const double m = log10(m_resistivityValuesPre[iBlk]);
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];

	assert( m >= a );
	assert( m <= b );
	assert( b > a );
	assert( !m_fixResistivityValues[iBlk] );

	return n * (b - m) * (m - a) / (b - a);
}

// Calculate derivative of transformed model parameter x with respect to logarithm of resistivity
double ResistivityBlock::calcDerivativeXWithRespectToLog10Resistivity( const int iBlk ) const{
	assert( m_typeBoundConstraints == ResistivityBlock::TRANSFORMING_METHOD );
	const double m = log10(m_resistivityValuesPre[iBlk]);
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];

	assert( !m_fixResistivityValues[iBlk] );

	const double val1 = b - m < m_minDistanceToBounds ? m_minDistanceToBounds : b - m;
	const double val2 = m - a < m_minDistanceToBounds ? m_minDistanceToBounds : m - a;

	return (b - a) / ( n * val1 * val2 );
}

//// Add vector resulting from the product of sensitivity matrix and model vector
//void ResistivityBlock::addProductOfSensitivityMatrixAndModelVector( const double senMat, const int numData, double* vecOut ) const{
//
//	// Matrix-vector product
//	for( int iRow = 0; iRow < numData; ++iRow ){
//
//		const int offset = ( iRow + 1 ) * iRow / 2;
//		double work(0.0);
//
//		// Lower triangle part of the transposed matrix
//		for( int iCol = 0; iCol < iRow; ++iCol ){
//			work += matrixToBeInverted[ offset + iCol ] * dataVectorThisPE[iCol];
//		}
//
//		// Upper triangle part of the transposed matrix
//		for( int iCol = iRow; iCol < numData; ++iCol ){
//			work += matrixToBeInverted[ iRow + ( iCol + 1 ) * iCol / 2 ] * dataVectorThisPE[iCol];
//		}
//
//		rhsVectorThisPE[iRow] = work;
//	}
//
//}

//// Get resistivity model
//const double* const ResistivityBlock::getResistivityValues() const{
//
//	return m_resistivityValues;
//
//}
