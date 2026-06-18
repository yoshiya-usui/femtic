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
#include "ResistivityBlockIsotropic.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"

#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <algorithm>

// Constructer
ResistivityBlockIsotropic::ResistivityBlockIsotropic():
	m_numResistivityBlockNotFixed(0),
	m_resistivityValues(NULL),
	m_resistivityValuesPre(NULL),
	m_resistivityValuesUpdatedFull(NULL),
	m_weightingConstants(NULL),
	m_fixResistivityValues(NULL),
	m_isolated(NULL),
	m_rougheningMatrix()
{}

// Destructer
ResistivityBlockIsotropic::~ResistivityBlockIsotropic(){

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

	m_rougheningMatrix.releaseMemory();

}

// Read isotropic reslstivity values from input file
void ResistivityBlockIsotropic::inputResistivityValues(const int nElem, std::ifstream& inFile){

	if (m_blockID2modelID != NULL) {
		delete[] m_blockID2modelID;
		m_blockID2modelID = NULL;
	}
	m_blockID2modelID = new int[m_numResistivityBlockTotal];

	if (m_resistivityValues != NULL) {
		delete[] m_resistivityValues;
		m_resistivityValues = NULL;
	}
	m_resistivityValues = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesPre != NULL) {
		delete[] m_resistivityValuesPre;
		m_resistivityValuesPre = NULL;
	}
	m_resistivityValuesPre = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesUpdatedFull != NULL) {
		delete[] m_resistivityValuesUpdatedFull;
		m_resistivityValuesUpdatedFull = NULL;
	}
	m_resistivityValuesUpdatedFull = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesMin != NULL) {
		delete[] m_resistivityValuesMin;
		m_resistivityValuesMin = NULL;
	}
	m_resistivityValuesMin = new double[m_numResistivityBlockTotal];

	if (m_resistivityValuesMax != NULL) {
		delete[] m_resistivityValuesMax;
		m_resistivityValuesMax = NULL;
	}
	m_resistivityValuesMax = new double[m_numResistivityBlockTotal];

	if (m_weightingConstants != NULL) {
		delete[] m_weightingConstants;
		m_weightingConstants = NULL;
	}
	m_weightingConstants = new double[m_numResistivityBlockTotal];

	if (m_fixResistivityValues != NULL) {
		delete[] m_fixResistivityValues;
		m_fixResistivityValues = NULL;
	}
	m_fixResistivityValues = new bool[m_numResistivityBlockTotal];

	if (m_isolated != NULL) {
		delete[] m_isolated;
		m_isolated = NULL;
	}
	m_isolated = new bool[m_numResistivityBlockTotal];

	for (int i = 0; i < m_numResistivityBlockTotal; ++i) {
		m_resistivityValues[i] = 0.0;
		m_resistivityValuesPre[i] = 0.0;
		m_resistivityValuesUpdatedFull[i] = 0.0;
		m_resistivityValuesMin[i] = 0.0;
		m_resistivityValuesMax[i] = 0.0;
		m_weightingConstants[i] = 1.0;
		m_fixResistivityValues[i] = false;
		m_isolated[i] = false;
	}

	m_numResistivityBlockNotFixed = 0;
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		int idum(0);
		int itype(0);
		inFile >> idum >> m_resistivityValues[iBlk] >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk] >> m_weightingConstants[iBlk] >> itype;
		if (idum != iBlk) {
			OutputFiles::m_logFile << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		if (m_resistivityValues[iBlk] <= 0.0) {
			OutputFiles::m_logFile << "Error : Resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValues[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMin[iBlk] <= 0.0) {
			OutputFiles::m_logFile << "Error : Minimum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMin[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMax[iBlk] <= 0.0) {
			OutputFiles::m_logFile << "Error : Maximum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMax[iBlk] << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMax[iBlk] < m_resistivityValues[iBlk]) {
			OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
		if (m_resistivityValuesMin[iBlk] > m_resistivityValues[iBlk]) {
			OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << m_resistivityValues[iBlk] << " )." << std::endl;
			exit(1);
		}
		if (m_weightingConstants[iBlk] <= 0.0) {
			OutputFiles::m_logFile << "Error : Weighting constant of block " << iBlk << " is less than or equal to zero !! : " << m_weightingConstants[iBlk] << std::endl;
			exit(1);
		}
		if ((AnalysisControl::getInstance())->getInversionMethod() == Inversion::GAUSS_NEWTON_DATA_SPECE && itype == FREE_AND_ISOLATED) {
			OutputFiles::m_logFile << "Error : Resistivity of isolated block must be fixed when data space inverson method is selected !!" << std::endl;
			exit(1);
		}
		switch (itype) {
		case ResistivityBlockIsotropic::FREE_AND_CONSTRAINED:// Go through
		case ResistivityBlockIsotropic::FREE_AND_ISOLATED:
			m_blockID2modelID[iBlk] = m_numResistivityBlockNotFixed++;
			break;
		case ResistivityBlockIsotropic::FIXED_AND_ISOLATED:// Go through
		case ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED:
			m_fixResistivityValues[iBlk] = true;
			m_blockID2modelID[iBlk] = -1;
			break;
		default:
			OutputFiles::m_logFile << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
			exit(1);
			break;
		}
		switch (itype) {
		case ResistivityBlockIsotropic::FREE_AND_CONSTRAINED:// Go through
		case ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED:
			m_isolated[iBlk] = false;
			break;
		case ResistivityBlockIsotropic::FIXED_AND_ISOLATED:// Go through
		case ResistivityBlockIsotropic::FREE_AND_ISOLATED:
			m_isolated[iBlk] = true;
			break;
		default:
			OutputFiles::m_logFile << "Error : Type of resistivity block is unknown !! : " << itype << std::endl;
			exit(1);
			break;
		}
#ifdef _DEBUG_WRITE
		std::cout << std::setw(5) << iBlk << std::setw(15) << m_resistivityValues[iBlk] << std::setw(15) << m_resistivityValuesMin[iBlk] << std::setw(15) << m_resistivityValuesMax[iBlk]
			<< std::setw(15) << m_weightingConstants[iBlk] << std::setw(5) << m_fixResistivityValues[iBlk] << std::setw(5) << m_blockID2modelID[iBlk] << std::endl; // For debug
#endif
	}

	if (!m_fixResistivityValues[0]) {
		OutputFiles::m_logFile << "Error : Resistivity block 0 must be the air. And, its resistivity must be fixed." << std::endl;
		exit(1);
	}

	inFile.close();

	memcpy(m_resistivityValuesPre, m_resistivityValues, sizeof(double) * (m_numResistivityBlockTotal));

	if (m_modelID2blockID != NULL) {
		delete[] m_modelID2blockID;
		m_modelID2blockID = NULL;
	}
	m_modelID2blockID = new int[m_numResistivityBlockNotFixed];

	int icount(0);
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {

		if (!m_fixResistivityValues[iBlk]) {
			m_modelID2blockID[icount] = iBlk;
			++icount;
		}

	}

	if (icount != m_numResistivityBlockNotFixed) {
		OutputFiles::m_logFile << "Error : icount is not equal to m_numResistivityBlockNotFixed. icount = " << icount << " m_numResistivityBlockNotFixed = " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}
	if (m_numResistivityBlockNotFixed <= 0) {
		OutputFiles::m_logFile << "Error : Total number of modifiable resisitivity value is zero or negative !! : " << m_numResistivityBlockNotFixed << std::endl;
		exit(1);
	}

#ifdef _DEBUG_WRITE
	for (int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl) {
		std::cout << " iMdl m_modelID2blockID[iMdl] : " << iMdl << " " << m_modelID2blockID[iMdl] << std::endl;
	}
#endif

}

// Get number of unfixed resistivity parameters
int ResistivityBlockIsotropic::getNumberOfUnfixedResistivityParameters() const {
	return m_numResistivityBlockNotFixed;
}

// Get resistivity values from resisitivity block ID
double ResistivityBlockIsotropic::getResistivityValuesFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValues[ iblk ];
}

// Get previous resistivity values from resisitivity block ID
double ResistivityBlockIsotropic::getResistivityValuesPreFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValuesPre[ iblk ] < 0 ){
		return 1.0e+20;
	}

	return m_resistivityValuesPre[ iblk ];

}

// Get conductivity values from resisitivity block ID
double ResistivityBlockIsotropic::getConductivityValuesFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	if( m_resistivityValues[ iblk ] < 0 ){
		return 0.0;
	}

	return 1.0/m_resistivityValues[ iblk ];
}

// Get resistivity values from element ID
double ResistivityBlockIsotropic::getResistivityValuesFromElemID( const int ielem ) const{

	const int iblk = getBlockIDFromElemID(ielem);
	return getResistivityValuesFromBlockID( iblk );
}

// Get conductivity values from element ID
double ResistivityBlockIsotropic::getConductivityValuesFromElemID( const int ielem ) const{

	const int iblk = getBlockIDFromElemID(ielem);
	return getConductivityValuesFromBlockID( iblk );
}

// Get model ID from block ID
int ResistivityBlockIsotropic::getModelIDFromBlockID( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_blockID2modelID[ iblk ];
}

// Get block ID from model ID
int ResistivityBlockIsotropic::getBlockIDFromModelID( const int imdl ) const{

	assert( imdl >= 0 );
	assert( imdl < m_numResistivityBlockNotFixed );

	return m_modelID2blockID[ imdl ];

}

// Get flag specifing whether resistivity value of each block is fixed or not
bool ResistivityBlockIsotropic::isFixedResistivityValue( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_fixResistivityValues[iblk];
}

// Get flag specifing whether resistivity block is excluded from roughing matrix
bool ResistivityBlockIsotropic::isolated( const int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	return m_isolated[iblk];

}

// Calculate volume of the specified resistivity block
double ResistivityBlockIsotropic::calcVolumeOfBlock( int iblk ) const{

	assert( iblk >= 0 );
	assert( iblk < m_numResistivityBlockTotal );

	const MeshData* const ptrMeshData = ( AnalysisControl::getInstance() )->getPointerOfMeshData();

	double volume(0.0);

	const std::vector<int>& blk2Elem = getBlockID2Elements(iblk);
	for(  std::vector<int>::const_iterator itr = blk2Elem.begin(); itr != blk2Elem.end(); ++itr ){
		const int elemID = *itr;
		volume += ptrMeshData->calcVolume( elemID );
	}

	return volume;
}

// Calculate pre-degenerated roughning matrix
void ResistivityBlockIsotropic::calcRougheningMatrix(){

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

	if((AnalysisControl::getInstance())->isBottomResistivityIncluded()){
		// Calculate roughning matrix from user-defined roughning factor
		addBottomResistivityContribution();
	}

	if((AnalysisControl::getInstance())->isSmallValueToRougheningMatrixDiagonals()){
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
void ResistivityBlockIsotropic::calcRougheningMatrixDegeneratedForLaplacianFilter( DoubleSparseMatrix& rougheningMatrixDegenerated, const double factor ) const{

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
void ResistivityBlockIsotropic::calcRougheningMatrixDegeneratedForDifferenceFilter( const double factor,
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
			const double log10rhoNeib = log10((AnalysisControl::getInstance())->getBottomResistivity());
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

// Calculate array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 ) from common logarithm
void ResistivityBlockIsotropic::calcResistivityUpdatedFullFromLog10ResistivityIncres( const double* const log10resistivity ){

	const int typeBoundConstraints = (AnalysisControl::getInstance())->getTypeBoundConstraints();

	int iMdl(0);
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		if( isFixedResistivityValue(iBlk) ){// If resistivity value is fixed
			m_resistivityValuesUpdatedFull[iBlk] = m_resistivityValuesPre[iBlk];
			continue;
		}
		assert( iMdl == getModelIDFromBlockID(iBlk) );
		if( typeBoundConstraints == ResistivityBlockIsotropic::TRANSFORMING_METHOD ){
			const double jacobian = calcDerivativeXWithRespectToLog10Resistivity(iBlk);// dx / dm
			const double xUpdatedFull = calcTransformedModelParameterFromResistivity(iBlk, m_resistivityValuesPre[iBlk]) + jacobian * log10resistivity[iMdl];
			m_resistivityValuesUpdatedFull[iBlk] = calcResistivityFromTransformedModelParameter(iBlk, xUpdatedFull);
		}
		else if( typeBoundConstraints == ResistivityBlockIsotropic::SIMPLE_BOUND_CONSTRAINING ){
			const double log10ResistivityIncreUpdatedFull = log10( m_resistivityValuesPre[iBlk] ) + log10resistivity[iMdl];
			m_resistivityValuesUpdatedFull[iBlk] = pow( 10.0, log10ResistivityIncreUpdatedFull );
		}
		else{
			OutputFiles::m_logFile << "Error : Wrong type of bound constraining metohd !! : " << typeBoundConstraints << " ." << std::endl;
			exit(1);
		}
		++iMdl;
	}
	
	assert( iMdl == m_numResistivityBlockNotFixed );
}

// Copy derivatives of logarithm of resistivities with respect to transformed model parameter x
void ResistivityBlockIsotropic::copyDerivativeLog10ResistivityWithRespectToX( double* derivs ) const{

	assert(!(AnalysisControl::getInstance())->isAnisotropyConsidered());

	const int numBlockNotFixed = getNumberOfUnfixedResistivityParameters();
	for( int iMdl = 0; iMdl < numBlockNotFixed; ++iMdl ){
		derivs[iMdl] = calcDerivativeLog10ResistivityWithRespectToX( getBlockIDFromModelID(iMdl) );
	}

}

// Update isotropic resistivity
void ResistivityBlockIsotropic::updateResistivityValues(){

	const double stepLengthDampingFactor = ( AnalysisControl::getInstance() )->getStepLengthDampingFactorCur();

	const int typeBoundConstraints = (AnalysisControl::getInstance())->getTypeBoundConstraints();
	if( typeBoundConstraints == ResistivityBlockIsotropic::TRANSFORMING_METHOD ){
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
	else if( typeBoundConstraints == ResistivityBlockIsotropic::SIMPLE_BOUND_CONSTRAINING ){
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
		OutputFiles::m_logFile << "Error : Wrong type of bound constraining metohd !! : " << typeBoundConstraints << " ." << std::endl;
		exit(1);
	}

}

// Output data of resisitivity block model to file
void ResistivityBlockIsotropic::outputResistivityValues(const int iterNum, FILE* fp) const{

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		fprintf(fp, "%10d%5s%15e%15e%15e%15e%10d\n", iBlk, "     ",
			m_resistivityValues[iBlk], m_resistivityValuesMin[iBlk], m_resistivityValuesMax[iBlk], m_weightingConstants[iBlk],
			getTypeOfResistivityBlock(m_fixResistivityValues[iBlk], m_isolated[iBlk]));
	}

}

// Output resistivity values to VTK file
void ResistivityBlockIsotropic::outputResistivityValuesToVTK() const {

	if (!OutputFiles::m_vtkFile.is_open()) {
		return;
	}

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->doesOutputToVTK(AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK)) {

		const int nElem = ((AnalysisControl::getInstance())->getPointerOfMeshData())->getNumElemTotal();

		OutputFiles::m_vtkFile << "SCALARS Resistivity[Ohm-m] float" << std::endl;
		OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
		for (int iElem = 0; iElem < nElem; ++iElem) {
			OutputFiles::m_vtkFile << m_resistivityValues[m_elementID2blockID[iElem]] << std::endl;
		}

	}

}

// Output resistivity values to binary file
void ResistivityBlockIsotropic::outputResistivityValuesToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();

	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Resistivity.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Resistivity[Ohm-m]";
	strcpy(line, ossTitle.str().c_str());
	fout.write(line, 80);

	strcpy(line, "part");
	fout.write(line, 80);

	int ibuf(1);
	fout.write((char*)&ibuf, sizeof(int));

	if (ptrAnalysisControl->getTypeOfMesh() == MeshData::TETRA) {
		strcpy(line, "tetra4");
	}
	else {
		strcpy(line, "hexa8");
	}
	fout.write(line, 80);

	const MeshData* const ptrMeshData = ptrAnalysisControl->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();

	for (int iElem = 0; iElem < nElem; ++iElem) {
		float dbuf = static_cast<float>(m_resistivityValues[m_elementID2blockID[iElem]]);
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();
}

// Copy current isotropic resistivity parameters to current ones
void ResistivityBlockIsotropic::copyUnfixedResistivityParametersCurToVector(double* vector) const {

	for (int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl) {
		vector[iMdl] = log10(m_resistivityValues[m_modelID2blockID[iMdl]]);
	}

}

// Copy previous isotropic resistivity parameters to current ones
void ResistivityBlockIsotropic::copyUnfixedResistivityParametersPreToVector(double* vector) const {

	for (int iMdl = 0; iMdl < m_numResistivityBlockNotFixed; ++iMdl) {
		vector[iMdl] = log10(m_resistivityValuesPre[m_modelID2blockID[iMdl]]);
	}

}

// Copy current isotropic resistivity parameters to previous ones
void ResistivityBlockIsotropic::copyUnfixedResistivityValuesCurToPre() {

	memcpy(m_resistivityValuesPre, m_resistivityValues, sizeof(double) * (m_numResistivityBlockTotal));

}

// Output resistivity data to VTK file
void ResistivityBlockIsotropic::outputResistivityDataToVTK() const {

	outputResistivityBlockIDsToVTK();

}

// Output resistivity data to binary file
void ResistivityBlockIsotropic::outputResistivityDataToBinary() const {

	outputResistivityBlockIDsToBinary();

}

// Calculate model roughness for laplacian filter
double ResistivityBlockIsotropic::calcModelRoughnessForLaplacianFilter() const{

	double* modelVec = new double[m_numResistivityBlockTotal];

	copyResistivityValuesToVectorLog10( modelVec );

	const double roughness = m_rougheningMatrix.calcModelRoughness( modelVec );

	delete[] modelVec;

	return roughness;

}

// Calculate model roughness for difference filter
double ResistivityBlockIsotropic::calcModelRoughnessForDifferenceFilter() const{

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
			const double log10rhoNeib = log10((AnalysisControl::getInstance())->getBottomResistivity());
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
double ResistivityBlockIsotropic::calcModelRoughnessAtBottom() const{

	if( !(AnalysisControl::getInstance())->isBottomResistivityIncluded() ){
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
		const double diff = (AnalysisControl::getInstance())->getRoughningFactorAtBottom() * (log10((AnalysisControl::getInstance())->getBottomResistivity()) - log10(m_resistivityValues[iBlk]));
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
void ResistivityBlockIsotropic::copyResistivityValuesToVectorLog10( double* vector ) const{

	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		vector[iBlk] = log10( m_resistivityValues[iBlk] );
	}

}

// Calculate weighting factor
double ResistivityBlockIsotropic::calWeightingFactor( const double alphaX, const double alphaY, const double alphaZ, const int iElem1, const int iElem2 ) const{
	const CommonParameters::locationXYZ diff = ( (AnalysisControl::getInstance())->getPointerOfMeshData() )->calDiffOfCenters(iElem1, iElem2);
	const double numer = sqrt( pow(alphaX*diff.X, 2) + pow(alphaY*diff.Y, 2) + pow(alphaZ*diff.Z, 2) ) * 1.0e-3;
	const double denom = sqrt( pow(diff.X, 2) + pow(diff.Y, 2) + pow(diff.Z, 2) ) * 1.0e-3;
	return numer / pow(denom, (AnalysisControl::getInstance())->getInverseDistanceWeightingFactor());
}

// Get type of resistivity block
int ResistivityBlockIsotropic::getTypeOfResistivityBlock( const bool fixed, const bool isolated ) const{
	if( fixed ){
		if( isolated ){
			return ResistivityBlockIsotropic::FIXED_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlockIsotropic::FIXED_AND_CONSTRAINED;
		}
	}
	else{// Free
		if( isolated ){
			return ResistivityBlockIsotropic::FREE_AND_ISOLATED;
		}
		else{// Constrained
			return ResistivityBlockIsotropic::FREE_AND_CONSTRAINED;
		}
	}
}

// Calculate roughning matrix with using elements share faces 
void ResistivityBlockIsotropic::calcRougheningMatrixUsingElementsShareFaces( const double factor ){

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
void ResistivityBlockIsotropic::calcRougheningMatrixUsingElementsShareFacesWeightingByAreaVolumeRatio( const double factor ){

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
	const double inverseDistanceWeightingFactor = (AnalysisControl::getInstance())->getInverseDistanceWeightingFactor();
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
void ResistivityBlockIsotropic::calcRougheningMatrixUsingResistivityBlocksShareFaces( const double factor ){

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
void ResistivityBlockIsotropic::calcRougheningMatrixUserDefined( const double factor ){

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
void ResistivityBlockIsotropic::addBottomResistivityContribution(){

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

	const double inverseDistanceWeightingFactor = (AnalysisControl::getInstance())->getInverseDistanceWeightingFactor();
	const int numElemBot = ptrMeshData->getNumElemOnBoundaryPlanes(MeshData::XYPlus);
	std::map< int, int> blockoCounter;
	std::map< int, double> blockToWeight;
	for( int iElemBot = 0; iElemBot < numElemBot; ++iElemBot ){
		const int iElem = ptrMeshData->getElemBoundaryPlanes(MeshData::XYPlus, iElemBot);
		const int iBlk = getBlockIDFromElemID( iElem );
		if( isolated(iBlk) ){
			continue;
		}
		double factor = (AnalysisControl::getInstance())->getRoughningFactorAtBottom();
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
			m_rougheningMatrix.addRightHandSideVector( iBlk, factor * log10((AnalysisControl::getInstance())->getBottomResistivity()) );
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
			m_rougheningMatrix.addRightHandSideVector( iBlk, weightedFactorAveraged * log10((AnalysisControl::getInstance())->getBottomResistivity()) );
		}
	}

	if( volumes != NULL ){
		delete [] volumes;
	}

}

// Add small value to diagonals
void ResistivityBlockIsotropic::addSmallValueToDiagonals(){

	const double value = (AnalysisControl::getInstance())->getSmallValueAddedToDiagonals();
	for( int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk ){
		m_rougheningMatrix.setStructureAndAddValueByTripletFormat( iBlk, iBlk, value );
	}

}

// Calculate transformed model parameter x from resistivity
double ResistivityBlockIsotropic::calcTransformedModelParameterFromResistivity( const int iBlk, const double resistivity ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD);
	return calcTransformedModelParameterFromLog10Resistivity(iBlk, log10(resistivity));
}

// Calculate transformed model parameter x from common logarithm of resistivity
double ResistivityBlockIsotropic::calcTransformedModelParameterFromLog10Resistivity( const int iBlk, const double log10Resistivity ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD );
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
double ResistivityBlockIsotropic::calcResistivityFromTransformedModelParameter( const int iBlk, const double x ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD );
	return pow( 10.0, calcLog10ResistivityFromTransformedModelParameter(iBlk, x) );
}

// Calculate common logarithm of resistivity from transformed model parameter x
double ResistivityBlockIsotropic::calcLog10ResistivityFromTransformedModelParameter( const int iBlk, const double x ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD );
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];
	return 0.5 * ( b - a ) * tanh( 0.5 * n * x ) + 0.5 * ( b + a );
}

// Calculate derivative of logarithm of resistivity with respect to transformed model parameter x
double ResistivityBlockIsotropic::calcDerivativeLog10ResistivityWithRespectToX( const int iBlk ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD );
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
double ResistivityBlockIsotropic::calcDerivativeXWithRespectToLog10Resistivity( const int iBlk ) const{
	assert((AnalysisControl::getInstance())->getTypeBoundConstraints() == ResistivityBlockIsotropic::TRANSFORMING_METHOD );
	const double m = log10(m_resistivityValuesPre[iBlk]);
	const double a = log10(m_resistivityValuesMin[iBlk]);
	const double b = log10(m_resistivityValuesMax[iBlk]);
	const double n = m_weightingConstants[iBlk];

	assert( !m_fixResistivityValues[iBlk] );

	const double minDistanceToBounds = (AnalysisControl::getInstance())->getMinDistanceToBounds();
	const double val1 = b - m < minDistanceToBounds ? minDistanceToBounds : b - m;
	const double val2 = m - a < minDistanceToBounds ? minDistanceToBounds : m - a;

	return (b - a) / ( n * val1 * val2 );
}
