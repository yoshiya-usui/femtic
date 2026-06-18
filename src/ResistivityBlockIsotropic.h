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
#ifndef DBLDEF_RESISTIVITY_BLOCK_ISOTROPIC
#define DBLDEF_RESISTIVITY_BLOCK_ISOTROPIC

#include "ResistivityBlock.h"
#include "RougheningSquareMatrix.h"

#include <iostream>

// Class of isotropic resistivity blocks
class ResistivityBlockIsotropic : public ResistivityBlock {

public:

	enum ResistivityBlockTypes{
		FREE_AND_CONSTRAINED = 0,
		FIXED_AND_ISOLATED,
		FIXED_AND_CONSTRAINED,
		FREE_AND_ISOLATED,
	};

	enum BoundconstrainingTypes{
		SIMPLE_BOUND_CONSTRAINING = 0,
		TRANSFORMING_METHOD,
	};

	// Constructer
	ResistivityBlockIsotropic();

	// Destructer
	virtual ~ResistivityBlockIsotropic();

	// Get number of unfixed resistivity parameters
	virtual int getNumberOfUnfixedResistivityParameters() const;

	// Change resistivity values
	virtual void updateResistivityValues();

	// Output resistivity values to VTK file
	virtual void outputResistivityValuesToVTK() const;

	// Output resistivity values to binary file
	virtual void outputResistivityValuesToBinary(const int iterNum) const;

	// Copy current resistivity parameters to  a vector
	virtual void copyUnfixedResistivityParametersCurToVector(double* vector) const;

	// Copy previous resistivity parameters to a vector
	virtual void copyUnfixedResistivityParametersPreToVector(double* vector) const;

	// Copy current resistivity parameters to previous ones
	virtual void copyUnfixedResistivityValuesCurToPre();

	// Output resistivity data to VTK file
	virtual void outputResistivityDataToVTK() const;

	// Output resistivity data to binary file
	virtual void outputResistivityDataToBinary() const;

	// Get resistivity values from resisitivity block ID
	double getResistivityValuesFromBlockID( const int iblk ) const;

	// Get previous resistivity values from resisitivity block ID
	double getResistivityValuesPreFromBlockID( const int iblk ) const;

	// Get conductivity values from resisitivity block ID
	double getConductivityValuesFromBlockID( const int iblk ) const;

	// Get resistivity values from element ID
	double getResistivityValuesFromElemID( const int ielem ) const;

	// Get conductivity values from element ID
	double getConductivityValuesFromElemID( const int ielem ) const;

	// Get model ID from block ID
	int getModelIDFromBlockID( const int iblk ) const;

	// Get block ID from model ID
	int getBlockIDFromModelID( const int imdl ) const;

	// Get flag specifing whether resistivity value of each block is fixed or not
	bool isFixedResistivityValue( const int iblk ) const;

	// Get flag specifing whether resistivity block is excluded from roughing matrix
	bool isolated( const int iblk ) const;

	// Calculate volume of the specified resistivity block
	double calcVolumeOfBlock( int iblk ) const;

	// Calculate pre-degenerated roughning matrix
	void calcRougheningMatrix();

	// Calculate roughning matrix degenerated for laplacian filter
	void calcRougheningMatrixDegeneratedForLaplacianFilter( DoubleSparseMatrix& rougheningMatrixDegenerated, const double factor ) const;

	// Calculate roughning matrix degenerated for difference filter
	void calcRougheningMatrixDegeneratedForDifferenceFilter( const double factor,
		std::vector< std::pair<int,int> >& nonZeroCols, std::vector<double>& matValues, std::vector<double>& rhsValues ) const;

	// Calculate array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 ) from common logarithm
	void calcResistivityUpdatedFullFromLog10ResistivityIncres( const double* const log10resistivity );

	// Copy derivatives of logarithm of resistivities with respect to transformed model parameter x
	void copyDerivativeLog10ResistivityWithRespectToX( double* derivs ) const;

	// Calculate model roughness for laplacian filter
	double calcModelRoughnessForLaplacianFilter() const;

	// Calculate model roughness for difference filter
	double calcModelRoughnessForDifferenceFilter() const;

	// Calculate model roughness at the bottom
	double calcModelRoughnessAtBottom() const;

private:

	// Copy constructer
	ResistivityBlockIsotropic(const ResistivityBlockIsotropic& rhs){
		std::cerr << "Error : Copy constructer of the class ResistivityBlockIsotropic is not implemented." << std::endl;
		exit(1);
	};

	// Assignment operator
	ResistivityBlockIsotropic& operator=(const ResistivityBlockIsotropic& rhs){
		std::cerr << "Error : Assignment operator of the class ResistivityBlockIsotropic is not implemented." << std::endl;
		exit(1);
	};

	// Number of resistivity blocks whose resistivity is NOT fixed
	int m_numResistivityBlockNotFixed;

	// Array of resistivity values of each block
	double* m_resistivityValues;

	// Array of previous resistivity values of each block
	double* m_resistivityValuesPre;

	// Array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 )
	double* m_resistivityValuesUpdatedFull;

	// Positive constant parameter n
	double* m_weightingConstants;

	// Flag specifing whether resistivity value of each block is fixed or not
	bool* m_fixResistivityValues;

	// Flag specifing whether resistivity block is excluded from roughing matrix
	bool* m_isolated;

	// Roughening matrix
	RougheningSquareMatrix m_rougheningMatrix;

	// Read reslstivity values from input file
	virtual void inputResistivityValues(const int nElem, std::ifstream& inFile);

	// Output resistivity values to resistivity_block_iter*.dat
	virtual void outputResistivityValues(const int iterNum, FILE* fp) const;

	// Calculate weighting factor
	double calWeightingFactor( const double alphaX, const double alphaY, const double alphaZ, const int iElem1, const int iElem2 ) const;

	// Get type of resistivity block
	int getTypeOfResistivityBlock( const bool fixed, const bool isolated ) const;

	// Calculate roughening matrix with using elements share faces 
	void calcRougheningMatrixUsingElementsShareFaces( const double factor );

	// Calculate roughening matrix with using elements share faces using area-volume ratio as weights
	void calcRougheningMatrixUsingElementsShareFacesWeightingByAreaVolumeRatio( const double factor );

	// Calculate roughening matrix with using resistivty blocks share faces
	void calcRougheningMatrixUsingResistivityBlocksShareFaces( const double factor );

	// Calculate roughning matrix from user-defined roughning factor
	void calcRougheningMatrixUserDefined( const double factor );

	// Add contribution of bottom resistivity to roughning matrix
	void addBottomResistivityContribution();

	// Add small value to diagonals
	void addSmallValueToDiagonals();

	// Calculate transformed model parameter x from resistivity
	double calcTransformedModelParameterFromResistivity( const int iBlk, const double resistivity ) const;

	// Calculate transformed model parameter x from common logarithm of resistivity
	double calcTransformedModelParameterFromLog10Resistivity( const int iBlk, const double log10Resistivity ) const;

	// Calculate resistivity from transformed model parameter x
	double calcResistivityFromTransformedModelParameter( const int iBlk, const double x ) const;

	// Calculate common logarithm of resistivity from transformed model parameter x
	double calcLog10ResistivityFromTransformedModelParameter( const int iBlk, const double x ) const;

	// Calculate derivative of logarithm of resistivity with respect to transformed model parameter x
	double calcDerivativeLog10ResistivityWithRespectToX( const int iBlk ) const;

	// Calculate derivative of transformed model parameter x with respect to logarithm of resistivity
	double calcDerivativeXWithRespectToLog10Resistivity( const int iBlk ) const;

	// Copy common logarithm of resistivity values
	void copyResistivityValuesToVectorLog10(double* vector) const;
};

#endif
