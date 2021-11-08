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
#ifndef DBLDEF_RESISTIVITY_BLOCK
#define DBLDEF_RESISTIVITY_BLOCK

#include "AnalysisControl.h"
#include "RougheningSquareMatrix.h"

#include <vector>

// Class of resistivity blocks
class ResistivityBlock{

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

	// Return the the instance of the class
    static ResistivityBlock* getInstance();

	// Read data of resisitivity block model from input file
	void inputResisitivityBlock();

	// Get resisitivity block ID from element ID
	inline int getBlockIDFromElemID( const int ielem ) const{
		return m_elementID2blockID[ ielem ];
	};

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

	// Get total number of resistivity blocks
	int getNumResistivityBlockTotal() const;

	// Get number of resistivity blocks whose resistivity values are not fixed
	int getNumResistivityBlockNotFixed() const;

	// Get flag specifing whether bottom resistivity is included in roughning
	bool includeBottomResistivity() const;

	// Get resistivity of the bottom of the model
	double getBottomResistivity() const;

	// Get roughning factor at the bottom of the model
	double getRoughningFactorAtBottom() const;

	// Get flag specifing whether small value is added to diagonals
	bool getFlagAddSmallValueToDiagonals() const;

	// Set small value added to diagonals
	double getSmallValueAddedToDiagonals() const;

	// Get minimum distance between current resistivity and resistivity bounds in common logarithm scale
	double getMinDistanceToBounds() const;

	// Get type of bound constraints
	int getTypeBoundConstraints() const;

	// Get positive real factor of inverse distance weighting
	double getInverseDistanceWeightingFactor() const;

	// Set flag specifing whether bottom resistivity is included in roughning
	void setFlagIncludeBottomResistivity( bool include );

	// Set resistivity of the bottom of the model
	void setBottomResistivity(  const double resistivity );

	// Set roughning factor at the bottom of the model
	void setRoughningFactorAtBottom( const double factor );

	// Set flag specifing whether small value is added to diagonals
	void setFlagAddSmallValueToDiagonals( const bool flag );

	// Set small value added to bottom
	void setSmallValueAddedToDiagonals( const double value );

	// Set minimum distance to resistivity bounds in common logarithm scale
	void setMinDistanceToBounds( const double distance );

	// Set type of bound constraints
	void setTypeBoundConstraints( const int type );

	// Set positive real factor of inverse distance weighting
	void setInverseDistanceWeightingFactor( const double factor );

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

	// Calculate array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 )
	// from common logarithm
	void calctResistivityUpdatedFullFromLog10ResistivityIncres( const double* const log10resistivity );

	// Copy derivatives of logarithm of resistivities with respect to transformed model parameter x
	void copyDerivativeLog10ResistivityWithRespectToX( double* derivs ) const;

	// Change resistivity values
	void updateResistivityValues();

	// Output resistivity values to VTK file
	void outputResistivityValuesToVTK() const;

	// Output resistivity values to binary file
	void outputResistivityValuesToBinary( const int iterNum ) const;

	// Output data of resisitivity block model to file
	void outputResisitivityBlock( const int iterNum ) const;

	// Copy common logarithm of free resistivity values to vector
	void copyResistivityValuesNotFixedToVectorLog10( double* vector ) const;

	// Copy common logarithm of previous free resistivity values to vector
	void copyResistivityValuesNotFixedPreToVectorLog10( double* vector ) const;

	// Copy common logarithm of current free resistivity values to previous ones
	void copyResistivityValuesNotFixedCurToPre() const;

	// Calculate model roughness for laplacian filter
	double calcModelRoughnessForLaplacianFilter() const;

	// Calculate model roughness for difference filter
	double calcModelRoughnessForDifferenceFilter() const;

	// Calculate model roughness at the bottom
	double calcModelRoughnessAtBottom() const;

	// Copy common logarithm of resistivity values
	void copyResistivityValuesToVectorLog10( double* vector ) const;

	//// Add vector resulting from the product of sensitivity matrix and model vector
	//void addProductOfSensitivityMatrixAndModelVector( const double senMat, const int numData, double* vecOut ) const;

	//// Get resistivity model
	//const double* const getResistivityValues() const;

	// Output resistivity data to VTK file
	void outputResistivityDataToVTK() const;

	// Output resistivity data to binary file
	void outputResistivityDataToBinary() const;

	//// Get arrays of elements belonging to each resistivity block
	//const std::vector<int>&  getBlockID2Elements( const int iBlk ) const;
	// Get arrays of elements belonging to each resistivity model
	//const std::vector<int>&  getBlockID2Elements( const int iBlk ) const;
	const std::vector< std::pair<int,double> >&  getBlockID2Elements( const int iBlk ) const;
	
#ifdef _ANISOTOROPY
	//// Calculate anisotropy coefficient rotated
	//void calcAisotropyCoefficientRotated( const int iBlk, CommonParameters::Vector3D& coeffX, CommonParameters::Vector3D& coeffY, CommonParameters::Vector3D& coeffZ ) const;

	// Calculate anisotropic conductivity tensor
	//void calcAisotropicConductivityTensor( const int blockID, CommonParameters::Vector3D& matX, CommonParameters::Vector3D& matY, CommonParameters::Vector3D& matZ ) const;
	void calcAisotropicConductivityTensor( const int blockID, double conductivityTensor[3][3] ) const;
#endif

private:

	// Constructer
	ResistivityBlock();

	// Destructer
	~ResistivityBlock();

	// Copy constructer
	ResistivityBlock(const ResistivityBlock& rhs){
		std::cerr << "Error : Copy constructer of the class ResistivityBlock is not implemented." << std::endl;
		exit(1);
	};

	// Assignment operator
	ResistivityBlock& operator=(const ResistivityBlock& rhs){
		std::cerr << "Error : Assignment operator of the class ResistivityBlock is not implemented." << std::endl;
		exit(1);
	};

	// Array of the resistivity block IDs of each element
	int* m_elementID2blockID;

	// Array convert IDs of resistivity block to model IDs 
	int* m_blockID2modelID;

	// Array convert model IDs to IDs of resistivity block  
	int* m_modelID2blockID;

	// Total number of resistivity blocks
	int m_numResistivityBlockTotal;

	// Number of resistivity blocks whose resistivity values are fixed
	int m_numResistivityBlockNotFixed;

	// Array of resistivity values of each block
	double* m_resistivityValues;

	// Array of previous resistivity values of each block
	double* m_resistivityValuesPre;

	// Array of resistivity values obtained by inversion which is the ones fully updated ( damping factor = 1 )
	double* m_resistivityValuesUpdatedFull;

	// Array of minimum resistivity values of each block
	double* m_resistivityValuesMin;

	// Array of maximum resistivity values of each block
	double* m_resistivityValuesMax;

	// Positive constant parameter n
	double* m_weightingConstants;

	// Flag specifing whether resistivity value of each block is fixed or not
	bool* m_fixResistivityValues;

	// Flag specifing whether resistivity block is excluded from roughing matrix
	bool* m_isolated;

#ifdef _ANISOTOROPY
	// Anisotoropy coefficient
	//CommonParameters::Vector3D* m_anisotoropyCoeff;

	// Array mapping block IDs to indexes of the blocks where anisotropic conductivity tensor components are specified
	std::map<int, int> m_mapBlockIDWithAnisotropyToIndex;

	// Array of the current components of resistivity tensor of axial anisotropy
	std::vector<CommonParameters::Vector3D> m_resistivityValuesAxialAnisotropy;

	// Array of the previous components of resistivity tensor of axial anisotropy
	std::vector<CommonParameters::Vector3D> m_resistivityValuesAxialAnisotropyPre;

	// Array of the current XX components of resistivity tensor obtained by inversion which is the ones fully updated ( damping factor = 1 )
	std::vector<CommonParameters::Vector3D> m_resistivityValuesAxialAnisotropyFull;

	// Axial anisotropy strike angles
	std::vector<double> m_axialAnisotropyStrileAngle;

	// Axial anisotropy dip angles
	std::vector<double> m_axialAnisotropyDipAngle;

	// Axial anisotropy slant angles
	std::vector<double> m_axialAnisotropySlantAngle;
#endif

	// Roughening matrix
	RougheningSquareMatrix m_rougheningMatrix;

	// Arrays of elements belonging to each resistivity model
	std::vector< std::pair<int,double> >* m_blockID2Elements;

	// Flag specifing whether bottom resistivity is included in roughning
	bool m_includeBottomResistivity;

	// Resistivity of the bottom of the model
	double m_bottomResistivity;

	// Roughning factor at the bottom of the model
	double m_roughningFactorAtBottom; 

	// Flag specifing whether small value is added to diagonals
	bool m_addSmallValueToDiagonals;

	// Small value added to bottom
	double m_smallValueAddedToDiagonals;

	// Type of bound constraints
	int m_typeBoundConstraints;

	// Minimum distance between current resistivity and resistivity bounds in common logarithm scale
	double m_minDistanceToBounds;

	// Positive real factor of inverse distance weighting
	double m_inverseDistanceWeightingFactor;

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

};

#endif
