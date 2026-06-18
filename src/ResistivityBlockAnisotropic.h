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
#ifndef DBLDEF_RESISTIVITY_BLOCK_ANISOTROPIC
#define DBLDEF_RESISTIVITY_BLOCK_ANISOTROPIC

#include "ResistivityBlock.h"
#include "CommonParameters.h"

#include <iostream>

// Class of anisotropic resistivity blocks
class ResistivityBlockAnisotropic : public ResistivityBlock {

public:

	enum AnisotropyTypes{
		ISOTROPY = 0,
		TRANSVERSE_ISOTROPY,
		GENERAL_ANISOTROPY,
		END_OF_ANISOTROPY_TYPE// THIS MUST BE WRITTEN AT THE END
	};

	enum AnisotropyParameters {
		RHO_XX = 0,
		RHO_YY,
		RHO_ZZ,
		STRIKE,
		DIP,
		SLANT,
		END_OF_ANISOTROPY_PARAMETERS// THIS MUST BE WRITTEN AT THE END
	};

	struct AnistropicResistivityParameters {
		double rhoXX;
		double rhoYY;
		double rhoZZ;
		double strike;
		double dip;
		double slant;
	};

	// Constructer
	ResistivityBlockAnisotropic();

	// Destructer
	virtual ~ResistivityBlockAnisotropic();

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

	// Get anisotropy type
	int getTypeOfAnisotropy(const int iblk) const;

	// Get flag specifing whether anisotrpic resistivity parameters is fixed or not
	bool isFixedAnisotropicResistivityParameters(const int iblk, const int iparam) const;

	// Get indexes of block and anistropic resistivity parameter from model index
	std::pair<int, int> getIndexesOfBlockAndAnisotropicResistivityParameterFromModelIndex(const int iMdl) const;

	// Calculate component-wise spatial difference filter for transversely isotropic conductivity
	void calcSpatialDifferenceFilterForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
		std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const;

	// Calculate model roughness for transversely isotropic conductivity
	double calcModelRoughnessForTransverselyIsotropicConductivity() const;

	// Calculate difference filter between sigma_x and sigma_y for transversely isotropic conductivity
	void calcDifferenceFilterBetweenComponentsForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
		std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const;

	// Calculate L2 norm of the differences between sigma_x and sigma_y for transversely isotropic conductivity
	double calcL2NormOfDifferencesBetweenResistivityComponentsForTransverselyIsotropicConductivity() const;

	// Calculate cross product for transversely isotropic conductivity
	void calcCrossProductForTransverselyIsotropicConductivity(const double factor, std::vector<double>& values) const;

	// Calculate L2 norm of cross product for transversely isotropic conductivity
	double calcL2NormOfCrossProductForTransverselyIsotropicConductivity() const;

	// Calculate derivative of the cross product for transversely isotropic conductivity
	void calcDerivativeOfCrossProductForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
		std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const;

	// Confirm that there is no element having general anisotropy
	void confirmNoElementHavingGeneralAnisotropy() const;

	// Calculate array of fully updated anisotropic resistivity parameters obtained by inversion
	void calcUnfixedAnisotropicResistivityParametersUpdatedFull(const double* const updatedModel);

	// Calculate anisotropic conductivity tensor
	void calcAnisotropicConductivityTensor( const int blockID, double conductivityTensor[3][3] ) const;

	// Calculate derivative of anisotropic conductivity tensor
	void calcDerivativeOfAnisotropicConductivityTensor(const int blockID, const int paramID, double derivative[3][3]) const;

	// Get model ID from block ID and anisotropic parameter
	int getModelIDFromBlockIDAndAnisotropicParameter(const int iblk, const int iparam) const;

private:

	// Copy constructer
	ResistivityBlockAnisotropic(const ResistivityBlockAnisotropic& rhs){
		std::cerr << "Error : Copy constructer of the class ResistivityBlockAnisotropic is not implemented." << std::endl;
		exit(1);
	};

	// Assignment operator
	ResistivityBlockAnisotropic& operator=(const ResistivityBlockAnisotropic& rhs){
		std::cerr << "Error : Assignment operator of the class ResistivityBlockAnisotropic is not implemented." << std::endl;
		exit(1);
	};

	// Number of anisotropic resistivity parameters whose values are NOT fixed
	int m_numAnisotropicResistivityParametersNotFixed;

	// Array convert index of anistropic resistivity patameter to model IDs 
	std::vector<int> m_anisotropicResistivityParameter2modelID[6];

	// Array convert model IDs to index of anistropic resistivity parameter  
	std::vector< std::pair<int, int> > m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes;

	// Array of the anisotropy type
	std::vector<int> m_typeOfAnisotropy;

	// Array of the current components of anistropic resistivity parameters
	std::vector<AnistropicResistivityParameters> m_anisotropicResistivityParameters;

	// Array of the previous components of anistropic resistivity parameters
	std::vector<AnistropicResistivityParameters> m_anisotropicResistivityParametersPre;

	// Array of the current components of resistivity tensor given that fully updates are perfomed in the inversion ( damping factor = 1 )
	std::vector<AnistropicResistivityParameters> m_anisotropicResistivityParametersUpdatedFull;

	// Array of the flag specifing whether anisotropic resistivity parameters are fixed or not
	std::vector<int> m_fixedAnisotropicResistivityParameters[6];

	// Calculate derivative of transverse isotropic conductivity tensor
	void calcDerivativeOfTransverseIsotropicConductivityTensor(const int paramID, const AnistropicResistivityParameters& anisotropicResistivityParam, double derivative[3][3]) const;

	// Calculate derivative of general anisotropic conductivity tensor
	void calcDerivativeOfGeneralAnisotropicConductivityTensor(const int paramID, const AnistropicResistivityParameters& anisotropicResistivityParam, double derivative[3][3]) const;

	// Auxiliary function for calculating derivative of anisotropic conductivity tensor
	void calcDerivativeOfAnisotropicConductivityTensorAux(const double lhsMatrix[3][3], double matrix[3][3], const double transposeOfRhsMatrix[3][3]) const;

	// Read reslstivity values from input file
	virtual void inputResistivityValues(const int nElem, std::ifstream& inFile);

	// Output resistivity values to resistivity_block_iter*.dat
	virtual void outputResistivityValues(const int iterNum, FILE* fp) const;

	// Output xx-component of the anisotropic resistivity tensor to binary file
	void outputRhoXXOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output yy-component of the anisotropic resistivity tensor to binary file
	void outputRhoYYOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output zz-component of the anisotropic resistivity tensor to binary file
	void outputRhoZZOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output strike of the anisotropic resistivity tensor to binary file
	void outputStrikeOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output dip of the anisotropic resistivity tensor to binary file
	void outputDipOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output slant of the anisotropic resistivity tensor to binary file
	void outputSlantOfAnisotropicResistivityTensorToBinary(const int iterNum) const;

	// Output abs(log10(rho_xx)-log10(rho_yy)) to binary file
	void outputfAnisotropyIndicatorToBinary(const int iterNum) const;

	// Output anisotropy types to VTK file
	void outputAnisotropyTypesToVTK() const;

	// Output anisotropy types to binary file
	void outputAnisotropyTypesToBinary() const;

	// Get anisotropic resistivity parameter from resisitivity block ID
	double getAnisotropicResistivityParameterFromBlockID(const int iblk, const int iparam) const;

	// Get anisotropic previous resistivity parameter from resisitivity block ID
	double getAnisotropicResistivityParameterPreFromBlockID(const int iblk, const int iparam) const;

	// Get anisotropic resistivity parameters from resisitivity block ID
	AnistropicResistivityParameters getAnisotropicResistivityParametersFromBlockID(const int iblk) const;

	// Get anisotropic previous resistivity parameters from resisitivity block ID
	AnistropicResistivityParameters getAnisotropicResistivityParametersPreFromBlockID(const int iblk) const;

	// Get anisotropic fully-updated resistivity parameters from resisitivity block ID
	AnistropicResistivityParameters getAnisotropicResistivityParametersUpdatedFullFromBlockID(const int iblk) const;

	// Calculate cross product for transversely isotropic conductivity
	CommonParameters::Vector3D calcCrossProductForTransverselyIsotropicConductivity(const AnistropicResistivityParameters& anisoParamsA, const AnistropicResistivityParameters& anisoParamsB) const;

	// Auxiliary function to calculate derivative of the cross product for transversely isotropic conductivity
	void calcDerivativeOfCrossProductForTransverselyIsotropicConductivityAux(const int iBlk, const int iBlkNeib,
		const AnistropicResistivityParameters& anisoParamsA, const AnistropicResistivityParameters& anisoParamsB, const double factor, std::vector<int> cols[3], std::vector<double> vals[3]) const;

};

#endif
