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
#include "ResistivityBlockAnisotropic.h"
#include "AnalysisControl.h"
#include "OutputFiles.h"
#include "Util.h"

#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <vector>
#include <algorithm>


// Constructer
ResistivityBlockAnisotropic::ResistivityBlockAnisotropic():
	m_numAnisotropicResistivityParametersNotFixed(0)
{}

// Destructer
ResistivityBlockAnisotropic::~ResistivityBlockAnisotropic(){
}

// Read reslstivity values from input file
void ResistivityBlockAnisotropic::inputResistivityValues(const int nElem, std::ifstream& inFile) {

	m_anisotropicResistivityParameters.reserve(m_numResistivityBlockTotal);
	m_anisotropicResistivityParametersPre.reserve(m_numResistivityBlockTotal);
	m_anisotropicResistivityParametersUpdatedFull.reserve(m_numResistivityBlockTotal);
	for (int iParam = 0; iParam < 6; ++iParam) {
		m_anisotropicResistivityParameter2modelID[iParam].reserve(m_numResistivityBlockTotal);
		m_fixedAnisotropicResistivityParameters[iParam].reserve(m_numResistivityBlockTotal);
	}

	if (m_resistivityValuesMin != NULL) {
		delete[] m_resistivityValuesMin;
	}
	m_resistivityValuesMin = new double[m_numResistivityBlockTotal];
	if (m_resistivityValuesMax != NULL) {
		delete[] m_resistivityValuesMax;
	}
	m_resistivityValuesMax = new double[m_numResistivityBlockTotal];
	for (int i = 0; i < m_numResistivityBlockTotal; ++i) {
		m_resistivityValuesMin[i] = 0.0;
		m_resistivityValuesMax[i] = 0.0;
	}

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		int idum(0);
		inFile >> idum;
		if (idum != iBlk) {
			OutputFiles::m_logFile << "Error : Block ID is wrong !!" << std::endl;
			exit(1);
		}
		inFile >> idum;
		const int typeOfAnistropy = idum;
		m_typeOfAnisotropy.push_back(typeOfAnistropy);
		AnistropicResistivityParameters anistropicResistivityParam = { -1.0, -1.0, -1.0, 0.0, 0.0, 0.0 };
		int isFixedArray[6] = { -1, -1, -1, -1, -1, -1 };
		if (typeOfAnistropy == ResistivityBlockAnisotropic::ISOTROPY) {
			inFile >> anistropicResistivityParam.rhoXX;
			if (anistropicResistivityParam.rhoXX <= 0.0) {
				OutputFiles::m_logFile << "Error : The resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoXX << std::endl;
				exit(1);
			}
			anistropicResistivityParam.rhoYY = anistropicResistivityParam.rhoXX;
			anistropicResistivityParam.rhoZZ = anistropicResistivityParam.rhoXX;
			inFile >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk];
			if (m_resistivityValuesMin[iBlk] <= 0.0) {
				OutputFiles::m_logFile << "Error : Minimum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMin[iBlk] << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] <= 0.0) {
				OutputFiles::m_logFile << "Error : Maximum resistivity value of block " << iBlk << " is less than or equal to zero !! : " << m_resistivityValuesMax[iBlk] << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoXX) {
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax << " ) is less than initial resistivity ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoXX) {
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax << " ) is greater than initial resistivity ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			int isFixed(0);
			inFile >> isFixed;
			isFixedArray[RHO_XX] = isFixed;
			isFixedArray[RHO_YY] = 1;
			isFixedArray[RHO_ZZ] = 1;
			isFixedArray[STRIKE] = 1;
			isFixedArray[DIP] = 1;
			isFixedArray[SLANT] = 1;
		}
		else if (typeOfAnistropy == TRANSVERSE_ISOTROPY || typeOfAnistropy == GENERAL_ANISOTROPY) {
			// Read principle resistivity values
			inFile >> anistropicResistivityParam.rhoXX >> anistropicResistivityParam.rhoYY;
			if (anistropicResistivityParam.rhoXX <= 0.0) {
				OutputFiles::m_logFile << "Error : The xx-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoXX << std::endl;
				exit(1);
			}
			if (anistropicResistivityParam.rhoYY <= 0.0) {
				OutputFiles::m_logFile << "Error : The yy-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoYY << std::endl;
				exit(1);
			}
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				anistropicResistivityParam.rhoZZ = anistropicResistivityParam.rhoXX;
			}
			else{
				inFile >> anistropicResistivityParam.rhoZZ;
				if (anistropicResistivityParam.rhoZZ <= 0.0) {
					OutputFiles::m_logFile << "Error : The zz-compoonent of the resistivity of block " << iBlk << " is less than or equal to zero !! : " << anistropicResistivityParam.rhoZZ << std::endl;
					exit(1);
				}
			}
			// Read orientation angles
			inFile >> anistropicResistivityParam.strike >> anistropicResistivityParam.dip;
			anistropicResistivityParam.strike *= CommonParameters::deg2rad;
			anistropicResistivityParam.dip *= CommonParameters::deg2rad;
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				anistropicResistivityParam.slant = 0.0;
			}
			else{
				inFile >> anistropicResistivityParam.slant;
				anistropicResistivityParam.slant *= CommonParameters::deg2rad;
			}
			// Read lower and upper bounds of resistivity values
			inFile >> m_resistivityValuesMin[iBlk] >> m_resistivityValuesMax[iBlk];
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoXX) {
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the xx-component ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoXX) {
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the xx-component ( " << anistropicResistivityParam.rhoXX << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoYY) {
				OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the yy-component ( " << anistropicResistivityParam.rhoYY << " )." << std::endl;
				exit(1);
			}
			if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoYY) {
				OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the yy-component ( " << anistropicResistivityParam.rhoYY << " )." << std::endl;
				exit(1);
			}
			if (typeOfAnistropy == GENERAL_ANISOTROPY) {
				if (m_resistivityValuesMax[iBlk] < anistropicResistivityParam.rhoZZ) {
					OutputFiles::m_logFile << "Error : Maximum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is less than initial resistivity of the zz-component ( " << anistropicResistivityParam.rhoZZ << " )." << std::endl;
					exit(1);
				}
				if (m_resistivityValuesMin[iBlk] > anistropicResistivityParam.rhoZZ) {
					OutputFiles::m_logFile << "Error : Minimum resistivity value ( " << m_resistivityValuesMax[iBlk] << " ) is greater than initial resistivity of the zz-component ( " << anistropicResistivityParam.rhoZZ << " )." << std::endl;
					exit(1);
				}
			}
			// Read fixing flags
			if (typeOfAnistropy == TRANSVERSE_ISOTROPY) {
				inFile >> isFixedArray[RHO_XX] >> isFixedArray[RHO_YY];
				isFixedArray[RHO_ZZ] = 1;
				inFile >> isFixedArray[STRIKE] >> isFixedArray[DIP];
				isFixedArray[SLANT] = 1;
			} 
			else {
				inFile >> isFixedArray[RHO_XX] >> isFixedArray[RHO_YY] >> isFixedArray[RHO_ZZ];
				inFile >> isFixedArray[STRIKE] >> isFixedArray[DIP] >> isFixedArray[SLANT];
			}
		}
		else {
			OutputFiles::m_logFile << "Error : Unsupported type of anisotropy (" << typeOfAnistropy << ")" << std::endl;
			exit(1);
		}
		m_anisotropicResistivityParameters.push_back(anistropicResistivityParam);
		m_anisotropicResistivityParametersPre.push_back(anistropicResistivityParam);
		m_anisotropicResistivityParametersUpdatedFull.push_back(anistropicResistivityParam);
		for (int iParam = 0; iParam < 6; ++iParam) {
			m_fixedAnisotropicResistivityParameters[iParam].push_back(isFixedArray[iParam]);
			if (isFixedArray[iParam]) {
				// Fixed
				m_anisotropicResistivityParameter2modelID[iParam].push_back(-1);
			}
			else {
				// Changable
				m_anisotropicResistivityParameter2modelID[iParam].push_back(m_numAnisotropicResistivityParametersNotFixed);
				++m_numAnisotropicResistivityParametersNotFixed;
			}
		}
	}

	inFile.close();

	if (m_fixedAnisotropicResistivityParameters[0][0] < 1) {
		OutputFiles::m_logFile << "Error : Resistivity block 0 must be the air, and its resistivity must be fixed." << std::endl;
		exit(1);
	}
	if (m_numAnisotropicResistivityParametersNotFixed <= 0) {
		OutputFiles::m_logFile << "Error : Total number of modifiable anistropic conductivity parameters is zero or negative !! : " << m_numAnisotropicResistivityParametersNotFixed << std::endl;
		exit(1);
	}
	m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.reserve(m_numAnisotropicResistivityParametersNotFixed);
	int icount(0);
	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		for (int iParam = 0; iParam < 6; ++iParam) {
			if (!m_fixedAnisotropicResistivityParameters[iParam][iBlk]) {
				m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.push_back(std::make_pair(iBlk, iParam));
				++icount;
			}
		}
	}
	assert(icount == m_numAnisotropicResistivityParametersNotFixed);
#ifdef _DEBUG_WRITE
	icount = 0;
	for (std::vector< std::pair<int, int> > ::const_iterator itr = m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.begin();
		itr != m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.end(); ++itr, ++icount) {
		std::cout << "i " << icount << " " << itr->first << " " << itr->second << std::endl;
	}
#endif

}


// Get number of anisotropic resistivity parameters whose values are NOT fixed
int ResistivityBlockAnisotropic::getNumberOfUnfixedResistivityParameters() const {

	return m_numAnisotropicResistivityParametersNotFixed;

}

// Get anisotropy type
int ResistivityBlockAnisotropic::getTypeOfAnisotropy(const int iblk) const {

	return m_typeOfAnisotropy[iblk];

}

// Get anisotropic resistivity parameter from resisitivity block ID
double ResistivityBlockAnisotropic::getAnisotropicResistivityParameterFromBlockID(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	const AnistropicResistivityParameters paramValues = getAnisotropicResistivityParametersFromBlockID(iblk);
	switch (iparam)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		return paramValues.rhoXX;
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		return paramValues.rhoYY;
		break;
	case ResistivityBlockAnisotropic::RHO_ZZ:
		return paramValues.rhoZZ;
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		return paramValues.strike;
		break;
	case ResistivityBlockAnisotropic::DIP:
		return paramValues.dip;
		break;
	case ResistivityBlockAnisotropic::SLANT:
		return paramValues.slant;
		break;
	default:
		OutputFiles::m_logFile << "Error : Unsupported type of anisotropy (" << iparam << ")" << std::endl;
		exit(1);
		break;
	}

}

// Get anisotropic previous resistivity parameter from resisitivity block ID
double ResistivityBlockAnisotropic::getAnisotropicResistivityParameterPreFromBlockID(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	const AnistropicResistivityParameters paramValues = getAnisotropicResistivityParametersPreFromBlockID(iblk);
	switch (iparam)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		return paramValues.rhoXX;
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		return paramValues.rhoYY;
		break;
	case ResistivityBlockAnisotropic::RHO_ZZ:
		return paramValues.rhoZZ;
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		return paramValues.strike;
		break;
	case ResistivityBlockAnisotropic::DIP:
		return paramValues.dip;
		break;
	case ResistivityBlockAnisotropic::SLANT:
		return paramValues.slant;
		break;
	default:
		OutputFiles::m_logFile << "Error : Unsupported type of anisotropy (" << iparam << ")" << std::endl;
		exit(1);
		break;
	}

}

// Get anisotropic resistivity from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParameters[iblk];


}

// Get anisotropic previous resistivity parameters from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersPreFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParametersPre[iblk];

}

// Get anisotropic fully-updated resistivity parameters from resisitivity block ID
ResistivityBlockAnisotropic::AnistropicResistivityParameters ResistivityBlockAnisotropic::getAnisotropicResistivityParametersUpdatedFullFromBlockID(const int iblk) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_anisotropicResistivityParametersUpdatedFull[iblk];

}

// Get model ID from block ID and anisotropic parameter
int ResistivityBlockAnisotropic::getModelIDFromBlockIDAndAnisotropicParameter(const int iblk, const int iparam) const {

	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);
	assert(iparam >= 0);
	assert(iparam < 6);

	return m_anisotropicResistivityParameter2modelID[iparam][iblk];

}

// Get flag specifing whether anisotrpic conductivity parameters is fixed or not
bool ResistivityBlockAnisotropic::isFixedAnisotropicResistivityParameters(const int iblk, const int iparam) const {

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());
	assert(iblk >= 0);
	assert(iblk < m_numResistivityBlockTotal);

	return m_fixedAnisotropicResistivityParameters[iparam][iblk];

}

// Get indexes of block and anistropic resistivity parameter from model index
std::pair<int, int> ResistivityBlockAnisotropic::getIndexesOfBlockAndAnisotropicResistivityParameterFromModelIndex(const int iMdl) const {

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	return m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes[iMdl];

}

// Calculate spatial difference filter for transversely isotropic conductivity
void ResistivityBlockAnisotropic::calcSpatialDifferenceFilterForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
	std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	const int paramIDs[2] = { ResistivityBlockAnisotropic::RHO_XX, ResistivityBlockAnisotropic::RHO_YY };

	if (ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		for (int i = 0; i < 2; ++i) {
			std::set< std::pair<int, int> > masterAndSlave;
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				for (int iFace = 0; iFace < 6; ++iFace) {
					if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
						continue;
					}
					if (!ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
						// First, only master elements are processed
						continue;
					}
					std::vector<int> cols;
					std::vector<double> vals;
					double rhsVectorVal(0.0);
					const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
					assert(numNeibs > 0);
					for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
						const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
						if (iElemNeib < 0) {
							continue;
						}
						const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
						if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
							continue;
						}
						const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
						masterAndSlave.insert(std::make_pair(iElem, iElemNeib));
						const double val = factor / static_cast<double>(numNeibs);
						if (isFixedAnisotropicResistivityParameters(iBlk, iParam) && isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
							continue;
						}
						else if (isFixedAnisotropicResistivityParameters(iBlk, iParam)) {
							cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
							vals.push_back(val);
							rhsVectorVal += val * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlk, iParam));
						}
						else if (isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
							cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
							vals.push_back(val);
							rhsVectorVal += val * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlkNeib, iParamNeib));
						}
						else {
							cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
							cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
							vals.push_back(val);
							vals.push_back(-val);
						}
					}
					if (cols.empty()) {
						continue;
					}
					nonZeroCols.push_back(cols);
					matValues.push_back(vals);
					rhsValues.push_back(rhsVectorVal);
				}
			}
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				for (int iFace = 0; iFace < 6; ++iFace) {
					if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
						continue;
					}
					if (ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
						// Slave elements and conforming elements are processed
						continue;
					}
					assert(ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace) == 1);
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, 0);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
						continue;
					}
					const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
					if (masterAndSlave.find(std::make_pair(iElemNeib, iElem)) != masterAndSlave.end()) {
						// This element face has already been processed
						continue;
					}
					if (iBlk >= iBlkNeib) {
						continue;
					}
					if (isFixedAnisotropicResistivityParameters(iBlk, iParam) && isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
						continue;
					}
					else if (isFixedAnisotropicResistivityParameters(iBlk, iParam)) {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						matValues.push_back(vals);
						rhsValues.push_back(factor * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlk, iParam)));
					}
					else if (isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						matValues.push_back(vals);
						rhsValues.push_back(factor * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlkNeib, iParamNeib)));
					}
					else {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						vals.push_back(-factor);
						matValues.push_back(vals);
						rhsValues.push_back(0.0);
					}
				}
			}
		}
	}
	else {
		for (int i = 0; i < 2; ++i) {
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				const int numNeibs = ptrMeshData->getNumNeighborElement();
				for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
					const int iElemNeib = ptrMeshData->getIDOfNeighborElement(iElem, iNeigh);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
						continue;
					}
					const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
					if (iBlk >= iBlkNeib) {
						continue;
					}
					if (isFixedAnisotropicResistivityParameters(iBlk, iParam) && isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
						continue;
					}
					else if (isFixedAnisotropicResistivityParameters(iBlk, iParam)) {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						matValues.push_back(vals);
						rhsValues.push_back(factor * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlk, iParam)));
					}
					else if (isFixedAnisotropicResistivityParameters(iBlkNeib, iParamNeib)) {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						matValues.push_back(vals);
						rhsValues.push_back(factor * log10(getAnisotropicResistivityParameterPreFromBlockID(iBlkNeib, iParamNeib)));
					}
					else {
						std::vector<int> cols;
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
						cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, iParamNeib));
						nonZeroCols.push_back(cols);
						std::vector<double> vals;
						vals.push_back(factor);
						vals.push_back(-factor);
						matValues.push_back(vals);
						rhsValues.push_back(0.0);
					}
				}
			}
		}
	}

	if ((AnalysisControl::getInstance())->isBottomResistivityIncluded()) {
		const double log10BottomResistivity = log10((AnalysisControl::getInstance())->getBottomResistivity());
		const double roughningFactorAtBottom = (AnalysisControl::getInstance())->getRoughningFactorAtBottom();
		for (int i = 0; i < 2; ++i) {
			const int numElemBot = ptrMeshData->getNumElemOnBoundaryPlanes(MeshData::XYPlus);
			for (int iElemBot = 0; iElemBot < numElemBot; ++iElemBot) {
				const int iElem = ptrMeshData->getElemBoundaryPlanes(MeshData::XYPlus, iElemBot);
				const int iBlk = getBlockIDFromElemID(iElem);
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				if (isFixedAnisotropicResistivityParameters(iBlk, iParam)) {
					continue;
				}
				std::vector<int> cols;
				cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, iParam));
				nonZeroCols.push_back(cols);
				std::vector<double> vals;
				vals.push_back(factor * roughningFactorAtBottom);
				matValues.push_back(vals);
				rhsValues.push_back(factor * roughningFactorAtBottom * log10BottomResistivity);
			}
		}
	}

}

// Calculate model roughness for transversely isotropic conductivity
double ResistivityBlockAnisotropic::calcModelRoughnessForTransverselyIsotropicConductivity() const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	const int paramIDs[2] = { ResistivityBlockAnisotropic::RHO_XX, ResistivityBlockAnisotropic::RHO_YY };

	double norm(0.0);
	if (ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		for (int i = 0; i < 2; ++i) {
			std::set< std::pair<int, int> > masterAndSlave;
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				const double log10Rho = log10(getAnisotropicResistivityParameterFromBlockID(iBlk, iParam));
				for (int iFace = 0; iFace < 6; ++iFace) {
					if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
						continue;
					}
					if (!ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
						// First, only master elements are processed
						continue;
					}
					const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
					assert(numNeibs > 0);
					double sumOfDiff(0.0);
					for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
						const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
						if (iElemNeib < 0) {
							continue;
						}
						const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
						if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
							continue;
						}
						const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
						const double log10RhoNeib = log10(getAnisotropicResistivityParameterFromBlockID(iBlkNeib, iParamNeib));
						sumOfDiff += (log10Rho - log10RhoNeib) / static_cast<double>(numNeibs);
						masterAndSlave.insert(std::make_pair(iElem, iElemNeib));
					}
					norm += pow(sumOfDiff, 2);
				}
			}
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				const double log10Rho = log10(getAnisotropicResistivityParameterFromBlockID(iBlk, iParam));
				for (int iFace = 0; iFace < 6; ++iFace) {
					if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
						continue;
					}
					if (ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
						// Slave elements and conforming elements are processed
						continue;
					}
					assert(ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace) == 1);
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, 0);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
						continue;
					}
					const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
					if (masterAndSlave.find(std::make_pair(iElemNeib, iElem)) != masterAndSlave.end()) {
						// This element face has already been processed
						continue;
					}
					if (iBlk >= iBlkNeib) {
						continue;
					}
					const double log10RhoNeib = log10(getAnisotropicResistivityParameterFromBlockID(iBlkNeib, iParamNeib));
					norm += pow(log10Rho - log10RhoNeib, 2);
				}
			}
 		}
	}
	else {
		for (int i = 0; i < 2; ++i) {
			for (int iElem = 0; iElem < nElem; ++iElem) {
				const int iBlk = getBlockIDFromElemID(iElem);
				if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
					continue;
				}
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				const double log10Rho = log10(getAnisotropicResistivityParameterFromBlockID(iBlk, iParam));
				const int numNeibs = ptrMeshData->getNumNeighborElement();
				for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
					const int iElemNeib = ptrMeshData->getIDOfNeighborElement(iElem, iNeigh);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY && isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::RHO_XX)) {
						continue;
					}
					const int iParamNeib = getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
					if (iBlk >= iBlkNeib) {
						continue;
					}
					const double log10RhoNeib = log10(getAnisotropicResistivityParameterFromBlockID(iBlkNeib, iParamNeib));
					norm += pow(log10Rho - log10RhoNeib, 2);
				}
			}
		}
	}

	if ((AnalysisControl::getInstance())->isBottomResistivityIncluded()) {
		const double log10BottomResistivity = log10((AnalysisControl::getInstance())->getBottomResistivity());
		const double roughningFactorAtBottom = (AnalysisControl::getInstance())->getRoughningFactorAtBottom();
		for (int i = 0; i < 2; ++i) {
			const int numElemBot = ptrMeshData->getNumElemOnBoundaryPlanes(MeshData::XYPlus);
			for (int iElemBot = 0; iElemBot < numElemBot; ++iElemBot) {
				const int iElem = ptrMeshData->getElemBoundaryPlanes(MeshData::XYPlus, iElemBot);
				const int iBlk = getBlockIDFromElemID(iElem);
				const int iParam = getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY ? ResistivityBlockAnisotropic::RHO_XX : paramIDs[i];
				if (isFixedAnisotropicResistivityParameters(iBlk, iParam)) {
					continue;
				}
				const double log10Rho = log10(getAnisotropicResistivityParameterFromBlockID(iBlk, iParam));
				norm += pow(roughningFactorAtBottom*(log10Rho - log10BottomResistivity), 2);
			}
		}
	}

	return norm;

}

// Calculate difference filter between sigma_x and sigma_y for transversely isotropic conductivity
void ResistivityBlockAnisotropic::calcDifferenceFilterBetweenComponentsForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
	std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	for (int iElem = 0; iElem < nElem; ++iElem) {
		const int iBlk = getBlockIDFromElemID(iElem);
		if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY){
			continue;
		}
		if (isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX) && isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY)) {
			continue;
		}
		else if (isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
			std::vector<int> cols;
			cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::RHO_YY));
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			vals.push_back(factor);
			matValues.push_back(vals);
			rhsValues.push_back(factor * log10(getAnisotropicResistivityParametersPreFromBlockID(iBlk).rhoXX));
		}
		else if (isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY)) {
			std::vector<int> cols;
			cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::RHO_XX));
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			vals.push_back(factor);
			matValues.push_back(vals);
			rhsValues.push_back(factor * log10(getAnisotropicResistivityParametersPreFromBlockID(iBlk).rhoYY));
		}
		else {
			std::vector<int> cols;
			cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::RHO_XX));
			cols.push_back(getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::RHO_YY));
			nonZeroCols.push_back(cols);
			std::vector<double> vals;
			vals.push_back(factor);
			vals.push_back(-factor);
			matValues.push_back(vals);
			rhsValues.push_back(0.0);
		}
	}

}

// Calculate L2 norm of the differences between sigma_x and sigma_y for transversely isotropic conductivity
double ResistivityBlockAnisotropic::calcL2NormOfDifferencesBetweenResistivityComponentsForTransverselyIsotropicConductivity() const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	double norm(0.0);
	const int nElem = ptrMeshData->getNumElemTotal();
	for (int iElem = 0; iElem < nElem; ++iElem) {
		const int iBlk = getBlockIDFromElemID(iElem);
		if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
			continue;
		}
		const AnistropicResistivityParameters paramValues = getAnisotropicResistivityParametersFromBlockID(iBlk);
		norm += pow(log10(paramValues.rhoXX) - log10(paramValues.rhoYY), 2);
	}
	return norm;

}

// Calculate cross product for transversely isotropic conductivity
void ResistivityBlockAnisotropic::calcCrossProductForTransverselyIsotropicConductivity(const double factor, std::vector<double>& values) const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	if (ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		std::set< std::pair<int, int> > masterAndSlave;
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (!ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// First, only master elements are processed
					continue;
				}
				double valueX(0.0);
				double valueY(0.0);
				double valueZ(0.0);
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				assert(numNeibs > 0);
				for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
						continue;
					}
					const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
					const double weightedFactor = factor / static_cast<double>(numNeibs);
					const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
					valueX += weightedFactor * crossProduct.X;
					valueY += weightedFactor * crossProduct.Y;
					valueZ += weightedFactor * crossProduct.Z;
					masterAndSlave.insert(std::make_pair(iElem, iElemNeib));
				}
				values.push_back(valueX);
				values.push_back(valueY);
				values.push_back(valueZ);
			}
		}
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// Slave elements and conforming elements are processed
					continue;
				}
				assert(ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace) == 1);
				const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, 0);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (masterAndSlave.find(std::make_pair(iElemNeib, iElem)) != masterAndSlave.end()) {
					// This element face has already been processed
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
				const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
				values.push_back(factor * crossProduct.X);
				values.push_back(factor * crossProduct.Y);
				values.push_back(factor * crossProduct.Z);
			}
		}
	}
	else {
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			const int numNeibs = ptrMeshData->getNumNeighborElement();
			for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement(iElem, iNeigh);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
				const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
				values.push_back(factor * crossProduct.X);
				values.push_back(factor * crossProduct.Y);
				values.push_back(factor * crossProduct.Z);
			}
		}
	}

}

// Calculate L2 norm of cross product for transversely isotropic conductivity
double ResistivityBlockAnisotropic::calcL2NormOfCrossProductForTransverselyIsotropicConductivity() const {

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	double norm(0.0);
	if (ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		std::set< std::pair<int, int> > masterAndSlave;
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (!ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// First, only master elements are processed
					continue;
				}
				double valueX(0.0);
				double valueY(0.0);
				double valueZ(0.0);
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				assert(numNeibs > 0);
				for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
						continue;
					}
					const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersFromBlockID(iBlkNeib);
					const double weight = 1.0 / static_cast<double>(numNeibs);
					const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
					valueX += weight * crossProduct.X;
					valueY += weight * crossProduct.Y;
					valueZ += weight * crossProduct.Z;
					masterAndSlave.insert(std::make_pair(iElem, iElemNeib));
				}
				norm += pow(valueX, 2);
				norm += pow(valueY, 2);
				norm += pow(valueZ, 2);
			}
		}
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// Slave elements and conforming elements are processed
					continue;
				}
				assert(ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace) == 1);
				const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, 0);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (masterAndSlave.find(std::make_pair(iElemNeib, iElem)) != masterAndSlave.end()) {
					// This element face has already been processed
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersFromBlockID(iBlkNeib);
				const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
				norm += pow(crossProduct.X, 2);
				norm += pow(crossProduct.Y, 2);
				norm += pow(crossProduct.Z, 2);
			}
		}
	}
	else {
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersFromBlockID(iBlk);
			const int numNeibs = ptrMeshData->getNumNeighborElement();
			for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement(iElem, iNeigh);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersFromBlockID(iBlkNeib);
				const CommonParameters::Vector3D crossProduct = calcCrossProductForTransverselyIsotropicConductivity(anisoParamsA, anisoParamsB);
				norm += pow(crossProduct.X, 2);
				norm += pow(crossProduct.Y, 2);
				norm += pow(crossProduct.Z, 2);
			}
		}
	}

	return norm;

}

// Calculate derivative of the cross product for transversely isotropic conductivity
void ResistivityBlockAnisotropic::calcDerivativeOfCrossProductForTransverselyIsotropicConductivity(const double factor, std::vector< std::vector<int> >& nonZeroCols,
	std::vector< std::vector<double> >& matValues, std::vector<double>& rhsValues) const{

	const MeshData* const ptrMeshData = (AnalysisControl::getInstance())->getPointerOfMeshData();
	const int nElem = ptrMeshData->getNumElemTotal();
	if (ptrMeshData->getMeshType() == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		std::set< std::pair<int, int> > masterAndSlave;
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (!ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// First, only master elements are processed
					continue;
				}
				std::vector<int> cols[3];
				std::vector<double> vals[3];
				const int numNeibs = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace);
				assert(numNeibs > 0);
				for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
					const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, iNeigh);
					if (iElemNeib < 0) {
						continue;
					}
					const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
					if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
						continue;
					}
					const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
					calcDerivativeOfCrossProductForTransverselyIsotropicConductivityAux(iBlk, iBlkNeib, anisoParamsA, anisoParamsB, factor / static_cast<double>(numNeibs), cols, vals);
					masterAndSlave.insert(std::make_pair(iElem, iElemNeib));
				}
				for (int i = 0; i < 3; ++i){
					nonZeroCols.push_back(cols[i]);
					matValues.push_back(vals[i]);
					rhsValues.push_back(0.0);
				}
			}
		}
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			for (int iFace = 0; iFace < 6; ++iFace) {
				if (ptrMeshDataNonConformingHexaElement->isOuterBoundary(iElem, iFace)) {
					continue;
				}
				if (ptrMeshDataNonConformingHexaElement->faceSlaveElements(iElem, iFace)) {
					// Slave elements and conforming elements are processed
					continue;
				}
				assert(ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, iFace) == 1);
				const int iElemNeib = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, iFace, 0);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (masterAndSlave.find(std::make_pair(iElemNeib, iElem)) != masterAndSlave.end()) {
					// This element face has already been processed
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
				std::vector<int> cols[3];
				std::vector<double> vals[3];
				calcDerivativeOfCrossProductForTransverselyIsotropicConductivityAux(iBlk, iBlkNeib, anisoParamsA, anisoParamsB, factor, cols, vals);
				for (int i = 0; i < 3; ++i) {
					nonZeroCols.push_back(cols[i]);
					matValues.push_back(vals[i]);
					rhsValues.push_back(0.0);
				}
			}
		}
	}
	else {
		for (int iElem = 0; iElem < nElem; ++iElem) {
			const int iBlk = getBlockIDFromElemID(iElem);
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				continue;
			}
			const AnistropicResistivityParameters anisoParamsA = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
			const int numNeibs = ptrMeshData->getNumNeighborElement();
			for (int iNeigh = 0; iNeigh < numNeibs; ++iNeigh) {
				const int iElemNeib = ptrMeshData->getIDOfNeighborElement(iElem, iNeigh);
				if (iElemNeib < 0) {
					continue;
				}
				const int iBlkNeib = getBlockIDFromElemID(iElemNeib);
				if (getTypeOfAnisotropy(iBlkNeib) == ResistivityBlockAnisotropic::ISOTROPY) {
					continue;
				}
				if (iBlk >= iBlkNeib) {
					continue;
				}
				const AnistropicResistivityParameters anisoParamsB = getAnisotropicResistivityParametersPreFromBlockID(iBlkNeib);
				std::vector<int> cols[3];
				std::vector<double> vals[3];
				calcDerivativeOfCrossProductForTransverselyIsotropicConductivityAux(iBlk, iBlkNeib, anisoParamsA, anisoParamsB, factor, cols, vals);
				for (int i = 0; i < 3; ++i) {
					nonZeroCols.push_back(cols[i]);
					matValues.push_back(vals[i]);
					rhsValues.push_back(0.0);
				}
			}
		}
	}

}

// Calculate cross product for transversely isotropic conductivity
CommonParameters::Vector3D ResistivityBlockAnisotropic::calcCrossProductForTransverselyIsotropicConductivity(const AnistropicResistivityParameters& anisoParamsA, const AnistropicResistivityParameters& anisoParamsB) const {

	const CommonParameters::Vector3D crossProduct = { 
		 cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * sin(anisoParamsB.dip) - sin(anisoParamsA.dip) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike),
		-sin(anisoParamsA.dip) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike) + cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * sin(anisoParamsB.dip),
		-cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike) + cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike)
	};
	return crossProduct;

}

// Auxiliary function to calculate derivative of the cross product for transversely isotropic conductivity
void ResistivityBlockAnisotropic::calcDerivativeOfCrossProductForTransverselyIsotropicConductivityAux(const int iBlk, const int iBlkNeib,
	const AnistropicResistivityParameters& anisoParamsA, const AnistropicResistivityParameters& anisoParamsB, const double factor, std::vector<int> cols[3], std::vector<double> vals[3]) const {

	if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP))
	{
		const int col = getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::DIP);
		cols[0].push_back(col);
		cols[1].push_back(col);
		cols[2].push_back(col);
		const CommonParameters::Vector3D derivativeOfCrossProduct = {
			-sin(anisoParamsA.dip) * cos(anisoParamsA.strike) * sin(anisoParamsB.dip) - cos(anisoParamsA.dip) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike),
			-cos(anisoParamsA.dip) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike) - sin(anisoParamsA.dip) * sin(anisoParamsA.strike) * sin(anisoParamsB.dip),
			 sin(anisoParamsA.dip) * sin(anisoParamsA.strike) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike) - sin(anisoParamsA.dip) * cos(anisoParamsA.strike) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike)
		};
		vals[0].push_back(factor * derivativeOfCrossProduct.X);
		vals[1].push_back(factor * derivativeOfCrossProduct.Y);
		vals[2].push_back(factor * derivativeOfCrossProduct.Z);
	}
	if (!isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::DIP))
	{
		const int col = getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, ResistivityBlockAnisotropic::DIP);
		cols[0].push_back(col);
		cols[1].push_back(col);
		cols[2].push_back(col);
		const CommonParameters::Vector3D derivativeOfCrossProduct = {
			cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * cos(anisoParamsB.dip) + sin(anisoParamsA.dip) * sin(anisoParamsB.dip) * cos(anisoParamsB.strike),
			sin(anisoParamsA.dip) * sin(anisoParamsB.dip) * sin(anisoParamsB.strike) + cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * cos(anisoParamsB.dip),
			cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * sin(anisoParamsB.dip) * cos(anisoParamsB.strike) - cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * sin(anisoParamsB.dip) * sin(anisoParamsB.strike)
		};
		vals[0].push_back(factor * derivativeOfCrossProduct.X);
		vals[1].push_back(factor * derivativeOfCrossProduct.Y);
		vals[2].push_back(factor * derivativeOfCrossProduct.Z);
	}
	if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::STRIKE))
	{
		const int col = getModelIDFromBlockIDAndAnisotropicParameter(iBlk, ResistivityBlockAnisotropic::STRIKE);
		cols[0].push_back(col);
		cols[1].push_back(col);
		cols[2].push_back(col);
		const CommonParameters::Vector3D derivativeOfCrossProduct = {
			-cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * sin(anisoParamsB.dip),
			 cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * sin(anisoParamsB.dip),
		    -cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * cos(anisoParamsB.dip)* cos(anisoParamsB.strike) - cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike)
		};
		vals[0].push_back(factor * derivativeOfCrossProduct.X);
		vals[1].push_back(factor * derivativeOfCrossProduct.Y);
		vals[2].push_back(factor * derivativeOfCrossProduct.Z);
	}
	if (!isFixedAnisotropicResistivityParameters(iBlkNeib, ResistivityBlockAnisotropic::STRIKE))
	{
		const int col = getModelIDFromBlockIDAndAnisotropicParameter(iBlkNeib, ResistivityBlockAnisotropic::STRIKE);
		cols[0].push_back(col);
		cols[1].push_back(col);
		cols[2].push_back(col);
		const CommonParameters::Vector3D derivativeOfCrossProduct = {
			 sin(anisoParamsA.dip) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike),
			-sin(anisoParamsA.dip) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike),
			 cos(anisoParamsA.dip) * sin(anisoParamsA.strike) * cos(anisoParamsB.dip) * sin(anisoParamsB.strike) + cos(anisoParamsA.dip) * cos(anisoParamsA.strike) * cos(anisoParamsB.dip) * cos(anisoParamsB.strike)
		};
		vals[0].push_back(factor * derivativeOfCrossProduct.X);
		vals[1].push_back(factor * derivativeOfCrossProduct.Y);
		vals[2].push_back(factor * derivativeOfCrossProduct.Z);
	}

}

// Confirm that there is no element having general anisotropy
void ResistivityBlockAnisotropic::confirmNoElementHavingGeneralAnisotropy() const {

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::GENERAL_ANISOTROPY) {
			OutputFiles::m_logFile << "Error : Resistivity block (" << iBlk << ") have general anisotropy !!" << std::endl;
			exit(1);
		}
	}

}

// Calculate array of fully updated anisotropic resistivity parameters obtained by inversion
void ResistivityBlockAnisotropic::calcUnfixedAnisotropicResistivityParametersUpdatedFull(const double* const updatedModel) {

	int index(0);
	for (std::vector< std::pair<int, int> >::const_iterator itr = m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.begin();
		itr != m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.end(); ++itr, ++index) {
		const int iBlk = itr->first;
		const AnistropicResistivityParameters anisoParamsPre = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
		const int iParam = itr->second;
		switch (iParam)
		{
		case ResistivityBlockAnisotropic::RHO_XX:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoXX = pow(10.0, log10(anisoParamsPre.rhoXX) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::RHO_YY:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoYY = pow(10.0, log10(anisoParamsPre.rhoYY) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::RHO_ZZ:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].rhoZZ = pow(10.0, log10(anisoParamsPre.rhoZZ) + updatedModel[index]);
			break;
		case ResistivityBlockAnisotropic::STRIKE:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].strike = anisoParamsPre.strike + updatedModel[index];
			break;
		case ResistivityBlockAnisotropic::DIP:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].dip = anisoParamsPre.dip + updatedModel[index];
			break;
		case ResistivityBlockAnisotropic::SLANT:
			m_anisotropicResistivityParametersUpdatedFull[iBlk].slant = anisoParamsPre.slant + updatedModel[index];
			break;
		default:
			OutputFiles::m_logFile << "Error : Unsupported type of anisotropy (" << iParam << ")" << std::endl;
			exit(1);
			break;
		}
	}

}

// Change anisotropic resistivity parameters
void ResistivityBlockAnisotropic::updateResistivityValues() {

	const double stepLengthDampingFactor = (AnalysisControl::getInstance())->getStepLengthDampingFactorCur();
	const double upperLimitOfAngleUpdates = (AnalysisControl::getInstance())->getUpperLimitOfAbsAngleUpdatesFoAnisotropy() * CommonParameters::deg2rad;

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		const AnistropicResistivityParameters anisoParamsUpdatedFull = getAnisotropicResistivityParametersUpdatedFullFromBlockID(iBlk);
		const AnistropicResistivityParameters anisoParamsPre = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
		AnistropicResistivityParameters anisoParams = anisoParamsPre;
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
			anisoParams.rhoXX = pow(10.0, stepLengthDampingFactor * log10(anisoParamsUpdatedFull.rhoXX) + (1.0 - stepLengthDampingFactor) * log10(anisoParamsPre.rhoXX));
			if (anisoParams.rhoXX < m_resistivityValuesMin[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_XX ( " << anisoParams.rhoXX << " [Ohm-m] ) of block " << iBlk << " is lower than the minimum value ( " << m_resistivityValuesMin[iBlk] << " [Ohm-m] ). Its resistivity is set to be the minimum value." << std::endl;
				anisoParams.rhoXX = m_resistivityValuesMin[iBlk];
			}
			else if (anisoParams.rhoXX > m_resistivityValuesMax[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_XX ( " << anisoParams.rhoXX << " [Ohm-m] ) of block " << iBlk << " is higher the maximum value ( " << m_resistivityValuesMax[iBlk] << " [Ohm-m] ). Its resistivity is set to be the maximum value." << std::endl;
				anisoParams.rhoXX = m_resistivityValuesMax[iBlk];
			}
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				anisoParams.rhoYY = anisoParams.rhoXX;
				anisoParams.rhoZZ = anisoParams.rhoXX;
			}
			else if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::TRANSVERSE_ISOTROPY) {
				anisoParams.rhoZZ = anisoParams.rhoXX;
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY)) {
			anisoParams.rhoYY = pow(10.0, stepLengthDampingFactor * log10(anisoParamsUpdatedFull.rhoYY) + (1.0 - stepLengthDampingFactor) * log10(anisoParamsPre.rhoYY));
			if (anisoParams.rhoYY < m_resistivityValuesMin[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_YY ( " << anisoParams.rhoYY << " [Ohm-m] ) of block " << iBlk << " is lower than the minimum value ( " << m_resistivityValuesMin[iBlk] << " [Ohm-m] ). Its resistivity is set to be the minimum value." << std::endl;
				anisoParams.rhoYY = m_resistivityValuesMin[iBlk];
			}
			else if (anisoParams.rhoYY > m_resistivityValuesMax[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_YY ( " << anisoParams.rhoYY << " [Ohm-m] ) of block " << iBlk << " is higher the maximum value ( " << m_resistivityValuesMax[iBlk] << " [Ohm-m] ). Its resistivity is set to be the maximum value." << std::endl;
				anisoParams.rhoYY = m_resistivityValuesMax[iBlk];
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_ZZ)) {
			anisoParams.rhoZZ = pow(10.0, stepLengthDampingFactor * log10(anisoParamsUpdatedFull.rhoZZ) + (1.0 - stepLengthDampingFactor) * log10(anisoParamsPre.rhoZZ));
			if (anisoParams.rhoZZ < m_resistivityValuesMin[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_YY ( " << anisoParams.rhoZZ << " [Ohm-m] ) of block " << iBlk << " is lower than the minimum value ( " << m_resistivityValuesMin[iBlk] << " [Ohm-m] ). Its resistivity is set to be the minimum value." << std::endl;
				anisoParams.rhoZZ = m_resistivityValuesMin[iBlk];
			}
			else if (anisoParams.rhoZZ > m_resistivityValuesMax[iBlk]) {
				OutputFiles::m_logFile << "Warning : Updated rho_YY ( " << anisoParams.rhoZZ << " [Ohm-m] ) of block " << iBlk << " is higher the maximum value ( " << m_resistivityValuesMax[iBlk] << " [Ohm-m] ). Its resistivity is set to be the maximum value." << std::endl;
				anisoParams.rhoZZ = m_resistivityValuesMax[iBlk];
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::STRIKE)) {
			double update = stepLengthDampingFactor * (anisoParamsUpdatedFull.strike - anisoParamsPre.strike);
			if (update > upperLimitOfAngleUpdates) {
				update = upperLimitOfAngleUpdates;
			}else if (update < -upperLimitOfAngleUpdates) {
				update = -upperLimitOfAngleUpdates;
			}
			anisoParams.strike = anisoParamsPre.strike + update;
			while (anisoParams.strike > CommonParameters::PI) {
				anisoParams.strike -= 2.0 * CommonParameters::PI;
			}
			while(anisoParams.strike < -CommonParameters::PI) {
				anisoParams.strike += 2.0 * CommonParameters::PI;
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP)) {
			double update = stepLengthDampingFactor * (anisoParamsUpdatedFull.dip - anisoParamsPre.dip);
			if (update > upperLimitOfAngleUpdates) {
				update = upperLimitOfAngleUpdates;
			}else if (update < -upperLimitOfAngleUpdates) {
				update = -upperLimitOfAngleUpdates;
			}
			anisoParams.dip = anisoParamsPre.dip + update;
			while (anisoParams.dip > 0.5 * CommonParameters::PI) {
				anisoParams.dip -= CommonParameters::PI;
			}
			while(anisoParams.dip < -0.5 * CommonParameters::PI) {
				anisoParams.dip += CommonParameters::PI;
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::SLANT)) {
			double update = stepLengthDampingFactor * (anisoParamsUpdatedFull.slant - anisoParamsPre.slant);
			if (update > upperLimitOfAngleUpdates) {
				update = upperLimitOfAngleUpdates;
			}else if (update < -upperLimitOfAngleUpdates) {
				update = -upperLimitOfAngleUpdates;
			}
			anisoParams.slant = anisoParamsPre.slant + update;
			while(anisoParams.slant > CommonParameters::PI) {
				anisoParams.slant -= 2.0 * CommonParameters::PI;
			}
			while(anisoParams.slant < -CommonParameters::PI) {
				anisoParams.slant += 2.0 * CommonParameters::PI;
			}
		}
		m_anisotropicResistivityParameters[iBlk] = anisoParams;
	}

}

// Output anisotropy types to VTK file
void ResistivityBlockAnisotropic::outputAnisotropyTypesToVTK() const {

	if (!OutputFiles::m_vtkFile.is_open()) {
		return;
	}

	const AnalysisControl* pAnalysisControl = AnalysisControl::getInstance();
	if (pAnalysisControl->doesOutputToVTK(AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK)) {
		const int nElem = ((AnalysisControl::getInstance())->getPointerOfMeshData())->getNumElemTotal();
		OutputFiles::m_vtkFile << "SCALARS AnisotropyTypes int" << std::endl;
		OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
		for (int iElem = 0; iElem < nElem; ++iElem) {
			OutputFiles::m_vtkFile << getTypeOfAnisotropy(getBlockIDFromElemID(iElem)) << std::endl;
		}
	}
}

// Output anisotropy types to binary file
void ResistivityBlockAnisotropic::outputAnisotropyTypesToBinary() const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}
	std::ofstream fout;
	fout.open("AnisotropyTypes", std::ios::out | std::ios::binary | std::ios::trunc);
	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Anisotropy types";
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
		float dbuf = static_cast<float>(getTypeOfAnisotropy(getBlockIDFromElemID(iElem)));
		fout.write((char*)&dbuf, sizeof(float));
	}
	fout.close();

}

// Calculate anisotropic conductivity tensor
void ResistivityBlockAnisotropic::calcAnisotropicConductivityTensor( const int blockID, double conductivityTensor[3][3] ) const{

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	// Zero clear
	for( int row = 0; row < 3; ++row ){
		for( int col = 0; col < 3; ++col ){
			conductivityTensor[row][col] = 0.0;
		}
	}

	const AnistropicResistivityParameters anisotropicResistivityParam = getAnisotropicResistivityParametersFromBlockID(blockID);
	conductivityTensor[0][0] = 1.0 / anisotropicResistivityParam.rhoXX;
	conductivityTensor[1][1] = 1.0 / anisotropicResistivityParam.rhoYY;
	conductivityTensor[2][2] = 1.0 / anisotropicResistivityParam.rhoZZ;

	if( m_typeOfAnisotropy[blockID] == ISOTROPY ) {
		// The conducitivity of the block is isotropic
		return;
	}

	// Rotation around z'-axis
	const double slant = anisotropicResistivityParam.slant;
	const double rotationTensorZ2[3][3] = {
		{  cos(slant), -sin(slant),  0.0 },
		{  sin(slant),  cos(slant),  0.0 },
		{  0.0,         0.0,         1.0 } };
	rotateTensor(conductivityTensor, rotationTensorZ2);

	// Rotation around x'-axis
	const double dip = anisotropicResistivityParam.dip;
	const double rotationTensorY[3][3] = {
		{  1.0, 0.0,       0.0      },
		{  0.0, cos(dip), -sin(dip) },
		{  0.0, sin(dip),  cos(dip) } };
	rotateTensor( conductivityTensor, rotationTensorY);

	// Rotation around z-axis
	const double strike = anisotropicResistivityParam.strike;
	const double rotationTensorZ[3][3] = {
		{  cos(strike), -sin(strike),  0.0 },
		{  sin(strike),  cos(strike),  0.0 },
		{  0.0,          0.0,          1.0 } };
	rotateTensor(conductivityTensor, rotationTensorZ);

}

// Calculate derivative of anisotropic conductivity tensor
void ResistivityBlockAnisotropic::calcDerivativeOfAnisotropicConductivityTensor(const int blockID, const int paramID, double derivative[3][3]) const {

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	// Zero clear
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col) {
			derivative[row][col] = 0.0;
		}
	}

	const AnistropicResistivityParameters anisotropicResistivityParam = getAnisotropicResistivityParametersFromBlockID(blockID);
	switch (getTypeOfAnisotropy(blockID))
	{
	case ResistivityBlockAnisotropic::ISOTROPY:
		derivative[0][0] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		derivative[1][1] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		derivative[2][2] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		break;
	case ResistivityBlockAnisotropic::TRANSVERSE_ISOTROPY:
		calcDerivativeOfTransverseIsotropicConductivityTensor(paramID, anisotropicResistivityParam, derivative);
		break;
	case ResistivityBlockAnisotropic::GENERAL_ANISOTROPY:
		calcDerivativeOfGeneralAnisotropicConductivityTensor(paramID, anisotropicResistivityParam, derivative);
		break;
	default:
		OutputFiles::m_logFile << "Error : Unsupported type of anisotropy (" << getTypeOfAnisotropy(blockID) << ")" << std::endl;
		exit(1);
		break;
	}

}

// Calculate derivative of transverse isotropic conductivity tensor
void ResistivityBlockAnisotropic::calcDerivativeOfTransverseIsotropicConductivityTensor(const int paramID, const AnistropicResistivityParameters& anisotropicResistivityParam,
	double derivative[3][3]) const {

	// Rotation around z-axis
	const double strike = anisotropicResistivityParam.strike;
	const double rotationTensorZ[3][3] = {
		{  cos(strike), -sin(strike),  0.0 },
		{  sin(strike),  cos(strike),  0.0 },
		{  0.0,          0.0,          1.0 } };
	const double derivRotationTensorZ[3][3] = {
		{ -sin(strike), -cos(strike),  0.0 },
		{  cos(strike), -sin(strike),  0.0 },
		{  0.0,          0.0,          0.0} };

	// Rotation around x'-axis
	const double dip = anisotropicResistivityParam.dip;
	const double rotationTensorX[3][3] = {
		{  1.0,  0.0,       0.0      },
		{  0.0,  cos(dip), -sin(dip) },
		{  0.0,  sin(dip),  cos(dip) } };
	const double derivRotationTensorX[3][3] = {
		{  0.0,  0.0,       0.0      },
		{  0.0, -sin(dip), -cos(dip) },
		{  0.0,  cos(dip), -sin(dip) } };

	double conductivityTensor1[3][3] = {
		{ 1.0 / anisotropicResistivityParam.rhoXX, 0.0, 0.0 },
		{ 0.0, 1.0 / anisotropicResistivityParam.rhoYY, 0.0 },
		{ 0.0, 0.0, 1.0 / anisotropicResistivityParam.rhoXX } };
	double conductivityTensor2[3][3] = {
		{ 1.0 / anisotropicResistivityParam.rhoXX, 0.0, 0.0 },
		{ 0.0, 1.0 / anisotropicResistivityParam.rhoYY, 0.0 },
		{ 0.0, 0.0, 1.0 / anisotropicResistivityParam.rhoXX } };

	switch (paramID)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		// Go through
	case ResistivityBlockAnisotropic::RHO_ZZ:
		derivative[0][0] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		derivative[2][2] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		// Rotation around x'-axis
		rotateTensor(derivative, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(derivative, rotationTensorZ);
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		derivative[1][1] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoYY;
		// Rotation around x'-axis
		rotateTensor(derivative, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(derivative, rotationTensorZ);
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		rotateTensor(conductivityTensor1, rotationTensorX);
		calcDerivativeOfAnisotropicConductivityTensorAux(derivRotationTensorZ, conductivityTensor1, rotationTensorZ);
		rotateTensor(conductivityTensor2, rotationTensorX);
		calcDerivativeOfAnisotropicConductivityTensorAux(rotationTensorZ, conductivityTensor2, derivRotationTensorZ);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				derivative[row][col] = conductivityTensor1[row][col] + conductivityTensor2[row][col];
			}
		}
		break;
	case ResistivityBlockAnisotropic::DIP:
		calcDerivativeOfAnisotropicConductivityTensorAux(derivRotationTensorX, conductivityTensor1, rotationTensorX);
		rotateTensor(conductivityTensor1, rotationTensorZ);
		calcDerivativeOfAnisotropicConductivityTensorAux(rotationTensorX, conductivityTensor2, derivRotationTensorX);
		rotateTensor(conductivityTensor2, rotationTensorZ);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				derivative[row][col] = conductivityTensor1[row][col] + conductivityTensor2[row][col];
			}
		}
		break;
	default:
		OutputFiles::m_logFile << "Error : Unsupported parameter for transverse isotropy (" << paramID << ")" << std::endl;
		exit(1);
		break;
	}

}

// Calculate derivative of general anisotropic conductivity tensor
void ResistivityBlockAnisotropic::calcDerivativeOfGeneralAnisotropicConductivityTensor(const int paramID, const AnistropicResistivityParameters& anisotropicResistivityParam, double derivative[3][3]) const {

	// Rotation around z-axis
	const double strike = anisotropicResistivityParam.strike;
	const double rotationTensorZ[3][3] = {
		{  cos(strike), -sin(strike),  0.0 },
		{  sin(strike),  cos(strike),  0.0 },
		{  0.0,          0.0,          1.0 } };
	const double derivRotationTensorZ[3][3] = {
		{ -sin(strike), -cos(strike),  0.0 },
		{  cos(strike), -sin(strike),  0.0 },
		{  0.0,          0.0,          0.0} };

	// Rotation around x'-axis
	const double dip = anisotropicResistivityParam.dip;
	const double rotationTensorX[3][3] = {
		{  1.0,  0.0,       0.0      },
		{  0.0,  cos(dip), -sin(dip) },
		{  0.0,  sin(dip),  cos(dip) } };
	const double derivRotationTensorX[3][3] = {
		{  0.0,  0.0,       0.0      },
		{  0.0, -sin(dip), -cos(dip) },
		{  0.0,  cos(dip), -sin(dip) } };

	// Rotation around z'-axis
	const double slant = anisotropicResistivityParam.slant;
	const double rotationTensorZ2[3][3] = {
		{  cos(slant), -sin(slant),  0.0 },
		{  sin(slant),  cos(slant),  0.0 },
		{  0.0,         0.0,         1.0 } };
	const double derivRotationTensorZ2[3][3] = {
		{ -sin(slant), -cos(slant),  0.0 },
		{  cos(slant), -sin(slant),  0.0 },
		{  0.0,         0.0,         0.0} };

	double conductivityTensor1[3][3] = {
		{ 1.0 / anisotropicResistivityParam.rhoXX, 0.0, 0.0 },
		{ 0.0, 1.0 / anisotropicResistivityParam.rhoYY, 0.0 },
		{ 0.0, 0.0, 1.0 / anisotropicResistivityParam.rhoZZ } };
	double conductivityTensor2[3][3] = {
		{ 1.0 / anisotropicResistivityParam.rhoXX, 0.0, 0.0 },
		{ 0.0, 1.0 / anisotropicResistivityParam.rhoYY, 0.0 },
		{ 0.0, 0.0, 1.0 / anisotropicResistivityParam.rhoZZ } };

	switch (paramID)
	{
	case ResistivityBlockAnisotropic::RHO_XX:
		derivative[0][0] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoXX;
		// Rotation around z'-axis
		rotateTensor(derivative, rotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(derivative, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(derivative, rotationTensorZ);
		break;
	case ResistivityBlockAnisotropic::RHO_YY:
		derivative[1][1] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoYY;
		// Rotation around z'-axis
		rotateTensor(derivative, rotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(derivative, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(derivative, rotationTensorZ);
		break;
	case ResistivityBlockAnisotropic::RHO_ZZ:
		derivative[2][2] = -CommonParameters::ln10 / anisotropicResistivityParam.rhoZZ;
		// Rotation around z'-axis
		rotateTensor(derivative, rotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(derivative, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(derivative, rotationTensorZ);
		break;
	case ResistivityBlockAnisotropic::STRIKE:
		// Rotation around z'-axis
		rotateTensor(conductivityTensor1, rotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(conductivityTensor1, rotationTensorX);
		// Rotation around z'-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(derivRotationTensorZ, conductivityTensor1, rotationTensorZ);
		// Rotation around z'-axis
		rotateTensor(conductivityTensor2, rotationTensorZ2);
		// Rotation around y'-axis
		rotateTensor(conductivityTensor2, rotationTensorX);
		// Rotation around z-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(rotationTensorZ, conductivityTensor2, derivRotationTensorZ);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				derivative[row][col] = conductivityTensor1[row][col] + conductivityTensor2[row][col];
			}
		}
		break;
	case ResistivityBlockAnisotropic::DIP:
		// Rotation around z'-axis
		rotateTensor(conductivityTensor1, rotationTensorZ2);
		// Rotation around x'-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(derivRotationTensorX, conductivityTensor1, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(conductivityTensor1, rotationTensorZ);
		// Rotation around z'-axis
		rotateTensor(conductivityTensor2, rotationTensorZ2);
		// Rotation around x'-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(rotationTensorX, conductivityTensor2, derivRotationTensorX);
		// Rotation around z-axis
		rotateTensor(conductivityTensor2, rotationTensorZ);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				derivative[row][col] = conductivityTensor1[row][col] + conductivityTensor2[row][col];
			}
		}
		break;
	case ResistivityBlockAnisotropic::SLANT:
		// Rotation around z'-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(derivRotationTensorZ2, conductivityTensor1, rotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(conductivityTensor1, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(conductivityTensor1, rotationTensorZ);
		// Rotation around z'-axis
		calcDerivativeOfAnisotropicConductivityTensorAux(rotationTensorZ2, conductivityTensor2, derivRotationTensorZ2);
		// Rotation around x'-axis
		rotateTensor(conductivityTensor2, rotationTensorX);
		// Rotation around z-axis
		rotateTensor(conductivityTensor2, rotationTensorZ);
		for (int row = 0; row < 3; ++row) {
			for (int col = 0; col < 3; ++col) {
				derivative[row][col] = conductivityTensor1[row][col] + conductivityTensor2[row][col];
			}
		}
		break;
	default:
		OutputFiles::m_logFile << "Error : Unsupported parameter for transverse isotropy (" << paramID << ")" << std::endl;
		exit(1);
		break;
	}

}

// Auxiliary function for calculating derivative of anisotropic conductivity tensor
void ResistivityBlockAnisotropic::calcDerivativeOfAnisotropicConductivityTensorAux(const double lhsMatrix[3][3], double matrix[3][3], const double transposeOfRhsMatrix[3][3]) const{

	double tempTensor[3][3];
	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col) {
			tempTensor[row][col] = 0.0;
			for (int i = 0; i < 3; ++i) {
				tempTensor[row][col] += lhsMatrix[row][i] * matrix[i][col];
			}
		}
	}

	for (int row = 0; row < 3; ++row) {
		for (int col = 0; col < 3; ++col) {
			matrix[row][col] = 0.0;
			for (int i = 0; i < 3; ++i) {
				matrix[row][col] += tempTensor[row][i] * transposeOfRhsMatrix[col][i];
			}
		}
	}

}

// Output resistivity values to VTK file
void ResistivityBlockAnisotropic::outputResistivityValuesToVTK() const {

	if (!OutputFiles::m_vtkFile.is_open()) {
		return;
	}
	if (!(AnalysisControl::getInstance())->doesOutputToVTK(AnalysisControl::OUTPUT_RESISTIVITY_VALUES_TO_VTK)) {
		return;
	}

	const int nElem = ((AnalysisControl::getInstance())->getPointerOfMeshData())->getNumElemTotal();
	OutputFiles::m_vtkFile << "SCALARS Rho_XX[Ohm-m] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_XX) << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Rho_YY[Ohm-m] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_YY) << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Rho_ZZ[Ohm-m] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_ZZ) << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Strike[deg.] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::STRIKE) * CommonParameters::rad2deg << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Dip[deg.] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::DIP) * CommonParameters::rad2deg << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Slant[deg.] float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		OutputFiles::m_vtkFile << getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::SLANT) * CommonParameters::rad2deg << std::endl;
	}
	OutputFiles::m_vtkFile << "SCALARS Anisotropy float" << std::endl;
	OutputFiles::m_vtkFile << "LOOKUP_TABLE default" << std::endl;
	for (int iElem = 0; iElem < nElem; ++iElem) {
		const double log10RhoXX = log10(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_XX));
		const double log10RhoYY = log10(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_YY));
		const double anisotropy = fabs(log10RhoXX - log10RhoYY);
		OutputFiles::m_vtkFile << anisotropy << std::endl;
	}

}

// Output resistivity values to binary file
void ResistivityBlockAnisotropic::outputResistivityValuesToBinary(const int iterNum) const {

	outputRhoXXOfAnisotropicResistivityTensorToBinary(iterNum);
	outputRhoYYOfAnisotropicResistivityTensorToBinary(iterNum);
	outputRhoZZOfAnisotropicResistivityTensorToBinary(iterNum);
	outputStrikeOfAnisotropicResistivityTensorToBinary(iterNum);
	outputDipOfAnisotropicResistivityTensorToBinary(iterNum);
	outputSlantOfAnisotropicResistivityTensorToBinary(iterNum);
	outputfAnisotropyIndicatorToBinary(iterNum);

}

// Output xx-component of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputRhoXXOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "RhoXX.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "RhoXX[Ohm-m]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_XX));
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output yy-component of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputRhoYYOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "RhoYY.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "RhoYY[Ohm-m]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_YY));
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output zz-component of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputRhoZZOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "RhoZZ.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "RhoZZ[Ohm-m]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_ZZ));
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output strike of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputStrikeOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Strike.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Strike[Deg.]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::STRIKE) * CommonParameters::rad2deg);
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();
}

// Output dip of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputDipOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Dip.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Dip[Deg.]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::DIP) * CommonParameters::rad2deg);
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output slant of the anisotropic resistivity tensor to binary file
void ResistivityBlockAnisotropic::outputSlantOfAnisotropicResistivityTensorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Slant.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Slant[Deg.]";
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
		float dbuf = static_cast<float>(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::SLANT) * CommonParameters::rad2deg);
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output abs(log10(rho_xx)-log10(rho_yy)) to binary file
void ResistivityBlockAnisotropic::outputfAnisotropyIndicatorToBinary(const int iterNum) const {

	const AnalysisControl* const ptrAnalysisControl = AnalysisControl::getInstance();
	if (ptrAnalysisControl->getMyPE() != 0) {// Only PE 0 output sensitivity
		return;
	}

	std::ostringstream oss;
	oss << "Anisotropy.iter" << iterNum;
	std::ofstream fout;
	fout.open(oss.str().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

	char line[80];
	std::ostringstream ossTitle;
	ossTitle << "Anisotropy";
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
		const double log10RhoXX = log10(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_XX));
		const double log10RhoYY = log10(getAnisotropicResistivityParameterFromBlockID(getBlockIDFromElemID(iElem), ResistivityBlockAnisotropic::RHO_YY));
		const float dbuf = static_cast<float>(fabs(log10RhoXX - log10RhoYY));
		fout.write((char*)&dbuf, sizeof(float));
	}

	fout.close();

}

// Output data of resisitivity block model to file (anisotropic case)
void ResistivityBlockAnisotropic::outputResistivityValues(const int iterNum, FILE* fp) const {

	for (int iBlk = 0; iBlk < m_numResistivityBlockTotal; ++iBlk) {
		const int typeOfAnisotropy = getTypeOfAnisotropy(iBlk);
		switch (typeOfAnisotropy) {
		case ResistivityBlockAnisotropic::ISOTROPY:
			fprintf(fp, "%10d%10d%5s%15e%15e%15e%10d\n", iBlk, typeOfAnisotropy, "     ",
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_XX),
				m_resistivityValuesMin[iBlk], m_resistivityValuesMax[iBlk],
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX) ? 1 : 0);
			break;
		case ResistivityBlockAnisotropic::TRANSVERSE_ISOTROPY:
			fprintf(fp, "%10d%10d%5s%15e%15e%15e%15e%15e%15e%10d%10d%10d%10d\n", iBlk, typeOfAnisotropy, "     ",
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_XX),
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_YY),
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::STRIKE) * CommonParameters::rad2deg,
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::DIP) * CommonParameters::rad2deg,
				m_resistivityValuesMin[iBlk], m_resistivityValuesMax[iBlk],
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::STRIKE) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP) ? 1 : 0);
			break;
		case ResistivityBlockAnisotropic::GENERAL_ANISOTROPY:
			fprintf(fp, "%10d%10d%5s%15e%15e%15e%15e%15e%15e%15e%15e%10d%10d%10d%10d%10d%10d\n", iBlk, typeOfAnisotropy, "     ",
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_XX),
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_YY),
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::RHO_ZZ),
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::STRIKE) * CommonParameters::rad2deg,
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::DIP) * CommonParameters::rad2deg,
				getAnisotropicResistivityParameterFromBlockID(iBlk, ResistivityBlockAnisotropic::SLANT) * CommonParameters::rad2deg,
				m_resistivityValuesMin[iBlk], m_resistivityValuesMax[iBlk],
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_ZZ) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::STRIKE) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP) ? 1 : 0,
				isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::SLANT) ? 1 : 0);
			break;
		}
	}

}

// Copy current anisotropic resistivity parameters to current ones
void ResistivityBlockAnisotropic::copyUnfixedResistivityParametersCurToVector(double* vector) const {

	int index(0);
	for (std::vector< std::pair<int, int> >::const_iterator itr = m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.begin();
		itr != m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.end(); ++itr, ++index) {
		const int iBlk = itr->first;
		const AnistropicResistivityParameters anisoParams = getAnisotropicResistivityParametersFromBlockID(iBlk);
		const int iParam = itr->second;
		switch (iParam)
		{
		case ResistivityBlockAnisotropic::RHO_XX:
			vector[index] = log10(anisoParams.rhoXX);
			break;
		case ResistivityBlockAnisotropic::RHO_YY:
			vector[index] = log10(anisoParams.rhoYY);
			break;
		case ResistivityBlockAnisotropic::RHO_ZZ:
			vector[index] = log10(anisoParams.rhoZZ);
			break;
		case ResistivityBlockAnisotropic::STRIKE:
			vector[index] = anisoParams.strike;
			break;
		case ResistivityBlockAnisotropic::DIP:
			vector[index] = anisoParams.dip;
			break;
		case ResistivityBlockAnisotropic::SLANT:
			vector[index] = anisoParams.slant;
			break;
		default:
			break;
		}
	}

}

// Copy previous anisotropic resistivity parameters to current ones
void ResistivityBlockAnisotropic::copyUnfixedResistivityParametersPreToVector(double* vector) const {

	int index(0);
	for (std::vector< std::pair<int, int> >::const_iterator itr = m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.begin();
		itr != m_modelIndexToBlockAndAnisotropicResistivityParameterIndexes.end(); ++itr, ++index) {
		const int iBlk = itr->first;
		const AnistropicResistivityParameters anisoParams = getAnisotropicResistivityParametersPreFromBlockID(iBlk);
		const int iParam = itr->second;
		switch (iParam)
		{
		case ResistivityBlockAnisotropic::RHO_XX:
			vector[index] = log10(anisoParams.rhoXX);
			break;
		case ResistivityBlockAnisotropic::RHO_YY:
			vector[index] = log10(anisoParams.rhoYY);
			break;
		case ResistivityBlockAnisotropic::RHO_ZZ:
			vector[index] = log10(anisoParams.rhoZZ);
			break;
		case ResistivityBlockAnisotropic::STRIKE:
			vector[index] = anisoParams.strike;
			break;
		case ResistivityBlockAnisotropic::DIP:
			vector[index] = anisoParams.dip;
			break;
		case ResistivityBlockAnisotropic::SLANT:
			vector[index] = anisoParams.slant;
			break;
		default:
			break;
		}
	}

}

// Copy current anisotropic resistivity parameters to previous ones
void ResistivityBlockAnisotropic::copyUnfixedResistivityValuesCurToPre(){

	const int nBlk = getNumResistivityBlockTotal();
	for (int iBlk = 0; iBlk < nBlk; ++iBlk) {
		const AnistropicResistivityParameters anisoParamsCur = getAnisotropicResistivityParametersFromBlockID(iBlk);
		AnistropicResistivityParameters anisoParamsPre = anisoParamsCur;
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_XX)) {
			anisoParamsPre.rhoXX = anisoParamsCur.rhoXX;
			if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::ISOTROPY) {
				anisoParamsPre.rhoYY = anisoParamsPre.rhoXX;
				anisoParamsPre.rhoZZ = anisoParamsPre.rhoXX;
			}
			else if (getTypeOfAnisotropy(iBlk) == ResistivityBlockAnisotropic::TRANSVERSE_ISOTROPY) {
				anisoParamsPre.rhoZZ = anisoParamsPre.rhoXX;
			}
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_YY)) {
			anisoParamsPre.rhoYY = anisoParamsCur.rhoYY;
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::RHO_ZZ)) {
			anisoParamsPre.rhoZZ = anisoParamsCur.rhoZZ;
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::STRIKE)) {
			anisoParamsPre.strike = anisoParamsCur.strike;
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP)) {
			anisoParamsPre.dip = anisoParamsCur.dip;
		}
		if (!isFixedAnisotropicResistivityParameters(iBlk, ResistivityBlockAnisotropic::DIP)) {
			anisoParamsPre.slant = anisoParamsCur.slant;
		}
		m_anisotropicResistivityParametersPre[iBlk] = anisoParamsPre;
	}
}

// Output resistivity data to VTK file
void ResistivityBlockAnisotropic::outputResistivityDataToVTK() const {

	outputResistivityBlockIDsToVTK();
	outputAnisotropyTypesToBinary();

}

// Output resistivity data to binary file
void ResistivityBlockAnisotropic::outputResistivityDataToBinary() const {

	outputResistivityBlockIDsToBinary();
	outputAnisotropyTypesToBinary();

}
