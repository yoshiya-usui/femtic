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
#include "Forward3DNonConformingHexaElement0thOrderAnisotropic.h"
#include "MeshDataNonConformingHexaElement.h"
#include "AnalysisControl.h"
#include "ResistivityBlockAnisotropic.h"
#include "CommonParameters.h"
#include "OutputFiles.h"
#include "ObservedData.h"
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <assert.h>

Forward3DNonConformingHexaElement0thOrderAnisotropic::Forward3DNonConformingHexaElement0thOrderAnisotropic():
	Forward3DNonConformingHexaElement0thOrder()
{
}

//Destructer
Forward3DNonConformingHexaElement0thOrderAnisotropic::~Forward3DNonConformingHexaElement0thOrderAnisotropic(){
}

// Calculate vector x of the reciprocity algorithm of Rodi (1976) for isotropic conductivity
void Forward3DNonConformingHexaElement0thOrderAnisotropic::calVectorXOfReciprocityAlgorithmForIsotropicConductivity(const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows) {

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);

}

// Is it used for calculating the EM field on a side boundary?
void Forward3DNonConformingHexaElement0thOrderAnisotropic::setNonZeroValues(ComplexSparseSquareSymmetricMatrix& matrix) {

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	const int nElem = m_MeshDataNonConformingHexaElement.getNumElemTotal();
	for (int iElem = 0; iElem < nElem; ++iElem) {
		const int elemID = iElem;
		//--- Calculate omega * mu
		const double factor1 = 1.0;
		const std::complex<double> factor2 = std::complex<double>(0.0, omega * CommonParameters::mu);// exp(-i*omega*t) form
		//--- Calculate conductivity tensor
		double conductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcAnisotropicConductivityTensor(pResistivityBlock->getBlockIDFromElemID(iElem), conductivityTensor);
		double length[12];
		for (int i = 0; i < 12; ++i) {
			length[i] = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge(elemID, i);
		}
		for (int iEdge1 = 0; iEdge1 < 12; ++iEdge1) {
			const int row = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge1]];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 12; ++iEdge2) {
				const int col = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge2]];
				if (col <= Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					const double xi = m_integralPointXi[ip];
					const double eta = m_integralPointEta[ip];
					const double zeta = m_integralPointZeta[ip];
					Forward3D::Matrix3x3 JacobMat;
					const double detJacob = calcJacobianMatrix(elemID, xi, eta, zeta, JacobMat);
					Forward3D::Matrix3x3 invJacobMat;
					calcInverseOfJacobianMatrix(JacobMat, detJacob, invJacobMat);
					integral1 += (getShapeFuncRotatedX(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedX(xi, eta, zeta, iEdge2, invJacobMat)
						+ getShapeFuncRotatedY(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedY(xi, eta, zeta, iEdge2, invJacobMat)
						+ getShapeFuncRotatedZ(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedZ(xi, eta, zeta, iEdge2, invJacobMat))
						* detJacob * m_weights[ip];
					const double Nx = getShapeFuncX(xi, eta, zeta, iEdge2, invJacobMat);
					const double Ny = getShapeFuncY(xi, eta, zeta, iEdge2, invJacobMat);
					const double Nz = getShapeFuncZ(xi, eta, zeta, iEdge2, invJacobMat);
					integral2 += (getShapeFuncX(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * conductivityTensor[0][0] + Ny * conductivityTensor[0][1] + Nz * conductivityTensor[0][2])
						+ getShapeFuncY(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * conductivityTensor[1][0] + Ny * conductivityTensor[1][1] + Nz * conductivityTensor[1][2])
						+ getShapeFuncZ(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * conductivityTensor[2][0] + Ny * conductivityTensor[2][1] + Nz * conductivityTensor[2][2])) * detJacob * m_weights[ip];
				}
				integral1 *= length[iEdge1] * length[iEdge2];
				integral2 *= length[iEdge1] * length[iEdge2];
				const std::complex<double> val = std::complex<double>(integral1 * factor1, 0.0) - std::complex<double>(integral2, 0.0) * factor2;// exp(-i*omega*t) form
				const std::vector< std::pair<int, double> >& rowMasters = m_slaveDofToMasterDofAndFactors[row];
				for (std::vector< std::pair<int, double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow) {
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
					if (col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE) {
						// Add to right hand side vector
						matrix.addRightHandSideVector(rowMod, -valMod * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]]);
					}
					else {
						// Add to right hand side vector corresponding to MPC constants
						matrix.addRightHandSideVector(rowMod, valMod * m_vectorMPCConstants[col]);
						const std::vector< std::pair<int, double> >& colMasters = m_slaveDofToMasterDofAndFactors[col];
						for (std::vector< std::pair<int, double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol) {
							const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
							const std::complex<double> valModMod = valMod * std::complex<double>(itrCol->second, 0.0);
							if (colMod >= rowMod) {// Store only upper triangle part
								const int loc = matrix.checkAndGetLocationNonZeroValue(rowMod, colMod);
								matrix.addNonZeroValuesWithoutSearchingLocation(loc, valModMod);// Add to matrix
							}
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}

#ifdef _DEBUG_WRITE
	std::cout << "matrix" << std::endl;
	matrix.debugWriteMatrix();
	matrix.debugWriteRightHandSide();
#endif

}

// Calculate vector x of the reciprocity algorithm of Rodi (1976) for anisotropic conductivity
void Forward3DNonConformingHexaElement0thOrderAnisotropic::calVectorXOfReciprocityAlgorithmForAnisotropicConductivity(const std::complex<double>* const vecIn, 
	const int blkID, const int paramID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows) {

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	assert(!pResistivityBlock->isFixedAnisotropicResistivityParameters(blkID, paramID));

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double ln10 = 2.30258509299405;
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

#ifdef _DEBUG_WRITE
	std::complex<double>** matrixTemp = new std::complex<double>*[getNumOfEquationFinallySolved()];
	std::complex<double>* rhsTemp = new std::complex<double>[getNumOfEquationFinallySolved()];
	for (int i = 0; i < getNumOfEquationFinallySolved(); ++i) {
		matrixTemp[i] = new std::complex<double>[getNumOfEquationFinallySolved()];
		rhsTemp[i] = 0.0;
		for (int j = 0; j < getNumOfEquationFinallySolved(); ++j) {
			matrixTemp[i][j] = 0.0;
		}
	}
#endif

	const std::vector<int>& mdl2Elem = pResistivityBlock->getBlockID2Elements(blkID);
	const int nElem = static_cast<int>(mdl2Elem.size());
	for (int iElem = 0; iElem < nElem; ++iElem) {
		// [Attention] : You must use elemID instead of iElem from this line
		const int elemID = mdl2Elem[iElem];
		//--- Calculate omega * mu
		const double factor1 = 0.0;
		const std::complex<double> factor2 = std::complex<double>(0.0, omega * CommonParameters::mu);// exp(-i*omega*t) form
		//--- Calculate derivative of the conductivity tensor
		double derivConductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcDerivativeOfAnisotropicConductivityTensor(blkID, paramID, derivConductivityTensor);
		double length[12];
		for (int i = 0; i < 12; ++i) {
			length[i] = m_MeshDataNonConformingHexaElement.calcEdgeLengthFromElementAndEdge(elemID, i);
		}
		for (int iEdge1 = 0; iEdge1 < 12; ++iEdge1) {
			const int row = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge1]];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 12; ++iEdge2) {
				const int col = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge2]];
				if (col <= Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					const double xi = m_integralPointXi[ip];
					const double eta = m_integralPointEta[ip];
					const double zeta = m_integralPointZeta[ip];
					Forward3D::Matrix3x3 JacobMat;
					const double detJacob = calcJacobianMatrix(elemID, xi, eta, zeta, JacobMat);
					Forward3D::Matrix3x3 invJacobMat;
					calcInverseOfJacobianMatrix(JacobMat, detJacob, invJacobMat);
					integral1 += (getShapeFuncRotatedX(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedX(xi, eta, zeta, iEdge2, invJacobMat)
						+ getShapeFuncRotatedY(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedY(xi, eta, zeta, iEdge2, invJacobMat)
						+ getShapeFuncRotatedZ(xi, eta, zeta, iEdge1, invJacobMat)
						* getShapeFuncRotatedZ(xi, eta, zeta, iEdge2, invJacobMat))
						* detJacob * m_weights[ip];
					const double Nx = getShapeFuncX(xi, eta, zeta, iEdge2, invJacobMat);
					const double Ny = getShapeFuncY(xi, eta, zeta, iEdge2, invJacobMat);
					const double Nz = getShapeFuncZ(xi, eta, zeta, iEdge2, invJacobMat);
					integral2 += (getShapeFuncX(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * derivConductivityTensor[0][0] + Ny * derivConductivityTensor[0][1] + Nz * derivConductivityTensor[0][2])
						+ getShapeFuncY(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * derivConductivityTensor[1][0] + Ny * derivConductivityTensor[1][1] + Nz * derivConductivityTensor[1][2])
						+ getShapeFuncZ(xi, eta, zeta, iEdge1, invJacobMat)
						* (Nx * derivConductivityTensor[2][0] + Ny * derivConductivityTensor[2][1] + Nz * derivConductivityTensor[2][2])) * detJacob * m_weights[ip];
				}
				integral1 *= length[iEdge1] * length[iEdge2];
				integral2 *= length[iEdge1] * length[iEdge2];
				const std::complex<double> val = std::complex<double>(integral1 * factor1, 0.0) - std::complex<double>(integral2, 0.0) * factor2;// exp(-i*omega*t) form
				const std::vector< std::pair<int, double> >& rowMasters = m_slaveDofToMasterDofAndFactors[row];
				for (std::vector< std::pair<int, double> >::const_iterator itrRow = rowMasters.begin(); itrRow != rowMasters.end(); ++itrRow) {
					const int rowMod = m_IDsAfterDegenerated2AfterConstrained[itrRow->first];
					const std::complex<double> valMod = val * std::complex<double>(itrRow->second, 0.0);
					if (col == Forward3DNonConformingHexaElement0thOrder::DIRICHLET_BOUNDARY_NONZERO_VALUE) {
						nonZeroRows.push_back(rowMod);
						// Add to right hand side vector
						vecOut[rowMod] -= valMod * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]];
#ifdef _DEBUG_WRITE
						rhsTemp[rowMod] -= valMod * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]];
#endif
					}
					else {
						nonZeroRows.push_back(rowMod);
						// Add to right hand side vector corresponding to MPC constants
						vecOut[rowMod] += valMod * m_vectorMPCConstants[col];
#ifdef _DEBUG_WRITE
						rhsTemp[rowMod] += valMod * m_vectorMPCConstants[col];
#endif
						const std::vector< std::pair<int, double> >& colMasters = m_slaveDofToMasterDofAndFactors[col];
						for (std::vector< std::pair<int, double> >::const_iterator itrCol = colMasters.begin(); itrCol != colMasters.end(); ++itrCol) {
							const int colMod = m_IDsAfterDegenerated2AfterConstrained[itrCol->first];
							const std::complex<double> valModMod = valMod * std::complex<double>(itrCol->second, 0.0);
							vecOut[rowMod] -= valModMod * vecIn[colMod];
#ifdef _DEBUG_WRITE
							matrixTemp[rowMod][colMod] += valModMod;
#endif
						}
					}
				}
			}// iEdge2
		}// iEdge1		
	}

	std::sort(nonZeroRows.begin(), nonZeroRows.end());
	nonZeroRows.erase(std::unique(nonZeroRows.begin(), nonZeroRows.end()), nonZeroRows.end());

#ifdef _DEBUG_WRITE
	std::cout << "derivatives" << std::endl;
	for (int i = 0; i < getNumOfEquationFinallySolved(); ++i) {
		for (int j = 0; j < getNumOfEquationFinallySolved(); ++j) {
			if (i <= j && std::abs(matrixTemp[i][j]) > CommonParameters::EPS) {
				std::cout << "row col val " << i << " " << j << " " << matrixTemp[i][j].real() << " " << matrixTemp[i][j].imag() << std::endl;
			}
		}
	}
	for (int i = 0; i < getNumOfEquationFinallySolved(); ++i) {
		if (std::abs(rhsTemp[i]) > CommonParameters::EPS) {
			std::cout << "row " << i << " " << rhsTemp[i].real() << " " << rhsTemp[i].imag() << std::endl;
		}
	}
	for (int i = 0; i < getNumOfEquationFinallySolved(); ++i) {
		delete[] matrixTemp[i];
	}
	delete[] matrixTemp;
	delete[] rhsTemp;
#endif

}

// Calculate electric current density vector
CommonParameters::ComplexValuedVector Forward3DNonConformingHexaElement0thOrderAnisotropic::calculateElectricCurrentDensityVector(const int iElem) const {

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	const CommonParameters::ComplexValuedVector electricField = {
		calcValueElectricFieldXDirection(iElem, 0.0, 0.0, 0.0),
		calcValueElectricFieldYDirection(iElem, 0.0, 0.0, 0.0),
		calcValueElectricFieldZDirection(iElem, 0.0, 0.0, 0.0)
	};
	double conductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
	pResistivityBlock->calcAnisotropicConductivityTensor(pResistivityBlock->getBlockIDFromElemID(iElem), conductivityTensor);
	const CommonParameters::ComplexValuedVector electricCurrentDensity = {
		conductivityTensor[0][0] * electricField.X + conductivityTensor[0][1] * electricField.Y + conductivityTensor[0][2] * electricField.Z,
		conductivityTensor[1][0] * electricField.X + conductivityTensor[1][1] * electricField.Y + conductivityTensor[1][2] * electricField.Z,
		conductivityTensor[2][0] * electricField.X + conductivityTensor[2][1] * electricField.Y + conductivityTensor[2][2] * electricField.Z
	};
	return electricCurrentDensity;

}
