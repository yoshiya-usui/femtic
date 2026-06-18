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
#include "Forward3DTetraElement0thOrderAnisotropic.h"
#include "OutputFiles.h"
#include "AnalysisControl.h"
#include "ResistivityBlockAnisotropic.h"

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <assert.h>

Forward3DTetraElement0thOrderAnisotropic::Forward3DTetraElement0thOrderAnisotropic():
	Forward3DTetraElement0thOrder()
{
}

//Destructer
Forward3DTetraElement0thOrderAnisotropic::~Forward3DTetraElement0thOrderAnisotropic(){
}

// Calculate vector x of the reciprocity algorithm of Rodi (1976) for isotropic conductivity
void Forward3DTetraElement0thOrderAnisotropic::calVectorXOfReciprocityAlgorithmForIsotropicConductivity(const std::complex<double>* const vecIn, const int blkID, std::complex<double>* const vecOut, std::vector<int>& nonZeroRows) {

	OutputFiles::m_logFile << "Error : " << __FUNCTION__ << " is not implemented" << std::endl;
	exit(1);

}

// Set non-zero values of matrix and right-hande side vector for forward calculation
void Forward3DTetraElement0thOrderAnisotropic::setNonZeroValues( ComplexSparseSquareSymmetricMatrix& matrix ){

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
	const double omega = 2.0 * CommonParameters::PI * freq;//Angular frequency

	//------------------------------------------
	//--- Components due to stiffness matrix ---
	//------------------------------------------
	const int nElem = m_MeshDataTetraElement.getNumElemTotal();
	for (int iElem = 0; iElem < nElem; ++iElem) {
		const int elemID = iElem;

		//--- Calculate Jacobian
		Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double detJacob(0.0);
		calcJacobianMatrix(elemID, jacobMat, detJacob);
		const double divDetJacob = 1.0 / detJacob;
		//--- Calculate inverse of Jacobian matrix multiplied by determinant
		Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		calcInverseOfJacobianMatrix(jacobMat, invJacobMat);

		double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
		for (int i = 0; i < 6; ++i) {
			length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge(elemID, i);
		}

		//--- Calculate omega * mu
		const double factor1 = 1.0;
		const std::complex<double> factor2 = std::complex<double>(0.0, omega * CommonParameters::mu);// exp(-i*omega*t) form

		//--- Calculate conductivity tensor
		double conductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcAnisotropicConductivityTensor(pResistivityBlock->getBlockIDFromElemID(iElem), conductivityTensor);

		for (int iEdge1 = 0; iEdge1 < 6; ++iEdge1) {
			const int row = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge1]];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 6; ++iEdge2) {
				const int col = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge2]];
				if (col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					integral1 += ((getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11
						+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21
						+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31)
						+ (getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12
							+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22
							+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32)
						+ (getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13
							+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23
							+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33)) * m_weights[ip];
					const double Nx = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat11
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat12
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat13;
					const double Ny = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat21
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat22
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat23;
					const double Nz = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat31
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat32
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat33;
					integral2 += ((getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat11
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat12
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat13)
						* (Nx * conductivityTensor[0][0] + Ny * conductivityTensor[0][1] + Nz * conductivityTensor[0][2])
						+ (getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat21
							+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat22
							+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat23)
						* (Nx * conductivityTensor[1][0] + Ny * conductivityTensor[1][1] + Nz * conductivityTensor[1][2])
						+ (getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat31
							+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat32
							+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat33)
						* (Nx * conductivityTensor[2][0] + Ny * conductivityTensor[2][1] + Nz * conductivityTensor[2][2])) * m_weights[ip];
				}

				if (m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2]) {
					integral1 *= -(divDetJacob * length[iEdge1] * length[iEdge2]);
					integral2 *= -(divDetJacob * length[iEdge1] * length[iEdge2]);
				}
				else {
					integral1 *= (divDetJacob * length[iEdge1] * length[iEdge2]);
					integral2 *= (divDetJacob * length[iEdge1] * length[iEdge2]);
				}

				const std::complex<double> val = std::complex<double>(integral1 * factor1, 0.0) - std::complex<double>(integral2, 0.0) * factor2;// exp(-i*omega*t) form

				if (col == DIRICHLET_BOUNDARY_NONZERO_VALUE) {
					matrix.addRightHandSideVector(row, -val * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]]);// Add to right hand side vector
				}
				else if (col >= row) {// Store only upper triangle part
					const int loc = matrix.checkAndGetLocationNonZeroValue(row, col);
					matrix.addNonZeroValuesWithoutSearchingLocation(loc, val);// Add to matrix
				}

			}// iEdge2
		}// iEdge1		
	}// iElem

#ifdef _DEBUG_WRITE
	std::cout << "matrix" << std::endl;
	matrix.debugWriteMatrix();
	matrix.debugWriteRightHandSide();
#endif

}

// Calculate vector x of the reciprocity algorithm of Rodi (1976) for anisotropic conductivity
void Forward3DTetraElement0thOrderAnisotropic::calVectorXOfReciprocityAlgorithmForAnisotropicConductivity(const std::complex<double>* const vecIn, const int blkID, const int paramID,
	std::complex<double>* const vecOut, std::vector<int>& nonZeroRows) {

	assert((AnalysisControl::getInstance())->isAnisotropyConsidered());

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();

	const int iPol = getPolarizationCurrent();
	const double freq = getFrequencyCurrent();
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
		//--- Calculate Jacobian
		Forward3DTetraElement0thOrder::Matrix3x3 jacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double detJacob(0.0);
		calcJacobianMatrix(elemID, jacobMat, detJacob);
		const double divDetJacob = 1.0 / detJacob;
		//--- Calculate inverse of Jacobian matrix multiplied by determinant
		Forward3DTetraElement0thOrder::Matrix3x3 invJacobMat = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		calcInverseOfJacobianMatrix(jacobMat, invJacobMat);

		double length[6] = { 0.0,0.0,0.0,0.0,0.0,0.0 };
		for (int i = 0; i < 6; ++i) {
			length[i] = m_MeshDataTetraElement.calcEdgeLengthFromElementAndEdge(elemID, i);
		}

		//--- Calculate omega * mu
		const double factor1 = 0.0;
		const std::complex<double> factor2 = std::complex<double>(0.0, omega * CommonParameters::mu);// exp(-i*omega*t) form
		//--- Calculate derivative of the conductivity tensor
		double derivConductivityTensor[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
		pResistivityBlock->calcDerivativeOfAnisotropicConductivityTensor(blkID, paramID, derivConductivityTensor);

		for (int iEdge1 = 0; iEdge1 < 6; ++iEdge1) {
			const int row = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge1]];
			if (row < 0) {
				continue;
			}
			for (int iEdge2 = 0; iEdge2 < 6; ++iEdge2) {
				const int col = m_IDsGlobal2AfterDegenerated[iPol][m_IDsLocal2Global[elemID][iEdge2]];
				if (col <= Forward3DTetraElement0thOrder::DIRICHLET_BOUNDARY_ZERO_VALUE) {
					continue;
				}
				double integral1 = 0.0;
				double integral2 = 0.0;
				for (int ip = 0; ip < m_numIntegralPoints; ++ip) {
					integral1 += ((getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat11
						+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat21
						+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat31)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat11
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat21
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat31)
						+ (getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat12
							+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat22
							+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat32)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat12
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat22
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat32)
						+ (getShapeFuncRotatedReferenceCoordU(iEdge1) * jacobMat.mat13
							+ getShapeFuncRotatedReferenceCoordV(iEdge1) * jacobMat.mat23
							+ getShapeFuncRotatedReferenceCoordW(iEdge1) * jacobMat.mat33)
						* (getShapeFuncRotatedReferenceCoordU(iEdge2) * jacobMat.mat13
							+ getShapeFuncRotatedReferenceCoordV(iEdge2) * jacobMat.mat23
							+ getShapeFuncRotatedReferenceCoordW(iEdge2) * jacobMat.mat33)) * m_weights[ip];
					const double Nx = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat11
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat12
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat13;
					const double Ny = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat21
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat22
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat23;
					const double Nz = getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat31
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat32
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge2) * invJacobMat.mat33;
					integral2 += ((getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat11
						+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat12
						+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat13)
						* (Nx * derivConductivityTensor[0][0] + Ny * derivConductivityTensor[0][1] + Nz * derivConductivityTensor[0][2])
						+ (getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat21
							+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat22
							+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat23)
						* (Nx * derivConductivityTensor[1][0] + Ny * derivConductivityTensor[1][1] + Nz * derivConductivityTensor[1][2])
						+ (getShapeFuncReferenceCoordU(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat31
							+ getShapeFuncReferenceCoordV(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat32
							+ getShapeFuncReferenceCoordW(m_uCoord[ip], m_vCoord[ip], m_wCoord[ip], iEdge1) * invJacobMat.mat33)
						* (Nx * derivConductivityTensor[2][0] + Ny * derivConductivityTensor[2][1] + Nz * derivConductivityTensor[2][2])) * m_weights[ip];
				}

				if (m_signInversion[elemID][iEdge1] != m_signInversion[elemID][iEdge2]) {
					integral1 *= -(divDetJacob * length[iEdge1] * length[iEdge2]);
					integral2 *= -(divDetJacob * length[iEdge1] * length[iEdge2]);
				}
				else {
					integral1 *= (divDetJacob * length[iEdge1] * length[iEdge2]);
					integral2 *= (divDetJacob * length[iEdge1] * length[iEdge2]);
				}

				const std::complex<double> val = std::complex<double>(integral1 * factor1, 0.0) - std::complex<double>(integral2, 0.0) * factor2;// exp(-i*omega*t) form

				if (col == DIRICHLET_BOUNDARY_NONZERO_VALUE) {
					vecOut[row] -= val * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]];
					nonZeroRows.push_back(row);
#ifdef _DEBUG_WRITE
					rhsTemp[row] -= val * m_globalID2NonZeroValues[m_IDsLocal2Global[elemID][iEdge2]];
#endif
				}
				else {
					vecOut[row] -= val * vecIn[col];
					nonZeroRows.push_back(row);
#ifdef _DEBUG_WRITE
					matrixTemp[row][col] += val;
#endif
				}
			}// iEdge2
		}// iEdge1		
	}// iElem

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
CommonParameters::ComplexValuedVector Forward3DTetraElement0thOrderAnisotropic::calculateElectricCurrentDensityVector(const int iElem) const {

	const ResistivityBlockAnisotropic* const pResistivityBlock = (AnalysisControl::getInstance())->getPointerOfResistivityBlockAnisotropic();
	const CommonParameters::ComplexValuedVector electricField = {
		calcValueElectricFieldXDirection(iElem, 0.25, 0.25, 0.25),
		calcValueElectricFieldYDirection(iElem, 0.25, 0.25, 0.25),
		calcValueElectricFieldZDirection(iElem, 0.25, 0.25, 0.25)
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
