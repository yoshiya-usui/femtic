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
#include <stddef.h> // For null pointer
#include <stdlib.h> // For exit
#include <string.h> // For memcpy
#include <iostream>
#include <assert.h>
#include "RougheningSquareMatrix.h"
#include "OutputFiles.h"
#include "AnalysisControl.h"

#include <iomanip>
#include <algorithm>
#include "mkl_lapacke.h"

// Default Constructer
RougheningSquareMatrix::RougheningSquareMatrix():
	DoubleSparseSquareUnsymmetricMatrix()
{}

// Constructer
RougheningSquareMatrix::RougheningSquareMatrix( const int nEq, const int nRhs ):
	DoubleSparseSquareUnsymmetricMatrix( nEq, nRhs )
{}

// Destructer
RougheningSquareMatrix::~RougheningSquareMatrix(){
}

// Set matrix structure ( locations of non-zero components ) by triplet format
// Note : This function must be called BEFORE the matrix is converted into CRS format
void RougheningSquareMatrix::setStructureByTripletFormat( const int row, const int col ){
	DoubleSparseMatrix::setStructureByTripletFormat( row, col );
}

// Set matrix structure ( locations of non-zero components ) and add values by triplet format
void RougheningSquareMatrix::setStructureAndAddValueByTripletFormat( const int row, const int col, const double val ){
	DoubleSparseMatrix::setStructureAndAddValueByTripletFormat( row, col, val );
}

// Add non-zero value to matrix
// Note : This function must be called AFTER the matrix is converted into CRS format
void RougheningSquareMatrix::addNonZeroValues( const int row, const int col, const double val ){
	DoubleSparseMatrix::addNonZeroValues( row, col, val );
}

//Make [R]T[R] matrix, where [R] is a constraining matrix
void RougheningSquareMatrix::makeRTRMatrix( DoubleSparseSquareSymmetricMatrix& RTRMatrix, const double smallValueOnDiagonals ) const{

	assert(!m_matrixTripletFormat);

	RTRMatrix.setDegreeOfEquation(m_numRows);

	for( int iRow = 0; iRow < m_numRows; ++iRow ){
		for( int iColRight = m_rowIndex[iRow]; iColRight < m_rowIndex[iRow+1]; ++iColRight ){
			const int col = m_columns[iColRight];
			for( int iColLeft = m_rowIndex[iRow]; iColLeft <= iColRight; ++iColLeft ){
				const int row = m_columns[iColLeft];
				const double value = m_values[iColRight] * m_values[iColLeft];
				RTRMatrix.setStructureAndAddValueByTripletFormat(row, col, value);
			}
		}
	}

	for( int iRow = 0; iRow < m_numRows; ++iRow ){
		RTRMatrix.setStructureAndAddValueByTripletFormat(iRow, iRow, smallValueOnDiagonals);
	}

	RTRMatrix.convertToCRSFormat();

}

// Calculate model roughness
double RougheningSquareMatrix::calcModelRoughness( const double* modelVec ) const{

	assert(!m_matrixTripletFormat);

	double modelRoughness(0.0);
	for( int i = 0; i < m_numRows; ++i ){
		double work(0.0);
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			 work += m_values[j] * modelVec[ m_columns[j] ];
		}
		modelRoughness += pow( m_rightHandSideVector[i] - work, 2 );
	}
	return modelRoughness;
	
}

// Calculate model roughness for difference filter
double RougheningSquareMatrix::calcModelRoughnessForDifferenceFilter( const double* modelVec ) const{

	assert(!m_matrixTripletFormat);

	double modelRoughness(0.0);
	for( int i = 0; i < m_numRows; ++i ){
		double work(0.0);
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			 work += m_values[j] * modelVec[ m_columns[j] ];
		}
		modelRoughness += pow( m_rightHandSideVector[i] - work, 2 );
	}
	return modelRoughness;

}

// Calculate vector of model roughness
void RougheningSquareMatrix::calcVectorOfModelRoughness( const double* modelVec, double* roughnessVec ) const{

	assert(!m_matrixTripletFormat);

	calcMatrixVectorProduct(modelVec, roughnessVec);

	for( int i = 0; i < m_numRows; ++i ){
		roughnessVec[i] = m_rightHandSideVector[i] - roughnessVec[i];
	}

}

// Postmultiply diagonal matrix
void RougheningSquareMatrix::postmultiplyDiagonalMatrix( const double* diagMatrix ){

	assert(!m_matrixTripletFormat);

	for( int i = 0; i < m_numRows; ++i ){
		for( int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			 m_values[j] *= diagMatrix[ m_columns[j] ];
		}
	}

}

//Copy constructer
RougheningSquareMatrix::RougheningSquareMatrix(const RougheningSquareMatrix &matrix ){
	std::cerr << "Error : Copy constructer of the class RougheningSquareMatrix is not implemented." << std::endl;
	exit(1);
}

// Assignment operator
RougheningSquareMatrix& RougheningSquareMatrix::operator=(const RougheningSquareMatrix& rhs){
	std::cerr << "Error : Assignment operator of the class RougheningSquareMatrix is not implemented." << std::endl;
	exit(1);
}

//Calculaste singular values
void RougheningSquareMatrix::calcSingularValues() const{

	const int myProcessID = (AnalysisControl::getInstance())->getMyPE();
	if( myProcessID != 0 ){
		return;
	}

	assert(!m_matrixTripletFormat);

	long long int lda = m_numRows;
	double* a = new double[lda*m_numRows];
	double* work = new double[m_numRows];
	long long int* lwork = new long long int[m_numRows];
	double* d = new double[m_numRows];
	double* e = new double[m_numRows];
	double* tauq = new double[m_numRows];
	double* taup = new double[m_numRows];
	for( long long int irow = 0; irow < m_numRows; ++irow ){
		for( long long int icol = 0; icol < m_numRows; ++icol ){
			a[ icol + irow * m_numRows ] = 0.0;
		}
	}
	for( long long int i = 0; i < m_numRows; ++i ){
		for( long long int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			a[ m_columns[j] + i * m_numRows ] = m_values[j];
		}
	}
	for( long long int irow = 0; irow < m_numRows; ++irow ){
		for( long long int icol = 0; icol < m_numRows; ++icol ){
			std::cout << " " << a[ icol + irow * m_numRows ];
		}
		std::cout << std::endl;
	}
	for( long long int irow = 0; irow < m_numRows; ++irow ){
		for( long long int icol = 0; icol < m_numRows; ++icol ){
			std::cout << "row col val : " << irow << " " << icol << " " << a[ icol + irow * m_numRows ] << std::endl;
		}
	}
	long long int m_numRows_64 = static_cast<long long int>(m_numRows);

	LAPACKE_dgebrd( LAPACK_ROW_MAJOR, m_numRows_64, m_numRows_64, a, lda, d, e, tauq, taup );
	
	long long int ncvt= 0;
	long long int nru = 0;
	long long int ncc = 0;
	long long int ldvt = m_numRows_64;
	long long int ldu = m_numRows_64;
	long long int ldc = m_numRows_64;
	double* vt = NULL;
	double* u = NULL;
	double* c = NULL;
	double* q = NULL;
	long long int* iq = NULL;
	//LAPACKE_dbdsqr( LAPACK_ROW_MAJOR, 'U', m_numRows_64, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc );
	LAPACKE_dbdsdc( LAPACK_ROW_MAJOR, 'U', 'N', m_numRows_64, d, e, u, ldu, vt, ldvt, q, iq );

	std::ofstream fout( "singularvalues.txt" );
	fout.precision(6);
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){		;
		fout << "row d : " << irow << std::setw(15) << std::scientific << d[irow] << std::endl;
	}	

	std::vector<double> valuesAbs;
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){
		valuesAbs.push_back(std::abs(d[irow]));
	}	

	std::sort(valuesAbs.begin(), valuesAbs.end());
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){		;
		fout << "row abs(d) : " << irow << std::setw(15) << std::scientific << valuesAbs[irow] << std::endl;
	}	

	delete [] a;
	delete [] work;
	delete [] lwork;
	delete [] d;
	delete [] e;
	delete [] tauq;
	delete [] taup;

}

//Calculaste eigen values
void RougheningSquareMatrix::calcEigenValues() const{

	const int myProcessID = (AnalysisControl::getInstance())->getMyPE();
	if( myProcessID != 0 ){
		return;
	}

	assert(!m_matrixTripletFormat);

	long long int m_numRows_64 = static_cast<long long int>(m_numRows);
	double* ap = NULL;
	double* d = NULL;
	double* e = NULL;
	double* tau = NULL;
	ap = new double[m_numRows*(m_numRows+1)/2];
	d = new double[m_numRows];
	e = new double[m_numRows-1];
	tau = new double[m_numRows-1];

	int icount(0);
	for( long long int irow = 0; irow < m_numRows; ++irow ){
		for( long long int icol = irow; icol < m_numRows; ++icol ){
			ap[icount] = 0.0;
			++icount;
		}
	}
	icount = 0;
	for( long long int i = 0; i < m_numRows; ++i ){
		for( long long int j = m_rowIndex[i]; j < m_rowIndex[i+1]; ++j ){		
			ap[ m_columns[j] + icount - i ] = m_values[j];
		}
		icount += m_numRows - i;
	}
	icount = 0;
	for( long long int irow = 0; irow < m_numRows; ++irow ){
		for( long long int icol = irow; icol < m_numRows; ++icol ){
			std::cout << "row col ap : " << irow << " " << icol << " " << ap[icount] << std::endl;
			++icount;
		}
	}	

	LAPACKE_dsptrd( LAPACK_ROW_MAJOR, 'U', m_numRows_64, ap, d, e, tau );
	LAPACKE_dsterf( m_numRows_64, d, e );

	std::ofstream fout( "eigenvalues.txt" );
	fout.precision(6);
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){		;
		fout << "row d : " << irow << std::setw(15) << std::scientific << d[irow] << std::endl;
	}	

	std::vector<double> valuesAbs;
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){
		valuesAbs.push_back(std::abs(d[irow]));
	}	

	std::sort(valuesAbs.begin(), valuesAbs.end());
	for( long long int irow = 0; irow < m_numRows_64; ++irow ){		;
		fout << "row abs(d) : " << irow << std::setw(15) << std::scientific << valuesAbs[irow] << std::endl;
	}	
	
	delete [] ap;
	delete [] d;
	delete [] e;
	delete [] tau;
}

// Output roughening matrix
void RougheningSquareMatrix::outputRougheningMatrix() const{

	assert( !m_hasConvertedToCRSFormat );

	std::string fileName="roughening_matrix.out";
	std::ofstream ofs( fileName.c_str() );

	if( ofs.fail() ){
		OutputFiles::m_logFile << "File open error : " << fileName.c_str() << " !!" << std::endl;
		exit(1);
	}

	ofs << std::setw(10) << m_numRows << std::endl;

	for( int irow = 0; irow < m_numRows; ++irow ){
		const std::map< int, double >::iterator itEnd = m_matrixTripletFormat[irow].end();
		const int numNonzeros = static_cast<int>( m_matrixTripletFormat[irow].size() );
		ofs << std::setw(10) << irow << std::endl;
		ofs << std::setw(10) << numNonzeros;
		for( std::map< int, double >::iterator it = m_matrixTripletFormat[irow].begin(); it != itEnd; ++it ){
			ofs << std::setw(15) << it->first;
		}
		ofs << std::endl;
		if( numNonzeros > 0 ){
			ofs << std::setw(10) << "";
			for( std::map< int, double >::iterator it = m_matrixTripletFormat[irow].begin(); it != itEnd; ++it ){
				ofs << std::setw(15) << std::scientific << std::setprecision(6) << it->second;
			}
			ofs << std::endl;
		}
	}

	ofs.close();

}