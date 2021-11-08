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
//-------------------------------------------------------------------------------------------------------.
#ifndef DBLDEF_DOUBLE_SPARSE_MATRIX
#define DBLDEF_DOUBLE_SPARSE_MATRIX

#include <map>

class DoubleSparseMatrix{

public:

	//Default Constructer
	explicit DoubleSparseMatrix();

	//Constructer
	explicit DoubleSparseMatrix( const int nrows, const int ncols, const int nrhs = 1 );

	//Destructer
	virtual ~DoubleSparseMatrix();

	// Set number of rows and columns
	virtual void setNumRowsAndColumns( const int nrows, const int ncols );

	// Set matrix structure ( locations of non-zero components ) by triplet format
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const double val );

	//Convert matrix from triplet format to CRS format
	void convertToCRSFormat();

	//Add non-zero value to matrix
	//virtual void addNonZeroValues( const int row, const int col, const std::complex<double> val );
	virtual void addNonZeroValues( const int row, const int col, const double val );

	//Zero clear non-zero values of matrix stored by CSR format
	void zeroClearNonZeroValues();

	//Add non-zero value to the right hand side vector
	void addRightHandSideVector( const int row, const double val, const int irhs = 0 );

	//Zero clear non-zero values of the right hand side vector
	void zeroClearRightHandSideVector();

	//Initialize matrix and right-hand side vectors
	virtual void initializeMatrixAndRhsVectors( const int nrows, const int ncols, const int nrhs );

	// Get number of rows
	int getNumRows() const;

	// Get number of columns
	int getNumColumns() const;

	//Get total number of right hand side vector 
	int getNumRightHandSideVectors() const;

	//Return whether matrix has already been converted to CRS format
	bool hasConvertedToCRSFormat() const;

	//Reallocate memory for right hand side vectors
	void reallocateMemoryForRightHandSideVectors( const int nrhs );

	//Release memory
	virtual void releaseMemory();

	////Multiply matrix by inputed vector and add calculated vector to another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndAddResult( const double* vecIn, const int* convertArray, double* vecOut ) const;

	////Multiply matrix by inputed vector and subtract calculated vector from another vector for the case the indexes of inputed vector does not correspond to colunm number
	//void multiplyMatrixByVectorAndSubtractResult( const double* vecIn, const int* convertArray, double* vecOut ) const;

	//Copy right-hand side vector to another vector
	void copyRhsVector( double* vecOut ) const;

	//Debug write the matrix componets
	void debugWriteMatrix() const;

	//Debug write the componets of right hand side vector 
	void debugWriteRightHandSide() const;

	//Get row indices of the compressed row storage format
	int getRowIndexCRS( const int iRow ) const;

	//Get column in which non-zero compnents exist
	int getColumnsCRS( const int iNonZero ) const;

	//Get value of non-zero compnents
	double getValueCRS( const int iNonZero ) const;

	//Get right hand side vector
	double getRightHandSideVector( const int row, const int irhs = 0 ) const;

	// Calculate matrix-vector product of coefficient matrix and inputted vector
	void calcMatrixVectorProduct( const double* invVec, double* outVec ) const;

	// Calculate matrix-vector product of transposed coefficient matrix and inputted vector
	void calcMatrixVectorProductUsingTransposedMatrix( const double* invVec, double* outVec ) const;

protected:
	//Delete the matrix of triplet ( Coordinate ) format
	void deleteTripletMatrix();

	//Total number of rows
	int m_numRows;

	//Total number of columns
	int m_numColumns;

	//Total number of non-zero elements in the matrix
	int m_numNonZeros;

	//Total number of right hand side vectors
	int m_numRightHandSideVectors;

	//Flag indicating whether matrix has already been converted to CRS format
	bool m_hasConvertedToCRSFormat;

	//Matrix stored by triplet ( Coordinate ) format
	std::map< int, double >* m_matrixTripletFormat;

	//Row indices of the compressed row storage format
	long long int* m_rowIndex;

	//Columns in which non-zero compnents exist
	long long int* m_columns;

	//Values of non-zero compnents
	double* m_values;

	//Right hand side vector
	double* m_rightHandSideVector;
	
private:
	//Copy constructer
	DoubleSparseMatrix(const DoubleSparseMatrix &matrix );

	// Assignment operator
	DoubleSparseMatrix& operator=(const DoubleSparseMatrix& rhs);

};

#endif
