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
#ifndef DBLDEF_OUTPUT_FILES
#define DBLDEF_OUTPUT_FILES

#include <fstream>
#include <stdio.h>
#include <stdlib.h>

// Class of output files
class OutputFiles{

public:
	// Log file
	static std::ofstream m_logFile;

	// VTK file
	static std::ofstream m_vtkFile;

	// VTK file for ploting observed station
	static std::ofstream m_vtkFileForObservedStation;

	// CSV file
	//static std::ofstream m_csvFile;
	static FILE* m_csvFile;

	// CSV file in which the results of 2D forward calculation is written
	//static std::ofstream m_csvFileFor2DFwd;
	static FILE* m_csvFileFor2DFwd;

	// CNV file
	static std::ofstream m_cnvFile;

	// Return the the instance of the class
	static OutputFiles* getInstance();

	// Close csv file in which the results of 3D forward calculation is written
	void openVTKFile( const int iterNum );

	// Open VTK file for ploting observed station
	void openVTKFileForObservedStation();

	// Open csv file in which the results of 2D forward calculation is written
	void openCsvFileFor2DFwd( const int iterNum );

	// Open csv file in which the results of 3D forward calculation is written
	void openCsvFileFor3DFwd( const int iterNum );

	// Open cnv file
	void openCnvFile();

	// Output case file
	void outputCaseFile() const;

private:
	// Constructer
	OutputFiles();

	// Destructer
	~OutputFiles();

	// Copy constructer
	OutputFiles(const OutputFiles& rhs);

	// Assignment operator
	OutputFiles& operator=(const OutputFiles& rhs);
		
};

#endif
