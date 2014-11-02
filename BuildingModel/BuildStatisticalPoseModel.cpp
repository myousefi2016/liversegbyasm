/*
* This file is part of the statismo library.
*
* Author: Marcel Luethi (marcel.luethi@unibas.ch)
*
* Copyright (c) 2011 University of Basel
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
*
* Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
*
* Neither the name of the project's author nor the names of its
* contributors may be used to endorse or promote products derived from
* this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "itkMatrix.h"
#include "itkEuler3DTransform.h"

#include "kmUtility.h"

//计算高斯拟合参数(mean,sigma)
template<class VectorType, class MatrixType>
void
calcuateGaussParameters( std::vector<VectorType> & samples, unsigned int N, VectorType& mean, MatrixType& sigma )
{
	mean.Fill( 0 );
	for (int i=0;i<N;i++)
	{
		mean += samples[i];
	}
	mean /= N;

	sigma.Fill( 0 );
	for (int i=0;i<N;i++)
	{
		// siama = ( xk-u )(xk-u)t
		VectorType sam = samples[i] - mean;
		MatrixType sig = sam.GetVnlMatrix() * sam.GetTranspose();

		sigma += sig;
	}
	sigma /= (N-1);
}

void buildPoseModel(const char* datalist, const char* modelname) 
{
	//////////////////////////////////////////////////////////////////////////
	//导入训练数据文件
	KM_DEBUG_INFO("Import Datalist");
	typedef std::vector<std::string>    StringVectorType;
	StringVectorType filenames;
	int numberOfData = km::getDataList(datalist, filenames);

	KM_DEBUG_PRINT( "Number of data: ", numberOfData );

	typedef itk::Matrix<double, 3, 3> Matrix3x3Type;
	typedef itk::Matrix<double, 3, 1> Matrix3x1Type;

	typedef std::vector<Matrix3x1Type> SampleVectorType;
	SampleVectorType posvector;  //存放所有肝脏质心位置样本
	SampleVectorType anglevector; //存放所有肝脏旋转角度样本
	
	typedef itk::Euler3DTransform<double> RigidTransformType;

	ofstream sampleListStream;
	sampleListStream.open( "posSamples.txt" );

	//高斯分布拟合
	for (unsigned int i=0;i<numberOfData;i++)
	{
		RigidTransformType::Pointer rigidTransform = km::readTransform<RigidTransformType>( filenames[i].c_str() );

		Matrix3x1Type possample;//质心位置样本
		RigidTransformType::OutputVectorType position = rigidTransform->GetTranslation();
		possample[0][0] = position[0];
		possample[1][0] = position[1];
		possample[2][0] = position[2];
		posvector.push_back( possample );

		Matrix3x1Type anglesample; //旋转角度样本
		double anglex = rigidTransform->GetAngleX();
		double angley = rigidTransform->GetAngleY();
		double anglez = rigidTransform->GetAngleZ();

		anglesample[0][0] = anglex;
		anglesample[1][0] = angley;
		anglesample[2][0] = anglez;
		anglevector.push_back( anglesample );

		sampleListStream << anglex << "\t" << angley << "\t" << anglez << "\t" << position[0] << "\t" <<position[1] << "\t" << position[2] << std::endl;
	}

	sampleListStream.close();

	ofstream posfile;
	posfile.open ( modelname );

	Matrix3x1Type posMean;
	Matrix3x3Type posSigma;
	calcuateGaussParameters<Matrix3x1Type, Matrix3x3Type>( posvector, numberOfData, posMean, posSigma  );
	posfile << "#pos-mean" << std::endl;
	posfile << posMean << std::endl;
	posfile << "#pos-sigma" << std::endl;
	posfile << posSigma << std::endl;
	KM_DEBUG_PRINT( "pos mean: ", posMean );
	KM_DEBUG_PRINT( "pos sigma: ", posSigma );

	Matrix3x1Type angleMean;
	Matrix3x3Type angleSigma;
	calcuateGaussParameters<Matrix3x1Type, Matrix3x3Type>( anglevector, numberOfData, angleMean, angleSigma );
	posfile << "#angle-mean" << std::endl;
	posfile << angleMean << std::endl;
	posfile << "#angle-sigma" << std::endl;
	posfile << angleSigma << std::endl;
	KM_DEBUG_PRINT( "angle mean: ", angleMean );
	KM_DEBUG_PRINT( "angle sigma: ", angleSigma );

	posfile.close();
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		std::cout << "usage " << argv[0] << " transformlist modelname(*.txt)" << std::endl;
		exit(-1);
	}

	const char* transformlist = argv[1];
	const char* modelname = argv[2];

	buildPoseModel(transformlist, modelname);

	std::cout << "Model building is completed successfully." << std::endl;
}

