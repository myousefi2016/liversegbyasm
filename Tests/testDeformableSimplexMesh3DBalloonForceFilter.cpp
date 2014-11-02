/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include <iostream>

#include "itkRegularSphereMeshSource.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkDeformableSimplexMesh3DBalloonForceFilter.h"

#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

#include "kmUtility.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

int main(int argc, char * argv[] )
{
	if(argc<8)
	{
		std::cout<<"Usage:"<<std::endl;
		std::cout<<"       ";
		return EXIT_FAILURE;
	}

	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;
	memorymeter.Start( "starting" );
	chronometer.Start( "starting" );

  // Declare the type of the input and output mesh
  typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
  typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double,double> SimplexMeshTraits;
  typedef itk::Mesh<double,3, TriangleMeshTraits> TriangleMeshType;
  typedef itk::SimplexMesh<double,3, SimplexMeshTraits> SimplexMeshType;

  // declare triangle mesh source
  typedef itk::RegularSphereMeshSource<TriangleMeshType>  SphereMeshSourceType;
  typedef SphereMeshSourceType::PointType PointType;
  typedef SphereMeshSourceType::VectorType VectorType;

  // declare the triangle to simplex mesh filter
  typedef itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType> SimplexFilterType;

  /*SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();
  PointType center;
  //center.Fill(50);
	center[0] = 95;
	center[1] = 145;
	center[2] = 195;
  PointType::ValueType scaleInit[3] = {10,10,10};
  VectorType scale = scaleInit;

  mySphereMeshSource->SetCenter(center);
  mySphereMeshSource->SetResolution( 5 );
  mySphereMeshSource->SetScale(scale);

	TriangleMeshType::Pointer triangleMesh = mySphereMeshSource->GetOutput();*/
	TriangleMeshType::Pointer triangleMesh = km::readMesh<TriangleMeshType>(argv[1]);

  std::cout << "Triangle mesh created. " << std::endl;

 // SimplexFilterType::Pointer simplexFilter = SimplexFilterType::New();
 // simplexFilter->SetInput( triangleMesh );
	//simplexFilter->Update();

	//SimplexMeshType::Pointer simplexMesh = simplexFilter->GetOutput();
	//simplexMesh->DisconnectPipeline();

	SimplexMeshType::Pointer simplexMesh = km::triangleMeshToSimplexMesh<TriangleMeshType, SimplexMeshType>(triangleMesh);
  
  km::writeMesh<SimplexMeshType>("initialSimplexMesh.vtk", simplexMesh);

	typedef itk::Image<float,3>                       OriginalImageType;
	
  typedef itk::DeformableSimplexMesh3DBalloonForceFilter<SimplexMeshType,SimplexMeshType> DeformFilterType;
	typedef DeformFilterType::GradientImageType       GradientImageType;
	
	/*
  std::cout << "Creating dummy image...";

  typedef OriginalImageType::PixelType              PixelType;
  typedef OriginalImageType::IndexType              IndexType;
  typedef OriginalImageType::SizeType               ImageSizeType;

	
  OriginalImageType::Pointer originalImage = OriginalImageType::New();

  ImageSizeType imageSize;
  //imageSize.Fill(100);
	imageSize[0] = 349;
	imageSize[1] = 349;
	imageSize[2] = 342;
  originalImage->SetRegions( imageSize );
  originalImage->Allocate();


  IndexType index;
  for (int x = 0; x < 349; x++)
  {
    for (int y = 0; y < 349; y++)
    {
      for (int z = 0; z < 342; z++)
      {
        index[0] = x;
        index[1] = y;
        index[2] = z;
        if ( ( (x == 65 || x == 125) && y >= 115 && y <= 175 && z >= 115 && z <= 175)  ||
             ( (y == 115 || y == 175) && x >= 65 && x <= 125 && z >= 115 && z <= 175)  ||
             ( (z == 115 || z == 175) && y >= 115 && y <= 175 && x >= 65 && x <= 125)
           )
        {
          originalImage->SetPixel(index, 1);
        }
        else
        {
          originalImage->SetPixel(index, 0);
        }
      }
    }
  }
  
  km::writeImage<OriginalImageType>("originalImage.mha", originalImage);*/

	/*OriginalImageType::Pointer originalImage = km::readImage<OriginalImageType>(argv[1]);

	std::cout << "detect edge image..." << std::endl;
  typedef itk::SobelEdgeDetectionImageFilter<OriginalImageType,OriginalImageType>   EdgeFilterType;
  EdgeFilterType::Pointer edgeFilter = EdgeFilterType::New();
  edgeFilter->SetInput( originalImage );
  edgeFilter->Update();
	OriginalImageType::Pointer edgeImage = edgeFilter->GetOutput();
	km::writeImage<OriginalImageType>( "edgeImage.mha", edgeFilter->GetOutput() );*/

	//OriginalImageType::Pointer edgeImage = km::readImage<OriginalImageType>(argv[2]);
	//edgeImage = km::calculateGradientMagnitudeImage<OriginalImageType, OriginalImageType>(edgeImage, 2);
	//edgeImage = km::thresholdImage<OriginalImageType>(edgeImage,10, 100, 0);
	//km::writeImage<OriginalImageType>( "edgeImage2.mha", edgeImage );

	//std::cout << "calculate gradient image..." << std::endl;
 // typedef DeformFilterType::GradientImageType       GradientImageType;
 // typedef itk::GradientRecursiveGaussianImageFilter<OriginalImageType,GradientImageType> GradientFilterType;

 // GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
 // gradientFilter->SetInput( edgeImage );
 // gradientFilter->SetSigma( 1.0 );
 // gradientFilter->Update();

	//GradientImageType::Pointer gradientImage = gradientFilter->GetOutput();

	//gradientImage = km::calcuateGVF<GradientImageType>(gradientImage, atoi(argv[4]), atoi(argv[5]));

	//km::writeImage<GradientImageType>( "gradientImage.mha", gradientImage );

	GradientImageType::Pointer gradientImage = km::readImage<GradientImageType>(argv[2]);

	std::cout << "start to fitting." << std::endl;
  DeformFilterType::Pointer deformFilter = DeformFilterType::New();
  deformFilter->SetInput( simplexMesh );
  deformFilter->SetGradient( gradientImage );
  deformFilter->SetAlpha(atof(argv[3])); //0.2
  deformFilter->SetBeta(atof(argv[4])); //0.1
  deformFilter->SetKappa(atof(argv[5])); //0.1
  deformFilter->SetIterations(atoi(argv[6])); //100
  deformFilter->SetRigidity(atoi(argv[7])); //1
  deformFilter->Update();

  SimplexMeshType::Pointer deformResult =  deformFilter->GetOutput();

  //std::cout << "Deformation Result: " << deformResult << std::endl;
  
  km::writeMesh<SimplexMeshType>("deformedSimplexMesh.vtk", deformResult);

  std::cout << "[TEST DONE]" << std::endl;

	chronometer.Stop( "starting" );
	memorymeter.Stop( "starting" );
	chronometer.Report( std::cout );
	memorymeter.Report( std::cout );

  return EXIT_SUCCESS;
}
