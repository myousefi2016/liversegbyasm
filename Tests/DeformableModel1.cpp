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

// Software Guide : BeginLatex
//
// This example illustrates the use of the \doxygen{DeformableMesh3DFilter}
// and \doxygen{BinaryMask3DMeshSource} in the hybrid segmentation framework.
//
// \begin{figure} \center
// \includegraphics[width=\textwidth]{DeformableModelCollaborationDiagram.eps}
// \itkcaption[Deformable model collaboration diagram]{Collaboration
// diagram for the DeformableMesh3DFilter applied to a segmentation task.}
// \label{fig:DeformableModelCollaborationDiagram}
// \end{figure}
//
// The purpose of the DeformableMesh3DFilter is to take an initial surface
// described by an \doxygen{Mesh} and deform it in order to adapt it to the
// shape of an anatomical structure in an image.
// Figure~\ref{fig:DeformableModelCollaborationDiagram} illustrates a typical
// setup for a segmentation method based on deformable models. First, an
// initial mesh is generated using a binary mask and an isocontouring
// algorithm (such as marching cubes) to produce an initial mesh. The binary
// mask used here contains a simple shape which vaguely resembles the
// anatomical structure that we want to segment. The application of the
// isocontouring algorithm produces a $3D$ mesh that has the shape of this
// initial structure. This initial mesh is passed as input to the deformable
// model which will apply forces to the mesh points in order to reshape the
// surface until make it fit to the anatomical structures in the image.
//
// The forces to be applied on the surface are computed from an approximate
// physical model that simulates an elastic deformation. Among the forces to
// be applied we need one that will pull the surface to the position of the
// edges in the anatomical structure. This force component is represented
// here in the form of a vector field and is computed as illustrated in the
// lower left of Figure~\ref{fig:DeformableModelCollaborationDiagram}. The
// input image is passed to a
// \doxygen{GradientMagnitudeRecursiveGaussianImageFilter}, which computes
// the magnitude of the image gradient. This scalar image is then passed to
// another gradient filter
// (\doxygen{GradientRecursiveGaussianImageFilter}). The output of this
// second gradient filter is a vector field in which every vector points to
// the closest edge in the image and has a magnitude proportional to the
// second derivative of the image intensity along the direction of the
// gradient. Since this vector field is computed using Gaussian derivatives,
// it is possible to regulate the smoothness of the vector field by playing
// with the value of sigma assigned to the Gaussian. Large values of sigma
// will result in a large capture radius, but will have poor precision in the
// location of the edges. A reasonable strategy may involve the use of large
// sigmas for the initial iterations of the model and small sigmas to refine
// the model when it is close to the edges. A similar effect could be
// achieved using multiresolution and taking advantage of the image pyramid
// structures already illustrated in the registration framework.
//
// \index{Deformable Models}
// \index{DeformableMesh3DFilter}
//
// Software Guide : EndLatex

#include <iostream>
#include "itkBinaryMask3DMeshSource.h"
#include "itkDeformableMesh3DFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkPointSetToImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//  Software Guide : EndCodeSnippet

#include "itkSimpleFilterWatcher.h"

#include "kmUtility.h"

int main( int argc, char *argv[] )
{

	if( argc < 4 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " TargetImage InputMesh DeformedMesh" << std::endl;
		return 1;
	}

	const     unsigned int    Dimension = 3;
	typedef   double                         PixelType;
	typedef itk::Image<PixelType, Dimension> ImageType;

	typedef itk::Image< unsigned char, Dimension >   BinaryImageType;

	typedef  itk::Mesh<double>     MeshType;

	typedef itk::CovariantVector< double, Dimension >  GradientPixelType;
	typedef itk::Image< GradientPixelType, Dimension > GradientImageType;

	typedef itk::GradientRecursiveGaussianImageFilter<ImageType, GradientImageType>
		GradientFilterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,ImageType>
		GradientMagnitudeFilterType;

	typedef itk::BinaryMask3DMeshSource< BinaryImageType, MeshType >  MeshSourceType;

	typedef itk::DeformableMesh3DFilter<MeshType,MeshType>  DeformableFilterType;

	typedef itk::ImageFileReader< ImageType       >  ReaderType;
	typedef itk::ImageFileReader< BinaryImageType >  BinaryReaderType;
	ReaderType::Pointer       imageReader   =  ReaderType::New();
	BinaryReaderType::Pointer maskReader    =  BinaryReaderType::New();

	//imageReader->SetFileName( argv[1] );
	ImageType::Pointer targetImage = km::readImage<ImageType>( argv[1] );

	targetImage = km::binaryThresholdImage<ImageType, ImageType>( targetImage, 1, 1, 200, 0 );


	GradientMagnitudeFilterType::Pointer  gradientMagnitudeFilter
		= GradientMagnitudeFilterType::New();

	gradientMagnitudeFilter->SetInput( targetImage/*imageReader->GetOutput()*/ );
	gradientMagnitudeFilter->SetSigma( 1.0 );

	GradientFilterType::Pointer gradientMapFilter = GradientFilterType::New();

	gradientMapFilter->SetInput( gradientMagnitudeFilter->GetOutput());
	gradientMapFilter->SetSigma( 1.0 );

	try
	{
		// Software Guide : BeginCodeSnippet
		gradientMapFilter->Update();
		// Software Guide : EndCodeSnippet
	}
	catch( itk::ExceptionObject & e )
	{
		std::cerr << "Exception caught when updating gradientMapFilter " << std::endl;
		std::cerr << e << std::endl;
		return -1;
	}

	std::cout << "The gradient map created!" << std::endl;

	km::writeImage<GradientImageType>( "gradient.mha", gradientMapFilter->GetOutput() );

	MeshSourceType::Pointer meshSource = MeshSourceType::New();

	DeformableFilterType::Pointer deformableModelFilter = DeformableFilterType::New();
	deformableModelFilter->SetGradient( gradientMapFilter->GetOutput() );

	MeshType::Pointer initMesh = km::readMesh<MeshType>( argv[2] );

	deformableModelFilter->SetInput(  initMesh/*meshSource->GetOutput()*/ );

	typedef itk::CovariantVector<double, 2>           double2DVector;
	typedef itk::CovariantVector<double, 3>           double3DVector;

	double2DVector stiffness;
	stiffness.Fill( 0.1 );
	//stiffness[0] = 0.0001;
	//stiffness[1] = 0.1;

	double3DVector scale;
	scale.Fill( 1.0 );
	//scale[0] = 1.0;
	//scale[1] = 1.0;
	//scale[2] = 1.0;

	deformableModelFilter->SetStiffness( stiffness );
	deformableModelFilter->SetScale( scale );
	deformableModelFilter->SetGradientMagnitude( 1.0 );
	deformableModelFilter->SetTimeStep( 0.1 );
	deformableModelFilter->SetStepThreshold( 100 );
	deformableModelFilter->SetObjectLabel( 200 );
	deformableModelFilter->SetPotentialOn( 0 );

	std::cout << "Deformable mesh fitting...";

	//std::cout<<initMesh->GetNumberOfPoints()<<std::endl;
	//std::cout<<initMesh->GetNumberOfCells()<<std::endl;

	km::assigneMesh<MeshType>( initMesh, 0 );

	km::writeMesh<MeshType>( "temp.vtk", initMesh );

	itk::SimpleFilterWatcher watcher(deformableModelFilter, "DeformableMeshFilter");

	try
	{
		deformableModelFilter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
		std::cerr << "Exception Caught !" << std::endl;
		std::cerr << excep << std::endl;
	}
	// Software Guide : EndCodeSnippet

	MeshType::Pointer outputMesh = deformableModelFilter->GetOutput();

	km::writeMesh<MeshType>( argv[3], outputMesh );

	/*

	typedef itk::PointSetToImageFilter<MeshType,ImageType> MeshFilterType;
	MeshFilterType::Pointer meshFilter = MeshFilterType::New();
	meshFilter->SetOrigin(mask->GetOrigin());
	meshFilter->SetSize(mask->GetLargestPossibleRegion().GetSize());
	meshFilter->SetSpacing(mask->GetSpacing());
	meshFilter->SetInput(meshSource->GetOutput());
	try
	{
	meshFilter->Update();
	}
	catch( itk::ExceptionObject & excep )
	{
	std::cerr << "Exception Caught !" << std::endl;
	std::cerr << excep << std::endl;
	}
	// Software Guide : EndCodeSnippet

	//  Software Guide : BeginLatex
	//
	//  The resulting deformed binary mask can be written on disk
	//  using the \doxygen{ImageFileWriter}.
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	typedef itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(meshFilter->GetOutput());
	writer->SetFileName(argv[3]);
	writer->Update();*/
	// Software Guide : EndCodeSnippet

	//  Software Guide : BeginLatex
	//
	//  Note that in order to successfully segment images, input
	//  parameters must be adjusted to reflect the characteristics of the
	//  data. The output of the filter is an Mesh.  Users can use
	//  their own visualization packages to see the segmentation results.
	//
	//  Software Guide : EndLatex

	return EXIT_SUCCESS;
}
