#ifndef __kmSegmentation_h
#define __kmSegmentation_h

#include "itkImage.h"
#include "itkRegularSphereMeshSource.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkSobelEdgeDetectionImageFilter.h"

#include "itkDeformableSimplexMesh3DICPForceFilter.h"
#include "itkDeformableSimplexMesh3DGradientConstraintForceFilter.h"
#include "itkDeformableSimplexMesh3DBalloonForceFilter.h"

#include "kmUtility.h"

namespace km 
{
	template<class InputImageType, class OutputImageType>
	typename OutputImageType::Pointer
		sobelEdge(typename InputImageType* inputImage)
	{
		typedef itk::SobelEdgeDetectionImageFilter<InputImageType,OutputImageType>   EdgeFilterType;

		EdgeFilterType::Pointer edgeFilter = EdgeFilterType::New();
		edgeFilter->SetInput( inputImage );
		edgeFilter->Update();

		return edgeFilter->GetOutput();
	}

	template<class MeshType>
	typename MeshType::Pointer
		createRegularSphereMesh( 
														typename MeshType::PointType center, 
														double scale, 
														unsigned int resolution )
	{
		typedef itk::RegularSphereMeshSource<MeshType>  SphereMeshSourceType;
		typedef SphereMeshSourceType::PointType PointType;
		typedef SphereMeshSourceType::VectorType VectorType;

		SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();
		VectorType scales;
		scales.Fill( scale );

		mySphereMeshSource->SetCenter(center);
		mySphereMeshSource->SetResolution( resolution );
		mySphereMeshSource->SetScale(scales);
		mySphereMeshSource->Update();

		return mySphereMeshSource->GetOutput();
	}

	template<class MeshType>
	void
		deformSegSimplexMeshByICP(  
																		typename MeshType::Pointer & targetMesh,
																		typename MeshType::Pointer & deformedMesh,
																		double alpha = 0.1,
																		double kappa = 1.0,
																		unsigned int iterations = 300,
																		unsigned int rigidity = 2)
	{
		typedef itk::DeformableSimplexMesh3DICPForceFilter<MeshType,MeshType>                 DeformFilterType;
		DeformFilterType::Pointer deformFilter = DeformFilterType::New();
		deformFilter->SetInput( deformedMesh );
		deformFilter->SetTargetMesh( targetMesh );
		deformFilter->SetAlpha( alpha );
		deformFilter->SetGamma( 0.01 );
		deformFilter->SetKappa( kappa );
		deformFilter->SetIterations( iterations );
		deformFilter->SetRigidity( rigidity );

		KM_DEBUG_INFO( "Start to deforming..." );

		itk::SimpleFilterWatcher watcher(deformFilter, "DeformableSimplexMesh3DICPForceFilter");

		deformFilter->Update();

		KM_DEBUG_INFO( "Now stop deformation fitting.." );
	}

	template<class ImageType, class MeshType>
	void
		deformSegSimplexMeshByGradient(  
		typename ImageType::Pointer & inputImage,
		typename MeshType::Pointer & deformedMesh,
		double gradientSigma = 1.0,
		double alpha = 0.2,
		double beta = 0.1,
		unsigned int iterations = 50,
		unsigned int rigidity = 1)
	{
		typedef itk::DeformableSimplexMesh3DBalloonForceFilter<MeshType,MeshType>             DeformFilterType;
		typedef DeformFilterType::GradientImageType                                           GradientImageType;
		typedef itk::Image<float, ImageType::ImageDimension>                                  EdgeImaegType;
		typedef itk::GradientRecursiveGaussianImageFilter<EdgeImaegType,GradientImageType>    GradientFilterType;

		EdgeImaegType::Pointer edgeMap = km::calculateGradientMagnitudeImage<EdgeImaegType, EdgeImaegType>( inputImage, gradientSigma );
		edgeMap = km::rescaleIntensity<EdgeImaegType, EdgeImaegType>( edgeMap, 0, 30.0 );

		//km::writeImage<EdgeImaegType>( "edgeMap.nii.gz", edgeMap );

		GradientImageType::Pointer  gradientMap = km::calculateRecursiveGradientImage<EdgeImaegType, GradientImageType>( edgeMap, gradientSigma );
		KM_DEBUG_INFO( "Generate gradient image done!" );

		//km::writeImage<GradientImageType>( "gradientMap.mha", gradientMap );

		DeformFilterType::Pointer deformFilter = DeformFilterType::New();
		//deformFilter->SetImage( inputImage );
		deformFilter->SetInput( deformedMesh );
		deformFilter->SetGradient( gradientMap );
		deformFilter->SetAlpha( alpha );
		deformFilter->SetBeta( beta );
		deformFilter->SetKappa( 0.2 );
		//deformFilter->SetGamma( 0.01 );
		//deformFilter->SetRange( 1 );
		deformFilter->SetIterations( iterations );
		deformFilter->SetRigidity( rigidity );

		KM_DEBUG_INFO( "Start to deforming..." );

		itk::SimpleFilterWatcher watcher(deformFilter, "DeformableSimplexMesh3DBalloonForceFilter");

		try
		{
			deformFilter->Update();
		}
		catch(itk::ExceptionObject e)
		{
			std::cout<<e<<std::endl;
		}

		KM_DEBUG_INFO( "Now stop deformation fitting.." );
	}

} //End of namespace km


#endif