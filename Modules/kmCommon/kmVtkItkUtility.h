#ifndef __kmVtkItkUtility_h
#define __kmVtkItkUtility_h

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include "vtkPolyData.h""
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkDecimatePro.h"

#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "itkImage.h"
#include "itkSmartPointer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "itkVTKPolyDataWriter.h"
#include "itkVTKPolyDataReader.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkPadImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkCastImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMedianImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkTriangleMeshToSimplexMeshFilter.h"
#include "itkSimplexMeshToTriangleMeshFilter2.h"
#include "itkSimplexMesh.h"
#include "itkWarpMeshFilter.h"
#include "itkTransformMeshFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"
#include "itkMinMaxCurvatureFlowImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientVectorFlowImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkShiftScaleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageDuplicator.h"
#include "itkWatershedImageFilter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkNormalizedCorrelationImageFilter.h"
#include "itkSimplexMeshAdaptTopologyFilter.h"
#include "itkSimplexMesh.h"
#include "itkIdentityTransform.h"
#include "itkTransformMeshFilter.h"
#include "itkAdaptSimplexMesh3DFilter.h"

#define kmStaticImageMacro(T) typedef typename T::PixelType PixelType;         \
	typedef typename T::IndexType IndexType;         \
	typedef typename T::PointType PointType;         \
	typedef typename T::SizeType  SizeType;          \
	typedef typename T::SpacingType SpacingType;     \
	typedef typename T::RegionType RegionType;       \
	typedef typename T::DirectionType DirectionType; \

using namespace itk;
using namespace std;

namespace km
{
	//-------------------------------------------------------//
	//              ITK image read/write                     //
	//-------------------------------------------------------//
	//Read ITK Image
	template<typename TImageType>
	typename TImageType::Pointer
		readImage(const char* filename)
	{
		typedef itk::ImageFileReader<TImageType> ImageFileReaderType;
		typedef typename ImageFileReaderType::Pointer ImageFileReaderPointer;
		ImageFileReaderPointer reader = ImageFileReaderType::New();
		reader->SetFileName( filename );
		try{
			reader->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL; // Since the goal of the example is to catch the exception, we declare this a success.
		}
		return reader->GetOutput();
	}
	
	//Read ITK Image
	template<typename TImageType>
	typename TImageType::Pointer
		readImage(const std::string filename)
	{
		return readImage<TImageType>(filename.c_str());
	}

	//Write ITK Image
	template<typename TImageType>
	void
		writeImage(const char* filename, const typename TImageType* image)
	{
		typedef itk::ImageFileWriter<TImageType> ImageFileWriterType;
		typedef typename ImageFileWriterType::Pointer ImageFileWriterPointer;
		ImageFileWriterPointer writer = ImageFileWriterType::New();
		writer->SetFileName( filename );
		writer->SetInput( image );
		try{
			writer->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}
	}

	//Write ITK Image
	template<typename TImageType>
	void writeImage(const std::string filename, const typename TImageType* image)
	{
		writeImage<TImageType>(filename.c_str(), image);
	}
	
	////Write ITK Image
	template<typename TImageType>
	void writeImage(
		const std::string filenameprefix, 
		int index, 
		const char* extension, 
		const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;
		km::writeImage<TImageType>( ss.str().c_str(), image );
	}
	

	template<typename TImageType>
	void
		writeImage(const char* dirname, 
		const char* filename, 
		const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filename;
		km::writeImage<TImageType>( ss.str().c_str(), image );
	}

	template<typename TImageType>
	void
		writeImage(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const typename TImageType* image)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;
		km::writeImage<TImageType>( ss.str().c_str(), image );
	}

	//-------------------------------------------------------//
	//                 ITK mesh read/write                   //
	//-------------------------------------------------------//	
	//Read Mesh
	template<typename TMeshType>
	typename TMeshType::Pointer
		readMesh( const char* filename )
	{
		//typedef itk::VTKPolyDataReader<TMeshType>        VTKPolyDataReaderType;
		//VTKPolyDataReaderType::Pointer reader = VTKPolyDataReaderType::New();
		typedef itk::MeshFileReader<TMeshType> MeshReaderType;
		MeshReaderType::Pointer reader = MeshReaderType::New();
		reader->SetFileName( filename );
		try{
			reader->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL;
		}
		return reader->GetOutput();
	}
	
	//Read Mesh
	template<typename TMeshType>
	typename TMeshType::Pointer
		readMesh( const std::string filename )
	{
		return readMesh<TMeshType>(filename.c_str());
	}

	//Write Mesh
	template<typename TMeshType>
	void
		writeMesh( const char* filename, const typename TMeshType* mesh )
	{
		typedef itk::MeshFileWriter<TMeshType> MeshWriterType;
		MeshWriterType::Pointer writer = MeshWriterType::New();
		writer->SetInput( mesh );
		writer->SetFileName( filename );
		try{
			writer->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
		}catch(...){
			std::cerr << "Shit happens!"<<std::endl;
		}
	}

	//Write Mesh
	template<typename TMeshType>
	void
		writeMesh( const std::string filename, const TMeshType* mesh )
	{
		writeMesh<TMeshType>(filename.c_str(), mesh);
	}

	template<typename TMeshType>
	void
		writeMesh(const char* dirname, 
		const std::string filename, 
		const TMeshType* mesh)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filename;
		km::writeMesh< TMeshType >( ss.str().c_str(), mesh );
	}

	template<typename TMeshType>
	void
		writeMesh(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const TMeshType* mesh)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;
		km::writeMesh< TMeshType >( ss.str().c_str(), mesh );
	}
	
	//-------------------------------------------------------//
	//             ITK transform read/write                  //
	//-------------------------------------------------------//
	template<typename TTransformType>
	void
		writeTransform(const char* dirname, 
		const std::string filenameprefix, 
		unsigned int index, 
		const char* extension, 
		const typename TTransformType* transform)
	{
		std::stringstream ss;
		ss  <<  dirname  <<  "/"  <<  filenameprefix   <<  "."  <<  index+1  <<  extension;
		km::writeTransform< TTransformType >( ss.str().c_str(), transform );
	}
	
	//-------------------------------------------------------//
	//            ITK image information calculation          //
	//-------------------------------------------------------//
	//Get image center
	template<typename TImageType>
	typename TImageType::PointType
		getCenter(const typename TImageType* inputImage)
	{
		TImageType::RegionType region = inputImage->GetLargestPossibleRegion();
		TImageType::IndexType  origin = region.GetIndex();
		TImageType::SizeType   size   = region.GetSize();
		TImageType::IndexType  centerIndex;
		for(int i=0;i<TImageType::ImageDimension;i++){
			centerIndex[i] = origin[i] + size[i]/2;
		}
		TImageType::PointType  centerPoint;
		inputImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);
		return centerPoint;
	}

	//Get image centroid
	template<typename TImageType>
	void getCentroid(
		const typename TImageType* inputImage, 
		typename TImageType::PointType &center_pt,
		typename TImageType::IndexType &center_idx)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		center_idx.Fill(0.0);
		int n = 0;
		RegionType requestedRegion = inputImage->GetLargestPossibleRegion();
		itk::ImageRegionConstIterator<TImageType> it(inputImage, requestedRegion);
		int count = 0;
		while( !it.IsAtEnd() ){
			PixelType val = it.Get();
			IndexType ind = it.GetIndex();
			if( val != 0 ){
				n++;
				for(int d = 0; d < Dimension; d++){
					center_idx[d] += ind[d];
				}
			}
			++it;
			++count;
		}
		for(int d = 0; d < Dimension; d++){
			center_idx[d] /= n;
		}
		inputImage->TransformIndexToPhysicalPoint(center_idx, center_pt); 
	}

	//Get image centroid
	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage)
	{
		TImageType::PointType center_pt;
		TImageType::IndexType center_idx;
		getCentroid<TImageType>( inputImage, center_pt, center_idx );
		return center_pt;
	}

	//Get image centroid
	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage, double lowThreshold, double highThreshold)
	{
		return getCentroid(inputImage, inputImage->GetLargestPossibleRegion(), lowThreshold, highThreshold);
	}

	//Get image centroid
	template<typename TImageType>
	typename TImageType::PointType
		getCentroid(const typename TImageType* inputImage, const typename TImageType::RegionType& searchRegion, double lowThreshold, double highThreshold)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);

		IndexType center_idx;
		center_idx.Fill(0.0);
		int n = 0;
		itk::ImageRegionConstIterator<TImageType> it(inputImage, searchRegion);
		int count = 0;
		while( !it.IsAtEnd() ){
			PixelType val = it.Get();
			IndexType ind = it.GetIndex();
			if( val>=lowThreshold && val<=highThreshold ){
				n++;
				for(int d = 0; d < Dimension; d++){
					center_idx[d] += ind[d];
				}
			}
			++it;
			++count;
		}
		for(int d = 0; d < Dimension; d++){
			center_idx[d] /= n;
		}
		PointType center_pt;
		inputImage->TransformIndexToPhysicalPoint(center_idx, center_pt); 
		return center_pt;
	}

	template<typename TImageType>
	int countPixels(const typename TImageType* inputImage, double lowThreshold, double highThreshold)
	{
		return countPixels<TImageType>( inputImage, inputImage->GetLargestPossibleRegion(), lowThreshold, highThreshold );
	}

	template<typename TImageType>
	int countPixels(const typename TImageType* inputImage, const typename TImageType::RegionType & searchRegion, double lowThreshold, double highThreshold)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		TImageType::RegionType croppedRegion(searchRegion);
		croppedRegion.Crop( inputImage->GetLargestPossibleRegion() );
		itk::ImageRegionConstIterator<TImageType> it(inputImage, croppedRegion);
		int count = 0;
		while( !it.IsAtEnd() ){
			PixelType val = it.Get();	
			if (val>=lowThreshold && val<=highThreshold){
				count++;
			}
			++it;
		}
		return count;
	}
	
	template<typename ImageType>
	typename ImageType::Pointer
		intensityWindow( typename ImageType* inputImage, double windowsmin, double windowmax, double outputmin, double outputmax )
	{
		typedef itk::IntensityWindowingImageFilter<ImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetWindowMinimum(windowsmin);
		filter->SetWindowMaximum(windowmax);
		filter->SetOutputMinimum(outputmin);
		filter->SetOutputMaximum(outputmax);
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TImageType>
	void calculateMinAndMax(const typename TImageType* inputImage, double& minimumValue, double& maximumValue)
	{
		typedef itk::MinimumMaximumImageCalculator<TImageType> MinimumMaximumImageCalculatorType;
		MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
		calculator->SetImage( inputImage );
		calculator->Compute();
		minimumValue = calculator->GetMinimum();
		maximumValue = calculator->GetMaximum();
	}

	template<typename TImageType>
	void calculateMinAndMax(
		const typename TImageType* inputImage, 
		double& minimumValue, 
		double& maximumValue,
		typename TImageType::IndexType & minIndex, 
		typename TImageType::IndexType & maxIndex)
	{
		typedef itk::MinimumMaximumImageCalculator<TImageType> MinimumMaximumImageCalculatorType;
		MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
		calculator->SetImage( inputImage );
		calculator->Compute();
		minimumValue = calculator->GetMinimum();
		maximumValue = calculator->GetMaximum();
		minIndex = calculator->GetIndexOfMinimum();
		maxIndex = calculator->GetIndexOfMaximum();
	}
	
	//Get upper/lower index of bound box
	template<typename TImageType>
	void getBoundSliceIndex(
		const typename TImageType* image, 
		int & lowestSliceIndex,
		int & highestSliceIndex,
		double backgroundValue = 0)
	{
		typedef itk::ImageSliceConstIteratorWithIndex< TImageType > ConstIteratorType;
		ConstIteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		it.SetFirstDirection( 0 );  // 0=x, 1=y, 2=z
		it.SetSecondDirection( 1 ); // 0=x, 1=y, 2=z
		lowestSliceIndex = 9999;
		highestSliceIndex = -9999;
		while( !it.IsAtEnd() ){
			while( !it.IsAtEndOfSlice() ){
				while( !it.IsAtEndOfLine() ){
					TImageType::IndexType ind = it.GetIndex();
					TImageType::PixelType val = it.Get();
					if ( val != backgroundValue ){
						if( ind[2]>highestSliceIndex ){
							highestSliceIndex = ind[2];
						}
						if( ind[2]<lowestSliceIndex ){
							lowestSliceIndex = ind[2];
						}
						it.NextSlice();
					}else{
						++it;
					}
				}
				it.NextLine();
			}
			it.NextSlice();
		}
	}

	//Get minimum bound box
	template<typename TImageType>
	void getBoundRegion(
		const typename TImageType* image, 
		typename TImageType::RegionType& region, 
		typename TImageType::PixelType backgroundValue = 0,
		double padding = 0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::ImageRegionConstIteratorWithIndex<TImageType> IteratorType;
		IteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		IndexType startIndex;
		startIndex.Fill(9999);
		IndexType endIndex;
		endIndex.Fill(0);
		while(!it.IsAtEnd()){
			PixelType val = it.Get();
			if(val !=backgroundValue ){
				IndexType idx = it.GetIndex();
				for(unsigned int i=0;i<Dimension;i++){
					if(idx[i]<startIndex[i]){
						startIndex[i] = idx[i];
					}		
					if(idx[i]>endIndex[i]){
						endIndex[i] = idx[i];
					}
				}
			}
			it++;
		}
		SpacingType spacing = image->GetSpacing();
		SizeType size;
		for(unsigned int i=0;i<Dimension;i++)
		{
			startIndex[i] -= static_cast<int>(padding/spacing[i]);
			endIndex[i] += static_cast<int>(padding/spacing[i]);
		}
		region.SetIndex( startIndex );
		region.SetUpperIndex( endIndex );
		region.Crop( image->GetLargestPossibleRegion() );
	}

	template<typename TMeshType, typename TImageType>
	void getBoundRegion(
		const typename TMeshType* mesh,
		const typename TImageType* referenceimage, 
		typename TImageType::RegionType& region,
		double padding = 0.0)
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef TMeshType::PointsContainerConstPointer     PointsContainer;
		typedef TMeshType::PointsContainerConstIterator    PointsIterator;
		TMeshType::PointType minPoint;
		minPoint.Fill(9999);
		TMeshType::PointType maxPoint;
		maxPoint.Fill(-9999);
		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();
		while ( pointItr != pointEnd ){
			TMeshType::PointType pt = pointItr.Value();
			for(int i=0;i<Dimension;i++){
				if(pt[i]<minPoint[i]){
					minPoint[i] = pt[i];
				}
				if(pt[i]>maxPoint[i]){
					maxPoint[i] = pt[i];
				}
			}
			++pointItr;
		}
		TImageType::PointType startPoint;
		TImageType::PointType endPoint;
		TImageType::IndexType startIndex;
		TImageType::IndexType endIndex;
		TImageType::SizeType  boundSize;
		for (int i=0;i<Dimension;i++){
			startPoint[i] = (minPoint[i]-padding);
			endPoint[i]   = (maxPoint[i]+padding);
		}
		referenceimage->TransformPhysicalPointToIndex( startPoint, startIndex );
		referenceimage->TransformPhysicalPointToIndex( endPoint, endIndex );
		region.SetIndex( startIndex );
		region.SetUpperIndex( endIndex );
	}
	
	//-------------------------------------------------------//
	//            ITK mesh information calculation          //
	//-------------------------------------------------------//
	//Get centroid of itk mesh
	template<class MeshType>
	typename MeshType::PointType
		getMeshCentroid(const MeshType* pointset)
	{
		typedef MeshType::PointType PointType;
		PointType centroid;
		centroid.Fill(0);
		unsigned int Dimension = MeshType::PointDimension;
		unsigned long numberOfPoints = pointset->GetNumberOfPoints();
		typedef MeshType::PointsContainer::ConstIterator PointsIterator;
		PointsIterator pointItr = pointset->GetPoints()->Begin();
		PointsIterator pointEnd = pointset->GetPoints()->End();
		while(pointItr!=pointEnd){
			PointType pt = pointItr.Value();
			for(unsigned int i=0;i<Dimension; i++){
				centroid[i] += pt[i];
			}
			pointItr++;
		}
		for(unsigned int i=0;i<Dimension; i++){
			centroid[i] /= numberOfPoints;
		}
		return centroid;
	}
	
	//-------------------------------------------------------//
	//            ITK binary image <==> ITK mesh             //
	//-------------------------------------------------------//
	//Generate itk mesh from binary image
	template<typename TImageType, typename TMeshType>
	typename TMeshType::Pointer
		generateMeshFromBinary( const typename TImageType* inputImage )
	{
		typedef itk::BinaryMask3DMeshSource<TImageType, TMeshType> BinaryMask3DMeshSourceType;
		typedef typename BinaryMask3DMeshSourceType::Pointer       BinaryMask3DMeshSourcePointer;
		BinaryMask3DMeshSourcePointer meshSource = BinaryMask3DMeshSourceType::New();
		meshSource->SetInput( inputImage );
		meshSource->SetObjectValue( 1 );
		meshSource->SetRegionOfInterest( inputImage->GetLargestPossibleRegion() );
		meshSource->Update();
		return meshSource->GetOutput();
	}

	//Generate binary image from itk mesh
	template<typename TMeshType, typename TImageType, typename TReferenceImageType>
	typename TImageType::Pointer
		generateBinaryFromMesh( 
		const typename TMeshType* inputMesh, 
		const typename TReferenceImageType* referenceImage )
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::TriangleMeshToBinaryImageFilter<TMeshType, TImageType> TriangleMeshToBinaryImageFilterType;
		typedef typename TriangleMeshToBinaryImageFilterType::Pointer       TriangleMeshToBinaryImageFilterPointer;
		TriangleMeshToBinaryImageFilterPointer filter = TriangleMeshToBinaryImageFilterType::New();
		filter->SetInput( const_cast<TMeshType*>(inputMesh) );
		filter->SetSize( referenceImage->GetLargestPossibleRegion().GetSize() );
		filter->SetSpacing( referenceImage->GetSpacing() );
		filter->SetIndex( referenceImage->GetLargestPossibleRegion().GetIndex() );
		filter->SetOrigin( referenceImage->GetOrigin() );
		filter->SetInsideValue( 1 );
		filter->SetOutsideValue( 0 );
		filter->Update();
		return filter->GetOutput();
	}

	//Generate binary image from itk mesh
	template<typename TMeshType, typename TImageType>
	typename TImageType::Pointer
		generateBinaryFromMesh( 
			const typename TMeshType* inputMesh, 
			typename TImageType::SizeType outputSize, 
			typename TImageType::SpacingType outputSpacing )
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::TriangleMeshToBinaryImageFilter<TMeshType, TImageType> TriangleMeshToBinaryImageFilterType;
		typedef typename TriangleMeshToBinaryImageFilterType::Pointer       TriangleMeshToBinaryImageFilterPointer;
		TriangleMeshToBinaryImageFilterPointer filter = TriangleMeshToBinaryImageFilterType::New();
		filter->SetInput( const_cast<TMeshType*>(inputMesh) );
		PointType origin;
		IndexType index;
		for( int i=0;i<Dimension;i++ ){
			origin[i] = 0.0;
			index[i] = 0;
		}
		filter->SetSize( outputSize );
		filter->SetSpacing( outputSpacing );
		filter->SetIndex( index );
		filter->SetOrigin( origin );
		filter->SetInsideValue( 1 );
		filter->SetOutsideValue( 0 );
		filter->Update();
		return filter->GetOutput();
	}
	
	//-------------------------------------------------------//
	//            ITK mesh transform/filtering               //
	//-------------------------------------------------------//
	//Warp/Transform a mesh
	template<typename TMeshType, typename TDisplacementType>
	typename TMeshType::Pointer
		warpMesh( const typename TMeshType* inputMesh, const typename TDisplacementType* displacement  )
	{
		typedef itk::WarpMeshFilter<TMeshType, TMeshType, TDisplacementType> WarpMeshFilterType;
		typedef typename WarpMeshFilterType::Pointer                         WarpMeshFilterPointer;
		WarpMeshFilterPointer meshWarper = WarpMeshFilterType::New();
		meshWarper->SetInput( inputMesh );
		meshWarper->SetDisplacementField( displacement );
		meshWarper->Update();
		return meshWarper->GetOutput();
	}

	template<typename TMeshType, typename TTransformType>
	typename TMeshType::Pointer
		transformMesh( const typename TMeshType* inputMesh, typename TTransformType* transform  )
	{
		typedef itk::TransformMeshFilter<TMeshType, TMeshType, TTransformType> TransformMeshFilterType;
		TransformMeshFilterType::Pointer filter = TransformMeshFilterType::New();
		filter->SetInput( inputMesh );
		filter->SetTransform( transform );
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TMeshType, typename TTransformType>
	void transformMesh( const typename TMeshType* inputMesh, typename TMeshType* outputMesh, typename TTransformType* transform  )
	{
		typedef TMeshType::PointType PointType;
		typedef TMeshType::PointsContainer::ConstPointer  PointsContainerConstPointer;
		typedef TMeshType::PointsContainer::Pointer       PointsContainerPointer;
		typedef TMeshType::PointsContainer::ConstIterator PointsContainerConstIterator;
		PointsContainerConstPointer inputPoints = inputMesh->GetPoints();
		PointsContainerPointer      outputPoints = outputMesh->GetPoints();
		if (outputPoints->Size() != inputPoints->Size()){
			outputPoints->Reserve( inputPoints->Size() );
		}
		PointsContainerConstIterator inputPointIt = inputPoints->Begin();
		PointsContainerConstIterator outputPointIt = outputPoints->Begin();
		while( (inputPointIt!=inputPoints->End()) && (outputPointIt!=outputPoints->End()) ){
			if ( inputPointIt->Index() != outputPointIt->Index() ){
				KM_DEBUG_ERROR( "Index of input mesh != index of output mesh" );
			}
			PointType p_old = inputPointIt->Value();
			PointType p_new = transform->TransformPoint( p_old );
			outputPoints->InsertElement( inputPointIt->Index(), p_new );
			inputPointIt++;
			outputPointIt++;
		}
	}
	
	//-------------------------------------------------------//
	//            Simplex mesh <==> triangle mesh            //
	//-------------------------------------------------------//
	//triangle Mesh To Simplex Mesh
	template<typename TriangleMeshType, typename SimplexMeshType>
	typename SimplexMeshType::Pointer
		triangleMeshToSimplexMesh(const typename TriangleMeshType* inputMesh)
	{
		typedef itk::TriangleMeshToSimplexMeshFilter<TriangleMeshType, SimplexMeshType> FilterType;
		typename FilterType::Pointer filter = FilterType::New();
		try{
			filter->SetInput( inputMesh );
			filter->Update();
		}catch ( itk::ExceptionObject & e  ){
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<e<<std::endl;
		}catch (...){
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
		}
		return filter->GetOutput();
	}

	//simplex Mesh To Triangle Mesh
	template<typename SimplexMeshType, typename TriangleMeshType>
	typename TriangleMeshType::Pointer
		simplexMeshToTriangleMesh(const typename SimplexMeshType* inputMesh)
	{
		typedef itk::SimplexMeshToTriangleMeshFilter<SimplexMeshType, TriangleMeshType> FilterType;
		typename FilterType::Pointer filter = FilterType::New();
		try{
			filter->SetInput( inputMesh );
			filter->Update();
		}catch ( itk::ExceptionObject & e  ){
			KM_PRINT_EXCEPTION(e);
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<e<<std::endl;
		}catch (...){
			std::cout<<"Exception occurs when converting triangle mesh to simplex mesh!"<<std::endl;
		}
		return filter->GetOutput();
	}
	
	//-------------------------------------------------------//
	//            ITK transform read/write                   //
	//-------------------------------------------------------//
	template<typename TTransform>
	void
		writeTransform(const char* filename, const typename TTransform* transform)
	{
		typedef itk::TransformFileWriter WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetInput(transform);
		writer->SetFileName(filename);
		writer->Update();
	}

	template<typename TTransform>
	typename TTransform::Pointer
		readTransform(const char* filename)
	{
		TTransform::Pointer transform = TTransform::New();
		const char* transformname = transform->GetNameOfClass();
		std::cout<<transformname<<std::endl;
		itk::TransformFactoryBase::RegisterDefaultTransforms();
		itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
		reader->SetFileName(filename);
		try{
			reader->Update();
		}catch( itk::ExceptionObject & excp ){
			std::cerr << "Error while reading the transform file" << std::endl;
			std::cerr << excp << std::endl;
			std::cerr << "[FAILED]" << std::endl;
			return transform;
		}
		typedef itk::TransformFileReader::TransformListType * TransformListType;
		TransformListType transforms = reader->GetTransformList();
		std::cout << "Number of transforms = " << transforms->size() << std::endl;
		itk::TransformFileReader::TransformListType::const_iterator it = transforms->begin();
		itk::TransformFileReader::TransformListType::const_iterator itEnd = transforms->end();
		while(it!=itEnd){
			if(!strcmp((*it)->GetNameOfClass(),transformname)){
				transform = static_cast<TTransform*>((*it).GetPointer());
			}
			it++;
		}
		return transform;
	}
	
	//-------------------------------------------------------//
	//                ITK image filtering                    //
	//-------------------------------------------------------//
	//Add images
	template<class ImageType1, class ImageType2, class OutputImageType>
	typename OutputImageType::Pointer
		addImage(const typename ImageType1* image1,
		typename ImageType2* image2)
	{
		typedef itk::AddImageFilter<ImageType1, ImageType2, OutputImageType> AddImageFilterType;
		AddImageFilterType::Pointer adder = AddImageFilterType::New();
		adder->SetInput1( image1 );
		adder->SetInput2( image2 );
		adder->Update();
		return adder->GetOutput();
	}

	//Cast image
	template<typename TImageType1, typename TImageType2>
	typename TImageType2::Pointer
		castImage( const typename TImageType1* inputImage )
	{
		typedef itk::CastImageFilter<TImageType1, TImageType2> CastImageFilterType;
		typename CastImageFilterType::Pointer caster = CastImageFilterType::New();
		caster->SetInput( inputImage );
		caster->Update();
		return caster->GetOutput();
	}

	//Resample Image
	template<typename TImageType>
	typename TImageType::Pointer
		resampleImage( 
		const typename TImageType* inputImage, 
		const typename TImageType::SpacingType spacing, 
		double defaultValue = 0,
		int interpolatorType = 0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		SizeType inputSize = inputImage->GetLargestPossibleRegion().GetSize();
		SpacingType inputSpacing = inputImage->GetSpacing();
		itk::Vector<double> inputLenth;
		for(int i = 0;i<Dimension; i++){
			inputLenth[i] = ((float)inputSize[i])*inputSpacing[i];
		}
		SizeType outputSize;
		SpacingType outputSpacing;
		for(int i = 0;i<Dimension; i++){
			outputSize[i] = inputLenth[i]/spacing[i];
			outputSpacing[i] = spacing[i];
		}
		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer        ResampleImageFilterPointer;
		ResampleImageFilterPointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputParametersFromImage( inputImage );
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetSize( outputSize );
		resampler->SetOutputSpacing( outputSpacing );
		resampler->SetDefaultPixelValue( defaultValue );

		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			{
				LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
				resampler->SetInterpolator( interpolator );
				break;
			}
			
		case 1:
			{
				NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
				resampler->SetInterpolator( interpolator2 );
				break;
			}
			
		default:
			{
				std::cout<<"Unknown interpolator type!"<<std::endl;
				break;
			}
			
		}
		resampler->Update();
		return resampler->GetOutput();
	}

	//Resample Image
	template<typename TImageType>
	typename TImageType::Pointer
		resampleImage( const typename TImageType* inputImage, const typename TImageType::SpacingType spacing, const typename TImageType::SizeType size, double defaultValue )
	{
		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer        ResampleImageFilterPointer;
		ResampleImageFilterPointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputParametersFromImage( inputImage );
		resampler->SetOutputOrigin( inputImage->GetOrigin() );
		resampler->SetSize( size );
		resampler->SetOutputSpacing( spacing );
		resampler->SetDefaultPixelValue( defaultValue );
		resampler->Update();

		return resampler->GetOutput();
	}

	//Padding Image
	template<typename TImageType>
	typename TImageType::Pointer
		padImage( const typename TImageType* inputImage, typename TImageType::SizeType lowerBound, typename TImageType::SizeType upperBound, double defaultValue=0)
	{
		kmStaticImageMacro(TImageType);
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::PadImageFilter<TImageType, TImageType> PadImageFilterType;
		typedef typename PadImageFilterType::Pointer        PadImageFilterPointer;
		PadImageFilterPointer filter = PadImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetPadLowerBound( lowerBound );
		filter->SetPadUpperBound( upperBound );
		itk::ConstantBoundaryCondition< TImageType > bc;
		bc.SetConstant( defaultValue );
		filter->SetBoundaryCondition( &bc );
		filter->Update();
		return filter->GetOutput();
	}

	//Translate Image
	template<typename TImageType>
	typename TImageType::Pointer
		translateImage( const typename TImageType* inputImage, itk::Vector<double> & offsetVector, const typename TImageType* referenceImage )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::TranslationTransform<double, Dimension> TranslationTransformType;
		typedef typename TranslationTransformType::Pointer TranslationTransformPointer;
		TranslationTransformPointer translationTransform = TranslationTransformType::New();
		translationTransform->Translate(offsetVector);
		typedef itk::ResampleImageFilter<TImageType, TImageType> ResampleImageFilterType;
		typedef typename ResampleImageFilterType::Pointer                 ResampleImageFilterPointer;
		ResampleImageFilterPointer resampleFilter = ResampleImageFilterType::New();
		resampleFilter->SetTransform(translationTransform.GetPointer());
		resampleFilter->SetInput(inputImage);
		resampleFilter->SetOutputParametersFromImage( referenceImage );
		resampleFilter->Update();
		return resampleFilter->GetOutput();
	}

	//Binary Close Image
	template<typename TImageType>
	typename TImageType::Pointer
		binaryClose( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryMorphologicalClosingImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> BinaryMorphologicalClosingImageFilterType;
		typedef typename BinaryMorphologicalClosingImageFilterType::Pointer BinaryMorphologicalClosingImageFilterPointer;
		BinaryMorphologicalClosingImageFilterPointer filter = BinaryMorphologicalClosingImageFilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++){
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();
		return filter->GetOutput();
	}

	//Binary Open Image
	template<typename TImageType>
	typename TImageType::Pointer
		binaryOpen( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryMorphologicalOpeningImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> BinaryMorphologicalClosingImageFilterType;
		typedef typename BinaryMorphologicalClosingImageFilterType::Pointer BinaryMorphologicalClosingImageFilterPointer;
		BinaryMorphologicalClosingImageFilterPointer filter = BinaryMorphologicalClosingImageFilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++){
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();
		return filter->GetOutput();
	}

	//Binary Dilate Image
	template<typename TImageType>
	typename TImageType::Pointer
		binaryDilate( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		typedef itk::BinaryBallStructuringElement<unsigned char, 3> BinaryBallStructuringElementType;
		typedef itk::BinaryDilateImageFilter<TImageType, TImageType, BinaryBallStructuringElementType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		BinaryBallStructuringElementType structuringElement;
		BinaryBallStructuringElementType::RadiusType radius;
		for(int i=0;i<Dimension;i++){
			radius[i] = kernelRadius[i];
		}
		structuringElement.SetRadius(radius);
		structuringElement.CreateStructuringElement();
		filter->SetInput( inputImage );
		filter->SetKernel( structuringElement );
		filter->SetForegroundValue( 1 );
		filter->Update();
		return filter->GetOutput();
	}

	//Binary median smooth image
	template<typename TImageType>
	typename TImageType::Pointer
		binaryMedianSmooth( const typename TImageType* inputImage, typename TImageType::SizeType kernelRadius )
	{
		typedef itk::BinaryMedianImageFilter<TImageType, TImageType> BinaryMedianImageFilterType;
		BinaryMedianImageFilterType::Pointer filter = BinaryMedianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetForegroundValue( 1 );
		filter->SetBackgroundValue( 0 );
		filter->SetRadius( kernelRadius );
		filter->Update();
		return filter->GetOutput();
	}

	//Calculate Mauer distance map of binary image
	template<typename TImageType, typename TDistanceMapType>
	typename TDistanceMapType::Pointer
		calculateDistanceMap( const typename TImageType* inputImage )
	{
		typedef itk::SignedMaurerDistanceMapImageFilter<TImageType, TDistanceMapType> SignedMaurerDistanceMapImageFilterType;
		typedef typename SignedMaurerDistanceMapImageFilterType::Pointer              SignedMaurerDistanceMapImageFilterPointer;
		SignedMaurerDistanceMapImageFilterPointer distanceFilter = SignedMaurerDistanceMapImageFilterType::New();
		distanceFilter->SetInput( inputImage );
		distanceFilter->Update();
		return distanceFilter->GetOutput();
	}

	template<typename TImageType, typename TGradientImageType>
	typename TGradientImageType::Pointer
		calculateRecursiveGradientImage( const typename TImageType* inputImage, double sigma )
	{
		typedef itk::GradientRecursiveGaussianImageFilter<TImageType, TGradientImageType> GradientRecursiveGaussianImageFilterType;
		GradientRecursiveGaussianImageFilterType::Pointer filter = GradientRecursiveGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetSigma( sigma );
		filter->SetNormalizeAcrossScale( true );
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TImageType, typename TGradientImageType>
	typename TGradientImageType::Pointer
		calculateGradientMagnitudeImage( const typename TImageType* inputImage, double sigma )
	{
		typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<TImageType, TGradientImageType> GradientMagnitudeRecursiveGaussianImageFilterType;
		GradientMagnitudeRecursiveGaussianImageFilterType::Pointer filter = GradientMagnitudeRecursiveGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetSigma( sigma );
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TGradientImageType>
	typename TGradientImageType::Pointer
		calcuateGVF(const typename TGradientImageType* inputGradient, int iterations, double noiseLevel)
	{
		typedef itk::GradientVectorFlowImageFilter< TGradientImageType, TGradientImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputGradient );
		filter->SetIterationNum( iterations );
		filter->SetNoiseLevel( noiseLevel );
		filter->Update();
		return filter->GetOutput();
	}



	//Warp/Transform a Image
	template<typename TImageType, typename TDisplacementType>
	typename TImageType::Pointer
		warpImage( const typename TImageType* inputImage, const typename TDisplacementType* displacement  )
	{
		typedef itk::WarpImageFilter<TImageType, TImageType, TDisplacementType> FilterType;
		FilterType::Pointer meshWarper = FilterType::New();
		meshWarper->SetInput( inputImage );
		meshWarper->SetDisplacementField( displacement );
		meshWarper->SetOutputSpacing( inputImage->GetSpacing() );
		meshWarper->SetOutputOrigin(  inputImage->GetOrigin() );
		meshWarper->Update();

		return meshWarper->GetOutput();
	}

	template<typename TImageType, typename TTransformType>
	typename TImageType::Pointer
		transformImage( const typename TImageType* inputImage, const typename TTransformType* transform  )
	{
		typedef itk::ResampleImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetTransform( transform );
		filter->SetSize(    inputImage->GetLargestPossibleRegion().GetSize() );
		filter->SetOutputOrigin(  inputImage->GetOrigin() );
		filter->SetOutputSpacing( inputImage->GetSpacing() );
		filter->SetOutputDirection( inputImage->GetDirection() );
		filter->Update();
		return filter->GetOutput();
	}

	//Theshold Image
	template<typename TImageType>
	typename TImageType::Pointer
		thresholdImage( const TImageType* inputImage, double lowValue, double highValue, double outsideValue = 0 )
	{
		typedef itk::ThresholdImageFilter<TImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->ThresholdOutside( static_cast<typename TImageType::PixelType>(lowValue), static_cast<typename TImageType::PixelType>(highValue) );
		thresholdFilter->SetOutsideValue( outsideValue );
		thresholdFilter->Update();
		return thresholdFilter->GetOutput();
	}

	//Theshold Image
	template<typename TImageType>
	typename TImageType::Pointer
		thresholdAboveImage( const TImageType* inputImage, double highValue, double outsideValue = highValue )
	{
		typedef itk::ThresholdImageFilter<TImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->ThresholdAbove ( static_cast<typename TImageType::PixelType>(highValue) );
		thresholdFilter->SetOutsideValue( outsideValue );
		thresholdFilter->Update();
		return thresholdFilter->GetOutput();
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TBinaryImageType::Pointer
		binaryThresholdImage( const TImageType* inputImage, double lowValue, double highValue, double insideValue, double outsideValue )
	{
		typedef itk::BinaryThresholdImageFilter<TImageType, TBinaryImageType> ThresholdImageFilterType;
		typename ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput( inputImage );
		thresholdFilter->SetLowerThreshold(lowValue);
		thresholdFilter->SetUpperThreshold(highValue);
		thresholdFilter->SetInsideValue(insideValue);
		thresholdFilter->SetOutsideValue(outsideValue);
		thresholdFilter->Update();
		return thresholdFilter->GetOutput();
	}
	
	template<class InputImageType, class OutputImageType>
	typename OutputImageType::Pointer
		binaryThresholdImage(const typename InputImageType* inputImage, const std::vector<std::pair<double, double>> & thresholds)
	{
		OutputImageType::Pointer mask = OutputImageType::New();
		mask->SetRegions( inputImage->GetLargestPossibleRegion() );
		mask->SetSpacing( inputImage->GetSpacing() );
		mask->SetOrigin( inputImage->GetOrigin() );
		mask->SetDirection( inputImage->GetDirection() );
		mask->Allocate();
		mask->FillBuffer(0);
		typedef itk::OrImageFilter <OutputImageType> OrImageFilterType;
		for(int i=0;i<thresholds.size();i++){
			OutputImageType::Pointer masktmp = km::binaryThresholdImage<InputImageType, OutputImageType>(inputImage, thresholds[i].first, thresholds[i].second, 1, 0);
			OrImageFilterType::Pointer orFilter = OrImageFilterType::New();
			orFilter->SetInput(0, mask);
			orFilter->SetInput(1, masktmp);
			orFilter->Update();
			mask = orFilter->GetOutput();
		}
		return mask;
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TBinaryImageType::Pointer
		binaryContourImage( const TImageType* inputImage, double foregroundValue, double backgourndValue )
	{
		typedef itk::BinaryContourImageFilter <TImageType, TBinaryImageType> BinaryContourImageFilterType;
		typename BinaryContourImageFilterType::Pointer filter = BinaryContourImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetForegroundValue(foregroundValue);
		filter->SetBackgroundValue(backgourndValue);
		filter->FullyConnectedOn();
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		gaussSmooth(const typename TImageType* inputImage, double sigma)
	{
		typedef itk::DiscreteGaussianImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter=FilterType::New();
		filter->SetInput(inputImage);
		filter->SetVariance(sigma);
		filter->Update();
		return filter->GetOutput();
	}


	template<typename TImageType>
	typename TImageType::Pointer
		minMaxSmooth(const typename TImageType* inputImage, int numberOfIterations, double timeStep, double radius)
	{
		typedef itk::Image<float, TImageType::ImageDimension> FloatImageType;
		FloatImageType::Pointer image = km::castImage<TImageType, FloatImageType>( inputImage );
		typedef itk::MinMaxCurvatureFlowImageFilter<FloatImageType, FloatImageType> MinMaxCurvatureFlowImageFilterType;
		typedef MinMaxCurvatureFlowImageFilterType::Pointer MinMaxCurvatureFlowImageFilterPointer;
		MinMaxCurvatureFlowImageFilterPointer minMaxFilter = MinMaxCurvatureFlowImageFilterType::New();
		minMaxFilter->SetInput( image );
		minMaxFilter->SetTimeStep( timeStep );
		minMaxFilter->SetNumberOfIterations( numberOfIterations );
		minMaxFilter->SetStencilRadius( radius );
		minMaxFilter->Update();
		FloatImageType::Pointer smoothed = minMaxFilter->GetOutput();
		TImageType::Pointer smoothedImage = km::castImage<FloatImageType, TImageType>( smoothed );
		return smoothedImage;
	}

	template<typename TImageType>
	typename TImageType::Pointer
		anisotropicSmooth(const typename TImageType* inputImage, int numberOfIterations, double timeStep, double conductance)
	{
		typedef itk::GradientAnisotropicDiffusionImageFilter<  TImageType,  TImageType >  AnisotropicSmoothFilterType;
		typedef AnisotropicSmoothFilterType::Pointer AnisotropicSmoothFilterPointer;
		AnisotropicSmoothFilterPointer anisoSmoothfilter = AnisotropicSmoothFilterType::New();
		anisoSmoothfilter->SetNumberOfIterations( numberOfIterations );
		anisoSmoothfilter->SetTimeStep( timeStep );
		anisoSmoothfilter->SetConductanceParameter( conductance );
		anisoSmoothfilter->SetInput( inputImage );
		anisoSmoothfilter->UseImageSpacingOn();
		anisoSmoothfilter->Update();
		return anisoSmoothfilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		medianSmooth(const typename TImageType* inputImage, int radiusValue)
	{
		typedef itk::MedianImageFilter<TImageType, TImageType > FilterType;
		FilterType::Pointer medianFilter = FilterType::New();
		FilterType::InputSizeType radius;
		radius.Fill(radiusValue);
		medianFilter->SetRadius(radius);
		medianFilter->SetInput( inputImage );
		medianFilter->Update();
		return medianFilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		medianSmooth(const typename TImageType* inputImage, typename TImageType::SizeType radius)
	{
		typedef itk::MedianImageFilter<TImageType, TImageType > FilterType;
		FilterType::Pointer medianFilter = FilterType::New();
		medianFilter->SetRadius(radius);
		medianFilter->SetInput( inputImage );
		medianFilter->Update();
		return medianFilter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		discreteGaussSmooth(const typename TImageType* inputImage, double variance)
	{
		typedef itk::DiscreteGaussianImageFilter<TImageType,TImageType> DiscreteGaussianImageFilterType;
		DiscreteGaussianImageFilterType::Pointer filter = DiscreteGaussianImageFilterType::New();
		filter->SetInput( inputImage );
		filter->SetVariance( variance );
		filter->Update();
		return filter->GetOutput();
	}

	template<class TImageType, class ReferenceImageType, class TTransformType>
	typename TImageType::Pointer
		transformImageByReference(
		const typename TImageType* inputimage, 
		const typename ReferenceImageType* referenceimage, 
		const typename TTransformType* transform,
		double defaultValue = 0,
		int   interpolatorType = 0) //0:liner interpolator 1:nearest...
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		resample->SetTransform( transform );
		resample->SetInput( inputimage );
		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			{
				LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
				resample->SetInterpolator( interpolator );
				break;
			}
			
		case 1:
			{
				NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
				resample->SetInterpolator( interpolator2 );
				break;
			}
			
		default:
			{
				std::cout<<"Unknown interpolator type!"<<std::endl;
				break;
			}
			
		}
		resample->Update();
		return resample->GetOutput();
	}

	template<class TImageType, class ReferenceImageType>
	typename TImageType::Pointer
		resampleImageByReference(
		const typename TImageType* inputimage, 
		const typename ReferenceImageType* referenceimage, 
		double defaultValue = 0,
		int interpolatorType = 0)
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		resample->SetInput( inputimage );
		typedef itk::LinearInterpolateImageFunction<TImageType> LinearInterpolateType;
		typedef itk::NearestNeighborInterpolateImageFunction<TImageType> NearestNeighborInterpolateType;
		switch (interpolatorType)
		{
		case 0:
			LinearInterpolateType::Pointer interpolator = LinearInterpolateType::New();
			resample->SetInterpolator( interpolator );
			break;
		case 1:
			NearestNeighborInterpolateType::Pointer interpolator2 = NearestNeighborInterpolateType::New();
			resample->SetInterpolator( interpolator2 );
			break;
		default:
			std::cout<<"Unknown interpolator type!"<<std::endl;
			break;
		}
		resample->Update();
		return resample->GetOutput();
	}

	template<class TImageType, class TTransformType>
	typename TImageType::Pointer
		transformImage(const typename TImageType* inputimage, const typename TImageType* referenceimage, const typename TTransformType* transform, double defaultValue = 0)
	{
		typedef itk::ResampleImageFilter< TImageType, TImageType,double >  ResampleFilterType;
		ResampleFilterType::Pointer resample = ResampleFilterType::New();
		resample->SetOutputOrigin( referenceimage->GetOrigin() );
		resample->SetOutputSpacing( referenceimage->GetSpacing() );
		resample->SetOutputDirection ( referenceimage->GetDirection() );
		resample->SetOutputStartIndex ( referenceimage->GetLargestPossibleRegion().GetIndex() );
		resample->SetSize ( referenceimage->GetLargestPossibleRegion().GetSize() );
		resample->SetDefaultPixelValue( defaultValue );
		resample->SetTransform( transform );
		resample->SetInput( inputimage );
		resample->Update();
		return resample->GetOutput();
	}

	template<typename TImageType, typename TBinaryImageType>
	typename TImageType::Pointer
		maskImage(const typename TImageType* inputimage, const typename TBinaryImageType* maskimage)
	{
		typedef itk::MultiplyImageFilter< TImageType, TBinaryImageType >  FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput1( inputimage );
		filter->SetInput2( maskimage );
		try{
			filter->Update();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return NULL; 
		}
		return filter->GetOutput();
	}

	template<typename ImageType1, typename ImageType2, typename OutputImageType>
	typename OutputImageType::Pointer
		multiplyImage( const ImageType1* image1, const ImageType2* image2 )
	{
		typedef itk::MultiplyImageFilter<ImageType1, ImageType2, OutputImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput1( image1 );
		filter->SetInput2( image2 );
		filter->Update();
		return filter->GetOutput();
	}

	template<typename TImageType>
	typename TImageType::Pointer
		shiftScale(const typename TImageType* inputImage, double shiftValue, double scaleValue)
	{
		typedef itk::ShiftScaleImageFilter<TImageType, TImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput( inputImage );
		filter->SetScale( scaleValue );
		filter->SetShift( shiftValue );
		filter->Update();
		return filter->GetOutput();
	}
	
	template<class ImageType>
	void
		shiftMinimum(typename ImageType::Pointer & image, double targetMin = -1024)
	{
		double minval, maxval;
		km::calculateMinAndMax<ShortImageType>( image, minval, maxval );

		image = km::shiftScale<ImageType>( image, targetMin - minval, 1.0 );
	}
	
	template<typename TImageType>
	typename TImageType::Pointer
		changeImage(
		const typename TImageType* inputImage, 
		const typename TImageType::SpacingType newSpacing,
		const typename TImageType::PointType   newOrigin)
	{
		typedef itk::ChangeInformationImageFilter<TImageType> ChangeInformationImageFilterType;
		ChangeInformationImageFilterType::Pointer changer = ChangeInformationImageFilterType::New();
		changer->SetInput( inputImage );
		changer->SetOutputSpacing( newSpacing );
		changer->ChangeSpacingOn();
		changer->SetOutputOrigin( newOrigin );
		changer->ChangeOriginOn();
		changer->Update();

		return changer->GetOutput();
	}

	template<typename InputImageType, typename OutputImageType>
	typename OutputImageType::Pointer
		rescaleIntensity(const typename InputImageType* inputImage, double outputMinimum, double outputMaximum)
	{
		typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> FilterType;
		FilterType::Pointer filter = FilterType::New();
		filter->SetInput(inputImage);
		filter->SetOutputMinimum(outputMinimum);
		filter->SetOutputMaximum(outputMaximum);
		filter->Update();
		return filter->GetOutput();
	}

	template<typename InputImageType, typename ReferenceImageType>
	typename InputImageType::Pointer
		changeByCopying(const typename InputImageType* inputImage, const typename ReferenceImageType* referenceImage)
	{
		typedef itk::ChangeInformationImageFilter<InputImageType> ChangeInformationImageFilterType;
		ChangeInformationImageFilterType::Pointer changer = ChangeInformationImageFilterType::New();
		changer->SetInput( inputImage );
		changer->SetOutputOrigin( referenceImage->GetOrigin() );
		changer->ChangeOriginOn();
		changer->Update();
		return changer->GetOutput();
	}

	template<typename InputImageType, typename ReferenceImageType>
	typename InputImageType::Pointer
		resampleByCopying(const typename InputImageType* inputImage, const typename ReferenceImageType* referenceImage)
	{
		typedef itk::ResampleImageFilter<InputImageType, InputImageType> ResampleImageFilterType;
		ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
		resampler->SetInput( inputImage );
		resampler->SetOutputOrigin ( referenceImage->GetOrigin() );
		resampler->SetOutputSpacing ( referenceImage->GetSpacing() );
		resampler->SetOutputDirection ( referenceImage->GetDirection() );
		resampler->SetOutputStartIndex ( referenceImage->GetLargestPossibleRegion().GetIndex() );
		resampler->SetSize ( referenceImage->GetLargestPossibleRegion().GetSize() );
		resampler->Update();
		return resampler->GetOutput();
	}

	template<typename TImageType>
	void boundRegions( 
		typename TImageType* image, 
		typename TImageType::RegionType region1,
		typename TImageType::RegionType region2,
		typename TImageType::RegionType regionResult)
	{
		itkStaticConstMacro(Dimension, unsigned int, TImageType::ImageDimension);
		TImageType::IndexType startIndex1 = region1.GetIndex();
		TImageType::IndexType endIndex1 = region1.GetUpperIndex();
		TImageType::IndexType startIndex2 = region2.GetIndex();
		TImageType::IndexType endIndex2 = region2.GetUpperIndex();
		for(int i=0;i<Dimension;i++){
			if(startIndex1[i]>startIndex2[i]){
				startIndex1[i] = startIndex2[i];
			}
			if(endIndex1[i]<endIndex2[i]){
				endIndex1[i] = endIndex2[i];
			}
		}
		regionResult.SetIndex( startIndex1 );
		regionResult.SetUpperIndex( endIndex1 );
	}


	template<typename TImageType>
	unsigned long calculateNumberOfPixels(
		const typename TImageType* image, 
		typename TImageType::PixelType lowerThreshold,
		typename TImageType::PixelType upperThreshold)
	{
		typedef itk::ImageRegionConstIteratorWithIndex<TImageType> IteratorType;
		IteratorType it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();

		unsigned long counts = 0;
		while(!it.IsAtEnd()){
			double val = it.Get();
			if(val>=lowerThreshold && val<=upperThreshold){
				counts++;
			}
			++it;
		}
		return counts;
	}

	template<typename ImageType>
	typename ImageType::Pointer
		extractByBoundRegion(const typename ImageType* inputImage, typename ImageType::RegionType roi)
	{
		typedef itk::RegionOfInterestImageFilter< ImageType, ImageType > ExtractFilterType;
		ExtractFilterType::Pointer extractfilter = ExtractFilterType::New();
		extractfilter->SetRegionOfInterest(roi);
		extractfilter->SetInput( inputImage );
		extractfilter->Update();
		return extractfilter->GetOutput();
	}

	template<class TImageType>
	void hightlight(typename TImageType* image, typename TImageType::PointType pt, double hvalue, int radius = 0 )
	{
		TImageType::IndexType index;
		image->TransformPhysicalPointToIndex( pt, index );
		TImageType::SizeType rad;
		rad.Fill( 2*radius + 1 );

		TImageType::IndexType idx;
		for(int i=0;i<TImageType::ImageDimension;i++){
			idx[i] = index[i] - radius;
		}

		TImageType::RegionType hRegion( idx, rad );
		hRegion.Crop( image->GetLargestPossibleRegion() );
		itk::ImageRegionIteratorWithIndex<TImageType> it( image, hRegion );
		it.GoToBegin();
		while(!it.IsAtEnd()){
			it.Set( hvalue );
			it++;
		}
	}

	template<typename ImageType>
	typename ImageType::Pointer
		cloneImage (const typename ImageType* image)
	{
		typedef itk::ImageDuplicator<ImageType> ImageDuplicatorType;
		ImageDuplicatorType::Pointer filter = ImageDuplicatorType::New();
		filter->SetInputImage( image );
		filter->Update();
		return filter->GetOutput();
	}

	template<typename ImageType, typename LabelType>
	typename LabelType::Pointer
		generateLabelImage( const ImageType* inputImage, unsigned int & numberOfComponents)
	{
		typedef itk::ConnectedComponentImageFilter<ImageType, LabelType> LabelFilterType;
		LabelFilterType::Pointer labelFilter = LabelFilterType::New();
		labelFilter->SetInput( inputImage );
		labelFilter->Update();
		numberOfComponents = labelFilter->GetObjectCount();
		return labelFilter->GetOutput();
	}

	template<typename ImageType>
	void
		thresholdOutImage(typename ImageType::Pointer & image, double thresholdLower, double thresholdUpper, double replaceValue)
	{
		itk::ImageRegionIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd()){
			ImageType::PixelType val = it.Get();
			if (val>=thresholdLower && val<=thresholdUpper){
				it.Set( replaceValue );
			}
			it++;
		}
	}

	template<typename ImageType>
	void thresholdLabels(
		typename ImageType::Pointer & image, 
		const std::vector<int> labels , 
		int insideValue, 
		int outsideValue )
	{
		if (labels.size() == 0){
			image->FillBuffer(0);
			return;
		}

		itk::ImageRegionIterator<ImageType> it( image, image->GetLargestPossibleRegion() );
		it.GoToBegin();
		while(!it.IsAtEnd()){
			ImageType::PixelType val = it.Get();
			for (int i=0;i<labels.size();i++){
				if (val == labels[i]){
					it.Set( insideValue );
				}else{
					it.Set( outsideValue );
				}
			}
			it++;
		}
	}

	template<class ImageType>
	typename ImageType::Pointer
		absImage(const typename ImageType* inputImage)
	{
		typedef itk::AbsImageFilter<ImageType, ImageType> AbsImageFilterType;
		AbsImageFilterType::Pointer filter = AbsImageFilterType::New();
		filter->SetInput( inputImage );
		filter->Update();
		return filter->GetOutput();
	}

	template<class ImageType>
	typename ImageType::Pointer
		powImage(const typename ImageType * inputImage, const double alpha)
	{
		typedef itk::PowImageFilter<ImageType> PowImageFilterType;
		PowImageFilterType::Pointer filter = PowImageFilterType::New();
		filter->SetInput1(inputImage);
		filter->SetConstant2(alpha);
		filter->Update();
		return filter->GetOutput();
	}

	template<class ImageType>
	typename ImageType::Pointer
		histogramMatch(const ImageType* inputImage, const ImageType* refImage)
	{
		typedef itk::HistogramMatchingImageFilter<ImageType,ImageType> HEFilterType;
		HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();
		IntensityEqualizeFilter->SetReferenceImage( refImage );
		IntensityEqualizeFilter->SetInput( inputImage );
		IntensityEqualizeFilter->SetNumberOfHistogramLevels( 300);
		IntensityEqualizeFilter->SetNumberOfMatchPoints( 15);
		IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
		IntensityEqualizeFilter->Update();
		return IntensityEqualizeFilter->GetOutput();
	}
	
	//-------------------------------------------------------//
	//            Simplex mesh utilities                     //
	//-------------------------------------------------------//
	typedef itk::CovariantVector<unsigned int, 3> GeometryVectorType;
	typedef itk::Image<GeometryVectorType, 1>     GeometryImageType;

	template< typename TInputMesh, typename TOutputMesh >
	void
		copyMeshToMeshGeometry(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		const unsigned int numberOfPoints = inputMesh->GetNumberOfPoints();
		typedef typename TInputMesh::GeometryMapType  InputGeometryMapType;
		typedef typename TOutputMesh::GeometryMapType OutputGeometryMapType;
		OutputGeometryMapType * outputGeometryData = outputMesh->GetGeometryData();
		const InputGeometryMapType * inputGeometryData = inputMesh->GetGeometryData();
		if (inputGeometryData)
		{
			if (outputGeometryData == NULL){
				outputGeometryData = OutputGeometryMapType::New();
			}
			outputGeometryData->Reserve( inputGeometryData->Size() );
			InputGeometryMapType::ConstIterator inputGeometryItr = inputGeometryData->Begin();
			InputGeometryMapType::ConstIterator inputGeometryEnd = inputGeometryData->End();
			OutputGeometryMapType::Iterator outputGeometryItr = outputGeometryData->Begin();
			while( inputGeometryItr!=inputGeometryEnd ){
				SimplexMeshGeometry * outputGeometryDataItem = new SimplexMeshGeometry;
				SimplexMeshGeometry * inputGeometryDataItem = inputGeometryItr.Value();
				outputGeometryDataItem->CopyFrom( *inputGeometryDataItem );
				outputGeometryItr.Value() = outputGeometryDataItem;
				++inputGeometryItr;
				++outputGeometryItr;
			}
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		copyMeshToMeshPoints(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::PointsContainer OutputPointsContainer;
		typedef typename TInputMesh::PointsContainer  InputPointsContainer;
		typename OutputPointsContainer::Pointer outputPoints = OutputPointsContainer::New();
		const InputPointsContainer *inputPoints = inputMesh->GetPoints();
		if ( inputPoints ){
			outputPoints->Reserve( inputPoints->Size() );
			typename InputPointsContainer::ConstIterator inputItr = inputPoints->Begin();
			typename InputPointsContainer::ConstIterator inputEnd = inputPoints->End();
			typename OutputPointsContainer::Iterator outputItr = outputPoints->Begin();
			while ( inputItr != inputEnd ){
				outputItr.Value() = inputItr.Value();
				++inputItr;
				++outputItr;
			}
			outputMesh->SetPoints(outputPoints);
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		copyMeshToMeshCells(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::CellsContainer  OutputCellsContainer;
		typedef typename TInputMesh::CellsContainer   InputCellsContainer;
		typedef typename TOutputMesh::CellAutoPointer CellAutoPointer;
		outputMesh->SetCellsAllocationMethod(TOutputMesh::CellsAllocatedDynamicallyCellByCell);
		typename OutputCellsContainer::Pointer outputCells = OutputCellsContainer::New();
		const InputCellsContainer *inputCells = inputMesh->GetCells();
		if ( inputCells ){
			outputCells->Reserve( inputCells->Size() );
			typename InputCellsContainer::ConstIterator inputItr = inputCells->Begin();
			typename InputCellsContainer::ConstIterator inputEnd = inputCells->End();
			typename OutputCellsContainer::Iterator outputItr = outputCells->Begin();
			CellAutoPointer clone;
			while ( inputItr != inputEnd ){
				inputItr.Value()->MakeCopy(clone);
				outputItr.Value() = clone.ReleaseOwnership();
				++inputItr;
				++outputItr;
			}
			outputMesh->SetCells(outputCells);
		}
	}

	template< typename TInputMesh, typename TOutputMesh >
	void
		copyMeshToMeshCellLinks(const typename TInputMesh *inputMesh, typename TOutputMesh *outputMesh)
	{
		typedef typename TOutputMesh::CellLinksContainer OutputCellLinksContainer;
		typedef typename TInputMesh::CellLinksContainer  InputCellLinksContainer;
		typename OutputCellLinksContainer::Pointer outputCellLinks = OutputCellLinksContainer::New();
		const InputCellLinksContainer *inputCellLinks = inputMesh->GetCellLinks();
		if ( inputCellLinks ){
			outputCellLinks->Reserve( inputCellLinks->Size() );
			typename InputCellLinksContainer::ConstIterator inputItr = inputCellLinks->Begin();
			typename InputCellLinksContainer::ConstIterator inputEnd = inputCellLinks->End();
			typename OutputCellLinksContainer::Iterator outputItr = outputCellLinks->Begin();
			while ( inputItr != inputEnd ){
				outputItr.Value() = inputItr.Value();
				++inputItr;
				++outputItr;
			}
			outputMesh->SetCellLinks(outputCellLinks);
		}
	}

	template<typename InputMeshType, typename OutputMeshType>
	void
		copyMeshToMeshPointData( typename InputMeshType* inputMesh, typename OutputMeshType* outputMesh )
	{
		unsigned int numberOfPointsInSource = inputMesh->GetNumberOfPoints();
		unsigned int numberOfPointsInTarget = outputMesh->GetNumberOfPoints();
		if (numberOfPointsInTarget!=numberOfPointsInSource){
			std::cerr<<"Number of points is different between source mesh and target mesh !"<<std::endl;
			return;
		}
		typedef typename InputMeshType::PointDataContainerPointer         InputPointDataContainerPointer;
		typedef typename InputMeshType::PointDataContainer::ConstIterator InputPointDataContainerConstIterator;
		typedef typename OutputMeshType::PointDataContainerPointer OutputPointDataContainerPointer;
		InputPointDataContainerPointer allInputPointData = inputMesh->GetPointData();
		InputPointDataContainerConstIterator inputPointDataIt = allInputPointData->Begin();
		InputPointDataContainerConstIterator inputPointDataItEnd = allInputPointData->End();
		OutputPointDataContainerPointer allOutputPointData = outputMesh->GetPointData();
		while( inputPointDataIt != inputPointDataItEnd ){
			allOutputPointData->InsertElement( inputPointDataIt.Index(), inputPointDataIt.Value() );
			inputPointDataIt++;
		}
	}

	template<class InputMeshType, class OutputMeshType>
	typename OutputMeshType::Pointer
		cloneMesh( const InputMeshType * inputmesh )
	{
		typedef itk::IdentityTransform<double, InputMeshType::PointDimension> TransformType;
		TransformType::Pointer transform = TransformType::New();
		typedef itk::TransformMeshFilter<InputMeshType, OutputMeshType, TransformType> TransformMeshFilterType;
		TransformMeshFilterType::Pointer transformFilter = TransformMeshFilterType::New();
		transformFilter->SetInput( inputmesh );
		transformFilter->SetTransform( transform );
		transformFilter->Update();
		return transformFilter->GetOutput();
	}

	template<class MeshType>
	GeometryImageType::Pointer
		generateGeoImage( const typename MeshType * mesh  )
	{
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();
		GeometryImageType::Pointer geoImage = GeometryImageType::New();
		GeometryImageType::IndexType start;
		GeometryImageType::SizeType  size;
		start.Fill( 0 );
		size.Fill( numberOfPoints );
		GeometryImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		geoImage->SetRegions( region );
		geoImage->Allocate();
		GeometryImageType::PixelType zeroVector;
		zeroVector.Fill( 0 );
		geoImage->FillBuffer( zeroVector );
		if( mesh->GetGeometryData().IsNull() ){
			std::cerr<< "Geometry data is NULL in mesh!" << std::endl;
			return geoImage;
		}
		typedef itk::ImageRegionIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int ppp = 0;
		while (!it.IsAtEnd()){
			MeshType::IndexArray neibours = mesh->GetNeighbors( ppp );
			if (neibours.Size()==3){
				GeometryImageType::PixelType pix;
				for (int k=0;k<3;k++){
					pix[k] = neibours[k];
				}
				it.Set( pix );
			}else{
				std::cout<<"No neibours!"<<std::endl;
			}
			++it;
			++ppp;
		}
		return geoImage;
	}

	template<class MeshType>
	void
		writeSimplexMeshGeometryData(const char* filename, const typename MeshType* mesh)
	{		
		GeometryImageType::Pointer geoImage = km::generateGeoImage<MeshType>( mesh );		
		km::writeImage<GeometryImageType>( filename, geoImage );
	}

	template<class MeshType>
	void
		readSimplexMeshGeometryData(const char* filename, typename MeshType* mesh)
	{
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();
		mesh->GetGeometryData()->Reserve( numberOfPoints );
		GeometryImageType::Pointer geoImage = km::readImage<GeometryImageType>( filename );
		typedef itk::ImageRegionConstIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int p = 0;
		while (!it.IsAtEnd()){
			GeometryImageType::PixelType pix = it.Get();
			SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
			data->pos = mesh->GetPoint( p );
			for (int k=0;k<3;k++){	
				data->neighborIndices[k] = pix[k];
			}
			mesh->SetGeometryData( p, data );

			++it;
			++p;
		}
		mesh->BuildCellLinks();
		ComputeGeometry<MeshType>( mesh );
	}

	template<class MeshType>
	void
		loadSimplexMeshGeometryData(const GeometryImageType* geoImage, typename MeshType* mesh)
	{
		unsigned int numberOfPoints = mesh->GetNumberOfPoints();
		mesh->GetGeometryData()->Reserve( numberOfPoints );
		typedef itk::ImageRegionConstIterator<GeometryImageType> IteratorType;
		IteratorType it( geoImage, geoImage->GetLargestPossibleRegion() );
		it.GoToBegin();
		unsigned int p = 0;
		while (!it.IsAtEnd()){
			GeometryImageType::PixelType pix = it.Get();
			SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
			data->pos = mesh->GetPoint( p );
			for (int k=0;k<3;k++){	
				data->neighborIndices[k] = pix[k];
			}
			mesh->SetGeometryData( p, data );

			++it;
			++p;
		}
		mesh->BuildCellLinks();
		ComputeGeometry<MeshType>( mesh );
	}

	template<typename MeshType>
	void
		copyPolyDataToMeshPoints(vtkPolyData* m_PolyData, typename MeshType::Pointer mesh)
	{
		const unsigned int numberOfPointsInPolydata = m_PolyData->GetNumberOfPoints();
		const unsigned int numberOfPointsInMesh = mesh->GetNumberOfPoints();
		if (numberOfPointsInPolydata!=numberOfPointsInMesh){
			std::cout<<"number of points in mesh is not equal to polydata!"<<std::endl;
			return;
		}
		vtkPoints * vtkpoints =  m_PolyData->GetPoints();
		for(unsigned int p =0; p < numberOfPointsInPolydata; p++){
			double* apoint = vtkpoints->GetPoint( p );
			MeshType::PointType pt;
			for(unsigned int i=0;i<3; i++){
				pt[i] = static_cast<vtkFloatingPointType>(apoint[i]);
			}
			mesh->SetPoint( p, pt);
		}
	}
	
	// Convert ITK Mesh to VTK PolyData
	template<typename ITKMeshType>
	vtkSmartPointer<vtkPolyData>
		mesh2PolyData( const typename ITKMeshType* mesh )
	{
		typedef typename ITKMeshType::CellType CellType;
		typedef itk::Point<vtkFloatingPointType, 3> ItkPoint;
		typedef typename CellType::PointIdIterator PointIdIterator;
		typedef typename ITKMeshType::CellsContainer::ConstIterator CellIterator;
		typedef typename ITKMeshType::PointsContainer::ConstIterator PointIterator;
		//typedef typename ITKMeshType::LinesContainer
		vtkSmartPointer<vtkPolyData> newPolyData = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		PointIterator pntIterator = mesh->GetPoints()->Begin();
		PointIterator pntItEnd = mesh->GetPoints()->End();
		for (int i = 0; pntIterator != pntItEnd; ++i, ++pntIterator){
			ItkPoint pnt = pntIterator.Value();
			points->InsertPoint(i, pnt[0], pnt[1], pnt[2]);
		}
		newPolyData->SetPoints(points);
		vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
		CellIterator cellIt = mesh->GetCells()->Begin();
		CellIterator cellItEnd = mesh->GetCells()->End();
		unsigned long numberoflines = 0;
		unsigned long numberoftriangles = 0;
		unsigned long numberofpolygons = 0;
		unsigned long numberofunknown = 0;
		for (int it = 0; cellIt != cellItEnd; ++it, ++cellIt)
		{
			CellType * cellptr = cellIt.Value();
			PointIdIterator pntIdIter = cellptr->PointIdsBegin();
			PointIdIterator pntIdEnd = cellptr->PointIdsEnd();
			vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
			for (; pntIdIter != pntIdEnd; ++pntIdIter){
				pts->InsertNextId( *pntIdIter );
			}
			switch( cellptr->GetType() )
			{
			case CellType::LINE_CELL:
				{
					lines->InsertNextCell(pts);
					numberoflines++;
					break;
				}	
			case CellType::TRIANGLE_CELL:
				{
					cells->InsertNextCell(pts);
					numberoftriangles++;
					break;
				}
			case CellType::POLYGON_CELL:
				{
					cells->InsertNextCell(pts);
					numberofpolygons++;
					break;
				}
			default:
				{
					cells->InsertNextCell(pts);
					numberofunknown++;
					break;
				}
			}
		}
		if (numberofpolygons>0 || numberoftriangles>0){
			newPolyData->SetPolys(cells);
		}
		if(numberoflines>0){
			newPolyData->SetLines(lines);
		}
		return newPolyData;
	}

	//Conver VTK PolyData to ITK Mesh
	template<typename MeshType>
	typename MeshType::Pointer
		polyData2Mesh( vtkPolyData* m_PolyData, typename MeshType::Pointer refmesh = NULL)
	{
		typedef typename MeshType::PixelType           PixelType;
		typedef typename MeshType::CellType            CellType;
		typedef itk::VertexCell< CellType >            VertexCellType;
		typedef itk::LineCell< CellType >              LineCellType;
		typedef itk::TriangleCell< CellType >          TriangleCellType;
		typedef itk::PolygonCell< CellType >           PolygonCellType;
		typedef typename MeshType::CellAutoPointer     CellAutoPointer;
		typedef typename MeshType::PointIdentifier     PointIdentifier;
		// Create a new mesh
		MeshType::Pointer m_itkMesh = MeshType::New();
		//Points and point data
		const unsigned int numberOfPointsInPolydata = m_PolyData->GetNumberOfPoints();
		unsigned int numberOfPointsInMesh = numberOfPointsInPolydata;
		vtkPoints * vtkpoints =  m_PolyData->GetPoints();
		KM_DEBUG_PRINT( "numberOfPoints: ", numberOfPointsInMesh );
		m_itkMesh->GetPoints()->Reserve( numberOfPointsInMesh );
		m_itkMesh->GetPointData()->Reserve( numberOfPointsInMesh );
		m_itkMesh->GetGeometryData()->Reserve( numberOfPointsInMesh );
		for(unsigned int p =0; p < numberOfPointsInPolydata; p++){
			double* apoint = vtkpoints->GetPoint( p );
			MeshType::PointType pt;
			for(unsigned int i=0;i<3; i++){
				pt[i] = static_cast<vtkFloatingPointType>(apoint[i]);
			}
			m_itkMesh->SetPoint( p, pt);
			m_itkMesh->SetPointData( p, static_cast<PixelType>( 0 ) );
			if (NULL!=refmesh){
				SimplexMeshGeometry *refdata = refmesh->GetGeometryData()->GetElement(p);
				SimplexMeshGeometry *data = new itk::SimplexMeshGeometry();
				data->CopyFrom( *refdata );
				m_itkMesh->SetGeometryData( p, data );
			}
		}
		//Polygons and cell data
		const unsigned int numberOfPolysInPolydata = m_PolyData->GetNumberOfPolys();
		const unsigned int numberOfLinesInPolydata = m_PolyData->GetNumberOfLines();

		//Count the number of all cells which will be inserted into mesh.
		unsigned int numberOfCellsInMesh = numberOfPolysInPolydata;
		vtkCellArray * lines = m_PolyData->GetLines();
		lines->InitTraversal();
		vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
		while( lines->GetNextCell( cellPoints ) ){
			numberOfCellsInMesh += ( cellPoints->GetNumberOfIds() - 1);
		}
		KM_DEBUG_PRINT( "numberOfCells: ", numberOfCellsInMesh );
		m_itkMesh->GetCells()->Reserve( numberOfCellsInMesh );
		m_itkMesh->GetCellData()->Reserve( numberOfCellsInMesh );

		vtkCellArray * polygons = m_PolyData->GetPolys();
		polygons->InitTraversal();
		vtkIdType cellId = 0;
		while( polygons->GetNextCell( cellPoints ) ){
			MeshType::CellAutoPointer c;
			PolygonCellType *newCell = new PolygonCellType;
			for ( unsigned int jj = 0; jj < cellPoints->GetNumberOfIds(); jj++ ){
				newCell->SetPointId( jj, cellPoints->GetId(jj) );
			}
			c.TakeOwnership(newCell);
			m_itkMesh->AddFace( c );
			m_itkMesh->SetCellData( cellId, static_cast<PixelType>( 0 ) );
			cellId++;
		}
		lines->InitTraversal();
		while( lines->GetNextCell( cellPoints ) ){
			MeshType::CellAutoPointer c;
			unsigned int numberOfIdsInLine = cellPoints->GetNumberOfIds();
			for (int i=0;i<numberOfIdsInLine-1;i++){
				vtkIdType pt1 = cellPoints->GetId( i );
				vtkIdType pt2 = cellPoints->GetId( i+1 );
				m_itkMesh->AddEdge( pt1, pt2 );
				m_itkMesh->SetCellData( cellId, static_cast<PixelType>( 0 ) );
				cellId++;
			}
		}
		m_itkMesh->BuildCellLinks();
		return m_itkMesh;
	}

	template<class MeshType>
	typename MeshType::Pointer
		decimateMesh( typename MeshType* inputMesh , unsigned int numberOfPoints = 3000)
	{
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		polydata = km::mesh2PolyData<MeshType>(inputMesh);
		polydata = km::smoothPolyData( polydata, 30 );
		polydata = km::decimatePolydata(polydata, static_cast<double>(numberOfPoints)/inputMesh->GetNumberOfPoints());
		MeshType::Pointer outputMesh = polyData2Mesh<MeshType>( polydata );
		return outputMesh;
	}

	//If computeAll=false, then only normal will be computed. Otherwise all geometry data like mean curve, sphere radius etc. will be computed.
	template<class MeshType>
	void
		ComputeGeometry(typename MeshType* inputMesh, bool computeAll = false)
	{
		typedef typename MeshType::GeometryMapType      GeometryMapType;
		typedef typename GeometryMapType::Pointer       GeometryMapPointer;
		typedef typename GeometryMapType::Iterator      GeometryMapIterator;
		typename MeshType::PointsContainerPointer allpoints = inputMesh->GetPoints();
		typename MeshType::GeometryMapPointer     allgeodata   = inputMesh->GetGeometryData();
		typename MeshType::PointDataContainerPointer allpointdata = inputMesh->GetPointData();
		if (allpointdata->Size() <= 0){
			allpointdata->Reserve( allpoints->Size() );
		}
		typedef typename SimplexMeshGeometry::PointType PointType;
		typedef typename PointType::VectorType          VectorType;
		typedef CovariantVector< typename VectorType::ValueType, 3 > CovariantVectorType;
		PointType           Foot;
		CovariantVectorType normal;
		CovariantVectorType z;
		VectorType          tmp;
		GeometryMapType::Iterator dataIt = allgeodata->Begin();
		SimplexMeshGeometry *data;
		while ( dataIt != allgeodata->End() )
		{
			data = dataIt.Value();
			data->neighbors[0] = allpoints->GetElement(data->neighborIndices[0]);
			data->neighbors[1] = allpoints->GetElement(data->neighborIndices[1]);
			data->neighbors[2] = allpoints->GetElement(data->neighborIndices[2]);

			// compute normal
			normal.Fill(0.0);
			z.Set_vnl_vector( itk_cross_3d( ( data->neighbors[1] - data->neighbors[0] ).Get_vnl_vector(),
				( data->neighbors[2] - data->neighbors[0] ).Get_vnl_vector() ) );
			z.Normalize();
			normal += z;
			// copy normal
			data->normal = normal;
			data->pos = allpoints->GetElement( dataIt.Index() );

			// compute the simplex angle
			data->ComputeGeometry();

			if ( computeAll ){
				tmp = data->neighbors[0] - data->pos;
				double D = 1.0 / ( 2 * data->sphereRadius ); /* */
				double tmpNormalProd = dot_product( tmp.GetVnlVector(), data->normal.GetVnlVector() );
				double sinphi =  2 *data->circleRadius *D *vnl_math_sgn(tmpNormalProd);
				double phi = vcl_asin(sinphi);
				data->phi = phi;
				data->meanCurvature = vcl_abs(sinphi / data->circleRadius);
				tmp = data->pos - data->neighbors[0];
				//compute the foot of p projection of p onto the triangle spanned by its
				// neighbors
				double distance = -tmpNormalProd;
				tmp.SetVnlVector( ( data->pos ).GetVnlVector() - distance * normal.GetVnlVector() );
				Foot.Fill(0.0);
				Foot += tmp;
				data->distance = ( ( data->circleCenter ) - Foot ).GetNorm();
				{
					PointType a, b, c;
					a = data->neighbors[0];
					b = data->neighbors[1];
					c = data->neighbors[2];
					VectorType n, na, nb, nc;
					n.SetVnlVector( itk_cross_3d( ( b - a ).GetVnlVector(), ( c - a ).GetVnlVector() ) );
					na.SetVnlVector( itk_cross_3d( ( c - b ).GetVnlVector(), ( Foot - b ).GetVnlVector() ) );
					nb.SetVnlVector( itk_cross_3d( ( a - c ).GetVnlVector(), ( Foot - c ).GetVnlVector() ) );
					nc.SetVnlVector( itk_cross_3d( ( b - a ).GetVnlVector(), ( Foot - a ).GetVnlVector() ) );
					PointType eps;
					eps[0] = dot_product( n.GetVnlVector(), na.GetVnlVector() ) / n.GetSquaredNorm();
					eps[1] = dot_product( n.GetVnlVector(), nb.GetVnlVector() ) / n.GetSquaredNorm();
					eps[2] = dot_product( n.GetVnlVector(), nc.GetVnlVector() ) / n.GetSquaredNorm();
					data->eps = eps;
				}
			}
			allpointdata->InsertElement( dataIt->Index(), data->phi );
			dataIt.Value() = data;
			dataIt++;
		}
	}

	template<class MeshType>
	void
		adaptMesh(typename MeshType* inputMesh, double m_Gamma, unsigned int m_Iterations)
	{
		typedef itk::AdaptSimplexMesh3DFilter<MeshType, MeshType> AdaptSimplexMesh3DFilterType;
		AdaptSimplexMesh3DFilterType::Pointer adaptorFilter = AdaptSimplexMesh3DFilterType::New();
		adaptorFilter->SetInput( inputMesh );
		adaptorFilter->SetGamma( m_Gamma );
		adaptorFilter->SetTangentFactor( 0.8 );
		adaptorFilter->SetNormalFactor( 0.2 );
		adaptorFilter->SetRigidity( 1 );
		adaptorFilter->SetIterations( m_Iterations );
		adaptorFilter->Update(); //This is just for initialization.
	}

	template<class MeshType>
	void
		smoothMesh(typename MeshType* inputMesh, double factor, unsigned int m_Iterations)
	{
		typedef itk::AdaptSimplexMesh3DFilter<MeshType, MeshType> AdaptSimplexMesh3DFilterType;
		AdaptSimplexMesh3DFilterType::Pointer adaptorFilter = AdaptSimplexMesh3DFilterType::New();
		adaptorFilter->SetInput( inputMesh );
		adaptorFilter->SetGamma( 0.0 );
		adaptorFilter->SetTangentFactor( 0.0 );
		adaptorFilter->SetNormalFactor( factor );
		adaptorFilter->SetRigidity( 1 );
		adaptorFilter->SetIterations( m_Iterations );
		adaptorFilter->Update(); //This is just for initialization.
	}

	template<class MeshType>
	void
		assigneMesh(typename MeshType::Pointer mesh, typename MeshType::PixelType value)
	{
		typedef MeshType::PointsContainer::Iterator      PointsIterator;
		mesh->GetPointData()->Reserve( mesh->GetNumberOfPoints() );
		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();
		while ( pointItr != pointEnd ){
			mesh->SetPointData( pointItr.Index(), value );	
			++pointItr;
		}
	}

	template<class MeshType>
	void
		smoothMeshData(typename MeshType::Pointer mesh, unsigned int radius)
	{
		MeshType::PointDataContainerPointer allpointdata = mesh->GetPointData();
		if (allpointdata.IsNull() || allpointdata->Size() <= 0){
			return;
		}
		typedef MeshType::PointDataContainer::ElementIdentifier ElementIdentifier;
		typedef MeshType::PointsContainer::Iterator PointsIterator;
		PointsIterator pointItr = mesh->GetPoints()->Begin();
		PointsIterator pointEnd = mesh->GetPoints()->End();
		while ( pointItr != pointEnd ){
			MeshType::PointIdentifier idx = pointItr.Index();
			double datasum = allpointdata->GetElement( idx );
			MeshType::NeighborListType * neighborlist = mesh->GetNeighbors(idx, radius);
			for (int i=0;i<neighborlist->size();i++){
				datasum += allpointdata->GetElement( (*neighborlist)[i] );
			}
			datasum /= (1+neighborlist->size());
			allpointdata->InsertElement( idx, datasum );
			++pointItr;
		}
	}
	
	//-------------------------------------------------------//
	//                      VTK related                      //
	//-------------------------------------------------------//
	//Calculate the normals of a PolyData
	vtkSmartPointer<vtkPolyData> calculateNormals( vtkPolyData* inputPolyData );
	
	//Smooth a PolyData 
	vtkSmartPointer<vtkPolyData> smoothPolyData( vtkPolyData * inputPolyData, int iterations, double relaxfactor = 0.025 );
	
	//Decimate a vtkpolydata
	vtkSmartPointer<vtkPolyData> decimatePolydata( vtkSmartPointer<vtkPolyData> inputPolyData, double newRatio );
	
	//Read VTKPolyData
	vtkSmartPointer<vtkPolyData> readPolyData( const std::string filename );
	
	//Write VTKPolyData
	void writePolyData( const std::string filename, vtkPolyData* polydata );
	void writePolyData( const std::string dirname, const std::string filenameprefix, unsigned int index, const char* extension, vtkPolyData* polydata);
}

#endif
