/*
Classification construction based on Hongzhi Wang's AdaBoost code.
*/
#include "util.h"
#include "AdaBoost.h"


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCommand.h"
#include "itkImageMaskSpatialObject.h"
#include "itkResampleImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkComposeImageFilter.h"

#include "itkDiscreteGaussianDerivativeImageFunction.h"
#include "itkDiscreteGradientMagnitudeGaussianImageFunction.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h" 
#include "itkLabelMapToLabelImageFilter.h"
#include "itkImageSample.h"
#include "itkVectorDataContainer.h"
#include "itkImageRandomSampler.h"
#include "itkImageSamplerBase.h"
#include "itkGaussianBlurImageFunction.h"
#include "itkImageMaskSpatialObject2.h"
#include "itkBinaryContourImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h" 
#include "itkLabelMapToLabelImageFilter.h"
#include "itkImageSample.h"
#include "itkVectorDataContainer.h"
#include "itkImageRandomSampler.h"
#include "itkImageSamplerBase.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkImageFullSampler.h"
#include "itkImageRandomSamplerSparseMask.h"

typedef  unsigned char                         MaskPixelType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
typedef itk::Image< PixelType, Dimension >     ImageType;
typedef itk::VectorImage<float, 3>             VectorImageType;
typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;
typedef itk::ImageFileReader< ImageType >      ReaderType;
typedef itk::ImageFileWriter< ImageType >      WriterType;
typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;
typedef itk::MinimumMaximumImageFilter<MaskImageType> MinimumMaximumImageFilterType;

typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>  GaussianFilterType;			//Í¼ÏñÆ½»¬£ºblurs an image by separable convolution with discrete guassian kernels
typedef itk::DerivativeImageFilter<ImageType,ImageType>         DerivativeFilterType;		//computer the directional derivative of an image.
typedef itk::ShiftScaleImageFilter<ImageType,ImageType>         ShiftScaleFilterType;	//shift and scale the pixel in an image
typedef itk::ImageRandomConstIteratorWithIndex<ImageType>       RandomIndexIteratorType;		//an multi-dimentsional image iterator that visits a random set of pixels within an image region
typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType>  GradientMagnitudeFilterType;		//computes the gradient magnitude of an image region at each pixel
typedef itk::ComposeImageFilter<ImageType>                      ImageToVectorImageFilterType;
typedef itk::ImageRandomConstIteratorWithIndex<VectorImageType> RandomVectorIteratorType;
typedef itk::BinaryImageToShapeLabelMapFilter< MaskImageType >  BinaryImageToShapeLabelMapFilterType;		//Converts a binary image to a label map and valuate the shape attributes

typedef itk::ImageRegionConstIterator< MaskImageType > ConstIteratorType2;

typedef itk::ImageRegionConstIteratorWithIndex<ImageType>                    RegionIteratorConstIteratorType;
typedef itk::DiscreteGaussianDerivativeImageFunction< ImageType, PixelType > GaussianDerivativeImageFunctionType;
typedef itk::CentralDifferenceImageFunction<ImageType>                       CentralDifferenceFunctionType;
typedef itk::DiscreteHessianGaussianImageFunction<ImageType>                 DiscreteHessianGaussianImageFunctionType;
typedef itk::DiscreteGradientMagnitudeGaussianImageFunction<ImageType>       DiscreteGradientMagnitudeFunctionType;
typedef itk::GaussianBlurImageFunction<ImageType,float>						 GaussianBlurFunctionType;
typedef itk::SignedMaurerDistanceMapImageFilter< MaskImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;
typedef itk::BinaryImageToLabelMapFilter< MaskImageType >  BinaryImageToLabelMapFilterType;		//Converts a binary image to a label map and valuate the shape attributes

typedef itk::BinaryThresholdImageFilter<MaskImageType,MaskImageType> BinaryThresholdFilterType;
typedef itk::ImageMaskSpatialObject2<Dimension>           ImageMaskSpatialObject2;
typedef itk::ImageSamplerBase< MaskImageType >            SamplerBaseType;
typedef itk::ImageRandomSampler< MaskImageType >          RandomSamplerType;
typedef SamplerBaseType::ImageSampleContainerType         SampleContainerType;
typedef SamplerBaseType::ImageSampleType                  SampleType;
typedef MaskImageType::IndexType IndexType;
typedef IndexType::IndexValueType  IndexValueType;
typedef itk::ContinuousIndex<double, Dimension> CIndexType;


typedef itk::ImageRandomSamplerSparseMask< MaskImageType >    RandomSparseMaskSamplerType;



typedef itk::ImageFullSampler< MaskImageType >                FullSamplerType;  


using namespace std;

int main( int argc, char ** argv )
{
	if ( argc != 4 )
	{
		std::cerr << "Missing parameters. " << std::endl;
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0]
		<< " manualSegmentationFile sampleRate outputImage" <<endl<<endl
			<< " Meanings of the parameters: "<<endl
			<< " manualSegmentationFile: labels "<<std::endl<<endl
			<< " output:                 output the random resampled voxel coordinate to a binary file"<<endl<<endl	
			<< " sampleNum:             number of sample (including background and foreground voxel"<<endl
			<< std::endl;
		return -1;
	}
	

	MaskReaderType::Pointer reader = MaskReaderType::New();

	reader->SetFileName(argv[1]);

	try
	{
		reader->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err.GetDescription() << std::endl;
		return -1;
	}



	int samplenum=atoi(argv[3]);
	
	MinimumMaximumImageFilterType::Pointer minmax = MinimumMaximumImageFilterType::New();

	minmax->SetInput(reader->GetOutput());
	minmax->Update();
	

	std::cout << "Number of labels: " << static_cast<unsigned int>( minmax->GetMaximum()) << std::endl;



	unsigned int labelnum = static_cast<int>(minmax->GetMaximum())+1;

	unsigned int *labels = new unsigned int [labelnum ];

	for(int i=0;i<labelnum;i++)
		labels[i]=0;

	ConstIteratorType2 it(reader->GetOutput(),reader->GetOutput()->GetLargestPossibleRegion());
	
	it.GoToBegin();
	while(!it.IsAtEnd())
	{
		labels[it.Get()]++;
		++it;
	}




	MaskImageType::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();

	MaskImageType::SizeType sizes = region.GetSize();

	unsigned int numbers = sizes[0]*sizes[1]*sizes[2];



	// location features
	MaskImageType::Pointer outputImage = MaskImageType::New();
	outputImage->SetRegions( reader->GetOutput()->GetRequestedRegion() );
	outputImage->CopyInformation( reader->GetOutput() );
	outputImage->Allocate();
	outputImage->FillBuffer(0);
	

	for(int i=0;i<labelnum;i++)
	{


		typedef itk::BinaryThresholdImageFilter <MaskImageType,MaskImageType> ThresholdImageFilterType;
		ThresholdImageFilterType::Pointer thresholdFilter = ThresholdImageFilterType::New();
		thresholdFilter->SetInput(reader->GetOutput());
		thresholdFilter->SetLowerThreshold(i);
		thresholdFilter->SetUpperThreshold(i);
		thresholdFilter->SetInsideValue(1);
		thresholdFilter->SetOutsideValue(0);
		thresholdFilter->Update();


		MaskWriterType::Pointer writer = MaskWriterType::New();
		writer->SetFileName("./tmp.mhd");
		writer->SetInput(thresholdFilter->GetOutput());
		writer->Update();
		

		ImageMaskSpatialObject2::Pointer maskSO = ImageMaskSpatialObject2::New();
		maskSO->SetImage(thresholdFilter->GetOutput());
		maskSO->Update();

		SamplerBaseType::Pointer sampler = 0;

		if(samplenum/labelnum < labels[i]/10 )
		{	
		//	std::cout<<"sample from label "<<i<<std::endl;
			RandomSamplerType::Pointer tempsampler1 = RandomSamplerType::New();
			tempsampler1->SetNumberOfSamples(samplenum/labelnum);
			sampler = tempsampler1;
		}
		else if(samplenum/labelnum > labels[i]/10 && samplenum/labelnum< labels[i])
		{
			RandomSparseMaskSamplerType::Pointer tempSampler = RandomSparseMaskSamplerType::New();    
			tempSampler->SetNumberOfSamples( samplenum/labelnum );
			sampler = tempSampler;
		}
		else
		{
			FullSamplerType::Pointer tempsampler1 = FullSamplerType::New();
		//	tempsampler1->SetNumberOfSamples(labels[i]);
			sampler = tempsampler1;
		//	std::cout<<"sample from label "<<i<<" Full sample!"<<std::endl;
		}


		sampler->SetInput( reader->GetOutput() );
		sampler->SetInputImageRegion( reader->GetOutput()->GetLargestPossibleRegion() );
		sampler->SetMask(maskSO);


		try
		{
			sampler->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "label "<<i<<std::endl;
			std::cerr << err << std::endl;
			std::cerr << "number of samples still obtained: " 
				<< sampler->GetOutput()->Size() << std::endl;
			return -1;
		}


		SampleContainerType::Pointer samples = sampler->GetOutput();

		IndexType index;
		CIndexType cindex;
		SampleType sample;


		for (unsigned int i = 0; i < samples->Size(); ++i)
		{

			sample = samples->ElementAt( i );
			reader->GetOutput()->TransformPhysicalPointToContinuousIndex( sample.m_ImageCoordinates, cindex);
			for ( unsigned int d = 0; d < Dimension; ++d )
			{
				index[d] = static_cast<IndexValueType>( vnl_math_rnd( cindex[d] ) );
			}

			outputImage->SetPixel(index,1);

		}

	}

		
	MaskWriterType::Pointer writer = MaskWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(outputImage);
	try
	{
		writer->Update();
	}
	catch ( itk::ExceptionObject &err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err.GetDescription() << std::endl;
		return -1;
	}


	return 0;
}
