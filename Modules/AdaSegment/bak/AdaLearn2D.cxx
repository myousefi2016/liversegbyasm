/*
Classification construction based on Hongzhi Wang's AdaBoost code.
*/

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
#include "itkImageRegionIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkComposeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkNumericTraits.h"
#include "itkContinuousIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageMaskSpatialObject2.h"
#include "itkBinaryContourImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkBinaryThresholdImageFilter.h" 
#include "itkLabelMapToLabelImageFilter.h"
#include "itkImageSample.h"
#include "itkVectorDataContainer.h"
#include "itkImageRandomSampler.h"
#include "itkImageSamplerBase.h"

const int Dimension = 2;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::ImageFileReader< ImageType >   ReaderType;
typedef itk::ImageFileWriter< ImageType >   WriterType;
typedef float WritePixelType;
typedef itk::Image< WritePixelType, Dimension > WriteImageType;
typedef itk::ImageLinearConstIteratorWithIndex< ImageType >  ConstIteratorType;
typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;
typedef itk::ImageRegionIterator< ImageType>        IteratorType;
typedef  unsigned char                         MaskPixelType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
typedef itk::Image< PixelType, Dimension >     ImageType;
typedef itk::VectorImage<float, Dimension>     VectorImageType;
typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;
typedef itk::ImageFileReader< ImageType >      ReaderType;
typedef itk::ImageFileWriter< ImageType >      WriterType;
typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;

typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>  GaussianFilterType;			//图像平滑：blurs an image by separable convolution with discrete guassian kernels
typedef itk::DerivativeImageFilter<ImageType,ImageType>         DerivativeFilterType;		//computer the directional derivative of an image.
typedef itk::ShiftScaleImageFilter<ImageType,ImageType>         ShiftScaleFilterType;		//shift and scale the pixel in an image
typedef itk::ImageRandomConstIteratorWithIndex<ImageType>       RandomIndexIteratorType;		//an multi-dimentsional image iterator that visits a random set of pixels within an image region
typedef itk::GradientMagnitudeImageFilter<ImageType,ImageType>  GradientMagnitudeFilterType;		//computes the gradient magnitude of an image region at each pixel
typedef itk::ComposeImageFilter<ImageType>                      ImageToVectorImageFilterType;
typedef itk::ImageRandomConstIteratorWithIndex<VectorImageType> RandomVectorIteratorType;
typedef itk::BinaryImageToLabelMapFilter< MaskImageType >  BinaryImageToLabelMapFilterType;		//Converts a binary image to a label map and valuate the shape attributes

typedef itk::BinaryThresholdImageFilter<MaskImageType,MaskImageType> BinaryThresholdFilterType;
typedef itk::ImageMaskSpatialObject2<Dimension>  ImageMaskSpatialObject2;
typedef itk::ImageSamplerBase< MaskImageType >            SamplerBaseType;
typedef itk::ImageRandomSampler< MaskImageType >          RandomSamplerType;
typedef SamplerBaseType::ImageSampleContainerType         SampleContainerType;
typedef SamplerBaseType::ImageSampleType                  SampleType;
typedef MaskImageType::IndexType IndexType;
typedef IndexType::IndexValueType  IndexValueType;
typedef itk::ContinuousIndex<double, Dimension> CIndexType;

typedef  itk::SignedMaurerDistanceMapImageFilter< MaskImageType, ImageType  > SignedMaurerDistanceMapImageFilterType;

using namespace std;

int main( int argc, char ** argv )
{
	if ( argc != 6 )
	{
		std::cerr << "Missing parameters. " << std::endl;
		std::cerr << "Usage: " << std::endl;
		std::cerr << argv[0]
		<< " inputImageFile manualSegmentationFile sampleRatio iteration AdaBoostOutputPrefix" <<endl<<endl
			<< " Learns to construct a Adaboost classifier based on multi-features."<<endl<<endl
			<< " Meanings of the parameters: "<<endl
			<< " inputImageFile:         A text file saves the locations of training images." <<endl
			<< "                         Each row of the file stores the location of one training image."<<endl
			<< "                         For example, a training set with two training images may have an"<<endl
			<< "                         inputImageFile with the following content:" <<endl 
			<< "                         /home/subject01.nii.gz" <<endl
			<< "                         /home/subject02.nii.gz" <<endl<<endl
			<< " manualSegmentationFile: A text file saves the locations of the corresponding manual segmentationds"<<endl
			<< "                         for the training images."<<endl
			<< "                         Each row of the file stores the location of one training image"<<endl
			<< "                         For example, the manualSegmentationFile for  the above example may have"<<endl
			<< "                         the following content:" <<endl 
			<< "                         /home/subject01_manualSeg.nii.gz" <<endl
			<< "                         /home/subject02_manulSeg.nii.gz" <<endl<<endl
			<< " sampleRate:             0<= sampleRate <=1. When the ROI is large, loading every voxel in ROI as a"<<endl
			<< "                         training sample may be impossible due to the memory limit. To address this "<<endl
			<< "                         problem, sampleRate specifies the percentage voxels from ROI that will be used"<<endl
			<< "                         as training samples. If sampleRate=0.01, 1 percent voxels will be used."<<endl<<endl
			<< " iteration:              We use AdaBoost learning. This parameter specifies the number of iterations"<<endl
			<< "                         for this learning task. Typically, hundreds of iterations are needed."<<endl<<endl
			<< " AdaBoostOutputPrefix:   Our program stores the learned results in two text files. This parameter"<<endl
			<< "                         specifies the path and the prefix of these two file names. For instance, if"<<endl
			<< "                         this parameter is given as /home/mytest, then the output files are"<<endl
			<< "                         /home/mytest-AdaBoostResults"<<endl
			<< "                         where ? represents the target label. For one host segmentation method, the"<<endl
			<< "                         learning task for each segmentation label should have the same AdaBoostOutputPrefix."<<endl
			<< std::endl;
		return -1;
	}


	ifstream imfile;
	string imstr;
	string manualstr;

	imfile.open(argv[1]);
	int subjectN=0;
	string t;

	while (!imfile.eof()){
		getline(imfile,t);
		if (t.length()>0)
			subjectN++;
	}

	imfile.clear();
	imfile.seekg(0, ios::beg);

	cout<<"subject #: "<<subjectN<<endl;

	ifstream manualfile;
	manualfile.open(argv[2]);

	vector<FLOATTYPE> X;
	vector<char> Y;
	int totalSample=0;


	FLOATTYPE sampleRatio=atof(argv[3]);

	cout<<"sampleRatio: "<<sampleRatio<<endl;
	if (sampleRatio>1)
	{
		std::cerr << "sampleRatio can not be greater than 1!"<<endl;
		return -1;  
	}
	int sampleRate=int(1.0/sampleRatio);
	int iterN=atoi(argv[4]);

	int fgcounter=0;
	int NFeature = 23;


	while (!imfile.eof())
	{
		getline(imfile,imstr);
		cout<<"image file: "<<imstr<<endl;
		getline(manualfile,manualstr);
		cout<<"manual seg: "<<manualstr<<endl;


		WriterType::Pointer writer = WriterType::New();
		ReaderType::Pointer im1 = ReaderType::New();
		im1->SetFileName(imstr.c_str());
		try
		{
			im1->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		MaskReaderType::Pointer manualSeg = MaskReaderType::New();
		manualSeg->SetFileName(manualstr.c_str());
		try
		{
			manualSeg->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}


		const int max_order = 2;
		const int no_scales = 3;
		float sigmas2[3]={0.25,1.0,4.0};

		GaussianFilterType::Pointer filter = GaussianFilterType::New();
		DerivativeFilterType::Pointer filterX = DerivativeFilterType::New();
		DerivativeFilterType::Pointer filterY = DerivativeFilterType::New();

		ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();
		filterX->SetDirection(0);
		filterY->SetDirection(1);


		//有的分割标记不一定是1，但肯定大于0，所以做一个阈值滤波
		BinaryThresholdFilterType::Pointer thresholder = BinaryThresholdFilterType::New();


		thresholder->SetInput(manualSeg->GetOutput());
		thresholder->SetLowerThreshold(1);
		thresholder->SetUpperThreshold(255);
		thresholder->SetInsideValue(1);
		thresholder->SetOutsideValue(0);
		thresholder->Update();

		BinaryImageToLabelMapFilterType::Pointer binaryer = BinaryImageToLabelMapFilterType::New();
		binaryer->SetInput(thresholder->GetOutput());
		binaryer->SetInputForegroundValue(1);
		binaryer->Update();


		typedef itk::LabelMapToLabelImageFilter<BinaryImageToLabelMapFilterType::OutputImageType, MaskImageType> LabelMapToLabelImageFilterType;
		LabelMapToLabelImageFilterType::Pointer labelMapToLabelImageFilter = LabelMapToLabelImageFilterType::New();
		labelMapToLabelImageFilter->SetInput(binaryer->GetOutput());
		labelMapToLabelImageFilter->Update();

		typedef itk::LabelStatisticsImageFilter< MaskImageType, MaskImageType > LabelStatisticsImageFilterType;
		LabelStatisticsImageFilterType::Pointer labelStatisticsImageFilter = LabelStatisticsImageFilterType::New();
		labelStatisticsImageFilter->SetLabelInput( labelMapToLabelImageFilter->GetOutput() );
		labelStatisticsImageFilter->SetInput(labelMapToLabelImageFilter->GetOutput());
		labelStatisticsImageFilter->Update();

		std::cout << "Number of labels: " << labelStatisticsImageFilter->GetNumberOfLabels() << std::endl;

		typedef LabelStatisticsImageFilterType::ValidLabelValuesContainerType ValidLabelValuesType;
		typedef LabelStatisticsImageFilterType::LabelPixelType                LabelPixelType;

		for(ValidLabelValuesType::const_iterator vIt=labelStatisticsImageFilter->GetValidLabelValues().begin();
			vIt != labelStatisticsImageFilter->GetValidLabelValues().end();
			++vIt)
		{
			if ( labelStatisticsImageFilter->HasLabel(*vIt) )
			{
				LabelPixelType labelValue = *vIt;
				std::cout << "min: " << labelStatisticsImageFilter->GetMinimum( labelValue ) << std::endl;
				std::cout << "max: " << labelStatisticsImageFilter->GetMaximum( labelValue ) << std::endl;
				std::cout << "median: " << labelStatisticsImageFilter->GetMedian( labelValue ) << std::endl;
				std::cout << "mean: " << labelStatisticsImageFilter->GetMean( labelValue ) << std::endl;
				std::cout << "sigma: " << labelStatisticsImageFilter->GetSigma( labelValue ) << std::endl;
				std::cout << "variance: " << labelStatisticsImageFilter->GetVariance( labelValue ) << std::endl;
				std::cout << "sum: " << labelStatisticsImageFilter->GetSum( labelValue ) << std::endl;
				std::cout << "count: " << labelStatisticsImageFilter->GetCount( labelValue ) << std::endl;
				//std::cout << "box: " << labelStatisticsImageFilter->GetBoundingBox( labelValue ) << std::endl; // can't output a box
				std::cout << "region: " << labelStatisticsImageFilter->GetRegion( labelValue ) << std::endl;
				std::cout << std::endl << std::endl;

			}
		}


		SignedMaurerDistanceMapImageFilterType::Pointer distanceMapImageFilter = SignedMaurerDistanceMapImageFilterType::New();

		typedef itk::BinaryContourImageFilter <MaskImageType, MaskImageType >  binaryContourImageFilterType;

		binaryContourImageFilterType::Pointer binaryContourFilter = binaryContourImageFilterType::New ();


		binaryContourFilter->SetInput(thresholder->GetOutput());
		binaryContourFilter->SetFullyConnected(1);
		binaryContourFilter->Update();

		distanceMapImageFilter->SetInput(binaryContourFilter->GetOutput());
		distanceMapImageFilter->Update();



		typedef itk::BinaryThresholdImageFilter <ImageType,MaskImageType> ThresholdImageFilterType;

	
		ThresholdImageFilterType::Pointer thresholdFilter_inside = ThresholdImageFilterType::New();
		ThresholdImageFilterType::Pointer thresholdFilter_outside = ThresholdImageFilterType::New();


		WriterType::Pointer dtwriter = WriterType::New();
		dtwriter->SetInput(distanceMapImageFilter->GetOutput());
		dtwriter->SetFileName("./dt.nii.gz");
		dtwriter->Update();


		thresholdFilter_inside->SetInput(distanceMapImageFilter->GetOutput());
		thresholdFilter_outside->SetInput(distanceMapImageFilter->GetOutput());


		thresholdFilter_inside->SetLowerThreshold(-10000.0);
		thresholdFilter_inside->SetUpperThreshold(0);
		thresholdFilter_inside->SetOutsideValue(0);
		thresholdFilter_inside->SetInsideValue(1);
		thresholdFilter_inside->Update();

		thresholdFilter_outside->SetLowerThreshold(1);
		thresholdFilter_outside->SetUpperThreshold(itk::NumericTraits<float>::max());
		thresholdFilter_outside->SetOutsideValue(0);
		thresholdFilter_outside->SetInsideValue(1);
		thresholdFilter_outside->Update();


		typedef itk::ImageFileWriter<MaskImageType> MaskWriterType;
		MaskWriterType::Pointer writer1 = MaskWriterType::New();

		writer1->SetInput(thresholdFilter_inside->GetOutput());
		writer1->SetFileName("./threshold_inside.nii.gz");
		writer1->Update();

		MaskWriterType::Pointer writer2 = MaskWriterType::New();
		writer2->SetInput(thresholdFilter_outside->GetOutput());
		writer2->SetFileName("./threshold_outside.nii.gz");
		writer2->Update();

		ImageMaskSpatialObject2::Pointer maskSO_inside = ImageMaskSpatialObject2::New();
		maskSO_inside->SetImage(thresholdFilter_inside->GetOutput());
		maskSO_inside->Update();

		ImageMaskSpatialObject2::Pointer maskSO_outside = ImageMaskSpatialObject2::New();
		maskSO_outside->SetImage(thresholdFilter_outside->GetOutput());
		maskSO_outside->Update();

		//BinaryImageToLabelMapFilterType::OutputImageType::LabelObjectType * labelObject = binaryer->GetOutput()->GetNthLabelObject(0);

		filter->SetInput(im1->GetOutput());
		filter->SetUseImageSpacingOn();

		char filename[256];
		char temp_dir[]="./";

		ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
		VectorImageType::Pointer vectorImage = imageToVectorImageFilter->GetOutput();


		int feat_n=0;
		for(int i=0; i<no_scales; i++)
		{
			for(int j=0; j<=max_order; j++)
			{
				// order j
				for(int ox=0; ox<=j; ox++)
				{
					for(int oy=0; oy<=j; oy++)
					{
						if(ox+oy==j)
						{
							filter->SetVariance(sigmas2[i]);
							filter->SetMaximumKernelWidth( 12*sigmas2[i] );
							filter->SetMaximumError(0.02);
							filter->Update();
							if(ox>0)
							{
								filterX->SetInput(filter->GetOutput());
								filterX->SetOrder(ox);
								filterX->Update();
							}
							if(oy>0)
							{
								if(ox>0)
								{
									filterY->SetInput(filterX->GetOutput());
								}
								else
								{
									filterY->SetInput(filter->GetOutput());
								}
								filterY->SetOrder(oy);
								filterY->Update();
							}

							ShiftScaleFilterType::Pointer shiftScale = ShiftScaleFilterType::New();

							ImageType::Pointer image = ImageType::New(); 


							if(oy>0)
							{
								filterY->Update();
								//		shiftScale->SetInput(filterY->GetOutput());
								image = filterY->GetOutput();
								image->DisconnectPipeline();
							} 
							else if(ox>0)
							{
								filterX->Update();
								//	shiftScale->SetInput(filterX->GetOutput());
								image = filterX->GetOutput();
								image->DisconnectPipeline();
							}
							else
							{
								filter->Update();
								//	shiftScale->SetInput(filter->GetOutput());
								image = filter->GetOutput();
								image->DisconnectPipeline();
							}

							imageToVectorImageFilter->SetInput(feat_n, image);

							feat_n++;

						}
					}
				}
			}
		}


		GaussianFilterType::Pointer filter1 = GaussianFilterType::New();
		GaussianFilterType::Pointer filter2 = GaussianFilterType::New();
		GaussianFilterType::Pointer filter3 = GaussianFilterType::New();

		filter1->SetInput(im1->GetOutput());
		filter2->SetInput(im1->GetOutput());
		filter3->SetInput(im1->GetOutput());

		filter1->SetVariance(sigmas2[0]);
		filter1->SetMaximumKernelWidth( 12*sigmas2[0] );
		filter1->SetMaximumError(0.02);

		filter2->SetVariance(sigmas2[1]);
		filter2->SetMaximumKernelWidth( 12*sigmas2[1] );
		filter2->SetMaximumError(0.02);

		filter3->SetVariance(sigmas2[2]);
		filter3->SetMaximumKernelWidth( 12*sigmas2[2] );
		filter3->SetMaximumError(0.02);


		GradientMagnitudeFilterType::Pointer gradientfilter1 = GradientMagnitudeFilterType::New();
		GradientMagnitudeFilterType::Pointer gradientfilter2 = GradientMagnitudeFilterType::New();
		GradientMagnitudeFilterType::Pointer gradientfilter3 = GradientMagnitudeFilterType::New();

		gradientfilter1->SetInput(filter1->GetOutput());
		gradientfilter1->Update();

		gradientfilter2->SetInput(filter2->GetOutput());
		gradientfilter2->Update();

		gradientfilter3->SetInput(filter3->GetOutput());
		gradientfilter3->Update();


		ImageType::Pointer gradient1 = gradientfilter1->GetOutput();
		ImageType::Pointer gradient2 = gradientfilter2->GetOutput();
		ImageType::Pointer gradient3 = gradientfilter3->GetOutput();


		gradient1->DisconnectPipeline();
		gradient2->DisconnectPipeline();
		gradient3->DisconnectPipeline();


		imageToVectorImageFilter->SetInput(feat_n++,gradient1);
		imageToVectorImageFilter->SetInput(feat_n++,gradient2);
		imageToVectorImageFilter->SetInput(feat_n++,gradient3);


		// location features
		ImageType::Pointer outputImageX = ImageType::New();
		ImageType::Pointer outputImageY = ImageType::New();
		outputImageX->SetRegions( im1->GetOutput()->GetRequestedRegion() );
		outputImageY->SetRegions( im1->GetOutput()->GetRequestedRegion() );
		outputImageX->CopyInformation( im1->GetOutput() );
		outputImageY->CopyInformation( im1->GetOutput() );
		outputImageX->Allocate();
		outputImageY->Allocate();

		IteratorType outputItX( outputImageX, outputImageX->GetRequestedRegion() );
		IteratorType outputItY( outputImageY, outputImageY->GetRequestedRegion() );
		ImageType::IndexType requestedIndexX =outputImageX->GetRequestedRegion().GetIndex();
		ImageType::IndexType requestedIndexY =outputImageY->GetRequestedRegion().GetIndex();
		ImageType::SizeType requestedSize = outputImageX->GetRequestedRegion().GetSize();
		for ( outputItX.GoToBegin(),outputItY.GoToBegin(); !outputItX.IsAtEnd(); ++outputItX,++outputItY)
		{
			ImageType::IndexType idx = outputItX.GetIndex();
			outputItX.Set(idx[0]+1);
			outputItY.Set(idx[1]+1);
		}


	//	imageToVectorImageFilter->SetInput(feat_n++,outputImageX);
	//	imageToVectorImageFilter->SetInput(feat_n++,outputImageY);
		imageToVectorImageFilter->Update();

		NFeature = feat_n;


		VectorImageType::Pointer vectorimage = imageToVectorImageFilter->GetOutput();
		vectorimage->DisconnectPipeline();

		RandomSamplerType::Pointer tempsampler1 = RandomSamplerType::New();
		RandomSamplerType::Pointer tempsampler2 = RandomSamplerType::New();

		tempsampler1->SetNumberOfSamples(2000);
		tempsampler2->SetNumberOfSamples(2000);

		tempsampler1->SetInput( manualSeg->GetOutput() );
		tempsampler1->SetInputImageRegion( manualSeg->GetOutput()->GetBufferedRegion() );
		tempsampler1->SetMask(maskSO_inside);				//.......

		tempsampler2->SetInput( manualSeg->GetOutput() );
		tempsampler2->SetInputImageRegion( manualSeg->GetOutput()->GetBufferedRegion() );
		tempsampler2->SetMask(maskSO_outside);		    

		std::cout << "Updating sampler..." << std::endl;
		try
		{
			tempsampler1->Update();
			tempsampler2->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << err << std::endl;
			std::cerr << "number of samples still obtained: " 
				<< tempsampler1->GetOutput()->Size() << std::endl;
			return -1;
		}

		/** Get sample container */
		std::cout << "Getting output of sampler..." << std::endl;
		SampleContainerType::Pointer samples = tempsampler1->GetOutput();


		/** Create image showing selected samples */
		std::cout << "Creating output image..." << std::endl;
		MaskImageType::Pointer outputImage = MaskImageType::New();
		outputImage->SetRegions( manualSeg->GetOutput()->GetLargestPossibleRegion().GetSize() );
		outputImage->CopyInformation( manualSeg->GetOutput() );
		outputImage->Allocate();
		outputImage->FillBuffer( 0 );  

		IndexType index;
		CIndexType cindex;
		SampleType sample;

		for (unsigned int i = 0; i < samples->Size(); ++i)
		{

			sample = samples->ElementAt( i );
			outputImage->TransformPhysicalPointToContinuousIndex( sample.m_ImageCoordinates, cindex);
			for ( unsigned int d = 0; d < Dimension; ++d )
			{
				index[d] = static_cast<IndexValueType>( vnl_math_rnd( cindex[d] ) );
			}
			//   outputImage->SetPixel( index, sample.m_ImageValue );
			//outputImage->SetPixel( index, 1);

			for(int j=0;j<NFeature;j++)
			{
				VectorImageType::PixelType pixel= vectorimage->GetPixel(index);
				X.push_back(pixel[j]);
			}

			Y.push_back(1);
			fgcounter++;

		}

		SampleContainerType::Pointer samples2 = tempsampler2->GetOutput();

		for (unsigned int i = 0; i < samples->Size(); ++i)
		{
			IndexType index;
			CIndexType cindex;
			sample = samples2->ElementAt( i );
			outputImage->TransformPhysicalPointToContinuousIndex( sample.m_ImageCoordinates, cindex);
			for ( unsigned int d = 0; d < Dimension; ++d )
			{
				index[d] = static_cast<IndexValueType>( vnl_math_rnd( cindex[d] ) );
			}
			//   outputImage->SetPixel( index, sample.m_ImageValue );
	//		outputImage->SetPixel( index, 1);

			for(int j=0;j<NFeature;j++)
			{
				VectorImageType::PixelType pixel= vectorimage->GetPixel(index);
				X.push_back(pixel[j]);
			}

			Y.push_back(-1);

		}

		totalSample+=4000;

	}

	std::cout<<fgcounter<<std::endl;


	char fileName[1024];
	sprintf (fileName, "%s-AdaBoostResults",argv[5]);
	ofstream AdaBoostParamFile(fileName, ios::out);

	AdaBoostParamFile.close();
	sprintf (fileName, "%s-AdaBoostResults",argv[5]);
	AdaBoostTrain(&X[0],&Y[0],totalSample,NFeature,iterN,fileName);
	return 0;
}
