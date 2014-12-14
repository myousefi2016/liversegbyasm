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
#include "itkImageRandomConstIteratorWithIndex.h"F
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


typedef  unsigned char                         MaskPixelType;
typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
typedef itk::Image< PixelType, Dimension >     ImageType;
typedef itk::VectorImage<float, 3>             VectorImageType;
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
typedef itk::BinaryImageToShapeLabelMapFilter< MaskImageType >  BinaryImageToShapeLabelMapFilterType;		//Converts a binary image to a label map and valuate the shape attributes

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
	WriterType::Pointer writer = WriterType::New();




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


	FLOATTYPE sampleRatio=atof(argv[3]);
	cout<<"sampleRatio: "<<sampleRatio<<endl;
	if (sampleRatio>1){
		std::cerr << "sampleRatio can not be greater than 1!"<<endl;
		return -1;  
	}
	int sampleRate=int(1.0/sampleRatio);
	int iterN=atoi(argv[4]);

	unsigned int maxKernelWidth = 100;
	unsigned int sample_num = 2000;
	float maxError = 0.02;
	float sigma[3]={0.7071,1.4142,2.8284};

	// Iterator all pixels

	const int max_order = 2;
	const int no_scales = 3;

	char filename[256];
	char temp_dir[]="./";


	vector<FLOATTYPE> X;
	int NFeature = 36;

	VectorImageType::PixelType pixel;
	VectorImageType::IndexType index;

	int subj = 0;

	float **features = new float *[subjectN*sample_num*2];

	for(int i=0;i<subjectN*sample_num*2;i++)
		features[i] = new float[NFeature];

    char *Y= new char[subjectN*sample_num*2];
	for(int i=0;i<subjectN*sample_num*2;i++)
		Y[i]=-1;


	while (!imfile.eof())
	{
		getline(imfile,imstr);
		cout<<"image file: "<<imstr<<endl;
		getline(manualfile,manualstr);
		cout<<"manual seg: "<<manualstr<<endl;

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

		ImageMaskSpatialObject2::Pointer maskSO_inside = ImageMaskSpatialObject2::New();
		maskSO_inside->SetImage(thresholdFilter_inside->GetOutput());
		maskSO_inside->Update();

		ImageMaskSpatialObject2::Pointer maskSO_outside = ImageMaskSpatialObject2::New();
		maskSO_outside->SetImage(thresholdFilter_outside->GetOutput());
		maskSO_outside->Update();


		RandomSamplerType::Pointer tempsampler1 = RandomSamplerType::New();
		RandomSamplerType::Pointer tempsampler2 = RandomSamplerType::New();

		tempsampler1->SetNumberOfSamples(sample_num);
		tempsampler2->SetNumberOfSamples(sample_num);

		tempsampler1->SetInput( manualSeg->GetOutput() );
		tempsampler1->SetInputImageRegion( manualSeg->GetOutput()->GetBufferedRegion() );
		tempsampler1->SetMask(maskSO_inside);				//.......

		tempsampler2->SetInput( manualSeg->GetOutput() );
		tempsampler2->SetInputImageRegion( manualSeg->GetOutput()->GetBufferedRegion() );
		tempsampler2->SetMask(maskSO_outside);		    

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

		
		SampleContainerType::Pointer samples = tempsampler1->GetOutput();

		IndexType index;
		CIndexType cindex;
		SampleType sample;

	

		std::vector<ImageType::IndexType> sample_array;

		for (unsigned int i = 0; i < samples->Size(); ++i)
		{

			sample = samples->ElementAt( i );
			im1->GetOutput()->TransformPhysicalPointToContinuousIndex( sample.m_ImageCoordinates, cindex);
			for ( unsigned int d = 0; d < Dimension; ++d )
			{
				index[d] = static_cast<IndexValueType>( vnl_math_rnd( cindex[d] ) );
			}

			std::cout<<index<<std::endl;
			sample_array.push_back(index);

		}

		SampleContainerType::Pointer samples2 = tempsampler2->GetOutput();

		for (unsigned int i = 0; i < samples->Size(); ++i)
		{
			IndexType index;
			CIndexType cindex;
			sample = samples2->ElementAt( i );
			im1->GetOutput()->TransformPhysicalPointToContinuousIndex( sample.m_ImageCoordinates, cindex);
			for ( unsigned int d = 0; d < Dimension; ++d )
			{
				index[d] = static_cast<IndexValueType>( vnl_math_rnd( cindex[d] ) );
			}
			std::cout<<index<<std::endl;

			sample_array.push_back(index);

		}




		unsigned int orderX[3],orderY[3],orderZ[3];
		orderX[0] = 1; orderX[1] = 0; orderX[2] = 0;
		orderY[0] = 0; orderY[1] = 1; orderY[2] = 0;
		orderZ[0] = 0; orderY[1] = 0; orderZ[2] = 1;

		GaussianDerivativeImageFunctionType::Pointer derivativeX =	GaussianDerivativeImageFunctionType::New();
		derivativeX->SetInputImage( im1->GetOutput() );
		derivativeX->SetMaximumError( maxError );
		derivativeX->SetMaximumKernelWidth( maxKernelWidth );
		derivativeX->SetOrder( orderX );
		derivativeX->SetNormalizeAcrossScale( false );
		derivativeX->SetUseImageSpacing( true );
		derivativeX->SetInterpolationMode( GaussianDerivativeImageFunctionType::NearestNeighbourInterpolation  );


		GaussianDerivativeImageFunctionType::Pointer derivativeY =	GaussianDerivativeImageFunctionType::New();
		derivativeY->SetInputImage( im1->GetOutput() );
		derivativeY->SetMaximumError( maxError );
		derivativeY->SetMaximumKernelWidth( maxKernelWidth );
		derivativeY->SetOrder( orderY );
		derivativeY->SetNormalizeAcrossScale( false );
		derivativeY->SetUseImageSpacing( true );
		derivativeY->SetInterpolationMode( GaussianDerivativeImageFunctionType::NearestNeighbourInterpolation );


		GaussianDerivativeImageFunctionType::Pointer derivativeZ =	GaussianDerivativeImageFunctionType::New();
		derivativeZ->SetInputImage( im1->GetOutput() );
		derivativeZ->SetMaximumError( maxError );
		derivativeZ->SetMaximumKernelWidth( maxKernelWidth );
		derivativeZ->SetOrder( orderZ );
		derivativeZ->SetNormalizeAcrossScale( false );
		derivativeZ->SetUseImageSpacing( true );
		derivativeZ->SetInterpolationMode( GaussianDerivativeImageFunctionType::NearestNeighbourInterpolation  );


		DiscreteHessianGaussianImageFunctionType::Pointer hessianfunction = DiscreteHessianGaussianImageFunctionType::New();
		DiscreteHessianGaussianImageFunctionType::TensorType hessian;
		hessianfunction->SetInputImage( im1->GetOutput() );
		hessianfunction->SetMaximumError( maxError );
		hessianfunction->SetMaximumKernelWidth( maxKernelWidth );
		hessianfunction->SetNormalizeAcrossScale( false );
		hessianfunction->SetUseImageSpacing( true );
		hessianfunction->SetInterpolationMode( DiscreteHessianGaussianImageFunctionType::NearestNeighbourInterpolation );



		typedef itk::DiscreteGradientMagnitudeGaussianImageFunction< ImageType, PixelType > DiscreteGradientMagnitudeGaussianFunctionType;
		DiscreteGradientMagnitudeGaussianFunctionType::Pointer gradientmagnitudefunction = DiscreteGradientMagnitudeGaussianFunctionType::New();
		gradientmagnitudefunction->SetInputImage( im1->GetOutput() );
		gradientmagnitudefunction->SetMaximumError( maxError );
		gradientmagnitudefunction->SetMaximumKernelWidth( maxKernelWidth );
		gradientmagnitudefunction->SetNormalizeAcrossScale( true );
		gradientmagnitudefunction->SetUseImageSpacing( true );
		gradientmagnitudefunction->SetInterpolationMode( DiscreteGradientMagnitudeGaussianFunctionType::NearestNeighbourInterpolation);

		GaussianBlurFunctionType::Pointer gaussianblur = GaussianBlurFunctionType::New();

		gaussianblur->SetInputImage(im1->GetOutput());
		gaussianblur->SetMaximumError(maxError);
		gaussianblur->SetMaximumKernelWidth(maxKernelWidth);
	//	gaussianblur->SetNormalizeAcrossScale( false);
		gaussianblur->SetUseImageSpacing(true);
	//	gaussianblur->SetInterpolationMode(GaussianBlurImageFunction::NearestNeighbourInterpolation);



			for(int scale=0;scale<3;scale++)
			{
				ImageType::IndexType index;

				gradientmagnitudefunction->SetSigma( sigma[scale] );
				derivativeX->SetSigma( sigma[scale] );
				derivativeY->SetSigma( sigma[scale] );
				derivativeZ->SetSigma( sigma[scale] );
				hessianfunction->SetSigma( sigma[scale] );
				gaussianblur->SetSigma(sigma[scale]);

				hessianfunction->Initialize( );
				derivativeX->Initialize();
				derivativeY->Initialize();
				derivativeZ->Initialize();
				gradientmagnitudefunction->Initialize();
		//		gaussianblur->Initialize();



				int arrayoffset=subj*sample_num*2;

				int samples2= sample_array.size();
				int sample_1= samples2 /2 ;

				for(int ii=0;ii<samples2;ii++)
				{
					float xx = derivativeX->EvaluateAtIndex(sample_array[ii]);
					float yy = derivativeY->EvaluateAtIndex(sample_array[ii]);
					float zz = derivativeZ->EvaluateAtIndex(sample_array[ii]);
					hessian = hessianfunction->EvaluateAtIndex( sample_array[ii] );
					float gm = gradientmagnitudefunction->EvaluateAtIndex( sample_array[ii] );

					int temp = arrayoffset+ii;
					int temp2= scale * 11;

					
					
					features[temp][temp2] = xx;
					features[temp][1+temp2] = yy;
					features[temp][2+temp2] = zz;
					
					for (int kk=0;kk<6;kk++)
						features[temp][3+kk+temp2] = hessian[kk];

					features[temp][9+temp2] = gm;
					features[temp][10+temp2] = gaussianblur->EvaluateAtIndex(sample_array[ii]);
				
					if(scale==0)
					{
						features[temp][33] = sample_array[ii][0];
						features[temp][34] = sample_array[ii][1];
						features[temp][35] = sample_array[ii][2];
						
						if(ii<sample_1)
							Y[temp]=1;
					}

				}
				
			}
	}

	X.reserve(subjectN*sample_num*2);


	for(int i=0;i<subjectN*sample_num*2;i++)
		for(int j=0;j<NFeature;j++)
			X.push_back(features[i][j]);

	for(int i=0;i<subjectN*sample_num*2;i++)
		delete []features[i];

	delete []features;

	char fileName[1024];
	sprintf (fileName, "%s-AdaBoostResults",argv[5]);
	ofstream AdaBoostParamFile(fileName, ios::out);

	AdaBoostParamFile.close();
	sprintf (fileName, "%s-AdaBoostResults",argv[5]);
	AdaBoostTrain(&X[0],&Y[0],subjectN*sample_num*2,NFeature,iterN,fileName);
	return 0;
}
