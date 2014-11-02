/*
AdaBoost classification based on Dr. Hongzhi Wang's AdaBoost code.
*/

//#include "util.h"
#include "AdaBoost.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCommand.h"
#include "itkImageMaskSpatialObject.h"
#include "itkResampleImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkComposeImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkStreamingImageFilter.h"
#include "itkImageRegionSplitterMultidimensional.h"
#include "itkPipelineMonitorImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMeanImageFunction.h"
#include "itkMedianImageFunction.h"
#include "itkVarianceImageFunction.h"
#include "itkSumOfSquaresImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGaussianBlurImageFunction.h"
#include "itkGaussianDerivativeImageFunction.h"
#include "itkDiscreteGaussianDerivativeImageFunction.h"
#include "itkDiscreteGradientMagnitudeGaussianImageFunction.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkSobelEdgeDetectionImageFilter.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkBoxMeanImageFilter.h"
#include "itkBoxSigmaImageFilter.h"
#include "itkBox2SumImageFilter.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkDifferenceOfGaussiansGradientImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkSubtractImageFilter.h"
#include "itkRankImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

using namespace std;

namespace AdaSegment
{

	typedef  unsigned char      MaskPixelType;
	typedef itk::Image<unsigned char, 3 > UCImageType;
	typedef itk::Image<signed short, 3> ImageType;
	typedef itk::Image<float, 3> FloatImageType;
	typedef itk::Image< MaskPixelType, 3 > MaskImageType;
	typedef itk::ImageFileReader< MaskImageType >  MaskReaderType;
	typedef itk::ImageFileReader< ImageType >      ReaderType;
	typedef itk::ImageFileWriter< ImageType >      WriterType;
	typedef itk::VectorImage<float, 3>             VectorImageType;
	typedef itk::ImageFileWriter< MaskImageType >  MaskWriterType;
	typedef itk::ImageFileWriter< FloatImageType > FloatWriterType;


	typedef itk::MeanImageFunction< ImageType >   MeanImageFunctionType;
	typedef itk::MedianImageFunction< ImageType > MedianImageFunctionType;
	typedef itk::VarianceImageFunction<ImageType> VarianceImageFunctionType;
	typedef itk::SumOfSquaresImageFunction<ImageType> SumOfSquaresImageFunctionType;
	typedef itk::CentralDifferenceImageFunction<ImageType> CenteralDifferenceImageFunctionType;
	typedef itk::DiscreteGaussianDerivativeImageFunction<ImageType> GaussianDerivativeImageFunctionType;
	typedef itk::DiscreteGradientMagnitudeGaussianImageFunction<ImageType> GradientMagnitudeImageFunctionType;
	typedef itk::DiscreteHessianGaussianImageFunction<ImageType> HessianGaussianImageFunctionType;
	typedef itk::BoxMeanImageFilter<ImageType,FloatImageType> MeanType;
	typedef itk::BoxSigmaImageFilter<ImageType,FloatImageType> SigmaType;
	typedef itk::Box2SumImageFilter<FloatImageType, FloatImageType> SumType;

	typedef itk::CastImageFilter<ImageType, FloatImageType> CastType;

	typedef itk::SobelEdgeDetectionImageFilter<ImageType, FloatImageType> SobelType;
	typedef itk::ZeroCrossingImageFilter<ImageType, FloatImageType> ZeroCrossingType;
	typedef itk::LaplacianRecursiveGaussianImageFilter<FloatImageType,FloatImageType >    LoGFilter;
	typedef itk::HessianRecursiveGaussianImageFilter< ImageType > HessianFilterType;
	typedef itk::Hessian3DToVesselnessMeasureImageFilter< float > VesselnessMeasureFilterType;
	typedef itk::DerivativeImageFilter< ImageType, FloatImageType >  derivativeFilterType;

	typedef itk::DifferenceOfGaussiansGradientImageFilter<FloatImageType,float> DOGFilterType;
	typedef itk::GradientMagnitudeImageFilter<ImageType,FloatImageType>  GradientMagnitudeFilterType;
	typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType,FloatImageType> GMRGFilterType;

	typedef itk::ImageRegionConstIterator< UCImageType > UCImageIteratorType;
	typedef itk::ImageRegionConstIterator< ImageType > ShortImageIteratorType;
	typedef itk::ImageRegionConstIterator< FloatImageType > FloatImageIteratorType;
	typedef itk::ImageRegionIterator< FloatImageType> FloatImageWriteIteratorType;
	typedef itk::ImageRegionIterator< UCImageType > UCImageWriteIteratorType;

	typedef itk::ImageRegionConstIteratorWithIndex< UCImageType > UCIndexIteratorType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType   > ShortIndexIteratorType;
	typedef itk::ImageRegionConstIteratorWithIndex< FloatImageType > FloatIndexIteratorType;
	typedef itk::ImageRegionIteratorWithIndex< UCImageType > UCIndexWriteIteratorType;
	typedef itk::RankImageFilter<ImageType, FloatImageType > RankType;


	typedef itk::ImageFileWriter< FloatImageType > FLOATWriterType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	typedef itk::ImageFileReader< UCImageType > UCReaderType;

	typedef itk::SubtractImageFilter<UCImageType,UCImageType,UCImageType> SubtractFilterType;
	typedef itk::ConstNeighborhoodIterator< FloatImageType > NeighborhoodIteratorType;

	typedef itk::SignedMaurerDistanceMapImageFilter<UCImageType, FloatImageType>  DistanceMapFilterType;
	typedef itk::ImageLinearIteratorWithIndex< FloatImageType > FloatLineIteratorType;
	typedef itk::ImageLinearConstIteratorWithIndex< UCImageType > UCLineIteratorType;



	template<class ProbabilityImage>
	typename ProbabilityImage::Pointer
		adaSegment(const char* inputImage, const char* abdominalMaskImage, const char* AdaBoostOutPutPrefix )
	{
		if(inputImage == NULL || abdominalMaskImage == NULL || AdaBoostOutPutPrefix == NULL)
		{
			std::cerr<<"Wrong input parameters!"<<std::endl;
			return NULL;
		}

		ReaderType::Pointer orgReader = ReaderType::New();
		orgReader->SetFileName(argv[1]);
		try
		{
			orgReader->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		ImageType::Pointer orgimage = orgReader->GetOutput();

		UCReaderType::Pointer abdominalmaskReader = UCReaderType::New();
		abdominalmaskReader->SetFileName(argv[2]);
		try
		{
			abdominalmaskReader->Update();
		}
		catch(itk::ExceptionObject &e)
		{
			std::cout<<e.GetDescription()<<std::endl;
			return -1;
		}

		UCImageType::Pointer abdominalmask = abdominalmaskReader->GetOutput();

		adaSegment<ProbabilityImage, ImageType, UCImageType>( orgimage, abdominalmask, AdaBoostOutPutPrefix );
	}

	template<class ProbabilityImage, class OrigImage, class MaskImage>
	typename ProbabilityImage::Pointer
		adaSegment(OrigImage* inputImage, MaskImage* abdominalMaskImage, const char* AdaBoostOutPutPrefix, const double thresholdLower = -1024, const double thresholdUpper = 1024 )
	{	
		vector<int> sign;
		vector<FLOATTYPE> alpha;
		vector<FLOATTYPE> threshold;
		vector<int> featureID;

		ifstream ifs; 
		ifs.open( AdaBoostOutPutPrefix , ifstream::in );

		FLOATTYPE t;
		int LC;

		while (ifs.good())
		{
			ifs >> t;
			if (!ifs.good())
				break;
			ifs >> t;
			alpha.push_back(t);
			ifs >> t;
			ifs >> t;
			featureID.push_back(int(t));
			ifs >> t;
			sign.push_back(int(t));
			ifs >> t;
			threshold.push_back(t);
			ifs >> t;
		}
		ifs.close();
		LC = alpha.size();
		cout<<"# weak learners: "<<LC<<endl;

		typedef itk::CastImageFilter<OrigImage, ImageType> InputCastImageFilterType;
		InputCastImageFilterType::Pointer inputCaster = InputCastImageFilterType::New();
		inputCaster->SetInput(inputImage);
		inputCaster->Update();	

		ImageType::Pointer image = inputCaster->GetOutput();

		typedef itk::CastImageFilter<MaskImage, UCImageType> MaskCastImageFilterType;
		MaskCastImageFilterType::Pointer maskCaster = MaskCastImageFilterType::New();
		maskCaster->SetInput(abdominalMaskImage);
		maskCaster->Update();	

		UCImageType::Pointer abmask = maskCaster->GetOutput();


		FloatImageType::Pointer nseg = FloatImageType::New();
		nseg->SetRegions(image->GetRequestedRegion());
		nseg->SetSpacing( image->GetSpacing() );
		nseg->SetOrigin( image->GetOrigin() );
		nseg->SetDirection(image->GetDirection());
		nseg->Allocate();

		MeanType::Pointer meanfilter = MeanType::New();
		SigmaType::Pointer sigmafilter = SigmaType::New();
		SumType::Pointer sumfilter = SumType::New();

		HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
		VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();

		CastType::Pointer caster = CastType::New();
		caster->SetInput(image);
		caster->Update();


		std::vector<FloatImageType::Pointer> imagefeatures;

		imagefeatures.push_back(caster->GetOutput());  // intensity 

		double sigma = 4.0;

		for(int ii=0;ii<3;ii++)
		{

			MeanType::Pointer meanfilter = MeanType::New();
			SigmaType::Pointer sigmafilter = SigmaType::New();
			RankType::Pointer rankfilter0 = RankType::New(); // minimum
			RankType::Pointer rankfilter1 = RankType::New(); // median
			RankType::Pointer rankfilter2 = RankType::New(); // maximum

			MeanType::RadiusType radius;
			radius.Fill( ii+1 );

			meanfilter->SetInput(image);
			meanfilter->SetRadius(radius);
			meanfilter->Update();
			imagefeatures.push_back(meanfilter->GetOutput());

			sigmafilter->SetInput(image);
			sigmafilter->SetRadius(radius);
			sigmafilter->Update();

			imagefeatures.push_back(sigmafilter->GetOutput());

			rankfilter0->SetInput(image);
			rankfilter0->SetRadius(radius);
			rankfilter0->SetRank(0);    // minimum value
			rankfilter0->Update();
			imagefeatures.push_back(rankfilter0->GetOutput());

			rankfilter1->SetInput(image);
			rankfilter1->SetRadius(radius);
			rankfilter1->SetRank(0.5);  // median value
			rankfilter1->Update();
			imagefeatures.push_back(rankfilter1->GetOutput());

			rankfilter2->SetInput(image);
			rankfilter2->SetRadius(radius);
			rankfilter2->SetRank(1);    // maximum value
			rankfilter2->Update();
			imagefeatures.push_back(rankfilter2->GetOutput());

		}

		sumfilter->SetInput(caster->GetOutput());
		sumfilter->Update();

		hessianFilter->SetInput( image );
		hessianFilter->SetSigma( sigma );



		// Create and setup a derivative filter
		derivativeFilterType::Pointer derivativeFilter_x = derivativeFilterType::New();
		derivativeFilterType::Pointer derivativeFilter_y = derivativeFilterType::New();
		derivativeFilterType::Pointer derivativeFilter_z = derivativeFilterType::New();

		derivativeFilter_x->SetInput( image );
		derivativeFilter_y->SetInput( image );
		derivativeFilter_z->SetInput( image );

		derivativeFilter_x->SetDirection(0); // "x" axis
		derivativeFilter_y->SetDirection(1); // "y" axis
		derivativeFilter_z->SetDirection(2); // "z" axis

		derivativeFilter_x->Update();
		derivativeFilter_y->Update();
		derivativeFilter_z->Update();

		imagefeatures.push_back(derivativeFilter_x->GetOutput());
		imagefeatures.push_back(derivativeFilter_y->GetOutput());
		imagefeatures.push_back(derivativeFilter_z->GetOutput());

		DOGFilterType::Pointer DOGFilter = DOGFilterType::New();
		DOGFilter->SetInput(caster->GetOutput());
		DOGFilter->SetWidth(4);
		DOGFilter->Update();

		//imagefeatures.push_back(DOGFilter->GetOutput());

		DOGFilterType::TOutputImage::Pointer gradientImage = DOGFilter->GetOutput();
		GradientMagnitudeFilterType::Pointer gradientmagnitudefilter = GradientMagnitudeFilterType::New();

		gradientmagnitudefilter->SetInput(image);
		gradientmagnitudefilter->Update();
		imagefeatures.push_back(gradientmagnitudefilter->GetOutput());

		GMRGFilterType::Pointer gmrgfilter = GMRGFilterType::New();
		gmrgfilter->SetInput(image);
		gmrgfilter->Update();

		imagefeatures.push_back(gmrgfilter->GetOutput());

		UCIndexIteratorType air_skin_it(abmask,abmask->GetLargestPossibleRegion());

		ShortImageIteratorType imageit(image,image->GetLargestPossibleRegion());


		UCImageType::Pointer im = UCImageType::New();

		im->SetRegions(image->GetLargestPossibleRegion());
		im->CopyInformation(image);
		im->Allocate();
		im->FillBuffer(1);

		//SubtractFilterType::Pointer subtracter = SubtractFilterType::New();

		//subtracter->SetInput(0,im);
		//subtracter->SetInput(1,abmask);
		//subtracter->Update();

		UCIndexIteratorType abdominalmask_it(abmask,abmask->GetLargestPossibleRegion());


		abdominalmask_it.GoToBegin();
		air_skin_it.GoToBegin();
		imageit.GoToBegin();


		// location features
		// location features
		FloatImageType::Pointer outputImageX = FloatImageType::New();
		FloatImageType::Pointer outputImageY = FloatImageType::New();
		FloatImageType::Pointer outputImageXdist = FloatImageType::New();
		FloatImageType::Pointer outputImageYdist = FloatImageType::New();
		outputImageX->SetRegions( image->GetLargestPossibleRegion() );
		outputImageY->SetRegions( image->GetLargestPossibleRegion() );
		outputImageXdist->SetRegions( image->GetLargestPossibleRegion() );
		outputImageYdist->SetRegions( image->GetLargestPossibleRegion() );

		outputImageX->CopyInformation( image );
		outputImageY->CopyInformation( image );
		outputImageXdist->CopyInformation( image );
		outputImageYdist->CopyInformation( image );
		outputImageX->Allocate();
		outputImageY->Allocate();
		outputImageXdist->Allocate();
		outputImageYdist->Allocate();

		outputImageXdist->FillBuffer(0);
		outputImageYdist->FillBuffer(0);


		FloatImageWriteIteratorType outputItX( outputImageX, outputImageX->GetLargestPossibleRegion() );
		FloatImageWriteIteratorType outputItY( outputImageY, outputImageY->GetLargestPossibleRegion() );
		FloatLineIteratorType outputItXdist( outputImageXdist, outputImageXdist->GetLargestPossibleRegion() );
		FloatLineIteratorType outputItYdist( outputImageYdist, outputImageYdist->GetLargestPossibleRegion() );

		FloatImageType::IndexType requestedIndexX =outputImageX->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexY =outputImageY->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexXdist = outputImageXdist->GetLargestPossibleRegion().GetIndex();
		FloatImageType::IndexType requestedIndexYdist = outputImageYdist->GetLargestPossibleRegion().GetIndex();

		FloatImageType::SizeType requestedSize = outputImageX->GetLargestPossibleRegion().GetSize();
		for ( outputItX.GoToBegin(),outputItY.GoToBegin() ;
			!outputItX.IsAtEnd(); ++outputItX,++outputItY)
		{
			FloatImageType::IndexType idx = outputItX.GetIndex();
			outputItX.Set(idx[0]);
			outputItY.Set(idx[1]);
		}

		imagefeatures.push_back(outputImageX);
		imagefeatures.push_back(outputImageY);

		DistanceMapFilterType::Pointer distancemapfilter = DistanceMapFilterType::New();
		distancemapfilter->SetInput( abmask );
		distancemapfilter->SetSquaredDistance( false );
		distancemapfilter->SetUseImageSpacing( true  );
		distancemapfilter->SetInsideIsPositive( true );
		distancemapfilter->Update();


		imagefeatures.push_back(distancemapfilter->GetOutput());


		UCLineIteratorType abdominallineit(abmask, abmask->GetLargestPossibleRegion());


		abdominallineit.SetDirection(0);
		outputItXdist.SetDirection(0);
		abdominallineit.GoToBegin();
		outputItXdist.GoToBegin();

		UCImageType::SpacingType spacing = abmask->GetSpacing();

		while(!abdominallineit.IsAtEnd())
		{
			int abdominalleft = 9999;
			abdominallineit.GoToBeginOfLine();

			if(abdominallineit.Get())
				abdominalleft = 0;
			else
				while(!abdominallineit.IsAtEndOfLine())
				{
					if(abdominallineit.Get())
					{
						abdominalleft = abdominallineit.GetIndex()[0];
						break;
					}
					++abdominallineit;
				}


				outputItXdist.GoToBeginOfLine();
				if(abdominalleft!=9999)
					while(!outputItXdist.IsAtEndOfLine())
					{
						float tmpdist = (outputItXdist.GetIndex()[0]-abdominalleft)*spacing[0];
						outputItXdist.Set(tmpdist);
						++outputItXdist;
					}

					abdominallineit.NextLine();
					outputItXdist.NextLine();
		}

		// Y dist

		abdominallineit.SetDirection(1);
		outputItYdist.SetDirection(1);
		abdominallineit.GoToBegin();
		outputItYdist.GoToBegin();

		while(!abdominallineit.IsAtEnd())
		{
			int abdominaltop = 9999;
			abdominallineit.GoToBeginOfLine();

			if(abdominallineit.Get())
			{
				abdominaltop = 0;
			}
			else
				while(!abdominallineit.IsAtEndOfLine())
				{
					if(abdominallineit.Get())
					{
						abdominaltop = abdominallineit.GetIndex()[1];
						break;
					}
					++abdominallineit;
				}


				outputItYdist.GoToBeginOfLine();
				if(abdominaltop!=9999)
					while(!outputItYdist.IsAtEndOfLine())
					{
						float tmpdist = (outputItYdist.GetIndex()[1]-abdominaltop)*spacing[1];
						outputItYdist.Set(tmpdist);
						++outputItYdist;
					}

					abdominallineit.NextLine();
					outputItYdist.NextLine();
		}


		imagefeatures.push_back(outputImageXdist);
		imagefeatures.push_back(outputImageYdist);


		std::vector<FloatImageIteratorType> imagefeaturesiterator;

		int NFeature = imagefeatures.size();

		for(int i=0;i<NFeature;i++)
		{
			FloatImageIteratorType it(imagefeatures[i],imagefeatures[i]->GetRequestedRegion());
			it.GoToBegin();
			imagefeaturesiterator.push_back(it);
		}



		NeighborhoodIteratorType::RadiusType neighborhood_radius;
		neighborhood_radius.Fill(3);

		//NeighborhoodIteratorType it(neighborhood_radius,mean,mean);

		std::vector<NeighborhoodIteratorType> mean_context_iterators;
		std::vector<NeighborhoodIteratorType> sigma_context_iterators;
		std::vector<NeighborhoodIteratorType> minimum_context_iterators;
		std::vector<NeighborhoodIteratorType> maximum_context_iterators;
		std::vector<NeighborhoodIteratorType> median_context_iterators;

		for(int i=0;i<3;i++)
		{
			int t = i*5;
			NeighborhoodIteratorType it1(neighborhood_radius,imagefeatures[t],imagefeatures[t]->GetLargestPossibleRegion());
			mean_context_iterators.push_back(it1);

			NeighborhoodIteratorType it2(neighborhood_radius,imagefeatures[t+1],imagefeatures[t+1]->GetLargestPossibleRegion());
			sigma_context_iterators.push_back(it2);

			NeighborhoodIteratorType it3(neighborhood_radius,imagefeatures[t+2],imagefeatures[t+2]->GetLargestPossibleRegion());
			minimum_context_iterators.push_back(it3);

			NeighborhoodIteratorType it4(neighborhood_radius,imagefeatures[t+3],imagefeatures[t+3]->GetLargestPossibleRegion());
			median_context_iterators.push_back(it4);

			NeighborhoodIteratorType it5(neighborhood_radius,imagefeatures[t+4],imagefeatures[t+4]->GetLargestPossibleRegion());
			maximum_context_iterators.push_back(it5);

		}
		NeighborhoodIteratorType::OffsetType offset1 = { {-3, 0, 0} };
		NeighborhoodIteratorType::OffsetType offset2 = { { 3, 0, 0} };
		NeighborhoodIteratorType::OffsetType offset3 = { { 0,-3, 0} };
		NeighborhoodIteratorType::OffsetType offset4 = { { 0, 3, 0} };
		NeighborhoodIteratorType::OffsetType offset5 = { { 0, 0,-3} };
		NeighborhoodIteratorType::OffsetType offset6 = { { 0, 0, 3} };


		FloatImageWriteIteratorType nsegit(nseg, nseg->GetRequestedRegion());


		for (nsegit.GoToBegin(); !nsegit.IsAtEnd(); ++nsegit)
		{
			nsegit.Set(0);
		}

		FLOATTYPE H,cH;


		std::cout<<"Start to classify by Adboost........."<<std::endl;

		nsegit.GoToBegin();
		while(!nsegit.IsAtEnd())
		{
			double val = inputImage->GetPixel( nsegit.GetIndex() );

			if(abdominalmask_it.Get() && (thresholdLower<=val && thresholdUpper>=val))
			{


				vector<FLOATTYPE> X;


				for(int i=0;i<3;i++)
				{
					X.push_back(mean_context_iterators[i].GetPixel(offset1));
					X.push_back(mean_context_iterators[i].GetPixel(offset2));
					X.push_back(mean_context_iterators[i].GetPixel(offset3));
					X.push_back(mean_context_iterators[i].GetPixel(offset4));
					X.push_back(mean_context_iterators[i].GetPixel(offset5));
					X.push_back(mean_context_iterators[i].GetPixel(offset6));

					X.push_back(sigma_context_iterators[i].GetPixel(offset1));
					X.push_back(sigma_context_iterators[i].GetPixel(offset2));
					X.push_back(sigma_context_iterators[i].GetPixel(offset3));
					X.push_back(sigma_context_iterators[i].GetPixel(offset4));
					X.push_back(sigma_context_iterators[i].GetPixel(offset5));
					X.push_back(sigma_context_iterators[i].GetPixel(offset6));

					X.push_back(minimum_context_iterators[i].GetPixel(offset1));
					X.push_back(minimum_context_iterators[i].GetPixel(offset2));
					X.push_back(minimum_context_iterators[i].GetPixel(offset3));
					X.push_back(minimum_context_iterators[i].GetPixel(offset4));
					X.push_back(minimum_context_iterators[i].GetPixel(offset5));
					X.push_back(minimum_context_iterators[i].GetPixel(offset6));

					X.push_back(median_context_iterators[i].GetPixel(offset1));
					X.push_back(median_context_iterators[i].GetPixel(offset2));
					X.push_back(median_context_iterators[i].GetPixel(offset3));
					X.push_back(median_context_iterators[i].GetPixel(offset4));
					X.push_back(median_context_iterators[i].GetPixel(offset5));
					X.push_back(median_context_iterators[i].GetPixel(offset6));


					X.push_back(maximum_context_iterators[i].GetPixel(offset1));
					X.push_back(maximum_context_iterators[i].GetPixel(offset2));
					X.push_back(maximum_context_iterators[i].GetPixel(offset3));
					X.push_back(maximum_context_iterators[i].GetPixel(offset4));
					X.push_back(maximum_context_iterators[i].GetPixel(offset5));
					X.push_back(maximum_context_iterators[i].GetPixel(offset6));

				}


				for(int i=0;i<NFeature;i++)
				{
					X.push_back(imagefeaturesiterator[i].Get());
				}

				//Apply AdaBoost classifier on the features to obtain classifier response in probability, stored in cH.
				AdaBoostClassify(&X[0], 1, NFeature, LC, &featureID[0], &alpha[0], &sign[0], &threshold[0], &H, &cH);

				nsegit.Set(static_cast<float>(H));
				//if (H>0.5)
				//	nsegit.Set(1);
			}

			++nsegit;	
			++abdominalmask_it;
			for(int i=0;i<NFeature;i++)
				++imagefeaturesiterator[i];

			for(int i=0;i<3;i++)
			{
				++mean_context_iterators[i];
				++sigma_context_iterators[i];
				++minimum_context_iterators[i];
				++maximum_context_iterators[i];
				++median_context_iterators[i];
			}


		}



		/*
		cout<<"saving the segmentation to "<<argv[4]<<endl;
		// write the segmentation into the output file
		FloatWriterType::Pointer writer = FloatWriterType::New();
		writer->SetFileName( argv[4] );
		writer->SetInput(nseg);
		try
		{
		writer->Update();
		}
		catch ( itk::ExceptionObject &err)
		{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return -1;
		}*/

		typedef itk::CastImageFilter<FloatImageType, ProbabilityImage> CastImageFilterType;
		CastImageFilterType::Pointer resultCaster = CastImageFilterType::New();
		resultCaster->SetInput(nseg);
		resultCaster->Update();

		return resultCaster->GetOutput();
	}

}
