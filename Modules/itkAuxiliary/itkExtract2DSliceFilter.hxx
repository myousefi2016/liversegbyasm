#ifndef __itkExtract2DSliceFilter_hxx
#define __itkExtract2DSliceFilter_hxx

#include "itkExtract2DSliceFilter.h"
#include "itkImageRegionConstIterator.h"

#define DEBUG 1

namespace itk
{
	/**
	* Constructor
	*/
	template< class TInputImage, class TOutputImage>
	Extract2DSliceFilter< TInputImage, TOutputImage>
		::Extract2DSliceFilter()
	{
		m_ExtractDimension = 0;
		m_ExtractSliceIndex = 0;
		m_DefaultOn = 1;
		Superclass::SetDirectionCollapseToIdentity();
	}

	template< class TInputImage, class TOutputImage>
	void
		Extract2DSliceFilter< TInputImage, TOutputImage>
		::SetDefaultOn()
	{
		typename Superclass::InputImageRegionType inputRegion;
		typename Superclass::InputImageIndexType inputIndex;
		typename Superclass::InputImageSizeType  inputSize;

		typename Superclass::InputImagePointer inputImage =
			const_cast< TInputImage * >( this->GetInput() );

		inputRegion = inputImage->GetLargestPossibleRegion();
		inputIndex = inputRegion.GetIndex();
		inputSize = inputRegion.GetSize();

		typedef itk::ImageRegionConstIterator<TInputImage> IteratorType;
		IteratorType it(inputImage, inputRegion);
		it.GoToEnd();

		int topSlice = 0;
		while(!it.IsAtBegin()){
			int val = it.Get();
			Superclass::InputImageIndexType ind = it.GetIndex();
			if(val != 0){
				//std::cout<<val<<std::endl;
				//std::cout<<ind<<std::endl;
				topSlice = ind[InputImageDimension - 1];
				break;
			}
			it--;
		}

		m_ExtractDimension = InputImageDimension - 1;
		m_ExtractSliceIndex = topSlice;

		if(DEBUG){
			std::cout<<"m_ExtractDimension = "<<m_ExtractDimension<<std::endl;
			std::cout<<"m_ExtractSliceIndex = "<<m_ExtractSliceIndex<<std::endl;
		}
	}

	//template< class TInputImage, class TOutputImage>
	//void
	//	Extract2DSliceFilter< TInputImage, TOutputImage>
	//	::SetInput(TInputImage *inputImage)
	//{
	//	std::cout<<"Here!"<<std::endl;
	//	TInputImage::Pointer tempInput = TInputImage::New();
	//	tempInput->Graft(inputImage);
	//	Superclass::SetInput(inputImage);
	//}

	template< class TInputImage, class TOutputImage >
	void
		Extract2DSliceFilter< TInputImage, TOutputImage >
		::PresetExtractRegion()
	{
		typename Superclass::InputImageRegionType extractRegion;
		typename Superclass::InputImageIndexType extractIndex;
		typename Superclass::InputImageSizeType  extractSize;

		typename Superclass::InputImagePointer inputImage =
			const_cast< TInputImage * >( this->GetInput() );
		extractIndex = inputImage->GetLargestPossibleRegion().GetIndex();
		extractSize = inputImage->GetLargestPossibleRegion().GetSize();

		extractIndex[m_ExtractDimension] = m_ExtractSliceIndex;
		extractSize[m_ExtractDimension] = 0;

		extractRegion.SetIndex(extractIndex);
		extractRegion.SetSize(extractSize);

		Superclass::SetExtractionRegion(extractRegion);
	}

	template< class TInputImage, class TOutputImage >
	void
		Extract2DSliceFilter< TInputImage, TOutputImage >
		::GenerateOutputInformation()
	{
		PresetExtractRegion();
		Superclass::GenerateOutputInformation();
	}

	template< class TInputImage, class TOutputImage >
	void
		Extract2DSliceFilter< TInputImage, TOutputImage >::PrintSelf(std::ostream & os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);

		os <<indent  << "ExtractDimension: " << m_ExtractDimension << std::endl;
		os <<indent  << "ExtractSliceIndex: " << m_ExtractSliceIndex << std::endl;
	}
} // end namespace km

#endif
