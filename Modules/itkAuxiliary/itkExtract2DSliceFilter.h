#ifndef __itkExtract2DSliceFilter_h
#define __itkExtract2DSliceFilter_h

#include "itkExtractImageFilter.h"

namespace itk
{
	template< typename TInputImage, typename TOutputImage >
	class ITK_EXPORT Extract2DSliceFilter:
		public ExtractImageFilter< TInputImage, TOutputImage >
	{
	public:
		/** Standard class typedefs. */
		typedef Extract2DSliceFilter                    Self;
		typedef ExtractImageFilter< TInputImage, TOutputImage > Superclass;
		typedef SmartPointer< Self >                            Pointer;
		typedef SmartPointer< const Self >                      ConstPointer;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods) */
		itkTypeMacro(Extract2DSliceFilter, ImageToImageFilter);

		//void SetInput(TInputImage *inputImage);

		itkSetMacro(ExtractDimension,  unsigned int);
		itkSetMacro(ExtractSliceIndex, int);
		
		itkGetConstReferenceMacro(ExtractDimension,  unsigned int);
		itkGetConstReferenceMacro(ExtractSliceIndex, int);
		
		void SetDefaultOn();

	protected:
		Extract2DSliceFilter();

		virtual ~Extract2DSliceFilter() {}

		virtual void PresetExtractRegion();

		virtual void GenerateOutputInformation();

		void PrintSelf(std::ostream & os, Indent indent) const;
	private:
		Extract2DSliceFilter(const Self &); //purposely not implemented
		void operator=(const Self &);               //purposely not implemented

		unsigned int m_ExtractDimension;
		int          m_ExtractSliceIndex;
		bool         m_DefaultOn;
	};  //end of class
} //end of km

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkExtract2DSliceFilter.hxx"
#endif

#endif