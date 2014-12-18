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
#ifndef __itkBox2SumImageFilter_h
#define __itkBox2SumImageFilter_h

#include "itkBoxImageFilter.h"

namespace itk
{
/**
 * \class BoxMeanImageFilter
 * \brief Implements a fast rectangular mean filter using the
 * accumulator approach
 *
 *
 * This code was contributed in the Insight Journal paper:
 * "Efficient implementation of kernel filtering"
 * by Beare R., Lehmann G
 * http://hdl.handle.net/1926/555
 * http://www.insight-journal.org/browse/publication/160
 *
 *
 * \author Richard Beare
 * \ingroup ITKReview
 */

template< class TInputImage, class TOutputImage = TInputImage >
class ITK_EXPORT Box2SumImageFilter:
  public BoxImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef Box2SumImageFilter                           Self;
  typedef BoxImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Standard New method. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(Box2SumImageFilter,
               BoxImageFilter);

  /** Image related typedefs. */
  typedef TInputImage                                InputImageType;
  typedef TOutputImage                               OutputImageType;
  typedef typename TInputImage::RegionType           RegionType;
  typedef typename TInputImage::SizeType             SizeType;
  typedef typename TInputImage::IndexType            IndexType;
  typedef typename TInputImage::PixelType            PixelType;
  typedef typename TInputImage::OffsetType           OffsetType;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename TOutputImage::PixelType           OutputPixelType;

  /** Image related typedefs. */
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( SameDimension,
                   ( Concept::SameDimension< itkGetStaticConstMacro(InputImageDimension),
                                             itkGetStaticConstMacro(OutputImageDimension) > ) );

  /** End concept checking */
#endif
protected:
  Box2SumImageFilter();
  ~Box2SumImageFilter() {}

  /** Multi-thread version GenerateData. */
  void  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId);

private:
  Box2SumImageFilter(const Self &);  //purposely not implemented
  void operator=(const Self &);     //purposely not implemented
};                                  // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBox2SumImageFilter.hxx"
#endif

#endif