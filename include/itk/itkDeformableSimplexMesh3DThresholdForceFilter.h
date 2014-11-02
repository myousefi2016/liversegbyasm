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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef __itkDeformableSimplexMesh3DThresholdForceFilter_h
#define __itkDeformableSimplexMesh3DThresholdForceFilter_h

#include "itkDeformableSimplexMesh3DFilter.h"
#include "itkMesh.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCovariantVector.h"
#include "itkLinearInterpolateImageFunction.h"

#include <set>

namespace itk
{
/** \class DeformableSimplexMesh3DThresholdForceFilter
  * \brief
  * Additional to its superclass this model adds an balloon force component to the
  * internal forces.
  *
  * The balloon force can be scaled, by setting the parameter kappa.
  *
  * \author Thomas Boettger. Division Medical and Biological Informatics, German Cancer Research Center, Heidelberg.
  *
  * \ingroup ITKDeformableMesh
  */
template< class TInputMesh, class TOutputMesh >
class ITK_EXPORT DeformableSimplexMesh3DThresholdForceFilter:public DeformableSimplexMesh3DFilter< TInputMesh,
                                                                                                 TOutputMesh >
{
public:
  /** Standard "Self" typedef. */
  typedef DeformableSimplexMesh3DThresholdForceFilter Self;

  /** Standard "Superclass" typedef. */
  typedef DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh > Superclass;

  /** Smart pointer typedef support */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method of creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DeformableSimplexMesh3DThresholdForceFilter, DeformableSimplexMesh3DFilter);

  /** Some typedefs. */
  typedef TInputMesh                                  InputMeshType;
  typedef TOutputMesh                                 OutputMeshType;
  typedef typename Superclass::PointType              PointType;
  typedef typename Superclass::GradientIndexType      GradientIndexType;
  typedef typename Superclass::GradientIndexValueType GradientIndexValueType;
  
  typedef typename Superclass::GradientImageType::PixelType GradientPixelType;

  /* Mesh pointer definition. */
  typedef typename InputMeshType::Pointer  InputMeshPointer;
  typedef typename OutputMeshType::Pointer OutputMeshPointer;

  typedef typename InputMeshType::PixelType PixelType;

	typedef Image< float, 3 >                  FeatureImageType;
	typedef typename FeatureImageType::Pointer FeatureImagePointer;

	typedef itk::LinearInterpolateImageFunction<GradientImageType, double> GradientInterpolatorType;
	typedef typename GradientInterpolatorType::Pointer                     GradientInterpolatorPointer;

	typedef itk::LinearInterpolateImageFunction<FeatureImageType, double>  FeatureInterpolatorType;
	typedef typename FeatureInterpolatorType::Pointer                      FeatureInterpolatorPointer;

    itkSetMacro(Kappa, double);
    itkGetConstMacro(Kappa, double);

	itkSetMacro(Threshold, double);
	itkGetConstMacro(Threshold, double);

	itkSetObjectMacro( FeatureImage, FeatureImageType );
	itkGetObjectMacro( FeatureImage, FeatureImageType );

protected:
  DeformableSimplexMesh3DThresholdForceFilter();
  ~DeformableSimplexMesh3DThresholdForceFilter();
  DeformableSimplexMesh3DThresholdForceFilter(const Self &)
  {}

  void operator=(const Self &)
  {}

  void PrintSelf(std::ostream & os, Indent indent) const;

  /**
   * Compute the external force component
   */
  virtual void ComputeExternalForce(SimplexMeshGeometry *data,const GradientImageType *gradientImage);

	virtual void Initialize();

  /** Parameters definitions. */

  /**
   * scalar for balloon force
   */
  double m_Kappa;

  FeatureImagePointer m_FeatureImage;

  double m_Threshold;

  GradientInterpolatorPointer m_GradientInterpolator;
  FeatureInterpolatorPointer  m_FeatureInterpolator;

}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformableSimplexMesh3DThresholdForceFilter.hxx"
#endif

#endif //__itkDeformableSimplexMesh3DThresholdForceFilter_H
