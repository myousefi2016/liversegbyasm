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
#ifndef __itkDeformableSimplexMesh3DFilterWithAdaptedSpeed_h
#define __itkDeformableSimplexMesh3DFilterWithAdaptedSpeed_h

#include "itkMeshToMeshFilter.h"
#include "itkDeformableSimplexMesh3DFilter.h"
#include "itkSimplexMesh.h"
#include "itkSphereSpatialFunction.h"
#include "itkFloodFilledSpatialFunctionConditionalIterator.h"
#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkArray.h"
#include "itkLinearInterpolateImageFunction.h"

#include <set>

namespace itk
{
/** \class DeformableSimplexMesh3DFilterWithAdaptedSpeed
  * \brief Three-dimensional deformable model for image segmentation
  *
  * DeformableSimplexMesh3DFilterWithAdaptedSpeed is a discrete three-dimensional deformable model, which
  * can be used to deform a 3-D SimplexMesh.
  *
  * The mesh deformation is constrained by internal forces. The interal force can be scaled
  * via SetAlpha (typical values are 0.01 < alpha < 0.3). The external force is derived from
  * the image one wants to delineate. Therefore an image of type GradientImageType needs to
  * be set by calling SetGradientImage(...). The external forces are scaled
  * via SetBeta (typical values are 0.01 < beta < 1). One still needs to play around with
  * these values.
  *
  * To control the smoothness of the mesh a rigidity parameter can be adjusted. Low values (1 or 0)
  * allow areas with high curvature. Higher values (around 7 or 8) will make the mesh smoother.
  *
  * By setting the gamma parameter the regularity of the mesh is controlled. Low values (< 0.03)
  * produce more regular mesh. Higher values ( 0.3 < gamma < 0.2) will allow to move the vertices to
  * regions of higher curvature.
  *
  * This approach for segmentation follows that of Delingette et al. (1997).
  *
  * This filter currently assumes that the spacing of the input image is 1.
  *
  * The user has to set the number of iterations for mesh evolution.
  *
  * \author Thomas Boettger. Division Medical and Biological Informatics, German Cancer Research Center, Heidelberg.
  *
  * \ingroup MeshFilters
  * \ingroup MeshSegmentation
  * \ingroup ITKDeformableMesh
  */
template< class TInputMesh, class TOutputMesh >
class ITK_EXPORT DeformableSimplexMesh3DFilterWithAdaptedConfidence:public DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh >
{
public:
  /** Standard "Self" typedef. */
  typedef DeformableSimplexMesh3DFilterWithAdaptedConfidence Self;

  /** Standard "Superclass" typedef. */
  typedef DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh > Superclass;

  /** Smart pointer typedef support */
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method of creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DeformableSimplexMesh3DFilterWithAdaptedConfidence, DeformableSimplexMesh3DFilter);

	typedef typename InputMeshType::PointDataContainer             InputPointDataContainer;
	typedef typename InputMeshType::PointDataContainer::iterator   InputPointDataContainerIterator;

	typedef itk::LinearInterpolateImageFunction<GradientImageType, double> GradientInterpolatorType;
	typedef typename GradientInterpolatorType::Pointer                     GradientInterpolatorPointer;

	itkSetMacro(ConfidenceThreshold, double);
	itkGetConstMacro(ConfidenceThreshold, double);

	GradientInterpolatorType * GetGradientIntepolator(){ return m_GradientInterpolator; }

	unsigned int error_count;

	void SetVarianceMap(InputMeshType* varianceMap)
	{
		this->m_VarianceMap = varianceMap;
	}

protected:
  DeformableSimplexMesh3DFilterWithAdaptedConfidence();
  ~DeformableSimplexMesh3DFilterWithAdaptedConfidence();
  DeformableSimplexMesh3DFilterWithAdaptedConfidence(const Self &); //purposely not implemented
  void operator=(const Self &);                //purposely not implemented

  void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void Initialize();

  virtual void ComputeDisplacement();

  virtual void GenerateData();

  /**
   * Compute the external force component
   * Pass in the gradient image, to avoid inner loop calls to GetGradient()
   */
  virtual void ComputeExternalForce(SimplexMeshGeometry *data, const GradientImageType *gradient, unsigned int indx, double & confidence);

  virtual void Intervene( );

  double m_ConfidenceThreshold;

  GradientInterpolatorPointer m_GradientInterpolator;

  InputMeshPointer m_VarianceMap;

}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDeformableSimplexMesh3DFilterWithAdaptedConfidence.hxx"
#endif

#endif //__itkDeformableSimplexMesh3DFilterWithAdaptedConfidence_h
