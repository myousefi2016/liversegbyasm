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
#ifndef __itkAdaptSimplexMesh3DFilter_h
#define __itkAdaptSimplexMesh3DFilter_h

#include "itkDeformableSimplexMesh3DFilter.h"

#include <set>

namespace itk
{
	template< typename TInputMesh, typename TOutputMesh >
	class AdaptSimplexMesh3DFilter:public DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh >
	{
	public:
		/** Standard "Self" typedef. */
		typedef AdaptSimplexMesh3DFilter Self;

		/** Standard "Superclass" typedef. */
		typedef DeformableSimplexMesh3DFilter< TInputMesh, TOutputMesh > Superclass;

		/** Smart pointer typedef support */
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		/** Method of creation through the object factory. */
		itkNewMacro(Self);

		/** Run-time type information (and related methods). */
		itkTypeMacro(AdaptSimplexMesh3DFilter, DeformableSimplexMesh3DFilter);

		itkSetMacro(NormalFactor, double);
		itkGetConstMacro(NormalFactor, double);

		itkSetMacro(TangentFactor, double);
		itkGetConstMacro(TangentFactor, double);

	protected:
		AdaptSimplexMesh3DFilter();
		~AdaptSimplexMesh3DFilter();
		AdaptSimplexMesh3DFilter(const Self &); //purposely not implemented
		void operator=(const Self &);                //purposely not implemented

		virtual void GenerateData();

		virtual void ComputeDisplacement();

		virtual void ComputeInternalForce(SimplexMeshGeometry *data);

		double m_NormalFactor;
		double m_TangentFactor;

	}; // end of class
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAdaptSimplexMesh3DFilter.hxx"
#endif

#endif //__itkAdaptSimplexMesh3DFilter_h
