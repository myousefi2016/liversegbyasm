/***************************************************************************
vtkIsotropicDiscreteRemeshing.h  -  description
-------------------
begin                : March 2006
copyright            : Sebastien Valette
email                : 
***************************************************************************/

/* ---------------------------------------------------------------------

* Copyright (c) CREATIS-LRMN (Centre de Recherche en Imagerie Medicale)
* Author : Sebastien Valette
*
*  This software is governed by the CeCILL-B license under French law and 
*  abiding by the rules of distribution of free software. You can  use, 
*  modify and/ or redistribute the software under the terms of the CeCILL-B 
*  license as circulated by CEA, CNRS and INRIA at the following URL 
*  http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html 
*  or in the file LICENSE.txt.
*
*  As a counterpart to the access to the source code and  rights to copy,
*  modify and redistribute granted by the license, users are provided only
*  with a limited warranty  and the software's author,  the holder of the
*  economic rights,  and the successive licensors  have only  limited
*  liability. 
*
*  The fact that you are presently reading this means that you have had
*  knowledge of the CeCILL-B license and that you accept its terms.
* ------------------------------------------------------------------------ */  

#ifndef _VTKANISOTROPICDISCRETEREMESHING_H_
#define _VTKANISOTROPICDISCRETEREMESHING_H_

#include <vtkObjectFactory.h>

#include "vtkDiscreteRemeshing.h"
#include "vtkVerticesProcessing.h"
#include "vtkQuadricAnisotropicMetricForClustering.h"

typedef vtkDiscreteRemeshing<vtkQuadricAnisotropicMetricForClustering> QuadricAnisotropicRemeshing;

/**
 * A Class to process anisotropic coarsening of vtkSurface PolyData.
 * Implemented from the paper: S. Valette, J. M. Chassery, and R. Prost, 
 * "Generic remeshing of 3D triangular meshes with metric-dependent discrete Voronoi Diagrams",
 * IEEE Trans Visu Comp Grap, vol. 14, no. 2, pp. 369–381, 2008
 */

class VTK_EXPORT vtkAnisotropicDiscreteRemeshing : public vtkVerticesProcessing<QuadricAnisotropicRemeshing>
{
public:

	static vtkAnisotropicDiscreteRemeshing* New()
	{
		// First try to	create the object from the vtkObjectFactory
		vtkObject* ret = vtkObjectFactory::CreateInstance("vtkAnisotropicDiscreteRemeshing");
		if(ret)
		{
			return (vtkAnisotropicDiscreteRemeshing*)ret;
		}
		// If the factory was unable to	create the object, then	create it here.
		return (new	vtkAnisotropicDiscreteRemeshing);
	}

protected:

	vtkAnisotropicDiscreteRemeshing()
	{
	}
	~vtkAnisotropicDiscreteRemeshing() {};
};

#endif
