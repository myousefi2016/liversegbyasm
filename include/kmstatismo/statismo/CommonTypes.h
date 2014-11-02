/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef __COMMON_TYPES_H
#define __COMMON_TYPES_H

#include "Config.h"
#include "Domain.h"
#include <iostream>

#include <exception>
#include <string>
#include <list>
#include <vector>

#include <Eigen/Dense>

namespace itk
{
	// helper function to compute the hash value of an itk point (needed by unorderd_map)
	template <typename PointType>
	int hash_value(const PointType& pt) {
		int hash_val = 1;
		for (unsigned i = 0; i < pt.GetPointDimension(); i++)
			hash_val *= pt[i];
		return hash_val;
	}
}

namespace statismo {
const double PI	=	3.14159265358979323846;

/// the type that is used for all vector and matrices throughout the library.
typedef float ScalarType;

/// Vector type used throughout the library
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 1> VectorType;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorTypeDoublePrecision;
typedef Eigen::Matrix<ScalarType, 1 , Eigen::Dynamic> RowVectorType;

/// Matrix type used throughout the library
/// Having RowMajor storage is important, as we want to be able to efficiently copy data to HDF5
/// and vnl, which are both storing the data as rowMajor
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixType;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixTypeDoublePrecision;

/// Diagonal matrix type used throughout the library
typedef Eigen::DiagonalMatrix<ScalarType, Eigen::Dynamic> DiagMatrixType;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatrixTypeDoublePrecision;

} //namespace statismo

#endif
