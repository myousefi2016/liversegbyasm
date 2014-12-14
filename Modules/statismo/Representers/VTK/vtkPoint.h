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


#ifndef VTKPOINT_H_
#define VTKPOINT_H_


/**
 * \brief Helper class that represents a vtkPoint.
 * In vtk a point is just of type  T*. The statismo library relies on a proper
 * copy semantics, and hence requires such a wrapper..
 */

class vtkPoint {
public:
	vtkPoint() {
		m_pt[0]  = 0;
		m_pt[1] = 0;
		m_pt[2] = 0;
	}

	vtkPoint(double x, double y, double z)
	{m_pt[0] = x; m_pt[1] = y; m_pt[2] = z; }

	vtkPoint(double* d) { m_pt[0] = d[0]; m_pt[1] = d[1]; m_pt[2] = d[2]; }

	double& operator[](unsigned i) {return m_pt[i];}
	const double& operator[](unsigned i) const {return m_pt[i];}

	const double* data() const {return m_pt; }


	vtkPoint& operator=(const vtkPoint& rhs) {
		if (this != &rhs) {
			m_pt[0] = rhs.m_pt[0];
			m_pt[1] = rhs.m_pt[1];
			m_pt[2] = rhs.m_pt[2];
		}
		return *this;
	}

	vtkPoint(const vtkPoint& orig) {
		operator=(orig);
	}


private:
	double m_pt[3];
};

#endif
