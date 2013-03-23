/* 
 * Pde.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of Moirai.
 *
 * Moirai is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Moirai is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Moirai.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Mar 22, 2013
 *      Author: mrobins
 */

#ifndef PDE_H_
#define PDE_H_

#include "Mesh.h"
#include "Types.h"

#include <vector>
#include <tuple>

namespace Moirai {

class Pde {
public:
	Pde(const double dx, const double dt);
	void timestep();
	void add_particle(const double x, const double y, const double z);
	void timestep_and_generate_particles(std::vector<std::tuple<double,double,double> > positions);
private:
	double dt;
	static const double omega = 1.0;
	Mesh mesh;

	RCP<vector_type> X,Y,u,lambda;
	RCP<sparse_matrix_type> LHS,RHS;

	void initialise_from_mesh();
};

} /* namespace Moirai */
#endif /* PDE_H_ */
