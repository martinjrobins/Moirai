/*
 * Integration.h
 * 
 * Copyright 2013 Martin Robinson
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
 *  Created on: 20 Mar 2013
 *      Author: robinsonm
 */

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include "Mesh.h"
#include "Types.h"
#include <functional>

namespace Moirai {

RCP<sparse_matrix_type> construct_stiffness_matrix(Mesh& mesh);

RCP<sparse_matrix_type> construct_stiffness_matrix(Mesh& mesh, std::function<ST(ST)> function);

RCP<sparse_matrix_type> construct_mass_matrix(Mesh& mesh);

RCP<sparse_matrix_type> construct_boundary_mass_matrix(Mesh& mesh);

} /* namespace Moirai */
#endif /* INTEGRATION_H_ */
