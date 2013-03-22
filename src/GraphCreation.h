/*
 * GraphCreation.h
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

#ifndef GRAPHCREATION_H_
#define GRAPHCREATION_H_

#include "Mesh.h"
#include "Types.h"

namespace Moirai {


RCP<sparse_graph_type>
create_sparse_graph_interior(Mesh& mesh);

RCP<sparse_graph_type>
create_sparse_graph_boundary_interior(Mesh& mesh);

RCP<sparse_graph_type>
create_sparse_graph_boundary_boundary(Mesh& mesh);



} /* namespace Moirai */
#endif /* GRAPHCREATION_H_ */
