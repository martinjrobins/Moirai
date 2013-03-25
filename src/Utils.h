/*
 * CreateMatricies.h
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

#ifndef CREATEMATRICIES_H_
#define CREATEMATRICIES_H_

#include "Mesh.h"
#include "Types.h"

namespace Moirai {


RCP<sparse_matrix_type> construct_saddle_point_matrix(RCP<const sparse_matrix_type> A, RCP<const sparse_matrix_type> B);

RCP<sparse_matrix_type> construct_anti_symetric_saddle_point_matrix(RCP<const sparse_matrix_type> A, RCP<const sparse_matrix_type> B);

RCP<sparse_matrix_type> get11(RCP<const sparse_matrix_type> A);

RCP<sparse_matrix_type> scatter_columns(RCP<const sparse_matrix_type> A,RCP<const map_type> col_map, const std::vector<LO>& new_columns);

void update(const double alpha, RCP<sparse_matrix_type> A, const double beta, RCP<const sparse_matrix_type> B);

}




#endif /* CREATEMATRICIES_H_ */
