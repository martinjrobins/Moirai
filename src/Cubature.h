/*
 * Cubature.h
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
 *  Created on: 22 Mar 2013
 *      Author: robinsonm
 */

#ifndef CUBATURE_H_
#define CUBATURE_H_

#include "Intrepid_Cubature.hpp"
#include "Intrepid_Basis.hpp"
#include "Types.h"

#include "Mesh.h"

namespace Moirai {

shards::CellTopology get_cell_type(Mesh& mesh);

shards::CellTopology get_face_type(Mesh& mesh);


RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > get_basis(Mesh& mesh);

RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > get_face_basis(Mesh& mesh);

} /* namespace Moirai */
#endif /* CUBATURE_H_ */
