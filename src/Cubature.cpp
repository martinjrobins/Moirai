/*
 * Cubature.cpp
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

#include "Cubature.h"
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>

namespace Moirai {

shards::CellTopology get_cell_type(Mesh& mesh) {
	shards::CellTopology cellType;

	switch (mesh.get_cell_type()) {
	case Mesh::Tet:
		cellType = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> > ());
		break;

	case Mesh::Hex:
		cellType = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> > ());
		break;

	default:
		ERROR("Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4.");
	}

	return cellType;
}

shards::CellTopology get_face_type(Mesh& mesh) {
	shards::CellTopology cellType;

	switch (mesh.get_cell_type()) {
	case Mesh::Tet:
		cellType = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> > ());
		break;

	case Mesh::Hex:
		cellType = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> > ());
		break;

	default:
		ERROR("Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4.");
	}

	return cellType;
}

RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > get_basis(Mesh& mesh) {
	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis;
	switch (mesh.get_cell_type()) {
	case Mesh::Tet:
		HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TET_C1_FEM<ST, Intrepid::FieldContainer<ST> > );
		break;

	case Mesh::Hex:
		HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_HEX_C1_FEM<ST, Intrepid::FieldContainer<ST> > );
		break;

	default:
		ERROR("Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4.");
	}

	return HGradBasis;
}

RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > get_face_basis(Mesh& mesh) {
	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis;
	switch (mesh.get_cell_type()) {
	case Mesh::Tet:
		HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_TRI_C1_FEM<ST, Intrepid::FieldContainer<ST> > );
		break;

	case Mesh::Hex:
		HGradBasis = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<ST, Intrepid::FieldContainer<ST> > );
		break;

	default:
		ERROR("Unknown cell topology for basis selction. Please use Hexahedron_8 or Tetrahedron_4.");
	}

	return HGradBasis;
}

} /* namespace Moirai */
