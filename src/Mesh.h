/*
 * Mesh.h
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

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include "Log.h"

namespace Moirai {

class Mesh {

public:

	typedef int LO;
	typedef int GO;
	typedef double ST;

	struct Node {

		Node() {}
		Node(const double _x, const double _y, const double _z) {
			x[0] = _x;
			x[1] = _y;
			x[2] = _z;
		}
		ST x[3];
		GO id;

	};

	struct Cell {
		std::vector<LO> node_ids;
	};

	struct SubCell {
		SubCell(LO cell_id,int ordinal):cell_id(cell_id),ordinal(ordinal) {}
		LO cell_id;
		int ordinal;
		std::vector<LO> face_node_ids;
		std::vector<LO> cell_node_ids;
	};

	enum CellType {Hex, Tet};

	Mesh();
	void initialise(const double dx);
	std::vector<LO>& get_local_boundary_node_ids();
	std::vector<LO>& get_global_boundary_node_ids();
	std::vector<LO>& get_global_interior_node_ids();
	std::vector<Cell>& get_cells();
	std::vector<Node>& get_nodes();
	std::vector<SubCell>& get_boundary_faces();
	LO get_num_interior_nodes();
	LO get_num_boundary_nodes();
	LO get_num_cells();
	LO get_num_boundary_faces();
	CellType get_cell_type();
	int get_nodes_per_cell();
	int get_nodes_per_face();
	int get_space_dim();
	bool is_in(double x, double y, double z);
	GO get_nearest_node(double x, double y, double z);

private:
	std::vector<Cell> cells;
	std::vector<SubCell> faces;
	std::vector<Node> nodes;
	std::vector<LO> local_boundary_ids;
	std::vector<GO> global_boundary_ids;
	std::vector<GO> global_interior_ids;

	int num_nodes_per_cell;
	int num_nodes_per_face;
	int space_dim;
	CellType cell_type;
	ST max[3];
	ST min[3];
};



} /* namespace Moirai */
#endif /* MESH_H_ */
