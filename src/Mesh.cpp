/*
 * Mesh.cpp
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

#include "Mesh.h"
#include <Shards_CellTopology.hpp>
// Pamgen includes
#include <create_inline_mesh.h>
#include <im_exodusII_l.h>
#include <im_ne_nemesisI_l.h>
#include <pamgen_extras.h>
#include <iostream>
#include <sstream>
#include <limits>
#include <math.h>


namespace Moirai {

Mesh::Mesh() {
}
void Mesh::initialise(const double dx) {
	shards::CellTopology cellType = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> > ());
	shards::CellTopology faceType = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> > ());
	cell_type = Hex;

	// Get dimensions
	num_nodes_per_cell = cellType.getNodeCount();
	num_nodes_per_face = 4;
	space_dim = cellType.getDimension();

	max[0] = 1.0;
	max[1] = 1.0;
	max[2] = 1.0;
	min[0] = 0.0;
	min[1] = 0.0;
	min[2] = 0.0;

	const int nx = (max[0]-min[0])/dx;
	const int ny = (max[1]-min[1])/dx;
	const int nz = (max[2]-min[2])/dx;

	/**********************************************************************************/
	/******************************* GENERATE MESH ************************************/
	/**********************************************************************************/

	LOG(2,"Generating mesh");
	using std::endl;
	std::ostringstream os;
	os << "mesh" << endl
			<< "\trectilinear" << endl
			<< "\t\tnx = " << nx << endl
			<< "\t\tny = " << ny << endl
			<< "\t\tnz = " << nz << endl
			<< "\t\tbx = 1" << endl
			<< "\t\tby = 1" << endl
			<< "\t\tbz = 1" << endl
			<< "\t\tgmin = 0 0 0" << endl
			<< "\t\tgmax = 1 1 1" << endl
			<< "\tend" << endl
			<< "\tset assign" << endl
			//     << "\t\tsideset, ilo, 1" << endl
			//     << "\t\tsideset, jlo, 2" << endl
			//     << "\t\tsideset, klo, 3" << endl
			<< "\t\tsideset, ihi, 4" << endl
			//     << "\t\tsideset, jhi, 5" << endl
			//     << "\t\tsideset, khi, 6" << endl
			<< "\tend" << endl
			<< "end";

	int error = 0; // Number of errors in generating the mesh

	// Generate mesh with Pamgen
	long long maxInt = 9223372036854775807LL;
	Create_Pamgen_Mesh (os.str().c_str(), space_dim, 0, 1, maxInt);

	// Get mesh size info
	char title[100];
	long long numDim;
	long long numNodes;
	long long numElems;
	long long numElemBlk;
	long long numNodeSets;
	long long numSideSets;
	int id = 0;

	im_ex_get_init_l (id, title, &numDim, &numNodes, &numElems, &numElemBlk,
			&numNodeSets, &numSideSets);

	ASSERT(numElems > 0,"The number of elements in the mesh is zero.");

	long long numNodesGlobal;
	long long numElemsGlobal;
	long long numElemBlkGlobal;
	long long numNodeSetsGlobal;
	long long numSideSetsGlobal;

	im_ne_get_init_global_l (id, &numNodesGlobal, &numElemsGlobal,
			&numElemBlkGlobal, &numNodeSetsGlobal,
			&numSideSetsGlobal);


	LOG(2,"Global number of elements:                     "
			<< numElemsGlobal << std::endl
			<< "Global number of nodes (incl. boundary nodes): "
			<< numNodesGlobal);


	long long * block_ids = new long long [numElemBlk];
	error += im_ex_get_elem_blk_ids_l(id, block_ids);

	long long  *nodes_per_element   = new long long [numElemBlk];
	long long  *element_attributes  = new long long [numElemBlk];
	long long  *elements            = new long long [numElemBlk];
	char      **element_types       = new char * [numElemBlk];
	long long **elmt_node_linkage   = new long long * [numElemBlk];

	for (long long i = 0; i < numElemBlk; ++i) {
		element_types[i] = new char [MAX_STR_LENGTH + 1];
		error += im_ex_get_elem_block_l (id,
				block_ids[i],
				element_types[i],
				(long long*)&(elements[i]),
				(long long*)&(nodes_per_element[i]),
				(long long*)&(element_attributes[i]));
	}

	// connectivity
	for (long long b = 0; b < numElemBlk; ++b) {
		elmt_node_linkage[b] =  new long long [nodes_per_element[b]*elements[b]];
		error += im_ex_get_elem_conn_l (id,block_ids[b], elmt_node_linkage[b]);
	}

	// Get node-element connectivity
	int telct = 0;
	cells.resize(numElems);
	for (long long b = 0; b < numElemBlk; b++) {
		for (long long el = 0; el < elements[b]; el++) {
			for (int j = 0; j < num_nodes_per_cell; ++j) {
				cells[telct].node_ids.push_back(elmt_node_linkage[b][el*num_nodes_per_cell + j]-1);
			}
			++telct;
		}
	}

	// Read node coordinates and place in field container
	nodes.resize(numNodes);
	global_interior_ids.resize(numNodes);
	ST * nodeCoordx = new ST [numNodes];
	ST * nodeCoordy = new ST [numNodes];
	ST * nodeCoordz = new ST [numNodes];
	im_ex_get_coord_l (id, nodeCoordx, nodeCoordy, nodeCoordz);
	for (int i=0; i<numNodes; i++) {
		nodes[i].id = i;
		global_interior_ids[i] = i;
		nodes[i].x[0] = nodeCoordx[i];
		nodes[i].x[1] = nodeCoordy[i];
		nodes[i].x[2] = nodeCoordz[i];
	}
	delete [] nodeCoordx;
	delete [] nodeCoordy;
	delete [] nodeCoordz;



	//
	// Mesh cleanup
	//
	delete [] block_ids;
	block_ids = NULL;
	delete [] nodes_per_element;
	nodes_per_element = NULL;
	delete [] element_attributes;
	element_attributes = NULL;
	for (long long b = 0; b < numElemBlk; ++b) {
		delete [] elmt_node_linkage[b];
		delete [] element_types[b];
	}
	delete [] element_types;
	element_types = NULL;
	delete [] elmt_node_linkage;
	elmt_node_linkage = NULL;
	delete [] elements;
	elements = NULL;

	// Container indicating whether a node is on the boundary (1-yes 0-no)


	//faceOnBoundary.resize(numFaces);

	// Get boundary (side set) information
	long long * sideSetIds = new long long [numSideSets];
	long long * numSidesInSet = new long long [numSideSets];
	long long numDFinSet;
	int numBndyFaces=0;
	im_ex_get_side_set_ids_l(id,sideSetIds);

	int i_boundary_face = 0;
	for (int i = 0; i < numSideSets; ++i) {
		im_ex_get_side_set_param_l (id,sideSetIds[i], &numSidesInSet[i], &numDFinSet);
		if (numSidesInSet[i] > 0) {
			long long * sideSetElemList = new long long [numSidesInSet[i]];
			long long * sideSetSideList = new long long [numSidesInSet[i]];
			im_ex_get_side_set_l (id, sideSetIds[i], sideSetElemList, sideSetSideList);

			i_boundary_face += numSidesInSet[i];
		}
	}

	std::vector<int> is_boundary(nodes.size(),0);
	for (int i = 0; i < numSideSets; ++i) {
		if (numSidesInSet[i] > 0){
			long long * sideSetElemList = new long long [numSidesInSet[i]];
			long long * sideSetSideList = new long long [numSidesInSet[i]];
			im_ex_get_side_set_l (id, sideSetIds[i], sideSetElemList, sideSetSideList);

			for (int j = 0; j < numSidesInSet[i]; ++j) {
				const int iface = sideSetSideList[j]-1;
				const int ielem = sideSetElemList[j]-1;
				faces.push_back(SubCell(ielem,iface));
				for (int ifacenode = 0; ifacenode < num_nodes_per_face; ++ifacenode) {
					const int sideNode = cellType.getNodeMap(2,iface,ifacenode);
					const int local_nodeid = cells[ielem].node_ids[sideNode];
					(faces.end()-1)->cell_node_ids.push_back(local_nodeid);
					is_boundary[local_nodeid] = 1;
				}

			}
			delete [] sideSetElemList;
			delete [] sideSetSideList;
		}
	}
	delete [] sideSetIds;
	delete [] numSidesInSet;


	Delete_Pamgen_Mesh ();

	std::vector<int> boundary_id(nodes.size(),0);
	for (int i = 0; i < is_boundary.size(); ++i) {
		if (is_boundary[i]) {
			boundary_id[i] = local_boundary_ids.size();
			local_boundary_ids.push_back(i);
			global_boundary_ids.push_back(boundary_id[i]);
			global_interior_ids_on_boundary.push_back(i);
		}
	}

	for (int i = 0; i < faces.size(); ++i) {
		for (int j = 0; j < faces[i].cell_node_ids.size(); ++j) {
			LO cell_node_id = faces[i].cell_node_ids[j];
			faces[i].face_node_ids.push_back(boundary_id[cell_node_id]);
		}
	}
}

std::vector<Mesh::LO>& Mesh::get_local_boundary_node_ids() {
	return local_boundary_ids;
}

std::vector<Mesh::GO>& Mesh::get_global_interior_ids_on_boundary() {
	return global_interior_ids_on_boundary;
}

std::vector<Mesh::LO>& Mesh::get_global_boundary_node_ids() {
	return global_boundary_ids;
}

std::vector<Mesh::LO>& Mesh::get_global_interior_node_ids() {
	return global_interior_ids;
}

Mesh::GO Mesh::get_num_interior_nodes() {
	return nodes.size();
}

Mesh::GO Mesh::get_num_boundary_nodes() {
	return local_boundary_ids.size();
}

Mesh::CellType Mesh::get_cell_type() {
	return cell_type;
}

int Mesh::get_nodes_per_cell() {
	return num_nodes_per_cell;
}

int Mesh::get_nodes_per_face() {
	return num_nodes_per_face;
}

bool Mesh::is_in(double x, double y, double z) {
	if (    (x >= min[0]) &&
			(x <= max[0]) &&
			(y <= min[1]) &&
			(y >= max[1]) &&
			(z <= min[2]) &&
			(z >= max[2])     ) {
		return true;
	} else {
		return false;
	}
}

Mesh::LO Mesh::get_nearest_node(double x, double y, double z) {
	const int n = nodes.size();
	int closest_node = -1;
	double closest_distance = std::numeric_limits<ST>::max();
	for (int i=0; i<n; i++) {
		const ST dx = nodes[i].x[0] - x;
		const ST dy = nodes[i].x[1] - y;
		const ST dz = nodes[i].x[2] - z;
		if ((dx < closest_distance) && (dy < closest_distance) && (dz < closest_distance)) {
			const ST r = sqrt(dx*dx + dy*dy + dz*dz);
			if (r < closest_distance) {
				closest_distance = r;
				closest_node = i;
				//					}
			} // if node within particle radius
		} // if node within the square
	} // loop through all nodes

	return closest_node;
}

Mesh::LO Mesh::get_num_cells() {
	return cells.size();
}

Mesh::LO Mesh::get_num_boundary_faces() {
	return faces.size();
}

std::vector<Mesh::Cell>& Mesh::get_cells() {
	return cells;
}

std::vector<Mesh::Node>& Mesh::get_nodes() {
	return nodes;
}

std::vector<Mesh::SubCell>& Mesh::get_boundary_faces() {
	return faces;
}

int Mesh::get_space_dim() {
	return space_dim;
}

} /* namespace Moirai */
