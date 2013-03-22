/*
 * GraphCreation.cpp
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

#include "GraphCreation.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"

#include "Teuchos_as.hpp"
#include "Teuchos_ArrayView.hpp"

namespace Moirai {

RCP<sparse_graph_type> create_sparse_graph_interior(Mesh& mesh) {
	const int numNodes = mesh.get_num_interior_nodes();
	const int numNodesPerElem = mesh.get_nodes_per_cell();
	const std::vector<Mesh::Cell> elems = mesh.get_cells();
	const std::vector<Mesh::GO> node_ids = mesh.get_global_interior_node_ids();
	const int numElems = elems.size();


	RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
	RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

	RCP<const map_type> interiorMap = Tpetra::createLocalMapWithNode<LO,GO,Node>(numNodes,comm,node);

	RCP<sparse_graph_type>	interior_graph = rcp (new sparse_graph_type (interiorMap, 0));

	for (int cell = 0; cell < numElems; ++cell) {
		for (int cellRow = 0; cellRow < numNodesPerElem; cellRow++) {
			LO localRow  = elems[cell].node_ids[cellRow];
			std::vector<LO> local_columns;
			for (int cellCol = 0; cellCol < numNodesPerElem; ++cellCol) {
				local_columns.push_back(elems[cell].node_ids[cellCol]);
			}// *** cell col loop ***
			interior_graph->insertGlobalIndices(localRow,Teuchos::arrayViewFromVector(local_columns));
		}// *** cell row loop ***
	}// *** cell loop **

	interior_graph->fillComplete ();
	return interior_graph;
}

RCP<sparse_graph_type> create_sparse_graph_boundary_interior(
		Mesh& mesh) {

	const int numNodes = mesh.get_num_interior_nodes();
	const int numBoundaryNodes = mesh.get_num_boundary_nodes();
	const int numNodesPerElem = mesh.get_nodes_per_cell();
	const int numNodesPerFace = mesh.get_nodes_per_face();
	const std::vector<Mesh::SubCell> faces = mesh.get_boundary_faces();
	const int numFaces = faces.size();


	RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
	RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

	RCP<const map_type> interiorMap = Tpetra::createLocalMapWithNode<LO, GO, Node>(numNodes,comm,node);
	RCP<const map_type> boundaryMap = Tpetra::createLocalMapWithNode<LO, GO, Node>(numBoundaryNodes,comm,node);

	RCP<sparse_graph_type>	boundary_interior_graph = rcp (new sparse_graph_type (boundaryMap, interiorMap, 0));

	for (int cell = 0; cell < numFaces; ++cell) {
		for (int cellRow = 0; cellRow < numNodesPerFace; cellRow++) {
			LO localRow  = faces[cell].face_node_ids[cellRow];
			std::vector<LO> local_columns;
			for (int cellCol = 0; cellCol < numNodesPerFace; ++cellCol) {
				local_columns.push_back(faces[cell].cell_node_ids[cellCol]);
			}// *** cell col loop ***
			boundary_interior_graph->insertGlobalIndices(localRow,Teuchos::arrayViewFromVector(local_columns));
		}// *** cell row loop ***
	}// *** cell loop **

	boundary_interior_graph->fillComplete(interiorMap,boundaryMap);
	return boundary_interior_graph;
}

RCP<sparse_graph_type> create_sparse_graph_boundary_boundary(
		Mesh& mesh) {
	const int numBoundaryNodes = mesh.get_num_boundary_nodes();
	const int numNodesPerFace = mesh.get_nodes_per_face();
	const std::vector<Mesh::SubCell> faces = mesh.get_boundary_faces();
	const int numFaces = faces.size();

	RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
	RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

	RCP<const map_type> boundaryMap = Tpetra::createLocalMapWithNode<LO, GO, Node>(numBoundaryNodes,comm,node);

	RCP<sparse_graph_type>	boundary_graph = rcp (new sparse_graph_type (boundaryMap, 0));


	for (int cell = 0; cell < numFaces; ++cell) {
		for (int cellRow = 0; cellRow < numNodesPerFace; cellRow++) {
			LO localRow  = faces[cell].face_node_ids[cellRow];
			std::vector<LO> local_columns;
			for (int cellCol = 0; cellCol < numNodesPerFace; ++cellCol) {
				local_columns.push_back(faces[cell].face_node_ids[cellCol]);
			}// *** cell col loop ***
			boundary_graph->insertGlobalIndices(localRow,Teuchos::arrayViewFromVector(local_columns));
		}// *** cell row loop ***
	}// *** cell loop **

	boundary_graph->fillComplete ();
	return boundary_graph;

}

}
