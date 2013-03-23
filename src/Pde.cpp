/* 
 * Pde.cpp
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

#include "Pde.h"
#include "zip.h"
#include "Integration.h"
#include "Utils.h"
#include <algorithm>



namespace Moirai {

Pde::Pde(const double dx, const double dt):dt(dt) {
	mesh.initialise(dx);

	initialise_from_mesh();
}

void Pde::timestep() {
}

void Pde::add_particle(const double x, const double y, const double z) {
}

void Pde::timestep_and_generate_particles(
		std::vector<std::tuple<double, double, double> > positions) {
}

void Pde::timestep_and_generate_particles(
		std::vector<std::tuple<double, double, double> > positions) {
}

void Pde::initialise_from_mesh() {
	std::vector<GO>& node_ids = mesh.get_global_interior_node_ids();
	std::vector<GO>& bnode_ids = mesh.get_global_boundary_node_ids();
	std::vector<GO> all_node_ids(node_ids.size()+bnode_ids.size());
	std::copy(node_ids.begin(),node_ids.end(),all_node_ids.begin());
	const int ni = node_ids.size();
	const int nb = bnode_ids.size();
	for (int i = 0; i < nb; ++i) {
		all_node_ids[i+ni] = bnode_ids[i]+ni;
	}

	RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
	RCP<Node> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
	RCP<const map_type> interiorMap = Tpetra::createNonContigMapWithNode(Teuchos::arrayViewFromVector(node_ids),comm,node);
	RCP<const map_type> boundaryMap = Tpetra::createNonContigMapWithNode(Teuchos::arrayViewFromVector(bnode_ids),comm,node);
	RCP<const map_type> allMap = Tpetra::createNonContigMapWithNode(Teuchos::arrayViewFromVector(all_node_ids),comm,node);

	X = Tpetra::createVector(allMap);
	X->putScalar(0);
	Y = Tpetra::createVector(allMap);
	u = X->offsetViewNonConst(interiorMap,0)->getVectorNonConst(0);
	lambda = X->offsetViewNonConst(boundaryMap,ni)->getVectorNonConst(0);

	RCP<sparse_matrix_type> K = construct_stiffness_matrix(mesh);
	RCP<sparse_matrix_type> Mi = construct_mass_matrix(mesh);
	RCP<sparse_matrix_type> A = rcp(new Tpetra::CrsMatrix(K->getCrsGraph()));
	const int nrow = K->getNodeNumRows();
	for (int row = 0; row < nrow; ++row) {
		Teuchos::ArrayView<const LO> indicies;
		Teuchos::ArrayView<const ST> values_K;
		Teuchos::ArrayView<const ST> values_Mi;


		K->getLocalRowView(row,indicies,values_K);
		Mi->getLocalRowView(row,indicies,values_Mi);

		const int n = indicies.size();
		Teuchos::Array<ST> new_values(n);
		for (int col = 0; col < n; ++col) {
			new_values[col] = values_K[col] + omega*values_Mi[col];
		}

		A->replaceLocalValues(row,indicies,new_values);
	}

	RCP<sparse_matrix_type> Mb = construct_boundary_mass_matrix(mesh);
	RCP<sparse_matrix_type> B = scatter_columns(Mb,mesh.get_boundary_to_interior_ids());
	B->scale(omega);

}


void Pde::timestep() {
}
 /* namespace Moirai */
