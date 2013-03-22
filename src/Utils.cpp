/*
 * Utils.cpp
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

#include "Utils.h"
#include <algorithm>

namespace Moirai {

RCP<sparse_matrix_type> construct_saddle_point_matrix(RCP<sparse_matrix_type> A, RCP<sparse_matrix_type> B) {

	RCP<const map_type> a_map = A->getRowMap();
	RCP<const map_type> b_map = B->getRowMap();

	const int num_a = a_map->getNodeNumElements();
	const int num_b = b_map->getNodeNumElements();

	Teuchos::ArrayView<const GO> a_elem_view = a_map->getNodeElementList();
	Teuchos::ArrayView<const GO> b_elem_view = b_map->getNodeElementList();

	std::vector<GO> new_map_elements(num_a+num_b);
	std::copy( a_elem_view.begin(), a_elem_view.begin()+num_a, new_map_elements.begin());
	std::copy( b_elem_view.begin(), b_elem_view.begin()+num_b, new_map_elements.begin() + num_a);

	RCP<const map_type> new_map = Tpetra::createNonContigMapWithNode<LO, GO, Node>(
				Teuchos::arrayViewFromVector(new_map_elements),a_map->getComm(),a_map->getNode());

	RCP<sparse_matrix_type> saddle =
			Tpetra::createCrsMatrix<ST,LO,GO,Node>(new_map,A->getCrsGraph()->getNodeMaxNumRowEntries() +
														   B->getCrsGraph()->getNodeMaxNumRowEntries());

	for (int row = 0; row < num_a; ++row) {
		/*
		 * copy A to saddle
		 */
		Teuchos::ArrayView<const LO> indicies;
		Teuchos::ArrayView<const ST> values;

		A->getLocalRowView(row,indicies,values);

		saddle->insertLocalValues(row,indicies,values);
	}

	for (int row = 0; row < num_b; ++row) {
		/*
		 * copy B to saddle
		 */
		Teuchos::ArrayView<const LO> indicies;
		Teuchos::ArrayView<const ST> values;
		size_t num;

		B->getLocalRowView(row,indicies,values);
		saddle->insertLocalValues(row+num_a,indicies,values);

		/*
		 * copy B^t to saddle
		 */
		const int numEntries = B->getNumEntriesInLocalRow(row);
		for (int col = 0; col < numEntries; ++col) {
			const LO index = indicies[col]+num_a;
			const ST value = values[col];
			saddle->insertLocalValues(row,Teuchos::arrayView(&index,1),Teuchos::arrayView(&value,1));
		}
	}

	saddle->fillComplete();
	return saddle;
}

RCP<sparse_matrix_type> get11(RCP<sparse_matrix_type> A) {

	RCP<const map_type> a_map = A->getRowMap();
	const int num_a = a_map->getNodeNumElements();

	RCP<sparse_matrix_type> A11 =
			Tpetra::createCrsMatrix<ST,LO,GO,Node>(a_map,A->getCrsGraph()->getNodeMaxNumRowEntries());

	for (int row = num_a/2; row < num_a; ++row) {
		/*
		 * copy lower right half of A to A11
		 */
		Teuchos::ArrayView<const LO> indicies;
		Teuchos::ArrayView<const ST> values;

		Teuchos::Array<LO> new_indicies;
		Teuchos::Array<ST> new_values;

		const int numEntries = A->getNumEntriesInLocalRow(row);
		for (int col = 0; col < numEntries; ++col) {
			if (indicies[col] > num_a/2) {
				new_indicies.push_back(indicies[col]);
				new_values.push_back(values[col]);
			}
		}

		A11->insertLocalValues(row,new_indicies(),new_values());
	}


	A11->fillComplete();
	return A11;
}

RCP<sparse_matrix_type> scatter_columns(RCP<sparse_matrix_type> A,std::vector<LO>& new_columns) {

}

}

