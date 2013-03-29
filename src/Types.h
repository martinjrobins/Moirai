/*
 * Types.h
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

#ifndef TYPES_H_
#define TYPES_H_

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Sacado_CacheFad_SFad.hpp>

namespace Moirai {

using Teuchos::RCP;
using Teuchos::rcp;

typedef double ST;
typedef int    LO;
typedef int    GO;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
typedef Tpetra::CrsMatrix<ST, LO, GO, Node>    sparse_matrix_type;
typedef Tpetra::Operator<ST, LO, GO, Node>     operator_type;
typedef Tpetra::MultiVector<ST, LO, GO, Node>  multivector_type;

typedef Tpetra::Vector<ST, LO, GO, Node>       vector_type;

typedef Tpetra::Map<LO, GO, Node>         map_type;
typedef Tpetra::Export<LO, GO, Node>      export_type;
typedef Tpetra::Import<LO, GO, Node>      import_type;
typedef Tpetra::CrsGraph<LO, GO, Node>    sparse_graph_type;
typedef Sacado::CacheFad::SFad<double,8> FadType;


}

#endif /* TYPES_H_ */
