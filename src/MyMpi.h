/* 
 * mpi.h
 *
 * Copyright 2012 Martin Robinson
 *
 * This file is part of PDE_BD.
 *
 * PDE_BD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PDE_BD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with PDE_BD.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Feb 16, 2013
 *      Author: mrobins
 */

#ifndef MYMPI_H_
#define MYMPI_H_

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

namespace Mpi {
	extern Teuchos::RCP<Teuchos::GlobalMPISession> mpiSession;
	void init(int argc, char *argv[]);
};

#endif /* MPI_H_ */
