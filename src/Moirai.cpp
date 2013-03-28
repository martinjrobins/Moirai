/*
 * Moirai.cpp
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
 *  Created on: 28 Mar 2013
 *      Author: robinsonm
 */

#include "Moirai.h"

namespace Moirai {
Teuchos::RCP<Teuchos::GlobalMPISession> mpiSession;

void init(int argc, char* argv[]) {
	mpiSession = Teuchos::rcp (new Teuchos::GlobalMPISession(&argc, &argv));
}

}
