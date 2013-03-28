/* 
 * io.h
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

#ifndef IO_MOIRAI_H_
#define IO_MOIRAI_H_

#include <string>
#include <vtkUnstructuredGrid.h>

namespace Moirai {
	void write_grid(std::string filename, vtkUnstructuredGrid* grid);
	void write_points(std::string filename,
			const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z);
	void write_column_vectors(std::string filename, std::string header,std::vector<std::vector<double>* >& columns);
}



#endif /* IO_H_ */
