/* 
 * test_pde_constructor.cpp
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

#include "Moirai.h"
#include "boost/format.hpp"
#include "zip.h"

#include <iostream>

int main(int argc, char **argv) {

	const double dt = 0.001;
	const double max_t = 20.0;
	const double dt_out = dt;
	const double dx = 0.1;

	Moirai::init(argc,argv);
	Moirai::Pde p(dx,dt);

	for (int i = 0; i < 1000; ++i) {
		p.add_particle(0,0,0);
	}

	std::vector<double> x,y,z;
	for (int i = 0; i < max_t/dt_out; ++i) {
		std::stringstream filename_grid,filename_boundary,filename_molecules;
		filename_grid <<boost::format("test%05d.pvtu")%i;
		filename_boundary<<boost::format("testBoundary%05d.pvtu")%i;
		filename_molecules<<boost::format("testMolecules%05d.pvtu")%i;


		const int iterations = int(dt_out/dt + 0.5);
		const double actual_dt = iterations*dt;
		for (int j = 0; j < iterations; ++j) {
			p.timestep();
			p.generate_particles(x,y,z);
			std::cout << "generated "<<x.size()<<" particles so far..."<<std::endl;
		}
	}
}


