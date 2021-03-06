/* 
 * Pde.h
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

#ifndef PDE_MOIRAI_H_
#define PDE_MOIRAI_H_

#include "Mesh.h"
#include "Types.h"

#include <vector>
#include <tuple>
#include <functional>

//#include <Thyra_TpetraThyraWrappers.hpp>
//
//#include <Thyra_LinearOpWithSolveBase.hpp>

#include <BelosSolverManager.hpp>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

namespace Moirai {

class Pde {
public:
	Pde(const double dx, const double dt);
	void timestep();
	void add_particle(const ST x, const ST y, const ST z);
	void generate_particles(std::vector<ST>& x,std::vector<ST>& y,std::vector<ST>& z);
	void react(const double ku, const double kb);
	vtkSmartPointer<vtkUnstructuredGrid> get_vtk_grid();
	std::string get_status_string();
	double get_number_of_particles();
private:
	double dt;
	int total_number_of_particles,number_of_particles_generated,number_of_particles_reacted;
	bool converged;
	int number_of_iterations;
	static constexpr double omega = 1.0;
	Mesh mesh;
	vtkSmartPointer<vtkUnstructuredGrid> vtk_grid;

	RCP<vector_type> volumes,areas;
	RCP<multivector_type> X,Y,u,f,lambda,flux,number_of_particles;
	RCP<sparse_matrix_type> LHS,LHS_prec,RHS,K,Mi,Mb;

//	RCP<Thyra::MultiVectorBase<ST> > X_w,Y_w;
//	RCP<const Thyra::LinearOpBase<ST> > RHS_w;
//	RCP<Thyra::LinearOpWithSolveBase<ST> > LHS_w;

	RCP<Belos::SolverManager<ST, multivector_type, operator_type> > solver;

	void initialise_from_mesh();
};

} /* namespace Moirai */
#endif /* PDE_H_ */
