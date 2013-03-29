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

#include <boost/random.hpp>

#include <Tpetra_RTI.hpp>
//#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosConfigDefs.hpp>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>


namespace Moirai {

constexpr double Pde::omega;

typedef boost::mt19937  base_generator_type;
static base_generator_type R;

Pde::Pde(const double dx, const double dt):
		dt(dt),
		total_number_of_particles(0),
		number_of_particles_generated(0) {

	mesh.initialise(dx);
	R.seed(time(NULL));
}

void Pde::timestep() {
	using Tpetra::RTI::reduce;
	using Tpetra::RTI::ZeroOp;
	using Tpetra::RTI::reductionGlob;

	if (RHS.is_null()) {
		initialise_from_mesh();
	}

	RHS->apply(*X.getConst(),*Y);


	solver->reset(Belos::Problem);
	Belos::ReturnType result = solver->solve ();

	converged = (result == Belos::Converged);

	Mi->apply(*u.getConst(),*number_of_particles);
	RCP<const vector_type> num_particles = number_of_particles->getVector(0);
	total_number_of_particles = TPETRA_REDUCE1(num_particles, num_particles, ZeroOp<ST>, std::plus<ST>());
}

void Pde::add_particle(const ST x, const ST y, const ST z) {
	const int local_index = mesh.get_nearest_node(x,y,z);
	X->getVectorNonConst(0)->sumIntoLocalValue(local_index, 1.0/volumes->get1dView()[local_index]);
	total_number_of_particles++;
}

void Pde::generate_particles(
		std::vector<ST>& x,std::vector<ST>& y,std::vector<ST>& z) {

	flux->elementWiseMultiply(1.0,*areas.getConst(),*lambda.getConst(),0.0);

	Teuchos::ArrayRCP<ST> flux_view = flux->get1dViewNonConst();
	//std::replace_if(bnv.begin(),bnv.end(),[](double d){return d<0;},0);
	std::vector<ST> flux_cumsum(flux_view.size());
	std::partial_sum(flux_view.begin(),flux_view.end(),flux_cumsum.begin());

	const double sum_values = *(flux_cumsum.end()-1);
	const double Ld = sum_values/dt;

	if (Ld <= 0) return;
	if (total_number_of_particles <= 0) return;

	boost::variate_generator<base_generator_type&, boost::uniform_real<> >
	U(R,boost::uniform_real<>(0,1));
	boost::variate_generator<base_generator_type&, boost::normal_distribution<> >
	N(R,boost::normal_distribution<>(0,1));
	boost::variate_generator<base_generator_type&, boost::poisson_distribution<> >
	P(R,boost::poisson_distribution<>(sum_values));

	//double tau = -log(U())/Ld;
	//number_of_particles_generated = 0;
	//while (tau < dt) {
	number_of_particles_generated = P();

	for (int i = 0; i < number_of_particles_generated; ++i) {
		const double r = U()*sum_values;
		std::vector<double>::iterator it = std::find_first_of(flux_cumsum.begin(),flux_cumsum.end(),&r,&r+1,
				[](double a, double b){return a >= b;});
		//		if ((it != bnv_cumsum.begin()) && (*it-r > r-*(it-1))) {
		//			it--;
		//		}
		const int r_find = mesh.get_local_boundary_node_ids()[it-flux_cumsum.begin()];
		//const int r_find = find_first(r,bnv_cumsum);
		//const double step_length = sqrt(2.0*D*(dt-tau));
		const double step_length = 0;

		x.push_back( mesh.get_nodes()[r_find].x[0] + step_length*N() );
		y.push_back( mesh.get_nodes()[r_find].x[1] + step_length*N() );
		z.push_back( mesh.get_nodes()[r_find].x[2] + step_length*N() );

		//tau = tau - log(U())/Ld;
		//number_of_particles_generated++;
	}

	u->scale(ST(total_number_of_particles-number_of_particles_generated+Ld*dt)/ST(total_number_of_particles));
	total_number_of_particles = total_number_of_particles-number_of_particles_generated+Ld*dt;
}


void Pde::set_reaction(std::function<ST(ST)> f) {
	function = f;
}

void Pde::initialise_from_mesh() {
	using Tpetra::RTI::reduce;
	using Tpetra::RTI::ZeroOp;
	using Tpetra::RTI::reductionGlob;

	/*
	 * Create vector and sparse matrix maps
	 */
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
	RCP<const map_type> interiorMap = Tpetra::createNonContigMapWithNode<LO,GO,Node>(
			Teuchos::arrayViewFromVector(node_ids),comm,node);
	RCP<const map_type> boundaryMap = Tpetra::createNonContigMapWithNode<LO,GO,Node>(
			Teuchos::arrayViewFromVector(bnode_ids),comm,node);
	RCP<const map_type> allMap = Tpetra::createNonContigMapWithNode<LO,GO,Node>(
			Teuchos::arrayViewFromVector(all_node_ids),comm,node);

	/*
	 * Create state vectors
	 */
	X = Tpetra::createMultiVector<ST,LO,GO,Node>(allMap,1);
	Y = Tpetra::createMultiVector<ST,LO,GO,Node>(allMap,1);
	u = X->offsetViewNonConst(interiorMap,0);
	lambda = X->offsetViewNonConst(boundaryMap,ni);
	flux = Tpetra::createMultiVector<ST,LO,GO,Node>(lambda->getMap(),1);
	number_of_particles = Tpetra::createMultiVector<ST,LO,GO,Node>(u->getMap(),1);
	volumes = Tpetra::createVector<ST,LO,GO,Node>(u->getMap());
	areas = Tpetra::createVector<ST,LO,GO,Node>(lambda->getMap());


	/*
	 * Create stiffness and mass matricies
	 */

	if (function) {
		K = construct_stiffness_matrix(mesh,function);
	} else {
		K = construct_stiffness_matrix(mesh);
	}
	Mi = construct_mass_matrix(mesh);
	Mb = construct_boundary_mass_matrix(mesh);

	/*
	 * calculate areas and volumes of nodes
	 */
	X->getVectorNonConst(0)->putScalar(1.0);
	Mi->apply(*u->getVector(0),*volumes);

	std::cout << "total volume is " <<
			TPETRA_REDUCE1(volumes, volumes, ZeroOp<ST>, std::plus<ST>()) << std::endl;

	Mb->apply(*lambda->getVector(0),*areas);

	std::cout << "total boundary area is " <<
			TPETRA_REDUCE1(areas, areas, ZeroOp<ST>, std::plus<ST>()) << std::endl;

	X->putScalar(0.0);

	RCP<sparse_matrix_type> A = Mi->convert<ST>();
	RCP<sparse_matrix_type> A_rhs = Mi->convert<ST>();

	/*
	 * A = Mi + omega*dt*K
	 */
	update(1.0,A,omega*dt,K);

	RCP<sparse_matrix_type> B = scatter_columns(Mb,interiorMap,mesh.get_global_interior_ids_on_boundary());

	RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(rcp(&std::cout,false));
//		Mb->describe(*out,Teuchos::VERB_EXTREME);
//		B->describe(*out,Teuchos::VERB_EXTREME);

	/*
	 * B = omega*B
	 */
	B->resumeFill();
	B->scale(omega);
	B->fillComplete();

//	B->describe(*out,Teuchos::VERB_EXTREME);


	/*
	 * LHS = [A,omegaB;omegaB^t,zero]
	 */
	LHS = construct_saddle_point_matrix(A,B);

//	A->describe(*out,Teuchos::VERB_EXTREME);
//	LHS->describe(*out,Teuchos::VERB_EXTREME);


	RCP<sparse_matrix_type> A11 = get11(A);
	RCP<sparse_matrix_type> Adiag = getDiag(A);

	/*
	 * LHS_prec = [A11,omegaB;omegaB^t,zero]
	 */
	LHS_prec = construct_saddle_point_matrix(Adiag,B);
//	const int n = LHS_prec->getNodeNumRows();
//	for (int row = 0; row < n; ++row) {
//		Teuchos::ArrayView<const LO> indicies;
//		Teuchos::ArrayView<const ST> values;
//		LHS_prec->getLocalRowView(row,indicies,values);
//		const int nc = indicies.size();
//		for (int j = 0; j < nc; ++j) {
//			if (values[j]==0) {
//				ERROR("found a zero!");
//			}
//		}
//	}

	/*
	 * A_rhs = Mi + (omega-1)*dt*K
	 */
	update(1.0,A_rhs,(omega-1)*dt,K);

	/*
	 * B = (omega-1)*B
	 */
	B->resumeFill();
	B->scale((omega-1.0)/omega);
	B->fillComplete();

	/*
	 * RHS = [A_rhs,(omega-1)B;(1-omega)B^t,zero]
	 */
	RHS = construct_anti_symetric_saddle_point_matrix(A_rhs,B);




	/*
	 * construct solver
	 */
	RCP<Teuchos::ParameterList> belosParams = Teuchos::parameterList ();
//	belosParams->set ("Block Size", 1);
//	belosParams->set ("Maximum Iterations", 1000);
//	belosParams->set ("Convergence Tolerance", sqrt(std::numeric_limits<ST>::min()));

	typedef Belos::LinearProblem<ST, multivector_type, operator_type > problem_type;
	RCP<problem_type> problem = rcp (new problem_type (LHS, X, Y));
//	problem->setLeftPrec (LHS_prec);
//	problem->setRightPrec (LHS_prec);

	const bool set = problem->setProblem ();
	TEUCHOS_TEST_FOR_EXCEPTION(! set, std::runtime_error, "solveWithBelos: The "
			"Belos::LinearProblem's setProblem() method returned false.  This probably "
			"indicates that there is something wrong with A, X, or B.");

	// Create the GMRES solver.
	Belos::SolverFactory<ST, multivector_type, operator_type> factory;
	solver = factory.create ("GMRES", belosParams);
//	Teuchos::Array<std::string> names = factory.supportedSolverNames();
//	std::cout << "list of supported solver names:" << std::endl;
//	for (int i = 0; i < names.size(); ++i) {
//		std::cout << names[i] << std::endl;
//	}
//	Teuchos::RCP<const Teuchos::ParameterList> params = solver->getCurrentParameters();
//	std::cout << "solver created with following parameters:" << std::endl;
//	for (Teuchos::ParameterList::ConstIterator i = params->begin(); i != params->end(); ++i) {
//		std::cout << i->first << ":  " << i->second << std::endl;
//	}
	// Tell the solver what problem you want to solve.
	solver->setProblem (problem);


//	/*
//	 * Create thyra objects to solve
//	 */
//	Stratimikos::DefaultLinearSolverBuilder strategy;
//	strategy.paramsXmlFileName("params.xml");
//	strategy.readParameters(&std::cout);
//	RCP<const Thyra::LinearOpWithSolveFactoryBase<ST> > lowsFactory = strategy.createLinearSolveStrategy("Belos");
//
//	X_w = Thyra::createMultiVector(X);
//	Y_w = Thyra::createMultiVector(Y);
//
//	RHS_w = Thyra::createConstLinearOp<ST, LO, GO, Node>(RHS);
//	RCP<const Thyra::LinearOpBase<ST> > LHS_prec_w = Thyra::createConstLinearOp<ST, LO, GO, Node>(LHS_prec);
//
//	RCP<const Thyra::LinearOpBase<ST> > LHS_op = Thyra::createConstLinearOp<ST, LO, GO, Node>(LHS);
//
//	LHS_w = lowsFactory->createOp();
//	//Thyra::initializeOp(*lowsFactory,LHSop,LHS.ptr());
//	//Thyra::initializeApproxPreconditionedOp(*lowsFactory,LHSop,LHS_prec,LHS.ptr());
//	//Thyra::initializePreconditionedOp<ST>(*lowsFactory,LHSop,Thyra::unspecifiedPrec<ST>(LHS_prec),LHS.ptr());
//	Thyra::initializePreconditionedOp<ST>(*lowsFactory,LHS_op,Thyra::rightPrec<ST>(LHS_prec_w),LHS_w.ptr());
}

vtkSmartPointer<vtkUnstructuredGrid> Pde::get_vtk_grid() {
	if (vtk_grid==NULL) {
		vtk_grid = mesh.get_vtk_grid();

		/*
		 * setup scalar data
		 */
		vtkSmartPointer<vtkDoubleArray> newScalars = vtkSmartPointer<vtkDoubleArray>::New();
		const int num_local_entries = u->getLocalLength();
		newScalars->SetArray(u->getDataNonConst(0).getRawPtr(),num_local_entries,1);
		newScalars->SetName("Concentration");

		vtk_grid->GetPointData()->SetScalars(newScalars);
	}
	return vtk_grid;
}

std::string Pde::get_status_string() {
	std::ostringstream ss;
	ss << "Pde Status:" << std::endl;
	ss << "\t" << total_number_of_particles << " particles." << std::endl;
	ss << "\t" << number_of_particles_generated << " particles generated." << std::endl;
	if (converged) {
		ss << "\tsolver CONVERGED after "<< solver->getNumIters() << " iterations " << std::endl;
	} else {
		ss << "\tsolver DID NOT CONVERGE after "<< solver->getNumIters() << " iterations " << std::endl;
	}
	return ss.str();
}

double Pde::get_number_of_particles() {
	return total_number_of_particles;
}

} /* namespace Moirai */
