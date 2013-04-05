/*
 * Integration.cpp
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

#include "Integration.h"
#include "GraphCreation.h"
#include "Cubature.h"

#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_ArrayTools.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_RealSpaceTools.hpp>
#include <Intrepid_DefaultCubatureFactory.hpp>
#include <Intrepid_Utils.hpp>
#include <Intrepid_FieldContainer.hpp>

// Sacado includes
#include "Sacado.hpp"

namespace Moirai {

template<typename Scalar>
void materialTensor (Scalar material[][3], const Scalar& x, const Scalar& y,const Scalar& z) {
	typedef Teuchos::ScalarTraits<Scalar> STS;

	material[0][0] = STS::one();
	material[0][1] = STS::zero();
	material[0][2] = STS::zero();

	material[1][0] = STS::zero();
	material[1][1] = STS::one();
	material[1][2] = STS::zero();

	material[2][0] = STS::zero();
	material[2][1] = STS::zero();
	material[2][2] = STS::one();

}


//! Compute the material tensor over a workset.
template<class ArrayOut, class ArrayIn>
void evaluateMaterialTensor (ArrayOut& matTensorValues, const ArrayIn& evaluationPoints) {
  typedef typename ArrayOut::scalar_type scalar_type;

  const int numWorksetCells  = evaluationPoints.dimension(0);
  const int numPoints        = evaluationPoints.dimension(1);
  const int spaceDim         = evaluationPoints.dimension(2);

  scalar_type material[3][3];

  for (int cell = 0; cell < numWorksetCells; ++cell) {
    for (int pt = 0; pt < numPoints; ++pt) {
      scalar_type x = evaluationPoints(cell, pt, 0);
      scalar_type y = evaluationPoints(cell, pt, 1);
      scalar_type z = evaluationPoints(cell, pt, 2);

      materialTensor<scalar_type> (material, x, y, z);

      for (int row = 0; row < spaceDim; ++row) {
        for(int col = 0; col < spaceDim; ++col) {
          matTensorValues(cell, pt, row, col) = material[row][col];
        }
      }
    }
  }
}

RCP<sparse_matrix_type> construct_stiffness_matrix(Mesh& mesh, std::function<FadType(FadType)> function) {

	using namespace Intrepid;

	typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
	typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools;
	typedef Intrepid::CellTools<ST>      IntrepidCTools;

	/**********************************************************************************/
	/********************************* GET CUBATURE For 3D cells***********************/
	/**********************************************************************************/
	shards::CellTopology cellType = get_cell_type(mesh);
	DefaultCubatureFactory<ST>  cubFactory;
	const double cubDegree = 2;
	RCP<Intrepid::Cubature<ST> > cubature = cubFactory.create (cellType, cubDegree);

	int cubDim       = cubature->getDimension ();
	int numCubPoints = cubature->getNumPoints ();

	Intrepid::FieldContainer<ST> cubPoints(numCubPoints, cubDim);
	Intrepid::FieldContainer<ST> cubWeights(numCubPoints);

	cubature->getCubature (cubPoints, cubWeights);


	/**********************************************************************************/
	/*********************************** GET BASIS ************************************/
	/**********************************************************************************/

	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis = get_basis(mesh);

	int numFieldsG = HGradBasis->getCardinality();
	Intrepid::FieldContainer<ST> HGBValues(numFieldsG, numCubPoints);
	Intrepid::FieldContainer<ST> HGBGrads(numFieldsG, numCubPoints, mesh.get_space_dim());

	// Evaluate basis values and gradients at cubature points
	HGradBasis->getValues(HGBValues, cubPoints, OPERATOR_VALUE);
	HGradBasis->getValues(HGBGrads, cubPoints, OPERATOR_GRAD);


	/**********************************************************************************/
	/*********************************** create and fill K*****************************/
	/**********************************************************************************/

	const int numElems = mesh.get_num_cells();
	const long long numNodes = mesh.get_num_interior_nodes();
	const int spaceDim = mesh.get_space_dim();

	RCP<sparse_graph_type> graph = create_sparse_graph_interior(mesh);

	RCP<sparse_matrix_type> K = rcp (new sparse_matrix_type (graph.getConst ()));


	K->setAllToScalar(0);


	/**********************************************************************************/
	/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
	/**********************************************************************************/


	// Define desired workset size and count how many worksets there are
	// on this processor's mesh block
	int desiredWorksetSize = numElems; // change to desired workset size!
	//int desiredWorksetSize = 100;    // change to desired workset size!
	int numWorksets        = numElems/desiredWorksetSize;

	// When numElems is not divisible by desiredWorksetSize, increase
	// workset count by 1
	if (numWorksets*desiredWorksetSize < numElems) {
		numWorksets += 1;
	}

	LOG(2,"Desired workset size:             " << desiredWorksetSize << std::endl
			<< "Number of worksets (per process): " << numWorksets);

	for (int workset = 0; workset < numWorksets; ++workset) {
		// Compute cell numbers where the workset starts and ends
		int worksetSize  = 0;
		int worksetBegin = (workset + 0)*desiredWorksetSize;
		int worksetEnd   = (workset + 1)*desiredWorksetSize;

		// When numElems is not divisible by desiredWorksetSize, the last
		// workset ends at numElems.
		worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

		// Now we know the actual workset size and can allocate the array
		// for the cell nodes.
		worksetSize = worksetEnd - worksetBegin;
		FieldContainer<ST> cellWorkset (worksetSize, numFieldsG, spaceDim);

		// array to contain boundary normals (=0 if not on boundary)
		FieldContainer<ST> boundary_normals(worksetSize, spaceDim);

		// Copy coordinates into cell workset
		int cellCounter = 0;
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
			for (int node = 0; node < numFieldsG; ++node) {
				const int node_num = mesh.get_cells()[cell].node_ids[node];
				for (int d = 0; d < spaceDim; ++d) {
					cellWorkset(cellCounter, node, d) = mesh.get_nodes()[node_num].x[d];
				}
			}
			++cellCounter;
		}

		/**********************************************************************************/
		/*                                Allocate arrays                                 */
		/**********************************************************************************/

		// Containers for Jacobians, integration measure & cubature points in workset cells
		FieldContainer<ST> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
		FieldContainer<ST> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
		FieldContainer<ST> worksetJacobDet  (worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubWeights(worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubPoints (worksetSize, numCubPoints, cubDim);

		// Containers for basis values transformed to workset cells and
		// them multiplied by cubature weights
		FieldContainer<ST> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
		FieldContainer<ST> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);
		FieldContainer<ST> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
		FieldContainer<ST> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

		FieldContainer<ST> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

		// Containers for material values and source term. Require
		// user-defined functions
		FieldContainer<ST> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);

		// Containers for workset contributions to the discretization
		// matrix and the right hand side
		FieldContainer<ST> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
		FieldContainer<ST> worksetGradOp (worksetSize, numFieldsG, numFieldsG);

		// Additional arrays used in AD-based assembly.
		FieldContainer<FadType> u_coeffsAD(worksetSize, numFieldsG);
		FieldContainer<FadType> u_FE_gradAD(worksetSize, numCubPoints, spaceDim);
		FieldContainer<FadType> u_FE_valAD(worksetSize, numCubPoints);
		FieldContainer<FadType> f_of_u_AD(worksetSize, numCubPoints);
		FieldContainer<FadType> cellResidualAD(worksetSize, numFieldsG);
		for (int c=0; c<worksetSize; c++) {
			for(int f=0; f<numFieldsG; f++) {
				u_coeffsAD(c,f) = FadType(numFieldsG, f, 1.0);
			}
		}

		/**********************************************************************************/
		/*                                Calculate Jacobians                             */
		/**********************************************************************************/

		IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
		IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
		IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

		/**********************************************************************************/
		/*          Cubature Points to Physical Frame and Compute Data                    */
		/**********************************************************************************/

		// Map cubature points to physical frame.
		IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

		// Evaluate the material tensor A at cubature points.
		evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

		/**********************************************************************************/
		/*                         Compute Stiffness Matrix                               */
		/**********************************************************************************/

		// Transform basis gradients to physical frame:
		IntrepidFSTools::HGRADtransformGRAD<ST> (worksetHGBGrads,   // DF^{-T}(grad u)
				worksetJacobInv,
				HGBGrads);
		// Compute integration measure for workset cells:
		IntrepidFSTools::computeCellMeasure<ST> (worksetCubWeights, // Det(DF)*w = J*w
				worksetJacobDet,
				cubWeights);
		// Multiply transformed (workset) gradients with weighted measure
		IntrepidFSTools::multiplyMeasure<ST> (worksetHGBGradsWeighted, // DF^{-T}(grad u)*J*w
				worksetCubWeights,
				worksetHGBGrads);
		// Compute the diffusive flux:
		IntrepidFSTools::tensorMultiplyDataField<ST> (worksetDiffusiveFlux, // A*(DF^{-T}(grad u)
				worksetMaterialVals,
				worksetHGBGrads);
		// Integrate to compute workset diffusion contribution to global matrix:
		IntrepidFSTools::integrate<ST> (worksetStiffMatrix, // (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
				worksetHGBGradsWeighted,
				worksetDiffusiveFlux,
				COMP_BLAS);

		/**********************************************************************************/
		/*                                   Compute Reaction                             */
		/**********************************************************************************/
		// Transform basis values to physical frame:
		IntrepidFSTools::HGRADtransformVALUE<ST> (worksetHGBValues, // clones basis values (u)
				HGBValues);
		// Multiply transformed (workset) values with weighted measure
		IntrepidFSTools::multiplyMeasure<ST> (worksetHGBValuesWeighted, // (u)*J*w
				worksetCubWeights,
				worksetHGBValues);

		// represent gradient of the current state (iterate) as a
		// linear combination of the gradients of basis functions
		// use AD arrays !
		u_FE_gradAD.initialize();
		IntrepidFSTools::evaluate<FadType>(u_FE_gradAD,
				u_coeffsAD,
				worksetHGBGrads);

		// represent value of the current state (iterate) as a
		// linear combination of the basis functions
		// use AD arrays !
		u_FE_valAD.initialize();
		IntrepidFSTools::evaluate<FadType>(u_FE_valAD,
				u_coeffsAD,
				worksetHGBValues);

		// compute nonlinear term
		for(int c=0; c<worksetSize; c++){
			for(int p=0; p<numCubPoints; p++){
				f_of_u_AD(c,p) = function(u_FE_valAD(c,p));
				//f_of_u_AD(c,p) = u_FE_valAD(c,p)*(1000.0-u_FE_valAD(c,p));
			}
		}

		// integrate to compute element residual
		IntrepidFSTools::integrate<FadType>(cellResidualAD,
				u_FE_gradAD,
				worksetHGBGradsWeighted, COMP_BLAS);
		IntrepidFSTools::integrate<FadType>(cellResidualAD,
				f_of_u_AD,
				worksetHGBValuesWeighted, COMP_BLAS, true);


		/**********************************************************************************/
		/*                         Assemble into Global Matrix                            */
		/**********************************************************************************/


		// "WORKSET CELL" loop: local cell ordinal is relative to numElems
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

			// Compute cell ordinal relative to the current workset
			const int worksetCellOrdinal = cell - worksetBegin;

			// "CELL EQUATION" loop for the workset cell: cellRow is
			// relative to the cell DoF numbering.
			for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
				LO localRow  = mesh.get_cells()[cell].node_ids[cellRow];
				std::vector<LO> local_columns;
				std::vector<ST> local_columns_contribution;
				for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
					local_columns.push_back(mesh.get_cells()[cell].node_ids[cellCol]);
					std::cout << 
					        *cellResidualAD (worksetCellOrdinal, cellRow).dx()	
						  << std::endl;
					local_columns_contribution.push_back(
					        *cellResidualAD (worksetCellOrdinal, cellRow).dx()	
					);
				}// *** cell col loop ***
				K->sumIntoGlobalValues (localRow, Teuchos::arrayViewFromVector(local_columns),
						Teuchos::arrayViewFromVector(local_columns_contribution));

			}// *** cell row loop ***
		}// *** workset cell loop **
	}// *** workset loop ***

	K->fillComplete();

	return K;
}

RCP<sparse_matrix_type> construct_stiffness_matrix(Mesh& mesh) {
	using namespace Intrepid;

	typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
	typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools;
	typedef Intrepid::CellTools<ST>      IntrepidCTools;

	/**********************************************************************************/
	/********************************* GET CUBATURE For 3D cells***********************/
	/**********************************************************************************/
	shards::CellTopology cellType = get_cell_type(mesh);
	DefaultCubatureFactory<ST>  cubFactory;
	const double cubDegree = 2;
	RCP<Intrepid::Cubature<ST> > cubature = cubFactory.create (cellType, cubDegree);

	int cubDim       = cubature->getDimension ();
	int numCubPoints = cubature->getNumPoints ();

	Intrepid::FieldContainer<ST> cubPoints(numCubPoints, cubDim);
	Intrepid::FieldContainer<ST> cubWeights(numCubPoints);

	cubature->getCubature (cubPoints, cubWeights);


	/**********************************************************************************/
	/*********************************** GET BASIS ************************************/
	/**********************************************************************************/

	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis = get_basis(mesh);

	int numFieldsG = HGradBasis->getCardinality();
	//Intrepid::FieldContainer<ST> HGBValues(numFieldsG, numCubPoints);
	Intrepid::FieldContainer<ST> HGBGrads(numFieldsG, numCubPoints, mesh.get_space_dim());

	// Evaluate basis values and gradients at cubature points
	//HGradBasis->getValues(HGBValues, cubPoints, OPERATOR_VALUE);
	HGradBasis->getValues(HGBGrads, cubPoints, OPERATOR_GRAD);


	/**********************************************************************************/
	/*********************************** create and fill K*****************************/
	/**********************************************************************************/

	const int numElems = mesh.get_num_cells();
	const long long numNodes = mesh.get_num_interior_nodes();
	const int spaceDim = mesh.get_space_dim();

	RCP<sparse_graph_type> graph = create_sparse_graph_interior(mesh);

	RCP<sparse_matrix_type> K = rcp (new sparse_matrix_type (graph.getConst ()));


	K->setAllToScalar(0);


	/**********************************************************************************/
	/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
	/**********************************************************************************/


	// Define desired workset size and count how many worksets there are
	// on this processor's mesh block
	int desiredWorksetSize = numElems; // change to desired workset size!
	//int desiredWorksetSize = 100;    // change to desired workset size!
	int numWorksets        = numElems/desiredWorksetSize;

	// When numElems is not divisible by desiredWorksetSize, increase
	// workset count by 1
	if (numWorksets*desiredWorksetSize < numElems) {
		numWorksets += 1;
	}

	LOG(2,"Desired workset size:             " << desiredWorksetSize << std::endl
			<< "Number of worksets (per process): " << numWorksets);

	for (int workset = 0; workset < numWorksets; ++workset) {
		// Compute cell numbers where the workset starts and ends
		int worksetSize  = 0;
		int worksetBegin = (workset + 0)*desiredWorksetSize;
		int worksetEnd   = (workset + 1)*desiredWorksetSize;

		// When numElems is not divisible by desiredWorksetSize, the last
		// workset ends at numElems.
		worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

		// Now we know the actual workset size and can allocate the array
		// for the cell nodes.
		worksetSize = worksetEnd - worksetBegin;
		FieldContainer<ST> cellWorkset (worksetSize, numFieldsG, spaceDim);

		// array to contain boundary normals (=0 if not on boundary)
		FieldContainer<ST> boundary_normals(worksetSize, spaceDim);

		// Copy coordinates into cell workset
		int cellCounter = 0;
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
			for (int node = 0; node < numFieldsG; ++node) {
				const int node_num = mesh.get_cells()[cell].node_ids[node];
				for (int d = 0; d < spaceDim; ++d) {
					cellWorkset(cellCounter, node, d) = mesh.get_nodes()[node_num].x[d];
				}
			}
			++cellCounter;
		}

		/**********************************************************************************/
		/*                                Allocate arrays                                 */
		/**********************************************************************************/

		// Containers for Jacobians, integration measure & cubature points in workset cells
		FieldContainer<ST> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
		FieldContainer<ST> worksetJacobInv  (worksetSize, numCubPoints, spaceDim, spaceDim);
		FieldContainer<ST> worksetJacobDet  (worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubWeights(worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubPoints (worksetSize, numCubPoints, cubDim);

		// Containers for basis values transformed to workset cells and
		// them multiplied by cubature weights
		FieldContainer<ST> worksetHGBGrads         (worksetSize, numFieldsG, numCubPoints, spaceDim);
		FieldContainer<ST> worksetHGBGradsWeighted (worksetSize, numFieldsG, numCubPoints, spaceDim);

		FieldContainer<ST> worksetDiffusiveFlux(worksetSize, numFieldsG, numCubPoints, spaceDim);

		// Containers for material values and source term. Require
		// user-defined functions
		FieldContainer<ST> worksetMaterialVals (worksetSize, numCubPoints, spaceDim, spaceDim);

		// Containers for workset contributions to the discretization
		// matrix and the right hand side
		FieldContainer<ST> worksetStiffMatrix (worksetSize, numFieldsG, numFieldsG);
		FieldContainer<ST> worksetGradOp (worksetSize, numFieldsG, numFieldsG);


		/**********************************************************************************/
		/*                                Calculate Jacobians                             */
		/**********************************************************************************/

		IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
		IntrepidCTools::setJacobianInv(worksetJacobInv, worksetJacobian );
		IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

		/**********************************************************************************/
		/*          Cubature Points to Physical Frame and Compute Data                    */
		/**********************************************************************************/

		// Map cubature points to physical frame.
		IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);

		// Evaluate the material tensor A at cubature points.
		evaluateMaterialTensor (worksetMaterialVals, worksetCubPoints);

		/**********************************************************************************/
		/*                         Compute Stiffness Matrix                               */
		/**********************************************************************************/

		// Transform basis gradients to physical frame:
		IntrepidFSTools::HGRADtransformGRAD<ST> (worksetHGBGrads,   // DF^{-T}(grad u)
				worksetJacobInv,
				HGBGrads);
		// Compute integration measure for workset cells:
		IntrepidFSTools::computeCellMeasure<ST> (worksetCubWeights, // Det(DF)*w = J*w
				worksetJacobDet,
				cubWeights);
		// Multiply transformed (workset) gradients with weighted measure
		IntrepidFSTools::multiplyMeasure<ST> (worksetHGBGradsWeighted, // DF^{-T}(grad u)*J*w
				worksetCubWeights,
				worksetHGBGrads);
		// Compute the diffusive flux:
		IntrepidFSTools::tensorMultiplyDataField<ST> (worksetDiffusiveFlux, // A*(DF^{-T}(grad u)
				worksetMaterialVals,
				worksetHGBGrads);
		// Integrate to compute workset diffusion contribution to global matrix:
		IntrepidFSTools::integrate<ST> (worksetStiffMatrix, // (DF^{-T}(grad u)*J*w)*(A*DF^{-T}(grad u))
				worksetHGBGradsWeighted,
				worksetDiffusiveFlux,
				COMP_BLAS);


		/**********************************************************************************/
		/*                         Assemble into Global Matrix                            */
		/**********************************************************************************/


		// "WORKSET CELL" loop: local cell ordinal is relative to numElems
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

			// Compute cell ordinal relative to the current workset
			const int worksetCellOrdinal = cell - worksetBegin;

			// "CELL EQUATION" loop for the workset cell: cellRow is
			// relative to the cell DoF numbering.
			for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
				LO localRow  = mesh.get_cells()[cell].node_ids[cellRow];
				std::vector<LO> local_columns;
				std::vector<ST> local_columns_contribution;
				for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
					local_columns.push_back(mesh.get_cells()[cell].node_ids[cellCol]);
					local_columns_contribution.push_back(worksetStiffMatrix (worksetCellOrdinal, cellRow, cellCol));
				}// *** cell col loop ***
				K->sumIntoGlobalValues (localRow, Teuchos::arrayViewFromVector(local_columns),
												  Teuchos::arrayViewFromVector(local_columns_contribution));

			}// *** cell row loop ***
		}// *** workset cell loop **
	}// *** workset loop ***

	K->fillComplete();

	return K;
}

RCP<sparse_matrix_type> construct_mass_matrix(Mesh& mesh) {
	using namespace Intrepid;

	typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
	typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools;
	typedef Intrepid::CellTools<ST>      IntrepidCTools;

	/**********************************************************************************/
	/********************************* GET CUBATURE For 3D cells***********************/
	/**********************************************************************************/
	shards::CellTopology cellType = get_cell_type(mesh);
	DefaultCubatureFactory<ST>  cubFactory;
	const double cubDegree = 2;
	RCP<Intrepid::Cubature<ST> > cubature = cubFactory.create (cellType, cubDegree);

	int cubDim       = cubature->getDimension ();
	int numCubPoints = cubature->getNumPoints ();

	Intrepid::FieldContainer<ST> cubPoints(numCubPoints, cubDim);
	Intrepid::FieldContainer<ST> cubWeights(numCubPoints);

	cubature->getCubature (cubPoints, cubWeights);


	/**********************************************************************************/
	/*********************************** GET BASIS ************************************/
	/**********************************************************************************/

	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis = get_basis(mesh);

	int numFieldsG = HGradBasis->getCardinality();
	Intrepid::FieldContainer<ST> HGBValues(numFieldsG, numCubPoints);
	//Intrepid::FieldContainer<ST> HGBGrads(numFieldsG, numCubPoints, mesh.get_space_dim());

	// Evaluate basis values and gradients at cubature points
	HGradBasis->getValues(HGBValues, cubPoints, OPERATOR_VALUE);
	//HGradBasis->getValues(HGBGrads, cubPoints, OPERATOR_GRAD);


	/**********************************************************************************/
	/*********************************** create and fill M*****************************/
	/**********************************************************************************/

	const int numElems = mesh.get_num_cells();
	const long long numNodes = mesh.get_num_interior_nodes();
	const int spaceDim = mesh.get_space_dim();

	RCP<sparse_graph_type> graph = create_sparse_graph_interior(mesh);

	RCP<sparse_matrix_type> M = rcp (new sparse_matrix_type (graph.getConst ()));


	M->setAllToScalar(0);

	/**********************************************************************************/
	/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
	/**********************************************************************************/


	// Define desired workset size and count how many worksets there are
	// on this processor's mesh block
	int desiredWorksetSize = numElems; // change to desired workset size!
	//int desiredWorksetSize = 100;    // change to desired workset size!
	int numWorksets        = numElems/desiredWorksetSize;

	// When numElems is not divisible by desiredWorksetSize, increase
	// workset count by 1
	if (numWorksets*desiredWorksetSize < numElems) {
		numWorksets += 1;
	}

	LOG(2,"Desired workset size:             " << desiredWorksetSize << std::endl
			<< "Number of worksets (per process): " << numWorksets);

	for (int workset = 0; workset < numWorksets; ++workset) {
		// Compute cell numbers where the workset starts and ends
		int worksetSize  = 0;
		int worksetBegin = (workset + 0)*desiredWorksetSize;
		int worksetEnd   = (workset + 1)*desiredWorksetSize;

		// When numElems is not divisible by desiredWorksetSize, the last
		// workset ends at numElems.
		worksetEnd = (worksetEnd <= numElems) ? worksetEnd : numElems;

		// Now we know the actual workset size and can allocate the array
		// for the cell nodes.
		worksetSize = worksetEnd - worksetBegin;
		FieldContainer<ST> cellWorkset (worksetSize, numFieldsG, spaceDim);

		// array to contain boundary normals (=0 if not on boundary)
		FieldContainer<ST> boundary_normals(worksetSize, spaceDim);

		// Copy coordinates into cell workset
		int cellCounter = 0;
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {
			for (int node = 0; node < numFieldsG; ++node) {
				const int node_num = mesh.get_cells()[cell].node_ids[node];
				for (int d = 0; d < spaceDim; ++d) {
					cellWorkset(cellCounter, node, d) = mesh.get_nodes()[node_num].x[d];
				}
			}
			++cellCounter;
		}

		/**********************************************************************************/
		/*                                Allocate arrays                                 */
		/**********************************************************************************/

		// Containers for Jacobians, integration measure & cubature points in workset cells
		FieldContainer<ST> worksetJacobian  (worksetSize, numCubPoints, spaceDim, spaceDim);
		FieldContainer<ST> worksetJacobDet  (worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubWeights(worksetSize, numCubPoints);
		FieldContainer<ST> worksetCubPoints (worksetSize, numCubPoints, cubDim);

		// Containers for basis values transformed to workset cells and
		// them multiplied by cubature weights
		FieldContainer<ST> worksetHGBValues        (worksetSize, numFieldsG, numCubPoints);
		FieldContainer<ST> worksetHGBValuesWeighted(worksetSize, numFieldsG, numCubPoints);


		// Containers for workset contributions to the discretization
		// matrix and the right hand side
		FieldContainer<ST> worksetMassMatrix (worksetSize, numFieldsG, numFieldsG);


		/**********************************************************************************/
		/*                                Calculate Jacobians                             */
		/**********************************************************************************/

		IntrepidCTools::setJacobian(worksetJacobian, cubPoints, cellWorkset, cellType);
		IntrepidCTools::setJacobianDet(worksetJacobDet, worksetJacobian );

		/**********************************************************************************/
		/*          Cubature Points to Physical Frame and Compute Data                    */
		/**********************************************************************************/

		// Map cubature points to physical frame.
		IntrepidCTools::mapToPhysicalFrame (worksetCubPoints, cubPoints, cellWorkset, cellType);


		/**********************************************************************************/
		/*                         Compute Mass Matrix                               */
		/**********************************************************************************/


		//Transform basis values to physical frame:
		IntrepidFSTools::HGRADtransformVALUE<ST> (worksetHGBValues, // clones basis values (u)
				HGBValues);
		// Compute integration measure for workset cells:
		IntrepidFSTools::computeCellMeasure<ST> (worksetCubWeights, // Det(DF)*w = J*w
				worksetJacobDet,
				cubWeights);

		// Multiply transformed (workset) gradients with weighted measure
		IntrepidFSTools::multiplyMeasure<ST> (worksetHGBValuesWeighted, // (u)*w
				worksetCubWeights,
				worksetHGBValues);
		// Integrate to compute workset contribution to global matrix:
		IntrepidFSTools::integrate<ST> (worksetMassMatrix, // (u)*(u)*w
				worksetHGBValues,
				worksetHGBValuesWeighted,
				COMP_BLAS);


		/**********************************************************************************/
		/*                         Assemble into Global Matrix                            */
		/**********************************************************************************/


		// "WORKSET CELL" loop: local cell ordinal is relative to numElems
		for (int cell = worksetBegin; cell < worksetEnd; ++cell) {

			// Compute cell ordinal relative to the current workset
			const int worksetCellOrdinal = cell - worksetBegin;

			// "CELL EQUATION" loop for the workset cell: cellRow is
			// relative to the cell DoF numbering.
			for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
				LO localRow  = mesh.get_cells()[cell].node_ids[cellRow];
				std::vector<LO> local_columns;
				std::vector<ST> local_columns_contribution;
				for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
					local_columns.push_back(mesh.get_cells()[cell].node_ids[cellCol]);
					local_columns_contribution.push_back(worksetMassMatrix (worksetCellOrdinal, cellRow, cellCol));
				}// *** cell col loop ***
				M->sumIntoGlobalValues (localRow, Teuchos::arrayViewFromVector(local_columns),
						Teuchos::arrayViewFromVector(local_columns_contribution));

			}// *** cell row loop ***
		}// *** workset cell loop **
	}// *** workset loop ***

	M->fillComplete();


	return M;
}

RCP<sparse_matrix_type> construct_boundary_mass_matrix(Mesh& mesh) {
	using namespace Intrepid;

	typedef Intrepid::FunctionSpaceTools IntrepidFSTools;
	typedef Intrepid::RealSpaceTools<ST> IntrepidRSTools;
	typedef Intrepid::CellTools<ST>      IntrepidCTools;

	/**********************************************************************************/
	/********************************* GET CUBATURE For 2D cells***********************/
	/**********************************************************************************/
	shards::CellTopology faceType = get_face_type(mesh);
	shards::CellTopology cellType = get_cell_type(mesh);
	DefaultCubatureFactory<ST>  cubFactory;
	const double cubDegree = 2;
	RCP<Intrepid::Cubature<ST> > cubature = cubFactory.create (faceType, cubDegree);

	int cubDim       = cubature->getDimension ();
	int numCubPoints = cubature->getNumPoints ();

	Intrepid::FieldContainer<ST> cubPoints(numCubPoints, cubDim);
	Intrepid::FieldContainer<ST> cubWeights(numCubPoints);

	cubature->getCubature (cubPoints, cubWeights);


	/**********************************************************************************/
	/*********************************** GET BASIS ************************************/
	/**********************************************************************************/

	// Define basis
	// select basis based on cell topology only for now, and assume first order basis
	RCP<Intrepid::Basis<ST, Intrepid::FieldContainer<ST> > > HGradBasis = get_face_basis(mesh);

	int numFieldsG = HGradBasis->getCardinality();
	int numNodesPerElem = mesh.get_nodes_per_cell();
	Intrepid::FieldContainer<ST> HGBValues(numFieldsG, numCubPoints);

	// Evaluate basis values and gradients at cubature points
	HGradBasis->getValues(HGBValues, cubPoints, OPERATOR_VALUE);


	/**********************************************************************************/
	/*********************************** create and fill M*****************************/
	/**********************************************************************************/

	const int numBoundaryFaces = mesh.get_num_boundary_faces();
	const int spaceDim = mesh.get_space_dim();


	RCP<sparse_graph_type> graph = create_sparse_graph_boundary_boundary(mesh);

	RCP<sparse_matrix_type> M = rcp (new sparse_matrix_type (graph.getConst ()));


	M->setAllToScalar(0);

	/**********************************************************************************/
	/******************** DEFINE WORKSETS AND LOOP OVER THEM **************************/
	/**********************************************************************************/
	if (numBoundaryFaces == 0) return M;
	// Define desired workset size and count how many worksets there are
	// on this processor's mesh block
	int desiredWorksetSize = numBoundaryFaces; // change to desired workset size!
	//int desiredWorksetSize = 100;    // change to desired workset size!
	int numWorksets        = numBoundaryFaces/desiredWorksetSize;

	// When numElems is not divisible by desiredWorksetSize, increase
	// workset count by 1
	if (numWorksets*desiredWorksetSize < numBoundaryFaces) {
		numWorksets += 1;
	}

	LOG(2,"Desired workset size:             " << desiredWorksetSize << std::endl
			<< "Number of worksets (per process): " << numWorksets);

	for (int workset = 0; workset < numWorksets; ++workset) {
		// Compute cell numbers where the workset starts and ends
		int worksetSize  = 0;
		int worksetBegin = (workset + 0)*desiredWorksetSize;
		int worksetEnd   = (workset + 1)*desiredWorksetSize;

		// When numElems is not divisible by desiredWorksetSize, the last
		// workset ends at numElems.
		worksetEnd = (worksetEnd <= numBoundaryFaces) ? worksetEnd : numBoundaryFaces;

		// Now we know the actual workset size and can allocate the array
		// for the cell nodes.
		worksetSize = worksetEnd - worksetBegin;
		FieldContainer<ST> cellWorkset (1, numNodesPerElem, spaceDim);
		FieldContainer<ST> worksetCubWeights(1, numCubPoints);

		FieldContainer<ST> worksetRefCubPoints (numCubPoints, spaceDim);
		// Copy coordinates and face cubature points (in the ref cell domain)
		// into cell workset
		FieldContainer<ST> worksetJacobian  (1, numCubPoints, spaceDim, spaceDim);

		FieldContainer<ST> worksetFaceValuesWeighted(1, numFieldsG, numCubPoints);
		FieldContainer<ST> worksetFaceValues        (1, numFieldsG, numCubPoints);

		// Containers for workset contributions to the boundary integral
		FieldContainer<ST> worksetWeakBC (1, numFieldsG, numFieldsG);

		for (int face = worksetBegin; face < worksetEnd; ++face) {
			const int iface = mesh.get_boundary_faces()[face].ordinal;
			const int ielem = mesh.get_boundary_faces()[face].cell_id;
			for (int node = 0; node < numNodesPerElem; ++node) {
				const int node_num = mesh.get_cells()[ielem].node_ids[node];
				for (int j = 0; j < spaceDim; ++j) {
					cellWorkset(0, node, j) = mesh.get_nodes()[node_num].x[j];
				}
			}


			IntrepidCTools::mapToReferenceSubcell(worksetRefCubPoints,
					cubPoints,
					2, iface, cellType);


			IntrepidCTools::setJacobian(worksetJacobian, worksetRefCubPoints,
					cellWorkset, cellType);

			IntrepidFSTools::computeFaceMeasure<ST> (worksetCubWeights, // Det(DF)*w = J*w
					worksetJacobian,
					cubWeights,
					iface,
					cellType);

			IntrepidFSTools::HGRADtransformVALUE<ST> (worksetFaceValues, // clones basis values (mu)
					HGBValues);

			IntrepidFSTools::multiplyMeasure<ST> (worksetFaceValuesWeighted, // (u)*w
					worksetCubWeights,
					worksetFaceValues);

			// Integrate to compute workset contribution to global matrix:
			IntrepidFSTools::integrate<ST> (worksetWeakBC, // (u)*(u)*w
					worksetFaceValues,
					worksetFaceValuesWeighted,
					COMP_BLAS);



			for (int cellRow = 0; cellRow < numFieldsG; ++cellRow) {
				LO localRow  = mesh.get_boundary_faces()[face].face_node_ids[cellRow];
				std::vector<LO> local_columns;
				std::vector<ST> local_columns_contribution;
				for (int cellCol = 0; cellCol < numFieldsG; ++cellCol) {
					local_columns.push_back(mesh.get_boundary_faces()[face].face_node_ids[cellCol]);
					local_columns_contribution.push_back(worksetWeakBC (0, cellRow, cellCol));
				}// *** cell col loop ***
				M->sumIntoGlobalValues (localRow, Teuchos::arrayViewFromVector(local_columns),
						Teuchos::arrayViewFromVector(local_columns_contribution));

			}// *** cell row loop ***
		}// *** workset cell loop **
	}// *** workset loop ***

	M->fillComplete();

	return M;
}

} /* namespace Moirai */
