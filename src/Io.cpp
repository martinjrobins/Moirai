/* 
 * io.cpp
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

#include "Io.h"
#include "Moirai.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "vtkSmartPointer.h"

#include <iostream>
#include <fstream>

namespace Moirai {

void write_grid(std::string filename, vtkUnstructuredGrid* grid) {
	const int my_rank = Moirai::mpiSession->getRank();
	const int num_procs = Moirai::mpiSession->getNProc();
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
	writer->SetNumberOfPieces(num_procs);
	writer->SetStartPiece(my_rank);
	writer->SetEndPiece(my_rank);
	writer->SetInput(grid);
	writer->SetDataModeToBinary();
	writer->SetFileName(filename.c_str());
	writer->Write();
}

void write_points(std::string filename, const std::vector<double>& x,
		const std::vector<double>& y, const std::vector<double>& z) {
	/*
	 * setup points
	 */
	vtkSmartPointer<vtkPoints> newPts = vtkSmartPointer<vtkPoints>::New();
	const int n = x.size();
	for (int i = 0; i < n; i++) {
		newPts->InsertNextPoint(x[i],y[i],z[i]);
	}

	vtkSmartPointer<vtkUnstructuredGrid> vtk_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtk_grid->SetPoints(newPts);

	const int my_rank = Moirai::mpiSession->getRank();
	const int num_procs = Moirai::mpiSession->getNProc();
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
	writer->SetNumberOfPieces(num_procs);
	writer->SetStartPiece(my_rank);
	writer->SetEndPiece(my_rank);
	writer->SetInput(vtk_grid);
	writer->SetDataModeToBinary();
	writer->SetFileName(filename.c_str());
	writer->Write();
}

void write_column_vectors(std::string filename, std::string header, std::vector<std::vector<double>* >& columns) {
	std::ofstream f;
	f.open (filename.c_str());
	f << header << std::endl;
	for (int i = 0; i < columns[0]->size(); ++i) {
		for (int j = 0; j < columns.size()-1; ++j) {
			f << columns[j]->at(i) << " ";
		}
		f << (*(columns.end()-1))->at(i) << std::endl;
	}
	f.close();
}

}
