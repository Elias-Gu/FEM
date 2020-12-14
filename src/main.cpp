#include <vector>
#include <iostream>

#include <Eigen/Dense>

#include "constants.h"
#include "helpers.h"
#include "heat_solver.h"
#include "poisson_solver.h"


/* Declarations */
HeatSolver* heat_solver;
PoissonSolver* poisson_solver;



/* -----------------------------------------------------------------------
|								MAIN									 |
----------------------------------------------------------------------- */


int main()
{	
	Output output;

	/* Run poisson solver */
	poisson_solver = new PoissonSolver();
	poisson_solver->FEMSolver();
	output.OutputNodeCoordinates(poisson_solver->nodes_coo, poisson_solver->sol);


	/* Run heat equation solver */
	heat_solver = new HeatSolver();
	heat_solver->Init();

	int step_counter = 0;
	double final_time = 8.0;
	output.OutputNodeCoordinates(heat_solver->nodes_coo, heat_solver->val, step_counter);
	while (heat_solver->tn < final_time)
	{
		heat_solver->Step();
		step_counter++;
		output.OutputNodeCoordinates(heat_solver->nodes_coo, heat_solver->val, step_counter);
	}

	std::cin.get();

	delete poisson_solver;
	delete heat_solver;

	return 0;
}