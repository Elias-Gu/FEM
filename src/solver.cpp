#include "solver.h"

/* Constructors */
Solver::Solver()
{
	Ne = 2 * (n - 1) * (n - 1);

	// node positions
	nodes_coo.resize(n * n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			nodes_coo[i * n + j] = Vector2d(i * dx, j * dx);
			nodes_coo[i * n + j][0] += 0.25;
			nodes_coo[i * n + j][1] += 0.25;
		}


	// node connectivity
	mesh.resize(Ne);
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			int q = (n - 1) * i + j;
			mesh[2 * q] = Eigen::Vector3i((i + 1) * n + j, (i + 1) * n + j + 1, i * n + j);
			mesh[2 * q + 1] = Eigen::Vector3i(i * n + j + 1, i * n + j, (i + 1) * n + j + 1);
		}
	}


	// Set boundary nodes/edges
	loads = Loads();
	neumann_edges = loads.NeumannEdges(nodes_coo);
	dirichlet_nodes = loads.DirichletNodes(nodes_coo);
}



/* -----------------------------------------------------------------------
|						 STIFFNESS MATRIX								 |
----------------------------------------------------------------------- */


Vector3d Solver::ShapeFunction(const Vector2d& local_coo)
{
	Vector3d N;
	N[0] = 1 - local_coo[0] - local_coo[1];
	N[1] = local_coo[0];
	N[2] = local_coo[1];
	return N;
}


std::vector<Vector2d> Solver::DShapeFunction()
{
	std::vector<Vector2d> dN(3);
	dN[0] = Vector2d(-1.0, -1.0);
	dN[1] = Vector2d(1.0, 0.0);
	dN[2] = Vector2d(0.0, 1.0);
	return dN;
}