#include "elastostatic_solver.h"

/* Constructors */
ElastostaticSolver::ElastostaticSolver()
{
	// Triangular mesh
	Nn = n * n;
	Ne = 2 * (n - 1) * (n - 1);
	Nv = 3;

	// Nodes positions
	nodes_coo.resize(Nn);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++)
		{
			nodes_coo[i * n + j] = Vector2d(i * dx, j * dx);
		}
	}


	// Nodes connectivity
	mesh.resize(Ne);
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++)
		{
			int q = (n - 1) * i + j;
			mesh[2 * q] = Eigen::Vector3i((i + 1) * n + j, (i + 1) * n + j + 1, i * n + j);
			mesh[2 * q + 1] = Eigen::Vector3i(i * n + j + 1, i * n + j, (i + 1) * n + j + 1);
		}
	}


	// Incident element
	incident_elements.resize(Nn);
	for (int i = 0; i < Ne; i++)
		for (int j = 0; j < 3; j++)
			incident_elements[mesh[i][j]].push_back(Vector2i(i, j));


	// Set boundary nodes/edges
	// loads = Loads();
	// neumann_edges = loads.NeumannEdges(nodes_coo);
	// dirichlet_nodes = loads.DirichletNodes(nodes_coo);


	// Resize global arrays
	// global_internal_force.setZero(Nn);
	// global_dirichlet_force.setZero(Nn);
	// global_neumann_force.setZero(Nn);
	// global_force.setZero(Nn);
	global_elastic_force.setZero(Nn);
	disp.setZero(Nn);

	verbose = true;
	if (verbose)
	{
		std::cout << "Number of nodes:    " << Nn << std::endl;
		std::cout << "Number of elements: " << Ne << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
	}
}



/* -----------------------------------------------------------------------
|							 SHAPE FUNCTIONS							 |
----------------------------------------------------------------------- */


Vector3d ElastostaticSolver::ShapeFunction(const Vector2d& local_coo)
{
	Vector3d N;
	N[0] = 1 - local_coo[0] - local_coo[1];
	N[1] = local_coo[0];
	N[2] = local_coo[1];
	return N;
}


Vector2d ElastostaticSolver::LineShapeFunction(const double local_coo)
{
	Vector2d Nl;
	Nl[0] = 1 - local_coo;
	Nl[1] = local_coo;
	return Nl;
}


std::vector<Vector2d> ElastostaticSolver::DShapeFunction()
{
	std::vector<Vector2d> dN(3);
	dN[0] = Vector2d(-1.0, -1.0);
	dN[1] = Vector2d(1.0, 0.0);
	dN[2] = Vector2d(0.0, 1.0);
	return dN;
}



/* -----------------------------------------------------------------------
|							  ELASTIC FORCES							 |
----------------------------------------------------------------------- */


Vector3d ElastostaticSolver::ElementElasticForce(const std::vector<Vector2d>& vertices_coo)
{
	std::vector<Vector2d> dN = DShapeFunction();

	Matrix2d Dm; Dm << vertices_coo[1][0] - vertices_coo[0][0], vertices_coo[2][0] - vertices_coo[0][0],
					   vertices_coo[1][1] - vertices_coo[0][1], vertices_coo[2][1] - vertices_coo[0][1];
	Matrix2d Dm_inverse = Dm.inverse();
	double area = Dm.determinant / 2.0;
}

