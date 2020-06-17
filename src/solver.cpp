#include "solver.h"

/* Constructors */
Solver::Solver()
{
	// Triangular mesh
	Nn = n * n;
	Ne = 2 * (n - 1) * (n - 1);
	Nv = 3;

	// Nodes positions
	nodes_coo.resize(Nn);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
		{
			nodes_coo[i * n + j] = Vector2d(i * dx, j * dx);
			nodes_coo[i * n + j][0] += 0.25;
			nodes_coo[i * n + j][1] += 0.25;
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


	// Set boundary nodes/edges
	loads = Loads();
	neumann_edges = loads.NeumannEdges(nodes_coo);
	dirichlet_nodes = loads.DirichletNodes(nodes_coo);


	// Resize global arrays
	global_stiffness.resize(Nn, Nn);
	global_internal_force.resize(Nn);
	global_dirichlet_force.resize(Nn);
}



/* -----------------------------------------------------------------------
|							 SHAPE FUNCTIONS							 |
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



/* -----------------------------------------------------------------------
|							STIFFNESS MATRIX							 |
----------------------------------------------------------------------- */


Matrix3d Solver::ElementStiffness(const std::vector<Vector2d>& vertices_coo)
{
	std::vector<Vector2d> dN = DShapeFunction();

	Matrix2d Dm; Dm << vertices_coo[1][0] - vertices_coo[0][0], vertices_coo[2][0] - vertices_coo[0][0],
					   vertices_coo[1][1] - vertices_coo[0][1], vertices_coo[2][1] - vertices_coo[0][1];
	Matrix2d Dm_inverse = Dm.inverse();
	Matrix2d B = Dm_inverse * Dm_inverse.transpose() * Dm.determinant() / 2.0;
	
	Matrix3d element_stiffness;
	for (int i = 0; i < Nv; i++) {
		for (int j = 0; j < Nv; j++)
		{
			element_stiffness(i, j) = dN[i].transpose() * B * dN[j];
		}
	}

	return element_stiffness;
}


void Solver::GlobalStiffness()
{
	std::vector<Triplet<double>> global_entries;

	for (size_t e = 0; e < Ne; e++)
	{
		// Get coordinates of element vertices
		std::vector<Vector2d> vertices_coo(Nv);
		for (size_t i = 0; i < Nv; i++)
			vertices_coo[i] = nodes_coo[mesh[e][i]];

		Matrix3d element_stiffness = ElementStiffness(vertices_coo);

		for (int i = 0; i < Nv; i++)
		{
			if (!dirichlet_nodes[mesh[e][i]])
			{
				for (int j = 0; j < Nv; j++)
				{
					if (!dirichlet_nodes[mesh[e][j]])	// Add stiffness entry for free nodes
						global_entries.push_back(Triplet<double>(mesh[e][i], mesh[e][j], element_stiffness(i, j)));
					else								// Add force contribution of dirichlet nodes
						global_dirichlet_force[mesh[e][i]] += element_stiffness(i, j) * loads.DirichletValue(nodes_coo[mesh[e][j]]);
				}
			}
		}
	}

	global_stiffness.setFromTriplets(global_entries.begin(), global_entries.end());
}