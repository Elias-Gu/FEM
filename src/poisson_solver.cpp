#include "poisson_solver.h"

/* Constructors */
PoissonSolver::PoissonSolver()
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
			//nodes_coo[i * n + j][0] += 0.25;
			//nodes_coo[i * n + j][1] += 0.25;
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
	loads = Loads();
	neumann_edges = loads.NeumannEdges(nodes_coo);
	dirichlet_nodes = loads.DirichletNodes(nodes_coo);


	// Resize global arrays
	global_stiffness.resize(Nn, Nn);
	global_internal_force.setZero(Nn);
	global_dirichlet_force.setZero(Nn);
	global_neumann_force.setZero(Nn);
	global_force.setZero(Nn);
	sol.resize(Nn);

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


Vector3d PoissonSolver::ShapeFunction(const Vector2d& local_coo)
{
	Vector3d N;
	N[0] = 1 - local_coo[0] - local_coo[1];
	N[1] = local_coo[0];
	N[2] = local_coo[1];
	return N;
}


Vector2d PoissonSolver::LineShapeFunction(const double local_coo)
{
	Vector2d Nl;
	Nl[0] = 1 - local_coo;
	Nl[1] = local_coo;
	return Nl;
}


std::vector<Vector2d> PoissonSolver::DShapeFunction()
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


Matrix3d PoissonSolver::ElementStiffness(const std::vector<Vector2d>& vertices_coo)
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


void PoissonSolver::GlobalStiffness()
{
	auto start = std::chrono::high_resolution_clock::now();
	if (verbose)
		std::cout << "GlobalStiffness ... ";

	std::vector<Triplet<double>> global_entries;
	std::vector<Vector3d> elements_dirichlet_force(Ne);

	// Each thread creates its own vector of triplets
	#pragma omp parallel
	{
		std::vector<Triplet<double>> global_entries_private;

		#pragma omp for nowait
		for (int e = 0; e < Ne; e++)
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
							global_entries_private.push_back(Triplet<double>(mesh[e][i], mesh[e][j], element_stiffness(i, j)));
						else								// Add force contribution of dirichlet nodes
							elements_dirichlet_force[e][i] += element_stiffness(i, j) * loads.DirichletValue(nodes_coo[mesh[e][j]]);
					}
				}
				else
					global_entries_private.push_back(Triplet<double>(mesh[e][i], mesh[e][i], 1.0 / incident_elements[mesh[e][i]].size()));
			}
		}

		// Combine all the vector of triplets
		#pragma omp critical
		global_entries.insert(global_entries.end(), global_entries_private.begin(), global_entries_private.end());
	}

	// Form global stiffness matrix from list of triplets
	global_stiffness.setFromTriplets(global_entries.begin(), global_entries.end());

	// Form global force created by dirichlet nodes
	#pragma omp parallel for
	for (int i = 0; i < Nn; i++)
		for (size_t j = 0; j < incident_elements[i].size(); j++)
			global_dirichlet_force[i] += elements_dirichlet_force[incident_elements[i][j][0]][incident_elements[i][j][1]];

	if (verbose)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		time_global_stiffness = (double)(elapsed.count() * 1000.0);
		std::cout << "      DONE";
		std::cout << " -- " << time_global_stiffness << " ms" << std::endl;
	}
}



/* -----------------------------------------------------------------------
|							 INTERNAL FORCE								 |
----------------------------------------------------------------------- */


Vector3d PoissonSolver::ElementInternalForce(const std::vector<Vector2d>& vertices_coo)
{
	// See 12.2 in Numerical Solution of Partial Differential Equations by the FEM, by C.Johnson
	// TODO: time and see if quadrature of degree 1 is much better
	double area = 0.5;
	std::vector<Vector2d> mids_iso = { Vector2d(0.5, 0.0), 
									   Vector2d(0.5, 0.5), 
									   Vector2d(0.0, 0.5) }; // Mid points in isoparametric element
	std::vector<Vector3d> N_iso = { ShapeFunction(mids_iso[0]),
									ShapeFunction(mids_iso[1]),
									ShapeFunction(mids_iso[2]) };
	// TODO: since we just need this, might be faster to get them from vertices_coo directly
	std::vector<Vector2d> mids_global(Nv);
	for (int i = 0; i < Nv; i++) 
	{	//TODO: Check that setZero is actually doing what it is supposed to do
		mids_global[i].setZero();

		for (int j = 0; j < Nv; j++)
			mids_global[i] += N_iso[i][j] * vertices_coo[j];
	}

	Matrix2d Dm; Dm << vertices_coo[1][0] - vertices_coo[0][0], vertices_coo[2][0] - vertices_coo[0][0],
		               vertices_coo[1][1] - vertices_coo[0][1], vertices_coo[2][1] - vertices_coo[0][1];
	double Dm_det = Dm.determinant();

	Vector3d element_internal_force; element_internal_force.setZero();
	for (int i = 0; i < Nv; i++)
		element_internal_force += N_iso[i] * loads.InternalForce(mids_global[i]) * area / 3.0;
	element_internal_force *= Dm_det;

	return element_internal_force;
}


void PoissonSolver::GlobalInternalForce()
{
	auto start = std::chrono::high_resolution_clock::now();
	if (verbose)
		std::cout << "GlobalInternalForce ... ";

	std::vector<Vector3d> elements_internal_force(Ne);

	#pragma omp parallel for
	for (int e = 0; e < Ne; e++)
	{
		// Get coordinates of element vertices
		std::vector<Vector2d> vertices_coo(Nv);
		for (size_t i = 0; i < Nv; i++)
			vertices_coo[i] = nodes_coo[mesh[e][i]];

		elements_internal_force[e] = ElementInternalForce(vertices_coo);
	}

	#pragma omp parallel for
	for (int i = 0; i < Nn; i++)
		for (size_t j = 0; j < incident_elements[i].size(); j++)
			global_internal_force[i] += elements_internal_force[incident_elements[i][j][0]][incident_elements[i][j][1]];

	if (verbose)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		time_global_internal_force = (double)(elapsed.count() * 1000.0);
		std::cout << "  DONE";
		std::cout << " -- " << time_global_internal_force << " ms" << std::endl;
	}
}



/* -----------------------------------------------------------------------
|							 BOUNDARY FORCES							 |
----------------------------------------------------------------------- */


Vector2d PoissonSolver::ElementNeumannForce(const std::vector<Vector2d>& vertices_coo)
{

	// TODO: time and see if quadrature of degree 1 is much better
	double length = Vector2d(vertices_coo[1] - vertices_coo[0]).norm();
	std::vector<double> pts_iso = { -1 / sqrt(3.0), 1 / sqrt(3.0) }; // Integration points on [-1,1]
	for (auto& pt : pts_iso) 
		pt = 0.5 * pt + 0.5;										 // Integration points on [0,1] (iso p element)
	std::vector<double> weights = { 0.5, 0.5 };
	std::vector<Vector2d> N_iso_line = { LineShapeFunction(pts_iso[0]), 
										 LineShapeFunction(pts_iso[1]) };
	std::vector<Vector2d> pts_global(2);
	for (int i = 0; i < 2; i++)
	{
		pts_global[i].setZero();

		for (int j = 0; j < 2; j++)
			pts_global[i] += N_iso_line[i][j] * vertices_coo[j];
	}

	Vector2d normal = Vector2d(-vertices_coo[1][1] + vertices_coo[0][1], vertices_coo[1][0] - vertices_coo[0][0]).normalized();

	Vector2d element_neumann_force; element_neumann_force.setZero();
	for (int i = 0; i < 2; i++)
		element_neumann_force += weights[i] * N_iso_line[i] * loads.NeumannForce(pts_global[i], normal);
	element_neumann_force *= length;

	return element_neumann_force;
}


void PoissonSolver::GlobalNeumannForce()
{
	auto start = std::chrono::high_resolution_clock::now();
	if (verbose)
		std::cout << "GlobalNeumannForce ... ";

	for (size_t e = 0; e < neumann_edges.size(); e++)
	{
		// Get coordinates of edge vertices
		std::vector<Vector2d> vertices_coo(2);
		for (size_t i = 0; i < 2; i++)
			vertices_coo[i] = nodes_coo[neumann_edges[e][i]];

		Vector2d element_neumann_force = ElementNeumannForce(vertices_coo);

		for (int i = 0; i < 2; i++)
			global_neumann_force[neumann_edges[e][i]] += element_neumann_force[i];
	}

	if (verbose)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		time_global_neumann_force = (double)(elapsed.count() * 1000.0);
		std::cout << "   DONE";
		std::cout << " -- " << time_global_neumann_force << " ms" << std::endl;
	}
}



/* -----------------------------------------------------------------------
|								SOLVE SYSTEM							 |
----------------------------------------------------------------------- */


void PoissonSolver::FEMSolver()
{
	GlobalStiffness();
	GlobalInternalForce();
	GlobalNeumannForce();

	auto start = std::chrono::high_resolution_clock::now();
	if (verbose)
		std::cout << "GlobalForce ... ";

	VectorXd global_force = global_internal_force - global_neumann_force;

	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < Nn; i++)
	{
		if (dirichlet_nodes[i])
		{
			global_dirichlet_force[i] = -loads.DirichletValue(nodes_coo[i]);
			global_force[i] = 0.0;
		}
	}
	global_force -= global_dirichlet_force;
	
	if (verbose)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		time_global_force = (double)(elapsed.count() * 1000.0);
		std::cout << "          DONE";
		std::cout << " -- " << time_global_force << " ms" << std::endl;
		std::cout << "Solving linear system ... ";
	}

	ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
	cg.compute(global_stiffness);
	sol = cg.solve(global_force);

	if (verbose)
	{
		auto stop = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = stop - start;
		time_solve = (double)(elapsed.count() * 1000.0);
		std::cout << "DONE";
		std::cout << " -- " << time_solve << " ms" << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		double total_time = time_global_stiffness + time_global_internal_force + time_global_neumann_force + time_global_force + time_solve;
		std::cout << std::fixed << std::setprecision(2);
		std::cout << "GlobalStiffness     -- " << time_global_stiffness / total_time * 100 << "%" << std::endl;
		std::cout << "GlobalInternalForce -- " << time_global_internal_force / total_time * 100 << "%" << std::endl;
		std::cout << "GlobalNeumannForce  -- " << time_global_neumann_force / total_time * 100 << "%" << std::endl;
		std::cout << "GlobalForce         -- " << time_global_force / total_time * 100 << "%" << std::endl;
		std::cout << "FEMSolver           -- " << time_solve / total_time * 100 << "%" << std::endl;
	}
}