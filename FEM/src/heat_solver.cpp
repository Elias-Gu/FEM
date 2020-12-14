#include "heat_solver.h"

/* Constructors */
HeatSolver::HeatSolver(const bool _verbose) : verbose(_verbose)
{
	// Initialize material properties
	rho = 1.0;
	capacity = 1.0;
	conduc = 1.0;


	// Initialize time & coefficient.
	alpha = 1.0;
	if (alpha < 0.5)
	{
		std::cout << "Scheme might not be stable for this value of alpha." << std::endl;
		std::cout << "Use alpha > 1/2 for unconditionally stable scheme." << std::endl;
		std::cout << "Press any key to continue." << std::endl;
		std::cin.get();
	}

	tn = 0;

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


	// Initial node value
	val.resize(Nn);
	for (int i = 0; i < Nn; i++)
		val[i] = 4.0;


	// Set boundary nodes/edges
	loads = Loads();
	neumann_edges = loads.NeumannEdges(nodes_coo);
	dirichlet_nodes = loads.DirichletNodes(nodes_coo);


	// Resize global arrays
	global_stiffness.resize(Nn, Nn);
	global_stiffness_reduced.resize(Nn, Nn);
	global_mass.resize(Nn, Nn);
	global_mass_reduced.resize(Nn, Nn);
	global_force.setZero(Nn);
	global_force_np1.setZero(Nn);

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


Vector3d HeatSolver::ShapeFunction(const Vector2d& local_coo)
{
	Vector3d N;
	N[0] = 1 - local_coo[0] - local_coo[1];
	N[1] = local_coo[0];
	N[2] = local_coo[1];
	return N;
}


Vector2d HeatSolver::LineShapeFunction(const double local_coo)
{
	Vector2d Nl;
	Nl[0] = 1 - local_coo;
	Nl[1] = local_coo;
	return Nl;
}


std::vector<Vector2d> HeatSolver::DShapeFunction()
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


Matrix3d HeatSolver::ElementStiffness(const std::vector<Vector2d>& vertices_coo)
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

	return element_stiffness * conduc;
}


void HeatSolver::GlobalStiffness()
{
	// ScopedTimer timer("", &time_global_stiffness, false);

	std::vector<Triplet<double>> global_entries;

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
						global_entries_private.push_back(Triplet<double>(mesh[e][i], mesh[e][j], element_stiffness(i, j)));
					}
				}
			}
		}

		// Combine all the vector of triplets
		#pragma omp critical
		global_entries.insert(global_entries.end(), global_entries_private.begin(), global_entries_private.end());
	}

	global_stiffness.setFromTriplets(global_entries.begin(), global_entries.end());
	global_stiffness.makeCompressed();
}



/* -----------------------------------------------------------------------
|							 MASS MATRIX								 |
----------------------------------------------------------------------- */

Matrix3d HeatSolver::ElementMass(const std::vector<Vector2d>& vertices_coo)
{
	double area = 0.5;
	std::vector<Vector2d> mids_iso = { Vector2d(0.5, 0.0),
									   Vector2d(0.5, 0.5),
									   Vector2d(0.0, 0.5) }; // Mid points in isoparametric element
	std::vector<Vector3d> N_iso = { ShapeFunction(mids_iso[0]),
									ShapeFunction(mids_iso[1]),
									ShapeFunction(mids_iso[2]) };

	Matrix2d Dm; Dm << vertices_coo[1][0] - vertices_coo[0][0], vertices_coo[2][0] - vertices_coo[0][0],
		vertices_coo[1][1] - vertices_coo[0][1], vertices_coo[2][1] - vertices_coo[0][1];
	double Dm_det = Dm.determinant();

	Matrix3d element_mass; element_mass.setZero();
	for (int i = 0; i < Nv; i++)
		for (int j = 0; j < Nv; j++)
			for (int k = 0; k < Nv; k++)
				element_mass(i, j) += N_iso[i][k] * N_iso[j][k];
	element_mass *= (Dm_det * area / 3.0) * rho * capacity;

	return element_mass;
}


void HeatSolver::GlobalMass()
{
	// ScopedTimer timer("", &time_global_mass, false);

	std::vector<Triplet<double>> global_entries;

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

			Matrix3d element_mass = ElementMass(vertices_coo);

			for (int i = 0; i < Nv; i++)
			{
				if (!dirichlet_nodes[mesh[e][i]])
				{
					for (int j = 0; j < Nv; j++)
					{
						global_entries_private.push_back(Triplet<double>(mesh[e][i], mesh[e][j], element_mass(i, j)));
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
	global_mass.setFromTriplets(global_entries.begin(), global_entries.end());
	global_mass.makeCompressed();
}


/* -----------------------------------------------------------------------
|							 INTERNAL FORCE								 |
----------------------------------------------------------------------- */


Vector3d HeatSolver::ElementInternalForce(const std::vector<Vector2d>& vertices_coo, const double _tn)
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

	std::vector<Vector2d> mids_global(Nv);
	for (int i = 0; i < Nv; i++) 
		for (int j = 0; j < Nv; j++)
			mids_global[i] += N_iso[i][j] * vertices_coo[j];

	Matrix2d Dm; Dm << vertices_coo[1][0] - vertices_coo[0][0], vertices_coo[2][0] - vertices_coo[0][0],
		               vertices_coo[1][1] - vertices_coo[0][1], vertices_coo[2][1] - vertices_coo[0][1];
	double Dm_det = Dm.determinant();

	Vector3d element_internal_force; element_internal_force.setZero();
	for (int i = 0; i < Nv; i++)
		element_internal_force += N_iso[i] * loads.InternalForce(mids_global[i], _tn) * area / 3.0;
	element_internal_force *= Dm_det;

	return element_internal_force;
}


VectorXd HeatSolver::GlobalInternalForce(const double _tn)
{
	VectorXd global_internal_force(Nn); global_internal_force.setZero();
	std::vector<Vector3d> elements_internal_force(Ne);

	#pragma omp parallel for
	for (int e = 0; e < Ne; e++)
	{
		// Get coordinates of element vertices
		std::vector<Vector2d> vertices_coo(Nv);
		for (size_t i = 0; i < Nv; i++)
			vertices_coo[i] = nodes_coo[mesh[e][i]];

		elements_internal_force[e] = ElementInternalForce(vertices_coo, _tn);
	}

	#pragma omp parallel for
	for (int i = 0; i < Nn; i++)
		for (size_t j = 0; j < incident_elements[i].size(); j++)
			global_internal_force[i] += elements_internal_force[incident_elements[i][j][0]][incident_elements[i][j][1]];

	return global_internal_force;
}



/* -----------------------------------------------------------------------
|							 BOUNDARY FORCES							 |
----------------------------------------------------------------------- */


Vector2d HeatSolver::ElementNeumannForce(const std::vector<Vector2d>& vertices_coo, const double _tn)
{

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
		element_neumann_force += weights[i] * N_iso_line[i] * loads.NeumannForce(pts_global[i], normal, _tn);
	element_neumann_force *= length;

	return element_neumann_force;
}


VectorXd HeatSolver::GlobalNeumannForce(const double _tn)
{
	VectorXd global_neumann_force(Nn); global_neumann_force.setZero();

	for (size_t e = 0; e < neumann_edges.size(); e++)
	{
		// Get coordinates of edge vertices
		std::vector<Vector2d> vertices_coo(2);
		for (size_t i = 0; i < 2; i++)
			vertices_coo[i] = nodes_coo[neumann_edges[e][i]];

		Vector2d element_neumann_force = ElementNeumannForce(vertices_coo, _tn);

		for (int i = 0; i < 2; i++)
			global_neumann_force[neumann_edges[e][i]] += element_neumann_force[i];
	}

	return global_neumann_force;
}



/* -----------------------------------------------------------------------
|		   					   TOTAL FORCE								 |
----------------------------------------------------------------------- */


void HeatSolver::TotalForce(const double _tn)
{
	// ScopedTimer timer("", &time_global_force, false);

	// Build internal and neumann force vectors at tn and tn+1
	VectorXd global_internal_force = GlobalInternalForce(_tn);
	VectorXd global_neumann_force = GlobalNeumannForce(_tn);

	// Build dirichlet force vectors at tn and tn+1
	VectorXd dir_values(Nn), dir_dtvalues(Nn);
	dir_values.setZero(), dir_dtvalues.setZero();

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Nn; i++)
	{
		if (dirichlet_nodes[i])
		{
			dir_values[i] = loads.DirichletValue(nodes_coo[i], _tn);
			dir_dtvalues[i] = loads.DtDirichletValue(nodes_coo[i], _tn);
		}
	}

	VectorXd global_dirichlet_force = global_stiffness * dir_values + global_mass * dir_dtvalues;

	// Update newest force vector
	global_force_np1 = global_internal_force + global_neumann_force - global_dirichlet_force;
}



/* -----------------------------------------------------------------------
|								SOLVE SYSTEM							 |
----------------------------------------------------------------------- */


void HeatSolver::Init()
{
	// Build semi-reduced matrices
	GlobalStiffness();
	GlobalMass();

	// Build reduced matrices
	global_stiffness_reduced = global_stiffness;
	global_mass_reduced = global_mass;
	for (int i = 0; i < Nn; i++)
		if (dirichlet_nodes[i])
		{
			global_stiffness_reduced.col(i) *= 0;
			global_mass_reduced.col(i) *= 0;
			global_mass_reduced.coeffRef(i, i) += 1.0;
		}
	global_stiffness_reduced.makeCompressed();
	global_mass_reduced.makeCompressed();

	// Build first force vector
	TotalForce(tn - dt); // -dt to be consistant with initial val[]
}


void HeatSolver::Step()
{
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

	// New force state
	global_force = global_force_np1;
	TotalForce(tn);					// Update global_force_np1 value

	// Build system

	SparseMatrix<double, RowMajor> LHS = global_mass_reduced + alpha * dt * global_stiffness_reduced;
	VectorXd RHS = (global_mass_reduced - (1 - alpha) * dt * global_stiffness_reduced) * val + dt * (alpha * global_force_np1 + (1.0 - alpha) * global_force);

	// Enforce dirichlet boundary conditions
	#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < Nn; i++)
		if (dirichlet_nodes[i])
			RHS[i] = loads.DirichletValue(nodes_coo[i], tn);

	// Solve linear system
	ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
	solver.setTolerance(1e-10);					// Can be lowered
	val = solver.compute(LHS).solve(RHS);
	if (solver.info() != Success)
	{
		std::cout << "Conjugate Gradient solver failed" << std::endl;
		return;
	}

	// Timer and output
	auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = stop - start;
	time_system = (double)(elapsed.count() * 1000.0);

	if (verbose)
	{
		std::cout << "Sim time: " << std::setfill(' ') << std::setw(9) << std::fixed << std::setprecision(3) << tn << " s     ||     Comp time" << std::setfill(' ') << std::setw(9) << std::fixed << std::setprecision(3) << time_system + time_global_force << " ms" << std::endl;
	}

	// Update time
	tn += dt;
}