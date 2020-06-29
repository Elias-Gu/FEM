#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "constants.h"
#include "helpers.h"
#include "loads.h"

/*
References:
"The Finite Element Method" by Thomas J.R. Hugues
"Numerical Solution of Partial Differential Equations by the Finite Element Method" by Claes Johnson
*/

class Solver
{
public:

	/* Data */
	double rho;												// Density
	double capacity;										// Heat capacity
	double conduc;											// Thermal conductivity

	double alpha;											// Coefficient for generalized trapezoidal method, in [0.1]
	double tn;
	int Nn;													// Number of nodes
	int Ne;													// Number of elements
	int Nv;													// Number of vertices per elements
	std::vector<Vector2d> nodes_coo;						// Node coordinates
	std::vector<Vector3i> mesh;								// Mesh connectivity
	std::vector<std::vector<Vector2i>> incident_elements;	// For each nodes, corresponding position in the incident element vectors

	Loads loads;
	std::vector<Vector2i> neumann_edges;
	std::vector<bool> dirichlet_nodes;

	SparseMatrix<double, RowMajor> global_stiffness;		// RowMajor for paralle matrix-vector product
	SparseMatrix<double, RowMajor> global_stiffness_reduced;// Might be faster with ColumnMajor though
	SparseMatrix<double, RowMajor> global_mass;
	SparseMatrix<double, RowMajor> global_mass_reduced;
	VectorXd global_force;
	VectorXd global_force_np1;

	VectorXd val;
	
	bool verbose;
	double time_global_stiffness;
	double time_global_mass;
	double time_global_force;
	double time_system;


	
	/* Constructors */
	Solver();
	~Solver() {}



	/* Functions */
	Vector3d ShapeFunction(const Vector2d& local_coo);
	Vector2d LineShapeFunction(const double local_coo);
	std::vector<Vector2d> DShapeFunction();

	Matrix3d ElementStiffness(const std::vector<Vector2d>& vertices_coo);
	void GlobalStiffness();

	Matrix3d ElementMass(const std::vector<Vector2d>& vertices_coo);
	void GlobalMass();

	Vector3d ElementInternalForce(const std::vector<Vector2d>& vertices_coo, const double _tn);
	VectorXd GlobalInternalForce(const double _tn);

	Vector2d ElementNeumannForce(const std::vector<Vector2d>& vertices_coo, const double _tn);
	VectorXd GlobalNeumannForce(const double _tn);

	void Init();
	void TotalForce(const double _tn);
	void Step();
};