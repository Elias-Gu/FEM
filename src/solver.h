#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

#include "constants.h"
#include "loads.h"


class Solver
{
public:

	/* Data */
	int Nn;													// Number of nodes
	int Ne;													// Number of elements
	int Nv;													// Number of vertices per elements
	std::vector<Vector2d> nodes_coo;						// Node coordinates
	std::vector<Vector3i> mesh;								// Mesh connectivity
	std::vector<std::vector<Vector2i>> incident_elements;	// For each nodes, corresponding position in the incident element vectors

	Loads loads;
	std::vector<Vector2i> neumann_edges;
	std::vector<bool> dirichlet_nodes;

	SparseMatrix<double> global_stiffness;
	SparseMatrix<double> global_mass;
	VectorXd global_internal_force;
	VectorXd global_dirichlet_force;
	VectorXd global_neumann_force;
	VectorXd global_force;

	VectorXd sol_prev;
	VectorXd sol_prev;
	
	bool verbose;
	double time_global_stiffness;
	double time_global_internal_force;
	double time_global_neumann_force;
	double time_global_force;
	double time_solve;


	
	/* Constructors */
	Solver();
	~Solver() {}



	/* Functions */
	Vector3d ShapeFunction(const Vector2d& local_coo);
	Vector2d LineShapeFunction(const double local_coo);
	std::vector<Vector2d> DShapeFunction();

	Matrix3d ElementStiffness(const std::vector<Vector2d>& vertices_coo);
	void GlobalStiffness();

	Vector3d ElementInternalForce(const std::vector<Vector2d>& vertices_coo);
	void GlobalInternalForce();

	Vector2d ElementNeumannForce(const std::vector<Vector2d>& vertices_coo);
	void GlobalNeumannForce();

	void FEMSolver();
};