#pragma once

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "constants.h"
#include "loads.h"


class Solver
{
public:

	/* Data */
	int Nn;											// Number of nodes
	std::vector<Vector2d> nodes_coo;				// Node coordinates
	int Ne;											// Number of elements
	std::vector<Vector3i> mesh;						// Mesh connectivity
	int Nv;											// Number of vertices per elements

	Loads loads;
	std::vector<Vector2i> neumann_edges;
	std::vector<bool> dirichlet_nodes;

	SparseMatrix<double> global_stiffness;
	std::vector<double> global_internal_force;
	std::vector<double> global_dirichlet_force;
	

	
	/* Constructors */
	Solver();
	~Solver() {}



	/* Functions */
	Vector3d ShapeFunction(const Vector2d& local_coo);
	std::vector<Vector2d> DShapeFunction();

	Matrix3d ElementStiffness(const std::vector<Vector2d>& vertices_coo);
	void GlobalStiffness();

	Vector3d ElementInternalForce(const std::vector<Vector2d>& vertices_coo);
	void GlobalInternalForce();
};