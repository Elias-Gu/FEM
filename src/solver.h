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
	std::vector<Vector2d> nodes_coo;			// Node coordinates
	std::vector<Vector3i> mesh;				// Mesh connectivity
	int Ne;									// Number of elements

	Loads loads;
	std::vector<Vector2i> neumann_edges;
	std::vector<bool> dirichlet_nodes;




	/* Constructors */
	Solver();
	~Solver() {}



	/* Functions */
	Vector3d ShapeFunction(const Vector2d& local_coo);
	std::vector<Vector2d> DShapeFunction();
	Matrix2d ElementStiffness(const std::vector<Vector2d>& vertices_coo);
	SparseMatrix<double> GlobalStiffness();
};