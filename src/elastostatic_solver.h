#pragma once

#include <iostream>

#include <Eigen/Dense>

#include "constants.h"
#include "loads.h"


class ElastostaticSolver
{
public:

	/* Data */
	int Nn;													// Number of nodes
	int Ne;													// Number of elements
	int Nv;													// Number of vertices per elements
	std::vector<Vector2d> nodes_coo;						// Node coordinates
	std::vector<Vector3i> mesh;								// Mesh connectivity
	std::vector<std::vector<Vector2i>> incident_elements;	// For each nodes, corresponding position in the incident element vectors

	// Loads loads;
	// std::vector<Vector2i> neumann_edges;
	// std::vector<bool> dirichlet_nodes;

	VectorXd global_elastic_force;

	VectorXd disp;

	bool verbose;
	//double time_global_stiffness;
	//double time_global_internal_force;
	//double time_global_neumann_force;
	//double time_global_force;
	//double time_solve;



	/* Constructors */
	ElastostaticSolver();
	~ElastostaticSolver() {}



	/* Functions */
	Vector3d ShapeFunction(const Vector2d& local_coo);
	Vector2d LineShapeFunction(const double local_coo);
	std::vector<Vector2d> DShapeFunction();

	Vector3d ElementElasticForce(const std::vector<Vector2d>& vertices_coo);
	void GlobalElasticForce();
	//Matrix3d ElementStiffness(const std::vector<Vector2d>& vertices_coo);
	//void GlobalStiffness();

	//Vector3d ElementInternalForce(const std::vector<Vector2d>& vertices_coo);
	//void GlobalInternalForce();

	//Vector2d ElementNeumannForce(const std::vector<Vector2d>& vertices_coo);
	//void GlobalNeumannForce();

	//void FEMSolver();
};