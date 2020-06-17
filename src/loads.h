#pragma once

#include <vector>

#include <Eigen/Dense>

#include "constants.h"


struct Loads
{

	/* Constructors */
	Loads() {};
	~Loads() {};



	/* Functions */
	double InternalForce(const Vector2d& coo);

	std::vector<Vector2i> NeumannEdges(const std::vector<Vector2d>& coo);
	double NeumannForce(const Vector2d& coo, const Vector2d& normal);

	std::vector<bool> DirichletNodes(const std::vector<Vector2d>& coo);
	double DirichletValue(const Vector2d& coo);
};