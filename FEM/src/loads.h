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
	double InternalForce(const Vector2d& coo, const double tn);

	std::vector<Vector2i> NeumannEdges(const std::vector<Vector2d>& coo);
	double NeumannForce(const Vector2d& coo, const Vector2d& normal);
	double NeumannForce(const Vector2d& coo, const Vector2d& normal, const double tn);

	std::vector<bool> DirichletNodes(const std::vector<Vector2d>& coo);
	double DirichletValue(const Vector2d& coo);
	double DirichletValue(const Vector2d& coo, const double tn);
	double DtDirichletValue(const Vector2d& coo, const double tn);
};