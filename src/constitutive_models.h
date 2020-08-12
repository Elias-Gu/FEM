#pragma once

#include <Eigen/Dense>

#include "constants.h"


struct LinearElasticity
{
	/* Data */
	double mu;
	double lambda;

	/* Constructors */
	LinearElasticity(const double _mu, const double _lambda) : mu(_mu), lambda(_lambda) {}
	~LinearElasticity() {}

	/* Functions */
	Matrix2d P(const Matrix2d& F)
	{
		Matrix2d id = Matrix2d::Identity();
		Matrix2d P = mu * (F + F.transpose() - 2.0 * id) + lambda * (F - id).trace() * id;

		return P;
	}
};